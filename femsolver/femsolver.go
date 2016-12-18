package femsolver

import (
  "fmt"
  "math"
  "github.com/gonum/matrix/mat64"
)

var DEBUG = true

type FEMsolver interface {
  AddBodyForce(func (float64) float64)
  CalcK()
  CalcLocK()
  Solve()
}

// Build a vector from nod and val
// nod: an array of node number which has value
// val: value of the node
func BuildVec(nod []int, val []float64, vec *mat64.Vector) {
  for i,n := range nod {
    vec.SetVec(n, val[i])
  }
}

// Calculate Gaussian Quadracture of f from a to b with num Gaussian point
func GausQuad(f func(float64) float64, a, b float64, num int) float64 {
  pt, w := getGausQuadPoint(num)
  f = changeInterval(f, a, b)
  sum := 0.0
  for i := 0; i < num; i++ {
    sum += w[i]*f(pt[i])
  }
  return sum
}

// Get the value of Gaussian points and the weight. Currently maximun num is 4.
func getGausQuadPoint(num int) (pt, w []float64) {
  pt = make([]float64, num)
  w = make([]float64, num)
  switch num {
  case 1:
    pt = []float64{0}
    w = []float64{0}
  case 2:
    pt = []float64{-0.5773502691896257, +0.5773502691896257}
    w = []float64{1.0, 1.0}
  case 3:
    pt = []float64{0, -0.7745966692414834, 0.7745966692414834}
    w = []float64{0.8888888888888888, 0.5555555555555556, 0.5555555555555556}
  case 4:
    pt = []float64{-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526}
    w = []float64{0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538}
  default:
    panic("Num of gaussian point should less than 4.")
  }
  return pt, w
}

// Change the interval from a, b to -1, 1
func changeInterval(f func (float64) float64, a, b float64) func (float64) float64 {
  return func (x float64) float64 {
    return (b-a)/2 * f( (b-a)*x/2 + (a+b)/2  )
  }
}

// Linear-element FEM solver for bar
// k: local stifness matrix
// K: global stifness matrix
// Ne: amount of element
// No: order of shape function
// Ng: amount of gaussian points
// Le: lenth of each element
// E: Young's modulus of each element
// A: area of each element
// u: displacement vector
// f: force vector
// uNod, uVal: displacement boundary condition
// fNod, fVal: force boundary condition
type FEMsolver1dBar struct {
  k, K *mat64.SymDense
  Ne, No, Ng int
  Le, E, A, u, f *mat64.Vector
  uNod, fNod []int
  uVal, fVal []float64
}

// Create a FEMsolver1dBar
func NewFEMsolver1dBar(No, Ne int, Le, E, A, u, f *mat64.Vector, uNod, fNod []int, uVal, fVal []float64) *FEMsolver1dBar {
  k := mat64.NewSymDense(No, nil)
  K := mat64.NewSymDense((No-1)*Ne+1, nil)
  BuildVec(uNod, uVal, u)
  fNod, fVal = InitF1dBar((No-1)*Ne+1, uNod, fNod, fVal)
  BuildVec(fNod, fVal, f)
  Ng := (No+1)/2 + 1
  return &FEMsolver1dBar{k, K, Ne, No, Ng, Le, E, A, u, f, uNod, fNod, uVal, fVal}
}

// Create a FEMsolver1dBar with const Le, E and A
func NewFEMsolver1dBarConstLeEA(No, Ne int, Le, E, A float64, u, f *mat64.Vector, uNod, fNod []int, uVal, fVal []float64) *FEMsolver1dBar {

  LeV := mat64.NewVector(Ne+1, nil)
  EV := mat64.NewVector(Ne+1, nil)
  AV := mat64.NewVector(Ne+1, nil)
  for i := 0; i < Ne+1; i++ {
    LeV.SetVec(i, Le)
    EV.SetVec(i, E)
    AV.SetVec(i, A)
  }

  return NewFEMsolver1dBar(No, Ne, LeV, EV, AV, u, f, uNod, fNod, uVal, fVal)
}

// Initialize the force boundary condition
func InitF1dBar (Nn int, uNod, fNod []int, fVal []float64) ([]int, []float64) {
  uNum := len(uNod)
  fNum := len(fNod)
  nfNum := Nn-uNum
  nfNod := make([]int, nfNum)
  nfVal := make([]float64, nfNum)
  c, j, k := 0, 0, 0
  for i := 0; i < Nn; i++ {
    if j < uNum && uNod[j] == i {
      j++
    } else if c < nfNum {
      nfNod[c] = i
      if k < fNum && fNod[k] == i {
        nfVal[c] = fVal[k]
        k++
      }
      c++
    }
  }
  return nfNod, nfVal
}

// Return shape function of e-th element
func (fem *FEMsolver1dBar) NElem(j, e int) func(float64) float64 {
  Le := fem.Le.At(e, 0)
  pLe := Le / float64(fem.No-1)
  xe := make([]float64, fem.No)
  for i := 0; i < fem.No; i++ {
    xe[i] = Le * float64(e) + pLe * float64(i)
  }
  return func(x float64) float64 {
    N := 1.0
    for k := 0; k < fem.No; k++ {
      if k != j {
        N *= (x-xe[k])/(xe[j]-xe[k])
      }
    }
    return N
  }
}

// Return derivative shape function of e-th element
func (fem *FEMsolver1dBar) BElem(j, e int) func(float64) float64 {
  Le := fem.Le.At(e, 0)
  pLe := Le / float64(fem.No-1)
  xe := make([]float64, fem.No)
  for i := 0; i < fem.No; i++ {
    xe[i] = Le * float64(e) + pLe * float64(i)
  }
  return func(x float64) float64 {
    N := 0.0
    for k := 0; k < fem.No; k++ {
      if k != j {
        Np := 1.0
        for l := 0; l < fem.No; l++ {
          if l != k && l != j {
            Np *= (x-xe[l])/(xe[j]-xe[l])
          }
        }
        N += Np * 1/(xe[j]-xe[k])
      }
    }
    return N
  }
}

// Add body force (or distributed force) b to the solver
func (fem *FEMsolver1dBar) AddBodyForce(b func(float64) float64) {
  fmt.Println("Add body force to force vector ...")
  for i := 0; i < fem.Ne; i++ {
    for j := 0; j < fem.No; j++ {
      Nj := fem.NElem(j, i)
      fj := GausQuad(func(x float64) float64 {return Nj(x)*b(x)},
              float64(i)*fem.Le.At(i, 0),
              float64(i+1)*fem.Le.At(i, 0), fem.Ng)
      for k, v := range fem.fNod {
        if v == (fem.No-1)*i+j {
          fem.fVal[k] += fj
        }
      }
    }
    BuildVec(fem.fNod, fem.fVal, fem.f)
  }
  if DEBUG {
    fmt.Println(fem.fNod)
    fmt.Println(fem.fVal)
    fmt.Printf("f = %0.4v\n", mat64.Formatted(fem.f, mat64.Prefix("    ")))
    fmt.Println()
  }
}

// Calculate local k matrix
func (fem *FEMsolver1dBar) CalcLocK() {
  fmt.Println("Calc local stifness matrix ...")
  for i := 0; i < fem.No; i++ {
    for j := 0; j < fem.No; j++ {
      if (i >= j) {
        ke := GausQuad(func(x float64) float64 { return fem.BElem(i, 0)(x)*fem.BElem(j, 0)(x) }, 0, fem.Le.At(0, 0), fem.Ng)
        fem.k.SetSym(i, j, ke)
      }
    }
  }
  if DEBUG {
    fmt.Printf("k = %0.4v\n", mat64.Formatted(fem.k, mat64.Prefix("    ")))
    fmt.Println()
  }
}

// Calculate the global K matrix
func (fem *FEMsolver1dBar) CalcK() {
  fmt.Println("Calc stifness matrix ...")
  for k := 0; k < fem.Ne; k++ {
    var kloc mat64.SymDense
    EAL := fem.E.At(k,0)*fem.A.At(k,0)/fem.Le.At(k,0)
    kp := k*(fem.No-1)
    kloc.ScaleSym(EAL, fem.k)
    klocR, klocC := kloc.Dims()
    for i := 0; i < klocR; i++ {
      for j := 0; j < klocC; j++ {
        if j >= i {
          glo := fem.K.At(kp+i, kp+j)
          loc := kloc.At(i, j)
          fem.K.SetSym(kp+i, kp+j, glo+loc)
        }
      }
    }
  }
  if DEBUG {
    fmt.Printf("K = %0.4v\n", mat64.Formatted(fem.K, mat64.Prefix("    ")))
    fmt.Println()
  }
}

// Solve the problem
func (fem *FEMsolver1dBar) Solve() {
  fmt.Println("Solving ...")

  fNum := len(fem.fNod)
  m := mat64.NewDense(fNum, fNum+1, nil)
  for i := 0; i < fNum; i++ {
    for j := 0; j < fNum; j++ {
      m.Set(i, j, fem.K.At(fem.fNod[i], fem.fNod[j]))
    }
    m.Set(i, fNum, fem.fVal[i])
  }
  if DEBUG {
    fmt.Printf("m = %0.4v\n", mat64.Formatted(m, mat64.Prefix("    ")))
    fmt.Println()
  }

  for i := 0; i < fNum; i++ {
    for j := i+1; j < fNum; j++ {
      c := m.At(j,i)/m.At(i,i)
      row := mat64.Row(nil, j, m)
      for k,v := range row {
        row[k] = v - c*m.At(i,k)
      }
      m.SetRow(j, row)
    }
  }
  if DEBUG {
    fmt.Printf("m = %0.4v\n", mat64.Formatted(m, mat64.Prefix("    ")))
    fmt.Println()
  }

  for i := fNum-1; i >= 0; i-- {
    sum := 0.0
    for j := i+1; j < fNum; j++ {
      sum += m.At(i,j)*fem.u.At(fem.fNod[j],0)
    }
    fem.u.SetVec(fem.fNod[i], (m.At(i, fNum) - sum)/m.At(i,i))
  }
  fmt.Printf("u = %0.4v\n", mat64.Formatted(fem.u, mat64.Prefix("    ")))
  fmt.Println()

  uNum := len(fem.uNod)
  for i := 0; i < uNum; i++ {
    fi := fem.uNod[i]
    for j := 0; j < fem.Ne+1; j++ {
      fv := fem.f.At(fi, 0)+fem.K.At(fi, j)*fem.u.At(j, 0)
      fem.f.SetVec(fi, fv)
    }
  }
  fmt.Printf("f = %0.4v\n", mat64.Formatted(fem.f, mat64.Prefix("    ")))
  fmt.Println()
}

// Linear-element FEM solver for beam
// k: local stifness matrix
// K: global stifness matrix
// Ne: amount of element
// No: order of shape function
// Ng: amount of gaussian points
// Le: lenth of each element
// E: Young's modulus of each element
// A: area of each element
// d: displacement vector
// f: force vector
// dNod, dVal: displacement boundary condition
// fNod, fVal: force boundary condition
type FEMsolver1dBeam struct {
  k, K *mat64.SymDense
  Ne, No, Ng int
  Le, E, I, d, f *mat64.Vector
  dNod, fNod []int
  dVal, fVal []float64
}

// Create a FEMsolver1dBeam
func NewFEMsolver1dBeam(No, Ne int, Le, E, I, d, f *mat64.Vector, dNod, fNod []int, dVal, fVal []float64) *FEMsolver1dBeam {
  k := mat64.NewSymDense(No*2, nil)
  K := mat64.NewSymDense((No-1)*2*Ne+2, nil)
  BuildVec(dNod, dVal, d)
  fNod, fVal = InitF1dBeam((No-1)*2*Ne+2, dNod, fNod, fVal)
  BuildVec(fNod, fVal, f)
  Ng := (No+1)/2 + 1
  return &FEMsolver1dBeam{k, K, Ne, No, Ng, Le, E, I, d, f, dNod, fNod, dVal, fVal}
}

// Create a FEMsolver1dBeam with const Le, E and I
func NewFEMsolver1dBeamConstLeEI(No, Ne int, Le, E, I float64, d, f *mat64.Vector, dNod, fNod []int, dVal, fVal []float64) *FEMsolver1dBeam {

  LeV := mat64.NewVector(Ne+1, nil)
  EV := mat64.NewVector(Ne+1, nil)
  IV := mat64.NewVector(Ne+1, nil)
  for i := 0; i < Ne+1; i++ {
    LeV.SetVec(i, Le)
    EV.SetVec(i, E)
    IV.SetVec(i, I)
  }

  return NewFEMsolver1dBeam(No, Ne, LeV, EV, IV, d, f, dNod, fNod, dVal, fVal)
}

// Initialize the force boundary condition
func InitF1dBeam (Nn int, dNod, fNod []int, fVal []float64) ([]int, []float64) {
  dNum := len(dNod)
  fNum := len(fNod)
  nfNum := Nn-dNum
  nfNod := make([]int, nfNum)
  nfVal := make([]float64, nfNum)
  c, j, k := 0, 0, 0
  for i := 0; i < Nn; i++ {
    if j < dNum && dNod[j] == i {
      j++
    } else if c < nfNum {
      nfNod[c] = i
      if k < fNum && fNod[k] == i {
        nfVal[c] = fVal[k]
        k++
      }
      c++
    }
  }
  return nfNod, nfVal
}

// Return shape function of e-th element
func (fem *FEMsolver1dBeam) NElem(i, e int) func(float64) float64 {
  Le := fem.Le.At(e, 0)
  xe := Le * float64(e)
  return func(x float64) float64 {
    N := 0.0
    switch(i) {
    case 1:
      N = 1 - 3*math.Pow((x-xe),2)/math.Pow(Le,2) + 2*math.Pow((x-xe),3)/math.Pow(Le,3)
    case 2:
      N = (x-xe) - 2*math.Pow((x-xe),2)/Le + math.Pow((x-xe),3)/math.Pow(Le,2)
    case 3:
      N = 3*math.Pow((x-xe),2)/math.Pow(Le,2) - 2*math.Pow((x-xe),3)/math.Pow(Le,3)
    case 4:
      N = -math.Pow((x-xe),2)/Le + math.Pow((x-xe),3)/math.Pow(Le,2)
    }
    return N
  }
}

// Calculate local k matrix
func (fem *FEMsolver1dBeam) CalcLocK() {
  fmt.Println("Calc local stifness matrix ...")
  Le := fem.Le.At(0, 0)
  var ke float64
  ke = 12
  fem.k.SetSym(0, 0, ke)
  fem.k.SetSym(2, 2, ke)
  fem.k.SetSym(0, 2, -ke)
  ke = 4*math.Pow(Le, 2)
  fem.k.SetSym(1, 1, ke)
  fem.k.SetSym(3, 3, ke)
  fem.k.SetSym(1, 3, ke/2)
  ke = 6*Le
  fem.k.SetSym(0, 1, ke)
  fem.k.SetSym(0, 3, ke)
  fem.k.SetSym(1, 2, -ke)
  fem.k.SetSym(2, 3, -ke)
  if DEBUG {
    fmt.Printf("k = %0.4v\n", mat64.Formatted(fem.k, mat64.Prefix("    ")))
    fmt.Println()
  }
}

// Calculate the global K matrix
func (fem *FEMsolver1dBeam) CalcK() {
  fmt.Println("Calc stifness matrix ...")
  for k := 0; k < fem.Ne; k++ {
    var kloc mat64.SymDense
    EIL := fem.E.At(k,0)*fem.I.At(k,0)/math.Pow(fem.Le.At(k,0),3)
    kp := k*(fem.No-1)*2
    kloc.ScaleSym(EIL, fem.k)
    klocR, klocC := kloc.Dims()
    for i := 0; i < klocR; i++ {
      for j := 0; j < klocC; j++ {
        if j >= i {
          glo := fem.K.At(kp+i, kp+j)
          loc := kloc.At(i, j)
          fem.K.SetSym(kp+i, kp+j, glo+loc)
        }
      }
    }
  }
  if DEBUG {
    fmt.Printf("K = %0.4v\n", mat64.Formatted(fem.K, mat64.Prefix("    ")))
    fmt.Println()
  }
}

// Add body force (or distributed force) b to the solver
func (fem *FEMsolver1dBeam) AddBodyForce(b func(float64) float64) {
}
  fmt.Println("Add body force to force vector ...")
  for i := 0; i < fem.Ne; i++ {
    for j := 0; j < fem.No; j++ {
      Nj := fem.NElem(j, i)
      fj := GausQuad(func(x float64) float64 {return Nj(x)*b(x)},
              float64(i)*fem.Le.At(i, 0),
              float64(i+1)*fem.Le.At(i, 0), fem.Ng)
      for k, v := range fem.fNod {
        if v == (fem.No-1)*i+j {
          fem.fVal[k] += fj
        }
      }
    }
    BuildVec(fem.fNod, fem.fVal, fem.f)
  }
  if DEBUG {
    fmt.Println(fem.fNod)
    fmt.Println(fem.fVal)
    fmt.Printf("f = %0.4v\n", mat64.Formatted(fem.f, mat64.Prefix("    ")))
    fmt.Println()
  }

// Solve the problem
func (fem *FEMsolver1dBeam) Solve() {
  fmt.Println("Solving ...")

  if DEBUG {
    fmt.Printf("f = %0.4v\n", mat64.Formatted(fem.f, mat64.Prefix("    ")))
    fmt.Println()
  }
  fNum := len(fem.fNod)
  m := mat64.NewDense(fNum, fNum+1, nil)
  for i := 0; i < fNum; i++ {
    for j := 0; j < fNum; j++ {
      m.Set(i, j, fem.K.At(fem.fNod[i], fem.fNod[j]))
    }
    m.Set(i, fNum, fem.fVal[i])
  }
  if DEBUG {
    fmt.Printf("m = %0.4v\n", mat64.Formatted(m, mat64.Prefix("    ")))
    fmt.Println()
  }

  for i := 0; i < fNum; i++ {
    for j := i+1; j < fNum; j++ {
      c := m.At(j,i)/m.At(i,i)
      row := mat64.Row(nil, j, m)
      for k,v := range row {
        row[k] = v - c*m.At(i,k)
      }
      m.SetRow(j, row)
    }
  }
  if DEBUG {
    fmt.Printf("m = %0.4v\n", mat64.Formatted(m, mat64.Prefix("    ")))
    fmt.Println()
  }

  for i := fNum-1; i >= 0; i-- {
    sum := 0.0
    for j := i+1; j < fNum; j++ {
      sum += m.At(i,j)*fem.d.At(fem.fNod[j],0)
    }
    fem.d.SetVec(fem.fNod[i], (m.At(i, fNum) - sum)/m.At(i,i))
  }
  fmt.Printf("d = %0.4v\n", mat64.Formatted(fem.d, mat64.Prefix("    ")))
  fmt.Println()

  dNum := len(fem.dNod)
  for i := 0; i < dNum; i++ {
    fi := fem.dNod[i]
    for j := 0; j < fem.No*2; j++ {
      fv := fem.f.At(fi, 0)+fem.K.At(fi, j)*fem.d.At(j, 0)
      fem.f.SetVec(fi, fv)
    }
  }
  fmt.Printf("f = %0.4v\n", mat64.Formatted(fem.f, mat64.Prefix("    ")))
  fmt.Println()
}
