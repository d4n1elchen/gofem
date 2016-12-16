package femsolver

import (
  "fmt"
  "github.com/gonum/matrix/mat64"
)

const FORCED = 3e+33
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

// Initialize the force boundary condition
func InitF (Nn int, uNod, fNod []int, fVal []float64) ([]int, []float64) {
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

// Linear-element FEM solver for column
// k: local stifness matrix
// K: global stifness matrix
// Ne: amount of element
// Le: lenth of each element
// E: Young's modulus of each element
// A: area of each element
// u: displacement vector
// f: force vector
// uNod, uVal: displacement boundary condition
// fNod, fVal: force boundary condition
type FEMsolver1d struct {
  k, K *mat64.SymDense
  Ne, No, Ng int
  Le, E, A, u, f *mat64.Vector
  uNod, fNod []int
  uVal, fVal []float64
}

// Create a FEMsolver1d
func NewFEMsolver1d(No, Ne int, Le, E, A, u, f *mat64.Vector, uNod, fNod []int, uVal, fVal []float64) *FEMsolver1d {
  k := mat64.NewSymDense(No, nil)
  K := mat64.NewSymDense((No-1)*Ne+1, nil)
  BuildVec(uNod, uVal, u)
  fNod, fVal = InitF((No-1)*Ne+1, uNod, fNod, fVal)
  BuildVec(fNod, fVal, f)
  Ng := (No+1)/2 + 1
  return &FEMsolver1d{k, K, Ne, No, Ng, Le, E, A, u, f, uNod, fNod, uVal, fVal}
}

// Create a FEMsolver1d with const Le, E and A
func NewFEMsolver1dConstLeEA(No, Ne int, Le, E, A float64, u, f *mat64.Vector, uNod, fNod []int, uVal, fVal []float64) *FEMsolver1d {

  LeV := mat64.NewVector(Ne+1, nil)
  EV := mat64.NewVector(Ne+1, nil)
  AV := mat64.NewVector(Ne+1, nil)
  for i := 0; i < Ne+1; i++ {
    LeV.SetVec(i, Le)
    EV.SetVec(i, E)
    AV.SetVec(i, A)
  }

  return NewFEMsolver1d(No, Ne, LeV, EV, AV, u, f, uNod, fNod, uVal, fVal)
}

// Return shape function of e-th element
func (fem *FEMsolver1d) NElem(j, e int) func(float64) float64 {
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
func (fem *FEMsolver1d) BElem(j, e int) func(float64) float64 {
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
func (fem *FEMsolver1d) AddBodyForce(b func(float64) float64) {
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
func (fem *FEMsolver1d) CalcLocK() {
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
func (fem *FEMsolver1d) CalcK() {
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
func (fem *FEMsolver1d) Solve() {
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
