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
  Ne int
  Le, E, A, u, f *mat64.Vector
  uNod, fNod []int
  uVal, fVal []float64
}

// Create a FEMsolver1d
func NewFEMsolver1d(Ne int, Le, E, A, u, f *mat64.Vector, uNod, fNod []int, uVal, fVal []float64) *FEMsolver1d {
  k := mat64.NewSymDense(2, []float64{
    1, -1,
    -1, 1,
  })
  K := mat64.NewSymDense(Ne+1, nil)
  BuildVec(uNod, uVal, u)
  fNod, fVal = InitF(Ne+1, uNod, fNod, fVal)
  BuildVec(fNod, fVal, f)
  return &FEMsolver1d{k, K, Ne, Le, E, A, u, f, uNod, fNod, uVal, fVal}
}

// Create a FEMsolver1d with const Le, E and A
func NewFEMsolver1dConstLeEA(Ne int, Le, E, A float64, u, f *mat64.Vector, uNod, fNod []int, uVal, fVal []float64) *FEMsolver1d {

  LeV := mat64.NewVector(Ne+1, nil)
  EV := mat64.NewVector(Ne+1, nil)
  AV := mat64.NewVector(Ne+1, nil)
  for i := 0; i < Ne+1; i++ {
    LeV.SetVec(i, Le)
    EV.SetVec(i, E)
    AV.SetVec(i, A)
  }

  return NewFEMsolver1d(Ne, LeV, EV, AV, u, f, uNod, fNod, uVal, fVal)
}

// Return shape function of e-th element
func (fem *FEMsolver1d) NElem(i, e int) func(float64) float64 {
  Le := fem.Le.At(e, 0)
  x1 := Le * float64(e)
  x2 := x1 + Le
  return func(x float64) float64 {
    return []float64{(x2-x)/Le, (x-x1)/Le}[i]
  }
}

// Return derivative shape function of e-th element
func (fem *FEMsolver1d) BElem(i, e int) func(float64) float64 {
  Le := fem.Le.At(e, 0)
  return func(x float64) float64 {
    return []float64{1/Le, -1/Le}[i]
  }
}

// Add body force (or distributed force) b to the solver
func (fem *FEMsolver1d) AddBodyForce(b func(float64) float64) {
  for i := 0; i < fem.Ne; i++ {

    N0 := fem.NElem(0, i)
    N1 := fem.NElem(1, i)
    f0 := GausQuad(func(x float64) float64 {return N0(x)*b(x)},
            float64(i)*fem.Le.At(i, 0),
            float64(i+1)*fem.Le.At(i, 0), 3)
    f1 := GausQuad(func(x float64) float64 {return N1(x)*b(x)},
            float64(i)*fem.Le.At(i, 0),
            float64(i+1)*fem.Le.At(i, 0), 3)

    for j, v := range fem.fNod {
      switch v {
      case i:
        fem.fVal[j] += f0
      case i+1:
        fem.fVal[j] += f1
      }
    }
    BuildVec(fem.fNod, fem.fVal, fem.f)
  }
}

// Calculate the global K matrix
func (fem *FEMsolver1d) CalcK() {
  fmt.Println("Calc stifness matrix ...")
  for k := 0; k < fem.Ne; k++ {
    var kloc mat64.SymDense
    EAL := fem.E.At(k,0)*fem.A.At(k,0)/fem.Le.At(k,0)
    kloc.ScaleSym(EAL, fem.k)
    klocR, klocC := kloc.Dims()
    for i := 0; i < klocR; i++ {
      for j := 0; j < klocC; j++ {
        if j>=i{
          glo := fem.K.At(k+i, k+j)
          loc := kloc.At(i, j)
          fem.K.SetSym(k+i, k+j, glo+loc)
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
