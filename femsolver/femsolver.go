package femsolver

import (
  "fmt"
  "github.com/gonum/matrix/mat64"
)

const FORCED = 3e+33
var DEBUG = true

type FEMsolver interface {
  CalcN()
  CalcK()
  Solve()
}

type FEMsolver1d struct {
  k, K *mat64.SymDense
  Ne int
  Le, E, A, u, f *mat64.Vector
  uNod, fNod []int
  uVal, fVal []float64
}

func BuildVec(nod []int, val []float64, vec *mat64.Vector) {
  for i,n := range nod {
    vec.SetVec(n, val[i])
  }
}

// TODO: Implement GausQuad
func GausQuad(f func(float64) float64, a, b float64, num int) float64 {
  pt, w := getGausQuadPoint(num)
  f = changeInterval(f, a, b)
  sum := 0.0
  for i := 0; i < num; i++ {
    sum += w[i]*f(pt[i])
  }
  return sum
}

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

func changeInterval(f func (float64) float64, a, b float64) func (float64) float64 {
  return func (x float64) float64 {
    return (b-a)/2 * f( (b-a)*x/2 + (a+b)/2  )
  }
}

func InitF (Ne int, uNod, fNod []int, fVal []float64) ([]int, []float64) {
  uNum := len(uNod)
  fNum := len(fNod)
  nfNum := Ne+1-uNum
  nfNod := make([]int, nfNum)
  nfVal := make([]float64, nfNum)
  c, j, k := 0, 0, 0
  for i := 0; i < Ne; i++ {
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

func NewFEMsolver1d(Ne int, Le, E, A, u, f *mat64.Vector, uNod, fNod []int, uVal, fVal []float64) *FEMsolver1d {
  k := mat64.NewSymDense(2, []float64{
    1, -1,
    -1, 1,
  })
  K := mat64.NewSymDense(Ne+1, nil)
  BuildVec(uNod, uVal, u)
  fNod, fVal = InitF(Ne, uNod, fNod, fVal)
  BuildVec(fNod, fVal, f)
  return &FEMsolver1d{k, K, Ne, Le, E, A, u, f, uNod, fNod, uVal, fVal}
}

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

// TODO: Implement CalcN()
func (fem *FEMsolver1d) CalcN() {
}

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
    for j := 0; j < fem.Ne; j++ {
      fem.f.SetVec(fi, fem.f.At(fi, 0)+fem.K.At(fi, j)*fem.u.At(j, 0))
    }
  }
  fmt.Printf("f = %0.4v\n", mat64.Formatted(fem.f, mat64.Prefix("    ")))
  fmt.Println()
}
