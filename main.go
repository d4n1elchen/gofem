package main

import (
  "fmt"
  "github.com/gonum/matrix/mat64"
  "github.com/team6612/gofem/femsolver"
  "math"
)

func main() {
  femsolver.DEBUG = false

  var fem femsolver.FEMsolver
  Ne := 12
  E := 70e+9
  A := 0.03
  L := 3.0

  uNod := []int{0, Ne}
  uVal := []float64{0, 0}
  u := mat64.NewVector(Ne+1, nil)

  fNod := []int{8}
  fVal := []float64{5000.0}
  f := mat64.NewVector(Ne+1, nil)

  fem = femsolver.NewFEMsolver1dConstLeEA(Ne, L/float64(Ne), E, A, u, f, uNod, fNod, uVal, fVal)
  fem.CalcK()
  fem.Solve()

  gausXSin := femsolver.GausQuad(fx, -5, 5, 4)
  analXSin := ff(5) - ff(-5)
  fmt.Printf("gaussian: inte x*Sin(x) from -5 to 5 = %v\n", gausXSin)
  fmt.Printf("analysis: inte x*Sin(x) from -5 to 5 = %v\n", analXSin)

  fmt.Println("Main end")
}

func fx(x float64) float64 {
  return math.Pow(x,2) + x + 1
}

func ff(x float64) float64 {
  return math.Pow(x,3)/3 + math.Pow(x,2)/2 + x
}
