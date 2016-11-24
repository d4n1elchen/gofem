package main

import (
  "fmt"
  "github.com/gonum/matrix/mat64"
  "github.com/team6612/gofem/femsolver"
  "math"
)

func main() {
  femsolver.DEBUG = true

  var fem femsolver.FEMsolver
  Ne := 2
  E := 100e+9
  A := 0.0001
  L := 2.0

  uNod := []int{0}
  uVal := []float64{0}
  u := mat64.NewVector(Ne+1, nil)

  fNod := []int{}
  fVal := []float64{}
  f := mat64.NewVector(Ne+1, nil)

  fem = femsolver.NewFEMsolver1dConstLeEA(Ne, L/float64(Ne), E, A, u, f, uNod, fNod, uVal, fVal)
  fem.AddBodyForce(b)
  fem.CalcK()
  fem.Solve()

  gausXSin := femsolver.GausQuad(fx, -5, 5, 3)
  analXSin := ff(5) - ff(-5)
  fmt.Printf("gaussian: inte x*Sin(x) from -5 to 5 = %v\n", gausXSin)
  fmt.Printf("analysis: inte x*Sin(x) from -5 to 5 = %v\n", analXSin)

  fmt.Println("Main end")
}

func fx(x float64) float64 {
  return 0.8*x + 1.2345
}

func ff(x float64) float64 {
  return math.Pow(x,2)*0.4 + 1.2345*x
}

func b(x float64) float64 {
  return 1000
}
