package main

import (
  "fmt"
  "image/color"

  "os"
  "strconv"

  "github.com/team6612/gofem/femsolver"

  "github.com/gonum/matrix/mat64"

  "github.com/gonum/plot"
  "github.com/gonum/plot/plotter"
  "github.com/gonum/plot/plotutil"
  "github.com/gonum/plot/vg"
)

func check(e error) {
  if e != nil {
    panic(e)
  }
}

func main() {
  femsolver.DEBUG = false

  var fem femsolver.FEMsolver
  Ne := 10
  if len(os.Args) > 1 {
    Ne, _ = strconv.Atoi(os.Args[1])
  }
  Nn := 2
  if len(os.Args) > 2 {
    Nn, _ = strconv.Atoi(os.Args[2])
  }
  E := 100e+9
  A := 0.0001
  L := 2.0
  Le := L/float64(Ne)

  uNod := []int{0}
  uVal := []float64{0}
  u := mat64.NewVector((Nn-1)*Ne+1, nil)

  fNod := []int{}
  fVal := []float64{}
  f := mat64.NewVector((Nn-1)*Ne+1, nil)

  fem = femsolver.NewFEMsolver1dBarConstLeEA(Nn, Ne, Le, E, A, u, f, uNod, fNod, uVal, fVal)
  fem.AddBodyForce(b, 3)
  fem.CalcLocK()
  fem.CalcK()
  fem.Solve()

  fmt.Println("Disp(0.5)=", fem.Disp(0.5))
  fmt.Println("Disp(1.5)=", fem.Disp(1.5))
  fmt.Println("Stress(0.5)=", fem.Stress(0.5))
  fmt.Println("Stress(1.5)=", fem.Stress(1.5))

  dc := int(L/0.01)
  dPts := make(plotter.XYs, dc)
  for i := 0; i < dc; i++ {
    dPts[i].X = float64(i)*0.01
    dPts[i].Y = fem.Disp(float64(i)*0.01)
  }
  plotDisp(dPts)

  sPts := make(plotter.XYs, dc)
  for i := 0; i < dc; i++ {
    sPts[i].X = float64(i)*0.01
    sPts[i].Y = fem.Stress(float64(i)*0.01)
  }
  plotStress(sPts)

  fmt.Println("Main end")
}

func b(x float64) float64 {
  return 1000
}

func plotDisp(pts plotter.XYs) {
  // Plot for diaplacement
  p, err := plot.New()
  if err != nil {
    panic(err)
  }

  p.Title.Text = "Position plot"
  p.X.Label.Text = "X"
  p.Y.Label.Text = "Y"

  err = plotutil.AddLines(p,
    "FEM", pts)
  if err != nil {
    panic(err)
  }

  exact := plotter.NewFunction(func(x float64) float64 { return 0.0002*x-0.00005*x*x })
  exact.Color = color.RGBA{B: 255, A: 255}

  p.Add(exact)
  p.Legend.Add("Exact", exact)

  // Save the plot to a PNG file.
  if err := p.Save(4*vg.Inch, 4*vg.Inch, "displacement.png"); err != nil {
    panic(err)
  }
}

func plotStress(pts plotter.XYs) {
  // Plot for stress
  p, err := plot.New()
  if err != nil {
    panic(err)
  }

  p.Title.Text = "Stress plot"
  p.X.Label.Text = "X"
  p.Y.Label.Text = "Y"

  err = plotutil.AddLines(p,
    "FEM", pts)
  if err != nil {
    panic(err)
  }

  exact := plotter.NewFunction(func(x float64) float64 { return 1000/0.0001*(2-x) })
  exact.Color = color.RGBA{B: 255, A: 255}

  p.Add(exact)
  p.Legend.Add("Exact", exact)

  p.X.Min = 0
  p.X.Max = 2
  p.Y.Min = 0
  p.Y.Max = 2e7

  // Save the plot to a PNG file.
  if err := p.Save(4*vg.Inch, 4*vg.Inch, "stress.png"); err != nil {
    panic(err)
  }
}
