package main

import (
  "fmt"
  "math"

  "bufio"
  "os"
  "strings"
  "regexp"
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
    fmt.Println(os.Args[1])
    Ne, _ = strconv.Atoi(os.Args[1])
  }
  Nn := 2
  E := 200e+9
  I := 5e-6
  L := 10.0
  Le := L/float64(Ne)

  dNod := []int{0, 1, Ne, 2*Ne}
  dVal := []float64{0, 0, 0, 0}
  // dNod := []int{0, 1}
  // dVal := []float64{0, 0}
  d := mat64.NewVector(2*(Nn-1)*Ne+2, nil)

  fNod := []int{}
  fVal := []float64{}
  f := mat64.NewVector(2*(Nn-1)*Ne+2, nil)

  fem = femsolver.NewFEMsolver1dBeamConstLeEI(Nn, Ne, Le, E, I, d, f, dNod, fNod, dVal, fVal)
  fem.AddBodyForce(q2, 4)
  fem.CalcLocK()
  fem.CalcK()
  fem.Solve()

  fmt.Println("Disp(1) = ", fem.Disp(1))
  fmt.Println("Disp(2) = ", fem.Disp(2))
  xitox := func (xi float64, e int) float64 {
    return 5.0*float64(e)+(xi+1)*Le/2
  }
  fmt.Println("Stress(gaus1) = ", fem.Stress(xitox(-1.0/math.Sqrt(3), 0)))
  fmt.Println("Stress(gaus2) = ", fem.Stress(xitox(+1.0/math.Sqrt(3), 0)))

  dc := int(L/0.01)
  dPts := make(plotter.XYs, dc)
  for i := 0; i < dc; i++ {
    dPts[i].X = float64(i)*0.01
    dPts[i].Y = -fem.Disp(float64(i)*0.01)
  }
  plotDisp(dPts)

  sc := int(L/0.01)
  sPts := make(plotter.XYs, sc)
  for i := 0; i < sc; i++ {
    sPts[i].X = float64(i)*0.01
    sPts[i].Y = -fem.Stress(float64(i)*0.01)
  }
  plotStress(sPts)

  fmt.Println("Main end")
}

func b(x float64) float64 {
  return 1000
}

func q(x float64) float64 {
  return 1000*(1-x/2)
}

func q1(x float64) float64 {
  if x < 5 {
    return 12.0
  } else {
    return 24.0
  }
}

func q2(x float64) float64 {
  if x < 5 {
    return 12.0+12.0/5.0*x
  } else {
    return 24.0
  }
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

  file, _ := os.Open("displacement_exact.txt")
  scanner := bufio.NewScanner(file)
  exact := make(plotter.XYs, 10001)
  i := 0
  for scanner.Scan() {
    lineStr := strings.TrimSpace(scanner.Text())
    s := regexp.MustCompile("\\s+").Split(lineStr, -1)
    exact[i].X, _ = strconv.ParseFloat(s[0], 64)
    exact[i].Y, _ = strconv.ParseFloat(s[1], 64)
    i += 1
  }
  err = plotutil.AddLines(p,
    "FEM", pts,
    "Exact", exact)
  if err != nil {
    panic(err)
  }

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

  file, _ := os.Open("stress_exact.txt")
  scanner := bufio.NewScanner(file)
  exact := make(plotter.XYs, 10001)
  i := 0
  for scanner.Scan() {
    lineStr := strings.TrimSpace(scanner.Text())
    s := regexp.MustCompile("\\s+").Split(lineStr, -1)
    exact[i].X, _ = strconv.ParseFloat(s[0], 64)
    exact[i].Y, _ = strconv.ParseFloat(s[1], 64)
    i += 1
  }
  err = plotutil.AddLines(p,
    "FEM", pts,
    "Exact", exact)
  if err != nil {
    panic(err)
  }

  // Save the plot to a PNG file.
  if err := p.Save(4*vg.Inch, 4*vg.Inch, "stress.png"); err != nil {
    panic(err)
  }
}
