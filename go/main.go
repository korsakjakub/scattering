package main

import (
	"fmt"
	"math"
	"os"
	"strconv"

	"go-hep.org/x/hep/fit"
	"gonum.org/v1/gonum/optimize"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
)

const De = 0.0011141
const Re = 10.98
const u = 1822.888
const m39k = 38.963703 * u
const m40k = 39.963999 * u

// const m41k = 40.961825 * u

var conf Config

func kSq(position, mu, collisionEnergy, angularMomentum float64) float64 {
	potential := func(position float64) float64 {
		return De * (math.Pow(Re/position, 12) - 2*math.Pow(Re/position, 6))
	}
	totalPotential := func(position float64) float64 {
		return potential(position) + angularMomentum*(angularMomentum+1)/(2*mu*math.Pow(position, 2))
	}
	return 2 * mu * (collisionEnergy - totalPotential(position))
}

func waveFunctionRatioGetN(prevRatio float64, index int, gridStep float64, kSqArray []float64) float64 {
	gridStepSq := math.Pow(gridStep, 2)
	return (2*(1-5*gridStepSq*kSqArray[index]/12) - (1+gridStepSq*kSqArray[index-1]/12)/prevRatio) / (1 + gridStepSq*kSqArray[index+1]/12)
}

func reducedMass(m1, m2 float64) float64 {
	return m1 * m2 / (m1 + m2)
}

func simulate(conf Config, collisionEnergy, angularMomentum float64) float64 {
	ratio, err := strconv.ParseFloat(conf.Ratio, 64)
	initialPosition, err := strconv.ParseFloat(conf.InitialPosition, 64)
	initialPhi, err := strconv.ParseFloat(conf.InitialPhi, 64)
	indexRangeFloat, err := strconv.ParseFloat(conf.IndexRange, 64)

	indexRange := int(indexRangeFloat)
	if err != nil {
		fmt.Println(err)
	}

	mu := reducedMass(m39k, m40k)
	gridStep := 2 * math.Pi / (50 * math.Sqrt(kSq(Re, mu, collisionEnergy, angularMomentum)))

	var k []float64
	for i := 0; i <= indexRange+1; i += 1 {
		k = append(k, kSq(initialPosition+float64(i)*gridStep, mu, collisionEnergy, angularMomentum))
	}

	var ratios []float64
	for i := 1; i <= indexRange; i += 1 {
		ratio = waveFunctionRatioGetN(ratio, i, gridStep, k)
		ratios = append(ratios, ratio)
	}

	var phi []float64
	var x []float64
	phi = append(phi, initialPhi)
	x = append(x, initialPosition)
	for i := 1; i < indexRange; i += 1 {
		phi = append(phi, phi[i-1]*ratios[i])
		x = append(x, initialPosition+float64(i)*gridStep)
	}

	/*
		Testing computing the K_l
	*/
	u := math.Sqrt(2.0 * mu * collisionEnergy)

	longRangeWFApproximation := func(position, scatteringPhase float64) float64 {
		return math.Sqrt((2 * mu / (math.Pi * u)) * (math.Sin(u*position-angularMomentum*math.Pi/2) + scatteringPhase*math.Cos(u*position-angularMomentum*math.Pi/2)))
	}

	res, err := fit.Curve1D(
		fit.Func1D{
			F: func(x float64, ps []float64) float64 {
				return longRangeWFApproximation(x, ps[0])
			},
			X: x,
			Y: phi,
			N: 1,
		},
		nil, &optimize.NelderMead{},
	)
	//if err := res.Status.Err(); err != nil {
	//fmt.Println(err)
	//}
	scatteringPhase := res.X[0]
	crossSection := 16 * math.Pi / math.Pow(u, 2) * math.Pow(scatteringPhase, 2) / (1 + math.Pow(scatteringPhase, 2)) * (2*angularMomentum + 1)
	return crossSection
}

func plotCrossSection(xys plotter.XYs, fileName string) {
	p := plot.New()
	s, err := plotter.NewScatter(xys)
	p.Add(s)
	wt, err := p.WriterTo(512, 512, "png")
	f, err := os.Create("go/out/" + fileName)
	_, err = wt.WriteTo(f)
	if err = f.Close(); err != nil {
		fmt.Printf("error: %v", err)
	}
}

func main() {
	conf = LoadConfig([]string{})
	collisionEnergy := 1e-10
	energies := make([]float64, 500)
	for i := 0; i < 500; i += 1 {
		energies[i] = collisionEnergy + float64(6*i)*collisionEnergy
	}
	angularMomentum := 0.0
	var xyss []plotter.XYs
	var x, y float64
	totalCrossSections := make([]float64, 500)
	for l := 0; l <= 10; l += 1 {
		angularMomentum = float64(l)
		var xys plotter.XYs
		for i := 0; i < 500; i += 1 {
			x = energies[i]
			y = simulate(conf, x, angularMomentum)
			totalCrossSections[i] += y
			xys = append(xys, struct{ X, Y float64 }{x, y})
		}
		xyss = append(xyss, xys)
		// plotCrossSection(xyss[len(xyss)-1], fmt.Sprintf("sigma%v.png", angularMomentum))
	}
	var totalXYs plotter.XYs
	for i := 0; i < 500; i += 1 {
		totalXYs = append(totalXYs, struct{ X, Y float64 }{energies[i], totalCrossSections[i]})
	}
	plotCrossSection(totalXYs, "totalSigma.png")

	/*
		x, phi := simulate(conf, collisionEnergy, angularMomentum)
		graph := chart.Chart{
			Width:  1000,
			Height: 1000,
			XAxis: chart.XAxis{
				Name: "R/bohr",
			},
			YAxis: chart.YAxis{
				Name: "wf",
			},
			Series: []chart.Series{
				chart.ContinuousSeries{
					Style: chart.Style{
						StrokeColor: chart.GetDefaultColor(0).WithAlpha(255),
					},
					XValues: x,
					YValues: phi,
				},
			},
		}
		f, _ := os.Create("go.png")
		defer func(f *os.File) {
			err := f.Close()
			if err != nil {
				fmt.Println(err)
			}
		}(f)
		err := graph.Render(chart.PNG, f)
		if err != nil {
			fmt.Println(err)
		}
	*/
}
