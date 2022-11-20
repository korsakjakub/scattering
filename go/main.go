package main

import (
	"fmt"
	"github.com/wcharczuk/go-chart/v2"
	"math"
	"os"
	"strconv"
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

func simulate(conf Config, collisionEnergy, angularMomentum float64) []float64 {
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
	phi = append(phi, initialPhi)
	for i := 0; i < indexRange; i += 1 {
		phi = append(phi, phi[i]*ratios[i])
	}

	return phi
}

func main() {
	conf = LoadConfig([]string{})
	collisionEnergy := 1e-8
	angularMomentum := 0.0
	phi := simulate(conf, collisionEnergy, angularMomentum)
	indexRangeFloat, err := strconv.ParseFloat(conf.IndexRange, 64)
	indexRange := int(indexRangeFloat)

	low, high := 1, indexRange
	x := make([]float64, high-low+1)
	for i := range x {
		x[i] = float64(i + low)
	}
	graph := chart.Chart{
		Width:  1000,
		Height: 1000,
		XAxis: chart.XAxis{
			Name: "R",
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
	err = graph.Render(chart.PNG, f)
	if err != nil {
		fmt.Println(err)
	}
}
