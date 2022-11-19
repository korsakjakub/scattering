package main

import (
	"fmt"
	"github.com/wcharczuk/go-chart"
	"math"
	"os"
)

const De = 0.0011141
const Re = 10.98
const angularMomentum = 3.0

func kSq(position, mu, collisionEnergy float64) float64 {
	potential := func(position float64) float64 {
		return De * (math.Pow(Re/position, 12) - 2*math.Pow(Re/position, 6))
	}
	totalPotential := func(position float64) float64 {
		return potential(position) + angularMomentum*(angularMomentum+1)/(2*mu*math.Pow(position, 2))
	}
	v := totalPotential(position)
	return 2 * mu * (collisionEnergy - v)
}

func waveFunctionRatioGetN(prevRatio float64, index int, gridStep float64, kSqArray []float64) float64 {
	gridStepSq := math.Pow(gridStep, 2)
	return (2*(1-5*gridStepSq*kSqArray[index]/12) - (1+gridStepSq*kSqArray[index-1]/12)/prevRatio) / (1 + gridStepSq*kSqArray[index+1]/12)
}

func reducedMass(m1, m2 float64) float64 {
	return m1 * m2 / (m1 + m2)
}

func main() {
	ratio := 10.0
	u := 1822.888
	m39k := 38.963703 * u
	m40k := 39.963999 * u
	// m41k := 40.961825 * 1822.888
	initialPosition := 10.0
	initialPhi := 1.0
	collisionEnergy := 1.0
	indexRange := 200

	mu := reducedMass(m39k, m40k)
	gridStep := 2 * math.Pi / (50 * math.Sqrt(kSq(Re, mu, collisionEnergy)))

	var k []float64
	for i := 0; i <= indexRange+1; i += 1 {
		k = append(k, kSq(initialPosition, mu, collisionEnergy))
	}

	var ratios []float64
	for i := 1; i <= indexRange; i += 1 {
		ratio = waveFunctionRatioGetN(ratio, i, gridStep, k)
		ratios = append(ratios, ratio)
		fmt.Println(ratio)
	}

	var phi []float64
	phi = append(phi, initialPhi)
	for i := 0; i < indexRange; i += 1 {
		phi = append(phi, phi[i]*ratios[i])
	}

	low, high := 1, indexRange
	x := make([]float64, high-low+1)
	for i := range x {
		x[i] = float64(i + low)
	}
	graph := chart.Chart{
		Series: []chart.Series{
			chart.ContinuousSeries{
				XValues: x,
				YValues: phi,
			},
		},
	}
	f, _ := os.Create("output.png")
	defer f.Close()
	err := graph.Render(chart.PNG, f)
	if err != nil {
		fmt.Println(err)
	}
}
