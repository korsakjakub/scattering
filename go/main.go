package main

import (
	"fmt"
	"image/color"
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
const m41k = 40.961825 * u
const hartreeEnergyInKelvins = 3.1668e-6

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

func simulate(conf Config, collisionEnergy, angularMomentum float64, plotWF ...bool) float64 {
	if plotWF == nil {
		plotWF[0] = false
	}
	ratio, err := strconv.ParseFloat(conf.Ratio, 64)
	if err != nil {
		fmt.Printf("error: %v", err)
	}
	initialPosition, err := strconv.ParseFloat(conf.InitialPosition, 64)
	if err != nil {
		fmt.Printf("error: %v", err)
	}
	initialPhi, err := strconv.ParseFloat(conf.InitialPhi, 64)
	if err != nil {
		fmt.Printf("error: %v", err)
	}
	indexRangeFloat, err := strconv.ParseFloat(conf.IndexRange, 64)
	if err != nil {
		fmt.Printf("error: %v", err)
	}

	indexRange := int(indexRangeFloat)
	if err != nil {
		fmt.Println(err)
	}

	mu := reducedMass(m39k, m39k)
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

	longRangeWFApproximation := func(position, scatteringPhase, a float64) float64 {
		return a * (math.Sin(u*position-angularMomentum*math.Pi/2) + scatteringPhase*math.Cos(u*position-angularMomentum*math.Pi/2))
		// return math.Sqrt(2*mu/(math.Pi*u)) * (math.Sin(u*position-angularMomentum*math.Pi/2) + scatteringPhase*math.Cos(u*position-angularMomentum*math.Pi/2))
	}

	res, err := fit.Curve1D(
		fit.Func1D{
			F: func(x float64, ps []float64) float64 {
				return longRangeWFApproximation(x, ps[0], ps[1])
			},
			X: x[200:],
			Y: phi[200:],
			N: 2,
		},
		nil, &optimize.NelderMead{},
	)
	//if err := res.Status.Err(); err != nil {
	//fmt.Println(err)
	//}
	scatteringPhase := res.X[0]
	a := res.X[1]
	if math.IsNaN(scatteringPhase) {
		scatteringPhase = 0.0
	}
	_ = a
	if plotWF[0] == true {
		var modelXYs plotter.XYs
		var fitXYs plotter.XYs
		var diffXYs plotter.XYs
		for i := 1; i < indexRange; i += 1 {
			modelXYs = append(modelXYs, struct{ X, Y float64 }{x[i], phi[i]})
			fitXYs = append(fitXYs, struct{ X, Y float64 }{x[i], longRangeWFApproximation(x[i], scatteringPhase, a)})
			diffXYs = append(diffXYs, struct{ X, Y float64 }{x[i-1], modelXYs[i-1].Y - fitXYs[i-1].Y})
		}
		plotWaveFunctions([]plotter.XYs{modelXYs, fitXYs, diffXYs}, fmt.Sprintf("model-vs-fit%v.png", collisionEnergy), "out/en/")
	}

	crossSection := 16 * math.Pi / math.Pow(u, 2) * math.Pow(scatteringPhase, 2) / (1 + math.Pow(scatteringPhase, 2)) * (2*angularMomentum + 1)
	return crossSection
}

func plotWaveFunctions(xyss []plotter.XYs, fileName, outDir string) {
	p := plot.New()
	for n, xys := range xyss {
		s, err := plotter.NewScatter(xys)
		if err != nil {
			fmt.Printf("error: %v", err)
		}
		s.GlyphStyle.Color = color.RGBA{R: uint8(50 + 200*n), G: uint8(50), B: uint8(50), A: 255}
		p.Add(s)
	}
	wt, err := p.WriterTo(512, 512, "png")
	if err != nil {
		fmt.Printf("error: %v", err)
	}
	f, err := os.Create(outDir + fileName)
	if err != nil {
		fmt.Printf("error: %v", err)
	}
	_, err = wt.WriteTo(f)
	if err != nil {
		fmt.Printf("error: %v", err)
	}
	if err = f.Close(); err != nil {
		fmt.Printf("error: %v", err)
	}
}

func plotCrossSection(xys plotter.XYs, fileName, outDir string) {
	p := plot.New()
	p.X.Label.Text = "T/K"
	p.Y.Label.Text = "cross section"
	p.X.Scale = plot.LogScale{}
	p.Y.Scale = plot.LogScale{}
	p.X.Tick.Marker = plot.LogTicks{}
	p.Y.Tick.Marker = plot.LogTicks{}
	s, err := plotter.NewScatter(xys)
	if err != nil {
		fmt.Printf("error: %v", err)
	}
	p.Add(s)
	wt, err := p.WriterTo(512, 512, "png")
	if err != nil {
		fmt.Printf("error: %v", err)
	}
	f, err := os.Create(outDir + fileName)
	if err != nil {
		fmt.Printf("error: %v", err)
	}
	_, err = wt.WriteTo(f)
	if err != nil {
		fmt.Printf("error: %v", err)
	}
	if err = f.Close(); err != nil {
		fmt.Printf("error: %v", err)
	}
}

func main() {
	conf = LoadConfig([]string{})

    // Prepare energy parameter array
	initColEnergy := hartreeEnergyInKelvins * 1e-5
    finColEnergy := hartreeEnergyInKelvins * 10
    energyCount, err := strconv.Atoi(conf.EnergyCount)
    if err != nil {
        fmt.Println(err)
    }
	energies := make([]float64, energyCount)
	for i := 0; i < energyCount; i += 1 {
		// energies[i] = initColEnergy + (finColEnergy - initColEnergy) * float64(i)/float64(energyCount-1)
		t := math.Pow(finColEnergy/initColEnergy, 1/float64(energyCount))
		energies[i] = initColEnergy * math.Pow(t, float64(i))
	}

    angularMomentumCount, err := strconv.Atoi(conf.AngularMomentumCount)
    if err != nil {
        fmt.Println(err)
    }
	angularMomenta := make([]float64, angularMomentumCount)
	for i := 0; i < angularMomentumCount; i += 1 {
		angularMomenta[i] = float64(i)
	}

    // Total Cross section
    totalCrossSections := make([]float64, energyCount)
    for e := range energies {
        totalCrossSections[e] = 0.0
    }

    // Calculate scatter points where x - energy, y - cross section, for some angularMomentum value
    for _, angularMomentum := range angularMomenta {
		xys := make(plotter.XYs, energyCount)
        for i, energy := range energies {
            crossSection := simulate(conf, energy, angularMomentum, true)
            xys[i] = struct{ X, Y float64 }{energy/hartreeEnergyInKelvins, crossSection}
            totalCrossSections[i] += crossSection
		}
        // for each angularMomentum plot scatter cross section
		plotCrossSection(xys, fmt.Sprintf("sigma%v.png", angularMomentum), conf.OutputDir)
	}

    // plot total cross section
    totalXYs := make(plotter.XYs, energyCount)
    for i, energy := range energies {
        totalXYs[i] = struct{ X, Y float64 }{energy, totalCrossSections[i]/hartreeEnergyInKelvins}
	}
	plotCrossSection(totalXYs, "totalSigma.png", conf.OutputDir)
}
