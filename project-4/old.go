package main

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

func inite(d, h, f mat.Matrix) (eelec float64) {
	r, c := d.Dims()
	for mu := 0; mu < r; mu++ {
		for nu := 0; nu < c; nu++ {
			eelec += d.At(mu, nu) * (h.At(mu, nu) + f.At(mu, nu))
		}
	}
	return
}

func fockmat(h, d mat.Matrix, t []float64) *mat.SymDense {
	r, _ := h.Dims()
	fock := mat.NewSymDense(r, nil)
	for m := 0; m < r; m++ {
		for n := 0; n < r; n++ {
			var sum float64
			for l := 0; l < r; l++ {
				for s := 0; s < r; s++ {
					sum += d.At(l, s) *
						(2*t[compind(m, n, l, s)] -
							t[compind(m, l, n, s)])
				}
			}
			fock.SetSym(m, n, h.At(m, n)+sum)
		}
	}
	return fock
}

func printmat(matrix mat.Matrix) {
	r, c := matrix.Dims()
	for i := 0; i < r; i++ {
		fmt.Printf("%5d", i+1)
		for j := 0; j < c; j++ {
			fmt.Printf("%12.8f", matrix.At(i, j))
		}
		fmt.Print("\n")
	}
}

func madd(a, b *mat.SymDense) *mat.SymDense {
	var sum mat.SymDense
	sum.AddSym(a, b)
	return &sum
}
