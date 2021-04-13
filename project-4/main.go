/*
   This is a script illustrating the intricacies of a quantum chemical suite's ability to calcualte the SCF
   converged energy specifically for H2O using STO-3G as a basis set and CCSD as the method. Previous lines
   of code in this script were pulled from a very similar code that computed the same energy at the same
   basis set level, but with HF as the method.
*/
package main

import (
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"strconv"
	"strings"

	"math"

	"gonum.org/v1/gonum/mat"
)

// This is for reading in the nuclear repulsion energy from a separate file and converting it to a float
func enuc(filename string) float64 {
	stuff, err := ioutil.ReadFile(filename)
	if err != nil {
		log.Fatal(err)
	}
	stuff2 := strings.TrimSpace(string(stuff))
	stuff3, _ := strconv.ParseFloat(stuff2, 64)
	// fmt.Println(stuff3)
	return stuff3
}

/*
   This takes an input file and splits the contents into a format that can be easily stored
   as a matrix using gonum/mat package.
*/
func maketrix(filename string) *mat.SymDense {
	input, _ := os.Open(filename)
	defer input.Close()
	scan := bufio.NewScanner(input)
	matrix := mat.NewSymDense(1, nil)
	for scan.Scan() {
		line := strings.Fields(scan.Text())
		i, _ := strconv.Atoi(line[0])
		j, _ := strconv.Atoi(line[1])
		val, _ := strconv.ParseFloat(line[2], 64)
		// fmt.Printf("%T %T %f\n", i, j, val)
		if i < j {
			i, j = j, i
		}
		n, _ := matrix.Dims()
		// fmt.Println(i, j, n)
		if n < i {
			matrix = (matrix.GrowSym(i - n)).(*mat.SymDense)
		}
		matrix.SetSym(j-1, i-1, val)
	}
	return matrix
}

/*
   This is a simple script that takes integers and stores them as a compound index that is used
   for pulling specific values from the two electron integrals
*/
func compind(m, n, l, s int) int {
	if m < n {
		m, n = n, m
	}
	mn := m*(m+1)/2 + n
	if l < s {
		l, s = s, l
	}
	ls := l*(l+1)/2 + s
	if mn < ls {
		mn, ls = ls, mn
	}
	mnls := mn*(mn+1)/2 + ls
	return mnls
}

/* This takes the two electron input file and, using the compind function, returns it as an indexable
   slice of floats.
*/
func twoelec(filename string) []float64 {
	input, _ := os.Open(filename)
	defer input.Close()
	scan := bufio.NewScanner(input)
	slice := make([]float64, 0)
	for scan.Scan() {
		line := strings.Fields(scan.Text())
		m, _ := strconv.Atoi(line[0])
		n, _ := strconv.Atoi(line[1])
		l, _ := strconv.Atoi(line[2])
		s, _ := strconv.Atoi(line[3])
		val, _ := strconv.ParseFloat(line[4], 64)
		// fmt.Printf("%T %T %T %T %f\n", m, n, l, s, val)
		mnls := compind(m-1, n-1, l-1, s-1)
		for mnls >= len(slice) {
			slice = append(slice, 0)
		}
		slice[mnls] = val
	}
	// fmt.Println(len(slice))
	return slice
}

/*
   This is a function that is used to make transform the AO-basis two electron integrals into the MO-basis
*/
func noddy(c *mat.Dense, t []float64) []float64 {
	row, _ := c.Dims()
	slice := make([]float64, 0)
	for p := 0; p < row; p++ {
		for q := 0; q < row; q++ {
			for r := 0; r < row; r++ {
				for s := 0; s < row; s++ {
					var sum float64
					for m := 0; m < row; m++ {
						for n := 0; n < row; n++ {
							for l := 0; l < row; l++ {
								for o := 0; o < row; o++ {
									sum += c.At(m, p) *
										c.At(n, q) *
										t[compind(m, n, l, o)] *
										c.At(l, r) *
										c.At(o, s)
								}
							}
						}
					}
					pqrs := compind(p, q, r, s)
					for pqrs >= len(slice) {
						slice = append(slice, 0)
					}
					slice[pqrs] = sum
				}
			}
		}
	}
	return slice
}

// A function that computes energy using MP2 as its method
func mp2energy(mo, noddy []float64, nocc int) float64 {
	var sum float64
	for i := 0; i < nocc; i++ {
		for j := 0; j < nocc; j++ {
			for a := nocc; a < len(mo); a++ {
				for b := nocc; b < len(mo); b++ {
					sum += noddy[compind(i, a, j, b)] *
						(2*noddy[compind(i, a, j, b)] -
							noddy[compind(i, b, j, a)]) /
						(mo[i] + mo[j] - mo[a] - mo[b])
				}
			}
		}
	}
	return sum
}

// A function that computes two electron integrals and is easier to store for the code.
func smart(c *mat.Dense, t []float64) (smarti []float64) {
	row, _ := c.Dims()
	map1 := make(map[[4]int]float64)
	map2 := make(map[[4]int]float64)
	map3 := make(map[[4]int]float64)
	for m := 0; m < row; m++ {
		for n := 0; n < row; n++ {
			for l := 0; l < row; l++ {
				for s := 0; s < row; s++ {
					var sumo float64
					for o := 0; o < row; o++ {
						sumo += c.At(o, s) * t[compind(m, n, l, o)]
					}
					map1[[4]int{m, n, l, s}] = sumo
				}
			}
		}
	}
	for m := 0; m < row; m++ {
		for n := 0; n < row; n++ {
			for r := 0; r < row; r++ {
				for s := 0; s < row; s++ {
					var suml float64
					for l := 0; l < row; l++ {
						suml += c.At(l, r) * map1[[4]int{m, n, l, s}]
					}
					map2[[4]int{m, n, r, s}] = suml
				}
			}
		}
	}
	for m := 0; m < row; m++ {
		for q := 0; q < row; q++ {
			for r := 0; r < row; r++ {
				for s := 0; s < row; s++ {
					var sumq float64
					for n := 0; n < row; n++ {
						sumq += c.At(n, q) * map2[[4]int{m, n, r, s}]
					}
					map3[[4]int{m, q, r, s}] = sumq
				}
			}
		}
	}
	for p := 0; p < row; p++ {
		for q := 0; q < row; q++ {
			for r := 0; r < row; r++ {
				for s := 0; s < row; s++ {
					var summ float64
					for m := 0; m < row; m++ {
						summ += c.At(m, p) * map3[[4]int{m, q, r, s}]
					}
					pqrs := compind(p, q, r, s)
					for pqrs >= len(smarti) {
						smarti = append(smarti, 0)
					}
					smarti[pqrs] = summ
				}
			}
		}
	}
	return
}

// This turns the spatial-basis two electron integrals into the spin-basis two electron integrals.
func spattospin(smarti []float64, nmo int) (spin [][][][]float64) {
	spin = make([][][][]float64, nmo)
	var v1 float64
	var v2 float64
	for p := 0; p < nmo; p++ {
		spin[p] = make([][][]float64, nmo)
		for q := 0; q < nmo; q++ {
			spin[p][q] = make([][]float64, nmo)
			for r := 0; r < nmo; r++ {
				spin[p][q][r] = make([]float64, nmo)
				for s := 0; s < nmo; s++ {
					if p%2 == r%2 && q%2 == s%2 {
						v1 = smarti[compind(p/2, r/2, q/2, s/2)]
					} else {
						v1 = 0
					}
					if p%2 == s%2 && q%2 == r%2 {
						v2 = smarti[compind(p/2, s/2, q/2, r/2)]
					} else {
						v2 = 0
					}
					spin[p][q][r][s] = v1 - v2
				}
			}
		}
	}
	return
}

// Transforms the Fock matrix into a spin basis version.
func spinfock(h *mat.SymDense, spin [][][][]float64, nmo int, nsocc int) (sofock *mat.Dense) {
	r, _ := h.Dims()
	sofock = mat.NewDense(2*r, 2*r, nil)
	for p := 0; p < nmo; p++ {
		for q := 0; q < nmo; q++ {
			var sum float64
			for m := 0; m < nsocc; m++ {
				sum += spin[p][m][q][m]
			}
			if p%2 == q%2 {
				sofock.Set(p, q, h.At(p/2, q/2)+sum)
			} else {
				sofock.Set(p, q, 0)
			}
		}
	}
	return
}

// Computes the enegrgy using coupled cluster amplitudes
func ccmp2(spinf *mat.Dense, spin [][][][]float64, nmo int, nsocc int, tijab [][][][]float64) float64 {
	// var tijab float64
	var sum float64
	var emp2 float64
	for i := 0; i < nsocc; i++ {
		for j := 0; j < nsocc; j++ {
			for a := nsocc; a < nmo; a++ {
				for b := nsocc; b < nmo; b++ {
					val := tijab[i][j][a][b]
					// tijab = spin[i][j][a][b] / (spinf.At(i, i) +
					// 	spinf.At(j, j) -
					// 	spinf.At(a, a) -
					// 	spinf.At(b, b))
					sum += spin[i][j][a][b] * val
				}
			}
		}
	}
	emp2 = (0.25) * sum
	return emp2
}

// Computes the single excitation contributions for CCSD
func tia(nsocc, nmo int) *mat.Dense {
	matrix := mat.NewDense(nmo, nmo, nil)
	for i := 0; i < nsocc; i++ {
		for a := nsocc; a < nmo; a++ {
			matrix.Set(i, a, 0.0)
		}
	}
	return matrix
}

// Computes the double excitation contributions for CCSD
func tijab(spin [][][][]float64, spinf *mat.Dense, nsocc, nmo int) [][][][]float64 {
	val := make([][][][]float64, nmo)
	for i := 0; i < nsocc; i++ {
		val[i] = make([][][]float64, nmo)
		for j := 0; j < nsocc; j++ {
			val[i][j] = make([][]float64, nmo)
			for a := nsocc; a < nmo; a++ {
				val[i][j][a] = make([]float64, nmo)
				for b := nsocc; b < nmo; b++ {
					val[i][j][a][b] = spin[i][j][a][b] / (spinf.At(i, i) +
						spinf.At(j, j) -
						spinf.At(a, a) -
						spinf.At(b, b))
				}
			}
		}
	}
	return val
}

// A simple function to represent the kronecker delta.
func kronk(i, j int) float64 {
	var val float64
	if i == j {
		val = 1.0
	} else {
		val = 0.0
	}
	return val
}

// Respresents the tilde tau function in the CCSD intermediate equations.
func ttau(spin [][][][]float64, tijab [][][][]float64, nmo, nsocc int, tia *mat.Dense) [][][][]float64 {
	val := make([][][][]float64, nmo)
	for i := 0; i < nsocc; i++ {
		val[i] = make([][][]float64, nmo)
		for j := 0; j < nsocc; j++ {
			val[i][j] = make([][]float64, nmo)
			for a := nsocc; a < nmo; a++ {
				val[i][j][a] = make([]float64, nmo)
				for b := nsocc; b < nmo; b++ {
					val[i][j][a][b] = tijab[i][j][a][b] +
						(0.5 * (tia.At(i, a)*tia.At(j, b) - tia.At(i, b)*tia.At(j, a)))
				}
			}
		}
	}
	return val
}

// Represens the tau function in the CCSD intermediate equations.
func tau(spin [][][][]float64, tijab [][][][]float64, nmo, nsocc int, tia *mat.Dense) [][][][]float64 {
	val := make([][][][]float64, nmo)
	for i := 0; i < nsocc; i++ {
		val[i] = make([][][]float64, nmo)
		for j := 0; j < nsocc; j++ {
			val[i][j] = make([][]float64, nmo)
			for a := nsocc; a < nmo; a++ {
				val[i][j][a] = make([]float64, nmo)
				for b := nsocc; b < nmo; b++ {
					val[i][j][a][b] = tijab[i][j][a][b] +
						(tia.At(i, a) * tia.At(j, b)) - (tia.At(i, b) * tia.At(j, a))
				}
			}
		}
	}
	return val
}

// F_ae intermediate equation for CCSD.
func faeint(spinf, tia *mat.Dense, spin [][][][]float64, ttau [][][][]float64, nsocc, nmo int) *mat.Dense {
	var val float64
	matrix := mat.NewDense(nmo, nmo, nil)
	for a := nsocc; a < nmo; a++ {
		for e := nsocc; e < nmo; e++ {
			val = (1 - kronk(a, e)) * spinf.At(a, e)
			for m := 0; m < nsocc; m++ {
				val -= 0.5 * spinf.At(m, e) * tia.At(m, a)
				for f := nsocc; f < nmo; f++ {
					val += tia.At(m, f) * spin[m][a][f][e]
					for n := 0; n < nsocc; n++ {
						val -= 0.5 * ttau[m][n][a][f] * spin[m][n][e][f]
					}
				}
			}
			matrix.Set(a, e, val)
		}
	}
	return matrix
}

// F_mi intermediate equation for CCSD.
func fmiint(spinf, tia *mat.Dense, spin [][][][]float64, ttau [][][][]float64, nsocc, nmo int) *mat.Dense {
	var val float64
	matrix := mat.NewDense(nmo, nmo, nil)
	for m := 0; m < nsocc; m++ {
		for i := 0; i < nsocc; i++ {
			val = (1 - kronk(m, i)) * spinf.At(m, i)
			for e := nsocc; e < nmo; e++ {
				val += 0.5 * tia.At(i, e) * spinf.At(m, e)
				for n := 0; n < nsocc; n++ {
					val += tia.At(n, e) * spin[m][n][i][e]
					for f := nsocc; f < nmo; f++ {
						val += 0.5 * ttau[i][n][e][f] * spin[m][n][e][f]
					}
				}
			}
			matrix.Set(m, i, val)
		}
	}
	return matrix
}

// F_me intermediate eqaution for CCSD.
func fmeint(spinf, tia *mat.Dense, spin [][][][]float64, nsocc, nmo int) *mat.Dense {
	var val float64
	matrix := mat.NewDense(nmo, nmo, nil)
	for m := 0; m < nsocc; m++ {
		for e := nsocc; e < nmo; e++ {
			val = spinf.At(m, e)
			for n := 0; n < nsocc; n++ {
				for f := nsocc; f < nmo; f++ {
					val += tia.At(n, f) * spin[m][n][e][f]
				}
			}
			matrix.Set(m, e, val)
		}
	}
	return matrix
}

// W_mnij intermediate equation.
func wmnij(spin [][][][]float64, nsocc, nmo int, tau [][][][]float64, tia *mat.Dense) [][][][]float64 {
	array := make([][][][]float64, nmo)
	for m := 0; m < nsocc; m++ {
		array[m] = make([][][]float64, nmo)
		for n := 0; n < nsocc; n++ {
			array[m][n] = make([][]float64, nmo)
			for i := 0; i < nsocc; i++ {
				array[m][n][i] = make([]float64, nmo)
				for j := 0; j < nsocc; j++ {
					array[m][n][i][j] = spin[m][n][i][j]
					for e := nsocc; e < nmo; e++ {
						array[m][n][i][j] += ((tia.At(j, e) * spin[m][n][i][e]) -
							(tia.At(i, e) * spin[m][n][j][e]))
						for f := nsocc; f < nmo; f++ {
							array[m][n][i][j] += 0.25 * tau[i][j][e][f] * spin[m][n][e][f]
						}
					}
				}
			}
		}
	}
	return array
}

// W_abef intermediate equation.
func wabef(spin [][][][]float64, nsocc, nmo int, tau [][][][]float64, tia *mat.Dense) [][][][]float64 {
	array := make([][][][]float64, nmo)
	for a := nsocc; a < nmo; a++ {
		array[a] = make([][][]float64, nmo)
		for b := nsocc; b < nmo; b++ {
			array[a][b] = make([][]float64, nmo)
			for e := nsocc; e < nmo; e++ {
				array[a][b][e] = make([]float64, nmo)
				for f := nsocc; f < nmo; f++ {
					array[a][b][e][f] = spin[a][b][e][f]
					for m := 0; m < nsocc; m++ {
						array[a][b][e][f] -= (tia.At(m, b) * spin[a][m][e][f]) -
							(tia.At(m, a) * spin[b][m][e][f])
						for n := 0; n < nsocc; n++ {
							array[a][b][e][f] += 0.25 * tau[m][n][a][b] * spin[m][n][e][f]
						}
					}
				}
			}
		}
	}
	return array
}

// W_mbej intermediate equation.
func wmbej(spin, tijab [][][][]float64, nsocc, nmo int, tia *mat.Dense) [][][][]float64 {
	array := make([][][][]float64, nmo)
	for m := 0; m < nsocc; m++ {
		array[m] = make([][][]float64, nmo)
		for b := nsocc; b < nmo; b++ {
			array[m][b] = make([][]float64, nmo)
			for e := nsocc; e < nmo; e++ {
				array[m][b][e] = make([]float64, nmo)
				for j := 0; j < nsocc; j++ {
					array[m][b][e][j] = spin[m][b][e][j]
					for f := nsocc; f < nmo; f++ {
						array[m][b][e][j] += tia.At(j, f) * spin[m][b][e][f]
						for n := 0; n < nsocc; n++ {
							array[m][b][e][j] -= ((0.5 * tijab[j][n][f][b]) + (tia.At(j, f) * tia.At(n, b))) * spin[m][n][e][f]
						}
					}
					for n := 0; n < nsocc; n++ {
						array[m][b][e][j] -= tia.At(n, b) * spin[m][n][e][j]
					}
				}
			}
		}
	}
	return array
}

// Singles denominator term for T1 amplitudes.
func dia(spinf *mat.Dense, nmo, nsocc int) *mat.Dense {
	var val float64
	r := nmo
	matrix := mat.NewDense(r, r, nil)
	for i := 0; i < nsocc; i++ {
		for a := nsocc; a < nmo; a++ {
			val = spinf.At(i, i) - spinf.At(a, a)
			matrix.Set(i, a, val)
		}
	}
	return matrix
}

// Doubles denominator term for T2 amplitudes.
func dijab(spinf *mat.Dense, nmo, nsocc int) [][][][]float64 {
	val := make([][][][]float64, nmo)
	for i := 0; i < nsocc; i++ {
		val[i] = make([][][]float64, nmo)
		for j := 0; j < nsocc; j++ {
			val[i][j] = make([][]float64, nmo)
			for a := nsocc; a < nmo; a++ {
				val[i][j][a] = make([]float64, nmo)
				for b := nsocc; b < nmo; b++ {
					val[i][j][a][b] = spinf.At(i, i) + spinf.At(j, j) - spinf.At(a, a) - spinf.At(b, b)
				}
			}
		}
	}
	return val
}

// T1 amplitude equation for CCSD.
func T1amp(spinf, fae, fmi, fme, dia, tia *mat.Dense, tijab, spin [][][][]float64, nsocc, nmo int) *mat.Dense {
	var val float64
	matrix := mat.NewDense(nmo, nmo, nil)
	for i := 0; i < nsocc; i++ {
		for a := nsocc; a < nmo; a++ {
			val = spinf.At(i, a)
			for m := 0; m < nsocc; m++ {
				val -= tia.At(m, a) * fmi.At(m, i)
			}
			for e := nsocc; e < nmo; e++ {
				val += tia.At(i, e) * fae.At(a, e)
				for m := 0; m < nsocc; m++ {
					val += tijab[i][m][a][e] * fme.At(m, e)
					for f := nsocc; f < nmo; f++ {
						val -= 0.5 * tijab[i][m][e][f] * spin[m][a][e][f]
					}
					for n := 0; n < nsocc; n++ {
						val -= 0.5 * tijab[m][n][a][e] * spin[n][m][e][i]
					}
				}
			}
			for n := 0; n < nsocc; n++ {
				for f := nsocc; f < nmo; f++ {
					val -= tia.At(n, f) * spin[n][a][i][f]
				}
			}
			val /= dia.At(i, a)
			matrix.Set(i, a, val)
		}
	}
	return matrix
}

// T2 amplitude equation for CCSD.
func T2amp(spinf, fae, fmi, fme, tia *mat.Dense, tijab, spin, wmnij, wabef, wmbej, dijab, tau [][][][]float64, nsocc, nmo int) [][][][]float64 {
	val := make([][][][]float64, nmo)
	for i := 0; i < nsocc; i++ {
		val[i] = make([][][]float64, nmo)
		for j := 0; j < nsocc; j++ {
			val[i][j] = make([][]float64, nmo)
			for a := nsocc; a < nmo; a++ {
				val[i][j][a] = make([]float64, nmo)
				for b := nsocc; b < nmo; b++ {
					val[i][j][a][b] = spin[i][j][a][b]
					for e := nsocc; e < nmo; e++ {
						summ1, summ2 := 0.0, 0.0
						for m := 0; m < nsocc; m++ {
							summ1 += tia.At(m, b) * fme.At(m, e)
							summ2 += tia.At(m, a) * fme.At(m, e)
						}
						val[i][j][a][b] += tijab[i][j][a][e] * (fae.At(b, e) - 0.5*summ1)
						val[i][j][a][b] -= tijab[i][j][b][e] * (fae.At(a, e) - 0.5*summ2)
					}
					for m := 0; m < nsocc; m++ {
						sume2, sume3 := 0.0, 0.0
						for e := nsocc; e < nmo; e++ {
							sume2 += tia.At(j, e) * fme.At(m, e)
							sume3 += tia.At(i, e) * fme.At(m, e)
						}
						val[i][j][a][b] -= (tijab[i][m][a][b] * (fmi.At(m, j) + (0.5 * sume2))) - (tijab[j][m][a][b] * (fmi.At(m, i) + (0.5 * sume3)))
						for n := 0; n < nsocc; n++ {
							val[i][j][a][b] += 0.5 * tau[m][n][a][b] * wmnij[m][n][i][j]
						}
					}
					for e := nsocc; e < nmo; e++ {
						for f := nsocc; f < nmo; f++ {
							val[i][j][a][b] += 0.5 * tau[i][j][e][f] * wabef[a][b][e][f]
						}
					}
					for m := 0; m < nsocc; m++ {
						for e := nsocc; e < nmo; e++ {
							val[i][j][a][b] += tijab[i][m][a][e]*wmbej[m][b][e][j] - tia.At(i, e)*tia.At(m, a)*spin[m][b][e][j]
							val[i][j][a][b] -= tijab[i][m][b][e]*wmbej[m][a][e][j] - tia.At(i, e)*tia.At(m, b)*spin[m][a][e][j]
							val[i][j][a][b] -= tijab[j][m][a][e]*wmbej[m][b][e][i] - tia.At(j, e)*tia.At(m, a)*spin[m][b][e][i]
							val[i][j][a][b] += tijab[j][m][b][e]*wmbej[m][a][e][i] - tia.At(j, e)*tia.At(m, b)*spin[m][a][e][i]
						}
					}
					for e := nsocc; e < nmo; e++ {
						val[i][j][a][b] += ((tia.At(i, e) * spin[a][b][e][j]) - (tia.At(j, e) * spin[a][b][e][i]))
					}
					for m := 0; m < nsocc; m++ {
						val[i][j][a][b] -= ((tia.At(m, a) * spin[m][b][i][j]) - (tia.At(m, b) * spin[m][a][i][j]))
					}
					val[i][j][a][b] /= dijab[i][j][a][b]
				}
			}
		}
	}
	return val
}

func ccenergy(sof, tia *mat.Dense, spin, tijab [][][][]float64, nsocc, nmo int) float64 {
	var val, sumia, sumijab1, sumijab2 float64
	for i := 0; i < nsocc; i++ {
		for a := nsocc; a < nmo; a++ {
			for j := 0; j < nsocc; j++ {
				for b := nsocc; b < nmo; b++ {
					sumijab2 += spin[i][j][a][b] * tia.At(i, a) * tia.At(j, b)
					sumijab1 += spin[i][j][a][b] * tijab[i][j][a][b]
				}
			}
			sumia += sof.At(i, a) * tia.At(i, a)
		}
	}
	val = sumia + (0.25 * sumijab1) + (0.5 * sumijab2)
	return val
}

func dijkabc(sof *mat.Dense, nmo, nsocc int) [][][][][][]float64 {
	val := make([][][][][][]float64, nmo)
	for i := 0; i < nsocc; i++ {
		val[i] = make([][][][][]float64, nmo)
		for j := 0; j < nsocc; j++ {
			val[i][j] = make([][][][]float64, nmo)
			for k := 0; k < nsocc; k++ {
				val[i][j][k] = make([][][]float64, nmo)
				for a := nsocc; a < nmo; a++ {
					val[i][j][k][a] = make([][]float64, nmo)
					for b := nsocc; b < nmo; b++ {
						val[i][j][k][a][b] = make([]float64, nmo)
						for c := nsocc; c < nmo; c++ {
							val[i][j][k][a][b][c] = sof.At(i, i) + sof.At(j, j) + sof.At(k, k) - sof.At(a, a) - sof.At(b, b) - sof.At(c, c)
						}
					}
				}
			}
		}
	}
	return val
}

func t3con(d3 [][][][][][]float64, spin, tijab [][][][]float64, nmo, nsocc int) [][][][][][]float64 {
	val := make([][][][][][]float64, nmo)
	for i := 0; i < nsocc; i++ {
		val[i] = make([][][][][]float64, nmo)
		for j := 0; j < nsocc; j++ {
			val[i][j] = make([][][][]float64, nmo)
			for k := 0; k < nsocc; k++ {
				val[i][j][k] = make([][][]float64, nmo)
				for a := nsocc; a < nmo; a++ {
					val[i][j][k][a] = make([][]float64, nmo)
					for b := nsocc; b < nmo; b++ {
						val[i][j][k][a][b] = make([]float64, nmo)
						for c := nsocc; c < nmo; c++ {
							for e := nsocc; e < nmo; e++ {
								val[i][j][k][a][b][c] += tijab[j][k][a][e] * spin[e][i][b][c]
								val[i][j][k][a][b][c] -= tijab[j][k][b][e] * spin[e][i][a][c]
								val[i][j][k][a][b][c] -= tijab[j][k][c][e] * spin[e][i][b][a]
								val[i][j][k][a][b][c] -= tijab[i][k][a][e] * spin[e][j][b][c]
								val[i][j][k][a][b][c] += tijab[i][k][b][e] * spin[e][j][a][c]
								val[i][j][k][a][b][c] += tijab[i][k][c][e] * spin[e][j][b][a]
								val[i][j][k][a][b][c] -= tijab[j][i][a][e] * spin[e][k][b][c]
								val[i][j][k][a][b][c] += tijab[j][i][b][e] * spin[e][k][a][c]
								val[i][j][k][a][b][c] += tijab[j][i][c][e] * spin[e][k][b][a]
							}
							for m := 0; m < nsocc; m++ {
								val[i][j][k][a][b][c] -= tijab[i][m][b][c] * spin[m][a][j][k]
								val[i][j][k][a][b][c] += tijab[i][m][a][c] * spin[m][b][j][k]
								val[i][j][k][a][b][c] += tijab[i][m][b][a] * spin[m][c][j][k]
								val[i][j][k][a][b][c] += tijab[j][m][b][c] * spin[m][a][i][k]
								val[i][j][k][a][b][c] -= tijab[j][m][a][c] * spin[m][b][i][k]
								val[i][j][k][a][b][c] -= tijab[j][m][b][a] * spin[m][c][i][k]
								val[i][j][k][a][b][c] += tijab[k][m][b][c] * spin[m][a][j][i]
								val[i][j][k][a][b][c] -= tijab[k][m][a][c] * spin[m][b][j][i]
								val[i][j][k][a][b][c] -= tijab[k][m][b][a] * spin[m][c][j][i]
							}
						}
					}
				}
			}
		}
	}
	return val
}

func t3dis(d3 [][][][][][]float64, spin [][][][]float64, tia *mat.Dense, nmo, nsocc int) [][][][][][]float64 {
	val := make([][][][][][]float64, nmo)
	for i := 0; i < nsocc; i++ {
		val[i] = make([][][][][]float64, nmo)
		for j := 0; j < nsocc; j++ {
			val[i][j] = make([][][][]float64, nmo)
			for k := 0; k < nsocc; k++ {
				val[i][j][k] = make([][][]float64, nmo)
				for a := nsocc; a < nmo; a++ {
					val[i][j][k][a] = make([][]float64, nmo)
					for b := nsocc; b < nmo; b++ {
						val[i][j][k][a][b] = make([]float64, nmo)
						for c := nsocc; c < nmo; c++ {
							val[i][j][k][a][b][c] += tia.At(i, a)*spin[j][k][b][c] - tia.At(i, b)*spin[j][k][a][c] - tia.At(i, c)*spin[j][k][b][a]
							val[i][j][k][a][b][c] -= tia.At(j, a)*spin[i][k][b][c] - tia.At(j, b)*spin[i][k][a][c] - tia.At(j, c)*spin[i][k][b][a]
							val[i][j][k][a][b][c] -= tia.At(k, a)*spin[j][i][b][c] - tia.At(k, b)*spin[j][i][a][c] - tia.At(k, c)*spin[j][i][b][a]
						}
					}
				}
			}
		}
	}
	return val
}

func ept(d3, t3discon, t3conn [][][][][][]float64, nmo, nsocc int) float64 {
	var val float64
	for i := 0; i < nsocc; i++ {
		for j := 0; j < nsocc; j++ {
			for k := 0; k < nsocc; k++ {
				for a := nsocc; a < nmo; a++ {
					for b := nsocc; b < nmo; b++ {
						for c := nsocc; c < nmo; c++ {
							val += t3conn[i][j][k][a][b][c] / d3[i][j][k][a][b][c] * (t3conn[i][j][k][a][b][c] + t3discon[i][j][k][a][b][c])
						}
					}
				}
			}
		}
	}
	return val / 36.0
}

func smat(matrix *mat.SymDense) *mat.Dense {
	// fmt.Println("SMAT")
	var eig mat.EigenSym
	var evec mat.Dense
	eig.Factorize(matrix, true)
	eval := eig.Values(nil)
	lambda := mat.NewDense(len(eval), len(eval), nil)
	for i := range eval {
		lambda.Set(i, i, 1/math.Sqrt(eval[i]))
	}
	eig.VectorsTo(&evec)
	var s12 mat.Dense
	s12.Mul(lambda, evec.T())
	s12.Mul(&evec, &s12)
	// printmat(&s12)
	return &s12
}

func ToSym(matrix *mat.Dense) *mat.SymDense {
	r, c := matrix.Dims()
	symmat := mat.NewSymDense(r, nil)
	for i := 0; i < r; i++ {
		for j := i; j < c; j++ {
			symmat.SetSym(i, j, matrix.At(i, j))
		}
	}
	return symmat
}

func cmat(h, s12 *mat.SymDense) *mat.Dense {
	var fmat mat.Dense
	fmat.Mul(h, s12)
	fmat.Mul(s12.T(), &fmat)
	fsmat := ToSym(&fmat)
	var eig mat.EigenSym
	var evec mat.Dense
	eig.Factorize(fsmat, true)
	eval := eig.Values(nil)
	energy := mat.NewDense(len(eval), len(eval), nil)
	for i := range eval {
		energy.Set(i, i, eval[i])
	}
	eig.VectorsTo(&evec)
	var c0 mat.Dense
	c0.Mul(s12, &evec)
	return &c0
}

func denmat(c0 *mat.Dense, nocc int) *mat.Dense {
	r, c := c0.Dims()
	d := mat.NewDense(r, c, nil)
	for mu := 0; mu < r; mu++ {
		for nu := 0; nu < c; nu++ {
			var sum float64
			for m := 0; m < nocc; m++ {
				sum += c0.At(mu, m) * c0.At(nu, m)
			}
			d.Set(mu, nu, sum)
		}
	}
	return d
}

func mofockmat(c0, m mat.Matrix) []float64 {
	r, _ := c0.Dims()
	mofock := mat.NewSymDense(r, nil)
	for i := 0; i < r; i++ {
		for j := 0; j < r; j++ {
			var sum float64
			for mu := 0; mu < r; mu++ {
				for nu := 0; nu < r; nu++ {
					sum += c0.At(mu, j) * c0.At(nu, i) * m.At(mu, nu)
				}
				mofock.SetSym(i, j, sum)
			}
		}
	}
	var eig mat.EigenSym
	eig.Factorize(mofock, true)
	eval := eig.Values(nil)
	// fmt.Println(eval)
	return eval
}

func mohmat(c0, m mat.Matrix) *mat.SymDense {
	r, _ := c0.Dims()
	moh := mat.NewSymDense(r, nil)
	for i := 0; i < r; i++ {
		for j := 0; j < r; j++ {
			var sum float64
			for mu := 0; mu < r; mu++ {
				for nu := 0; nu < r; nu++ {
					sum += c0.At(mu, j) * c0.At(nu, i) * m.At(mu, nu)
				}
				moh.SetSym(i, j, sum)
			}
		}
	}
	return moh
}

func CISHam(sof *mat.Dense, spin [][][][]float64, nsocc, nmo int) *mat.SymDense {
	r := 40
	var val float64
	cham := mat.NewSymDense(r, nil)
	for i := 0; i < nsocc; i++ {
		for a := nsocc; a < nmo; a++ {
			for j := 0; j < nsocc; j++ {
				for b := nsocc; b < nmo; b++ {
					val = sof.At(a, b)*kronk(i, j) - sof.At(i, j)*kronk(a, b) + spin[a][j][i][b]
					cham.SetSym(i, j, val)
				}
			}
		}
	}
	return cham
}

func main() {
	nucrep := enuc("h2o_sto3g/enuc.dat")
	s := maketrix("h2o_sto3g/s.dat")
	t := maketrix("h2o_sto3g/t.dat")
	v := maketrix("h2o_sto3g/v.dat")
	nocc := 5
	nmo := 14
	nsocc := 10
	h := madd(t, v)
	s12 := ToSym(smat(s))
	c := cmat(h, s12)
	d := denmat(c, nocc)
	eri := twoelec("h2o_sto3g/eri.dat")
	f := fockmat(h, d, eri)
	// fmt.Println(eri)
	// // fmt.Println("Nuclear Repulsion Energy = ", nucrep, "\n")
	// // fmt.Println("AO Overlap Integrals:", "\n")
	// // printmat(s)
	// fmt.Println("\n")
	// fmt.Println("Kinetic Energy Integrals:", "\n")
	// // printmat(t)
	// fmt.Println("\n")
	// fmt.Println("Nuclear Attraction Integrals:", "\n")
	// // printmat(v)
	// fmt.Println("\n")
	// fmt.Println("Core Hamiltonian:", "\n")
	// // printmat(h)
	// twoelec("h2o_sto3g/eri.dat")
	// fmt.Println("\n")
	// fmt.Println("Orthogonalization Matrix:", "\n")
	// // printmat(s12)
	// fmt.Println("\n")
	// fmt.Println("Initial Guess Density Matrix:", "\n")
	// // printmat(d)
	// fmt.Println("\n")
	scfinit := inite(d, h, h)
	// fmt.Println("Initial SCF Energy = ", scfinit, "\n")
	// fmt.Println("New Fock Matrix:", "\n")
	// // printmat(f)
	// fmt.Println("\n")
	fmt.Println("Iter:", " E(elec):", "            E(tot):", "              Delta(E):", "            RMS(D):")
	iter, de, rmsd := 1, 1.0, 1.0
	var newc *mat.Dense
	var newe float64
	var etot float64
	for de > 1e-12 && rmsd > 1e-12 {
		newc = cmat(f, s12)
		newden := denmat(newc, nocc)
		newe = inite(newden, h, f)
		f = fockmat(h, newden, eri)
		de = math.Abs(newe - scfinit)
		scfinit = newe
		var diff mat.Dense
		diff.Sub(newden, d)
		rmsd = mat.Norm(&diff, 2)
		d = newden
		etot = newe + nucrep
		fmt.Printf("%02d %20.12f %20.12f %20.12f %20.12f\n", iter, newe, etot, de, rmsd)
		iter++
	}
	fmt.Println("\n")
	EHF := etot
	fmt.Println("HF/STO-3G Energy = ", EHF, "\n")
	mof := mofockmat(newc, f)
	// nod := noddy(newc, eri)
	// fmt.Println(mp2energy(mo, nod, nocc, nao))
	newi := smart(newc, eri)
	mp2 := mp2energy(mof, newi, nocc)
	Emp2 := EHF + mp2
	fmt.Println("MP2 Energy Correction = ", mp2, "\n")
	fmt.Println("MP2/STO-3G Energy = ", Emp2, "\n")
	spin := spattospin(newi, nmo)
	moh := mohmat(newc, h)
	sof := spinfock(moh, spin, nmo, nsocc)
	singlestia := tia(nsocc, nmo)
	doublestijab := tijab(spin, sof, nsocc, nmo)
	CCmp2 := ccmp2(sof, spin, nmo, nsocc, doublestijab)
	fmt.Println("MP2 Energy Correction from CC Initial Amplitudes = ", CCmp2, "\n")
	tildetau := ttau(spin, doublestijab, nmo, nsocc, singlestia)
	regtau := tau(spin, doublestijab, nmo, nsocc, singlestia)
	fae := faeint(sof, singlestia, spin, tildetau, nsocc, nmo)
	fmi := fmiint(sof, singlestia, spin, tildetau, nsocc, nmo)
	fme := fmeint(sof, singlestia, spin, nsocc, nmo)
	Wmnij := wmnij(spin, nsocc, nmo, regtau, singlestia)
	Wabef := wabef(spin, nsocc, nmo, regtau, singlestia)
	Wmbej := wmbej(spin, doublestijab, nsocc, nmo, singlestia)
	Dia := dia(sof, nmo, nsocc)
	Dijab := dijab(sof, nmo, nsocc)
	T1 := T1amp(sof, fae, fmi, fme, Dia, singlestia, doublestijab, spin, nsocc, nmo)
	T2 := T2amp(sof, fae, fmi, fme, singlestia, doublestijab, spin, Wmnij, Wabef, Wmbej, Dijab, regtau, nsocc, nmo)
	fmt.Println("Iter:", "  E(elec):", "           CC-E(tot):", "            Delta(E):", "           T1 RMS(D):")
	iter, de, rmsd = 1, 1.0, 1.0
	// rmsd2 := 1.0
	var firstEcc float64
	var CCEtot float64
	var newt1 *mat.Dense
	var newt2 [][][][]float64
	for de > 1e-12 && rmsd > 1e-12 {
		tildetau := ttau(spin, T2, nmo, nsocc, T1)
		regtau := tau(spin, T2, nmo, nsocc, T1)
		fae := faeint(sof, T1, spin, tildetau, nsocc, nmo)
		fmi := fmiint(sof, T1, spin, tildetau, nsocc, nmo)
		fme := fmeint(sof, T1, spin, nsocc, nmo)
		Wmnij := wmnij(spin, nsocc, nmo, regtau, T1)
		Wabef := wabef(spin, nsocc, nmo, regtau, T1)
		Wmbej := wmbej(spin, T2, nsocc, nmo, T1)
		Dia := dia(sof, nmo, nsocc)
		Dijab := dijab(sof, nmo, nsocc)
		newt1 = T1amp(sof, fae, fmi, fme, Dia, T1, T2, spin, nsocc, nmo)
		newt2 = T2amp(sof, fae, fmi, fme, T1, T2, spin, Wmnij, Wabef, Wmbej, Dijab, regtau, nsocc, nmo)
		newEcc := ccenergy(sof, T1, spin, T2, nsocc, nmo)
		de = math.Abs(newEcc - firstEcc)
		firstEcc = newEcc
		var diff mat.Dense
		diff.Sub(T1, newt1)
		rmsd = mat.Norm(&diff, 2)
		T1 = newt1
		T2 = newt2
		CCEtot = newEcc + EHF
		fmt.Printf("%02d %20.12f %20.12f %20.12f %20.12f\n", iter, newEcc, CCEtot, de, rmsd)
		iter++
	}
	fmt.Println("\n")
	fmt.Println("HF/STO-3G Energy = ", EHF, "\n")
	fmt.Println("CCSD Energy Correction = ", firstEcc, "\n")
	fmt.Println("CCSD/STO-3G Energy = ", CCEtot, "\n")
	d3 := dijkabc(sof, nmo, nsocc)
	t3conn := t3con(d3, spin, newt2, nmo, nsocc)
	t3disc := t3dis(d3, spin, newt1, nmo, nsocc)
	eparant := ept(d3, t3disc, t3conn, nmo, nsocc)
	fmt.Println("(T) Energy Correction = ", eparant, "\n")
	fmt.Println("CCSD(T)/STO-3G Energy = ", EHF+firstEcc+eparant, "\n")
	cham := CISHam(sof, spin, nsocc, nmo)
	printmat(cham)
}

// (progn (setq compile-command "time go run .") (recompile))
