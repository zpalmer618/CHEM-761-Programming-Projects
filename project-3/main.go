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

func denmat(h, s12 *mat.SymDense, nocc int) *mat.Dense {
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

func main() {
	nucrep := enuc("h2o_sto3g/enuc.dat")
	s := maketrix("h2o_sto3g/s.dat")
	t := maketrix("h2o_sto3g/t.dat")
	v := maketrix("h2o_sto3g/v.dat")
	nocc := 5
	h := madd(t, v)
	s12 := ToSym(smat(s))
	d := denmat(h, s12, nocc)
	eri := twoelec("h2o_sto3g/eri.dat")
	f := fockmat(h, d, eri)
	fmt.Println("Nuclear Repulsion Energy = ", nucrep, "\n")
	fmt.Println("AO Overlap Integrals:", "\n")
	printmat(s)
	fmt.Println("\n")
	fmt.Println("Kinetic Energy Integrals:", "\n")
	printmat(t)
	fmt.Println("\n")
	fmt.Println("Nuclear Attraction Integrals:", "\n")
	printmat(v)
	fmt.Println("\n")
	fmt.Println("Core Hamiltonian:", "\n")
	printmat(h)
	twoelec("h2o_sto3g/eri.dat")
	fmt.Println("\n")
	fmt.Println("Orthogonalization Matrix:", "\n")
	printmat(s12)
	fmt.Println("\n")
	fmt.Println("Initial Guess Density Matrix:", "\n")
	printmat(d)
	fmt.Println("\n")
	scfinit := inite(d, h, h)
	fmt.Println("Initial SCF Energy = ", scfinit, "\n")
	fmt.Println("New Fock Matrix:", "\n")
	printmat(f)
	fmt.Println("\n")
	fmt.Println("Iter:", " E(elec):", "            E(tot):", "              Delta(E):", "            RMS(D):")
	iter, de, rmsd := 1, 1.0, 1.0
	for de > 1e-12 && rmsd > 1e-12 {
		newden := denmat(f, s12, nocc)
		newe := inite(newden, h, f)
		f = fockmat(h, newden, eri)
		de = math.Abs(newe - scfinit)
		scfinit = newe
		var diff mat.Dense
		diff.Sub(newden, d)
		rmsd = mat.Norm(&diff, 2)
		d = newden
		etot := newe + nucrep
		fmt.Printf("%02d %20.12f %20.12f %20.12f %20.12f\n", iter, newe, etot, de, rmsd)
		iter++
	}
}
