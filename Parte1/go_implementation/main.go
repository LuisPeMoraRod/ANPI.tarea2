package main

import (
	"fmt"
	"log"

	"gonum.org/v1/gonum/mat"
)

func main() {

	m := 10 // size of matrix
	step := 1

	//create p and q vectors
	p_values := make([]float64, m-1)
	for i := 0; i < m-1; i++ {
		p_values[i] = float64(1 + float64(i)*float64(step))
	}
	p := mat.NewVecDense(m-1, p_values)
	q := mat.NewVecDense(m-1, p_values)

	tridiagonal(*p, *q, m)

}

func matPrint(X mat.Matrix) {
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%v\n", fa)
}

// Checks if size of vector is m-1
func correctLen(v mat.VecDense, m int) bool {
	if v.Len() == m-1 {
		return true
	}
	return false
}

//Sets diagonal of the matrix
func setDiag(p mat.VecDense, q mat.VecDense, m int) mat.Matrix {
	d_values := make([]float64, m)
	d_values[0] = float64(2) * q.AtVec(0)
	for i := 1; i <= m-2; i++ {
		d_values[i] = float64(2) * (p.AtVec(i-1) + q.AtVec(i))
	}
	d_values[m-1] = p.AtVec(m - 2)

	D := mat.NewDiagDense(m, d_values) //Diagonal matrix
	return D
}

//Returns square matrix of size @m with first sub diagonal setted by @v vector
func firstSubDiag(m int, v mat.VecDense) *mat.Dense {

	M := zeroMat(m)
	for i := 1; i < m; i++ {
		M.Set(i, i-1, v.At(i-1, 0))
	}
	return M
}

//Returns square matrix of size @m with first super diagonal setted by @v vector
func firstSuperDiag(m int, v mat.VecDense) *mat.Dense {

	M := zeroMat(m)
	for i := 0; i < m-1; i++ {
		M.Set(i, i+1, v.At(i, 0))
	}
	return M
}

//Returns square zero matrix of size @m
func zeroMat(m int) *mat.Dense {
	zero_v := make([]float64, m*m)
	for i := 0; i < m*m; i++ {
		zero_v[i] = 0
	}
	Z := mat.NewDense(m, m, zero_v)
	return Z
}

//Generates a tridiagonal mxm matrix, based on the requirements specifed
//parameters:
//   @p = vector of length m
//   @p = vector of length m
//   @m = size of the square matrix
//output:
//   resultant tridiagonal matrix
func tridiagonal(p mat.VecDense, q mat.VecDense, m int) mat.Matrix {
	if m >= 3 {
		if correctLen(p, m) && correctLen(q, m) {

			D := setDiag(p, q, m)       //create diagonal
			Dm1 := firstSubDiag(m, p)   //create first sub diagonal
			Dp1 := firstSuperDiag(m, q) //create first super diagonal

			M := zeroMat(m)
			M.Add(Dm1, D)
			M.Add(M, Dp1)

			return M
		} else {
			log.Fatal("p and q vector must match the size: m-1")
			return nil
		}

	} else {
		log.Fatal("m must be greater or equal than 3")
		return nil
	}

}
