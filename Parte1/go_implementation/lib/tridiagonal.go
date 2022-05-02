package tridiagonal

import (
	"log"
	"math"

	"gonum.org/v1/gonum/mat"
)

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
	d_values[m-1] = float64(2) * p.AtVec(m-2)

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
//   @q = vector of length m
//   @m = size of the square matrix
//output:
//   resultant tridiagonal matrix
func Tridiagonal(p mat.VecDense, q mat.VecDense, m int) *mat.Dense {
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

//Checks if matrix is diagonally dominant
func IsDiagDominant(A *mat.Dense) bool {
	return isDiagDomByRows(A) && isDiagDomByCols(A)
}

//Checks if matrix is diagonally dominant by rows
func isDiagDomByRows(A *mat.Dense) bool {
	m, n := A.Caps()
	var sumRow float64
	isValid := true

	for i := 0; i < m; i++ {
		diagVal := math.Abs(A.At(i, i))
		sumRow = 0
		for j := 0; j < n; j++ {
			if i != j {
				sumRow += A.At(i, j)
			}
		}
		if diagVal <= sumRow { //Abs(A[i,i]) must be greater than sum of the rest of values in row i
			isValid = false
		}
	}

	return isValid
}

//Checks if matrix is diagonally dominant by collumns
func isDiagDomByCols(A *mat.Dense) bool {
	m, n := A.Caps()
	var sumCol float64
	isValid := true

	for j := 0; j < m; j++ {
		diagVal := math.Abs(A.At(j, j))
		sumCol = 0
		for i := 0; i < n; i++ {
			if i != j {
				sumCol += A.At(i, j)
			}
		}
		if diagVal <= sumCol { //Abs(A[j,j]) must be greater than sum of the rest of values in collumn j
			isValid = false
		}
	}

	return isValid
}
