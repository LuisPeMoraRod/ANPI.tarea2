package main

import (
	"fmt"
	"log"

	tridiagonal "./lib"
	"gonum.org/v1/gonum/mat"
)

func main() {

	m := 10   // size of matrix
	step := 1 //step of q and p vectors elements

	//create p and q vectors
	p_values := make([]float64, m-1)
	for i := 0; i < m-1; i++ {
		p_values[i] = float64(1 + float64(i)*float64(step))
	}
	p := mat.NewVecDense(m-1, p_values)
	q := mat.NewVecDense(m-1, p_values)

	A := tridiagonal.Tridiagonal(*p, *q, m)
	matPrint(A)

	b := allEqualVec(m, 1)
	x_0 := allEqualVec(m, 0)
	tol := 1e-5
	iterMax := 1000

	x_k, k, err := secJacobi(A, b, x_0, iterMax, tol)
	matPrint(x_k)
	fmt.Println(k)
	fmt.Println(err)
}

//Returns pointer to a vector of size @m and all values equal to @value
func allEqualVec(m int, value float64) *mat.VecDense {
	v_values := make([]float64, m)
	for i := 0; i < m; i++ {
		v_values[i] = value
	}
	v := mat.NewVecDense(m, v_values)
	return v
}

//Formats matrix to print it in console
func matPrint(X mat.Matrix) {
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%v\n", fa)
}

//Secuential implementation of Jacobi's method
//Parameters:
//@A : coefficient matrix
//@b : constant terms vector
//@x_0 : inital solution vector
//@iterMax : max number of possible iterations
//@tol : error tolerance
//Output
//1) Solution vector
//2) Number of iterations taken
//3) Gotten error
func secJacobi(A *mat.Dense, b *mat.VecDense, x_0 *mat.VecDense, iterMax int, tol float64) (*mat.VecDense, int, float64) {

	m, n := A.Caps() //get number of rows and collumns of the matrix

	//Returns x_i for k iteration
	getX_i := func(i int, x_k *mat.VecDense) float64 {
		var sum_i float64
		sum_i = 0
		for j := 0; j < m; j++ {
			if j != i {
				sum_i += A.At(i, j) * x_k.AtVec(j)
			}
		}
		x_i := 1 / A.At(i, i) * (b.AtVec(i) - sum_i)
		return x_i
	}

	if m == n {
		if tridiagonal.IsDiagDominant(A) {
			x_k := x_0 //initialize solution vector

			var x_i float64
			var err float64
			sol := mat.NewDense(m, 1, nil)
			for k := 0; k < iterMax; k++ {
				x_kp1 := allEqualVec(m, 0) //intialize x_(k+1) with zero values
				for i := 0; i < m; i++ {   //set all values of x_(k+1)
					x_i = getX_i(i, x_k)
					x_kp1.SetVec(i, x_i)
				}
				x_k = x_kp1
				sol.Product(A, x_k)
				sol.Sub(sol, b)
				err = sol.Norm(2)
				// fmt.Println(err)
				if err < tol {
					return x_k, k, err
				}
			}
			return x_k, iterMax, err

		} else {
			log.Fatal("Matrix must be diagonally dominant")
			return nil, 0, 0
		}
	} else {
		log.Fatal("Matrix must be square")
		return nil, 0, 0
	}

}
