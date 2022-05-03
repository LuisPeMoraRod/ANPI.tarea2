package main

import (
	"fmt"

	tridiagonal "../lib"
	secuential "./secuential"
	"gonum.org/v1/gonum/mat"
)

func main() {

	m := 242    // size of matrix
	step := 0.1 //step of q and p vectors elements

	//create p and q vectors
	p_values := make([]float64, m-1)
	for i := 0; i < m-1; i++ {
		p_values[i] = float64(1 + float64(i)*float64(step))
	}
	p := mat.NewVecDense(m-1, p_values)
	q := mat.NewVecDense(m-1, p_values)

	A := tridiagonal.Tridiagonal(*p, *q, m) //set coefficient matrix (tridiagonally dominant)

	b := secuential.AllEqualVec(m, 1)   //constant terms vector
	x_0 := secuential.AllEqualVec(m, 0) //init solution vector
	tol := 1e-5
	iterMax := 1000

	x_k, k, err := secuential.SecJacobi(A, b, x_0, iterMax, tol)
	fmt.Println("x_k = ")
	secuential.MatPrint(x_k)
	fmt.Println("k = ", k)
	fmt.Println("error = ", err)

}
