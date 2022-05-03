package main

import (
	"fmt"
	"runtime"

	tridiagonal "../lib"
	parallel "./parallel"
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

	b := parallel.AllEqualVec(m, 1)   //constant terms vector
	x_0 := parallel.AllEqualVec(m, 0) //init solution vector
	tol := 1e-5
	iterMax := 1000
	processors := runtime.NumCPU() //number of available processors

	x_k, k, err := parallel.ParallelJacobi(A, b, x_0, iterMax, tol, processors)
	fmt.Println("x_k = ")
	parallel.MatPrint(x_k)
	fmt.Println("k = ", k)
	fmt.Println("error = ", err)

}
