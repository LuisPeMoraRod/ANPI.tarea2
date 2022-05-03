package main

import (
	"fmt"
	"runtime"
	"time"

	tridiagonal "../lib"
	secuential "../parte1_p2/secuential"
	parallel "../parte1_p3/parallel"
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

	b := allEqualVec(m, 1)   //constant terms vector
	x_0 := allEqualVec(m, 0) //init solution vector
	tol := 1e-5
	iterMax := 1000
	processors := runtime.NumCPU() //number of available processors

	var time_sec time.Duration
	time_sec = getTs(A, b, x_0, iterMax, tol) //get time for secuential implementation

	var time_par time.Duration
	accels := make([]float64, processors) //create vector to store acceleration ratios
	for i := 0; i < processors; i++ {
		accels[i] = 0
	}

	for i := 1; i <= processors; i++ {
		x_0 := allEqualVec(m, 0) //init solution vector
		time_par = getTp(A, b, x_0, iterMax, tol, i)
		accel := float64(time_sec) / float64(time_par)
		accels[i-1] = accel
	}

	for i := 0; i < processors; i++ {
		fmt.Println("S(", i+1, ") = ", accels[i])
	}

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

//Calculates the time taken by parallel implementation with  @p processes
//Parameters:
//@A : coefficient matrix
//@b : constant terms vector
//@x_0 : inital solution vector
//@iterMax : max number of possible iterations
//@tol : error tolerance
//@p : number of processes for parallel implementation
//Output : time taken to get the solution
func getTp(A *mat.Dense, b *mat.VecDense, x_0 *mat.VecDense, iterMax int, tol float64, p int) time.Duration {
	var time_par time.Duration

	start := time.Now() //start timer for secuential
	x_k, k, err := parallel.ParallelJacobi(A, b, x_0, iterMax, tol, p)
	fmt.Println("x_k = ")
	parallel.MatPrint(x_k)
	fmt.Println("k = ", k)
	fmt.Println("error = ", err)
	time_par = time.Since(start) //stop timer for secuential
	return time_par
}

//Calculates the time taken by secuential implementation
//Parameters:
//@A : coefficient matrix
//@b : constant terms vector
//@x_0 : inital solution vector
//@iterMax : max number of possible iterations
//@tol : error tolerance
//Output : time taken to get the solution
func getTs(A *mat.Dense, b *mat.VecDense, x_0 *mat.VecDense, iterMax int, tol float64) time.Duration {
	var time_sec time.Duration

	start := time.Now() //start timer for secuential
	x_k, k, err := secuential.SecJacobi(A, b, x_0, iterMax, tol)
	fmt.Println("x_k = ")
	secuential.MatPrint(x_k)
	fmt.Println("k = ", k)
	fmt.Println("error = ", err)
	time_sec = time.Since(start) //stop timer for secuential
	return time_sec
}
