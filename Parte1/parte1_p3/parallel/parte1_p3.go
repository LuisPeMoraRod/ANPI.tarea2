package parallel

import (
	"fmt"
	"log"
	"sync"

	tridiagonal "../../lib"
	"gonum.org/v1/gonum/mat"
)

//Returns pointer to a vector of size @m and all values equal to @value
func AllEqualVec(m int, value float64) *mat.VecDense {
	v_values := make([]float64, m)
	for i := 0; i < m; i++ {
		v_values[i] = value
	}
	v := mat.NewVecDense(m, v_values)
	return v
}

//Formats matrix to print it in console
func MatPrint(X mat.Matrix) {
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%v\n", fa)
}

//Computes X_kp1 using a parallel aproach
//Parameters:
//@processes : amount of goroutines to be executed
//@A : coefficient matrix
//@b : constant terms vector
//@x_k : last solution vector
//@m : size of matrix
//Output: next solution vector
func getX_kp1(processes int, A *mat.Dense, b *mat.VecDense, x_k *mat.VecDense, m int) *mat.VecDense {

	x_kp1 := AllEqualVec(m, 0) //intialize x_(k+1) with zero values
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

	//Sets 1/p parts of the total of values in x_i, where p is the number of processors
	//This function will be called concurrently for parallel processing
	getSubX_n := func(initial int, wg *sync.WaitGroup) {
		defer wg.Done()
		var x_i float64
		for i := initial; i < m; i += processes {
			x_i = getX_i(i, x_k)
			x_kp1.SetVec(i, x_i)
		}
	}

	var wg sync.WaitGroup
	for i := 0; i < processes; i++ {
		wg.Add(1)            //add new goroutine to wait-group
		go getSubX_n(i, &wg) //start goroutine
	}

	wg.Wait() //wait until all goroutines finish

	return x_kp1

}

//Parallel implementation of Jacobi's method
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
func ParallelJacobi(A *mat.Dense, b *mat.VecDense, x_0 *mat.VecDense, iterMax int, tol float64, processes int) (*mat.VecDense, int, float64) {

	m, n := A.Caps() //get number of rows and collumns of the matrix

	if m == n {
		if tridiagonal.IsDiagDominant(A) {
			x_k := x_0 //initialize solution vector

			var err float64
			sol := mat.NewDense(m, 1, nil)
			for k := 0; k < iterMax; k++ {

				x_kp1 := getX_kp1(processes, A, b, x_k, m) //computes x_kp1 with parallel processing
				*x_k = *x_kp1
				sol.Product(A, x_k)
				sol.Sub(sol, b)
				err = sol.Norm(2)

				if err < tol { //stopping condition
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
