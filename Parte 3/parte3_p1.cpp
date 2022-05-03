/** Estudientes Grupo 2:
 *  Fabián Crawford Barquero
 *  Irene Muñoz Castro
 *  Luis Morales Rodríguez
 *  Steven Badilla Soto
 *  Adrián Trejos Salazar
 * 
 * 
 *                  Parte 3 Tarea 2
 * Se quiere crear un metodo iterativo en C++ que calcule la aproximacion
 * de la pseudoinversa de una matriz A ∈ R^m×n. Este metodo es la generalizacion de la 
 * ecuacion de Newton-Schulz, que se presenta en el articulo "A family of 
 * iterative methods for computing the approximate inverse of a square matrix and inner
 * inverse of a non-square matrix desarrollado" por W. Li y Z. Li.
 * 
 * Se utilizara:
 * Una tolerancia de 10^−5
 * Iteraciones maximas 500
 * un valor de q de 
 * y un valor de p de
 * 
 */


#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

double tolerancia = 0.00001;
int iter_max = 500;
int q = 1;
int p = 1;


 /** Funcion Factorial, Esta funcion retorna el numero de entrada en factorial
    * @param n  es el numero natural positivo al que se le aplica factorial
    * @return resultado de la operacion
    */
double factorial(double n){
    if ((n==0)||(n==1)){  //si el numero ingresado es cero o uno
       return 1;         //retorna 1 
    }else{  
       return n*factorial(n-1); // retorna el factorial del numero menos 1
    }
}

/** Metodo iterativo para aproximar la pseudoinversa de una Matriz.
  * @param A una matriz A ∈ R^m×n
  * @return resultado de la aproximacion de la pseudoinversa
  */
mat pseudoinversa(mat A){
    // Calculos necesarios para obtener normsq (||A||2)
    mat mat_aux= A * A.t();
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, mat_aux);
    int normsq= max(eigval); // Calculo de normsq
    
    // Calcular el valor inicial
    mat Xk = ((1+0.0)/normsq) * A.t(); 
    mat b = ones(A.n_rows);
    mat xk = Xk*b;
  
    //Generar matriz identidad y valorenes extra necesarios
    int n = A.n_cols;
    mat I = eye(n,n);
    mat multiplicador;
    mat numerador;
    mat denominador;
    mat xk_n;
    double error;

    // Inicia la iteraciones
    for (int k=1; k<iter_max; k++){
        multiplicador = pow( (A * xk), (q-1) );
        numerador = pow(-1, (q-1) ) * factorial(p);
        denominador = factorial(q) * factorial(p - q);

        xk_n = xk * ((numerador / denominador) * multiplicador);

        error = norm(xk_n*A - I,"fro"); // Se calcula el error
        if (error < tolerancia){        // Se revisa si se cumple la condicion de parada
            cout << "Numero de iteraciones: "<< k << endl;
            break;
        }
        xk = xk_n; // Se actualiza xk
    }
    cout << "Error: "<< error << endl;
    cout << "xk: "<< xk << endl;
    return xk;
}

/** Funcion que genera una matriz A ∈ R^45×30
 * donde A i,j = i 2 + j 2 , para i = 1, 2, ..., 45 y j = 1, 2, ..., 30
 * 
 */
mat generar_Matriz(){
    mat A = zeros(45,30);
    for (int i=0; i<A.n_rows; i++){
        for (int j=0; j<A.n_cols; j++){
            A(i,j) = pow((i+1),2) + pow((j+1),2);
        }
    }
    return A;
}



//g++ parte3_p1.cpp -o parte3_p1 -O2 -larmadillo -llapack -lblas

int main(){
    cout << "    --- Aproximacion de pseudoinversa:    " << endl;
    cout << "- Se genera la matriz" << endl;
    mat A = generar_Matriz();
    cout << "- Se hace el calculo: " << endl;
    pseudoinversa(A);
    cout << "//" << endl;
    return 0;
}