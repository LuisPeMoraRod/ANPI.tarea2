/** Estudientes Grupo 2:
 *  Fabián Crawford Barquero
 *  Irene Muñoz Castro
 *  Luis Morales Rodríguez
 *  Steven Badilla Soto
 *  Adrián Trejos Salazar
 * 
 * 
 *                  Parte 3 Tarea 2 pregunta 1
 * Se quiere crear un metodo iterativo en C++ que calcule la aproximacion
 * de la pseudoinversa de una matriz A ∈ R^m×n. Este metodo es la generalizacion de la 
 * ecuacion de Newton-Schulz, que se presenta en el articulo "A family of 
 * iterative methods for computing the approximate inverse of a square matrix and inner
 * inverse of a non-square matrix desarrollado" por W. Li y Z. Li.
 * 
 * Se utilizara:
 * Una tolerancia de 10^−5
 * y un valor de p de 7
 * 
 */


#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

double tolerancia = 0.00001;
int p = 7;


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
    //Se definen las constantes para las condiciones iniciales de la iteración
    double alpha_0 = 1e-6; //Constante que acompaña el valor de x0
    mat x0 = (alpha_0)*A.t(); //Se define el valor inicial para x0
    mat xk = x0;
    //Generar matriz identidad y valorenes extra necesarios
    int n = A.n_cols;
    mat I = eye(n,n);
    mat multiplicador;
    double numerador;
    double denominador;
    double division;
    mat sumatoria = zeros(A.n_rows,A.n_rows);
    mat xk_n;
    double error;
    double p_Factorial = factorial(p);
    int iteraciones;
    
    // Inicia la iteraciones
    for (int q = 1; q < p; q++){
        multiplicador = pow( (A * xk), (q-1) );
        numerador = pow(-1, (q-1) );
        denominador = factorial(q) * factorial((p - q)  ); 
        division = numerador / denominador;
        sumatoria = sumatoria + ((division)* multiplicador);
        xk_n = xk * sumatoria;
        error = norm(xk_n*A - I,"for");       // Se calcula el error
        iteraciones = q;
        if (error < tolerancia){        // Se revisa si se cumple la condicion de parada
            break;
        }
        xk = xk_n; // Se actualiza xk
    }
    cout << "Numero de iteraciones: "<< iteraciones << endl;
    cout << "Error: "<< error << endl;
    cout << "xk: "<< xk << " //"<< endl;
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
            A.row(i).col(j).fill((pow(i+1,2)+ pow(j+1,2)));
        }
    }
    return A;
}

int main(){
    cout << "    --- Aproximacion de pseudoinversa:    " << endl;
    cout << "- Se genera la matriz" << endl;
    mat A = generar_Matriz();
    cout << "- Se hace el calculo: " << endl;
    pseudoinversa(A);
    return 0;
}