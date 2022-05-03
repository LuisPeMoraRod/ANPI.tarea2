/** Estudientes Grupo 2:
 *  Fabián Crawford Barquero
 *  Irene Muñoz Castro
 *  Luis Morales Rodríguez
 *  Steven Badilla Soto
 *  Adrián Trejos Salazar
 * 
 * 
 *                  Parte 3 Tarea 2 Pregunta 2
 * 
 * Se utilizara:
 * Una tolerancia de 10^−5
 * y un valor de p de 8
 * 
 */


#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

double tolerancia = 0.00001;
int p = 8;


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
    double alpha_0 = 1.7e-3; //Constante que acompaña el valor de x0
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



int main(){
    cout << "    --- Aproximacion de pseudoinversa:    " << endl;
    cout << "- Se genera la matriz" << endl;
    mat A = {{1777,0},{1838,12},{1752,0},{1826,15},{1862,2},{1854,5},{1882,0},{1815,0},{1835,2},{1843,20}};
    cout << "- Se hace el calculo: " << endl;
    pseudoinversa(A);
    return 0;
}