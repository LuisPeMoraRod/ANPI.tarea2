import sympy as sy
import numpy as ny
import matplotlib.pyplot as plt


##Metodo para calcular la solucion de un sistema de ecuaciones no lineales mediante newthon-raphson
#Entradas: x0 -> vector inicial, f -> vector de funciones, x -> verctor de variables, tol-> tolerancia, iterMax -> iteraciones maximas
#Salidas: xk -> aproximacion de la solucion, k -> numero de iteraciones, error -> error de la aproximacion
def newton_raphson(x0,f,x,tol,iterMax):
    xk = x0
    fx = ny.zeros(len(f)) #vector para almacenar los resultados de las funcionas evaluadas en xk
    k = 0
    error = tol + 1
    eaux = []
    kaux = []
    j = 0
    i = 0
    #Se ejecuta hasta tener un error < que la tolerancia o superar las iteraciones maximas
    while k < iterMax and error > tol:
        #Iteraciones para evaluar cada funcion del vector f
        while i < len(f):
            faux = sy.sympify(f[i]) #Almacena la funcion de f actual para poder ir sustituyendo cada variable
            #Iteraciones para sustituir cada variable en el vector x por su valor actual para xk en la funcion
            while j < len(x0):
                faux = sy.N(faux.subs(x[j],xk[j]))
                j += 1
            fx[i] = faux
            j = 0
            i += 1
        jac = jacobiano(f,x,xk,len(f)) #Se obtiene el jacobiano evaluado en xk
        y = ny.linalg.solve(jac,fx) #Se resuelve Jf(xk)y=f(xk) para evitar el calculo de la inversa
        xk = xk - y #Se aplica la formula para obtener xk+1 mediante el metodo newton-raphson
        error = ny.linalg.norm(fx,2) #Se calcula el error como ||f(xk)||
        eaux.append(error)
        kaux.append(k)
        i = 0
        k += 1

    #Grafica error vs iteraciones
    plt.plot(kaux,eaux) #Se grafican las listas con los errores para cada iteracion
    plt.scatter(kaux, eaux)
    plt.xlabel("Iteraciones")
    plt.ylabel("Error")
    plt.title("Error vs Iteraciones")
    plt.show()

    return [xk,k,error]

#Metodo para calcular el Jacobiano de un vector de funciones evaluado en los valores dados
#Entradas: x0 -> vector con los valores a evaluar, f -> vector de funciones, x -> verctor de variables, m -> tamaño dado
#Salidas: jac -> l la matriz jacobiana evaluada en los valores dados
def jacobiano(f,x,x0,m):
    jac = ny.zeros((m,m)) #Matriz de ceros para almacenar los resultados
    i = 0
    j = 0
    n = 0
    #Iteraciones para calcular cada funcion de f
    while i < m:
        #Iteraciones para calcular la derivada parcial para cada variable
        while j < m:
            xaux = sy.symbols(x[j])#Se crea una variable auxiliar para almacenar como simbolo la variable que se va a evaluar de el vector x
            df = sy.diff(f[i],xaux)#Se calcula la derivada parcial de la funcion actual con la variable actual
            while n < m: #Se evaluan todas las variables en la derivada parcial
                df = sy.N(df.subs(x[n],x0[n]))
                n += 1
            jac[i,j] = df #Se almacena el valor evaluado
            j += 1
            n = 0
        j = 0
        i += 1
    return jac







