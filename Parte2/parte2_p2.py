import sympy as sy
import numpy as ny

def newton_raphson(x0,f,x,tol,iterMax):
    xk = x0
    fx = ny.zeros(len(f))
    k = 0
    error = tol + 1
    j = 0
    i = 0
    while k < iterMax and error > tol:
        while i < len(f):
            faux = sy.sympify(f[i])
            while j < len(x0):
                faux = sy.N(faux.subs(x[j],xk[j]))
                j += 1
            fx[i] = faux
            j = 0
            i += 1
        jac = jacobiano(f,x,xk,len(f))
        y = ny.linalg.solve(jac,fx)
        xk = xk - y
        error = ny.linalg.norm(fx,2)
        i = 0
        k += 1
    print("xk:",xk)
    print("k:",k)
    print("error:",error)
    return [xk,k,error]

def jacobiano(f,x,x0,m):
    jac = ny.zeros((m,m))
    i = 0
    j = 0
    while i < m:
        while j < m:
            xaux = sy.symbols(x[j])
            df = sy.diff(f[i],xaux)
            jac[i,j] = sy.N(df.subs(x[j],x0[j]))
            j += 1
        j = 0
        i += 1
    return jac

f = ["x**2+y**2+z**2-1","2*x**2+y**2-4*z","3*x**2-4*y+z**2"]
x = ["x","y","z"]
x0 = [0.5,0.5,0.5]


print(newton_raphson(x0,f,x,1e-10,500))