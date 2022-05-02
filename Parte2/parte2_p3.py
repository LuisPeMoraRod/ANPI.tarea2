# -*- coding: utf-8 -*-
"""
Ejemplos numéricos de sección 4 del documento "ejemplos.pdf"
"""
from parte2_p2 import *
    
def p_newton_raph_a():
    print('\n')
    print('Ejemplo numérico a)\n')
    f = ["exp(x1**2)-exp(sqrt(2)*x1)","x1-x2"]
    x = ["x1","x2"]
    x0 = [2.3, 2.3]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    x0 = [1.8, 1.8]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    x0 = [0.8, 0.8]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    
def p_newton_raph_b():
    print('\n')
    print('Ejemplo numérico b)\n')
    f = ["x1+exp(x2)-cos(x2)","3*x1-x2-sin(x2)"]
    x = ["x1","x2"]
    x0 = [1.5, 2]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    x0 = [0.3, 0.5]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    
def p_newton_raph_c():
    print('\n')
    print('Ejemplo numérico c)\n')
    f = ["x1**2-2*x1-x2+0.5","x1**2+4*x2**2-4"]
    x = ["x1","x2"]
    x0 = [3, 2]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    x0 = [1.6, 0]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    
def p_newton_raph_d():
    print('\n')
    print('Ejemplo numérico d)\n')
    f = ["x1**2+x2**2-1","x1**2-x2**2+0.5"]
    x = ["x1","x2"]
    x0 = [0.7, 1.2]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    x0 = [-1, -2]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    
    """
def p_newton_raph_e():
    print('\n')
    print('Ejemplo numérico e)\n')
    f = ["sin(x1)+(x2*cos(x1))","x1-x2"]
    x = ["x1","x2"]
    x0 = [1.2, -1.5]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    x0 = [-0.6, 0.6]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    """
    """
def p_newton_raph_g():
    print('\n')
    print('Ejemplo numérico g)\n')
    f = ["(x2*x3)+(x4*(x2+x3))", "(x1*x3)+(x4*(x1+x3))", "(x1*x2)+(x4*(x1+x2))", "(x1*x2)+(x1*x3)+(x2*x3)-1"]
    x = ["x1", "x2", "x3", "x4"]
    x0 = [-1, -1, -1, -1]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
    x0 = [2, 2, 2, 0]
    print("Usando x0:", x0)
    sol = newton_raphson(x0,f,x,1e-10,500)
    print("Aproximacion de la solucion: xk =", sol[0])
    print("En", sol[1], "iteraciones")
    print("Con un error de", sol[2])
    print('\n')
"""


p_newton_raph_a()
p_newton_raph_b()
p_newton_raph_c()
p_newton_raph_d()
#p_newton_raph_e()
#p_newton_raph_g()
