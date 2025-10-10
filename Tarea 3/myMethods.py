import sys
import myFunctions
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, c, m_p, m_n, eV

EPSILON = sys.float_info.epsilon

#Definimos el método de la Bisección
def mBis(f, a, b, maxIter = 1000, TOL = EPSILON):
    history = []
    fn = []

    try:

        if f(a) * f(b) > 0:
            raise ValueError(f"No se puede asegurar que la función tenga raíz en [{a}, {b}]")
        
        count = 0

        #Verifica si la pendiente de la recta que une ambos extremos es positivo o negativo
        #Si es positivo, entonces, si f(a) * f(m) < 0, debe enfocarse en la valores a su derecha para hallar el cero
        #Si es negativo, entonces, si f(a) * f(m) < 0, debe enfocarse en los valores a su izquierda para hallar el cero
        if f(a) < f(b):
            while abs(a - b) > TOL and count < maxIter:
                m = a + 0.5 * (a - b)
                history.append(m)
                fn.append(abs(f(m)))

                if f(a) * f(m) < 0:
                    a = m
                else:
                    b = m

                count += 1

            return m, count, history, fn
        
        else:
            while abs(a - b) > TOL and count < maxIter:
               m = a + 0.5*(b - a)
               history.append(m)
               fn.append(abs(f(m)))

               if f(a) * f(m) < 0:
                    b = m     
               else:
                    a = m
               
               count += 1
          
            return m, count, history, fn
        
    except Exception as e:
          print(f"[Error en Bisección] {e}")
          return None, 0, []
    
#Definimos el método de Newton-Raphson    
def mNR(f, df, x0, maxIter = 1000, TOL = EPSILON):
    history = []
    fn = []

    try:

        for i in range(maxIter):
            if df(x0) == 0:
                raise ZeroDivisionError("Se está dividiendo entre cero")
            
            x = x0 - f(x0) / df(x0)
            history.append(x)
            fn.append(abs(f(x)))

            if abs(x - x0) < TOL * max(1.0, abs(x)):
                return x, i+1, history, fn

            x0 = x
        
        print("Error: Newton-Raphson no llegó a EPSILON; devolviendo última aproximación")
        return x, maxIter, history, fn
    
    except Exception as e:
               print(f"[Error en Newton-Raphson] {e}")
               return None, 0, []
    
#Definimos el método de la Secante
def mSec(f, x1, x0, maxIter = 1000, TOL = EPSILON):
    history = []
    fn = []

    try:
    
        for i in range(maxIter):
            if f(x1) - f(x0) == 0:
                raise ZeroDivisionError("Se está dividiendo entre cero")
                pass
        
            x2 = x1 - (((x1 - x0) / (f(x1) - f(x0))) * f(x1))
            history.append(x2)
            fn.append(abs(f(x2)))

            if abs(x2 - x1) < TOL * max(1.0, abs(x2)):
                return x2, i+1, history, fn

            x0, x1 = x1, x2

        print("Error: Secante no llegó a EPSILON; devolviendo última aproximación")
        return x1, maxIter, history, fn
    
    except Exception as e:
        print(f"[Error en Secante] {e}")
        return None, 0, []
    







