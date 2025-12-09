import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.constants import hbar, c, m_p, m_n

# =============================================================================
# MÓDULO myFunctions.py ACTUALIZADO
# =============================================================================

TOL = sys.float_info.epsilon
overflow = sys.float_info.max

def myFactorial(n):
    """Calcula el factorial de n"""
    try:
        if n < 0:
            raise ValueError("Factorial no definido para números negativos")
        if n == 0:
            return 1
        factorial = 1
        for i in range(1, n+1):
            factorial *= i
        return factorial
    except Exception as e:
        print(f"Error en myFactorial: {e}")
        return None

def Cos(x):
    """Calcula cos(x) usando series de Taylor"""
    try:
        pi = np.pi
        if x > 2 * pi:
            x = 2 * pi * abs(round(x / (2 * pi)) - x / (2 * pi)) 
            cos = Cos(x)
        else:
            cos = 0
            newTerm = 1
            n = 1
            while True:
                cos += newTerm
                newTerm = (((-1)**n)/(myFactorial(2*n)))*(x**(2*n))
                n += 1
                if abs(newTerm) < TOL:
                    break
        return cos
    except Exception as e:
        print(f"Error en Cos: {e}")
        return None

def Sin(x):
    """Calcula sin(x) usando series de Taylor"""
    try:
        pi = np.pi
        if x > 2 * pi:
            x = -2 * pi * abs(round(x / (2 * pi)) - x / (2 * pi)) 
            sin = Sin(x)
        else:
            sin = 0
            newTerm = x
            n = 1
            while True:
                sin += newTerm
                newTerm = (((-1)**n)/myFactorial((2*n) + 1))*(x**((2*n) + 1))
                n += 1
                if abs(newTerm) < TOL:
                    break
        return sin
    except Exception as e:
        print(f"Error en Sin: {e}")
        return None

def E(x):
    """Calcula e^x usando series de Taylor"""
    try:
        e = 0
        newTerm = 1
        n = 1
        while abs(newTerm) > TOL:
            e += newTerm
            newTerm = (x**n)/(myFactorial(n))
            n += 1
        return e
    except Exception as e:
        print(f"Error en E: {e}")
        return None

def Cosh(x):
    """Calcula cosh(x)"""
    try:
        return (1/2) * (E(x) + (1/E(x)))
    except Exception as e:
        print(f"Error en Cosh: {e}")
        return None

def Sinh(x):
    """Calcula sinh(x)"""
    try:
        return (1/2) * (E(x) - (1/E(x)))
    except Exception as e:
        print(f"Error en Sinh: {e}")
        return None

def ln(x):
    """Calcula ln(x) usando series"""
    try:
        if x <= 0:
            raise ValueError("ln no definido para números no positivos")
        
        ln_val = 0
        newTerm = 0
        n = 1

        while True:
            newTerm = (1/n) * (((x - 1)/(x + 1))**n)
            ln_val += newTerm
            n += 2
            if abs(newTerm) < TOL:
                break
        
        return 2*ln_val
    except Exception as e:
        print(f"Error en ln: {e}")
        return None

def sqrt(x, tol=None):
    """Calcula raíz cuadrada usando Newton-Raphson"""
    try:
        if not isinstance(x, (int, float)):
            raise TypeError(f"x debe ser numérico, se recibió: {type(x).__name__}")
        if x < 0:
            raise ValueError("No se puede calcular la raíz cuadrada de un número negativo.")
        if tol is not None and not isinstance(tol, float):
            raise TypeError(f"tol debe ser numérico o None, se recibió: {type(tol).__name__}")
        if tol is None:
            tol = TOL

        # Casos triviales
        if x == 0:
            return 0.0
        if x == 1:
            return 1.0

        # Método de Newton-Raphson
        y = x / 2.0
        for _ in range(1000):
            y_new = 0.5 * (y + x / y)
            if abs(y_new - y) < tol * abs(y_new):
                return y_new
            y = y_new

        raise RuntimeError("No se alcanzó convergencia en 1000 iteraciones.")
    except Exception as e:
        print(f"Error en sqrt: {e}")
        return None

def cot(x):
    """Función cotangente"""
    try:
        return Cos(x) / Sin(x)
    except ZeroDivisionError:
        return float('inf')
    except Exception as e:
        print(f"Error en cot: {e}")
        return None


