import sys
import numpy as np

TOL = sys.float_info.epsilon
overflow = sys.float_info.max

def myFactorial(n):
    factorial = 1
    for i in range(1, n+1):
        factorial *= i
    
    return factorial


def Cos(x):
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

def Sin(x):
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
 

def E(x):
    e = 0
    newTerm = 1
    n = 1
    while abs(newTerm) > TOL:
        e += newTerm
        newTerm = (x**n)/(myFactorial(n))
        n += 1

    return e

def E2(x):
    e = 0
    newTerm = 1
    n = 1
    while abs(newTerm) > TOL:
        e += newTerm
        newTerm = (x**(2*n))/(myFactorial(n))
        n += 1

    return e


def Cosh(x):
    return (1/2) * (E(x) + (1/E(x)))


def ln(x):
    ln = 0
    newTerm = 0
    n = 1

    if x == 0:
        pass

    else:
        while True:
            newTerm = (1/n) * (((x - 1)/(x + 1))**n)
            ln += newTerm
            n += 2

            if abs(newTerm) < TOL:
                break
        
    return 2*ln

def sqrt(x: float, tol: float = None) -> float:
    """
    Calcula la raíz cuadrada de un número positivo usando el método de Newton-Raphson.

    Args:
        x (float): Número del que se quiere calcular la raíz cuadrada.
        tol (float, opcional): Tolerancia de convergencia. Si no se da, se usa el épsilon de máquina.

    Returns:
        float: Aproximación de sqrt(x).

    Raises:
        TypeError: Si x o tol no son del tipo correcto.
        ValueError: Si x es negativo.
        RuntimeError: Si no converge.
    """
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

    except (TypeError, ValueError, RuntimeError) as e:
        ###print("Error en mi_raiz_cuadrada:", e)
        return None
    ##finally:
        ###print(f"Se intentó calcular sqrt({x}) con tolerancia {tol}.")

