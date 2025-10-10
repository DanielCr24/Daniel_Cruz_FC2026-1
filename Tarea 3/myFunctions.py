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

