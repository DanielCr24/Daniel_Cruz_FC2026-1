#Tarea 1. Física Computacional#
#Ángel Daniel Cruz Flores#

#Problema 1.#

#a. Al ser distancias y tiempos muy grandes es necesario usar Unidades Astronómicas y años para evitar trabajar con números muy grandes en la simulación.

#b. 

distM = float(input("Ingrese la distancia en metros: "))
velM_S = float(input("Ingrese la velocidad en m/s: "))

def convertDist(d: float) -> float:
    '''
    Recibe la distancia en metros y la convierte a unidades astronómicas UA.

    Args:
        d (float): Distancia en metros a convertir a UA.
    
    Returns:
        (float): Distancia en UA.
    '''

    if not isinstance(d, float):
        raise TypeError(f"Se esperaba un float, pero se recibió {type(d).__name__}")

    return d/1.496e11

def convertVel(v: float) -> float:
    '''
    Recibe la velocidad en metros sobre segundos y la convierte a unidades astronómicas UA por año.

    Args:
        v (float): Velocidad en metros sobre segundos a convertir a UA por año.
    
    Returns:
        (float): Velocidad en UA por años.
    '''

    if not isinstance(v, float):
        raise TypeError(f"Se esperaba un float, pero se recibió {type(v).__name__}")
    
    return (v*(365.25 * 24 * 3600)) / 1.496e11

#b.

distAfelio, velAfelio = convertDist(5.28e12), convertVel(9.12e2)
print(distAfelio, velAfelio)

#c.



#Problema 2.#

#a.

from misFunciones import isPrime

def nPrimes(n: int) -> list:
    '''
    Recibe un entero positivo n y genera los primeros n primos.

    Args:
        n (int): Cantidad de primos a generar.

    Returns:
        list: Lista con los primeros n primos.
    '''

    if not isinstance(n, int):
        raise TypeError(f"Se esoeraba un int, pero se recibió {type(n).__name__}")
    
    primes = []
    i = 0

    while len(primes) < n:
        if isPrime(i) == True:
            primes.append(i)
        i += 1
    
    return primes

#b.

def twinPrimes(n: int) -> list:
    '''
    Recibe un entero positivo n y genera una lista de tuplas de todos los pares primos iguales o menores a n.

    Args:
        n (int): Primos menores o iguales a n.

    Returns:
        list: Lista de tuplas de todos los pares de primos gemelos menores o iguales a n.
    '''

    if not isinstance(n, int):
        raise TypeError(f"Se esperaba un int, pero se recibió {type(n).__name__}")
    
    twinPrimes = []
    primes = []
    i = 0
    j = 1
    
    while i < n:
        if isPrime(i) == True:
            primes.append(i)
        i += 1

    while j < n:
        if primes[j-1] - primes[j] == 2:
            twinPrimes.append((primes[j], primes[j-1]))
        j += 1

    return twinPrimes



#Problema 4.#

def myAverage(x: list) -> float:
    '''
    Recibe una lista finita y genera la media.

    Args:
        x (list): Secuencia finita de datos.

    Returns:
        average (float): Media de la secuencia de datos.
    '''

    if not isinstance(x, int):
        raise TypeError(f"Se esperaba una lista, pero se recibió {type(x).__name__}")
    
    average = 0

    for i in x:
        average += i
    
    average *= 1/len(x)

    return average

def myDesvest2s(x: list) -> float:
    '''
    Recibe una lista finita y genera la desviación estandar de dos pasos.

    Args:
        x (list): Secuencia finita de datos.

    Returns:
        desvest (float): Desviación estándar de dos pasos.
    '''

    if not isinstance(x, int):
        raise TypeError(f"Se esperaba una lista, pero se recibió {type(x).__name__}")

    sum = 0

    for i in x:
        sum += (i - myAverage(x))**2
    
    desvest = ((1/len(x)-1)*(sum))**(1/2)

    return desvest

def myDesvest1s(x: list) -> float:
    '''
    Recibe una lista finita y genera la desviación estandar de un paso.

    Args:
        x (list): Secuencia finita de datos.

    Returns:
        average (float): Desviación estandat de un paso.
    '''

    if not isinstance(x, int):
        raise TypeError(f"Se esperaba una lista, pero se recibió {type(x).__name__}")
    
    sum = 0

    for i in x:
        sum += i**2

    desvest = (1/(len(x)-1))*(sum - (len(x)*(myAverage(x))**2))

    return desvest
