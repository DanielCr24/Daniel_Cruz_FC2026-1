#Tarea 2. Física Computacional#
#Ángel Daniel Cruz Flores#

#Problema 1.#

#a. El tamaño en memoria de la grilla inicial es (2^2)*(2^3)*8 bytes = 256 bytes.
#   La memoria RAM =  16 GB, por cada iteración, se duplican las columnas o filas, 
#               2*..*2(nm)*8 = 2^k(nm) = 2^k(2^5)*8 = 2^(k+5)*8 = 1.6x10^{10} bytes
#   despejando k de la ecuación, 
#                                   2^{k} = (1.6x10^{10})/(2^{8})
#                                   2^{k} = 6.25x10^{7}
#                                       k = 25.897
#   por lo tanto, se pueden hacer aproximadamente 25 iteraciones antes de llenar la memoria RAM.

#b. Si el programa debe almacenar la grilla actual y la anterior, entonces, 
#                  2^{8}(2^{k} + 2^{k-1}) = 1.6x10^{10}
#                                2^{2k-1} = (1x10^{10})/(2^{8})
#                                2^{2k-1} = 6.25x10^{7}
#                                  2k - 1 = 25.897
#                                      2k = 26.897
#                                       k = 13.449
#   por lo tanto, si el programa guarda la grilla actual y la de la iteración anterior solo se pueden realizar 13 iteraciones.

#c. Si el sistema operativo y las librerias ocupan 2GB de memoria RAM y la computadora reserva el 10% como buffer, entonces,
#            M_D = (16 - 2 - 1.6) GB = 12.4 GB = 1.24x10^{10} bytes
#   Caso 1. Solo se guarda la grilla actual.
#                                   2^{k} = (1.24x10^{10})/(2^{8}) 
#                                   2^{k} = 48,437,500
#                                       k = 25.53
#   Caso 2. Se almacena la grilla actual y la anterior.
#                                2^{2k-1} = 48,437,500  
#                                  2k - 1 = 25.53
#                                       k = 13.265
# Notamos que se pueden realizar los mismos pasos que usando la totalidad de la memoria RAM.

#d. Escribimos una función que simule el crecimiento de la grilla.
def simGrid(maxMemory, keepPrevious = False):
    n = 0
    currentMemory = 2**8

    while currentMemory < maxMemory:
        n += 1
        newMemory = currentMemory*2 

        if keepPrevious:
            currentMemory = newMemory + currentMemory
        else:
            currentMemory = newMemory
    
    currentMemory /= 2

    currentMemoryB = currentMemory
    currentMemoryMB = currentMemory*(1e-6)
    currentMemoryGB = currentMemory*(1e-9)
    
    print(f"El límite de memoria se alcanzó en la iteración {n-1} y el tamaño de la grilla fue {currentMemoryGB} GB.")
    return currentMemoryB, currentMemoryMB, currentMemoryGB


#Problema 2.#

import matplotlib.pyplot as plt
import numpy as np

#Definimos una función que evalue un polinomio de grado n usando la regla de Horner#
def HornerR(a: list, x: int) -> float:
    result = 0
    for i in reversed(a):
        result = (x*result) + i
    
    return result

#Definimos una lista con los coeficientes del polinomio
y = [-0.981, 35.355, 19.62, -4.905]

#Definimos una lista con los valores del tiempo t a evaluar en el polinomio y(t)
t = [0.5, 1.0, 1.5, 2.0, 2.5]

#Creamos un ciclo for para calcular los valores del polinomio evaluado en cada t
y_t = []
for i in t:
    y_t.append(HornerR(y, i))

#Graficamos la trayectoria completa del proyectil

x = np.linspace(0, 5.4, 100)
h = HornerR(y, x)

#plt.plot(x, h)

plt.title('Gráfica de la función y = x^2 + 2x + 1')
plt.xlabel('Eje X')
plt.ylabel('Eje Y')

#plt.show()

#Encontramos la derivada
def rootAprox(coeficientesD, intervalo, tolerancia=0.0001):
    """
    Encuentra una raíz aproximada de un polinomio usando el método de bisección.
    """
    a, b = intervalo
    fa = HornerR(coeficientesD, a)
    fb = HornerR(coeficientesD, b)
    
    # Verificar que hay cambio de signo
    if fa * fb > 0:
        return None
    
    # Método de bisección
    while (b - a) > tolerancia:
        medio = (a + b) / 2
        f_medio = HornerR(coeficientesD, medio)
        
        if abs(f_medio) < tolerancia:
            return medio
        
        if fa * f_medio < 0:
            b = medio
            fb = f_medio
        else:
            a = medio
            fa = f_medio
    
    return (a + b) / 2

# Derivada del polinomio: y'(t) = -19.62t³ + 58.86t² + 70.71t - 0.981
coeficientesD = [-0.981, 70.71, 58.86, -19.62]

# Encontrar raíz en el intervalo [0, 5] 
tiempo_maximo = rootAprox(coeficientesD, [0, 5])

#Calculand0 la altura máxima
if tiempo_maximo is not None:
    altura_maxima = HornerR(y, tiempo_maximo)
    print(f"\nTiempo de altura máxima: {tiempo_maximo:.3f} s")
    print(f"Altura máxima: {altura_maxima:.3f} m")
else:
    print("No se pudo encontrar el tiempo de altura máxima")

#Problema 3.#

import sys

def machine_epsilon():
    """Calcula la precisión de la máquina (épsilon)"""
    epsilon = 1.0
    while 1.0 + epsilon > 1.0:
        epsilon /= 2.0
    return epsilon * 2.0

EPSILON = machine_epsilon()

def sin_taylor(x, tol=EPSILON):
    """
    Calcula sin(x) usando serie de Taylor
    """
    result = 0.0
    term = x
    n = 1
    
    while abs(term) > tol:
        result += term
        n += 2
        term = (-term * x * x) / ((n - 1) * n)
    
    return result

def cos_taylor(x, tol=EPSILON):
    """
    Calcula cos(x) usando serie de Taylor
    """
    result = 0.0
    term = 1.0
    n = 0
    
    while abs(term) > tol:
        result += term
        n += 2
        term = (-term * x * x) / ((n - 1) * n)
    
    return result

def tan_taylor(x, tol=EPSILON):
    """
    Calcula tan(x) = sin(x)/cos(x)
    """
    cos_val = cos_taylor(x, tol)
    if abs(cos_val) < tol:
        raise ValueError("tan(x) no está definido para este valor de x")
    
    return sin_taylor(x, tol) / cos_val

def log_taylor(x, tol=EPSILON):
    """
    Calcula log(x) usando serie de Taylor centrada en 1
    """
    if x <= 0:
        raise ValueError("log(x) solo está definido para x > 0")
    
    # Si x > 2, usar propiedad log(x) = -log(1/x)
    if x > 2:
        return -log_taylor(1/x, tol)
    
    # Si x < 0.5, usar propiedad log(x) = -log(1/x)
    if x < 0.5:
        return -log_taylor(1/x, tol)
    
    # Transformar para usar serie centrada en 1
    u = x - 1
    if abs(u) >= 1:
        raise ValueError("No se puede calcular con la serie centrada en 1")
    
    result = 0.0
    term = u
    n = 1
    
    while abs(term) > tol:
        result += term
        n += 1
        term = (-term * u * (n - 1)) / n
    
    return result


#Problema 4.#

def eNorm(x: list) -> float:
    norm = 0
    
    for i in x:
        norm += i**2

    norm **= 1/2

    return norm

#Problema 5.#

import decimal
from decimal import Decimal, getcontext

class Interval:
    def __init__(self, a, b=None):

        if b is None:
            b = a
        
        # Configurar el contexto decimal para redondeo
        getcontext().prec = 28  # Precisión suficiente
        
        # Redondear a hacia abajo y b hacia arriba
        self.a = Decimal(a).to_integral_value(rounding=decimal.ROUND_FLOOR)
        self.b = Decimal(b).to_integral_value(rounding=decimal.ROUND_CEILING)
        
        if self.a > self.b:
            raise ValueError("a debe ser menor o igual que b")
    
    def __repr__(self):
        return f"{float(self.a)}, {float(self.b)}"
    
    def __str__(self):
        return f"[{float(self.a)}, {float(self.b)}]"
    
    def __add__(self, other):
        """Suma de intervalos: [a + c, b + d]"""
        if not isinstance(other, Interval):
            other = Interval(other)
        
        # Redondear extremos izquierdo hacia abajo y derecho hacia arriba
        left = (self.a + other.a).to_integral_value(rounding=decimal.ROUND_FLOOR)
        right = (self.b + other.b).to_integral_value(rounding=decimal.ROUND_CEILING)
        
        return Interval(left, right)
    
    def __sub__(self, other):
        """Resta de intervalos: [a - d, b - c]"""
        if not isinstance(other, Interval):
            other = Interval(other)
        
        left = (self.a - other.b).to_integral_value(rounding=decimal.ROUND_FLOOR)
        right = (self.b - other.a).to_integral_value(rounding=decimal.ROUND_CEILING)
        
        return Interval(left, right)
    
    def __mul__(self, other):
        """Multiplicación de intervalos: [min(ac, ad, bc, bd), max(ac, ad, bc, bd)]"""
        if not isinstance(other, Interval):
            other = Interval(other)
        
        # Calcular todos los productos posibles
        products = [
            self.a * other.a,
            self.a * other.b,
            self.b * other.a,
            self.b * other.b
        ]
        
        # Encontrar mínimo y máximo con redondeo
        min_val = min(products).to_integral_value(rounding=decimal.ROUND_FLOOR)
        max_val = max(products).to_integral_value(rounding=decimal.ROUND_CEILING)
        
        return Interval(min_val, max_val)
    
    def __truediv__(self, other):
        """División de intervalos: [min(a/c, a/d, b/c, b/d), max(a/c, a/d, b/c, b/d)]"""
        if not isinstance(other, Interval):
            other = Interval(other)
        
        # Verificar que 0 no está en el denominador
        if other.a <= 0 <= other.b:
            raise ValueError("División por intervalo que contiene cero")
        
        # Calcular todos los cocientes posibles
        quotients = [
            self.a / other.a,
            self.a / other.b,
            self.b / other.a,
            self.b / other.b
        ]
        
        # Encontrar mínimo y máximo con redondeo apropiado
        min_val = min(quotients).to_integral_value(rounding=decimal.ROUND_FLOOR)
        max_val = max(quotients).to_integral_value(rounding=decimal.ROUND_CEILING)
        
        return Interval(min_val, max_val)
    
    def __radd__(self, other):
        return self + other
    
    def __rsub__(self, other):
        return Interval(other) - self
    
    def __rmul__(self, other):
        return self * other
    
    def __rtruediv__(self, other):
        return Interval(other) / self

# Pruebas de la implementación
if __name__ == "__main__":
    print("=== Pruebas de la aritmética de intervalos ===\n")
    
    # Prueba 1: Sumar diez veces el número 0.1
    print("1. Sumar diez veces 0.1:")
    cero_punto_uno = Interval(0.1)
    resultado = Interval(0.0)
    
    for _ in range(10):
        resultado += cero_punto_uno
    
    print(f"   0.1 + 0.1 + ... + 0.1 (10 veces) = {resultado}")
    print(f"   Valor exacto debería ser: 1.0")
    print(f"   ¿Contiene el valor exacto? {1.0 >= float(resultado.a) and 1.0 <= float(resultado.b)}")
    print()
    
    # Prueba 2: Calcular el cociente 1/3
    print("2. Calcular 1/3:")
    uno = Interval(1.0)
    tres = Interval(3.0)
    division = uno / tres
    
    print(f"   1/3 = {division}")
    print(f"   Valor exacto: 0.3333333333333333...")
    print(f"   ¿Contiene el valor exacto? {1/3 >= float(division.a) and 1/3 <= float(division.b)}")
    print()
    
    # Prueba 3: Operaciones básicas
    print("3. Otras operaciones básicas:")
    x = Interval(1, 2)
    y = Interval(3, 4)
    
    print(f"   x = {x}")
    print(f"   y = {y}")
    print(f"   x + y = {x + y}")
    print(f"   x - y = {x - y}")
    print(f"   x * y = {x * y}")
    print(f"   x / y = {x / y}")
    print()
    
    # Prueba 4: Mostrar el efecto del redondeo
    print("4. Efecto del redondeo controlado:")
    # Sin control de redondeo
    a_simple = 0.1
    suma_simple = sum([a_simple for _ in range(10)])
    
    # Con aritmética de intervalos
    a_interval = Interval(0.1)
    suma_interval = Interval(0.0)
    for _ in range(10):
        suma_interval += a_interval
    
    print(f"   Suma simple (sin control): {suma_simple}")
    print(f"   Suma con intervalos: {suma_interval}")
    print(f"   Error absoluto simple: {abs(suma_simple - 1.0)}")
    print(f"   El intervalo garantiza contener el valor verdadero: {1.0 >= float(suma_interval.a) and 1.0 <= float(suma_interval.b)}")


#Problema 6.#

#Definimos una función que recorra elemento por elemento hasta encontrar el valor buscado.
def myFindLin(x: list, n: int):
    count = 0
    for i in range(len(x)):
        count += 1
        if x[i] == n:
            return i
            break

    return -1

#a. En el peor de los casos realiza len(x) operaciones (recorre cada elemento de la lista), que es cuando no existe el valor buscado en la lista dada.
#b. O(n): En el peor de los casos, realiza n comparaciones.
#   Ω: En el mejor caso, encuentra el elemento en la primera posición.
#   Θ: En el caso promedio, realiza aproximadamente n/2 comparaciones.     
        
def myFindBin(x: list, n: int):
    x.sort()
    
    i = 0
    left = 0
    right = len(x) - 1

    while left <= right:
        i += 1
        mid = (right - left)//2

        if x[mid] == n:
            return mid
        elif x[mid] < n:
            left = mid + 1
        else:
            right = mid - 1
    
    return -1

#a. Notamos que la lista se ordena, en cada paso, se divide el espacio de busqueda a la mitad,
#   por lo tanto, en el pero de los casos, el número máximo de divisiones es log_2(n).
#b. O(log n): En el peor de los casos crece logarítmicamente con n.
#   Ω(1): En el mejor de los casos, el elemento buscado está en el medio, solo se hace una comparación.
#   Θ(log n): En el caso promedio, también es logarítmico.

#Experimentos computacionales

#Importamos la libreria random para crear números aleatorios y la libreria time para medir tiempos.
import random
import time
#Genera listas de diferentes tamaños
size = [10**3, 10**4, 10**5, 10**6]
linealTime = []
binaryTime = []
ordenTime = []

for i in size:
    #Generamos lista aleatoria
    list = random.sample(range(i*2), i)
    n = -1

    #Medimos tiempo de búsqueda lineal
    startLineal = time.time()
    myFindLin(list, n)
    finishLineal = time.time()
    linealTime.append(finishLineal)

    #Medimos tiempo de busqueda binaria incluyendo el ordenamiento
    startOrder = time.time()
    listOrder = sorted(list)
    finishOrder = time.time() - startOrder
    ordenTime.append(finishOrder)

    startBin = time.time()
    myFindBin(listOrder, n)
    finishBin = time.time() - startBin
    binaryTime.append(finishBin)

    print(f"  Lineal: {finishLineal:.6f}s")
    print(f"  Binaria (búsqueda + ordenamiento): {finishOrder + finishBin:.6f}s")
    print(f"  - Ordenamiento: {finishOrder:.6f}s")
    print(f"  - Búsqueda: {finishBin:.6f}s")

#Graficamos los resultados
#Tiempos completos
plt.subplot(2, 1, 1)
plt.plot(size, linealTime, 'o-', label='Búsqueda Lineal')
plt.plot(size, binaryTime, 's-', label='Búsqueda Binaria (con ordenamiento)')
plt.plot(size, ordenTime, '^--', label='Solo ordenamiento')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Tamaño de la lista (n)')
plt.ylabel('Tiempo de ejecución (s)')
plt.title('Comparación de algoritmos de búsqueda (listas aleatorias)')
plt.legend()
plt.grid(True)

#Solo búsqueda sin ordenamientos
plt.subplot(2, 1, 2)
plt.plot(size, linealTime, 'o-', label='Búsqueda Lineal')
plt.plot(size, [binaryTime[k] - ordenTime[k] for k in range(len(size))], 
            's-', label='Búsqueda Binaria (solo búsqueda)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Tamaño de la lista (n)')
plt.ylabel('Tiempo de ejecución (s)')
plt.title('Comparación de la operación de búsqueda solamente')
plt.legend()
plt.grid(True)

plt.tight_layout()
#plt.show()


