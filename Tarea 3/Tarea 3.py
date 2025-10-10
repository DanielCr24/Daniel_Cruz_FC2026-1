import myFunctions
import myMethods
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import hbar, c, m_p, m_n, eV

EPSILON = sys.float_info.epsilon

#Tarea 3.#
#Ángel Daniel Cruz Flores

#Problema 1.#

#Definimos la ecuación de vibraciones de una viga en voladizo f(ω)
def ecVib(w): 
    return (myFunctions.Cos(w) * myFunctions.Cosh(w)) + 1

#Importamos el método de Bisección de la libreria myMethods
root, count, history, fn = myMethods.mBis(ecVib, 0, 4)

#a. 
#Imprimimos la menor raíz positiva de f(ω)
print(f"La menor raíz positiva de la ecuación de vibraciones es: {root} y la alcanza en {count} pasos.")

#b. 
#Graficamos f(ω) en ω ∈ [0, 4]

w = np.linspace(0, 4, 100)
f = []
for i in w:
    f.append(ecVib(i))

plt.plot(w, f, color = "orange")
plt.axhline(y = 0, linestyle = ":")
plt.axvline(0, linestyle = ":")
plt.axvline(4, linestyle = ":")
plt.plot(root, ecVib(root), 'ro', markersize = 6)
plt.xlabel("ω")
plt.ylabel("f(ω)")
plt.title("Vibraciones de una viga en voladizo")
plt.grid(True)
plt.show()

#c.
#Graficamos |f(ωn)| vs. iteración en escala semilog 
n = []
for i in range(1, count + 1):
    n.append(i)

plt.semilogy(n, fn)
plt.xlabel("n")
plt.ylabel("|f(ωn)|")
plt.title("Convergencia del método de la Bisección")
plt.grid(True)
plt.show()

#Notamos que el método converge de manera exponencial, sin embargo, no es uniforme debido a los picos que se observan.

#Problema 2.#

#Definimos el modelo de enfriamiento de Newton
def modelNewton(t, Ta, T0, k):
    return Ta + ((T0 - Ta) * myFunctions.E(-k * t))

#Definimos el modelo de enfriamiento de Newton para f(t) = 0
def Termic(t):
    return 20 + (90 - 20) * myFunctions.E(-0.07 * t) - 50
#Definimos la derivada de f(t)
def DTermic(t):
     return (-0.07)*(90 - 20) * myFunctions.E(0.07 * t)

#a.
#Importamos el método de Newton-Raphson de la libreria myMethods
root, count, history, fn = myMethods.mNR(Termic, DTermic, 10)
#Imprimimos t tal que T(t) = 50◦C empleando el método de Newton-Raphson usando t0 = 10
print(f"El valor de t tal que T(t) = 50◦C empleando el método de Newton-Rapshon es: {root} y lo alcanza en {count} pasos.")

#b.
#Graficamos T(t) en t ∈ [0, 30] junto con la l ́ınea horizontal T = 50◦C
t = np.linspace(0, 30, 100)
T = []
for i in t:
    T.append(modelNewton(i, 20, 90, 0.07))

plt.plot(t, T, color = "orange")
plt.axhline(y = 50, linestyle = ":")
plt.plot(root, modelNewton(root, 20, 90, 0.07), 'ro', markersize = 6)
plt.xlabel("t")
plt.ylabel("T(t)")
plt.title("Modelo de enfriamiento de Newton")
plt.grid(True)
plt.show()

#c.
#Graficamos |f(tn)| vs. iteración n en escala semilog
n = []
for i in range(1, count + 1):
    n.append(i)

plt.semilogy(n, fn)
plt.xlabel("n")
plt.ylabel("|T(tn)|")
plt.title("Convergencia del método Newton-Raphson")
plt.grid(True)
plt.show()

#Notamos que, al igual que el método de la bisección converge de manera exponencial, pero mucho más lento y de manera más uniforme.

#Problema 3.#

#Definimos la función f(x)
def f(x):
    return (myFunctions.E2(x) * (2 * myFunctions.ln(abs(x)))) - x

#a.
#Importamos el método de la secante de la librería myMethods
root, count, history, fn = myMethods.mSec(f, 1.5, 0.5)
#Imprimimos la raíz encontrada con el método de la secante
print(f"La raíz encontrada con el método de la secante es: {root} y lo alcanza en {count} pasos.")

#b.
#Graficamos f(x) en x ∈ [−2, 2]
x = np.linspace(-2, 2, 100)
y = []
for i in x:
    y.append(f(i))

plt.plot(x, y, color = "orange")
plt.plot(0.5, f(0.5), 'bo', markersize = 6)
plt.plot(1.5, f(1.5), 'bo', markersize = 6)
plt.plot(root, f(root), 'ro', markersize = 6)
plt.axhline(y = 0, linestyle = ":")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Gráfica de la función f(x)")
plt.grid(True)
plt.show()

#c.
#Grafique |f(xn)| vs. iteración n en escala semilog
n = []
for i in range(1, count + 1):
    n.append(i)

plt.semilogy(n, fn)
plt.xlabel("n")
plt.ylabel("|f(xn)|")
plt.title("Convergencia del método de la Secante")
plt.grid(True)
plt.show()

#El método converge de manera logarítmica, sin embargo, notamos que lo hace mucho más rápido que las otras.

#Problema 4.#

def mBN(f, df, a, b, maxIter = 1000, TOL = EPSILON):
    #Definimos la lista history para guardar el historial
    history = []
    #Tomamos un intervalo de encuadre y calculamos el valor medio
    m = a + 0.5 * (a - b)
    try:

        if f(a) * f(b) > 0:
            raise ValueError(f"No se puede asegurar que la función tenga raíz en [{a}, {b}]")
        
        for i in range(maxIter):
            #Realizamos un paso de Newton
            x = m - f(m) / df(m)
            history.append(x)

            #Si el valor converge, sale de la iteración y regresa los valores
            if abs(x - m) < TOL * max(1.0, abs(x)) or abs(a - b) < TOL:
                    return x, i+1, history

            #Si el valor de x está dentro del intervalo entonces realiza otro paso de Newton
            if a <= x <= b:
                m = x
            #Si el valor de x quiere ''salir'' del intervalo [a,b], entonces se realiza un paso de bisección
            else:
                if f(a) * f(m) < 0:
                    b = m
                else:
                    a = m
                #Calculamos el nuevo valor intermedio para realizar otro un paso de Newton
                m = a + 0.5 * (a - b)
        
        #Si no converge devuelve el último valor calculado
        print("Error: El método de Newton-Bisección no llegó a EPSILON; devolviendo última aproximación")
        return x, maxIter, history

    except Exception as e:
          print(f"[Error en el método Bisección-Newton] {e}")
          return None, 0, []
    

#Problema 5.#

#Definimos la ecuación trascendental
def h(E):
    V0 = 60
    a = 1.45
    h = (hbar*c*(1/eV*10**6)*10**15) / 1e12
    m = (2 * ((m_p * m_n) / (m_p + m_n))) * 5.61e29
    k = ((m * (V0 - abs(E))))**(1/2) / h
    b = ((m * abs(E)))**(1/2) / h

    return (k * (myFunctions.Cos(k * a) / myFunctions.Sin(k * a))) + b

def Dh(E):
    V0 = 60
    a = 1.54
    h = (hbar*c*(1/eV*10**6)*10**15) / 1e12
    m = (2 * ((m_p * m_n) / (m_p + m_n))) * 5.61e29
    k = ((m * (V0 - abs(E))) / h**2)**(1/2)
    b = ((m * abs(E)) / h**2)**(1/2)
    cot = myFunctions.Cos(k * a) / myFunctions.Sin(k * a)
    csc = 1 / myFunctions.Sin(k * a)

    try:
        E = ((- m / (2 * h**2)) * (E / abs(E))) * ((a * csc**2) + (cot / k) - (1 / b))
    except ZeroDivisionError:
        pass

    return E

#Hallamos E usando el método de Bisección
root, count, history, fn = myMethods.mBis(h, -2, 0)
print(f"La Energía que satisface la ecuación trascendental es {root} y el método de la Bisección llegó a el en {count} pasos.")

#Hallamos E usando el método de Newton-Raphson
root, count, history, fn, = myMethods.mNR(h, Dh, -1.3)
print(f"La Energía que satisface la ecuación trascendental es {root} y el método de Newton-Raphson llegó a el en {count} pasos.")

#Hallamos E usando el método de la Secante
root, count, history, fn = myMethods.mSec(h, -1.3, -1)
print(f"La Energía que satisface la ecuación trascendental es {root} y el método de la Secante llegó a el en {count} pasos.")

#Graficamos la función de la energía E
E = np.linspace(-2, 0, 200)
hE = []
for i in E:
    hE.append(h(i))

plt.plot(E, hE, color = "orange")
plt.plot(root, h(root), 'ro', markersize = "6")
plt.axhline(y = 0, linestyle = ":")
plt.xlabel("E [MeV]")
plt.ylabel("h(E)")
plt.title("Ecuación trascendental para las enegías del sistema")
plt.grid(True)
plt.show()

#De todos los métodos para calcular las raíces de la ecuación trascendental notamos que el que lo hace en un número menor de pasos, 
#   es decir, es más eficiente es el método de la secante, seguida del método de la Bisección y por último el método de Newton-Raphson
#   pues este no convergió en 1000 iteraciones como los otros dos métodos. 
#   Para este último método se cambió el número máximo de iteraciones a 10000 y se usó un valor muy cercano a la raíz como x0, sin embargo, 
#   no convergió este método.

#Notamos que V0 no puede ser menor a |E|, pues,
#                       V0 < |E|
#                     V0 - |E| < 0
#                 (V0 - |E|)^{1/2} está en C 
# 
#   así mismo, V0 no puede ser igual a |E|, pues, 
#                       cot(ka) = cot(0)
#   y cot(0) no está definido.
#   Por lo tanto, V0 > E.