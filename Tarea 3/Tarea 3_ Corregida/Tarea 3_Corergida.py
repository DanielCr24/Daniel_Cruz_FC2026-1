import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, c, eV
import myFunctions as mf
import myMethods as mm 

# Configuración general
TOL = sys.float_info.epsilon

# =============================================================================
# MÉTODOS NUMÉRICOS
# =============================================================================

def bisection(f, a, b, tol=TOL, max_iter=1000):
    """Método de bisección"""
    try:
        raiz, iteraciones, historia, fvals = mm.mBis(f, a, b, max_iter, tol)
        return raiz, list(range(iteraciones)), fvals
    except Exception as e:
        print(f"Error en bisección: {e}")
        return None, [], []

def newton_raphson(f, df, x0, tol=TOL, max_iter=1000):
    """Método de Newton-Raphson"""
    try:
        raiz, iteraciones, historia, fvals = mm.mNR1(f, df, x0, max_iter, tol)
        return raiz, list(range(iteraciones)), fvals
    except Exception as e:
        print(f"Error en Newton-Raphson: {e}")
        return None, [], []

def secante(f, x0, x1, tol=TOL, max_iter=1000):
    """Método de la secante"""
    try:
        raiz, iteraciones, historia, fvals = mm.mSec(f, x1, x0, max_iter, tol)
        return raiz, list(range(iteraciones)), fvals
    except Exception as e:
        print(f"Error en método de la secante: {e}")
        return None, [], []

def newton_seguro(f, df, a, b, tol=TOL, max_iter=1000):
    """
    Método de Newton seguro que combina con bisección
    """
    try:
        if f(a) * f(b) >= 0:
            raise ValueError("La función debe cambiar de signo en el intervalo [a,b]")
        
        x = (a + b) / 2  # Punto inicial en el medio del intervalo
        iterations = []
        function_values = []
        
        for i in range(max_iter):
            fx = f(x)
            iterations.append(i)
            function_values.append(abs(fx))
            
            # Verificar convergencia
            if abs(fx) < tol or (b - a) < tol * max(1.0, abs(x)):
                return x, iterations, function_values
            
            # Intentar paso de Newton
            try:
                dfx = df(x)
                
                if abs(dfx) > tol:  # Solo usar Newton si la derivada es significativa
                    x_newton = x - fx / dfx
                    
                    # Verificar si x_newton está dentro del intervalo y mejora la solución
                    if a < x_newton < b and abs(f(x_newton)) < abs(fx):
                        x_new = x_newton
                    else:
                        # Si Newton no es bueno, usar bisección
                        x_new = (a + b) / 2
                else:
                    # Derivada muy pequeña, usar bisección
                    x_new = (a + b) / 2
            except:
                # En caso de error, usar bisección
                x_new = (a + b) / 2
            
            # Actualizar intervalo
            if f(a) * f(x_new) < 0:
                b = x_new
            else:
                a = x_new
            
            x = x_new
        
        raise RuntimeError("No se alcanzó convergencia en el número máximo de iteraciones")
    
    except Exception as e:
        print(f"Error en Newton seguro: {e}")
        return None, [], []

# =============================================================================
# EJERCICIO 1: Vibraciones de viga en voladizo
# =============================================================================

def ejercicio1():
    try:
        # Definir la función f(ω)
        def f(omega):
            return mf.Cos(omega) * mf.Cosh(omega) + 1
        
        # Encontrar la raíz
        a, b = 0, 2  # Intervalo inicial
        raiz, iter_bisec, fvals_bisec = bisection(f, a, b)
        
        if raiz is not None:
            print(f"Ejercicio 1a - Raíz encontrada: ω = {raiz:.6f}")
            print(f"Número de iteraciones: {len(iter_bisec)}")
            
            # b) Gráfica de f(ω)
            omega_vals = np.linspace(0, 4, 1000)
            f_vals = [f(omega) for omega in omega_vals]
            
            plt.figure(figsize=(12, 4))
            
            plt.subplot(1, 2, 1)
            plt.plot(omega_vals, f_vals, 'b-', linewidth=2, label='f(ω)')
            plt.axhline(y=0, color='k', linestyle='--', alpha=0.7)
            plt.axvspan(a, b, alpha=0.3, color='red', label='Intervalo inicial')
            plt.axvline(x=raiz, color='green', linestyle='--', label=f'Raíz ≈ {raiz:.4f}')
            plt.xlabel('ω')
            plt.ylabel('f(ω)')
            plt.title('Vibraciones de viga en voladizo')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            # c) Gráfica de convergencia
            plt.subplot(1, 2, 2)
            plt.semilogy(iter_bisec, fvals_bisec, 'ro-', linewidth=2, markersize=4, label='|f(ωₙ)|')
            plt.axhline(y=TOL, color='k', linestyle='--', alpha=0.7, label=f'TOL = {TOL:.1e}')
            plt.xlabel('Iteración (n)')
            plt.ylabel('|f(ωₙ)|')
            plt.title('Convergencia del método de bisección')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.show()
        
    except Exception as e:
        print(f"Error en ejercicio 1: {e}")
    finally:
        print("\nEl método de bisección exhibe un régimen de convergencia lineal garantizado, siendo el más robusto pero también el más lento entre los métodos estudiados. Su convergencia es monótona y predecible, reduciendo el error a la mitad en cada iteración, lo que se refleja en las gráficas semilogarítmicas como una disminución constante del error con pendiente uniforme. A diferencia de Newton-Raphson y la secante, la bisección no acelera su convergencia cerca de la raíz, manteniendo siempre la misma tasa de convergencia independientemente de la suavidad de la función o la proximidad a la solución. Esta característica lo hace extremadamente confiable - siempre converge cuando la función cambia de signo en el intervalo inicial - pero computacionalmente costoso para alcanzar alta precisión, requiriendo aproximadamente 3-4 iteraciones adicionales por cada dígito decimal de precisión deseado.")
        print("\nEjercicio 1 completado")

# =============================================================================
# EJERCICIO 2: Modelo de enfriamiento de Newton
# =============================================================================

def ejercicio2():
    try:
        # Parámetros del problema
        T_a = 20  # °C
        T_0 = 90  # °C  
        k = 0.07  # min⁻¹
        
        # Definir T(t) y f(t)
        def T(t):
            return T_a + (T_0 - T_a) * mf.E(-k * t)
        
        def f(t):
            return T(t) - 50
        
        # Derivada de f(t)
        def df(t):
            return -k * (T_0 - T_a) * mf.E(-k * t)
        
        # Encontrar la raíz
        t0 = 10
        t_raiz, iter_newton, fvals_newton = newton_raphson(f, df, t0)
        
        if t_raiz is not None:
            print(f"\nEjercicio 2a - Tiempo encontrado: t = {t_raiz:.6f} min")
            print(f"Temperatura en t = {t_raiz:.4f}: T = {T(t_raiz):.4f} °C")
            print(f"Número de iteraciones: {len(iter_newton)}")
            
            # b) Gráfica de T(t)
            t_vals = np.linspace(0, 30, 1000)
            T_vals = [T(t) for t in t_vals]
            
            plt.figure(figsize=(12, 4))
            
            plt.subplot(1, 2, 1)
            plt.plot(t_vals, T_vals, 'b-', linewidth=2, label='T(t)')
            plt.axhline(y=50, color='r', linestyle='--', label='T = 50°C')
            plt.axvline(x=t_raiz, color='green', linestyle='--', label=f't* ≈ {t_raiz:.2f} min')
            plt.xlabel('Tiempo (min)')
            plt.ylabel('Temperatura (°C)')
            plt.title('Enfriamiento de Newton')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            # c) Gráfica de convergencia
            plt.subplot(1, 2, 2)
            plt.semilogy(iter_newton, fvals_newton, 'ro-', linewidth=2, markersize=4, label='|f(tₙ)|')
            plt.axhline(y=TOL, color='k', linestyle='--', alpha=0.7, label=f'TOL = {TOL:.1e}')
            plt.xlabel('Iteración (n)')
            plt.ylabel('|f(tₙ)|')
            plt.title('Convergencia del método de Newton-Raphson')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.show()
        
    except Exception as e:
        print(f"Error en ejercicio 2: {e}")
    finally:
        print("\n El método de Newton-Raphson exhibe un régimen de convergencia cuadrático cuando se cumplen las condiciones ideales: la función es suficientemente diferenciable, la derivada no se anula en la vecindad de la raíz, y el punto inicial está lo suficientemente cerca de la solución. Esta convergencia cuadrática se manifiesta en las gráficas semilogarítmicas como una disminución extremadamente rápida del error, donde el número de dígitos correctos se duplica aproximadamente en cada iteración. Sin embargo, cuando el punto inicial está lejos de la raíz o la función presenta regiones de alta no linealidad, el método puede exhibir convergencia lineal inicial antes de transitionar al régimen cuadrático, o incluso divergir si la derivada se aproxima a cero o si la función tiene puntos de inflexión cerca del punto inicial. La robustez del método depende críticamente de la elección del punto inicial y del comportamiento local de la función y su derivada en la región de búsqueda.")
        print("\nEjercicio 2 completado")

# =============================================================================
# EJERCICIO 3: Método de la secante (ambas raíces)
# =============================================================================

def ejercicio3():
    try:
        # Definir la función f(x)
        def f(x):
            return (mf.E2(x) * (2 * mf.ln(abs(x)))) - x
        
        # Encontrar la raíz positiva
        x0_pos, x1_pos = 0.5, 1.5
        raiz_pos, iter_secante_pos, fvals_secante_pos = secante(f, x0_pos, x1_pos)
        
        # Encontrar la raíz negativa
        x0_neg, x1_neg = -1.5, -0.5
        raiz_neg, iter_secante_neg, fvals_secante_neg = secante(f, x0_neg, x1_neg)
        
        print(f"\nEjercicio 3 - Raíces encontradas:")
        if raiz_pos is not None:
            print(f"Raíz positiva: x = {raiz_pos:.6f}, f(x) = {f(raiz_pos):.2e}")
            print(f"Iteraciones raíz positiva: {len(iter_secante_pos)}")
        
        if raiz_neg is not None:
            print(f"Raíz negativa: x = {raiz_neg:.6f}, f(x) = {f(raiz_neg):.2e}")
            print(f"Iteraciones raíz negativa: {len(iter_secante_neg)}")
        
        # a) Gráfica de f(x)
        x_vals = np.linspace(-2, 2, 1000)
        f_vals = [f(x) for x in x_vals]
        
        plt.figure(figsize=(15, 5))
        
        # Gráfica de la función
        plt.subplot(1, 3, 1)
        plt.plot(x_vals, f_vals, 'b-', linewidth=2, label='f(x)')
        plt.axhline(y=0, color='k', linestyle='--', alpha=0.7)
        
        if raiz_pos is not None:
            plt.axvline(x=raiz_pos, color='green', linestyle='--', label=f'Raíz + ≈ {raiz_pos:.4f}')
            plt.plot([x0_pos, x1_pos], [f(x0_pos), f(x1_pos)], 'go', markersize=6, label='Puntos iniciales +')
        
        if raiz_neg is not None:
            plt.axvline(x=raiz_neg, color='red', linestyle='--', label=f'Raíz - ≈ {raiz_neg:.4f}')
            plt.plot([x0_neg, x1_neg], [f(x0_neg), f(x1_neg)], 'ro', markersize=6, label='Puntos iniciales -')
        
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.title('Función f(x) = exp(x²)ln(x²) - x')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Gráfica de convergencia raíz positiva
        plt.subplot(1, 3, 2)
        if raiz_pos is not None:
            plt.semilogy(iter_secante_pos, fvals_secante_pos, 'go-', linewidth=2, markersize=4, label='Raíz positiva')
            plt.axhline(y=TOL, color='k', linestyle='--', alpha=0.7, label=f'TOL = {TOL:.1e}')
            plt.xlabel('Iteración (n)')
            plt.ylabel('|f(xₙ)|')
            plt.title('Convergencia - Raíz positiva')
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        # Gráfica de convergencia raíz negativa
        plt.subplot(1, 3, 3)
        if raiz_neg is not None:
            plt.semilogy(iter_secante_neg, fvals_secante_neg, 'ro-', linewidth=2, markersize=4, label='Raíz negativa')
            plt.axhline(y=TOL, color='k', linestyle='--', alpha=0.7, label=f'TOL = {TOL:.1e}')
            plt.xlabel('Iteración (n)')
            plt.ylabel('|f(xₙ)|')
            plt.title('Convergencia - Raíz negativa')
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print(f"Error en ejercicio 3: {e}")
    finally:
        print("\nEl método de la secante muestra una sensibilidad moderada a la elección de los puntos iniciales, siendo más robusto que Newton-Raphson pero menos que la bisección. La convergencia óptima se obtiene cuando los puntos iniciales se encuentran en lados opuestos de la raíz, similar al concepto de encerrar la raíz en la bisección, ya que esto garantiza que la secante cruce el eje x cerca de la solución verdadera. Puntos muy cercanos entre sí pueden causar inestabilidad numérica debido a divisiones por valores pequeños en el cálculo de la pendiente, mientras que puntos muy alejados pueden converger más lentamente si la función es altamente no lineal. La peor situación ocurre cuando ambos puntos iniciales se localizan del mismo lado de la raíz, particularmente en regiones de baja curvatura, donde el método puede divergir o requerir significativamente más iteraciones. La efectividad del método depende críticamente de que los puntos iniciales proporcionen una aproximación razonable de la comportamiento local de la función cerca de la raíz.")
        print("\nEjercicio 3 completado")

# =============================================================================
# EJERCICIO 4: Método de Newton seguro mejorado
# =============================================================================

def ejercicio4():
    try:
        # Probar con una función ''fácil''
        def f_ejemplo(x):
            return x**3 - 2*x - 5
        
        def df_ejemplo(x):
            return 3*x**2 - 2
        
        # También probar con una función ''difícil''
        def f_dificil(x):
            return x**4 - 3*x**3 + x**2 - 5*x + 2
        
        def df_dificil(x):
            return 4*x**3 - 9*x**2 + 2*x - 5
        
        print("\nEjercicio 4 - Método de Newton seguro mejorado")
        print("=" * 50)
        
        # Prueba con función fácil
        a1, b1 = 1, 3
        raiz1, iter_newton_seg1, fvals_newton_seg1 = newton_seguro(f_ejemplo, df_ejemplo, a1, b1)
        
        if raiz1 is not None:
            print(f"Función fácil: raíz = {raiz1:.6f}, f(x) = {f_ejemplo(raiz1):.2e}")
            print(f"Iteraciones: {len(iter_newton_seg1)}")
        
        # Prueba con función difícil
        a2, b2 = 0, 2
        raiz2, iter_newton_seg2, fvals_newton_seg2 = newton_seguro(f_dificil, df_dificil, a2, b2)
        
        if raiz2 is not None:
            print(f"Función difícil: raíz = {raiz2:.6f}, f(x) = {f_dificil(raiz2):.2e}")
            print(f"Iteraciones: {len(iter_newton_seg2)}")
        
        # Gráfica de convergencia comparativa
        plt.figure(figsize=(12, 5))
        
        plt.subplot(1, 2, 1)
        if raiz1 is not None:
            plt.semilogy(iter_newton_seg1, fvals_newton_seg1, 'bo-', linewidth=2, markersize=4, label='Función fácil')
            plt.axhline(y=TOL, color='k', linestyle='--', alpha=0.7, label=f'TOL = {TOL:.1e}')
            plt.xlabel('Iteración (n)')
            plt.ylabel('|f(xₙ)|')
            plt.title('Newton seguro - Función fácil')
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 2, 2)
        if raiz2 is not None:
            plt.semilogy(iter_newton_seg2, fvals_newton_seg2, 'ro-', linewidth=2, markersize=4, label='Función difícil')
            plt.axhline(y=TOL, color='k', linestyle='--', alpha=0.7, label=f'TOL = {TOL:.1e}')
            plt.xlabel('Iteración (n)')
            plt.ylabel('|f(xₙ)|')
            plt.title('Newton seguro - Función difícil')
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print(f"Error en ejercicio 4: {e}")
    finally:
        print("Ejercicio 4 completado")

# =============================================================================
# EJERCICIO 5: Ecuación de Schrödinger para el deuterón
# =============================================================================

def ejercicio5():
    try:
        # Energía del estado ligado conocido
        E_known = -2.223  # MeV
        
        def h_func(E):
            V0 = 60
            a = 1.45
            m = 938
            h = (hbar*c*(1/eV*10**6)*10**15) / 1e12 #= 197.3269804 MeV·fm

            k = mf.sqrt((m * (V0 - abs(E))) / (h**2)) 
            B = mf.sqrt((m * abs(E)) / (h**2))    

            return k * (mf.Cos(k * a) / mf.Sin(k * a)) + B

        def h_derivative(E):

            """
            Derivada analítica de h(E) para E < 0
            """
            V0 = 60
            a = 1.45
            m = 938
            h = 197.3269804  # MeV·fm
            h2 = h**2  # ħ²
            
            # Para E < 0, |E| = -E
            k = mf.sqrt((m * (V0 - abs(E))) / h2) 
            B = mf.sqrt((m * abs(E)) / h2)
            
            ka = k * a
            
            # Términos trigonométricos
            cos_ka = mf.Cos(ka)
            sin_ka = mf.Sin(ka)
            cot_ka = cos_ka / sin_ka
            csc2_ka = 1.0 / (sin_ka**2)  # cosec²(ka) = 1/sin²(ka)
            
            # Término constante
            constant = m / (2 * h2)
            
            return constant * (cot_ka / k - a * csc2_ka - 1.0 / B)
        
        # Gráfica para encontrar autovalores
        print("Generando gráfica de la ecuación trascendental...")
        E = np.linspace(-2, 0, 200)
        hE = []
        for i in E:
            hE.append(h_func(i))

        plt.plot(E, hE, color = "orange")
        plt.plot(-1.380574, h_func(-1.380574), 'ro', markersize = "6")
        plt.axhline(y = 0, linestyle = ":")
        plt.xlabel("E [MeV]")
        plt.ylabel("h(E)")
        plt.title("Ecuación trascendental para las enegías del sistema")
        plt.grid(True)
        plt.show()
        
        # Encontrar raíz usando los tres métodos
        print(f"\nEjercicio 5 - Búsqueda del estado ligado del deuterón")
        print("=" * 60)
        
        # Intervalo para bisección
        E_a, E_b = -10.0, -1.0
        
        # Método de bisección
        print("\n1. MÉTODO DE BISECCIÓN:")
        raiz_bisec, iter_bisec, fvals_bisec = bisection(h_func, E_a, E_b)
        if raiz_bisec is not None:
            print(f"   Energía encontrada: E = {raiz_bisec:.6f} MeV")
            print(f"   Iteraciones: {len(iter_bisec)}")
            print(f"   Error respecto al valor conocido: {abs(raiz_bisec - E_known):.6f} MeV")
            print(f"   h(E) = {h_func(raiz_bisec):.2e}")
        else:
            print("   No convergió")
            raiz_bisec, iter_bisec, fvals_bisec = None, [], []
        
        # Método de Newton-Raphson
        print("\n2. MÉTODO DE NEWTON-RAPHSON:")
        E0_newton = -1
        raiz_newton, iter_newton, fvals_newton = newton_raphson(h_func, h_derivative, E0_newton)
        
        if raiz_newton is not None:
            print(f"   Energía encontrada: E = {raiz_newton:.6f} MeV")
            print(f"   Iteraciones: {len(iter_newton)}")
            print(f"   Error respecto al valor conocido: {abs(raiz_newton - E_known):.6f} MeV")
            print(f"   h(E) = {h_func(raiz_newton):.2e}")
        else:
            print("   No convergió")
            raiz_newton, iter_newton, fvals_newton = None, [], []
        
        # Método de la secante
        print("\n3. MÉTODO DE LA SECANTE:")
        E0_secante, E1_secante = -5.0, -2.0
        raiz_secante, iter_secante, fvals_secante = secante(h_func, E0_secante, E1_secante)
        
        if raiz_secante is not None:
            print(f"   Energía encontrada: E = {raiz_secante:.6f} MeV")
            print(f"   Iteraciones: {len(iter_secante)}")
            print(f"   Error respecto al valor conocido: {abs(raiz_secante - E_known):.6f} MeV")
            print(f"   h(E) = {h_func(raiz_secante):.2e}")
        else:
            print("   No convergió")
            raiz_secante, iter_secante, fvals_secante = None, [], []
        
        # Gráfica comparativa de convergencia
        plt.figure(figsize=(12, 8))
        
        if raiz_bisec is not None and len(iter_bisec) > 0:
            plt.semilogy(iter_bisec, fvals_bisec, 'ro-', label='Bisección', linewidth=2, markersize=6)
        
        if raiz_newton is not None and len(iter_newton) > 0:
            plt.semilogy(iter_newton, fvals_newton, 'go-', label='Newton-Raphson', linewidth=2, markersize=6)
        
        if raiz_secante is not None and len(iter_secante) > 0:
            plt.semilogy(iter_secante, fvals_secante, 'bo-', label='Secante', linewidth=2, markersize=6)
        
        plt.axhline(y=TOL, color='k', linestyle='--', alpha=0.7, label=f'TOL = {TOL:.1e}')
        plt.xlabel('Iteración (n)')
        plt.ylabel('|h(Eₙ)|')
        plt.title('Comparación de convergencia para la ecuación del deuterón')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

        # CÁLCULO DEL VALOR MÍNIMO DE V0

        def h_func(E, V0=60):
            a = 1.45
            m = 938
            h = (hbar*c*(1/eV*10**6)*10**15) / 1e12  # = 197.3269804 MeV·fm

            k = mf.sqrt((m * (V0 - abs(E))) / (h**2)) 
            B = mf.sqrt((m * abs(E)) / (h**2))    

            return k * (mf.Cos(k * a) / mf.Sin(k * a)) + B
        
        print("\n" + "="*60)
        print("CÁLCULO DEL VALOR MÍNIMO DE V0")
        print("="*60)
        
        def tiene_estado_ligado(V0_test):
            """Verifica si existe estado ligado para un V0 dado"""
            # Para V0 mínimo, la energía tiende a 0
            E_test = -1e-14
            h_valor = h_func(E_test, V0_test)
            return h_valor < 0  # Si h(E) < 0, hay estado ligado
        
        # Búsqueda simple del V0 mínimo
        V0_min = 40.0
        V0_max = 60.0
        
        for i in range(2000):
            V0_medio = (V0_min + V0_max) / 2
            if tiene_estado_ligado(V0_medio):
                V0_max = V0_medio
            else:
                V0_min = V0_medio
        
        V0_min_calculado = (V0_min + V0_max) / 2
        print(f"Valor mínimo de V0: {V0_min_calculado:.2f} MeV")
        
    except Exception as e:
        print(f"Error en ejercicio 5: {e}")
        import traceback
        traceback.print_exc()
    finally:
        print("\nEjercicio 5 completado")

# =============================================================================
# EJECUCIÓN PRINCIPAL
# =============================================================================

if __name__ == "__main__":
    print("\nTAREA 3 - FÍSICA COMPUTACIONAL")
    print("Ángel Daniel Cruz Flores")
    print("=" * 50)
    
    # Ejecutar todos los ejercicios
    ejercicio1()
    ejercicio2() 
    ejercicio3()
    ejercicio4()
    ejercicio5()