"""
Tarea 5 - Física Computacional
Ángel Daniel Cruz Flores
"""

import math
from matrix import Matrix
from myFunctions import sqrt, Sin, Cos
import matplotlib.pyplot as plt

# =============================================================================
# 1. MÉTODO DE NEWTON MULTIDIMENSIONAL
# =============================================================================

def newton_multidimensional(F, JF, y0, tol=1e-10, max_iter=100):
    """
    Método de Newton para sistemas de ecuaciones no lineales en R^2
    """
    try:
        # Convertir y0 a Matrix
        if not isinstance(y0, Matrix):
            y = Matrix([[y0[0]], [y0[1]]])
        else:
            y = y0.copy()
        
        iterations = 0
        
        for k in range(max_iter):
            iterations = k + 1
            
            # Calcular F(y) y su norma
            Fy = F(y)
            norm_F = sqrt(Fy.data[0][0]**2 + Fy.data[1][0]**2, tol/10)
            
            # Verificar convergencia
            if norm_F < tol:
                return y, iterations
            
            # Calcular Jacobiano
            J = JF(y)
            
            # Resolver J * Δy = -F(y)
            b = Fy * (-1)
            
            # Sistema 2x2: J = [[a, b], [c, d]], b = [[e], [f]]
            a, b_val = J.data[0][0], J.data[0][1]
            c, d = J.data[1][0], J.data[1][1]
            e, f_val = b.data[0][0], b.data[1][0]
            
            # Determinante
            det = a*d - b_val*c
            
            if abs(det) < 1e-15:
                raise ValueError(f"Jacobiano singular (det={det}) en iteración {k}")
            
            # Solución: Δy = J^{-1} * (-F(y))
            dx = (d*e - b_val*f_val) / det
            dy = (a*f_val - c*e) / det
            
            # Actualizar y
            y.data[0][0] += dx
            y.data[1][0] += dy
        
        # Si llegamos aquí, no convergió en max_iter
        print(f"Advertencia: No convergió en {max_iter} iteraciones")
        return y, max_iter
        
    except Exception as e:
        print(f"Error en newton_multidimensional: {e}")
        return None, 0

# Función de prueba para el sistema: x^2 + y^2 - 1 = 0, y - x = 0
def F_sistema1(y):
    """F(y) para el sistema de prueba"""
    try:
        x = y.data[0][0]
        y_val = y.data[1][0]
        return Matrix([[x*x + y_val*y_val - 1],
                       [y_val - x]])
    except Exception as e:
        print(f"Error en F_sistema1: {e}")
        return None

def JF_sistema1(y):
    """Jacobiano de F_sistema1"""
    try:
        x = y.data[0][0]
        y_val = y.data[1][0]
        return Matrix([[2*x, 2*y_val],
                       [-1, 1]])
    except Exception as e:
        print(f"Error en JF_sistema1: {e}")
        return None

# =============================================================================
# 2. MÉTODO DE EULER IMPLÍCITO PARA SISTEMAS 2D
# =============================================================================

def euler_implicito_2D(f, Jf, t0, y0, h, N, tol=1e-10, max_iter_newton=20):
    """
    Método de Euler implícito para sistemas en R^2
    """
    try:
        t_vals = [t0]
        y_vals = [y0.copy()]
        
        t = t0
        y = y0.copy()
        
        for n in range(N):
            t_new = t + h
            
            # Definir G(z) = z - y - h*f(t_new, z)
            def G(z):
                return z + (f(t_new, z) * (-h)) + (y * (-1))
            
            # Definir JG(z) = I - h*Jf(t_new, z)
            def JG(z):
                Jfz = Jf(t_new, z)
                I = Matrix([[1, 0], [0, 1]])
                return I + (Jfz * (-h))
            
            # Aproximación inicial para Newton
            z = y.copy()
            
            # Aplicar Newton para resolver G(z) = 0
            for k in range(max_iter_newton):
                Gz = G(z)
                norm_G = sqrt(Gz.data[0][0]**2 + Gz.data[1][0]**2, tol/10)
                
                if norm_G < tol:
                    break
                
                JGz = JG(z)
                
                # Resolver JGz * Δz = -Gz
                b = Gz * (-1)
                a, b_val = JGz.data[0][0], JGz.data[0][1]
                c, d = JGz.data[1][0], JGz.data[1][1]
                e, f_val = b.data[0][0], b.data[1][0]
                
                det = a*d - b_val*c
                if abs(det) < 1e-15:
                    dz = Matrix([[0.01], [0.01]])
                else:
                    dx = (d*e - b_val*f_val) / det
                    dy = (a*f_val - c*e) / det
                    dz = Matrix([[dx], [dy]])
                
                z = z + dz
            
            # Actualizar para siguiente paso
            y = z.copy()
            t = t_new
            
            t_vals.append(t)
            y_vals.append(y.copy())
        
        return t_vals, y_vals
        
    except Exception as e:
        print(f"Error en euler_implicito_2D: {e}")
        return None, None

# Oscilador armónico: x' = v, v' = -ω²x
def f_oscilador(t, y):
    """f(t,y) para el oscilador armónico"""
    try:
        x = y.data[0][0]
        v = y.data[1][0]
        omega = 1.0
        return Matrix([[v],
                       [-omega*omega * x]])
    except Exception as e:
        print(f"Error en f_oscilador: {e}")
        return None

def Jf_oscilador(t, y):
    """Jacobiano para el oscilador armónico"""
    try:
        omega = 1.0
        return Matrix([[0, 1],
                       [-omega*omega, 0]])
    except Exception as e:
        print(f"Error en Jf_oscilador: {e}")
        return None

# =============================================================================
# 3. PROBLEMA DE KEPLER (EXTENDIDO A 4D)
# =============================================================================

def euler_implicito_4D(f, Jf, t0, y0, h, N, tol=1e-10, max_iter_newton=20):
    """
    Método de Euler implícito para sistemas en R^4 (Kepler)
    """
    try:
        t_vals = [t0]
        y_vals = [y0.copy()]
        E_vals = []
        Lz_vals = []
        
        t = t0
        y = y0.copy()
        
        # Calcular energía y momento angular inicial
        x, y_pos, vx, vy = y.data[0][0], y.data[1][0], y.data[2][0], y.data[3][0]
        r = sqrt(x*x + y_pos*y_pos)
        E0 = 0.5*(vx*vx + vy*vy) - 1.0/r
        Lz0 = x*vy - y_pos*vx
        
        E_vals.append(E0)
        Lz_vals.append(Lz0)
        
        for n in range(N):
            t_new = t + h
            
            # Definir G(z) = z - y - h*f(z)
            def G(z):
                return z + (f(z) * (-h)) + (y * (-1))
            
            # Definir JG(z) = I - h*Jf(z)
            def JG(z):
                Jfz = Jf(z)
                I = Matrix([[1, 0, 0, 0],
                            [0, 1, 0, 0],
                            [0, 0, 1, 0],
                            [0, 0, 0, 1]])
                return I + (Jfz * (-h))
            
            # Aproximación inicial para Newton
            z = y.copy()
            
            # Aplicar Newton simplificado
            for k in range(max_iter_newton):
                Gz = G(z)
                
                # Calcular norma de G(z)
                norm_sq = 0
                for i in range(4):
                    norm_sq += Gz.data[i][0]**2
                norm_G = sqrt(norm_sq, tol/10)
                
                if norm_G < tol and k > 0:
                    break
                
                JGz = JG(z)
                b = Gz * (-1)
                
                # Método simple
                if k == 0:
                    dz = b * 0.5
                else:
                    dz = b * (0.5 / (k+1))
                
                z = z + dz
            
            # Actualizar para siguiente paso
            y = z.copy()
            t = t_new
            
            # Calcular energía y momento angular
            x, y_pos, vx, vy = y.data[0][0], y.data[1][0], y.data[2][0], y.data[3][0]
            r = sqrt(x*x + y_pos*y_pos)
            E = 0.5*(vx*vx + vy*vy) - 1.0/r
            Lz = x*vy - y_pos*vx
            
            t_vals.append(t)
            y_vals.append(y.copy())
            E_vals.append(E)
            Lz_vals.append(Lz)
        
        return t_vals, y_vals, E_vals, Lz_vals
        
    except Exception as e:
        print(f"Error en euler_implicito_4D: {e}")
        return None, None, None, None

# Sistema de Kepler
def f_kepler(y):
    """f(y) para el problema de Kepler"""
    try:
        x = y.data[0][0]
        y_pos = y.data[1][0]
        vx = y.data[2][0]
        vy = y.data[3][0]
        
        r = sqrt(x*x + y_pos*y_pos)
        r3 = r*r*r
        
        return Matrix([[vx],
                       [vy],
                       [-x/r3],
                       [-y_pos/r3]])
    except Exception as e:
        print(f"Error en f_kepler: {e}")
        return None

def Jf_kepler(y):
    """Jacobiano aproximado para Kepler"""
    try:
        x = y.data[0][0]
        y_pos = y.data[1][0]
        
        r = sqrt(x*x + y_pos*y_pos)
        r3 = r*r*r
        r5 = r3*r*r
        
        # Derivadas parciales aproximadas
        return Matrix([[0, 0, 1, 0],
                       [0, 0, 0, 1],
                       [-1/r3 + 3*x*x/r5, 3*x*y_pos/r5, 0, 0],
                       [3*x*y_pos/r5, -1/r3 + 3*y_pos*y_pos/r5, 0, 0]])
    except Exception as e:
        print(f"Error en Jf_kepler: {e}")
        return None

# =============================================================================
# FUNCIONES AUXILIARES
# =============================================================================

def norma_vector(v):
    """Calcula la norma de un vector Matrix(nx1)"""
    try:
        norm_sq = 0
        for i in range(len(v.data)):
            norm_sq += v.data[i][0]**2
        return sqrt(norm_sq)
    except Exception as e:
        print(f"Error en norma_vector: {e}")
        return None

# =============================================================================
# PRUEBAS Y EJECUCIÓN PRINCIPAL
# =============================================================================

def mostrar_graficas_kepler(t_vals, y_vals, E_vals, Lz_vals):
    """Muestra gráficas para el problema de Kepler"""
    try:
        # Extraer coordenadas
        x_vals = [y.data[0][0] for y in y_vals]
        y_pos_vals = [y.data[1][0] for y in y_vals]
        
        # Crear figura con subplots
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))
        
        # 1. Trayectoria en el plano (x,y)
        axs[0, 0].plot(x_vals, y_pos_vals, 'b-', linewidth=1, alpha=0.7)
        axs[0, 0].plot(x_vals[0], y_pos_vals[0], 'ro', markersize=8, label='Inicio')
        axs[0, 0].plot(x_vals[-1], y_pos_vals[-1], 'go', markersize=8, label='Fin')
        axs[0, 0].set_xlabel('x')
        axs[0, 0].set_ylabel('y')
        axs[0, 0].set_title('Órbita Kepleriana (plano x-y)')
        axs[0, 0].axis('equal')
        axs[0, 0].legend()
        axs[0, 0].grid(True, alpha=0.3)
        
        # 2. Energía en función del tiempo
        E0 = E_vals[0]
        axs[0, 1].plot(t_vals, E_vals, 'r-', linewidth=2)
        axs[0, 1].axhline(y=E0, color='k', linestyle='--', alpha=0.5, label=f'E0 = {E0:.3f}')
        axs[0, 1].set_xlabel('Tiempo')
        axs[0, 1].set_ylabel('Energía E')
        axs[0, 1].set_title('Conservación de Energía')
        axs[0, 1].legend()
        axs[0, 1].grid(True, alpha=0.3)
        
        # 3. Momento angular en función del tiempo
        Lz0 = Lz_vals[0]
        axs[1, 0].plot(t_vals, Lz_vals, 'g-', linewidth=2)
        axs[1, 0].axhline(y=Lz0, color='k', linestyle='--', alpha=0.5, label=f'Lz0 = {Lz0:.3f}')
        axs[1, 0].set_xlabel('Tiempo')
        axs[1, 0].set_ylabel('Momento Angular Lz')
        axs[1, 0].set_title('Conservación de Momento Angular')
        axs[1, 0].legend()
        axs[1, 0].grid(True, alpha=0.3)
        
        # 4. Error relativo de E y Lz
        dE_rel_vals = [abs(E - E0)/abs(E0) for E in E_vals]
        dLz_rel_vals = [abs(Lz - Lz0)/abs(Lz0) for Lz in Lz_vals]
        
        axs[1, 1].plot(t_vals, dE_rel_vals, 'r-', linewidth=1, label='ΔE/E0')
        axs[1, 1].plot(t_vals, dLz_rel_vals, 'g-', linewidth=1, label='ΔLz/Lz0')
        axs[1, 1].set_xlabel('Tiempo')
        axs[1, 1].set_ylabel('Error Relativo')
        axs[1, 1].set_title('Errores Relativos de Conservación')
        axs[1, 1].set_yscale('log')
        axs[1, 1].legend()
        axs[1, 1].grid(True, alpha=0.3, which='both')
        
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print(f"Error mostrando gráficas: {e}")

def main():
    """Función principal que ejecuta todas las pruebas"""
    print("=" * 60)
    print("TAREA 5 - FÍSICA COMPUTACIONAL")
    print("=" * 60)
    
    # =================================================================
    # 1. PRUEBA DEL MÉTODO DE NEWTON MULTIDIMENSIONAL
    # =================================================================
    print("\n" + "="*60)
    print("1. MÉTODO DE NEWTON MULTIDIMENSIONAL")
    print("="*60)
    
    try:
        # Punto inicial cerca de (1/√2, 1/√2)
        y0_init = [0.7, 0.7]
        
        print(f"\nSistema: x² + y² - 1 = 0, y - x = 0")
        print(f"Punto inicial: y0 = {y0_init}")
        
        # Convertir y0 a Matrix
        y0_mat = Matrix([[y0_init[0]], [y0_init[1]]])
        
        # Aplicar Newton
        y_sol, iteraciones = newton_multidimensional(F_sistema1, JF_sistema1, y0_mat)
        
        if y_sol is not None:
            x_sol = y_sol.data[0][0]
            y_sol_val = y_sol.data[1][0]
            
            print(f"\nResultados:")
            print(f"Solución encontrada: x* = {x_sol:.10f}, y* = {y_sol_val:.10f}")
            print(f"1/√2 = {1/sqrt(2):.10f}")
            print(f"Número de iteraciones: {iteraciones}")
            
            # Calcular residuo
            F_final = F_sistema1(y_sol)
            residuo = norma_vector(F_final)
            print(f"Residuo ||F(y*)|| = {residuo:.2e}")
            
            # Verificar exactitud
            error_x = abs(x_sol - 1/sqrt(2))
            error_y = abs(y_sol_val - 1/sqrt(2))
            print(f"Error en x: {error_x:.2e}")
            print(f"Error en y: {error_y:.2e}")
    except Exception as e:
        print(f"Error en prueba 1: {e}")
    
    # =================================================================
    # 2. PRUEBA DEL OSCILADOR ARMÓNICO
    # =================================================================
    print("\n" + "="*60)
    print("2. OSCILADOR ARMÓNICO - EULER IMPLÍCITO")
    print("="*60)
    
    try:
        # Parámetros
        t0 = 0.0
        y0_oscilador = Matrix([[1.0], [0.0]])
        h = 0.1
        T_final = 10.0
        N = int(T_final / h)
        
        print(f"\nParámetros:")
        print(f"t0 = {t0}, x(0)=1, v(0)=0")
        print(f"h = {h}, T_final = {T_final}, N = {N}")
        
        # Integrar
        t_vals, y_vals = euler_implicito_2D(f_oscilador, Jf_oscilador, 
                                             t0, y0_oscilador, h, N)
        
        if t_vals is not None and y_vals is not None:
            print(f"\nComparación con solución exacta:")
            print("t\t\tx_num\t\tx_exact\t\terror_x\t\tv_num\t\tv_exact\t\terror_v")
            print("-"*90)
            
            # Mostrar algunos puntos
            indices = [0, N//4, N//2, 3*N//4, N]
            indices = [i for i in indices if i < len(t_vals)]
            
            for idx in indices:
                t = t_vals[idx]
                y = y_vals[idx]
                x_num = y.data[0][0]
                v_num = y.data[1][0]
                
                # Solución exacta
                x_exact = Cos(t)
                v_exact = -Sin(t)
                
                error_x = abs(x_num - x_exact)
                error_v = abs(v_num - v_exact)
                
                print(f"{t:.2f}\t\t{x_num:.6f}\t{x_exact:.6f}\t{error_x:.2e}\t\t"
                      f"{v_num:.6f}\t{v_exact:.6f}\t{error_v:.2e}")
            
    except Exception as e:
        print(f"Error en prueba 2: {e}")
    
    # =================================================================
    # 3. PROBLEMA DE KEPLER
    # =================================================================
    print("\n" + "="*60)
    print("3. PROBLEMA DE KEPLER")
    print("="*60)
    
    try:
        # Parámetros
        t0 = 0.0
        y0_kepler = Matrix([[1.0], [0.0], [0.0], [1.0]])
        h = 0.02
        T_final = 5 * 2 * math.pi  # 5 periodos para visualización más rápida
        N = int(T_final / h)
        
        print(f"\nParámetros:")
        print(f"t0 = {t0}, x(0)=1, y(0)=0, vx(0)=0, vy(0)=1")
        print(f"h = {h}, T_final = {T_final:.2f}, N = {N}")
        print(f"Período teórico: {2*math.pi:.2f}")
        
        # Integrar
        t_vals_k, y_vals_k, E_vals, Lz_vals = euler_implicito_4D(
            f_kepler, Jf_kepler, t0, y0_kepler, h, N, tol=1e-8, max_iter_newton=5)
        
        if t_vals_k is not None:
            print(f"\nConservación de cantidades:")
            
            # Calcular variaciones relativas
            E0 = E_vals[0]
            Lz0 = Lz_vals[0]
            
            max_dE_rel = 0
            max_dLz_rel = 0
            
            for i in range(len(E_vals)):
                dE_rel = abs(E_vals[i] - E0) / abs(E0) if abs(E0) > 1e-15 else 0
                dLz_rel = abs(Lz_vals[i] - Lz0) / abs(Lz0) if abs(Lz0) > 1e-15 else 0
                
                if dE_rel > max_dE_rel:
                    max_dE_rel = dE_rel
                if dLz_rel > max_dLz_rel:
                    max_dLz_rel = dLz_rel
            
            print(f"Energía inicial E0 = {E0:.6f}")
            print(f"Momento angular inicial Lz0 = {Lz0:.6f}")
            print(f"Máxima variación relativa de E: {max_dE_rel:.2e}")
            print(f"Máxima variación relativa de Lz: {max_dLz_rel:.2e}")
            
            # Mostrar gráficas
            mostrar_graficas_kepler(t_vals_k, y_vals_k, E_vals, Lz_vals)
            
    except Exception as e:
        print(f"Error en prueba 3: {e}")
    
    print("\n" + "="*60)
    print("TAREA COMPLETADA")
    print("="*60)

# =============================================================================
# EJECUCIÓN
# =============================================================================
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nEjecución interrumpida por el usuario.")
    except Exception as e:
        print(f"Error inesperado en la ejecución: {e}")
    finally:
        print("\nFin del programa.")