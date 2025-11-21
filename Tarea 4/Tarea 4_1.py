import numpy as np
import time
import matplotlib.pyplot as plt
from matrix import Matrix
from myLinearMethods import lu_solve_pp
import sys

class HeatConduction1D:
    """
    Resuelve la ecuación de conducción de calor 1D estacionaria:
    -κ d²T/dx² = q(x), T(0)=100, T(1)=0
    """
    
    def __init__(self, L=1.0, kappa=1.0, T0=100.0, T1=0.0):
        self.L = L
        self.kappa = kappa
        self.T0 = T0
        self.T1 = T1
    
    def q_source(self, x):
        """Fuente volumétrica q(x) = 50 sin(πx)"""
        return 50 * np.sin(np.pi * x)
    
    def assemble_system(self, N):
        """
        Ensambla el sistema lineal AT = b para N nodos interiores
        
        Parameters
        ----------
        N : int
            Número de nodos interiores
            
        Returns
        -------
        A : Matrix
            Matriz del sistema (tridiagonal)
        b : list
            Vector del lado derecho
        x_nodes : numpy.array
            Posiciones de los nodos (incluyendo fronteras)
        """
        h = self.L / (N + 1)  # Paso espacial
        x_nodes = np.linspace(0, self.L, N + 2)  # Nodos incluyendo fronteras
        
        # Inicializar matriz A (tridiagonal)
        A_data = [[0.0] * N for _ in range(N)]
        
        # Llenar matriz A usando esquema de tres puntos
        for i in range(N):
            if i > 0:
                A_data[i][i-1] = -self.kappa / h**2  # Subdiagonal
            A_data[i][i] = 2 * self.kappa / h**2     # Diagonal
            if i < N - 1:
                A_data[i][i+1] = -self.kappa / h**2  # Superdiagonal
        
        A = Matrix(A_data)
        
        # Construir vector b
        b = [0.0] * N
        for i in range(N):
            x = x_nodes[i+1]  # Nodos interiores
            b[i] = self.q_source(x)
            
            # Añadir contribución de condiciones de frontera
            if i == 0:
                b[i] += (self.kappa / h**2) * self.T0
            if i == N - 1:
                b[i] += (self.kappa / h**2) * self.T1
        
        return A, b, x_nodes
    
    def solve_with_lu(self, N):
        """
        Resuelve el sistema usando factorización LUP
        
        Returns
        -------
        T_lu : list
            Solución con LU (solo nodos interiores)
        x_nodes : numpy.array
            Posiciones de los nodos
        time_lu : float
            Tiempo de ejecución
        operations_estimate : int
            Estimación de operaciones
        """
        A, b, x_nodes = self.assemble_system(N)
        
        # Resolver con LU
        start_time = time.time()
        T_lu = lu_solve_pp(A, b)
        time_lu = time.time() - start_time
        
        # Estimación de operaciones: O(N³) para LU
        operations_estimate = N**3
        
        return T_lu, x_nodes, time_lu, operations_estimate
    
    def solve_with_numpy(self, N):
        """
        Resuelve el sistema usando numpy.linalg.solve como referencia
        """
        A, b, x_nodes = self.assemble_system(N)
        
        # Convertir a arrays de numpy
        A_np = np.array(A.data)
        b_np = np.array(b)
        
        T_np = np.linalg.solve(A_np, b_np)
        
        return T_np.tolist()
    
    def calculate_error(self, T_lu, T_ref):
        """Calcula el error con la norma infinito ∥T_LU - T_ref∥∞"""
        error = max(abs(T_lu[i] - T_ref[i]) for i in range(len(T_lu)))
        return error
    
    def plot_solutions(self, N_values):
        """
        Grafica las soluciones para diferentes valores de N
        """
        plt.figure(figsize=(12, 8))
        
        for N in N_values:
            # Obtener solución
            T_lu, x_nodes, time_lu, _ = self.solve_with_lu(N)
            
            # Construir solución completa (incluyendo fronteras)
            T_full = [self.T0] + T_lu + [self.T1]
            
            # Graficar
            plt.plot(x_nodes, T_full, 'o-', label=f'N = {N}', markersize=3)
        
        plt.xlabel('Posición x')
        plt.ylabel('Temperatura T(x)')
        plt.title('Distribución de temperatura - Conducción de calor 1D')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

def main():
    """Función principal para ejecutar el análisis"""
    
    # Crear instancia del problema
    heat_problem = HeatConduction1D()
    
    # Valores de N a analizar
    N_values = [50, 100, 500]
    
    print("=" * 60)
    print("ANÁLISIS DE CONDUCCIÓN DE CALOR 1D")
    print("=" * 60)
    
    # Resultados para cada N
    results = []
    
    for N in N_values:
        print(f"\n--- N = {N} ---")
        print(f"Paso espacial h = {1/(N+1):.6f}")
        
        # Solución con LU
        T_lu, x_nodes, time_lu, ops_est = heat_problem.solve_with_lu(N)
        
        # Solución de referencia con numpy
        T_ref = heat_problem.solve_with_numpy(N)
        
        # Calcular error
        error = heat_problem.calculate_error(T_lu, T_ref)
        
        # Almacenar resultados
        results.append({
            'N': N,
            'time_lu': time_lu,
            'operations': ops_est,
            'error': error,
            'T_lu': T_lu,
            'x_nodes': x_nodes
        })
        
        print(f"Tiempo LU: {time_lu:.6f} segundos")
        print(f"Operaciones estimadas: {ops_est}")
        print(f"Error ∥T_LU - T_ref∥∞: {error:.2e}")
    
    # Análisis de complejidad
    print("\n" + "=" * 60)
    print("ANÁLISIS DE COMPLEJIDAD")
    print("=" * 60)
    
    for i in range(len(results)):
        N = results[i]['N']
        time_lu = results[i]['time_lu']
        ops_est = results[i]['operations']
        
        print(f"N = {N}:")
        print(f"  Tiempo: {time_lu:.6f} s")
        print(f"  Operaciones: {ops_est}")
        print(f"  Tiempo/Operación: {time_lu/ops_est:.2e} s/op")
        
        if i > 0:
            prev_N = results[i-1]['N']
            prev_time = results[i-1]['time_lu']
            ratio = time_lu / prev_time if prev_time > 0 else 0
            expected_ratio = (N/prev_N)**3
            print(f"  Ratio tiempo: {ratio:.2f} (esperado: {expected_ratio:.2f})")
    
    # Graficar soluciones
    print("\nGenerando gráficas...")
    heat_problem.plot_solutions(N_values)
    
    # Graficar error vs N
    plt.figure(figsize=(10, 6))
    N_list = [r['N'] for r in results]
    errors = [r['error'] for r in results]
    
    plt.loglog(N_list, errors, 'bo-', linewidth=2, markersize=8)
    plt.xlabel('N (número de nodos interiores)')
    plt.ylabel('Error ∥T_LU - T_ref∥∞')
    plt.title('Error vs Resolución espacial')
    plt.grid(True, alpha=0.3)
    plt.show()
    
    # Graficar tiempo vs N
    plt.figure(figsize=(10, 6))
    times = [r['time_lu'] for r in results]
    
    plt.loglog(N_list, times, 'ro-', linewidth=2, markersize=8, label='Tiempo medido')
    
    # Añadir línea de referencia O(N³)
    ref_N = np.array(N_list)
    ref_time = times[0] * (ref_N / N_list[0])**3
    plt.loglog(ref_N, ref_time, 'k--', label='O(N³)')
    
    plt.xlabel('N (número de nodos interiores)')
    plt.ylabel('Tiempo de ejecución (s)')
    plt.title('Tiempo de ejecución vs N')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

if __name__ == "__main__":
    main()