import numpy as np
import time
import matplotlib.pyplot as plt
from matrix import Matrix
from myLinearMethods import qr_decomposition
import myFunctions as myF

class QREigenSolver:
    """
    Resuelve el problema de autovalores usando el método QR
    """
    
    def __init__(self, A):
        self.A = A
        self.n = A.rows
        
    def qr_algorithm_no_shift(self, max_iter=50, tol=1e-12):
        """
        Aplica el método QR sin desplazamiento
        
        Parameters
        ----------
        max_iter : int
            Número máximo de iteraciones
        tol : float
            Tolerancia para convergencia
            
        Returns
        -------
        eigenvalues : list
            Autovalores aproximados
        eigenvectors : Matrix
            Matriz de autovectores aproximados
        convergence_history : list
            Historia de la convergencia
        """
        print("=" * 70)
        print("MÉTODO QR PARA AUTOVALORES - MATRIZ LAPLACIANA 1D")
        print("=" * 70)
        
        # Matriz inicial
        Ak = self.A.copy()
        Q_total = Matrix([[1.0 if i == j else 0.0 for j in range(self.n)] 
                         for i in range(self.n)])
        
        convergence_history = []
        
        print(f"\nMatriz inicial A (3×3 Laplaciana discreta):")
        self.print_matrix(Ak)
        
        for k in range(max_iter):
            # Factorización QR usando Gram-Schmidt modificado
            Qk, Rk = qr_decomposition(Ak)
            
            # Reconstrucción: A_{k+1} = R_k Q_k
            Ak_plus_1 = Rk * Qk
            
            # Acumular transformaciones ortogonales
            Q_total = Q_total * Qk
            
            # Almacenar elementos diagonales para seguimiento
            diag_elements = [Ak_plus_1.data[i][i] for i in range(self.n)]
            convergence_history.append(diag_elements)
            
            # Verificar convergencia (matriz triangular superior)
            off_diag_norm = self.off_diagonal_norm(Ak_plus_1)
            
            if k < 3 or k % 10 == 0 or off_diag_norm < tol:
                print(f"\nIteración {k+1}:")
                print(f"Elementos diagonales: {[f'{x:.8f}' for x in diag_elements]}")
                print(f"Norma elementos fuera de diagonal: {off_diag_norm:.2e}")
            
            # Actualizar para siguiente iteración
            Ak = Ak_plus_1
            
            if off_diag_norm < tol:
                print(f"\nConvergencia alcanzada en {k+1} iteraciones")
                break
        
        # Autovalores aproximados (elementos diagonales)
        eigenvalues = [Ak.data[i][i] for i in range(self.n)]
        
        return eigenvalues, Q_total, convergence_history
    
    def off_diagonal_norm(self, A):
        """
        Calcula la norma de los elementos fuera de la diagonal
        """
        norm = 0.0
        for i in range(A.rows):
            for j in range(A.cols):
                if i != j:
                    norm += A.data[i][j] ** 2
        return myF.sqrt(norm)
    
    def print_matrix(self, A, precision=6):
        """
        Imprime una matriz de forma legible
        """
        for i in range(A.rows):
            row_str = "  ["
            for j in range(A.cols):
                row_str += f"{A.data[i][j]:.{precision}f}"
                if j < A.cols - 1:
                    row_str += ", "
            row_str += "]"
            print(row_str)
    
    def compare_with_numpy(self, qr_eigenvalues, qr_eigenvectors):
        """
        Compara resultados con numpy.linalg.eig
        """
        print("\n" + "=" * 70)
        print("COMPARACIÓN CON NUMPY.LINALG.EIG")
        print("=" * 70)
        
        # Convertir a numpy para comparación
        A_np = np.array(self.A.data)
        
        # Calcular autovalores exactos con numpy
        np_eigenvalues, np_eigenvectors = np.linalg.eig(A_np)
        
        # Ordenar autovalores (ambos en orden descendente para comparar)
        qr_sorted = sorted(qr_eigenvalues, reverse=True)
        np_sorted = sorted(np_eigenvalues, reverse=True)
        
        print(f"\nAutovalores QR:    {[f'{x:.8f}' for x in qr_sorted]}")
        print(f"Autovalores NumPy: {[f'{x:.8f}' for x in np_sorted]}")
        
        # Calcular errores relativos
        errors = []
        for i in range(self.n):
            error = abs(qr_sorted[i] - np_sorted[i]) / abs(np_sorted[i])
            errors.append(error)
            print(f"λ{i+1}: Error relativo = {error:.2e}")
        
        return np_sorted, np_eigenvectors, errors
    
    def verify_orthogonality(self, Q):
        """
        Verifica que QᵀQ ≈ I
        """
        print("\n" + "=" * 70)
        print("VERIFICACIÓN DE ORTOGONALIDAD QᵀQ ≈ I")
        print("=" * 70)
        
        Q_T = Q.transpose()
        Q_T_Q = Q_T * Q
        
        print(f"\nMatriz QᵀQ:")
        self.print_matrix(Q_T_Q)
        
        # Calcular residuo ∥QᵀQ - I∥
        identity_diff = Matrix([[Q_T_Q.data[i][j] - (1.0 if i == j else 0.0) 
                              for j in range(self.n)] for i in range(self.n)])
        
        residual_norm = 0.0
        for i in range(self.n):
            for j in range(self.n):
                residual_norm += identity_diff.data[i][j] ** 2
        residual_norm = myF.sqrt(residual_norm)
        
        print(f"\n∥QᵀQ - I∥₂ = {residual_norm:.2e}")
        
        # Verificar elementos diagonales ≈ 1 y fuera de diagonal ≈ 0
        max_off_diag = 0.0
        min_diag = float('inf')
        max_diag = 0.0
        
        for i in range(self.n):
            for j in range(self.n):
                if i == j:
                    min_diag = min(min_diag, abs(Q_T_Q.data[i][j]))
                    max_diag = max(max_diag, abs(Q_T_Q.data[i][j]))
                else:
                    max_off_diag = max(max_off_diag, abs(Q_T_Q.data[i][j]))
        
        print(f"Elementos diagonales: min = {min_diag:.2e}, max = {max_diag:.2e}")
        print(f"Máximo elemento fuera de diagonal: {max_off_diag:.2e}")
        
        return residual_norm
    
    def analyze_convergence(self, convergence_history, np_eigenvalues):
        """
        Analiza la convergencia de los autovalores
        """
        print("\n" + "=" * 70)
        print("ANÁLISIS DE CONVERGENCIA")
        print("=" * 70)
        
        np_sorted = sorted(np_eigenvalues, reverse=True)
        
        # Calcular errores en cada iteración
        errors_history = []
        print("\nEvolución de errores relativos:")
        print("Iter | Error λ₁ | Error λ₂ | Error λ₃")
        print("-" * 40)
        
        for iter_idx, diag_elems in enumerate(convergence_history):
            diag_sorted = sorted(diag_elems, reverse=True)
            iter_errors = [abs(diag_sorted[i] - np_sorted[i]) / abs(np_sorted[i]) 
                          for i in range(self.n)]
            errors_history.append(iter_errors)
            
            if iter_idx < 5 or iter_idx % 10 == 9 or iter_idx == len(convergence_history)-1:
                error_str = " | ".join([f"{e:.2e}" for e in iter_errors])
                print(f"{iter_idx+1:4d} | {error_str}")
        
        return errors_history

def create_laplacian_matrix(n=3):
    """
    Crea la matriz Laplaciana discreta 1D de tamaño n×n
    
    A = [[2, -1, 0, ..., 0],
         [-1, 2, -1, ..., 0],
         [0, -1, 2, ..., 0],
         ...
         [0, 0, 0, ..., 2]]
    """
    data = [[0.0] * n for _ in range(n)]
    
    for i in range(n):
        data[i][i] = 2.0
        if i > 0:
            data[i][i-1] = -1.0
        if i < n-1:
            data[i][i+1] = -1.0
    
    return Matrix(data)

def plot_convergence(convergence_history, np_eigenvalues):
    """
    Grafica la convergencia de los autovalores
    """
    plt.figure(figsize=(12, 8))
    
    np_sorted = sorted(np_eigenvalues, reverse=True)
    iterations = range(1, len(convergence_history) + 1)
    
    # Para cada autovalor
    for i in range(3):
        lambda_errors = []
        for iter_idx, diag_elems in enumerate(convergence_history):
            diag_sorted = sorted(diag_elems, reverse=True)
            error = abs(diag_sorted[i] - np_sorted[i]) / abs(np_sorted[i])
            lambda_errors.append(error)
        
        plt.semilogy(iterations, lambda_errors, 'o-', 
                    label=f'λ{i+1} = {np_sorted[i]:.3f}', markersize=4)
    
    plt.xlabel('Iteración')
    plt.ylabel('Error relativo')
    plt.title('Convergencia de Autovalores - Método QR')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

def main():
    """Función principal"""
    
    # b) Crear matriz Laplaciana discreta 1D de 3×3
    A = create_laplacian_matrix(3)
    
    print("MATRIZ LAPLACIANA DISCRETA 1D (3×3)")
    print("Representa vibraciones de cadena de 3 masas acopladas")
    print("A = [[2, -1, 0],")
    print("     [-1, 2, -1],") 
    print("     [0, -1, 2]]")
    
    # Crear solver QR
    solver = QREigenSolver(A)
    
    # c) Aplicar método QR iterativo (50 iteraciones)
    start_time = time.time()
    qr_eigenvalues, qr_eigenvectors, convergence_history = solver.qr_algorithm_no_shift(max_iter=50)
    qr_time = time.time() - start_time
    
    # d) Comparar con numpy
    np_eigenvalues, np_eigenvectors, errors = solver.compare_with_numpy(
        qr_eigenvalues, qr_eigenvectors)
    
    # e) Verificar ortogonalidad
    ortho_error = solver.verify_orthogonality(qr_eigenvectors)
    
    # Análisis de convergencia
    errors_history = solver.analyze_convergence(convergence_history, np_eigenvalues)
    
    # Resultados teóricos esperados
    print("\n" + "=" * 70)
    print("RESULTADOS TEÓRICOS ESPERADOS")
    print("=" * 70)
    
    # Para matriz Laplaciana 1D de tamaño n, los autovalores son:
    # λ_k = 2 - 2cos(kπ/(n+1)), k = 1,2,...,n
    n = 3
    expected_eigenvalues = []
    for k in range(1, n+1):
        lambda_k = 2 - 2 * myF.Cos(k * np.pi / (n + 1))
        expected_eigenvalues.append(lambda_k)
    
    expected_eigenvalues_sorted = sorted(expected_eigenvalues, reverse=True)
    
    print(f"Autovalores teóricos: {[f'{x:.8f}' for x in expected_eigenvalues_sorted]}")
    print(f"Autovalores NumPy:    {[f'{x:.8f}' for x in sorted(np_eigenvalues, reverse=True)]}")
    
    # Gráfica de convergencia
    print("\nGenerando gráfica de convergencia...")
    plot_convergence(convergence_history, np_eigenvalues)
    
    # Interpretación física
    print("\n" + "=" * 70)
    print("INTERPRETACIÓN FÍSICA")
    print("=" * 70)
    
    print("La matriz representa el operador Laplaciano discreto para")
    print("una cadena de 3 masas idénticas acopladas por resortes idénticos.")
    print("\nLos autovalores representan frecuencias al cuadrado:")
    print(f"λ₁ ≈ {expected_eigenvalues_sorted[0]:.3f}: Modo de alta frecuencia (masas en contrafase)")
    print(f"λ₂ ≈ {expected_eigenvalues_sorted[1]:.3f}: Modo de frecuencia intermedia") 
    print(f"λ₃ ≈ {expected_eigenvalues_sorted[2]:.3f}: Modo de baja frecuencia (masas en fase)")
    
    if max(errors) < 1e-6 and ortho_error < 1e-6:
        print("RESULTADOS SATISFACTORIOS: Alta precisión alcanzada")
    elif max(errors) < 1e-3 and ortho_error < 1e-3:
        print("RESULTADOS ACEPTABLES: Precisión moderada")
    else:
        print("PRECISIÓN LIMITADA: Considerar más iteraciones o método con desplazamiento")

if __name__ == "__main__":
    main()