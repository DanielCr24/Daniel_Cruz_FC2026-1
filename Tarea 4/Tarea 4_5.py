import numpy as np
import matplotlib.pyplot as plt
from matrix import Matrix
from myLinearMethods import householder_qr
import myFunctions as myF

def least_squares_householder(A: Matrix, y: list):
    """
    Resuelve el problema de mínimos cuadrados A β ≈ y usando Householder QR
    
    Parameters
    ----------
    A : Matrix
        Matriz de mxn
    y : list
        Vector de observaciones
    
    Returns
    -------
    beta : list
        Vector de coeficientes β
    residual_norm : float
        Norma del residuo
    """
    m, n = A.shape()
    
    # Factorización QR de A usando Householder
    Q, R = householder_qr(A)
    
    # Calcular Q^T y
    Q_T = Q.transpose()
    QTy = [0.0] * m
    for i in range(m):
        QTy[i] = sum(Q_T.data[i][j] * y[j] for j in range(m))
    
    # Extraer la parte triangular superior de R y el vector correspondiente de Q^T y
    R_upper = Matrix([[R.data[i][j] for j in range(n)] for i in range(n)])
    QTy_upper = QTy[:n]
    
    # Resolver R β = Q^T y usando sustitución hacia atrás
    beta = [0.0] * n
    for i in range(n-1, -1, -1):
        beta[i] = QTy_upper[i]
        for j in range(i+1, n):
            beta[i] -= R_upper.data[i][j] * beta[j]
        beta[i] /= R_upper.data[i][i]
    
    # Calcular norma del residuo
    residual = [0.0] * m
    for i in range(m):
        Ax_i = sum(A.data[i][j] * beta[j] for j in range(n))
        residual[i] = y[i] - Ax_i
    
    residual_norm = myF.sqrt(sum(r**2 for r in residual))
    
    return beta, residual_norm

def generate_parabolic_data(v0=20.0, theta_deg=45.0, g=9.81, sigma=0.5, num_points=20):
    """
    Genera datos simulados de movimiento parabólico con ruido
    
    Parameters
    ----------
    v0 : float
        Velocidad inicial (m/s)
    theta_deg : float
        Ángulo de lanzamiento (grados)
    g : float
        Gravedad (m/s²)
    sigma : float
        Desviación estándar del ruido
    num_points : int
        Número de puntos de datos
    
    Returns
    -------
    x_data : list
        Posiciones horizontales
    y_data : list
        Posiciones verticales con ruido
    y_theoretical : list
        Posiciones verticales teóricas sin ruido
    """
    theta_rad = np.radians(theta_deg)
    
    # Rango máximo teórico
    x_max = (v0**2 * myF.Sin(2*theta_rad)) / g
    x_data = np.linspace(0, x_max * 0.9, num_points).tolist()
    
    # Ecuación teórica del movimiento parabólico
    y_theoretical = []
    for x in x_data:
        y = x * np.tan(theta_rad) - (g * x**2) / (2 * v0**2 * myF.Cos(theta_rad)**2)
        y_theoretical.append(y)
    
    # Añadir ruido gaussiano
    y_data = [y + np.random.normal(0, sigma) for y in y_theoretical]
    
    return x_data, y_data, y_theoretical

def create_design_matrix(x_data, degree=2):
    """
    Crea la matriz de diseño para ajuste polinomial
    
    Parameters
    ----------
    x_data : list
        Datos de entrada
    degree : int
        Grado del polinomio
    
    Returns
    -------
    A : Matrix
        Matriz de diseño
    """
    m = len(x_data)
    n = degree + 1
    
    A_data = []
    for i in range(m):
        row = []
        for j in range(n):
            row.append(x_data[i]**j)
        A_data.append(row)
    
    return Matrix(A_data)

def parabolic_trajectory(x, v0=20.0, theta_deg=45.0, g=9.81):
    """
    Calcula la trayectoria parabólica teórica
    """
    theta_rad = np.radians(theta_deg)
    return x * np.tan(theta_rad) - (g * x**2) / (2 * v0**2 * myF.Cos(theta_rad)**2)

def verify_orthogonality(Q: Matrix):
    """
    Verifica que QᵀQ ≈ I
    """
    n = Q.rows
    
    Q_T = Q.transpose()
    Q_T_Q = Q_T * Q
    
    # Calcular residuo ∥QᵀQ - I∥
    identity_diff = Matrix([[Q_T_Q.data[i][j] - (1.0 if i == j else 0.0) 
                          for j in range(n)] for i in range(n)])
    
    residual_norm = 0.0
    for i in range(n):
        for j in range(n):
            residual_norm += identity_diff.data[i][j] ** 2
    residual_norm = myF.sqrt(residual_norm)
    
    return residual_norm, Q_T_Q

def main():
    """
    Función principal para el ajuste por mínimos cuadrados
    """
    print("=" * 70)
    print("AJUSTE POR MÍNIMOS CUADRADOS - MOVIMIENTO PARABÓLICO")
    print("=" * 70)
    
    # Parámetros del problema
    v0 = 20.0  # m/s
    theta_deg = 45.0  # grados
    g = 9.81  # m/s²
    sigma = 0.5  # desviación estándar del ruido
    
    # c) Generar datos simulados
    print("\nGenerando datos simulados...")
    print(f"Parámetros: v0 = {v0} m/s, θ = {theta_deg}°, g = {g} m/s², σ = {sigma} m")
    
    x_data, y_data, y_theoretical = generate_parabolic_data(v0, theta_deg, g, sigma)
    
    print(f"Número de puntos: {len(x_data)}")
    print(f"Rango de x: {min(x_data):.2f} a {max(x_data):.2f} m")
    print(f"Rango de y: {min(y_data):.2f} a {max(y_data):.2f} m")
    
    # d) Crear matriz de diseño para ajuste cuadrático
    print("\nCreando matriz de diseño para ajuste cuadrático...")
    A = create_design_matrix(x_data, degree=2)
    
    print(f"Dimensión de A: {A.rows} × {A.cols}")
    print("Primeras filas de la matriz de diseño A:")
    for i in range(min(3, A.rows)):
        print(f"  {[f'{val:.3f}' for val in A.data[i]]}")
    
    # Verificar ortogonalidad de Q
    print("\nVerificando factorización QR...")
    Q, R = householder_qr(A)
    ortho_error, Q_T_Q = verify_orthogonality(Q)
    print(f"Error de ortogonalidad ∥QᵀQ - I∥ = {ortho_error:.2e}")
    
    # Resolver usando Householder QR
    print("\nResolviendo el problema de mínimos cuadrados...")
    beta, residual_norm = least_squares_householder(A, y_data)
    
    print(f"\nCoeficientes del ajuste cuadrático:")
    print(f"β₀ (término constante) = {beta[0]:.6f}")
    print(f"β₁ (coeficiente lineal) = {beta[1]:.6f}")
    print(f"β₂ (coeficiente cuadrático) = {beta[2]:.6f}")
    
    print(f"Norma del residuo = {residual_norm:.6f}")
    
    # Coeficientes teóricos esperados
    theta_rad = np.radians(theta_deg)
    beta0_theoretical = 0.0  # y(0) = 0
    beta1_theoretical = np.tan(theta_rad)  # dy/dx en x=0
    beta2_theoretical = -g / (2 * v0**2 * np.cos(theta_rad)**2)  # coeficiente de x²
    
    print(f"\nCoeficientes teóricos:")
    print(f"β₀ teórico = {beta0_theoretical:.6f}")
    print(f"β₁ teórico = {beta1_theoretical:.6f}")
    print(f"β₂ teórico = {beta2_theoretical:.6f}")
    
    # Errores relativos
    error_beta0 = abs(beta[0] - beta0_theoretical) / (abs(beta0_theoretical) + 1e-12)
    error_beta1 = abs(beta[1] - beta1_theoretical) / abs(beta1_theoretical)
    error_beta2 = abs(beta[2] - beta2_theoretical) / abs(beta2_theoretical)
    
    print(f"\nErrores relativos:")
    print(f"Error en β₀ = {error_beta0:.2e}")
    print(f"Error en β₁ = {error_beta1:.2e}")
    print(f"Error en β₂ = {error_beta2:.2e}")
    
    # e) Graficar resultados
    print("\nGenerando gráficas...")
    
    # Crear puntos para la curva ajustada
    x_fine = np.linspace(min(x_data), max(x_data), 100)
    y_fitted = [beta[0] + beta[1]*x + beta[2]*x**2 for x in x_fine]
    y_theoretical_fine = [parabolic_trajectory(x, v0, theta_deg, g) for x in x_fine]
    
    plt.figure(figsize=(12, 8))
    
    # Datos experimentales con ruido
    plt.scatter(x_data, y_data, color='red', s=50, alpha=0.7, 
                label='Datos experimentales con ruido', zorder=5)
    
    # Curva ajustada
    plt.plot(x_fine, y_fitted, 'b-', linewidth=2, 
             label=f'Ajuste cuadrático: y = {beta[0]:.3f} + {beta[1]:.3f}x + {beta[2]:.3f}x²')
    
    # Trayectoria teórica sin ruido
    plt.plot(x_fine, y_theoretical_fine, 'g--', linewidth=2, 
             label='Trayectoria teórica sin ruido')
    
    plt.xlabel('Distancia horizontal x (m)', fontsize=12)
    plt.ylabel('Altura y (m)', fontsize=12)
    plt.title('Ajuste por Mínimos Cuadrados - Movimiento Parabólico\n' +
              f'(Householder QR, v₀ = {v0} m/s, θ = {theta_deg}°, σ = {sigma} m)', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    
    # Añadir información del ajuste en el gráfico
    textstr = '\n'.join([
        f'Coeficientes del ajuste:',
        f'β₀ = {beta[0]:.4f}',
        f'β₁ = {beta[1]:.4f}',
        f'β₂ = {beta[2]:.4f}',
        f'Norma del residuo = {residual_norm:.4f}'
    ])
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.show()
    
    # Gráfica de residuos
    plt.figure(figsize=(12, 6))
    
    # Calcular residuos
    residuals = []
    for i in range(len(x_data)):
        y_pred = beta[0] + beta[1]*x_data[i] + beta[2]*x_data[i]**2
        residuals.append(y_data[i] - y_pred)
    
    plt.scatter(x_data, residuals, color='purple', s=60, alpha=0.7)
    plt.axhline(y=0, color='black', linestyle='-', alpha=0.8)
    plt.xlabel('Distancia horizontal x (m)', fontsize=12)
    plt.ylabel('Residuos (m)', fontsize=12)
    plt.title('Análisis de Residuos - Ajuste Cuadrático', fontsize=14)
    plt.grid(True, alpha=0.3)
    
    # Añadir estadísticas de residuos
    mean_residual = np.mean(residuals)
    std_residual = np.std(residuals)
    
    textstr = '\n'.join([
        f'Estadísticas de residuos:',
        f'Media = {mean_residual:.4f} m',
        f'Desv. estándar = {std_residual:.4f} m'
    ])
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.show()
        
    if ortho_error < 1e-10:
        print("Excelente ortogonalidad de Q")
    elif ortho_error < 1e-6:
        print("Buena ortogonalidad de Q")
    else:
        print("Ortogonalidad de Q podría mejorarse")

if __name__ == "__main__":
    # Configurar semilla para reproducibilidad
    np.random.seed(42)
    main()