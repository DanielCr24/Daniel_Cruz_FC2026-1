from matrix import Matrix
import myFunctions as myF
from linear_systems import forward_substitution, backward_substitution
import numpy as np

#Implementamos la factorización LUP
def lu_factorization_pp(A: Matrix):
    """
    Factoriza PA = LU con pivoteo parcial (Doolittle con permutaciones)
    
    Parameters
    ----------
    A : Matrix
        Matriz cuadrada a factorizar
        
    Returns
    -------
    P : Matrix
        Matriz de permutación
    L : Matrix
        Matriz triangular inferior con diagonal unitaria
    U : Matrix
        Matriz triangular superior
        
    Raises
    ------
    ValueError
        Si la matriz no es cuadrada
    """
    if A.rows != A.cols:
        raise ValueError("La matriz debe ser cuadrada para la factorización LU")
    
    n = A.rows
    
    # Inicializar matrices
    U = A.copy()  # Trabajaremos sobre una copia
    L = Matrix([[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)])
    P = Matrix([[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)])
    
    # Vector para seguir los intercambios de filas
    pivot_history = []
    
    for k in range(n - 1):
        # Encontrar el pivote máximo en la columna k
        max_val = abs(U.data[k][k])
        pivot_row = k
        
        for i in range(k + 1, n):
            if abs(U.data[i][k]) > max_val:
                max_val = abs(U.data[i][k])
                pivot_row = i
        
        pivot_history.append((k, pivot_row))
        
        # Intercambiar filas si es necesario
        if pivot_row != k:
            # Intercambiar filas en U
            U.data[k], U.data[pivot_row] = U.data[pivot_row], U.data[k]
            
            # Intercambiar filas en P
            P.data[k], P.data[pivot_row] = P.data[pivot_row], P.data[k]
            
            # Intercambiar elementos en L (solo los ya calculados)
            for j in range(k):
                L.data[k][j], L.data[pivot_row][j] = L.data[pivot_row][j], L.data[k][j]
        
        # Verificar que el pivote no sea cero
        if abs(U.data[k][k]) < 1e-12:
            raise ZeroDivisionError("Matriz singular o casi singular")
        
        # Calcular multiplicadores y eliminar
        for i in range(k + 1, n):
            L.data[i][k] = U.data[i][k] / U.data[k][k]
            for j in range(k, n):
                U.data[i][j] -= L.data[i][k] * U.data[k][j]
    
    return P, L, U

def lu_solve_pp(A: Matrix, b: list) -> list:
    """
    Resuelve Ax = b usando factorización LU con pivoteo parcial
    
    Parameters
    ----------
    A : Matrix
        Matriz del sistema
    b : list
        Vector del lado derecho
        
    Returns
    -------
    list
        Solución x del sistema
    """
    # Factorización PA = LU
    P, L, U = lu_factorization_pp(A)
    
    # Aplicar permutación al vector b: Pb
    n = len(b)
    b_perm = [0.0] * n
    for i in range(n):
        b_perm[i] = sum(P.data[i][j] * b[j] for j in range(n))
    
    # Resolver Ly = Pb
    y = forward_substitution(L, b_perm)
    
    # Resolver Ux = y
    x = backward_substitution(U, y)
    
    return x

#Podemos usar esta función para verificar el resultado, pues verifica que Ax ≈ b
def verify_solution(A: Matrix, b: list, x: list, tol=1e-10) -> bool:
    """
    Verifica que Ax ≈ b
    
    Parameters
    ----------
    A : Matrix
        Matriz del sistema
    b : list
        Vector del lado derecho
    x : list
        Solución calculada
    tol : float
        Tolerancia para la verificación
        
    Returns
    -------
    bool
        True si la solución es correcta dentro de la tolerancia
    """
    n = len(b)
    residual = 0.0
    
    for i in range(n):
        ax_i = sum(A.data[i][j] * x[j] for j in range(n))
        residual += (ax_i - b[i]) ** 2
    
    return (residual ** 0.5) < tol

#Implementamos el algoritmo RREF
def rref(A: Matrix, tol=1e-12):
    """
    Reduce una matriz a su forma escalonada reducida por filas (RREF)
    con pivoteo parcial
    """
    A_rref = A.copy()
    rows, cols = A_rref.rows, A_rref.cols
    
    pivot_row = 0
    pivot_cols = []  # Columnas con pivote
    
    for col in range(cols):
        # Encontrar el pivote máximo en la columna actual
        max_row = pivot_row
        max_val = abs(A_rref.data[pivot_row][col])
        
        for row in range(pivot_row + 1, rows):
            if abs(A_rref.data[row][col]) > max_val:
                max_val = abs(A_rref.data[row][col])
                max_row = row
        
        # Si el elemento máximo es cero, pasar a la siguiente columna
        if max_val < tol:
            continue
        
        # Intercambiar filas si es necesario
        if max_row != pivot_row:
            A_rref.data[pivot_row], A_rref.data[max_row] = A_rref.data[max_row], A_rref.data[pivot_row]
        
        # Normalizar la fila del pivote
        pivot_val = A_rref.data[pivot_row][col]
        for j in range(col, cols):
            A_rref.data[pivot_row][j] /= pivot_val
        
        # Eliminar elementos arriba y abajo del pivote
        for i in range(rows):
            if i != pivot_row:
                factor = A_rref.data[i][col]
                for j in range(col, cols):
                    A_rref.data[i][j] -= factor * A_rref.data[pivot_row][j]
        
        pivot_cols.append(col)
        pivot_row += 1
        
        if pivot_row == rows:
            break
    
    return A_rref, pivot_cols

#Implementamos el espacio nulo de una matriz
def nullspace_basis(A: Matrix, tol=1e-12):
    """
    Encuentra una base para el espacio nulo de A
    """
    A_rref, pivot_cols = rref(A)
    rows, cols = A_rref.rows, A_rref.cols
    
    # Identificar columnas pivote y libres
    pivot_columns = set(pivot_cols)
    free_columns = [j for j in range(cols) if j not in pivot_columns]
    
    basis = []
    
    for free_col in free_columns:
        # Construir vector del espacio nulo
        null_vector = [0.0] * cols
        
        # La variable libre corresponde a 1
        null_vector[free_col] = 1.0
        
        # Resolver para las variables dependientes
        for i in range(min(rows, cols)):
            if i < len(pivot_cols):
                pivot_col = pivot_cols[i]
                if pivot_col < free_col:
                    # Calcular variable dependiente
                    suma = 0.0
                    for j in range(pivot_col + 1, cols):
                        suma += A_rref.data[i][j] * null_vector[j]
                    null_vector[pivot_col] = -suma
        
        basis.append(null_vector)
    
    return basis

#Implementamos la función qr_decomposition(A) que realice la composición A = QR utilizando el método de Gram-Schimdt modificado
def qr_decomposition(A: Matrix):
    """
    Descomposición QR usando Gram-Schmidt modificado (MGS).
    
    Parameters
    ----------
    A : Matrix
        Matriz de entrada m×n (m ≥ n)
    
    Returns
    -------
    Q : Matrix
        Matriz con columnas ortonormales
    R : Matrix
        Matriz triangular superior
    """
    m, n = A.shape()
    
    # Copias de las columnas de A (editables)
    Acols = [A.get_col(j)[:] for j in range(n)]
    
    # Inicializar Q y R
    Q_data = [[0.0] * n for _ in range(m)]
    R_data = [[0.0] * n for _ in range(n)]
    
    for j in range(n):
        v = Acols[j][:]  # Copia de la columna j
        
        # Ortogonalización respecto a todas las q_i anteriores
        for i in range(j):
            # r_ij = q_i^T v
            r_ij = sum(Q_data[k][i] * v[k] for k in range(m))
            R_data[i][j] = r_ij
            
            # v = v - r_ij * q_i
            for k in range(m):
                v[k] -= r_ij * Q_data[k][i]
        
        # r_jj = ||v||_2
        r_jj = myF.sqrt(sum(v[k]**2 for k in range(m)))
        R_data[j][j] = r_jj
        
        if abs(r_jj) < 1e-12:  # Columna casi dependiente
            # q_j = vector cero
            for k in range(m):
                Q_data[k][j] = 0.0
        else:
            # q_j = v / r_jj
            for k in range(m):
                Q_data[k][j] = v[k] / r_jj
    
    return Matrix(Q_data), Matrix(R_data)

#Implementamos una función householder_qr(A) que construya matrices Q y R mediante reflexiones ortogonales.
def householder_qr(A: Matrix):
    """
    Factorización QR usando transformaciones de Householder.
    
    Parameters
    ----------
    A : Matrix
        Matriz de entrada m×n (m ≥ n)
    
    Returns
    -------
    Q : Matrix
        Matriz ortogonal m×m
    R : Matrix
        Matriz triangular superior m×n
    """
    m, n = A.shape()
    
    # Copia de A que se convertirá en R
    R = A.copy()
    
    # Inicializar Q como matriz identidad m×m
    Q = Matrix([[1.0 if i == j else 0.0 for j in range(m)] for i in range(m)])
    
    # Lista para almacenar los vectores v de Householder
    householder_vectors = []
    
    for k in range(min(m, n)):
        # Vector x: columna k de R desde la fila k hasta el final
        x = [R.data[i][k] for i in range(k, m)]
        
        # Calcular norma del vector x
        norm_x = myF.sqrt(sum(xi**2 for xi in x))
        
        if norm_x < 1e-12:
            # Columna casi cero, saltar esta transformación
            householder_vectors.append(None)
            continue
        
        # Calcular vector v = x + sign(x0)*||x||*e1
        sign = 1.0 if x[0] >= 0 else -1.0
        v = [xi for xi in x]
        v[0] += sign * norm_x
        
        # Normalizar v
        norm_v = myF.sqrt(sum(vi**2 for vi in v))
        if norm_v < 1e-12:
            householder_vectors.append(None)
            continue
            
        v = [vi / norm_v for vi in v]
        householder_vectors.append(v)
        
        # Aplicar transformación de Householder a R (submatriz desde fila k, columna k)
        for j in range(k, n):
            # Calcular producto punto: v^T * columna j de R
            col_j = [R.data[i][j] for i in range(k, m)]
            dot_product = sum(v[r] * col_j[r] for r in range(len(v)))
            
            # Actualizar columna j: R_col = R_col - 2 * v * (v^T * R_col)
            for i in range(len(v)):
                R.data[k + i][j] -= 2 * v[i] * dot_product
        
        # Acumular transformación en Q: Q = Q * H_k
        # Aplicar H_k a Q (a todas las columnas de Q)
        for j in range(m):
            # Columna j de Q desde fila k
            q_col = [Q.data[i][j] for i in range(k, m)]
            dot_product = sum(v[r] * q_col[r] for r in range(len(v)))
            
            for i in range(len(v)):
                Q.data[k + i][j] -= 2 * v[i] * dot_product
    
    return Q, R