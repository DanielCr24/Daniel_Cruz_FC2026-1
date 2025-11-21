"""
triangular.py
-------------

Este módulo define la clase TriangularMatrix que extiende Matrix para representar
matrices triangulares superiores e inferiores de forma eficiente.

Almacena solo los elementos no cero del triángulo correspondiente en una lista
unidimensional por filas.
"""

from matrix import Matrix
import math

class TriangularMatrix(Matrix):
    """
    Clase para representar matrices triangulares de forma eficiente.

    Attributes
    ----------
    n : int
        Tamaño de la matriz (n x n).
    kind : str
        Tipo de triángulo: "upper" o "lower".
    store : list
        Lista unidimensional que almacena los elementos del triángulo por filas.
    """

    def __init__(self, n, kind="upper", data=None):
        """
        Inicializa una matriz triangular.

        Parameters
        ----------
        n : int
            Tamaño de la matriz (n x n).
        kind : str, optional
            Tipo de triángulo: "upper" (superior) o "lower" (inferior).
        data : list, optional
            Lista con los elementos del triángulo en orden por filas.
            Si es None, se inicializa con ceros.
        """
        if kind not in ["upper", "lower"]:
            raise ValueError("El tipo debe ser 'upper' o 'lower'")
        
        self.n = n
        self.kind = kind
        self.rows = n
        self.cols = n
        
        # Número de elementos a almacenar: n(n+1)/2
        total_elements = n * (n + 1) // 2
        
        if data is None:
            self.store = [0.0] * total_elements
        else:
            if len(data) != total_elements:
                raise ValueError(f"Se esperaban {total_elements} elementos, pero se recibieron {len(data)}")
            self.store = list(data)
        
        # Inicializamos data para compatibilidad con Matrix
        # pero no la usaremos para almacenamiento principal
        dense_data = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if self._is_stored(i, j):
                    dense_data[i][j] = self._get_from_store(i, j)
        super().__init__(dense_data)

    def _is_stored(self, i, j):
        """Determina si el elemento (i,j) está almacenado."""
        if self.kind == "upper":
            return i <= j
        else:  # lower
            return i >= j

    def _get_index(self, i, j):
        """
        Calcula el índice en store para el elemento (i,j).
        
        Almacenamiento por filas:
        - Upper: fila 0: n elementos, fila 1: n-1 elementos, etc.
        - Lower: fila 0: 1 elemento, fila 1: 2 elementos, etc.
        """
        if not self._is_stored(i, j):
            raise IndexError(f"Elemento ({i},{j}) no está almacenado")
        
        if self.kind == "upper":
            # Para upper triangular: índice = i*n - i(i-1)/2 + (j-i)
            return i * self.n - i * (i - 1) // 2 + (j - i)
        else:  # lower
            # Para lower triangular: índice = i(i+1)/2 + j
            return i * (i + 1) // 2 + j

    def _get_from_store(self, i, j):
        """Obtiene el elemento (i,j) desde store."""
        if not self._is_stored(i, j):
            return 0.0
        return self.store[self._get_index(i, j)]

    def get(self, i, j):
        """Devuelve el elemento (i,j)."""
        if i < 0 or i >= self.n or j < 0 or j >= self.n:
            raise IndexError("Índices fuera de rango")
        return self._get_from_store(i, j)

    def set(self, i, j, val):
        """Establece el elemento (i,j). Lanza excepción si está fuera del triángulo."""
        if i < 0 or i >= self.n or j < 0 or j >= self.n:
            raise IndexError("Índices fuera de rango")
        
        if not self._is_stored(i, j):
            raise ValueError(f"No se puede modificar elemento ({i},{j}) fuera del triángulo {self.kind}")
        
        self.store[self._get_index(i, j)] = val
        # Actualizar también data para compatibilidad
        self.data[i][j] = val

    def toDense(self):
        """Convierte la matriz triangular a formato denso."""
        dense_data = [[0.0] * self.n for _ in range(self.n)]
        for i in range(self.n):
            for j in range(self.n):
                dense_data[i][j] = self.get(i, j)
        return Matrix(dense_data)

    def transpose(self):
        """Devuelve la transpuesta, intercambiando upper con lower."""
        # Crear nueva matriz triangular del tipo opuesto
        new_kind = "lower" if self.kind == "upper" else "upper"
        result = TriangularMatrix(self.n, new_kind)
        
        # Copiar elementos transpuestos
        for i in range(self.n):
            for j in range(self.n):
                if self._is_stored(i, j):
                    result.set(j, i, self.get(i, j))
        
        return result

    def matvec(self, x):
        """
        Producto matriz-vector optimizado.
        
        Parameters
        ----------
        x : list
            Vector de tamaño n.
        
        Returns
        -------
        list
            Resultado del producto A*x.
        """
        if len(x) != self.n:
            raise ValueError("El vector debe tener tamaño n")
        
        result = [0.0] * self.n
        
        if self.kind == "upper":
            for i in range(self.n):
                for j in range(i, self.n):
                    result[i] += self.get(i, j) * x[j]
        else:  # lower
            for i in range(self.n):
                for j in range(0, i + 1):
                    result[i] += self.get(i, j) * x[j]
        
        return result

    def __add__(self, other):
        """Suma de matrices triangulares o con Matrix densa."""
        if isinstance(other, TriangularMatrix):
            if self.n != other.n:
                raise ValueError("Las matrices deben tener el mismo tamaño")
            if self.kind != other.kind:
                # Convertir a denso si los tipos son diferentes
                return self.toDense() + other.toDense()
            
            # Suma eficiente de triangulares del mismo tipo
            result = TriangularMatrix(self.n, self.kind)
            for i in range(len(self.store)):
                result.store[i] = self.store[i] + other.store[i]
            return result
        
        elif isinstance(other, Matrix):
            # Convertir a denso y delegar en la clase base
            return self.toDense() + other
        
        else:
            raise TypeError("Operación no soportada")

    def __mul__(self, other):
        """Multiplicación por escalar."""
        if isinstance(other, (int, float)):
            result = TriangularMatrix(self.n, self.kind)
            for i in range(len(self.store)):
                result.store[i] = self.store[i] * other
            return result
        else:
            raise TypeError("Solo se permite multiplicación por escalar")

    def __eq__(self, other, tol=1e-10):
        """Comparación con tolerancia numérica."""
        if not isinstance(other, TriangularMatrix):
            return False
        
        if self.n != other.n or self.kind != other.kind:
            return False
        
        for i in range(len(self.store)):
            if abs(self.store[i] - other.store[i]) > tol:
                return False
        
        return True

    def __str__(self):
        """Representación para impresión."""
        return self.toDense().__str__()

    def shape(self):
        """Devuelve la dimensión de la matriz."""
        return (self.n, self.n)

    def get_row(self, i):
        """Devuelve la fila i."""
        return [self.get(i, j) for j in range(self.n)]

    def get_col(self, j):
        """Devuelve la columna j."""
        return [self.get(i, j) for i in range(self.n)]


# ======================
# Pruebas unitarias
# ======================
if __name__ == "__main__":
    import time
    
    print("=== PRUEBAS DE TRIANGULARMATRIX ===\n")
    
    # Prueba 1: Construcción y acceso
    print("1. Construcción y acceso:")
    n = 3
    upper_data = [1, 2, 3, 4, 5, 6]  # [[1,2,3],[0,4,5],[0,0,6]]
    upper = TriangularMatrix(n, "upper", upper_data)
    print("Upper matrix:")
    print(upper)
    
    lower_data = [1, 2, 3, 4, 5, 6]  # [[1,0,0],[2,3,0],[4,5,6]]
    lower = TriangularMatrix(n, "lower", lower_data)
    print("\nLower matrix:")
    print(lower)
    
    # Prueba 2: Suma de triangulares del mismo tipo
    print("\n2. Suma upper + upper:")
    upper2 = TriangularMatrix(n, "upper", [2, 3, 4, 5, 6, 7])
    result_upper = upper + upper2
    print(result_upper)
    
    print("\nSuma lower + lower:")
    lower2 = TriangularMatrix(n, "lower", [2, 3, 4, 5, 6, 7])
    result_lower = lower + lower2
    print(result_lower)
    
    # Prueba 3: Producto matriz-vector
    print("\n3. Producto matriz-vector:")
    x = [1, 2, 3]
    result_upper_vec = upper.matvec(x)
    result_lower_vec = lower.matvec(x)
    print(f"Upper * {x} = {result_upper_vec}")
    print(f"Lower * {x} = {result_lower_vec}")
    
    # Verificación contra matriz densa
    upper_dense = upper.toDense()
    expected_upper = [sum(upper_dense.data[i][j] * x[j] for j in range(n)) for i in range(n)]
    print(f"Verificación upper: {result_upper_vec == expected_upper}")
    
    # Prueba 4: Transpuesta
    print("\n4. Transpuesta:")
    upper_t = upper.transpose()
    print("Transpuesta de upper (debería ser lower):")
    print(upper_t)
    
    upper_tt = upper_t.transpose()
    print("\nTranspuesta de la transpuesta (debería ser original):")
    print(upper_tt)
    print(f"¿Transpuesta doble igual a original? {upper == upper_tt}")
    
    # Prueba 5: Multiplicación por escalar
    print("\n5. Multiplicación por escalar:")
    upper_scaled = upper * 2
    print(upper_scaled)
    
    # Prueba 6: Comparación
    print("\n6. Comparación:")
    upper_copy = TriangularMatrix(n, "upper", upper_data)
    print(f"¿Upper igual a copia? {upper == upper_copy}")
    
    # Prueba 7: Escalabilidad
    print("\n7. Prueba de escalabilidad:")
    sizes = [256, 512, 1024]
    
    for size in sizes:
        print(f"\nTamaño: {size}x{size}")
        
        # Crear matrices triangulares grandes
        start = time.time()
        upper_large = TriangularMatrix(size, "upper")
        for i in range(size):
            for j in range(i, size):
                upper_large.set(i, j, i + j)
        
        upper_large2 = TriangularMatrix(size, "upper")
        for i in range(size):
            for j in range(i, size):
                upper_large2.set(i, j, i * j)
        
        # Medir tiempo de suma
        start_time = time.time()
        result = upper_large + upper_large2
        end_time = time.time()
        
        print(f"Tiempo de suma triangular: {end_time - start_time:.4f} segundos")
        
        # Comparar con suma densa
        dense1 = upper_large.toDense()
        dense2 = upper_large2.toDense()
        
        start_time = time.time()
        dense_result = dense1 + dense2
        end_time = time.time()
        
        print(f"Tiempo de suma densa: {end_time - start_time:.4f} segundos")