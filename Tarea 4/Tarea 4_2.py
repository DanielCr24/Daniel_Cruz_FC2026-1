import numpy as np
from matrix import Matrix
from myLinearMethods import rref, nullspace_basis
import myFunctions as myF

class TwoMassVibration:
    """
    Sistema de dos masas acopladas por resortes
    m1, m2: masas
    k1, k2: resortes laterales  
    kc: resorte de acoplamiento
    """
    
    def __init__(self, m1=1, m2=2, k1=3, k2=4, kc=1):
        self.m1 = m1
        self.m2 = m2
        self.k1 = k1
        self.k2 = k2
        self.kc = kc
        
        # Matrices de masa y rigidez
        self.M = Matrix([[m1, 0], [0, m2]])
        self.K = Matrix([[k1 + kc, -kc], [-kc, k2 + kc]])
    
    def solve_eigenvalues_determinant(self):
        """
        a) Calcula autovalores resolviendo det(K - ω²M) = 0
        """
        print("=" * 60)
        print("CÁLCULO DE AUTOVALORES POR DETERMINANTE")
        print("=" * 60)
        
        # Para matriz 2x2: det([[a-ω²m1, b], [c, d-ω²m2]]) = 0
        # (a-ω²m1)(d-ω²m2) - bc = 0
        a, b, c, d = self.K.data[0][0], self.K.data[0][1], self.K.data[1][0], self.K.data[1][1]
        m1, m2 = self.M.data[0][0], self.M.data[1][1]
        
        print(f"Matriz K = [[{a}, {b}], [{c}, {d}]]")
        print(f"Matriz M = diag([{m1}, {m2}])")
        
        # Coeficientes de la ecuación cuadrática: A(ω²)² + B(ω²) + C = 0
        A = m1 * m2
        B = -(a * m2 + d * m1)
        C = a * d - b * c
        
        print(f"\nEcuación característica: {A}(ω²)² + {B}(ω²) + {C} = 0")
        
        # Resolver ecuación cuadrática
        discriminante = B**2 - 4*A*C
        
        if discriminante < 0:
            raise ValueError("Discriminante negativo, autovalores complejos")
        
        omega2_1 = (-B + np.sqrt(discriminante)) / (2*A)
        omega2_2 = (-B - np.sqrt(discriminante)) / (2*A)
        
        # Ordenar de menor a mayor
        omega2_sorted = sorted([omega2_1, omega2_2])
        omega_sorted = [np.sqrt(w2) for w2 in omega2_sorted]
        
        print(f"\nDiscriminante: {discriminante}")
        print(f"ω₁² = {omega2_sorted[0]:.6f}, ω₁ = {omega_sorted[0]:.6f}")
        print(f"ω₂² = {omega2_sorted[1]:.6f}, ω₂ = {omega_sorted[1]:.6f}")
        
        return omega2_sorted, omega_sorted
    
    def solve_eigenvectors_gauss_jordan(self, omega2):
        """
        b) Calcula autovectores usando eliminación Gauss-Jordan
        para (K - ω²M)v = 0
        """
        print(f"\n" + "=" * 60)
        print(f"CÁLCULO DE AUTOVECTOR PARA ω² = {omega2:.6f}")
        print("=" * 60)
        
        # Construir matriz (K - ω²M)
        K_minus_omega2M = Matrix([
            [self.K.data[0][0] - omega2 * self.M.data[0][0], 
             self.K.data[0][1] - omega2 * self.M.data[0][1]],
            [self.K.data[1][0] - omega2 * self.M.data[1][0],
             self.K.data[1][1] - omega2 * self.M.data[1][1]]
        ])
        
        print(f"Matriz (K - ω²M):")
        print(f"  [{K_minus_omega2M.data[0][0]:.6f}, {K_minus_omega2M.data[0][1]:.6f}]")
        print(f"  [{K_minus_omega2M.data[1][0]:.6f}, {K_minus_omega2M.data[1][1]:.6f}]")
        
        # Aplicar RREF
        K_rref, pivot_cols = rref(K_minus_omega2M)
        
        print(f"\nMatriz en RREF:")
        print(f"  [{K_rref.data[0][0]:.6f}, {K_rref.data[0][1]:.6f}]")
        print(f"  [{K_rref.data[1][0]:.6f}, {K_rref.data[1][1]:.6f}]")
        
        # Encontrar base del espacio nulo
        null_basis = nullspace_basis(K_minus_omega2M)
        
        if not null_basis:
            raise ValueError("No se encontró autovector para este ω²")
        
        eigenvector = null_basis[0]  # Tomar el primer vector de la base
        
        print(f"\nAutovector crudo: [{eigenvector[0]:.6f}, {eigenvector[1]:.6f}]")
        
        return eigenvector
    
    def normalize_eigenvector(self, eigenvector, omega2):
        """
        Normaliza el autovector con producto M-ponderado: vᵢMvᵢ = 1
        """
        # Calcular norma M-ponderada
        v0, v1 = eigenvector[0], eigenvector[1]
        m_norm = v0 * self.M.data[0][0] * v0 + v1 * self.M.data[1][1] * v1
        
        if m_norm <= 0:
            raise ValueError("Norma M-ponderada no positiva")
        
        normalization_factor = 1.0 / np.sqrt(m_norm)
        normalized_vector = [v0 * normalization_factor, v1 * normalization_factor]
        
        # Verificar normalización
        v0_norm, v1_norm = normalized_vector
        check_norm = v0_norm * self.M.data[0][0] * v0_norm + v1_norm * self.M.data[1][1] * v1_norm
        
        print(f"Factor de normalización: {normalization_factor:.6f}")
        print(f"Autovector normalizado: [{v0_norm:.6f}, {v1_norm:.6f}]")
        print(f"Verificación vᵢMvᵢ = {check_norm:.6f}")
        
        return normalized_vector
    
    def modal_diagonalization(self, omega2_list, eigenvectors):
        """
        c) Diagonalización modal: construye Φ y verifica ΦᵀMΦ ≈ I, ΦᵀKΦ ≈ Ω²
        """
        print("\n" + "=" * 60)
        print("DIAGONALIZACIÓN MODAL")
        print("=" * 60)
        
        # Construir matriz de modos Φ = [v1 v2]
        Phi_data = [
            [eigenvectors[0][0], eigenvectors[1][0]],
            [eigenvectors[0][1], eigenvectors[1][1]]
        ]
        Phi = Matrix(Phi_data)
        
        print("Matriz de modos Φ:")
        print(f"  [{Phi.data[0][0]:.6f}, {Phi.data[0][1]:.6f}]")
        print(f"  [{Phi.data[1][0]:.6f}, {Phi.data[1][1]:.6f}]")
        
        # Calcular ΦᵀMΦ
        Phi_T = Phi.transpose()
        M_Phi = self.M * Phi
        PhiT_M_Phi = Phi_T * M_Phi
        
        print(f"\nΦᵀMΦ (debería ser ≈ I):")
        print(f"  [{PhiT_M_Phi.data[0][0]:.6f}, {PhiT_M_Phi.data[0][1]:.6f}]")
        print(f"  [{PhiT_M_Phi.data[1][0]:.6f}, {PhiT_M_Phi.data[1][1]:.6f}]")
        
        # Calcular ΦᵀKΦ
        K_Phi = self.K * Phi
        PhiT_K_Phi = Phi_T * K_Phi
        
        Omega2_expected = Matrix([[omega2_list[0], 0], [0, omega2_list[1]]])
        
        print(f"\nΦᵀKΦ (debería ser ≈ Ω²):")
        print(f"  [{PhiT_K_Phi.data[0][0]:.6f}, {PhiT_K_Phi.data[0][1]:.6f}]")
        print(f"  [{PhiT_K_Phi.data[1][0]:.6f}, {PhiT_K_Phi.data[1][1]:.6f}]")
        print(f"Ω² esperado: diag([{omega2_list[0]:.6f}, {omega2_list[1]:.6f}])")
        
        return Phi, PhiT_M_Phi, PhiT_K_Phi
    
    def calculate_residuals(self, PhiT_M_Phi, PhiT_K_Phi, omega2_list):
        """
        d) Calcula residuos de las condiciones de diagonalización
        """
        print("\n" + "=" * 60)
        print("CÁLCULO DE RESIDUOS")
        print("=" * 60)
        
        # Residuo para ΦᵀMΦ ≈ I
        identity_diff = Matrix([
            [PhiT_M_Phi.data[0][0] - 1, PhiT_M_Phi.data[0][1]],
            [PhiT_M_Phi.data[1][0], PhiT_M_Phi.data[1][1] - 1]
        ])
        
        residual_M = max(
            abs(identity_diff.data[0][0]), abs(identity_diff.data[0][1]),
            abs(identity_diff.data[1][0]), abs(identity_diff.data[1][1])
        )
        
        # Residuo para ΦᵀKΦ ≈ Ω²
        omega2_diff = Matrix([
            [PhiT_K_Phi.data[0][0] - omega2_list[0], PhiT_K_Phi.data[0][1]],
            [PhiT_K_Phi.data[1][0], PhiT_K_Phi.data[1][1] - omega2_list[1]]
        ])
        
        residual_K = max(
            abs(omega2_diff.data[0][0]), abs(omega2_diff.data[0][1]),
            abs(omega2_diff.data[1][0]), abs(omega2_diff.data[1][1])
        )
        
        print(f"Residuo ∥ΦᵀMΦ - I∥∞ = {residual_M:.2e}")
        print(f"Residuo ∥ΦᵀKΦ - Ω²∥∞ = {residual_K:.2e}")
        
        return residual_M, residual_K

def main():
    """Función principal"""
    
    # Crear sistema
    system = TwoMassVibration(m1=1, m2=2, k1=3, k2=4, kc=1)
    
    # a) Calcular autovalores
    omega2_list, omega_list = system.solve_eigenvalues_determinant()
    
    eigenvectors = []
    normalized_eigenvectors = []
    
    # b) Calcular autovectores para cada autovalor
    for i, omega2 in enumerate(omega2_list):
        print(f"\n>>> Procesando modo {i+1}:")
        
        # Calcular autovector
        eigenvector = system.solve_eigenvectors_gauss_jordan(omega2)
        
        # Normalizar autovector
        normalized_vector = system.normalize_eigenvector(eigenvector, omega2)
        
        eigenvectors.append(eigenvector)
        normalized_eigenvectors.append(normalized_vector)
    
    # c) Diagonalización modal
    Phi, PhiT_M_Phi, PhiT_K_Phi = system.modal_diagonalization(omega2_list, normalized_eigenvectors)
    
    # d) Cálculo de residuos
    residual_M, residual_K = system.calculate_residuals(PhiT_M_Phi, PhiT_K_Phi, omega2_list)
    
    # Reporte final
    print("\n" + "=" * 60)
    print("REPORTE FINAL")
    print("=" * 60)
    
    print("\nFRECUENCIAS NATURALES:")
    for i, (omega2, omega) in enumerate(zip(omega2_list, omega_list)):
        print(f"Modo {i+1}: ω² = {omega2:.6f}, ω = {omega:.6f}")
    
    print("\nMODOS NORMALIZADOS (M-ortonormales):")
    for i, vec in enumerate(normalized_eigenvectors):
        print(f"v{i+1} = [{vec[0]:.6f}, {vec[1]:.6f}]")
    
    print("\nPRECISIÓN DE LA DIAGONALIZACIÓN:")
    print(f"∥ΦᵀMΦ - I∥∞ = {residual_M:.2e}")
    print(f"∥ΦᵀKΦ - Ω²∥∞ = {residual_K:.2e}")
    
    # Interpretación física
    print("\nINTERPRETACIÓN FÍSICA:")
    v1 = normalized_eigenvectors[0]
    v2 = normalized_eigenvectors[1]
    
    print(f"Modo 1 (ω₁ = {omega_list[0]:.3f}):")
    if abs(v1[0]) > abs(v1[1]):
        ratio1 = v1[1] / v1[0] if v1[0] != 0 else float('inf')
        print(f"  Las masas se mueven en { 'la misma' if ratio1 > 0 else 'opuesta' } dirección")
        print(f"  Amplitud relativa m2/m1 = {ratio1:.3f}")
    
    print(f"Modo 2 (ω₂ = {omega_list[1]:.3f}):")
    if abs(v2[0]) > abs(v2[1]):
        ratio2 = v2[1] / v2[0] if v2[0] != 0 else float('inf')
        print(f"  Las masas se mueven en { 'la misma' if ratio2 > 0 else 'opuesta' } dirección")
        print(f"  Amplitud relativa m2/m1 = {ratio2:.3f}")

if __name__ == "__main__":
    main()