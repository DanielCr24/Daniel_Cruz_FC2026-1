# Rúbrica — Daniel Cruz — Tarea 4



## 1) Factorización LU con pivoteo parcial — 19/20
- Ensamble A tridiagonal y b para calor 1D (q=50 sin(πx), T0=100, T1=0) correcto; usa `Matrix` y `lu_solve_pp` propios.
- Compara contra `numpy.linalg.solve`, reporta error ∥·∥∞, tiempos y curva T(x) para N={50,100,500}.
- Mejora: el conteo de operaciones es solo O(N^3) genérico; podrías registrar flops/ops reales o tiempos repetidos para escalar mejor.

## 2) Gauss–Jordan en vibraciones pequeñas — 19/20
- Autovalores por determinante cerrada (2x2) correctos y ordenados; imprime ω y ω².
- Autovectores vía RREF y `nullspace_basis`; normalización M-ponderada, Φ, verificas ΦᵀMΦ≈I y ΦᵀKΦ≈Ω² con residuos.
- Mejora: añade una tabla compacta de residuos/valores y comenta sensibilidad a parámetros de masas/resortes.

## 3) Clase TriangularMatrix — 20/20
- Almacena n(n+1)/2, getters/setters con validación, `toDense`, `transpose` upper↔lower, `matvec`, suma entre triangulares y con densas, escalares y comparación con tolerancia.
- Pruebas incluidas (upper/lower, suma, matvec vs denso, transpuestas, escalabilidad n∈{256,512,1024}).

## 4) QR (Gram–Schmidt modificado) — 19/20
- `qr_algorithm_no_shift` con MGS, acumula Q_total, monitorea norma fuera de diagonal, compara autovalores con NumPy y verifica ortogonalidad de QᵀQ.
- Analiza convergencia y grafica errores; matriz Laplaciana 3×3 como caso base.
- Mejora: usar desplazamientos aceleraría convergencia; incluye tiempo y más iteraciones para ver la tasa.

## 5) Mínimos cuadrados con Householder — 20/20
- `householder_qr` aplicado a matriz de diseño cuadrática, resuelve Rβ=Qᵀy por back-substitution; verifica ortogonalidad de Q.
- Datos parabólicos con ruido, coeficientes β comparados con teóricos, errores relativos y gráficas de ajuste y residuos con estadísticas.

## Calificación final
- Total: **97/100**.
- Sugerencias generales: tabular tiempos/errores para LU y QR; añade comentario breve sobre estabilidad numérica en QR sin desplazamiento y sobre cómo el ruido afecta β en mínimos cuadrados.
