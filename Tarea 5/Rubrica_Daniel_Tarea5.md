# Rúbrica — Ángel Daniel Cruz — Tarea 5

Evaluación concisa por ejercicio.

## 1) Newton multidimensional (R²) — 30/30
- Implementaste Jacobiano analítico 2×2, resolución cerrada y criterio de paro por norma.
- Sistema círculo–recta bien definido, punto inicial cercano a 1/√2; converges y reportas residuo y errores vs exacta.

## 2) Euler implícito 2D con Newton — 28/30
- Reutilizas Newton, armas G y J_G=I−hJf; usas h=0.1 en [0,10].
- Comparas con solución exacta en varios tiempos y reportas errores; buen uso de Jacobiano constante.
- Falta gráfica de comparación/errores (solo tabla) y el Newton interno es básico (sin control de fallos/iteraciones reportados), por eso −2 pts.

## 3) Kepler (E y L_z) — 35/40
- Implementas modelo 4D y Jacobiano analítico; integras con h=0.02, ~5 periodos.
- Calculas E y L_z en cada paso, graficas trayectoria y conservación; reportas variaciones relativas máximas.
- El Newton interno es muy simplificado (paso proporcional al residual) y el horizonte es menor a lo solicitado (~10 periodos), por eso −5 pts.

## Calificación final
- Total: **93/100**.
- Sugerencias: añade gráfica de errores para el oscilador; en Kepler usa Newton completo (resolver J_G Δz = −G) y extiende a ~10 periodos con un h pequeño para evaluar mejor la conservación.
