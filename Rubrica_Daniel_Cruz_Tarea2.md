## Rúbrica (simplificada): Daniel_Cruz_FC2026-1 — Tarea 2

- Archivo revisado: `Tarea 2.py`
- Evaluación por problema (1.0 c/u):
  - P1 (Grid y memoria): 0.6 / 1.0
  - P2 (Horner y proyectil): 0.6 / 1.0
  - P3 (Series de Taylor): 0.5 / 1.0
  - P4 (Norma euclidiana robusta): 0.2 / 1.0
  - P5 (Aritmética de intervalos): 0.6 / 1.0
  - P6 (Búsqueda y análisis): 0.3 / 1.0

Total: 2.8 / 6.0

Solo incisos con observaciones

- P1
  - Derivación (b) incorrecta (pasa a 13 iteraciones; no cuadra con el modelo de memoria nueva + anterior).
  - Conclusión “mismos pasos” en (c) es dudosa con 12.4 GB efectivos.
  - `simGrid` simula, pero no se muestran corridas para (a)(b)(c).

- P2
  - Polinomio evaluado es cúbico; el enunciado pide el de grado 4. Falta el término constante 0 y el coeficiente de t^4.
  - La gráfica está comentada y los valores de y(t) no se reportan.

- P3
  - `log_taylor(x)`: para x>2 cae en recursión infinita (usa `-log(1/x)` y luego vuelve a `-log(1/(1/x))`).
  - `sin_taylor`, `cos_taylor`, `tan_taylor` calculan, pero sin reducción de rango (convergen lento en ángulos grandes).

- P4
  - Solo se implementa norma ingenua; falta rutina robusta y la comparación de precisión y rendimiento solicitada.

- P5
  - Se acepta el redondeo a enteros con `to_integral_value` (FLOOR/CEILING) según el criterio de referencia; por ello no se penaliza.
  - Las salidas quedan anchas (p.ej., 10×0.1 → [0, 10], 1/3 → [0, 1]), lo cual es esperable con redondeo entero; faltan pruebas más informativas o validación adicional (p.ej., divisor que toque 0) si aplica.

- P6
  - Búsqueda binaria: `mid = (right - left)//2` ignora el desplazamiento; debe incluir `left`.
  - Medición de tiempo lineal: se guarda el timestamp final en lugar del delta.
  - En “solo búsqueda” se resta el tiempo de ordenamiento al de búsqueda, dando resultados inconsistentes.

Estado: rúbrica generada.

Evidencia de ejecución (resumen)
- P2: “No se pudo encontrar el tiempo de altura máxima” (la bisección no halló raíz en [0,5] con la derivada usada/intervalo elegido).
- P5: Suma 10×0.1 → [0.0, 10.0]; 1/3 → [0.0, 1.0]. Intervalos excesivamente anchos por usar `to_integral_value` (trunca a enteros).
- P6: Tiempos lineal reportan valores tipo 1761025046.s (timestamp, no delta). Desglose binaria muestra ordenamiento + búsqueda; las cifras confirman los problemas anotados.
