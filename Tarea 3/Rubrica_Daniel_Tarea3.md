# Rúbrica — Tarea 3 (Ángel Daniel Cruz Flores) — Actualizada 2025-11-14



Total: 50/100

Desglose por problema:

- P1 (Bisección, viga en voladizo): 19/20
- P2 (Enfriamiento de Newton, Newton–Raphson): 10/20
- P3 (Secante: f(x) = e^{x^2} ln(x^2) − x): 16/20
- P4 (Híbrido Newton–Bisección): 0/20
- P5 (Deuterón): 5/20

---

## P1 — Bisección: viga en voladizo f(ω) = cos(ω) cosh(ω) + 1

- Evidencia encontrada:
  - Implementación de bisección `mBis(ecVib, 0, 4)` y graficado de f en [0, 4]; semilog de |f(ω_n)| vs iteración.
- Referencias/canon:
  - Menor raíz positiva: ω* ≈ 1.875104068711961 (tolerancia 1e−5 o mejor).
  - Semilog(|f|) debe decaer de forma lineal en escala log (bisección reduce intervalo a la mitad por iteración).
- Observaciones y retroalimentación:
  - Bien el encuadre [0, 4] y la gráfica; faltó incluir línea de tolerancia en el semilog y un criterio explícito (por ejemplo, |Δω| < TOL·max(1,|ω|) y/o |f| < TOL).
- Calificación: 19/20

## P2 — Enfriamiento de Newton (T(t) = 20 + (90−20)e^{−0.07 t}) y Newton–Raphson

- Evidencia encontrada:
  - Definición correcta de f(t) para T(t)=50°C; uso de Newton–Raphson con t0=10; gráficas de T(t) y semilog(|f|).
-- Observaciones (opcionales de mejora):
  - Derivada implementada con el signo del exponente incorrecto:
    - Código actual: `DTermic(t) = (-0.07)*(90 - 20) * e^{+0.07 t}`.
    - Debe ser: `DTermic(t) = (-0.07)*(90 - 20) * e^{−0.07 t}`.
  - Esto puede afectar severamente la convergencia de Newton y los pasos calculados.
- Referencias/canon:
  - Solución: t* ≈ 12.10425515 min (convergencia en ~4 iteraciones con |error| ~ 1e−8).
- Sugerencia de corrección (código):
  - `def DTermic(t): return (-0.07)*(90-20) * myFunctions.E(-0.07*t)`
- Calificación: 10/20

## P3 — Secante: f(x) = e^{x^2} ln(x^2) − x

- Evidencia encontrada:
  - Implementación de la secante con x0=1.5, x1=0.5; gráfica de f en [−2,2] y semilog(|f|).
  - La forma `E2(x)*(2*ln(|x|)) − x` es equivalente a `e^{x^2} ln(x^2) − x` (para x ≠ 0), correcto.
-- Observaciones (opcionales):
  - No se muestra búsqueda de la raíz negativa ni análisis de sensibilidad a condiciones iniciales (ambos solicitados en la guía estricta).
  - No se incluye línea de tolerancia en el semilog ni criterio explícito de paro.
- Referencias/canon:
  - Raíz positiva: x+ ≈ 1.1624063469; raíz negativa: x− ≈ −0.8105015386.
- Calificación: 16/20

## P4 — Método híbrido Newton–Bisección

- Evidencia encontrada:
  - Se implementa una función `mBN`, pero no se prueba en funciones de test (f1, f2) ni se incluyen resultados o gráficas.
-- Observaciones (necesita corrección para subir puntaje):
  - El “punto medio” se calcula como `m = a + 0.5*(a-b)`; esto no es el punto medio. Debe ser `m = 0.5*(a+b)`.
  - El encuadre [a,b] solo se actualiza cuando el paso de Newton sale del intervalo; si Newton cae dentro de [a,b], no se reduce el intervalo, rompiendo la garantía de convergencia por bisección.
  - Criterio de paro incompleto (no combina residuo |f| con |Δx| ni tolerancias relativas de forma clara).
- Referencias/canon (lo esperado):
  - Probar con dos funciones (p.ej., f1 y f2) y mostrar que el híbrido converge en ~4–6 pasos con buena robustez; semilog(|f|) con línea de tolerancia.
- Calificación: 0/20

## P5 — Deuterón (unidades MeV–fm)

- Evidencia encontrada:
  - Se implementan `h(E)` y `Dh(E)`, se prueban bisección, Newton y secante; se grafica h(E) en E ∈ [−2,0]. Reportas que Newton no converge.
-- Observaciones (opcionales / correcciones sugeridas):
  - Inconsistencias de unidades/constantes y en el parámetro a: usas a=1.45 en h y a=1.54 en Dh.
  - Definiciones de k y β con factores de conversión no estándar; esto impide reproducir resultados canónicos y arruina Newton.
  - La derivada Dh(E) no corresponde a la derivada de g(E)=k cot(ka)+β en MeV–fm; manejo de signos y singularidades poco claro.
- Referencias/canon (objetivos):
  - Unidades: ħc ≈ 197.3269804 MeV·fm; masa m ≈ 938 MeV (masa reducida en MeV/c² equivalente, en el formalismo estándar).
  - Ecuación: g(E) = k cot(ka) + β = 0, con k = sqrt(2m(V0−|E|))/ħ y β = sqrt(2m|E|)/ħ, en MeV–fm.
  - Energía ligada: E* ≈ −1.3805743 MeV (Newton converge en ~6 iter.), bisección ~42 iter.; residuo hasta ~1e−15.
  - Profundidad mínima del pozo para estado ligado: V0_min ≈ 48.72 MeV.
- Sugerencias de corrección:
  - Fijar a=1.45 fm en todo el código; usar m=938 MeV y ħc=197.32698 MeV·fm; reescalar g(E) exclusivamente en MeV–fm.
  - Implementar g(E), g′(E) consistentes; comparar Newton/Bisección/Secante en semilog(|g|) con línea de tolerancia; reportar V0_min.
- Calificación: 5/20

---

## Sugerencias opcionales de mejora

1) P2 (derivada correcta):

```python
# Reemplaza por:
def DTermic(t):
    return (-0.07)*(90 - 20) * myFunctions.E(-0.07 * t)
```

2) P4 (núcleo del híbrido):
- Punto medio correcto: `m = 0.5*(a + b)`.
- Lógica típica: intentar paso de Newton desde m; si x_new ∉ [a,b] o no reduce |f| de forma aceptable, hacer bisección y actualizar el encuadre; combinar criterios de paro: |Δx| < TOL·max(1,|x|) y |f(x)| < TOL.
- Probar con dos funciones de referencia (f1, f2) y reportar iteraciones, residuo y gráfica semilog(|f|) con TOL.

3) P5 (física y unidades):
- Usa ħc = 197.32698 MeV·fm y m = 938 MeV.
- Define k = sqrt(2 m (V0 − |E|)) / ħ y β = sqrt(2 m |E|) / ħ en el mismo sistema de unidades (MeV–fm), y g(E) = k cot(ka) + β.
- Objetivos a verificar: E* ≈ −1.3805743 MeV y V0_min ≈ 48.72 MeV; Newton debe converger en ~6 pasos si g′(E) es correcto.

---

## Criterios generales (aplican a todos los problemas)
- Usa criterios de paro combinados: |Δx| < TOL·max(1,|x|) y |f(x)| < TOL.
- En las gráficas semilog(|f|) añade una línea horizontal de TOL para evidenciar convergencia.
- Reporta tablas/resúmenes con iteraciones, aproximaciones, |f| y error relativo.

## Nota final
Entrega con buen intento general en P1–P3, pero P2 tiene un bug crítico en la derivada, P4 no está operativo ni probado y P5 requiere rehacer las unidades y la derivada para reproducir los valores canónicos.