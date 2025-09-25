# Rúbrica: Daniel_Cruz_FC2026-1 - Tarea 1

Referencia: ejercicios 1,2,3,5 del enunciado (Problema 4 descartado). Total = 4.0 puntos (P1, P2, P3, P5).

Resumen de comprobaciones

- Archivo revisado: `Tarea_1.py`
- Evaluación por problema:
  - Problema 1 (Unidad y conversión): 1.0 / 1.0
  - Problema 2 (Primos): 0.9 / 1.0
  - Problema 3 (nameFile): 1.0 / 1.0
  - Problema 5 (Media y desviación): 0.5 / 1.0

Total: 3.4 / 4.0

Evidencia y comentarios

Problema 1:
- `convertDist(d)` y `convertVel(v)` implementadas y se prueban con los valores del cometa Halley.
- Observación: las funciones validan que la entrada sea `float` (estricto); es aceptable pero limitaría entrada de `int`.

Problema 2:
- `nPrimes(n)` genera correctamente la lista de los primeros n primos.
- `twinPrimes(n)` crea la lista de primos pero el bucle final itera `while j < n` y usa índices `primes[j-1]`/`primes[j]`, lo cual es incorrecto: debe iterar sobre el rango de índices de `primes` (p.ej. 1..len(primes)-1) y comprobar `primes[i+1] - primes[i] == 2`; además el par debe agregarse en orden `(primes[i], primes[i+1])`. Penalización menor por ese bug lógico.

Problema 3:
- `nameFile(simName, *args, **kwargs)` construye el nombre con el prefijo, reemplaza puntos por guiones bajos y añade `.dat`. Cumple el formato requerido y maneja listas/tuplas en `args`.

Problema 5:
- `myAverage`, `myDesvest2s`, `myDesvest1s` están presentes pero contienen errores graves:
  - Las validaciones usan `isinstance(x, int)` en lugar de comprobar que `x` sea una secuencia (lista/tuple) — esto produce TypeError para entradas válidas.
  - En `myDesvest2s` la expresión para la varianza usa `(1/len(x)-1)` en vez de `1/(len(x)-1)`, lo que produce un factor incorrecto.
  - En `myDesvest1s` falta la raíz cuadrada final (devuelve la varianza sin sqrt) y hay la misma validación errónea de tipo.
- Por estas fallas, la implementación de las desviaciones no es correcta y requiere corrección; otorgo 0.5 por haber entendido las fórmulas pero implementar mal.

Recomendaciones:
- Corregir las comprobaciones de tipo: usar `isinstance(x, (list, tuple))` o simplemente validar iterabilidad y longitud.
- En `twinPrimes` iterar sobre `range(1, len(primes))` y comprobar `primes[i] - primes[i-1] == 2` y añadir `(primes[i-1], primes[i])`.
- Corregir las fórmulas de desviación:
  - `desvest_dos = sqrt(sum((xi - mu)**2) / (n - 1))`
  - `desvest_uno = sqrt((sum(xi**2) - n*mu**2) / (n - 1))`
- Considerar aceptar `int` y `float` en las funciones de conversión (usar `isinstance(d, (int, float))`) para mayor flexibilidad.

Archivos modificados:
- Ninguno (solo se creó esta rúbrica).

Estado: Hecho.
