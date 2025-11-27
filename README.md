# Simulador de Epidemias SIRV en Redes

Este proyecto es una aplicaci0n web interactiva creada con R y Shiny que permite simular y visualizar la propagaci0n de una enfermedad utilizando modelos epidemiol0gicos. La aplicaci0n ofrece dos enfopatques complementarios: un **modelo determinista (SIRV)** basado en ecuaciones diferenciales y un **modelo estoc0stico basado en eventos** que se ejecuta sobre diferentes tipos de redes complejas.

## Caracter0sticas

- **Doble Enfoque de Simulaci0n**:
  - **Modelo Determinista**: Simulaci0n de campo medio (mean-field) que muestra la evoluci0n promedio de la epidemia.
  - **Modelo Estoc0stico**: Simulaci0n basada en el algoritmo de Gillespie sobre una red, capturando la aleatoriedad inherente a las interacciones individuales.
- **Parametrizaci0n Interactiva**: Permite al usuario ajustar en tiempo real los par0metros epidemiol0gicos (tasas de infecci0n, recuperaci0n, vacunaci0n) y de la red (tama0o, topolog0a).
- **Topolog0as de Red**: Soporta la generaci0n de redes Regulares, Erd0os-R0nyi y Barab0asi-Albert.
- **Visualizaciones Detalladas**:
  - **Pesta0a "Simulaciones"**: Muestra las trayectorias individuales de m0ltiples simulaciones estoc0sticas.
  - **Pesta0a "Modelo Determinista"**: Presenta la curva suave y predecible del modelo de ecuaciones diferenciales.
  - **Pesta0a "Red"**: Visualiza la estructura de la red generada y el estado de cada nodo.
  - **Pesta0a "Estad0sticas"**: Agrega los resultados de las simulaciones estoc0sticas para mostrar la media y la desviaci0n est0ndar.

---

## C0mo Funciona

### Modelo Determinista (Ecuaciones Diferenciales)

Este modelo asume que la poblaci0n est0 bien mezclada y describe los cambios en la proporci0n de individuos en cada compartimento (Susceptible, Infectado, Recuperado, Vacunado) a trav0s de un sistema de Ecuaciones Diferenciales Ordinarias (EDO).

El sistema de ecuaciones para el modelo SIRV implementado es:

-   $$ \frac{dS}{dt} = - \frac{\beta S I}{N} - \nu S $$
-   $$ \frac{dI}{dt} = \frac{\beta S I}{N} - \mu I $$
-   $$ \frac{dR}{dt} = \mu I $$
-   $$ \frac{dV}{dt} = \nu S $$

Donde:
-   `S`, `I`, `R`, `V` son los n0meros de individuos en cada compartimento.
-   `N` es la poblaci0n total (`S+I+R+V`).
-   `\beta` (beta) es la tasa de transmisi0n.
-   `\mu` (miu) es la tasa de recuperaci0n.
-   `\nu` (nu) es la tasa de vacunaci0n.

En la aplicaci0n, este sistema se resuelve num0ricamente utilizando la funci0n `ode` del paquete de R `deSolve`, generando curvas suaves que representan la tendencia general de la epidemia.

### Modelo Estoc0stico (Algoritmo de Gillespie)

A diferencia del modelo determinista, el modelo estoc0stico simula la epidemia como una secuencia de **eventos discretos y aleatorios** que ocurren en nodos individuales de una red. Para determinar qu0 evento ocurre y cu0ndo, se utiliza el **Algoritmo de Gillespie** (tambi0n conocido como M0todo Directo o SSA).

Este algoritmo es ideal para sistemas donde el azar juega un papel crucial, como en poblaciones peque0as o en las fases iniciales de un brote.

#### Conceptos Clave

1.  **Estado del Sistema**: El estado de cada nodo en la red (`S`, `I`, `R` o `V`).
2.  **Propensi0n (`\lambda_i`)**: Es la probabilidad por unidad de tiempo de que un evento le ocurra a un nodo `i`. Se calcula para cada nodo:
    -   **Para un nodo Susceptible (S)**: Puede infectarse o vacunarse. Su propensi0n es la suma de las tasas de ambos eventos. La tasa de infecci0n depende de cu0ntos de sus vecinos est0n infectados.
        $$ \lambda_i = (\beta \times k_{infectados}) + \nu $$
        donde `k_infectados` es el n0mero de vecinos infectados del nodo `i`.
    -   **Para un nodo Infectado (I)**: Solo puede recuperarse.
        $$ \lambda_i = \mu $$
    -   **Para un nodo Vacunado (V)**: Puede perder la inmunidad y volver a ser susceptible (si `\omega > 0`).
        $$ \lambda_i = \omega $$
3.  **Propensi0n Total (`\Lambda`)**: Es la suma de las propensiones de todos los nodos de la red. Representa la tasa a la que ocurre *cualquier* evento en el sistema.
    $$ \Lambda = \sum_{i=1}^{N} \lambda_i $$

#### Pasos del Algoritmo de Gillespie

La simulaci0n avanza de evento en evento siguiendo estos pasos:

1.  **Inicializaci0n (t = 0)**:
    -   Se define el estado inicial de cada nodo en la red.
    -   Se calculan las propensiones iniciales `\lambda_i` para cada nodo y la propensi0n total `\Lambda`.

2.  **Paso 1: Calcular el tiempo hasta el pr0ximo evento (`\tau`)**:
    -   Se genera un n0mero aleatorio `u₁` de una distribuci0n uniforme U(0, 1).
    -   El tiempo `\tau` hasta que ocurra el *siguiente* evento se calcula con la f0rmula:
        $$ \tau = -\frac{\ln(u₁)}{\Lambda} $$9
    -   El tiempo de la simulaci0n se actualiza: `t = t + \tau`. Este m0todo genera saltos de tiempo variables: cortos cuando los eventos son muy probables y largos cuando son improbables.

3.  **Paso 2: Decidir qu0 evento ocurre**:
    -   Se genera un segundo n0mero aleatorio `u₂` de U(0, 1).
    -   Se elige el nodo `j` que sufrir0 el cambio de estado. La probabilidad de que se elija un nodo `j` es proporcional a su propensi0n `\lambda_j`. Esto se logra encontrando el primer nodo `j` que satisface la condici0n:
        $$ \sum_{k=1}^{j} \lambda_k \ge u_2 \times \Lambda $$
    -   En la aplicaci0n, esto se implementa eficientemente con `which(cumsum(lambda) >= u2 * Lambda)[1]`.

4.  **Paso 3: Actualizar el estado del sistema**:
    -   Se actualiza el estado del nodo `j` seleccionado. Por ejemplo, si era `S`, podr0a cambiar a `I` o a `V` (determinado por otra tirada aleatoria que compara las tasas `\beta*k_inf` vs. `\nu`).
    -   **Paso crucial**: Se recalculan las propensiones del nodo `j` y de todos sus vecinos, ya que su cambio de estado afecta sus probabilidades de cambiar. Por ejemplo, si `j` pasa de `I` a `R`, sus vecinos susceptibles ahora tienen un menor riesgo de infecci0n.
    -   Se actualiza la propensi0n total `\Lambda`.

5.  **Repetici0n**:
    -   Se vuelve al Paso 1 y se repite el ciclo hasta que se alcanza el tiempo m0ximo de simulaci0n `T` o hasta que `\Lambda = 0` (lo que significa que no pueden ocurrir m0s eventos).

#### Implementaci0n Optimizada

En `app.r`, las funciones `calculateNodePropensities`, `draw_next_event_direct`, y `update_states` implementan este algoritmo de forma vectorizada para mejorar el rendimiento, evitando bucles lentos en R y aprovechando operaciones matriciales y de vectores.

---

## C0mo Ejecutar la Aplicaci0n

1.  **Prerrequisitos**:
    -   Tener instalado R y preferiblemente RStudio.
2.  **Librer0as Requeridas**:
    -   Instala las siguientes librer0as desde la consola de R:
    ```R
    install.packages(c("shiny", "bslib", "igraph", "ggplot2", "dplyr", "tidyr", "deSolve"))
    ```
3.  **Ejecutar**:
    -   Abre el archivo `app.r` en RStudio.
    -   Haz clic en el bot0n **"Run App"** que aparece en la parte superior del editor de c0digo.
