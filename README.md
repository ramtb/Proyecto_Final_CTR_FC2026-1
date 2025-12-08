# Simulación Estocástica de Crecimiento Tumoral (Modelo de Monte Carlo)



 



Este repositorio contiene el proyecto final para la materia de **Física Computacional**. Implementa una simulación de crecimiento tumoral avascular utilizando **Autómatas Celulares** y el **Método de Monte Carlo**.



El proyecto modela la dinámica de proliferación, migración y muerte celular (apoptosis y necrosis) en una malla bidimensional, validando los resultados microscópicos con modelos macroscópicos de crecimiento poblacional (Curva Logística/Gompertz).



-----



## 📋 Descripción del Proyecto



El objetivo principal es simular la evolución temporal de un tumor basándose en reglas probabilísticas locales. El sistema evoluciona a través de **Pasos de Monte Carlo**, donde cada célula activa tiene la oportunidad de realizar una acción basada en su entorno y sus parámetros biológicos intrínsecos.



### Estados de la Célula



La simulación considera tres estados posibles para cada celda de la malla:



1.  ⬜ **Sana / Vacía (0):** Espacio disponible o tejido sano.

2.  🟥 **Cancerosa (1):** Célula activa capaz de reproducirse, migrar o morir.

3.  ⬛ **Necrótica (2):** Célula muerta por hipoxia/densidad, capaz de liberar toxicidad antes de ser limpiada.



-----



## 🚀 Características Principales



  * **Modelo Híbrido Estocástico:** Combina autómatas celulares con probabilidades de transición (Monte Carlo).

  * **Dinámica de Necrosis Avanzada ($\alpha_t$):** Implementación de un factor `alpha` dinámico. La probabilidad de necrosis depende de la **densidad global** del tumor y de un índice de **resistencia celular** intrínseco.

      * Formula: $\alpha_t = (1 - R) \times D_{global}$

  * **Ciclo de Vida Completo:**

      * **Inhibición por Contacto:** Las células rodeadas no pueden reproducirse.

      * **Necrosis vs. Apoptosis:** Diferenciación entre muerte limpia (espacio libre) y muerte tóxica (necrosis) según la disponibilidad de espacio.

      * **Limpieza:** El sistema simula la respuesta inmune/limpieza eliminando células necróticas con una probabilidad $P_{clean}$.

  * **Optimización Computacional:** Uso de estructuras de datos tipo `Set` para el seguimiento de células activas, reduciendo la complejidad de búsqueda de $O(N \times M)$ a $O(N_{activas})$.

  * **Visualización en Tiempo Real:** Generación de mapas de calor de la malla y gráficas de métricas poblacionales.



-----



## 📊 Validación Física y Matemática



El proyecto incluye un módulo de validación teórica ("El Regreso") que compara los datos de la simulación con la solución analítica de la Ecuación Logística de crecimiento:



$$\frac{dN}{dt} = rN \left( 1 - \frac{N}{K} \right)$$



Se calcula automáticamente la tasa de crecimiento efectiva ($r$) basada en las probabilidades microscópicas ($P_{repro}$, $P_{dead}$) ajustada por un factor geométrico de superficie, demostrando que el comportamiento emergente del autómata celular coincide con la teoría macroscópica.



-----



## 🛠️ Tecnologías Utilizadas



El proyecto fue construido utilizando **Python** puro para la lógica de simulación (sin librerías de autómatas externas) para demostrar el dominio de los algoritmos.



  * **NumPy:** Exclusivamente para la generación de números aleatorios (`random`).

  * **Matplotlib:** Para la visualización de la malla y graficación de datos (Crecimiento, Densidad, Derivada $dN/dt$).



-----



## 🖼️ Visualizaciones



El código genera cuatro tipos de salidas visuales:



1.  **Grid de Evolución:** Visualización paso a paso del tumor (Blanco/Rojo/Negro).

2.  **Curva de Crecimiento:** Número de células vs. Tiempo (Curva Sigmoidal).

3.  **Densidad Global:** Evolución de la saturación del tejido.

4.  **Tasa de Crecimiento:** Derivada numérica ($dN/dt$) suavizada con media móvil.



-----



## 🔧 Instalación y Uso



1.  Clonar el repositorio:



    ```bash

    git clone https://github.com/tu-usuario/simulacion-tumor-montecarlo.git

    ```



2.  Instalar dependencias:



    ```bash

    pip install numpy matplotlib

    ```



3.  Ejecutar la simulación en el documento Simulation.ipynb dentro de la carpeta Funcional.

-----

## 🛠️ Código Adicionales

1.  **Aproximación matemática:** El cual permite observar el fenomeno con operaciones mátematicas, siguiendo la función de Gompertz.



2.  **Obtención de la probabilidad:** El cual permite obtener los valores para las probabilidades que se utilizan durante la simulación.

-----



## ✒️ Autores

Tellez Becerra Angel Ramses

Canela Lupercio Paola

Robles Ávila Fernando Daniel



-----



### Notas Adicionales



El código incluye comentarios detallados sobre la lógica de recanalización de probabilidades y el manejo de condiciones de frontera. Se recomienda revisar el archivo principal para ajustar parámetros como `RESISTANCE_SCORE` o `ROWS/COLS`.

