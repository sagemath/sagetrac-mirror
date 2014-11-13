Introducción
============

¿Que es Sage?
-------------

Sage (see http://sagemath.org) és un sistema de software matemático
extenso para cómputos en muchas áreas de la matemática pura y aplicada.
Programamos a Sage utilizando Python (see http://python.org),
el lenguaje de programación de uso corriente, o su variante compilada Cython.
És tambien muy fácil usar eficientemente el código escrito en C/C++ desde Sage.

El autor de éste artículo inició el proyecto Sage en el 2005.

Sage es libre y de código abierto, dando a entender que puedes cambiar cualquier parte
de Sage y redistribuir el resultado sin tener que pagar por alguna licencia,
y Sage tambien puede aventajar al poder del software matemático comercial
como el de Magma y Mathematica, si sucede que
tienes acceso a esos sistemas comerciales de código cerrado.

Éste documento no asume algún conocimiento previo ya sea de Python o Sage.
Nuestra meta es ayudar a los que estudian teoría de números a hacer cómputos que involucren
cámpos de números y formas modulares usando Sage.

TODO: Overview of Article

A medida que léas este artículo, por favor prueba cada ejemplo en Sage,
y asegurate que las cosas funcionen como yo pretendo, y ház todos los ejercicios.
Además, debes experimentar introduciendo ejemplos similares y
verificando que el resultado que obtengas está de acuerdo con lo que esperas.

Usando Sage
-----------

Para usar Sage, instálalo en tu computador, y usa ya séa la línea de comandos
o inicia el notebook de Sage tecleando ``notebook()`` en la línea de comandos.

Mostramos algunas sesiones de Sage a continuación::

    sage: factor(123456)
    2^6 * 3 * 643


Ésto significa que si tecléas ``factor(123456)`` como entrada de datos en Sage, entonces
obtendrás ``2^6 * 3 * 643`` como resultado. Si estás usando la línea de comandos de Sage,
tecléa ``factor(123456)`` y presiona enter; Si estás usando el notebook de Sage a traves
de tu web browser, tecléa ``factor(123456)`` en una celda de entrada de datos y presiona shift-enter;
en la celda de resultados verás ``2^6 * 3 * 643``.

Después de probar el comando ``factor`` en el párrafo anterior
(¡hazlo ahora!), debes tratar de factorizar otros números.

.. note::

    ¿Qué pasa si factorizas un número negativo? ¿Un número racional?


También puedes dibujar gráficos tanto en 2D como en 3D usando Sage.
Por ejemplo, las siguiente entrada de datos traza el número de divisores primos
de cada entero positivo hasta :math:`500`.

::

    sage: line([(n, len(factor(n))) for n in [1..500]])


Y este ejemplo dibuja un trazo similar en 3D::

    sage: v = [[len(factor(n*m)) for n in [1..15]] for m in [1..15]]
    sage: list_plot3d(v, interpolation_type='nn')


El Ecosistema Sage-Pari-Magma
-----------------------------

* *La diferencia principal entre Sage y Pari es que Sage es mucho más
  grande que Pari con un rango mucho más amplio de funcionalidad, y tiene
  muchos más tipos de datos y objetos mucho más estructurados.* Sage de hecho
  incluye a Pari, y una típica instalación de Sage se lleva cási un gigabyte de
  espacio en disco, mientras que una una típica instalación de Pari és mucho más liviana, usando
  sólo unos cuantos megabytes. Hay muchos algorítmos de teoría de números
  que están incluídos en Sage, los cuales nunca se han implementado en Pari, y
  Sage posee gráficos en 2D y 3D que pueden ser útiles para visualizar
  ideas en teoría de números además de un interfáz gráfico del usuario. Tanto Pari como
  Sage son libres y de código abierto, lo cual significa que cualquier persona puede leer o cambiar
  cualquier cosa en cualquiera de los dos programas y el software és grátis.

* *La más grande diferencia entre Sage y Magma es de que Magma és
  de código cerrado, no es libre y és dificil de extender por los usuarios.* Esto
  significa que la mayor parte de Magma no puede cambiarse excepto por los desarrolladores
  del núcleo de Magma, ya que que la misma Magma tiene arriba de dos millones de líneas
  de código C compilado, combinados con cási medio millón de líneas de
  código interpretado de Magma (que nadie puede leer ni modificar).
  Al diseñar a Sage, tomamos algunas de las excelentes ideas de diseño
  de Magma, tales como el padre, el elemento y la jerarquía de categorías.

* *Cualquier matemático que seriamente quiere realizar trabajos de cómputo extensos
  en teoría algebraica de números y geometría aritmética se le urge enérgicamente
  a que se familiarize con los tres sistemas*, ya que todos ellos tienen sus pros
  y sus contras. Pari es pulcro y pequeño, Magma tiene funcionalidad única
  para cómputos en geometría aritmética, y Sage tiene un amplio rango
  de funcionalidad en la mayoría de las areas de la matemática, una grán
  comunidad de desarrolladores y un código único y nuevo.
