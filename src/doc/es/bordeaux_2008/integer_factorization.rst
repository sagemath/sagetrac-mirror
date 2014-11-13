Factorización De Enteros
========================


Criba Cuadrática
----------------

La criba cuadrática de Bill Hart está incluída en Sage. La criba cuadrática
és el mejor algorítmo para factorizar números de la forma :math:`pq` de 
hasta alrededor de 100 dígitos. Involucra la búsqueda de relaciones, resolviendo
un problema de álgebra lineal modulo :math:`2`, luego factorizando :math:`n`
utilizando una relación :math:`x^2 \equiv y^2 \mod n`.

::

    sage: qsieve(next_prime(2^90)*next_prime(2^91), time=True)   # aún sin probar
    ([1237940039285380274899124357, 2475880078570760549798248507],
     '14.94user 0.53system 0:15.72elapsed 98%CPU (0avgtext+0avgdata 0maxresident)k')

El uso de ``qsieve`` és dos veces más rápido que el comando general ``factor`` de Sage en
este ejemplo. Nótese que el comando general ``factor`` de Sage no hace nada más que
invocar a la función ``factor`` de la biblioteca en C de Pari.

::

    sage: time factor(next_prime(2^90)*next_prime(2^91))     # aún sin probar
    CPU times: user 28.71 s, sys: 0.28 s, total: 28.98 s
    Wall time: 29.38 s
    1237940039285380274899124357 * 2475880078570760549798248507

Obviamente, el comando ``factor`` de Sage no debiera solamente llamar a Pari, pero
todavía nadie se ha acercado a reescribirla.

GMP-ECM
-------

El método GMP-ECM de Paul Zimmerman está incluído en Sage. El algorítmo
de factorización para curvas elípticas (ECM) es el mejor algorítmo para factorizar
números de la forma :math:`n=pm`, donde :math:`p` no es "demasiado
grande". ECM es un algorítmo que se debe a Hendrik Lenstra, que funciona
"fingiendo" que :math:`n` és primo, escogiendo una curva elíptica al azar
sobre :math:`\ZZ/n\ZZ` y efectuando aritmética sobre esa curva
--si algo sale mal cuando se efectúa la aritmética, factorizamos a
:math:`n`.

En el siguiente ejemplo, GMP-ECM está por encima de 10 veces más rápido que
la función genérica ``factor`` de Sage. De nuevo, ésto enfatiza que el
comando genérico ``factor`` de Sage se beneficiaría de una reescritura
que utilice a GMP-ECM y a qsieve.

::

    sage: time ecm.factor(next_prime(2^40) * next_prime(2^300))    # not tested
    CPU times: user 0.85 s, sys: 0.01 s, total: 0.86 s
    Wall time: 1.73 s
    [1099511627791,
     2037035976334486086268445688409378161051468393665936250636140449354381299763336706183397533]
    sage: time factor(next_prime(2^40) * next_prime(2^300))        # not tested
    CPU times: user 23.82 s, sys: 0.04 s, total: 23.86 s
    Wall time: 24.35 s
    1099511627791 * 2037035976334486086268445688409378161051468393665936250636140449354381299763336706183397533
