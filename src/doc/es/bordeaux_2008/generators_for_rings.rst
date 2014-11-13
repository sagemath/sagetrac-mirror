Generadores Para Anillos De Formas Modulares
============================================

Calculando Generadores
----------------------

Para cualquier subgrupo de congruencia :math:`\Gamma`, la suma directa

.. math:: 

   M(\Gamma) =  \bigoplus_{k\geq 0} M_k(\Gamma)

és un anillo, ya que el producto de formas modulares
:math:`f\in M_k(\Gamma)` y :math:`g \in M_{k'}(\Gamma)` es
un elemento :math:`fg \in M_{k+k'}(\Gamma)`. Sage puede calcular
posibles generadores para anillos de formas modulares, pero corrientemente no
demuestra ninguno de estos resultados.

Verificamos el enunciado que se demostró en "A Course in Arithmetic" (Un Curso En Aritmética) de Serre
de que :math:`E_4` y :math:`E_6` generan el espacio de formas modulares de nivel uno.

::

    sage: ModularFormsRing(SL2Z).generators(prec=4)
    [(4, 1 + 240*q + 2160*q^2 + 6720*q^3 + O(q^4)), 
     (6, 1 - 504*q - 16632*q^2 - 122976*q^3 + O(q^4))]

¿Alguna vez te has preguntado qué formas generan al anillo :math:`M(\Gamma_0(2))`?
Resulta ser que una forma de peso 2 y una forma de peso 4 son suficientes.

::

    sage: ModularFormsRing(Gamma0(2)).generators(prec=12)
    [(2, 1 + 24*q + 24*q^2 + 96*q^3 + 24*q^4 + 144*q^5 + 96*q^6 + 192*q^7 + 24*q^8 + 312*q^9 + 144*q^10 + 288*q^11 + O(q^12)), 
    (4, 1 + 240*q^2 + 2160*q^4 + 6720*q^6 + 17520*q^8 + 30240*q^10 + O(q^12))]

Aqui estan los generadores para :math:`M(\Gamma_0(3))`. Nótese que los
elementos de peso :math:`6` son ahora requeridos, ademas de los
pesos :math:`2` y :math:`4`.

::

    sage: ModularFormsRing(Gamma0(3)).generators()
    [(2, 1 + 12*q + 36*q^2 + 12*q^3 + 84*q^4 + 72*q^5 + 36*q^6 + 96*q^7 + 180*q^8 + 12*q^9 + O(q^10)), 
    (4, 1 + 240*q^3 + 2160*q^6 + 6720*q^9 + O(q^10)), 
    (6, q - 6*q^2 + 9*q^3 + 4*q^4 + 6*q^5 - 54*q^6 - 40*q^7 + 168*q^8 + 81*q^9 + O(q^10))]

Tambien podemos manejar anillos de formas modulares para subgrupos de congruencia impar, pero
con la acostumbrada advertencia de que no podemos calcular formas de peso 1. Asi que éstos son los
elementos que generan el anillo graduado de formas de peso 0 o :math:`\ge 2`.

::

    sage: ModularFormsRing(Gamma1(3)).generators()
    [(2, 1 + 12*q + 36*q^2 + 12*q^3 + 84*q^4 + 72*q^5 + 36*q^6 + 96*q^7 + 180*q^8 + 12*q^9 + O(q^10)), 
    (3, 1 + 54*q^2 + 72*q^3 + 432*q^5 + 270*q^6 + 918*q^8 + 720*q^9 + O(q^10)), 
    (3, q + 3*q^2 + 9*q^3 + 13*q^4 + 24*q^5 + 27*q^6 + 50*q^7 + 51*q^8 + 81*q^9 + O(q^10)), 
    (4, 1 + 240*q^3 + 2160*q^6 + 6720*q^9 + O(q^10))]


