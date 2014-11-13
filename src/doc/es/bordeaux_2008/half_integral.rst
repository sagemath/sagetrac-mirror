Formas Con Peso De Media Integral
================================

Algorítmo De Basmaji
--------------------

Basmaji (página 55 de su tésis de Essen, "Ein Algorithmus zur Berechnung von
Hecke-Operatoren und Anwendungen auf modulare Kurven" (Un algorítmo para calcular operadores de Hecke y aplicaciones a curvas modulares),
http://wstein.org/scans/papers/basmaji/).

Séa :math:`S = S_{k+1}(\varepsilon)` el espacio de formas cuspidales de peso
entero pár :math:`k+1` y caracter :math:`\varepsilon = \chi
\psi^{(k+1)/2}`, donde :math:`\psi` és el caracter no trivial módulo 4 de Dirichlet.
Séa :math:`U` el subespacio de :math:`S \times S` de
elementos :math:`(a,b)` tales que :math:`\Theta_2 a = \Theta_3 b`. Entonces
:math:`U` es isomorfa a :math:`S_{k/2}(\chi)` a traves de el mapa
:math:`(a,b) \mapsto a/\Theta_3`.

Éste algorítmo está implementado en Sage. Estoy seguro de que podría
implementarse de manera que séa mucho más rápido que la presente
implementación...

::

    sage: half_integral_weight_modform_basis(DirichletGroup(16,QQ).1, 3, 10)
    []
    sage: half_integral_weight_modform_basis(DirichletGroup(16,QQ).1, 5, 10)
    [q - 2*q^3 - 2*q^5 + 4*q^7 - q^9 + O(q^10)]
    sage: half_integral_weight_modform_basis(DirichletGroup(16*7).0^2,3,30)
    [q - 2*q^2 - q^9 + 2*q^14 + 6*q^18 - 2*q^21 - 4*q^22 - q^25 + O(q^30), 
     q^2 - q^14 - 3*q^18 + 2*q^22 + O(q^30), 
     q^4 - q^8 - q^16 + q^28 + O(q^30), q^7 - 2*q^15 + O(q^30)]
