.. escape-backslashes

Démo Sage: Introduction
=======================

`Sage <http://sagemath.org>`_ est un logiciel libre (licence GPL) qui aggrège
d'autres logiciels libres et librairies:
    - `pari <http://pari.math.u-bordeaux.fr/>`_  (théorie des nombres),
    - `gap  <http://www.gap-system.org/>`_ (théorie des groupes),
    - `maxima <http://maxima.sourceforge.net/>`_ (calcul symbolique),
    - et plein d'autres: `numpy <www.numpy.org/>`_, `scipy <www.scipy.org/>`_,
      `matplotlib <http://matplotlib.org/>`_, `m4ri <http://m4ri.sagemath.org/>`_,
      `mpfr <http://www.mpfr.org/>`_, `gmp <https://gmplib.org/>`_, ...

Interface en Python/Ipython
---------------------------

Sage est basé sur le langage `Python <http://www.python.org>`_ qui est très
populaire (programmation web, interfaces graphiques, script, ...) et facile
d'accès.

Python est un langage expressif
+++++++++++++++++++++++++++++++


$\Big\{17n\ \Big|\ n \in \{0,1,\ldots, 9\}\text{ and }n\text{ is odd}\Big\}$

::

    sage: S = {17*n for n in range(10) if n%2 == 1}
    sage: S

    sage: 124 in S

    sage: sum(S)

    sage: {3*i for i in S}


Sage lui ajoute des objets mathématiques
++++++++++++++++++++++++++++++++++++++++

::

    sage: 8324074213.factor()

::

    sage: m = matrix(ZZ,3,3,[0,3,-2,1,4,3,0,0,1])

    sage: m.eigenvalues()

    sage: m.inverse()

::

    sage: R.<x> = PolynomialRing(ZZ,'x')

    sage: R

    sage: P = 6*x^4 + 6*x^3 - 6*x^2 - 12*x - 12
    sage: P.factor()

    sage: P2 = P.change_ring(QQ)
    
    sage: P2.factor()

    sage: P3 = P.change_ring(QQbar)

    sage: P3.factor()


Orienté objet, autocompletion, doc, sources
+++++++++++++++++++++++++++++++++++++++++++

Python est un langage orienté objet. Le notebook de Sage a certaines facilités d'utilisation:

    - autocomplétion avec la touche <TAB>
    - accès à la documentation avec "?"
    - accès au code source avec "??" 

Calculer l'intégrale d'une fonction symbolique.

::

    sage: f(x) = sin(x)^2 -sin(x)

    sage: f

    sage: f.in


Exercice: Dessiner le graphe de Petersen. Quel est l'algorithme utilisé pour
trouver un vertex cover de ce graphe ?

::

    sage: G = grap

    sage: # edit here

    sage: G.vertex

    sage: # edit here



Ensembles, itérateurs
+++++++++++++++++++++

Les itérateurs sont des "générateurs à la demande". Avantages par rapport à un
tableau ou une liste: peu de consomation en mémoire.

::

    sage: Primes()

    sage: 13 in Primes()

    sage: Jumeaux = ((p,p+2) for p in Primes() if (p+2).is_prime())

    sage: [Jumeaux.next() for i in range(10)]

    sage: [Jumeaux.next() for i in range(20)]

Calculette
----------

Symbolique vs algébrique

::

    sage: # edit here

Calcul intégral

::

    sage: integral(e^(-x^2),x,-Infinity,Infinity)


Racines::

    sage: f(x) = x^5 - 1/3*x^2 - 5*x + 1

    sage: plot(f, xmin=-2, xmax=2)

    sage: r1 = find_root(f,-2,-1)
    sage: r1

    sage: r2 = find_root(f,0,1)
    sage: r2

    sage: r3 = find_root(f,1,2)
    sage: r3

    sage: plot(f, xmin=-2, xmax=2) + point2d([(r1,0),(r2,0),(r3,0)], pointsize=50, color='red')


Latex::

    sage: M = Matrix(QQ,[[1,2,3],[4,5,6],[7,8,9]]); M

    sage: latex(M)

    sage: M.parent()

    sage: latex(M.parent())


Graphiques::

    sage: x, y = var('x,y')
    sage: plot3d(sin(x-y)*y*cos(x),(x,-3,3),(y,-3,3))


Interaction::

    sage: var('x')
    sage: @interact
    sage: def g(f=sin(x), c=0, n=(1..30),
    ...         xinterval=range_slider(-10, 10, 1, default=(-8,8), label="x-interval"),
    ...         yinterval=range_slider(-50, 50, 1, default=(-3,3), label="y-interval")):
    ...     x0 = c
    ...     degree = n
    ...     xmin,xmax = xinterval
    ...     ymin,ymax = yinterval
    ...     p   = plot(f, xmin, xmax, thickness=4)
    ...     dot = point((x0,f(x=x0)),pointsize=80,rgbcolor=(1,0,0))
    ...     ft = f.taylor(x,x0,degree)
    ...     pt = plot(ft, xmin, xmax, color='red', thickness=2, fill=f)
    ...     show(dot + p + pt, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
    ...     html('$f(x)\;=\;%s$'%latex(f))
    ...     html('$P_{%s}(x)\;=\;%s+R_{%s}(x)$'%(degree,latex(ft),degree))


Programmation linéaire (réelle ou entière)
------------------------------------------

(voir: http://fr.wikipedia.org/wiki/Optimisation_lin%C3%A9aire)

Demandons à Sage de résoudre le problème linéaire suivant:

Maximiser:

    `2x_0 + x_1 + 3x_2`

Sous les contraintes:

    `x_0 + 2x_1 \leq 4`

    `5x_2 - x_1 \leq 8`

    `x_0, x_1, x_2 \geq 0`

::

    sage: p = MixedIntegerLinearProgram()
    sage: x = p.new_variable(nonnegative=True)

::

    sage: p.set_objective(2*x[0] + x[1] + 3*x[2])

::

    sage: p.add_constraint(x[0] + 2*x[1] <= 4)
    sage: p.add_constraint(5*x[2] - 3*x[1] <= 8)

::

    sage: p.solve()

::

    sage: p.get_values(x)

::

    sage: P = p.polyhedron()
    sage: P.plot()

::

    sage: P.plot(fill=False) + points(P.integral_points())


Quelques liens
--------------

Les indispensables
++++++++++++++++++

- Un forum pour poser des questions : http://ask.sagemath.org
- Calcul mathématique avec Sage, un livre sur Sage en français, mis a jour sur http://sagebook.gforge.inria.fr/ 

TP introductifs
+++++++++++++++

Voici une sélection de TP élémentaires pour vous apprendre à utiliser Sage et
Python. Tous ces documents font partie de Sage. Vous pouvez les retrouver dans
la documentation sur votre ordinateur (depuis le notebook, cliquez sur "Help"
en haut à droite puis sur "Thematic Tutorials") ou bien sur
http://sagemath.org/doc

Si vous voulez explorer la partie calculette de Sage et faire des dessins:

- Calcul: faire des fonctions, des intégrales, des graphiques élémentaires : http://sagemath.org/doc/prep/Calculus.html
- Graphiques : http://sagemath.org/doc/prep/Advanced-2DPlotting.html 

Pour faire de la programmation:

- Introduction à Sage avec un peu de programmation (fonctions, boucles, ...) : http://sagemath.org/doc/prep/Programming.html
- Apprendre à faire des itérateurs : http://sagemath.org/doc/thematic_tutorials/tutorial-comprehensions.html
- Introduction assez complète à la programmation en Python et en Sage : http://sagemath.org/doc/thematic_tutorials/tutorial-programming-python.html
- Des problèmes de mathématiques à résoudre avec des programmes : http://projecteuler.net

Programmation linéaire
++++++++++++++++++++++

- La programmation linéaire dans Sage : http://www.steinertriples.fr/ncohen/tut/LP/|Utiliser
- Exemples de problèmes : http://www.steinertriples.fr/ncohen/tut/LP_examples/
- Worksheet de la doc. officielle de Sage(en anglais) : http://sagemath.org/doc/thematic_tutorials/linear_programming.html


Python Scientifique
+++++++++++++++++++

Tout paquet Python s'installe facilement avec la commande (depuis un
terminal)::

    sage -pip install <paquet>

  Une liste assez exhaustive de paquets Python pour le calcul scientifique : http://scipy.org/topical-software.html


----

Un peu de combinatoire
----------------------

Partitions d'un entier
++++++++++++++++++++++

::

    sage: P = Partitions(12)
    sage: P

    sage: [5,4,1,1,1] in P

    sage: [5,4,2,1,1] in P

    sage: P.list()

    sage: Partitions(100000).cardinality()

    sage: Permutations(20).random_element()


Permutations
++++++++++++

::

    sage: s = Permutation([5,3,2,6,4,8,9,7,1])
    sage: s

    sage: (p,q) = s.robinson_schensted()

    sage: p.pp()

    sage: q.pp()

    sage: p.row_stabilizer()


Points entiers d'un polytope
++++++++++++++++++++++++++++

::

    sage: V = ZZ^3
    sage: vectors = [V.random_element(x=0,y=10) for _ in range(6)]
    sage: L = LatticePolytope(vectors)
    sage: L.plot3d()

    sage: L.npoints()


Graphes à isomorphisme près
+++++++++++++++++++++++++++

::

    sage: show(graphs(5, lambda G: G.size() <= 4))


Une géodésique sur une surface dans `R^3`
+++++++++++++++++++++++++++++++++++++++++

::

    sage: E = surfaces.Ellipsoid(axes=(1,3,2))

    sage: E_plot = E.plot()

    sage: E_plot.show(aspect_ratio=1)

    sage: xy = (1.,1.)        # un point

    sage: E.point(xy)         # ses coordonnees en 3d
    (1.29192658172643, 3.45464871341284, 2.84147098480790)

    sage: v = (0.245, 0.312)   # un vecteur

    sage: E.tangent_vector(xy,v)
    (-0.253239333370952, -0.448190681935137, 0.337148638861719)    

    sage: pts = [E.point(x[1]) for x in E.geodesics_numerical(xy, v, (0,10000,1000))]

    sage: pts[0]
    (0.291926581726429, 1.36394614023852, 1.68294196961579)
    sage: pts[1]
    (0.281703243506732, 1.34596947848564, 1.69629081709544)

    sage: G = E.plot(opacity=0.5)
    sage: G += line3d(pts, color='red')
    sage: G += arrow3d(start=E.point(xy), end=E.point(xy)+E.tangent_vector(xy,v), color='red')

    sage: G.show(aspect_ratio=1)
