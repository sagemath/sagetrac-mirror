.. -*- coding: utf-8 -*-
.. _agregation.algebre_lineaire:

==================================================================================
Option Algèbre et Calcul Formel de l'Agrégation de Mathématiques: Algèbre linéaire
==================================================================================

.. MODULEAUTHOR:: `Nicolas M. Thiéry <http://Nicolas.Thiery.name/>`_ <Nicolas.Thiery at u-psud.fr>

*******************************
Formes normales et applications
*******************************

L'algorithme de Gauss revisité
==============================

On se place dans un corps `K` quelconque.

.. TOPIC:: Exercice: matrices à deux lignes

    Soit `M` une matrice générique à deux lignes. Écrire sous forme de
    multiplication à gauche par une matrice `2\times 2` le pivot de
    Gauss appliqué à `M`.


.. TOPIC:: Remarque

    Si `M` est obtenue de `N` par l'algorithme du pivot de Gauß, alors
    `M=PM` où `P` est une matrice inversible (éventuellement de
    déterminant `1`).

Disons ici que deux matrices `M` et `N` de `M_{n,m}(K)` sont
*équivalentes* (modulo l'action de `GL_n(K)` à gauche) s'il existe une
matrice inversible `P` telle que `M=PN`.

.. TOPIC:: Exercice:

    Vérifier que cela définit une relation d'équivalence!

.. TOPIC:: Question

    La remarque précédente dit que si deux matrices `M` et `N` donnent
    la même forme échelon par Gauß, alors elles sont équivalentes.

    Réciproque?

.. TOPIC:: Théorème

   On considère les matrices `n\times m` à coefficients dans un corps
   `K`. La forme échelon réduite donne une *forme normale* pour les
   matrices modulo l'action de `GL_n(K)` à gauche.

Voir par exemple [Storjohan.2004]_ pour une présentation d'ensemble
des formes normales pour les différentes relations d'équivalences sur
les matrices (conjugaison, similitude, ...) sur les différents types
d'anneaux, et les algorithmes associés.


Interprétation géométrique
--------------------------

Décrire un objet comme étant le résultat d'un algorithme est
*opératoire*, mais pas très *conceptuel*. Peut-on faire mieux?

.. TOPIC:: Exercice

    Soient `M` et `N` deux matrices de `M_{n,m}(K)`, que l'on voit
    comme deux paquets de `n` vecteurs de `K^m`. Montrer que `M` et
    `N` sont équivalentes (modulo l'action de `GL_n(K)` à gauche) si
    et seulement si les vecteurs engendrent le même sous-espace
    vectoriel de `K^m`.

    Autrement dit, l'ensemble quotient `GL_n(K) \backslash M_{n,m}(K)`
    représente l'ensemble des sous-espaces vectoriels de dimension au
    plus `n` dans `K^m`. Cet ensemble est naturellement muni d'une
    structure de variété appelée variété Grassmanienne.

.. TOPIC:: Rappel: groupes de permutations

    Pour manipuler un sous-groupe `G` du groupe symétrique `S_n`, on
    avait considéré le sous-groupe `G_{n-1}` des éléments fixant `n`,
    puis ceux fixant `n` et `n-1`, et ainsi de suite récursivement.

    Formellement, on avait considéré la suite des groupes symétriques
    emboîtés:

    .. MATH::

        \{id\} = S_0\subsetneq S_1 \subsetneq \cdots \subsetneq S_n

    et la suite induite des groupes emboîtés `G_i:=G \cap S_i`:

    .. MATH::

        \{id\} = G_0\subset G_1 \subset \cdots \subset G_n=G

    L'étude de `G` se ramenait alors à l'étude des quotients
    successifs `G_i/G_{i-1}`.

Appliquons le même programme.

.. TOPIC:: Définition: Drapeau

    Un drapeau complet d'un espace vectoriel `V` de dimension `n` est
    une suite maximale de sous-espaces strictement emboîtés:

    .. MATH::

        \{0\} = V_0 \subsetneq V_1 \subsetneq \cdots \subsetneq V_n=V

.. TOPIC:: Définition: Drapeau canonique

    À chaque base ordonnée, on peut associer naturellement un drapeau
    complet.  Ici on considérera principalement le drapeau canonique
    associé à la base canonique `e_1,\cdots,e_m` de `V=K^m`:

    .. MATH::

        V_i:=\langle e_{m-i+1} \cdots e_m \rangle

    Note: on prend les éléments dans cet ordre pour que cela colle
    avec nos petites habitudes de calcul du pivot de Gauß. Et pour
    alléger les notations, on utilisera plutôt:

    .. MATH::

        \overline V_i:=\langle e_i \cdots e_m \rangle=V_{n-i+1}

.. TOPIC:: Formes échelon et bases adaptées

    Dans ce formalisme, qu'est-ce qu'une matrice sous forme échelon?

    C'est une base d'un espace vectoriel `E` adaptée à un drapeau
    complet donné. C'est-à-dire une base sur laquelle on peut lire
    immédiatement les sous espaces `E_i:=E\cap \overline V_i`.

    Le pivot de Gauß est un algorithme de calcul de base adaptée.

.. TOPIC:: Définition intrinsèque des colonnes caractéristiques

    Remarque: en passant de `E_{i+1}` à `E_i`, la dimension croît de
    `0` ou de `1`.

    Cela permet de donner une définition intrinsèque de la notion de
    colonnes caractéristiques d'un sous espace vectoriel `E`: les `i`
    tels que la dimension de `E_i` croît strictement. Cela décrit la
    position de `E` par rapport à un drapeau complet fixé.

    Évidemment, sur une forme échelon pour `E`, cela correspond aux
    colonnes `i` pour lesquelles on a un vecteur de la forme
    `e_i+\cdots`.


.. TOPIC:: Formes échelon réduites

    Considérons deux bases adaptées d'un même espace vectoriel
    `E`. Pour `i` une colonne caractéristique, on note `a_i` et `b_i`
    les vecteurs de la forme `a_i=e_i+\cdots` et `b_i=e_i+\cdots`.

    Alors `a_i-b_i\in V_{i+1}`; autrement dit `a_i=b_i` dans le
    quotient `E_i/E_{i+1}`.

    Prendre une forme échelon réduite, c'est faire un choix d'un
    représentant (relativement canonique) `a_i` dans chaque quotient
    `E_i/E_{i+1}`: celui qui a des zéros aux autres colonnes
    caractéristiques.

    Ce formalisme montre que le vecteur `a_i` est intrinsèque à `E`
    (et au choix du drapeau complet). En particulier il est clair
    qu'il est complètement indépendant des autres coefficients de la
    forme échelon réduite, même si opératoirement le calcul de `a_i`
    par Gauß passe par ceux-ci.

.. TOPIC:: Remarque

    La permutation `P` apparaissant dans le calcul de l'algorithme de
    Gauß a une interprétation géométrique naturelle.

    Les variétés Grassmaniennes et ses variantes (variétés de
    drapeaux, ...) et leur multiples généralisations sont l'objet
    d'études approfondies en géométrie. La combinatoire y joue un rôle
    important: l'apparition d'une permutation `P` dans le pivot de
    Gauß est le prototype du type de lien.

Applications des formes échelon
-------------------------------

.. TOPIC:: Exercice: résolution d'équations linéaires

    Soit `E` un ensemble d'équations linéaires/affines. Retrouver les
    algorithmes usuels de résolution: existence de solution,
    dimension, base et paramétrisation de l'espace des solutions.

.. TOPIC:: Exercice: calcul avec les sous espaces vectoriels

    On considère des sous espaces `E`, `F`, ... de `V=K^n` donnés par
    des générateurs ou des équations. Donner des algorithmes pour:

    #.  Déterminer une base de `E`.

    #.  Tester si un vecteur appartient à `E`.

    #.  Tester si `E=F`.

    #.  Tester si deux vecteurs `x` et `y` de `V` sont égaux modulo `E`

    #.  Calculer l'orthogonal d'un sous-espace vectoriel

    #.  Calculer la somme `E+F` et l'intersection `E\cap F` de deux espaces vectoriels

    #.  Calculer la sous-algèbre de `V` engendrée par `E`
	(en supposant `V` muni d'une structure d'algèbre `(V,+,.,*)`)

        Plus généralement: clôture de `E` sous des opérations linéaires

    #.  Calculer dans l'espace quotient `E/F`

    #.  Cas de la dimension infinie?


.. TOPIC:: Exercice: calcul avec les morphismes

    Soit `\phi` une application linéaire entre deux espaces vectoriels
    `E` et `F` de dimension fini. Donner des algorithmes pour:

    #.  Calculer le noyau de `\phi`

    #.  Calculer l'image de `\phi`

    #.  Calculer l'image réciproque par `\phi` d'un vecteur `f` de `F`

    #.  Arithmétique: composition, combinaison linéaires, inverse

    #.  Calculer le polynôme caractéristique

    #.  Calculer les valeurs propres de `\phi`

    #.  Calculer les espaces propres de `\phi`


Forme de Hermite
================

On considère maintenant l'anneau `\ZZ`. On est maintenant en train de
travailler avec des `\ZZ`-modules au lieu d'espaces vectoriels.
Peut-on procéder comme précédemment?

.. TOPIC:: Exercice: matrices à deux lignes

    Exemple::

        sage: M = matrix([[10,1,2], [6,2,-1]]); M
	[10  1  2]
	[ 6  2 -1]

	sage:

    #.  Quel candidat pour une forme échelon?

    #.  Interprétation en terme de multiplication par une matrice?

    #.  Interprétation en terme de sous-espace engendré?

    #.  Cette forme échelon est elle réduite?

	::

            sage: M.echelon_form()
	    [  2   3  -4]
	    [  0   7 -11]

    #.  Description du quotient?

.. TOPIC:: Remarque clef

    Soit `\begin{pmatrix}a\\b\end{pmatrix}` un vecteur de `\ZZ^2`, et
    `r, u,v` les résultats du pgcd étendu de `a` et `b`: `r = a\wedge
    b = ua+bv`. Posons:

    .. MATH::
        M :=
        \begin{pmatrix}
	    u        & v        \\
	    \frac br & \frac ar \\
	\end{pmatrix}

    alors: `M\begin{pmatrix}a\\b\end{pmatrix} = \begin{pmatrix}r\\0\end{pmatrix}` et `M\in GL(\ZZ)`!


Moralité: la majeure partie de ce que l'on a vu précédemment
s'applique mutatis-mutandis. L'algèbre linéaire sur `\ZZ` n'est pas
foncièrement plus compliquée ou coûteuse que sur un corps.

Il y a juste quelques points techniques à traiter, qui apparaissent
déjà en dimension `1`:

.. TOPIC:: Exercice: Résolution

    Déterminer l'ensemble des solutions entières de l'équation
    `6x+4y+10z=18`.

.. TOPIC:: Exercice: Torsion

    #.  Donner un exemple de quotient d'un module libre `\ZZ^n` qui
        n'est pas isomorphe à un module libre.

    #.  Donner un exemple de drapeau infini

    #.  Existe-t'il des drapeaux croissants infinis?

	.. TODO:: Donner un exemple illustrant le phénomène en
		  ajoutant des vecteurs au hasard les uns après les
		  autres

.. TOPIC:: Généralisations

    Tout ce que l'on vient de dire se généralise immédiatement pour un
    anneau principal quelconque comme `A=\QQ[x]`; à condition bien
    entendu que `A` soit *constructif*, et en particulier, qu'il y ait
    un algorithme pour calculer le PGCD étendu.

Application: classification des groupes abéliens de type fini
-------------------------------------------------------------

.. TOPIC:: Exercice

    Soit `G` un groupe additif abélien engendré par un nombre fini `n`
    d'éléments, mettons `a`, `b`, `c`, avec `n=3`.

    #.  Que peut-on dire sur l'ensemble des relations entre ces
        éléments?

    #.  En déduire la structure de `G`


Gauß sans fractions et Gauß-Bareiss
===================================

Soit `A` un anneau. Par exemple un anneau de polynômes multivariés
`A=\QQ[x,y]`. Qu'est-ce qui subsiste de tout cela?

Exemple: dimension 1
--------------------

Un `A`-sous-module de `A^1` est juste un idéal de `A`.

Exemple: `\langle x^2y, xy^2\rangle`

Calcul avec les idéaux et sous-modules: bases de Gröbner

Combinatoire sous-jacente: idéaux monomiaux!

Algorithme de Gauß sans fraction
--------------------------------

Explorons un exemple::

    sage: pretty_print_default()

    sage: A = QQ['a']
    sage: a = A.gen()
    sage: M = matrix(A, random_matrix(ZZ, 3, 8)); M[0,0] = a; M
    [ a  2  0  0  0  2  2  2]
    [ 1  2  0 -2  0  0  1  2]
    [ 1 -2  1 -2 -1 -2 -2 -2]
    sage: N = copy(M)

    sage: M[1] = M[0,0] * M[1] - M[1,0] * M[0]
    sage: M[2] = M[0,0] * M[2] - M[2,0] * M[0]
    sage: M

    sage: M[2] = M[1,1] * M[2] - M[2,1] * M[1]
    sage: M

Algorithme de Gauß-Bareiss
--------------------------

Revenons sur notre exemple::

    sage: M

    sage: M[2] = M[2] / a
    sage: M

    sage: N[:,:3]
    sage: N[:,:3].det()

    sage: 

.. TOPIC:: Proposition

    Soit `M` une matrice sur un anneau intègre.

    Après avoir traité de les `i` premières colonnes, `M_{i,i}` est le
    déterminant du mineur de la matrice d'origine, et ce déterminant
    divise toutes les coefficients `M_{i',j'}` avec `i',j'>i`.

    Par récurrence, on obtient l'algorithme de Gauß-Bareiss qui permet
    de calculer le déterminant.

.. TOPIC:: Exercice

    Vérifier la proposition dans le cas d'une matrice triangulaire
    supérieure.


Conclusion
==========

.. TODO:: !!!

TP
==

.. TOPIC:: Exercice: algorithme de Gauß-Bareiss

    Dans tout cet exercice, on pourra supposer que la matrice d'entrée
    est inversible, voire que ses `n` premiers mineurs sont non nuls
    (pas de permutation des lignes nécessaire).

    #.  (Échauffement) Écrire une fonction qui met une matrice à
	coefficients dans un corps sous forme échelon à l'aide de
	l'algorithme de Gauß. Vérifier votre programme pour::

	    sage: M = matrix([[2, 1, 3], [1, 4, 9], [1, 8, 27]]); M

    #.  Écrire une fonction qui met une matrice à coefficients entiers
	sous forme échelon à l'aide de l'algorithme de
	Gauß-Bareiss. Vérifier votre programme pour la matrice
	ci-dessus, puis sur une matrice aléatoire de grande taille.

	Évaluer la complexité pratique en prenant des matrices
	aléatoire de taille `n=2^k`. Comparer avec ce que l'on obtient
	avec Gauss, et avec Gauss sur un corps fini.

	Qu'en pensez-vous?

    #.  En déduire une fonction qui calcule le déterminant d'une
        matrice à coefficients entiers.

    #.  Faire la même chose pour des matrices à coefficients
	polynomiaux univariés.

    #.  En déduire une fonction qui calcule le polynôme
        caractéristique d'une matrice.


.. TOPIC:: Algèbre linéaire, représentations des monoïdes et Chaînes de Markov

    Voir: :ref:`agregation.bibliotheque_tsetlin`.

    Ce texte est à approcher comme les textes de l'agrégation: il
    s'agit d'un menu à la carte; vous pouvez choisir d'étudier
    certains points, pas tous, pas nécessairement dans l'ordre, et de
    façon plus ou moins fouillée. Vous pouvez aussi vous poser
    d'autres questions que celles indiquées plus bas. L'objectif final
    est de concevoir un mini-développement de 5 minutes comportant une
    partie traitée sur ordinateur et, si possible, des représentations
    graphiques de vos résultats.

Textes connexes
===============

- `Algorithme Page Rank de Google <http://nicolas.thiery.name/Enseignement/Agregation/Textes/PageRankGoogle.pdf>`_

- `Résolution de systèmes linéaires en entiers <http://nicolas.thiery.name/Enseignement/Agregation/Textes/560-ResolutionDeSystemesLineairesEnEntiers.pdf>`_

- `Pseudo inverses de matrices <http://nicolas.thiery.name/Enseignement/Agregation/Textes/PseudoInverseMatrice.pdf>`_

**********************
Programmation linéaire
**********************

Simplexe, dualité, applications

Voir Cours/TP dans OldNotes

*******************
Quelques références
*******************

.. [Storjohan.2004] `Algorithms for Matrix Canonical Forms <https://cs.uwaterloo.ca/~astorjoh/diss2up.pdf>`_,
   Arne Storjohan, PhD Thesis,
   Department of Computer Science,
   Swiss Federal Institute of Technology -- ETH, 2000
