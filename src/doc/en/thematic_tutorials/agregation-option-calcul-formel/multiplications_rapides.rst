.. -*- coding: utf-8 -*-
.. _agregation.multiplications_rapides:

==================================================================================
Option Algèbre et Calcul Formel de l'Agrégation de Mathématiques: Produits rapides
==================================================================================

.. MODULEAUTHOR:: `Nicolas M. Thiéry <http://Nicolas.Thiery.name/>`_ <Nicolas.Thiery at u-psud.fr>

.. linkall

*****************************************
Motivation: «Tout» se ramène aux produits
*****************************************

Inversion de matrices, pivot de Gauss
=====================================

.. TOPIC:: Exercice: matrices `2\times 2` génériques

    Soit `M=\begin{pmatrix}a&b\\c&d\end{pmatrix}`

    #. Calculer `M^{-1}`

    #. Effectuer pas à pas le pivot de Gauss sur `M`.

    #. Généralisation à une matrice par blocs `M=\begin{pmatrix}A&B\\C&D\end{pmatrix}`?

    .. TODO:: Peaufiner

Méthode de Newton
=================

Approximation numérique de solution de `f(x) = 0`
-------------------------------------------------

Approximation en série d'une équation implicite `F(G(z)) = 0`
-------------------------------------------------------------

Inversion de séries
-------------------

Division Euclidienne de polynômes
---------------------------------

*****************************
Produits rapides de polynômes
*****************************

Algorithme naïf
===============

Dans la suite, on considère un anneau `K` et deux polynômes dans `K[z]`:

.. MATH::

    A = A(z) = a_0 + a_1 z + \cdots + a_n z^n

.. MATH::

    B = B(z) = b_0 + b_1 z + \cdots + b_n z^m

L'objectif est de calculer les coefficients `c_k` du polynôme `C(z) = A(z)B(z)`

.. TOPIC:: Algorithme naïf

    On se contente d'utiliser la formule `c_k = \sum_{i+j=k} a_i b_j`.

.. TOPIC:: Exercice

    Quelle est la complexité du calcul du produit des polynômes `A(z)`
    et `B(z)` par l'algorithme naïf?


Karatsuba
=========

.. TOPIC:: Exercice

    Donner des formules pour calculer les coefficients du polynôme
    `(a_0+a_1z)(b_0+b_1z)` en fonction de `a_0,a_1,b_0,b_1` utilisant
    un nombre minimal de produits.

Nous allons maintenant appliquer les deux principes suivants:

- «Si vous avez une bonne idée, appliquez la par récurrence, vous
  obtiendrez une meilleure idée.»

- Diviser pour régner!

.. TOPIC:: Étape de récurrence

    Supposons pour que `n=m=2l`, et écrivons

    .. MATH:: A = A_0 + A_1 z^l

    .. MATH:: B = B_0 + B_1 z^l

    où `A_0=A_0(z)` `A_1=A_1(z)`, `B_0=B_0(z)`, `B_1=B_1(z)` sont de degré `\leq l`.

    On peut calculer `AB` en calculant récursivement quatre produits
    de polynômes de degré `l`:

    .. MATH::

        AB = A_0B_0 + ( A_0B_1 + A_1B_0 ) z^l + (A_1B_1)z^{2l}

    Ou seulement avec trois:

    .. MATH::

        AB = A_0B_0 + ( (A_0+A_1)(B_0+B_1) - A_0B_0 - A_1B_1 ) z^l + (A_1B_1)z^{2l}

L'algorithme de Karatsuba consiste à calculer le produits de polynômes
de degré `2^r` en appliquant récursivement l'étape précédente.


.. TOPIC:: Complexité

    L'algorithme de multiplication de Karatsuba est de complexité
    `O(n^{\log_2(3)})\approx O(n^{1.59})`.

.. TOPIC:: Démonstration

    On supposes d'abord que `n=2^r`, et on ne compte que le nombre
    `f(r)` de multiplications requises dans `K`.
    Clairement:

    .. MATH::

        f(r)=3f(r-1)=3^rf(0)=3^r

    Pour calculer le produit de deux polynômes de degré `n`, on les
    complète en polynômes de degré `2^{\lceil \log_2(n)\rceil}`.  Le
    nombre de multiplications dans `K` est alors borné par:

    .. MATH::

        3^{\lceil \log_2 n \rceil} \leq 3. 3^{\log_2 n} = 3.2^{\log_2 3 . \log_2 n} = 3 n^{\log_2 3}

    Il est clair que le nombre d'additions est négligeable (de l'ordre
    de `O(4n\log_2 n)`).

.. TOPIC:: En pratique: implantation

    L'algorithme de Karatsuba, étant plus compliqué en particulier à
    cause de la récursion, est moins performant en petit degré que
    l'algorithme naïf. Aussi les implantations utilisent l'étape de
    récurrence en haut degré, et basculent sur un produit naïf en deçà
    d'un certain seuil.

    Ce seuil est déterminé expérimentalement par bancs d'essais. Dans
    certains cas la détermination du seuil optimal pour une
    architecture donnée est effectuée automatiquement à la
    compilation.

    C'est un principe très général. On l'avait déjà vu avec les tris,
    et on le retrouve par exemple en algèbre linéaire avec la
    bibliothèque ATLAS:

    ATLAS: Automatically Tuned Linear Algebra Software

.. TOPIC:: En pratique: usage

    L'algorithme de Karatsuba requière des soustractions:

    - Il ne s'applique pas aux polynômes sur des semi-anneaux (par
      exemple `\NN[x]`, algèbre tropicale, ...)

    - Il peut poser des problèmes de stabilité numérique en calcul
      approché (flottants, ...)

Produit par évaluation
======================

.. TOPIC:: Remarque stupide

    Si `x_0` est un élément de `K`, et `C(z) = A(z)B(z)` alors:

    .. MATH::

        C(x_0) = A(x_0) B(x_0)

.. TOPIC:: Corollaire

    Soient `x_1,\dots,x_n` des éléments de `K` et munissons `K^n` de
    l'addition et de la multiplication point à point.

    L'application d'évaluation:

    .. MATH::

        \Phi:
	\begin{cases}
	  K[z] &\mapsto (K^n,+,.)\\
	  P(z)  &\mapsto ( P(x_1), \ldots, P(x_n) )
	\end{cases}

    est un morphisme d'algèbre.

    C'est même un isomorphisme si on se restreint à l'ensemble
    `K[z]_n` des polynômes de degré `<n`.

Le produit dans `(K^n,+,.)` est de complexité `n`. Donc il est tentant
d'utiliser cet isomorphisme pour calculer les produits:

.. MATH::

    A(z)B(z) = \Phi^{-1} ( \Phi(A) \Phi(B) )

.. TOPIC:: Problème

    Rentable si le calcul de `\Phi` (évaluation) et de `\Phi^{-1}`
    (par ex. interpolation) est peu coûteux. Pour des points
    quelconques, c'est au moins du `O(n^2)`.

    Comment choisir de bons points d'évaluation?


Produit par transformée de Fourier Discrète
===========================================

Transformée de Fourier Discrète
-------------------------------

.. TOPIC:: Proposition

    Supposons que l'anneau `K` contienne une racine primitive `\omega`
    de l'unité. Alors le morphisme d'algèbre:

    .. MATH::

        DFT_\omega:
	\begin{cases}
	  K[z] &\mapsto (K^n,+,.)\\
	  P(z)  &\mapsto ( P(1), P(w), \ldots, P(w^{n-1}) )
	\end{cases}

    induit un isomorphisme d'algèbre de `K[z] / (z^n-1)`.

.. TOPIC:: Démonstration

    Regarder le noyau + dimension.

.. TOPIC:: Remarque

    On retrouve la même algèbre que dans les codes cycliques; entre
    autres, la multiplication par `x` donne une action du groupe
    cyclique `C_n`.

.. TOPIC:: Exercice

    #. `DFT_\omega` est une application linéaire. Donner sa matrice.

    #. Donner la matrice inverse.


    Indication: `\sum_{k=0}^{n-1} \omega^{ik} = \begin{cases}n&\text{si $j\equiv 0[n]$}\\0&\text{sinon}\end{cases}`

.. TOPIC:: Proposition

    La transformée de Fourier discrète inverse est encore une
    transformée de Fourier discrète:

    .. MATH::

        DFT_\omega^{-1} = \frac 1n DFT_{\omega^{-1}}


.. TOPIC:: Remarque: lien avec la théorie des représentations

    La matrice de `DFT_\omega` est aussi la table des caractères du
    groupe cyclique `C_n`. Le fait qu'elle soit hermitienne à un
    scalaire près est un cas particulier d'une proposition générale
    sur les tables de caractères. L'espace `K[z]/(z^n-1)` se
    décompose en `n` modules simples de dimension `1`, et la
    transformation `DFT_\omega` correspond à la décomposition d'un
    polynôme dans ces modules simples.

    Il existe des notions de transformées de Fourier discrète pour
    d'autres groupes.

Il reste à calculer efficacement la transformée de Fourier discrète.

Transformée de Fourier rapide (FFT: Fast Fourier Transform)
-----------------------------------------------------------

.. TOPIC:: Diviser pour régner

    Supposons que `P` soit un polynôme de degré `n=2k`.

    Noter que `z^{2k} - 1 = (z^k-1) (z^k+1)`.

    Les racines de `z^k-1` sont les racines `k`-ème de l'unité; on
    peut donc utiliser la transformée de Fourier discrète pour évaluer
    `P(z)` dessus. Plus précisément, on calcule

    .. MATH::

        P_+(z) = P(z) [ z^k - 1 ]

    (ce calcul est léger!) et on utilise `DFT_{\omega^2}(P_+(z))`.

    Les racines de `z^k-1` sont les racines `k`-ème de l'unité,
    décalées par un facteur `\omega`. on calcule alors

    .. MATH::

        P_-(z) = P(z) [ z^k + 1 ]

    et on utilise `DFT_{\omega^2}(P_-(\omega z))`.

.. TOPIC:: Algorithme de multiplication par FFT

    On considère une racine `2^k`-ème de l'unité, et on applique
    récursivement l'idée précédente.

    Complexité: O(n\log n), comme pour les tris.


.. TOPIC:: Problème

    Et s'il n'y a pas de racine primitive de l'unité dans `K`?

    On la rajoute!

    Exemple: les corps cyclotomiques obtenus par extension algébrique
    de `\QQ` par un polynôme cyclotomique::

        sage: K = CyclotomicFields(6)
	sage: omega = K.gen()
	sage: omega^6

    Souci: ces corps cyclotomiques nécessitent de calculer dans des
    extensions de corps de haut degré; donc un bon produit; cela
    pourrait boucler!

    Algorithme de Schönhage et Strassen: `O(n\log n\log\log n)`

    Autre souci: on a divisé par `n=2^k`; ce n'est pas forcément
    possible, par exemple en caractéristique `2`!


**************************
Produits rapides d'entiers
**************************

Même principe que pour les polynômes; juste plus technique à cause de
la gestion des retenues. On retrouve le produit par Karatsuba, par
FFT, ...

Ce que l'on a remarqué pour les séries s'applique aux calculs sur les
nombres réels à précision arbitraire.

****************************
Produits rapides de matrices
****************************

.. TOPIC:: Algorithme de Strassen

    Même principe que Karatsuba!

    - Pour multiplier deux matrices `2\times 2`, il existe des formules
      n'utilisant que 7 produits au lieu de 8.

    - On découpe les matrices de taille `2^k` en 4 blocs de taille
      `2^{k-1}` et on utilise les formules ci-dessus récursivement.


    Complexité: `O(n^{\log_2 7}) \approx O(n^{2,8})`

.. TOPIC:: Algorithme de Coppersmith-Winograd, ...

    Complexité: `O(n^{2,3755\cdots})`, `O(n^{2,3736\cdots})`, `O(n^{2,3727\cdots})`

    Inutilisable en pratique ...

    .. TODO: donner une intuition du principe de fonctionnement


**
TP
**

Parcourir les exercices suivants et en piocher un pour préparer une
démonstration courte (5 minutes). Ensuite, jouez avec les exercices
que vous voulez. En fin de séance (mettons vers 11h45), chacun d'entre
vous présentera sa démonstration aux autres.

.. TOPIC:: Exercice: Karatsuba

    #.  Implanter l'algorithme naïf pour multiplier deux polynômes

    #.  Implanter l'algorithme de Karatsuba pour multiplier deux
	polynômes

    #.  Faire un banc d'essai pour ces deux algorithmes, et tracer un
	graphe permettant de comparer simultanément leur complexité
	pratique entre elles et avec leur complexité théorique.

    #.  Avec votre implantation, à partir de quel seuil est-il
	préférable d'utiliser l'algorithme de Karatsuba?

    Prolongements possibles:

    #.  Implanter un algorithme mixte Karatsuba/naïf qui tienne compte
	du seuil obtenu. Comparer.

    #.  Comparer la complexité pratique de votre implantation du
	produit avec celle de la bibliothèque de Sage.

    #.  Deviner, d'après sa complexité pratique, le ou les algorithmes
        utilisés par Sage.

    #.  Implanter le produit de deux entiers par Karatsuba; comparer
        avec l'implantation pour les polynômes.

.. TOPIC:: Transformée de Fourier rapide

    Voir le `sujet de TP <../Jouve/FFT1213.pdf>`_ de l'année dernière.

.. TOPIC:: Exercice: Illustration de Newton numérique

    Réaliser une animation similaire à celle de
    http://en.wikipedia.org/wiki/Newton%27s_method#Newton-Fourier_method

.. TOPIC:: Exercice: Convergence de Newton numérique

    #.  Choisir une équation de la forme `f(x) = 0` et calculer des
        approximations successives `x_0, x_1,\dots,` de l'une de ses
        solutions à l'aide d'une itération de Newton.

    #.  Tracer le graphe du nombre de décimales correctes en fonction
        du nombre d'itérations.

.. TOPIC:: Exercice: Inversion de séries formelle par itération de Newton

    Soit `A(z)` une série formelle dans `K[[z]]` dont on veut calculer
    l'inverse `B(z)=A^{-1}(z)`. En particulier, on supposera que son
    terme constant `a_0=A(0)` est inversible dans `K`.

    On pose la fonction `F(X,z) = A(z) - 1/X`, de sorte que `B(z)`
    satisfait l'équation fonctionnelle implicite `F(B(z), z)=0`.

    #.  Choisir `B_0(z)` tel que `B_0(z)\equiv B(z) [z]`

    #.  Supposer que l'on ait trouvé `B_i(z)` tel que
        `B_i(z)\equiv B(z)[z^k]`.  Appliquer une itération de Newton pour
        retrouver l'expression de `B_{i+1}` vue en cours, et donner sa
        précision (i.e. combien de termes de `B(z)` sont obtenus).

.. TOPIC:: Exercice: Comptage des arbres par itération de Newton

    Cet exercice est un complément pour la section 15.1.2
    «Dénombrement d'arbres par séries génératrices» du livre «Calcul
    Mathématique avec Sage».

    On rappelle que l'ensemble `C` des arbres binaires complets est
    défini récursivement en spécifiant qu'un arbre binaire complet est
    soit une feuille, soit consiste en une racine à laquelle sont
    attachés un sous-arbre gauche et un sous-arbre droit. Soit `C(z)`
    la série génératrices des arbres binaires complets comptés par
    nombres de feuilles.

    #.  Écrire l'équation ensembliste satisfaite par `C`.

    #.  La traduire en équation algébrique satisfaite par `C(z)`.

    #.  Choisir `C_0(z)` tel que `C_0(z)\equiv C(z) [z]`.

    #.  Par itération de Newton, calculer successivement `C_1(z)`,
        `C_2(z)`, ... et indiquer le nombre de termes de `C(z)`
        obtenus à chaque étape.

	Indication: on pourra au choix représenter les `C_i(z)` par:

	- Des fractions rationnelles, en utilisant la commande
	  :func:`taylor` pour les développer en série entière.

	- Des séries tronquées à l'ordre approprié (éventuellement
          représentée par un simple polynôme), en utilisant l'exercice
          précédent pour les calculs d'inverse.


*******************
Quelques références
*******************

.. [MCA 2013] Modern Computer Algebra, Joachim von zur Gathen, Jürgen Gerhard

.. [Salvy 2013] `Newton iteration in Computer Algebra and Combinatorics <http://www.liafa.univ-paris-diderot.fr/jifp/salvy.pdf>`_

.. [Riou] `Notes de cours de Joël Riou <http://www.math.u-psud.fr/~riou/enseignement/2012-2013/mao/cours.pdf>`_

.. [Roblot 2013] `Mini-cours sur l'arithmétique algorithmique <http://math.univ-lyon1.fr/~roblot/ens.html>`_
