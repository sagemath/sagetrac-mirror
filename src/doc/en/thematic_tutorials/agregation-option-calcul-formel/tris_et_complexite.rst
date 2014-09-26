.. -*- coding: utf-8 -*-
.. _agregation.tris_et_complexite:

====================================================================================
Option Algèbre et Calcul Formel de l'Agrégation de Mathématiques: Tris et complexité
====================================================================================

.. MODULEAUTHOR:: `Nicolas M. Thiéry <http://Nicolas.Thiery.name/>`_ <Nicolas.Thiery at u-psud.fr>

****************************
Introduction à la complexité
****************************

Quelques problèmes
==================

- Quelle est le meilleur algorithme pour trouver un nom dans un
  annuaire ?

- Quel est le meilleur méthode pour calculer le déterminant d’une
  matrice ?

- Comment prédire le temps que va mettre un programme pour s’exécuter?

- Comment savoir, entre deux algorithmes, lequel est le plus efficace?

- Comment savoir si un algorithme est optimal?

- Comment déterminer si un problème est insoluble en pratique?

Complexité d'un algorithme
==========================

Exemple: recherche naïve dans une liste
---------------------------------------

Je recherche le nom «Zorro» dans un annuaire en utilisant la méthode
suivante:

#. Je pars du début de l’annuaire;

#. Je compare le nom avec «Zorro»;

#. Si oui, j’ai terminé;

#. Si non, je recommence en 2 avec le nom suivant.

.. TOPIC:: Problème

   Combien est-ce que cela va me prendre de temps?

Synthèse
--------

On s’est donné un problème (rechercher un mot dans un dictionnaire),
un algorithme pour le résoudre (recherche exhaustive). Puis on a
introduit un *modèle de calcul*:

#. Choix de la mesure de *la taille d’une instance problème* (le nombre
   de mots d’un dictionnaire donné)

#. Choix des *opérations élémentaires* (comparer deux mots)

Dans ce modèle, on a cherché le nombre d’opérations élémentaires
effectuées par l’algorithme pour un problème de taille `n`.
C’est ce que l’on appelle la *complexité de l’algorithme*.

En fait, on a vu deux variations:

#. Complexité au pire (`n` opérations)

#. Complexité en moyenne (`\frac{n}{2}` opérations)

À partir de cette information, et en connaissant le temps nécessaire
pour de petites instances du problème on peut évaluer le temps
nécessaire pour résoudre n’importe quelle instance du problème.

.. TODO:: Mentionner la complexité en mémoire

Exercices
---------

.. TOPIC:: Exercice

    Donner des algorithmes et leur complexité au pire et en moyenne
    pour les problèmes suivants:

    #. Calculer la somme de deux matrices carrées

    #. Calculer le produit de deux matrices carrées

    #. Rechercher un élément dans une liste

    #. Calculer le plus grand élément d'une liste

        ::

            sage: %hide
            sage: def plus_grand_element(liste):
            ...       """
            ...       Renvoie le plus grand élément de la liste
            ...       EXAMPLES::
            ...           sage: plus_grand_element([7,3,1,10,4,10,2,9])
            ...           10
            ...           sage: plus_grand_element([7])
            ...           7
            ...       """
            ...       resultat = liste[0]
            ...       for i in range(1, len(liste)-1):
            ...           # Invariant: resultat est le plus grand element de liste[:i]
            ...           assert resultat in liste[:i]
            ...           assert all(resultat >= x for x in liste[:i])
            ...           if liste[i] > resultat:
            ...               resultat = liste[i]
            ...       return resultat
	    ...
	    sage: plus_grand_element([7,3,1,10,4,10,2,9])
	    10

        .. NOTE:: Digression: invariants, preuve et test

    #. Rechercher un élément dans une liste triée

    #. Insérer un élément dans une liste triée

Ordres de grandeur
==================

Exemple: recherche dichotomique
-------------------------------

Quelques courbes de complexité
------------------------------

::

    sage: %hide
    sage: var('n')
    sage: funs = [n^0, log(n, 10), sqrt(n), n, 50*n, n*(log(n,10)), n^2, n^2.3727, n^log(2,7), n^3, 2^n, 10^n, factorial(n)]
    sage: colors = rainbow(len(funs))
    sage: def time_label(s, t): return text(s, (1,t), horizontal_alignment = "left")
    sage: time_labels = sum(time_label(t,s)
    ...                     for t,s in [["seconde", 1], ["minute", 60], ["jour",24*3600], 
    ...                                 ["annee",365*24*3600], ["siecle",100*365*24*3600],["age de l'univers",14*10^9*365*24*3600]])
    sage: sum( plot(f/10^9,
    ...             xmin=1, xmax=(100 if f(n=100)>10^19 else 10^10),
    ...             ymax=10^20,
    ...             scale="loglog",
    ...             color=color, legend_label=repr(f))
    ...        for f,color in zip(funs, colors)) + time_labels

    sage: %hide

Synthèse
--------

La plupart du temps, il suffit d’avoir un ordre de grandeur du nombre
d’opérations: les constantes sont sans grande importance. Un
algorithme en `1000\log_{2}n+50` sera meilleur qu’un algorithme en
`\frac{n}{1000}` dès que l’on s’intéressera à des instances
suffisamment grandes.

Mais voir aussi [CTFM1993]_

.. TOPIC:: Définition

    Soient `f` et `g` deux fonctions de `\mathbb{N}` dans `\mathbb{N}`
    (par exemple les complexités de deux algorithmes).

    On note `f=O(g)` si, asymptotiquement, `f` est au plus du même
    ordre de grandeur que `g`; formellement: il existe une constante
    `a` et un entier `N` tels que `f(n)\leq ag(n)` pour `n\geq N`.

    On note `f=o(g)` si, assymptotiquement, `f` est négligeable devant
    `g`; formellement: pour toute constante `a` il existe `N` tel que
    `f(n)\leq ag(n)` pour `n\geq N`.

.. TOPIC:: Proposition

    Quelques règles de calculs sur les `O()`:

    #. `O(4n+3)=O(n)`

    #. `O(\log n)+O(\log n)=O(\log n)`

    #. `O(n^{2})+O(n)=O(n^{2})`

    #. `O(n^{3})O(n^{2}\log n)=O(n^{5}\log n)`

Exercices
---------

.. TOPIC:: Exercice

    Retrouver les règles de calcul analogues pour les `o()`, et pour
    les opérations mixtes avec `o()` et `O()`.

    .. TODO:: Donner à la place des calculs intéressants

.. TOPIC:: Exercice

    Donner des algorithmes et leur complexité au pire et en moyenne
    pour les problèmes suivants:

    #. Effectuer un pivot de Gauss sur une matrice

        .. NOTE:: Digression: Complexité arithmétique versus complexité binaire

    #. Calculer le déterminant d'une matrice


Complexité d'un problème
========================

.. TOPIC:: Exercice

    #. Donner un algorithme pour rechercher le plus grand élément d’une liste de nombres.
    #. Évaluer la complexité de cet algorithme.
    #. Existe-t’il un meilleur algorithme ?

.. TOPIC:: Définition

    La *complexité d’un problème* est la complexité du meilleur
    algorithme pour le résoudre.

    On dit qu’un algorithme est *optimal* si sa complexité coïncide
    avec celle du problème.

.. TOPIC:: Exercices

    #. Les algorithmes vus précédemment sont-ils optimaux?

    #. Démontrer que la recherche d'un élément dans une liste triée de taille `n` est un problème de complexité `O(\log n)`.

	#. On dispose d’un ordinateur pouvant exécuter `10^{9}` opérations élémentaires par seconde (1GHz). On a un problème (par exemple, chercher un mot dans une liste, calculer le déterminant d’une matrice), et des instances de taille `1,10,100,1000` de ce problème. Enfin, on a plusieurs algorithmes pour résoudre ce problème, dont on connaît les complexités respectives: `O(\log n)`, `O(n)`, `O(n\log n)`, `O(n^{2})`, `O(n^{3})`, `O(n^{10})`, `O(2^{n})`, `O(n!)`, `O(n^{n})`. Évaluer dans chacun des cas le temps nécessaire.

***********************************************************
Comparaison de la complexité de quelques algorithmes de tri
***********************************************************

On a une liste que l’on veut trier, mettons `[7,8,4,2,5,9,3,5]`.

Quelques algorithmes de tri
===========================

Tri sélection
-------------

#. On échange le premier élément avec le plus petit des
   éléments: `2,8,4,7,5,9,3,5`

#. On échange le deuxième élément avec le plus petit des
   éléments restants: `2,3,4,7,5,9,8,5`

#. Etc.

#. Au bout de `k` étapes, les `k` premiers
   éléments sont triés; on échange alors le `k+1`-ième
   élément avec le plus petit des éléments restants.

#. À la fin, la liste est triée: `2,3,4,5,5,7,8,9`.

Tri fusion
----------

#. On groupe les éléments par paquets de deux, et on trie chacun
   de ces paquets: `(7,8),(2,4),(5,9),(3,5)`.

#. On groupe les éléments par paquets de quatres, et on trie
   chacun de ces paquets: `(2,4,7,8),(3,5,5,9)`.

#. ...

#. Au bout de `k` étapes, les paquets de `2^{k}`
   éléments sont triés; on les regroupe par paquets de
   `2^{k+1}` que l’on trie.

#. À la fin, tous les éléments sont dans le même paquet et sont
   triés: `(2,3,4,5,5,7,8,9)`.

Tri rapide
----------

#. On choisit une valeur `p` dans la liste que l'on appelle pivot

#. On fait des échanges judicieux jusqu'à ce que toutes les valeurs
   strictement plus petites que `p` soient placées avant `p`, et les
   valeurs plus grandes soient placées après.

#. On applique récursivement l'algorithme sur les éléments avant et
   après `p`.

Tri insertion, tri par arbre binaire de recherche
-------------------------------------------------

Analyse de complexité
=====================

.. TOPIC:: Problèmes

    Quelle est le meilleur algorithme de tri?

    Les algorithmes de tris en `O(n\log n)` sont ils optimaux?

.. TOPIC:: Théorème

    Le tri d’une liste de taille `n` est un problème de complexité `O(n\log n)`.

.. TOPIC:: Exercices

    Évaluer au mieux la complexité des problèmes suivants:

    #. Calcul du `n`-ième nombre de Fibonacci;

    #. Calcul du déterminant d’une matrice;

    #. Calcul du rang d’une matrice;

    #. Calcul de l’inverse d’une matrice;

    #. Calcul d’un vecteur `x` solution de `Ax=b`, où
       `A` est une matrice et `b` un vecteur;

    #. Calcul du pgcd de deux nombres;

    #. Test de primalité de `n`;

    #. Recherche du plus court chemin entre deux stations de métro à Paris;

    #. Calcul de la `n`-ième décimale de `\sqrt{2}`;

    #. Calcul de l’inverse d’un nombre modulo `3`;

    #. Recherche d’un échec et mat en `4` coups à partir d’une
       position donnée aux échecs.

    #. Problème du sac-à-dos: étant donné un ensemble d’objets de hauteur et
       de poids variables, et un sac à dos de hauteur donnée, charger au
       maximum le sac-à-dos?

*****************
Travaux pratiques
*****************

Première étude pratique de complexité
=====================================

.. TOPIC:: Exercice

    Implantez une fonction ``recherche(liste, valeur)`` renvoyant la
    première position de ``valeur`` dans la ``liste``. Par exemple::

	sage: recherche([9,20,3,40,37,42,69,65,21,66,1,74,50], 21)
	9
	sage: recherche([9,20,3,40,37,42,69,65,21,66,1,74,50], 69)
	7
	sage: recherche([9,20,3,40,37,42,69,65,21,66,1,74,50], 5)
	None

    Indications: utilisez les tests suivants::

	sage: recherche([],1)
	None
	sage: recherche([2],1)
	None
	sage: recherche([2],2)
	1
	sage: recherche([9,20,3,40,37,42,69,65,21,66,1,74,50],21)
	9
	sage: recherche([9,20,3,40,37,42,69,65,21,66,1,74,50],69)
	7
	sage: recherche([9,20,3,40,37,42,69,65,21,66,1,74,50],5)
	None
	sage: recherche([1,3,9,20,21,37,40,42,50,65,66,69,74],21)
	5
	sage: recherche([1,3,9,20,21,37,40,42,50,65,66,69,74],69)
	13
	sage: recherche([1,3,9,20,21,37,40,42,50,65,66,69,74],5)
	None

    Insérer un compteur dans la fonction pour compter le nombre de
    comparaisons effectuées lors d’un appel. Faire quelques statistiques
    et tracer une courbe donnant le nombre de comparaisons en moyenne et
    au pire en fonction de la taille de la liste.

    .. TODO:: C'est difficile. Donner plus d'indications!

    Indications:

    #. Voir :func:`randint` pour créer une liste aléatoire

    #. Voir :func:`point` pour afficher un nuage de points

	Que fait l'exemple suivant?::

	    sage: point( [ [i, i^2] for i in range(10) ] )


.. TOPIC:: Exercice:

    Même exercice précédement, mais en supposant que les listes sont
    triées et utilisant une recherche dichotomique.

    Indications: Voir ``while``, :func:`sort`, et utiliser deux bornes
    ``inf`` et ``sup``, vérifiant à chaque étape l’invariant ``inf <=
    i < sup``, où ``i`` est la première position (éventuelle) de
    ``valeur`` dans la ``liste``.

    Comparer les deux courbes. Évaluer la taille maximale d’une liste
    dans laquelle on peut faire une recherche en moins d’une heure et
    d’une semaine.

Exercice: complexité de l'algorithme de tri de Python
=====================================================

#. Estimez la complexité de l'algorithme de tri de Python (:func:`sort`)

#. Estimez la complexité de la fonction suivante::

       sage: def fusion(l1, l2):
       ...       sort(l1+l2)

   lorsque elle est appliquée à des listes aléatoires, respectivement triées.

   Qu'en déduisez vous?

   Pour en savoir plus: [TimSort]_

Exercice: Implantation de quelques algorithmes de tri
=====================================================

Écrire quelques tests pour une fonction de tri de liste d’entiers.

Implanter des fonctions de tri utilisant chacun des algorithmes
suivants. Pour chacune tracer des courbes statistiques de complexité
au pire et en moyenne. Comparer avec les courbes théoriques.

#. Tri à bulle en place.

   Indication: *choisir au préalable le bon invariant!*

#. Tri fusion.

   Indication: utiliser une fonction récursive; le cas échéant, s’entraîner en implantant au préalable une fonction récursive pour calculer `n!`

#. Tri rapide en place.

   Indication: *choisir au préalable le bon invariant!*

#. Tri par tas.

   Indication: utiliser le module `heapq <http://docs.python.org/library/heapq.html>`_

#. Tri insertion dans un arbre binaire de recherche (équilibré ou non).

   Indication: définir une fonction récursive ``insert(arbre, i)`` qui
   insère un nombre ``i`` dans un arbre binaire de recherche.

   .. NOTE::

       Sage n'a pas encore de bonne structure de donnée pour les
       arbres binaires; on y travaille! En attendant, on peut
       représenter un arbre binaire récursivement par une liste ``[l,
       t1, t2]`` où ``l`` est le label de la racine, ``t1`` est le
       sous-arbre gauche et ``t2`` le sous-arbre droit, en utilisant
       ``None`` pour dénoter un sous-arbre vide.

*******************
Quelques références
*******************

-  Wikipédia Française: `Complexité algorithmique <http://fr.wikipedia.org/wiki/Complexité_algorithmique>`_

.. -  `Un support de cours sur les tris <http://dept-info.labri.u-bordeaux.fr/~lachaud/IUT/ASD-Prog-E1-2000/planning-prof.html>`_

-  `Une fiche de TP sur les tris <http://www.lri.fr/~denise/M2Spec/97-98.1/TDSpec6.ps>`_

.. -  `Démonstration de bubble sort et quicksort <http://jade.lim.univ-mrs.fr/~vancan/mait/demo/SortDemo/example1.html>`_

.. [TimSort] `Tim sort <http://en.wikipedia.org/wiki/Timsort>`_
.. [CTFM1993] `Constant Time Factor do Matter <http://scholar.google.fr/scholar?hl=fr&q=constant+time+factor+do+matter>`_

