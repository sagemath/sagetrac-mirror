# -*- coding: utf-8 -*-
r'''
Tutorial of Combinatorial Hopf algebras

Abstract
--------

In this module we implements several ``combinatorial Hopf algebras`` (CHA)
(like ``FQSym`` (the Malvenuto-Reutenauer Hopf algebra, [MalReut]_, [NCSF-VI]_,
[NCSF-VII]_), ``WQSym`` ([PhdHivert]_, [CDendTri]_), ``PQSym`` ([PQS-FPSAC]_,
[PQS-Comp]_), ``CQSym`` ([PQS-Comp]_), ``PBT`` (the Loday-Ronco Hopf algebra,
[LodRon]_, [HNT04]_, [HNT05]_), ``FSym`` (Standard Tableaux
Hopf algebra of Poirier_Reutenauer, [FSym-PR]_, [NCSF-VI]_), etc).

Roughly, we could define a CHA as graded and connexe Hopf algebra with basis
indexed by a combinatorial class. For a formal and concrete definition of a CHA
we refer to [CHA-def]_ or [CHA-def2]_ (or for french people [CHA-def3]_).

Introduction
------------

This tutorial presents how to use theses CHAs and several tools implemented
with likes:

    - product, coproduct,
    - antipode,
    - internal product,
    - scalar product,
    - #-product ([AvaVien]_ and [AvNoThi]_),
    - (bi or tri)-dendriform operations ([CHA-def]_),
    - expand the polynomial realization.

Basic use
---------

This first part of the tutorial, we quickly expose how to use these CHA. If the
reader is familiar with how to use NCSF/QSym (.. SEE ALSO:
:mod:`sage.combinat.ncsf_qsym.tutorial`) this part is useless.

Each of our CHAs are described on differents bases which are `Realizations` in
Sage. So there is the implementation of the abstract CHA (for example `FQSym`)
and the concrete realization on a specific basis on this CHA (for example the
fundamental basis of FQSym).

To call the abstract CHA over a specific ring::

    sage: FQS = FreeQuasiSymmetricFunctions(QQ); FQS
    The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field

And sometime we have shorthands::

    sage: FQS = FQSym(QQ); FQS
    The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field

(It is possible the CHA is not load with Sage so you could have to load
manually::

    sage: from sage.combinat.cha.fqsym import FreeQuasiSymmetricFunctions
    sage: FQS = FreeQuasiSymmetricFunctions(QQ)

)

To call the concrete realization of a CHA in a specific basis, we call the
abstract CHA and select the realization::

    sage: F = FQS.Fundamental()       # FQS.F()
    sage: G = FQS.FundamentalDual()   # FQS.G()
    sage: E = FQS.Elementary()        # FQS.E()
    sage: H = FQS.Homogene()          # FQS.H()
    sage: n = FQS.ElementaryDual()    # FQS.n()
    sage: m = FQS.HomogeneDual()      # FQS.m()

And generally if you use often one CHA for research or just for fun you could
use shorthands::

    sage: FQSym(QQ).inject_shorthands()
    Injecting F as shorthand for The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the Fundamental basis
    Injecting G as shorthand for The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the FundamentalDual basis
    Injecting E as shorthand for The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the Elementary basis
    Injecting n as shorthand for The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the ElementaryDual basis
    ...
    Injecting H as shorthand for The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the Homogene basis
    Injecting m as shorthand for The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the HomogeneDual basis
    sage: F
    The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the Fundamental basis
    sage: n
    The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the ElementaryDual basis

Furthermore quickly you don't support to read the longest description of your
favorite CHA, so you could rename it::

    sage: FQS = FQSym(ZZ); FQS
    The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field
    sage: FQS.rename("FQSym"); FQS
    FQSym
    sage: FQS.reset_name(); FQS
    The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field

Now we have one concrete realization of a basis of our CHA, how to construct
an element::

    sage: F(Permutation([3,1,2]))
    F[3, 1, 2]

As this is rather cumbersome, the following abuses of notation are (often)
allowed::

    sage: F([3,1,2])
    F[3, 1, 2]
    sage: F[[3,1,2]]
    F[3, 1, 2]
    sage: F[3,1,2]
    F[3, 1, 2]

Unfortunately, due to a limitation in Python syntax, one cannot use::

    sage: F[]       # not implemented

Instead, you can use::

    sage: F[[]]
    F[]

Now, we can construct linear combinations of basis elements and compute
product, coproduct, antipode, unit and counit::

    sage: F[2,1,3] + 2 * (F[1,3,4,2] + F[2,1])
    2*F[1, 3, 4, 2] + 2*F[2, 1] + F[2, 1, 3]
    sage: F[3,1,2] * F[2,1]
    F[3, 1, 2, 5, 4] + F[3, 1, 5, 2, 4] + F[3, 1, 5, 4, 2] +
    F[3, 5, 1, 2, 4] + F[3, 5, 1, 4, 2] + F[3, 5, 4, 1, 2] +
    F[5, 3, 1, 2, 4] + F[5, 3, 1, 4, 2] + F[5, 3, 4, 1, 2] +
    F[5, 4, 3, 1, 2]
    sage: F[3,1,2].coproduct()
    F[] # F[3, 1, 2] + F[1] # F[1, 2] + F[2, 1] # F[1] + F[3, 1, 2] # F[]
    sage: F[1,2,3].antipode()
    -F[3, 2, 1]
    sage: F[1,3,2].antipode()
    -F[2, 1, 3] - F[2, 3, 1] + F[3, 1, 2]
    sage: F.one()
    F[]
    sage: (3 + F[3,1,2] + 2*F[1]).counit()
    3

We also could use morphism with other realization of our CHA or with an other::

    sage: G = FreeQuasiSymmetricFunctions.FundamentalDual()   # FQSym(QQ).G()
    sage: M = QuasiSymmetricFunctions(QQ).M()
    sage: G(F[3,1,2])
    G[2, 3, 1]
    sage: M(G[3,1,2])
    M[1, 1, 1] + M[2, 1]

And it works with tensor too::

    sage: F.rename("F-basis")
    sage: G.rename("G-basis")
    sage: FxF = F.tensor_square(); FxF
    F-basis # F-basis
    sage: FxF(G[3,1,4,2].coproduct())
    F[] # F[2, 4, 1, 3] + F[1] # F[3, 1, 2] + F[1, 2] # F[1, 2] + F[2, 3, 1] # F[1] + F[2, 4, 1, 3] # F[]
    sage: FxG = tensor((F,G)); FxG
    F-basis # G-basis
    sage: FxG(G[3,1,4,2].coproduct())
    F[] # G[3, 1, 4, 2] + F[1] # G[2, 3, 1] + F[1, 2] # G[1, 2] + F[2, 3, 1] # G[1] + F[2, 4, 1, 3] # G[]

In the following parts, we expose several tools often present in these CHAs.

Internal Product
----------------

.. TODO

Scalar Product
--------------

.. TODO

#-product
---------

.. TODO

Dendriform Structures
---------------------

.. TODO

Polynomial realizations
-----------------------

.. TODO

References
----------

.. [MalReut] Duality between quasi-symmetrical functions and the solomon
    descent algebra,
    Claudia Malvenuto and
    Christophe Reutenauer

.. [NCSF-VI] Noncommutative Symmetric Function VI: Free Quasi-Symmetric
    Functions and Related Algebras,
    Gérard Duchamp,
    Florent Hivert and
    Jean-Yves Thibon

.. [NCSF-VII] Noncommutative symmetric functions VII : free quasi-symmetric
    functions revisited,
    Gérard Duchamp,
    Florent Hivert,
    Jean-Christophe Novelli and
    Jean-Yves Thibon

.. [PhdHivert] Combinatoire des fonctions quasi-symétriques,
    Florent Hivert

.. [CDendTri] Construction of dendriform trialgebras,
    Jean-Christophe Novelli,
    Jean-Yves Thibon, and al

.. [PQS-FPSAC] A hopf algebra of parking functions,
    Jean-Christophe Novelli and
    Jean-Yves Thibon

.. [PQS-Comp] Parking functions and descent algebras,
    Jean-Christophe Novelli and
    Jean-Yves Thibon

.. [LodRon] Hopf algebra of the planar binary trees,
    Jean-Louis Loday and
    Maria Ronco

.. [HNT04] An analogue of the plactic monoid for binary search trees,
    Florent Hivert,
    Jean-Christophe Novelli and
    Jean-Yves Thibon

.. [HNT05] The algebra of binary search trees,
    Florent Hivert,
    Jean-Christophe Novelli and
    Jean-Yves Thibon

.. [FSym-PR] Algebres de hopf de tableaux,
    Stéphane Poirier and
    Christophe Reutenauer

.. [CHA-def] Combinatorial Hopf algebras,
    Jean-Louis Loday and
    Maria Ronco

.. [CHA-def2] Hopf algebras in combinatorics,
    Victor Reiner

.. [CHA-def3] Algèbres de Hopf combinatoires,
    Loïc Foissy

.. [AvaVien] The product of trees in the Loday-Ronco algebra through Catalan
    alternative tableaux,
    Jean-Christophe Aval and
    Xavier Viennot,

.. [AvNoThi] The # product in combinatorial Hopf algebras,
    Jean-Christophe Aval,
    Jean-Christophe Novelli and
    Jean-Yves Thibon
'''
