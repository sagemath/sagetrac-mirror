r"""
Combinatorics quickref
----------------------

Integer Sequences::

    sage: s = oeis([1,3,19,211]); s                  # optional - internet
    0: A000275: Coefficients of a Bessel function (reciprocal of J_0(z)); also pairs of permutations with rise/rise forbidden.
    sage: s[0].programs() # optional - internet
    0: (PARI) {a(n) = if( n<0, 0, n!^2 * 4^n * polcoeff( 1 / besselj(0, x + x * O(x^(2*n))), 2*n))}; /* _Michael Somos_, May 17 2004 */

Combinatorial objects::

    sage: S = Subsets([1,2,3,4]); S.list(); S.<tab>        # not tested
    sage: P = Partitions(10000); P.cardinality()
    3616...315650422081868605887952568754066420592310556052906916435144
    sage: Combinations([1,3,7]).random_element()           # random
    sage: Compositions(5, max_part = 3).unrank(3)
    [2, 2, 1]

    sage: DyckWord([1,0,1,0,1,1,0,0]).to_binary_tree()
    [., [., [[., .], .]]]
    sage: Permutation([3,1,4,2]).robinson_schensted()
    [[[1, 2], [3, 4]], [[1, 3], [2, 4]]]
    sage: StandardTableau([[1, 4], [2, 5], [3]]).schuetzenberger_involution()
    [[1, 3], [2, 4], [5]]

Constructions and Species::

    sage: for (p, s) in cartesian_product([P,S]): print((p, s)) # not tested
    sage: DisjointUnionEnumeratedSets(Family(lambda n: IntegerVectors(n, 3), NonNegativeIntegers))  # not tested

Words::

    sage: Words('abc', 4).list()
    [word: aaaa, ..., word: cccc]

    sage: Word('aabcacbaa').is_palindrome()
    True
    sage: WordMorphism('a->ab,b->a').fixed_point('a')
    word: abaababaabaababaababaabaababaabaababaaba...

Polytopes::

    sage: points = random_matrix(ZZ, 6, 3, x=7).rows()
    sage: L = LatticePolytope(points)
    sage: L.npoints(); L.plot3d()                         # random

:ref:`Root systems, Coxeter and Weyl groups <sage.combinat.root_system>`::

    sage: WeylGroup(["B",3]).bruhat_poset()
    Finite poset containing 48 elements
    sage: RootSystem(["A",2,1]).weight_lattice().plot()   # not tested

:ref:`Crystals <sage.combinat.crystals>`::

    sage: CrystalOfTableaux(["A",3], shape = [3,2]).some_flashy_feature() # not tested

:mod:`Symmetric functions and combinatorial Hopf algebras <sage.combinat.algebraic_combinatorics>`::

    sage: Sym = SymmetricFunctions(QQ); Sym.inject_shorthands()
    ...
    doctest:...: RuntimeWarning: redefining global value `e`
    sage: m( ( h[2,1] * (1 + 3 * p[2,1]) ) + s[2](s[3]) )
    3*m[1, 1, 1] + ... + 10*m[5, 1] + 4*m[6]

:ref:`Discrete groups, Permutation groups <sage.groups.groups_catalog>`::

    sage: S = SymmetricGroup(4)
    sage: M = PolynomialRing(QQ, 'x0,x1,x2,x3')
    sage: M.an_element() * S.an_element()
    x1

Graph theory, posets, lattices (:ref:`sage.graphs`, :ref:`sage.combinat.posets`)::

    sage: Poset({1: [2,3], 2: [4], 3: [4]}).linear_extensions().cardinality()
    2
"""
