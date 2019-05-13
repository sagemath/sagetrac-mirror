##      -*-   coding: utf-8   -*-     ##
##          Sage Doctest File         ##
#**************************************#
#*    Generated from PreTeXt source   *#
#*    on 2017-08-24T11:43:34-07:00    *#
#*                                    *#
#*   http://mathbook.pugetsound.edu   *#
#*                                    *#
#**************************************#
##
"""
Please contact Rob Beezer (beezer@ups.edu) with
any test failures here that need to be changed
as a result of changes accepted into Sage.  You
may edit/change this file in any sensible way, so
that development work may procede.  Your changes
may later be replaced by the authors of "Abstract
Algebra: Theory and Applications" when the text is
updated, and a replacement of this file is proposed
for review.
"""
##
## To execute doctests in these files, run
##   $ $SAGE_ROOT/sage -t <directory-of-these-files>
## or
##   $ $SAGE_ROOT/sage -t <a-single-file>
##
## Replace -t by "-tp n" for parallel testing,
##   "-tp 0" will use a sensible number of threads
##
## See: http://www.sagemath.org/doc/developer/doctesting.html
##   or run  $ $SAGE_ROOT/sage --advanced  for brief help
##
## Generated at 2017-08-24T11:43:34-07:00
## From "Abstract Algebra"
## At commit 26d3cac0b4047f4b8d6f737542be455606e2c4b4
##
## Section 11.5 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C12 = CyclicPermutationGroup(12)
    sage: C20 = CyclicPermutationGroup(20)
    sage: domain_gens = C12.gens()
    sage: [g.order() for g in domain_gens]
    [12]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: x = C20.gen(0)
    sage: y = x^5
    sage: y.order()
    4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: phi = PermutationGroupMorphism(C12, C20, [y])
    sage: phi
    Permutation group morphism:
      From: Cyclic group of order 12 as a permutation group
      To:   Cyclic group of order 20 as a permutation group
      Defn: [(1,2,3,4,5,6,7,8,9,10,11,12)] ->
            [(1,6,11,16)(2,7,12,17)(3,8,13,18)(4,9,14,19)(5,10,15,20)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = C12("(1,6,11,4,9,2,7,12,5,10,3,8)")
    sage: phi(a)
    (1,6,11,16)(2,7,12,17)(3,8,13,18)(4,9,14,19)(5,10,15,20)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = C12("(1,3,5,7,9,11)(2,4,6,8,10,12)")
    sage: phi(b)
    (1,11)(2,12)(3,13)(4,14)(5,15)(6,16)(7,17)(8,18)(9,19)(10,20)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c = C12("(1,9,5)(2,10,6)(3,11,7)(4,12,8)")
    sage: phi(c)
    ()

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K = phi.kernel(); K
    Subgroup generated by [(1,5,9)(2,6,10)(3,7,11)(4,8,12)]
    of (Cyclic group of order 12 as a permutation group)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Im = phi.image(C12); Im
    Subgroup generated by
    [(1,6,11,16)(2,7,12,17)(3,8,13,18)(4,9,14,19)(5,10,15,20)]
    of (Cyclic group of order 20 as a permutation group)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Im.is_isomorphic(C12.quotient(K))
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = DihedralGroup(5)
    sage: H = DihedralGroup(20)
    sage: G.gens()
    [(1,2,3,4,5), (1,5)(2,4)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.gens()
    [(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),
     (1,20)(2,19)(3,18)(4,17)(5,16)(6,15)(7,14)(8,13)(9,12)(10,11)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: x = H.gen(0)^4
    sage: y = H.gen(1)
    sage: rho = PermutationGroupMorphism(G, H, [x, y])
    sage: rho.kernel()
    Subgroup generated by
    [()] of (Dihedral group of order 10 as a permutation group)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Im = rho.image(G); Im
    Subgroup generated by
    [(1,5,9,13,17)(2,6,10,14,18)(3,7,11,15,19)(4,8,12,16,20),
     (1,20)(2,19)(3,18)(4,17)(5,16)(6,15)(7,14)(8,13)(9,12)(10,11)]
    of (Dihedral group of order 40 as a permutation group)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Im.is_subgroup(H)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Im.is_isomorphic(G)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = CyclicPermutationGroup(7)
    sage: H = CyclicPermutationGroup(4)
    sage: tau = PermutationGroupMorphism_im_gens(G, H, H.gens())
    sage: tau
    Permutation group morphism:
      From: Cyclic group of order 7 as a permutation group
      To:   Cyclic group of order 4 as a permutation group
      Defn: [(1,2,3,4,5,6,7)] -> [(1,2,3,4)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tau.kernel()
    Traceback (most recent call last):
    ...
    RuntimeError: Gap produced error output
    ...

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = CyclicPermutationGroup(3)
    sage: H = DihedralGroup(4)
    sage: results = G.direct_product(H)
    sage: results[0]
    Permutation Group with generators [(4,5,6,7), (4,7)(5,6), (1,2,3)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: results[1]
    Permutation group morphism:
      From: Cyclic group of order 3 as a permutation group
      To:   Permutation Group with generators
            [(4,5,6,7), (4,7)(5,6), (1,2,3)]
      Defn: Embedding( Group( [ (1,2,3), (4,5,6,7), (4,7)(5,6) ] ), 1 )

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: results[2]
    Permutation group morphism:
      From: Dihedral group of order 8 as a permutation group
      To:   Permutation Group with generators
            [(4,5,6,7), (4,7)(5,6), (1,2,3)]
      Defn: Embedding( Group( [ (1,2,3), (4,5,6,7), (4,7)(5,6) ] ), 2 )

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: results[3]
    Permutation group morphism:
      From: Permutation Group with generators
            [(4,5,6,7), (4,7)(5,6), (1,2,3)]
      To:   Cyclic group of order 3 as a permutation group
      Defn: Projection( Group( [ (1,2,3), (4,5,6,7), (4,7)(5,6) ] ), 1 )

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: results[4]
    Permutation group morphism:
      From: Permutation Group with generators
            [(4,5,6,7), (4,7)(5,6), (1,2,3)]
      To:   Dihedral group of order 8 as a permutation group
      Defn: Projection( Group( [ (1,2,3), (4,5,6,7), (4,7)(5,6) ] ), 2 )

"""
