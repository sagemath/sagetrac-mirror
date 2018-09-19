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
## Section 5.4 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = SymmetricGroup(5)
    sage: sigma = G("(1,3)(2,5,4)")
    sage: sigma*sigma
    (2,4,5)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: rho = G("(2,4)(1,5)")
    sage: rho^3
    (1,5)(2,4)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigma*rho
    (1,3,5,2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: rho*sigma
    (1,4,5,3)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: rho^-1*sigma*rho
    (1,2,4)(3,5)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigma1 = G("(1,3)(2,5,4)")
    sage: sigma1
    (1,3)(2,5,4)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigma2 = G([(1,3),(2,5,4)])
    sage: sigma2
    (1,3)(2,5,4)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigma3 = G([3,5,1,2,4])
    sage: sigma3
    (1,3)(2,5,4)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigma1 == sigma2
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigma2 == sigma3
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigma2.cycle_tuples()
    [(1, 3), (2, 5, 4)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [sigma3(x) for x in G.domain()]
    [3, 5, 1, 2, 4]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H = SymmetricGroup(4)
    sage: sigma = H("(1,2,3,4)")
    sage: G = SymmetricGroup(6)
    sage: tau = G("(1,2,3,4,5,6)")
    sage: rho = tau * sigma
    sage: rho
    (1,3)(2,4,5,6)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigma.parent()
    Symmetric group of order 4! as a permutation group

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tau.parent()
    Symmetric group of order 6! as a permutation group

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: rho.parent()
    Symmetric group of order 6! as a permutation group

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tau.parent() == rho.parent()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigmaG = G(sigma)
    sage: sigmaG.parent()
    Symmetric group of order 6! as a permutation group

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tauH = H(tau)
    Traceback (most recent call last):
    ...
    ValueError: Invalid permutation vector: (1,2,3,4,5,6)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: groups.permutation.   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D = DihedralGroup(5)
    sage: elements = D.list(); elements
    [(),
     (1,5,4,3,2),
     (1,4,2,5,3),
     (1,3,5,2,4),
     (1,2,3,4,5),
     (2,5)(3,4),
     (1,5)(2,4),
     (1,4)(2,3),
     (1,3)(4,5),
     (1,2)(3,5)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: rotate = elements[4]
    sage: flip = elements[5]
    sage: flip*rotate == rotate* flip
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D = DihedralGroup(5)
    sage: D.is_abelian()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A4 = AlternatingGroup(4)
    sage: A4.order()
    12

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A4.is_finite()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A4.is_abelian()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A4.is_cyclic()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigma = A4("(1,2,4)")
    sage: sigma^-1
    (1,4,2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sigma.order()
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = SymmetricGroup(3)
    sage: sigma = G("(1,2)")
    sage: tau = G("(1,3)")
    sage: rho = sigma*tau
    sage: sigma.sign()
    -1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: rho.sign()
    1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A4 = AlternatingGroup(4)
    sage: sigma = A4("(1,2,4)")
    sage: sg = A4.subgroup([sigma])
    sage: sg
    Subgroup of (Alternating group of order 4!/2 as a permutation group) 
    generated by [(1,2,4)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sg.order()
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sg.list()
    [(), (1,4,2), (1,2,4)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sg.is_abelian()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sg.is_cyclic()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sg.is_subgroup(A4)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = SymmetricGroup(5)
    sage: sigma = G("(4,5)")
    sage: tau = G("(1,3)")
    sage: H = G.subgroup([sigma, tau])
    sage: H.list()
    [(), (4,5), (1,3), (1,3)(4,5)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: text_names = ['id', 'sigma', 'tau', 'mu']
    sage: H.cayley_table(names=text_names)
        *     id sigma   tau    mu
         +------------------------
       id|    id sigma   tau    mu
    sigma| sigma    id    mu   tau
      tau|   tau    mu    id sigma
       mu|    mu   tau sigma    id

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = SymmetricGroup(8)
    sage: above = G("(1,2,3,4)(5,6,7,8)")
    sage: front = G("(1,4,8,5)(2,3,7,6)")
    sage: right = G("(1,2,6,5)(3,7,8,4)")
    sage: cube = G.subgroup([above, front, right])
    sage: cube.order()
    24

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: cube.list()
    [(),
     (1,3)(2,4)(5,7)(6,8),
     (1,6)(2,5)(3,8)(4,7),
     (1,8)(2,7)(3,6)(4,5),
     (1,4,3,2)(5,8,7,6),
     (1,2,3,4)(5,6,7,8),
     (1,5)(2,8)(3,7)(4,6),
     (1,7)(2,6)(3,5)(4,8),
     (2,5,4)(3,6,8),
     (1,3,8)(2,7,5),
     (1,6,3)(4,5,7),
     (1,8,6)(2,4,7),
     (1,4)(2,8)(3,5)(6,7),
     (1,2,6,5)(3,7,8,4),
     (1,5,6,2)(3,4,8,7),
     (1,7)(2,3)(4,6)(5,8),
     (2,4,5)(3,8,6),
     (1,3,6)(4,7,5),
     (1,6,8)(2,7,4),
     (1,8,3)(2,5,7),
     (1,4,8,5)(2,3,7,6),
     (1,2)(3,5)(4,6)(7,8),
     (1,5,8,4)(2,6,7,3),
     (1,7)(2,8)(3,4)(5,6)]
"""
