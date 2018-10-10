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
## Section 4.7 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = 3*ZZ
    sage: -12 in G
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: 37 in G
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.gen()
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = AdditiveAbelianGroup([14])
    sage: G.order()
    14

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.list()
    [(0), (1), (2), (3), (4), (5), (6), (7),
     (8), (9), (10), (11), (12), (13)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = G.gen(0)
    sage: a
    (1)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a + a
    (2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a + a + a + a
    (4)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: 4*a
    (4)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: 37*a
    (9)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G([2])
    doctest:...: DeprecationWarning: The default behaviour changed! If you
    *really* want a linear combination of smith generators, use
    .linear_combination_of_smith_form_gens.
    See http://trac.sagemath.org/16261 for details.
    (2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = G([2]); b
    (2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b + b
    (4)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: 2*b == 4*a
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: 7*b
    (0)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b.order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c = a - 6*b; c
    (3)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c + c + c + c
    (12)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c.order()
    14

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H = G.submodule([b]); H
    Additive abelian group isomorphic to Z/7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.list()
    [(0), (2), (4), (6), (8), (10), (12)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: e = H.gen(0); e
    (2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: 3*e
    (6)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: e.order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: f = 12*a; f
    (12)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: f.order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K = G.submodule([f]); K
    Additive abelian group isomorphic to Z/7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K.order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K.list()
    [(0), (2), (4), (6), (8), (10), (12)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K.gen(0)
    (2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H == K
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.<a> = AbelianGroup([14])
    sage: G.order()
    14

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.list()
    (1, a, a^2, a^3, a^4, a^5, a^6, a^7, a^8, a^9, a^10, a^11, a^12, a^13)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a.order()
    14

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = a^2
    sage: b.order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b*b*b
    a^6

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c = a^7
    sage: c.order()
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c^2
    1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b*c
    a^9

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b^37*c^42
    a^4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H = G.subgroup([a^2])
    sage: H.order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K = G.subgroup([a^12])
    sage: K.order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: allsg = G.subgroups(); allsg
    [Multiplicative Abelian subgroup isomorphic to C2 x C7 generated by {a},
     Multiplicative Abelian subgroup isomorphic to C7 generated by {a^2},
     Multiplicative Abelian subgroup isomorphic to C2 generated by {a^7},
     Trivial Abelian subgroup]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sub = allsg[2]
    sage: sub.order()
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G=CyclicPermutationGroup(14)
    sage: a = G.gen(0); a
    (1,2,3,4,5,6,7,8,9,10,11,12,13,14)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = a^2
    sage: b = a^2; b
    (1,3,5,7,9,11,13)(2,4,6,8,10,12,14)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b.order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a*a*b*b*b
    (1,9,3,11,5,13,7)(2,10,4,12,6,14,8)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c = a^37*b^26; c
    (1,6,11,2,7,12,3,8,13,4,9,14,5,10)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c.order()
    14

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H = G.subgroup([a^2])
    sage: H.order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.gen(0)
    (1,3,5,7,9,11,13)(2,4,6,8,10,12,14)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.list()
    [(),
     (1,13,11,9,7,5,3)(2,14,12,10,8,6,4),
     (1,11,7,3,13,9,5)(2,12,8,4,14,10,6),
     (1,9,3,11,5,13,7)(2,10,4,12,6,14,8),
     (1,7,13,5,11,3,9)(2,8,14,6,12,4,10),
     (1,5,9,13,3,7,11)(2,6,10,14,4,8,12),
     (1,3,5,7,9,11,13)(2,4,6,8,10,12,14)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.<a> = AbelianGroup([14])
    sage: G.cayley_table()
    *  a b c d e f g h i j k l m n
     +----------------------------
    a| a b c d e f g h i j k l m n
    b| b c d e f g h i j k l m n a
    c| c d e f g h i j k l m n a b
    d| d e f g h i j k l m n a b c
    e| e f g h i j k l m n a b c d
    f| f g h i j k l m n a b c d e
    g| g h i j k l m n a b c d e f
    h| h i j k l m n a b c d e f g
    i| i j k l m n a b c d e f g h
    j| j k l m n a b c d e f g h i
    k| k l m n a b c d e f g h i j
    l| l m n a b c d e f g h i j k
    m| m n a b c d e f g h i j k l
    n| n a b c d e f g h i j k l m

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K.<b> = AbelianGroup([10])
    sage: K.cayley_table(names='elements')
      *    1   b b^2 b^3 b^4 b^5 b^6 b^7 b^8 b^9
       +----------------------------------------
      1|   1   b b^2 b^3 b^4 b^5 b^6 b^7 b^8 b^9
      b|   b b^2 b^3 b^4 b^5 b^6 b^7 b^8 b^9   1
    b^2| b^2 b^3 b^4 b^5 b^6 b^7 b^8 b^9   1   b
    b^3| b^3 b^4 b^5 b^6 b^7 b^8 b^9   1   b b^2
    b^4| b^4 b^5 b^6 b^7 b^8 b^9   1   b b^2 b^3
    b^5| b^5 b^6 b^7 b^8 b^9   1   b b^2 b^3 b^4
    b^6| b^6 b^7 b^8 b^9   1   b b^2 b^3 b^4 b^5
    b^7| b^7 b^8 b^9   1   b b^2 b^3 b^4 b^5 b^6
    b^8| b^8 b^9   1   b b^2 b^3 b^4 b^5 b^6 b^7
    b^9| b^9   1   b b^2 b^3 b^4 b^5 b^6 b^7 b^8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = CyclotomicField(14)
    sage: w = G.gen(0); w
    zeta14

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: wc = CDF(w)
    sage: wc.abs()
    1.0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: wc.arg()/N(2*pi/14)
    1.0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = w^2
    sage: b.multiplicative_order()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: bc = CDF(b); bc
    0.62348980185... + 0.781831482468...*I

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: bc.abs()
    1.0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: bc.arg()/N(2*pi/14)
    2.0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sg = [b^i for i in range(7)]; sg
    [1, zeta14^2, zeta14^4,
    zeta14^5 - zeta14^4 + zeta14^3 - zeta14^2 + zeta14 - 1,
    -zeta14, -zeta14^3, -zeta14^5]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c = sg[3]; d = sg[5]
    sage: c*d
    zeta14^2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c = sg[3]; d = sg[6]
    sage: c*d in sg
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c*d == sg[2]
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sg[5]*sg[6] == sg[4]
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.multiplication_table(elements=sg)
    *  a b c d e f g
     +--------------
    a| a b c d e f g
    b| b c d e f g a
    c| c d e f g a b
    d| d e f g a b c
    e| e f g a b c d
    f| f g a b c d e
    g| g a b c d e f

"""
