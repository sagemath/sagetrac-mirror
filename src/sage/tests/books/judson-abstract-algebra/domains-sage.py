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
## Section 18.5 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q = ZZ.fraction_field(); Q
    Rational Field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q == QQ
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: R.<x> = ZZ[]
    sage: P = R.fraction_field();P
    Fraction Field of Univariate Polynomial Ring in x over Integer Ring

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: f = P((x^2+3)/(7*x+4))
    sage: g = P((4*x^2)/(3*x^2-5*x+4))
    sage: h = P((-2*x^3+4*x^2+3)/(x^2+1))
    sage: ((f+g)/h).numerator()
    3*x^6 + 23*x^5 + 32*x^4 + 8*x^3 + 41*x^2 - 15*x + 12

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: ((f+g)/h).denominator()
    -42*x^6 + 130*x^5 - 108*x^4 + 63*x^3 - 5*x^2 + 24*x + 48

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.<c> = FiniteField(3^5)
    sage: F.characteristic()
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = F.prime_subfield(); G
    Finite Field of size 3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.list()
    [0, 1, 2]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K.<y>=QuadraticField(-7); K
    Number Field in y with defining polynomial x^2 + 7 with y = 2.645751311064591?*I

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K.prime_subfield()
    Rational Field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K.<x> = ZZ[sqrt(-3)]; K
    Order in Number Field in a with defining polynomial x^2 + 3 with a = 0.?e-18 + 1.732050807568878?*I

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K.is_integral_domain()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K.basis()
    [1, a]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: x
    a

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (1+x)*(1-x) == 2*2
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: four = K(4)
    sage: four.is_unit()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: four^-1
    1/4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: T.<x>=ZZ[]
    sage: T.is_integral_domain()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: J = T.ideal(5, x); J
    Ideal (5, x) of Univariate Polynomial Ring in x over Integer Ring

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q = T.quotient(J); Q
    Quotient of Univariate Polynomial Ring in x over
    Integer Ring by the ideal (5, x)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: J.is_principal()
    Traceback (most recent call last):
    ...
    NotImplementedError

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.is_field()
    Traceback (most recent call last):
    ...
    NotImplementedError

"""
