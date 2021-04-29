## -*- encoding: utf-8 -*-
"""
This file (./domaines_doctest.sage) was *autogenerated* from ./domaines.tex,
with sagetex.sty version 2011/05/27 v2.3.1.
It contains the contents of all the sageexample environments from this file.
You should be able to doctest this file with:
sage -t ./domaines_doctest.sage
It is always safe to delete this file; it is not used in typesetting your
document.

Sage example in ./domaines.tex, line 10::

  sage: x = var('x')

Sage example in ./domaines.tex, line 69::

  sage: o = 12/35
  sage: type(o)
  <... 'sage.rings.rational.Rational'>

Sage example in ./domaines.tex, line 82::

  sage: type(12/35)
  <... 'sage.rings.rational.Rational'>

Sage example in ./domaines.tex, line 131::

  sage: o = 720
  sage: o.factor()
  2^4 * 3^2 * 5

Sage example in ./domaines.tex, line 142::

  sage: type(o).factor(o)
  2^4 * 3^2 * 5

Sage example in ./domaines.tex, line 157::

  sage: 720.factor()
  2^4 * 3^2 * 5

Sage example in ./domaines.tex, line 166::

  sage: o = 720 / 133
  sage: o.numerator().factor()
  2^4 * 3^2 * 5

Sage example in ./domaines.tex, line 253::

  sage: 3 * 7
  21

Sage example in ./domaines.tex, line 261::

  sage: (2/3) * (6/5)
  4/5

Sage example in ./domaines.tex, line 267::

  sage: (1 + I)  *  (1 - I)
  2

Sage example in ./domaines.tex, line 274::

  sage: (x + 2) * (x + 1)
  (x + 2)*(x + 1)
  sage: (x + 1) * (x + 2)
  (x + 2)*(x + 1)

Sage example in ./domaines.tex, line 308::

  sage: def fourth_power(a):
  ....:      a = a * a
  ....:      a = a * a
  ....:      return a

Sage example in ./domaines.tex, line 330::

  sage: fourth_power(2)
  16
  sage: fourth_power(3/2)
  81/16
  sage: fourth_power(I)
  1
  sage: fourth_power(x+1)
  (x + 1)^4
  sage: M = matrix([[0,-1],[1,0]]); M
  [ 0 -1]
  [ 1  0]
  sage: fourth_power(M)
  [1 0]
  [0 1]

Sage example in ./domaines.tex, line 375::

  sage: t = type(5/1); t
  <... 'sage.rings.rational.Rational'>
  sage: t == type(5)
  False

Sage example in ./domaines.tex, line 476::

  sage: a = 5; a
  5
  sage: a.is_unit()
  False

Sage example in ./domaines.tex, line 484::

  sage: a = 5/1; a
  5
  sage: a.is_unit()
  True

Sage example in ./domaines.tex, line 507::

  sage: parent(5)
  Integer Ring
  sage: parent(5/1)
  Rational Field

Sage example in ./domaines.tex, line 515::

  sage: ZZ
  Integer Ring
  sage: QQ
  Rational Field

Sage example in ./domaines.tex, line 525::

  sage: QQ(5).parent()
  Rational Field
  sage: ZZ(5/1).parent()
  Integer Ring
  sage: ZZ(1/5)
  Traceback (most recent call last):
    ...
  TypeError: no conversion of this rational to integer

Sage example in ./domaines.tex, line 543::

  sage: ZZ(1), QQ(1), RR(1), CC(1)
  (1, 1, 1.00000000000000, 1.00000000000000)

Sage example in ./domaines.tex, line 568::

  sage: cartesian_product([QQ, QQ])
  The Cartesian product of (Rational Field, Rational Field)

Sage example in ./domaines.tex, line 574::

  sage: ZZ.fraction_field()
  Rational Field

Sage example in ./domaines.tex, line 580::

  sage: ZZ['x']
  Univariate Polynomial Ring in x over Integer Ring

Sage example in ./domaines.tex, line 591::

  sage: Z5 = GF(5); Z5
  Finite Field of size 5
  sage: P = Z5['x']; P
  Univariate Polynomial Ring in x over Finite Field of size 5
  sage: M = MatrixSpace(P, 3, 3); M
  Full MatrixSpace of 3 by 3 dense matrices over
  Univariate Polynomial Ring in x over Finite Field of size 5

Sage example in ./domaines.tex, line 602::

  sage: M.random_element()                           # random
  [2*x^2 + 3*x + 4 4*x^2 + 2*x + 2     4*x^2 + 2*x]
  [            3*x   2*x^2 + x + 3     3*x^2 + 4*x]
  [      4*x^2 + 3 3*x^2 + 2*x + 4         2*x + 4]

Sage example in ./domaines.tex, line 697::

  sage: QQ.category()
  Join of Category of number fields and Category of quotient fields and Category of metric spaces

Sage example in ./domaines.tex, line 704::

  sage: QQ in Fields()
  True

Sage example in ./domaines.tex, line 712::

  sage: QQ in CommutativeAdditiveGroups()
  True

Sage example in ./domaines.tex, line 718::

  sage: QQ['x'] in EuclideanDomains()
  True

Sage example in ./domaines.tex, line 859::

  sage: 5.parent()
  Integer Ring

Sage example in ./domaines.tex, line 872::

  sage: type(factor(4))
  <class 'sage.structure.factorization_integer.IntegerFactorization'>

Sage example in ./domaines.tex, line 895::

  sage: int(5)
  5
  sage: type(int(5))
  <... 'int'>

Sage example in ./domaines.tex, line 909::

  sage: Integer(5)
  5
  sage: type(Integer(5))
  <... 'sage.rings.integer.Integer'>

Sage example in ./domaines.tex, line 926::

  sage: factorial(99) / factorial(100) - 1 / 50
  -1/100

Sage example in ./domaines.tex, line 974::

  sage: 72/53 - 5/3 * 2.7
  -3.14150943396227

Sage example in ./domaines.tex, line 982::

  sage: cos(1), cos(1.)
  (cos(1), 0.540302305868140)

Sage example in ./domaines.tex, line 1000::

  sage: pi.n(digits=50)    # variant: n(pi,digits=50)
  3.1415926535897932384626433832795028841971693993751

Sage example in ./domaines.tex, line 1020::

  sage: z = CC(1,2); z.arg()
  1.10714871779409

Sage example in ./domaines.tex, line 1036::

  sage: I.parent()
  Number Field in I with defining polynomial x^2 + 1 with I = 1*I

Sage example in ./domaines.tex, line 1043::

  sage: (1.+2.*I).parent()
  Complex Field with 53 bits of precision
  sage: (1.+2.*SR(I)).parent()
  Symbolic Ring

Sage example in ./domaines.tex, line 1064::

  sage: z = 3 * exp(I*pi/4)
  sage: z.real(), z.imag(), z.abs().canonicalize_radical()
  (3/2*sqrt(2), 3/2*sqrt(2), 3)

Sage example in ./domaines.tex, line 1094::

  sage: a, b, c = 0, 2, 3
  sage: a == 1 or (b == 2 and c == 3)
  True

Sage example in ./domaines.tex, line 1147::

  sage: x, y = var('x, y')
  sage: bool( (x-y)*(x+y) == x^2-y^2 )
  True

Sage example in ./domaines.tex, line 1171::

  sage: Z4 = IntegerModRing(4); Z4
  Ring of integers modulo 4
  sage: m = Z4(7); m
  3

Sage example in ./domaines.tex, line 1184::

  sage: 3 * m + 1
  2

Sage example in ./domaines.tex, line 1191::

  sage: Z3 = GF(3); Z3
  Finite Field of size 3

Sage example in ./domaines.tex, line 1243::

  sage: a = matrix(QQ, [[1,2,3],[2,4,8],[3,9,27]])
  sage: (a^2 + 1) * a^(-1)
  [  -5 13/2  7/3]
  [   7    1 25/3]
  [   2 19/2   27]

Sage example in ./domaines.tex, line 1259::

  sage: M = MatrixSpace(QQ,3,3); M
  Full MatrixSpace of 3 by 3 dense matrices over Rational Field
  sage: a = M([[1,2,3],[2,4,8],[3,9,27]])
  sage: (a^2 + 1) * a^(-1)
  [  -5 13/2  7/3]
  [   7    1 25/3]
  [   2 19/2   27]

Sage example in ./domaines.tex, line 1283::

  sage: P = ZZ['x']; P
  Univariate Polynomial Ring in x over Integer Ring
  sage: F = P.fraction_field(); F
  Fraction Field of Univariate Polynomial Ring in x over Integer Ring
  sage: p = P(x+1) * P(x); p
  x^2 + x
  sage: p + 1/p
  (x^4 + 2*x^3 + x^2 + 1)/(x^2 + x)
  sage: parent(p + 1/p)
  Fraction Field of Univariate Polynomial Ring in x over Integer Ring

Sage example in ./domaines.tex, line 1382::

  sage: k.<a> = NumberField(x^3 + x + 1); a^3; a^4+3*a
  -a - 1
  -a^2 + 2*a

Sage example in ./domaines.tex, line 1416::

  sage: parent(sin(x))
  Symbolic Ring

Sage example in ./domaines.tex, line 1422::

  sage: SR
  Symbolic Ring

Sage example in ./domaines.tex, line 1428::

  sage: SR.category()
  Category of fields

Sage example in ./domaines.tex, line 1482::

  sage: R = QQ['x1,x2,x3,x4']; R
  Multivariate Polynomial Ring in x1, x2, x3, x4 over Rational Field
  sage: x1, x2, x3, x4 = R.gens()

Sage example in ./domaines.tex, line 1489::

  sage: x1 * (x2 - x3)
  x1*x2 - x1*x3

Sage example in ./domaines.tex, line 1496::

  sage: (x1+x2)*(x1-x2) - (x1^2 - x2^2)
  0

Sage example in ./domaines.tex, line 1509::

  sage: P = prod( (a-b) for (a,b) in Subsets([x1,x2,x3,x4],2) ); P * P.lc()
  x1^3*x2^2*x3 - x1^2*x2^3*x3 - x1^3*x2*x3^2 + x1*x2^3*x3^2
  + x1^2*x2*x3^3 - x1*x2^2*x3^3 - x1^3*x2^2*x4 + x1^2*x2^3*x4
  + x1^3*x3^2*x4 - x2^3*x3^2*x4 - x1^2*x3^3*x4 + x2^2*x3^3*x4
  + x1^3*x2*x4^2 - x1*x2^3*x4^2 - x1^3*x3*x4^2 + x2^3*x3*x4^2
  + x1*x3^3*x4^2 - x2*x3^3*x4^2 - x1^2*x2*x4^3 + x1*x2^2*x4^3
  + x1^2*x3*x4^3 - x2^2*x3*x4^3 - x1*x3^2*x4^3 + x2*x3^2*x4^3

Sage example in ./domaines.tex, line 1531::

  sage: x1, x2, x3, x4 = SR.var('x1, x2, x3, x4')
  sage: got = prod( (a-b) for (a,b) in Subsets([x1,x2,x3,x4],2) )
  sage: expected1 = -(x1 - x2)*(x1 - x3)*(x1 - x4)*(x2 - x3)*(x2 - x4)*(x3 - x4)
  sage: expected2 = (x1 - x2)*(x1 - x3)*(x1 - x4)*(x2 - x3)*(x2 - x4)*(x3 - x4)
  sage: bool(got == expected1 or got == expected2)
  True

Sage example in ./domaines.tex, line 1581::

  sage: x = var('x')
  sage: p = 54*x^4+36*x^3-102*x^2-72*x-12
  sage: factor(p)
  6*(x^2 - 2)*(3*x + 1)^2

Sage example in ./domaines.tex, line 1616::

  sage: R = ZZ['x']; R
  Univariate Polynomial Ring in x over Integer Ring

Sage example in ./domaines.tex, line 1622::

  sage: q = R(p); q
  54*x^4 + 36*x^3 - 102*x^2 - 72*x - 12

Sage example in ./domaines.tex, line 1629::

  sage: parent(q)
  Univariate Polynomial Ring in x over Integer Ring

Sage example in ./domaines.tex, line 1635::

  sage: factor(q)
  2 * 3 * (3*x + 1)^2 * (x^2 - 2)

Sage example in ./domaines.tex, line 1642::

  sage: R = QQ['x']; R
  Univariate Polynomial Ring in x over Rational Field
  sage: q = R(p); q
  54*x^4 + 36*x^3 - 102*x^2 - 72*x - 12
  sage: factor(q)
  (54) * (x + 1/3)^2 * (x^2 - 2)

Sage example in ./domaines.tex, line 1665::

  sage: R = ComplexField(16)['x']; R
  Univariate Polynomial Ring in x over Complex Field
  with 16 bits of precision
  sage: q = R(p); q
  54.00*x^4 + 36.00*x^3 - 102.0*x^2 - 72.00*x - 12.00
  sage: factor(q)
  (54.00) * (x - 1.414) * (x + 0.3333)^2 * (x + 1.414)

Sage example in ./domaines.tex, line 1685::

  sage: R = QQ[sqrt(2)]['x']; R
    Univariate Polynomial Ring in x over Number Field in sqrt2 with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?
  sage: q = R(p); q
  54*x^4 + 36*x^3 - 102*x^2 - 72*x - 12
  sage: factor(q)
  (54) * (x - sqrt2) * (x + sqrt2) * (x + 1/3)^2

Sage example in ./domaines.tex, line 1698::

  sage: R = GF(5)['x']; R
  Univariate Polynomial Ring in x over Finite Field of size 5
  sage: q = R(p); q
  4*x^4 + x^3 + 3*x^2 + 3*x + 3
  sage: factor(q)
  (4) * (x + 2)^2 * (x^2 + 3)

"""

