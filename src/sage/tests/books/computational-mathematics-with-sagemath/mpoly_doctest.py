## -*- encoding: utf-8 -*-
"""
This file (./mpoly_doctest.sage) was *autogenerated* from ./mpoly.tex,
with sagetex.sty version 2011/05/27 v2.3.1.
It contains the contents of all the sageexample environments from this file.
You should be able to doctest this file with:
sage -t ./mpoly_doctest.sage
It is always safe to delete this file; it is not used in typesetting your
document.

Sage example in ./mpoly.tex, line 65::

  sage: R = PolynomialRing(QQ, 'x,y,z')
  sage: x,y,z = R.gens() # gives the tuples of indeterminates

Sage example in ./mpoly.tex, line 78::

  sage: R = PolynomialRing(QQ, 'x', 10)

Sage example in ./mpoly.tex, line 84::

  sage: x = R.gens()
  sage: sum(x[i] for i in range(5))
  x0 + x1 + x2 + x3 + x4

Sage example in ./mpoly.tex, line 96::

  sage: def test_poly(ring, deg=3):
  ....:     monomials = Subsets(
  ....:         flatten([(x,)*deg for x in (1,) + ring.gens()]),
  ....:         deg, submultiset=True)
  ....:     return add(mul(m) for m in monomials)

Sage example in ./mpoly.tex, line 103::

  sage: test_poly(QQ['x,y']) # py2
  x^3 + x^2*y + x*y^2 + y^3 + x^2 + x*y + y^2 + x + y + 1
  sage: test_poly(QQ['y,x']) # py2
  y^3 + y^2*x + y*x^2 + x^3 + y^2 + y*x + x^2 + y + x + 1
  sage: test_poly(QQ['x,y']) == test_poly(QQ['y,x']) # py2
  True

Sage example in ./mpoly.tex, line 127::

  sage: test_poly(PolynomialRing(QQ, 'x,y', order='deglex')) # py2
  x^3 + x^2*y + x*y^2 + y^3 + x^2 + x*y + y^2 + x + y + 1

Sage example in ./mpoly.tex, line 230::

  sage: R.<x,y> = InfinitePolynomialRing(ZZ, order='lex')
  sage: p = mul(x[k] - y[k] for k in range(2)); p
  x_1*x_0 - x_1*y_0 - x_0*y_1 + y_1*y_0
  sage: p + x[100]
  x_100 + x_1*x_0 - x_1*y_0 - x_0*y_1 + y_1*y_0

Sage example in ./mpoly.tex, line 343::

  sage: R.<x,y,z> = QQ[]
  sage: p = 7*y^2*x^2 + 3*y*x^2 + 2*y*z + x^3 + 6
  sage: p.lt()
  7*x^2*y^2

Sage example in ./mpoly.tex, line 358::

  sage: p[x^2*y] == p[(2,1,0)] == p[2,1,0] == 3
  True

Sage example in ./mpoly.tex, line 366::

  sage: p(0, 3, -1)
  0
  sage: p.subs(x = 1, z = x^2+1)
  2*x^2*y + 7*y^2 + 5*y + 7

Sage example in ./mpoly.tex, line 381::

  sage: print("total={d}    (in x)={dx}    partial={ds}"\
  ....:   .format(d=p.degree(), dx=p.degree(x), ds=p.degrees()))
  total=4    (in x)=3    partial=(3, 2, 1)

Sage example in ./mpoly.tex, line 441::

  sage: R.<x,y> = QQ[]; p = x^2 + y^2; q = x + y
  sage: print("({quo})*({q}) + ({rem}) == {p}".format( \
  ....:         quo=p//q, q=q, rem=p%q, p=p//q*q+p%q))
  (-x + y)*(x + y) + (2*x^2) == x^2 + y^2
  sage: p.mod(q)  # is NOT equivalent to p%q
  2*y^2

Sage example in ./mpoly.tex, line 459::

  sage: R.<x,y> = QQ[exp(2*I*pi/5)][]
  sage: (x^10 + y^5).gcd(x^4 - y^2)
  x^2 + y
  sage: (x^10 + y^5).factor()
  (x^2 + y) * (x^2 + (a^3)*y) * (x^2 + (a^2)*y) * (x^2 + (a)*y) * (x^2 + (-a^3 - a^2 - a - 1)*y)

Sage example in ./mpoly.tex, line 564::

  sage: R.<x,y,z> = QQ[]
  sage: J = R.ideal(x^2 * y * z - 18,
  ....:             x * y^3 * z - 24,
  ....:             x * y * z^4 - 6)

Sage example in ./mpoly.tex, line 575::

  sage: J.dimension()
  0

Sage example in ./mpoly.tex, line 584::

  sage: J.variety() # py2
  [{y: 2, z: 1, x: 3}]
  sage: J.variety() # py3
  [{z: 1, y: 2, x: 3}]

Sage example in ./mpoly.tex, line 596::

  sage: V = J.variety(QQbar)
  sage: len(V)
  17

Sage example in ./mpoly.tex, line 603::

  sage: sorted(V, key=str)[-3:]
  [{z: 0.9324722294043558? + 0.3612416661871530?*I,
   y: -1.700434271459229? - 1.052864325754712?*I,
   x: 1.337215067329615? + 2.685489874065187?*I},
  {z: 0.9324722294043558? - 0.3612416661871530?*I,
   y: -1.700434271459229? + 1.052864325754712?*I,
   x: 1.337215067329615? - 2.685489874065187?*I},
   {z: 1, y: 2, x: 3}]

Sage example in ./mpoly.tex, line 619::

  sage: (xx, yy, zz) = QQbar['x,y,z'].gens()
  sage: sorted([pt[xx].degree() for pt in V])
  [1, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
  16, 16]

Sage example in ./mpoly.tex, line 640::

  sage: Set(tuple(abs(pt[i]) for i in (xx,yy,zz)) for pt in V)
  {(3, 2, 1)}

Sage example in ./mpoly.tex, line 662::

  sage: w = QQbar.zeta(17); w  # primitive root of 1
  0.9324722294043558? + 0.3612416661871530?*I
  sage: Set(pt[zz] for pt in V) == Set(w^i for i in range(17))
  True

Sage example in ./mpoly.tex, line 690::

  sage: set(pt[zz].minpoly() for pt in sorted(V, key=str)[:-1])
    {x^16 + x^15 + x^14 + x^13 + x^12 + x^11 + x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1}

Sage example in ./mpoly.tex, line 706::

  sage: def polar_form(z):
  ....:     rho = z.abs(); rho.simplify()
  ....:     theta = 2 * pi * z.rational_argument()
  ....:     return (SR(rho) * exp(I*theta))
  sage: sorted([tuple(polar_form(pt[i]) for i in [xx,yy,zz])
  ....:  for pt in sorted(V, key=str)[-3:]])
  [(3*e^(6/17*I*pi), 2*e^(-14/17*I*pi), e^(2/17*I*pi)),
   (3*e^(-6/17*I*pi), 2*e^(14/17*I*pi), e^(-2/17*I*pi)),
   (3, 2, 1)]

Sage example in ./mpoly.tex, line 733::

  sage: J.triangular_decomposition()
  [Ideal (z^17 - 1, y - 2*z^10, x - 3*z^3) of Multivariate
  Polynomial Ring in x, y, z over Rational Field]
  sage: J.transformed_basis()
  [z^17 - 1, -2*z^10 + y, -3/4*y^2 + x]

Sage example in ./mpoly.tex, line 909::

  sage: R.<x,y> = QQ[]
  sage: J = R.ideal(x^2 + y^2 - 1, 16*x^2*y^2 - 1)

Sage example in ./mpoly.tex, line 917::

  sage: ybar2 = R.quo(J)(y^2)
  sage: [ybar2^i for i in range(3)]
  [1, ybar^2, ybar^2 - 1/16]
  sage: ((ybar2 + 1)^2).lift()
  3*y^2 + 15/16

Sage example in ./mpoly.tex, line 958::

  sage: u = (16*y^4 - 16*y^2 + 1).lift(J); u
  [16*y^2, -1]
  sage: u[0]*J.0 + u[1]*J.1
  16*y^4 - 16*y^2 + 1

Sage example in ./mpoly.tex, line 968::

  sage: (y^4).mod(J)
  y^2 - 1/16

Sage example in ./mpoly.tex, line 978::

  sage: (y^4).reduce([x^2 + y^2 - 1, 16*x^2*y^2 - 1])
  y^4

Sage example in ./mpoly.tex, line 1047::

  sage: 1 in ideal(x^2+y^2-1, (x-4)^2+y^2-1)
  False

Sage example in ./mpoly.tex, line 1056::

  sage: R(1).lift(ideal(x^2+y^2-1, (x-4)^2+y^2-1, x-y))
  [-1/28*y + 1/14, 1/28*y + 1/14, -1/7*x + 1/7*y + 4/7]

Sage example in ./mpoly.tex, line 1079::

  sage: J1 = (x^2 + y^2 - 1, 16*x^2*y^2 - 1)*R
  sage: J2 = (x^2 + y^2 - 1,  4*x^2*y^2 - 1)*R
  sage: J1.radical() == J1
  True
  sage: J2.radical()
  Ideal (2*y^2 - 1, 2*x^2 - 1) of Multivariate Polynomial
  Ring in x, y over Rational Field
  sage: 2*y^2 - 1 in J2
  False

Sage example in ./mpoly.tex, line 1122::

  sage: C = ideal(x^2 + y^2 - 1); H = ideal(16*x^2*y^2 - 1)
  sage: C + H == J1
  True

Sage example in ./mpoly.tex, line 1151::

  sage: CH = C.intersection(H).quotient(ideal(4*x*y-1)); CH
  Ideal (4*x^3*y + 4*x*y^3 + x^2 - 4*x*y + y^2 - 1) of
  Multivariate Polynomial Ring in x, y over Rational Field
  sage: CH.gen(0).factor()
  (4*x*y + 1) * (x^2 + y^2 - 1)

Sage example in ./mpoly.tex, line 1161::

  sage: H.quotient(C) == H
  True

Sage example in ./mpoly.tex, line 1184::

  sage: [J.dimension() for J in [J1, J2, C, H, H*J2, J1+J2]]
  [0, 0, 1, 1, 1, -1]

Sage example in ./mpoly.tex, line 1285::

  sage: R.<x,y,z> = QQ[]
  sage: J = ideal(2*x+y-2*z, 2*x+2*y+z-1)
  sage: J.elimination_ideal(x)
  Ideal (y + 3*z - 1) of Multivariate Polynomial Ring in x, y, z
  over Rational Field
  sage: J.elimination_ideal([x,y])
  Ideal (0) of Multivariate Polynomial Ring in x, y, z over Rational Field

Sage example in ./mpoly.tex, line 1312::

  sage: R.<x,y> = QQ[]
  sage: J1 = ideal(x^2 + y^2 - 1, 16*x^2*y^2 - 1)

Sage example in ./mpoly.tex, line 1328::

  sage: g = J1.elimination_ideal(y).gens(); g
  [16*x^4 - 16*x^2 + 1]
  sage: SR(g[0]).solve(SR(x)) # solves by radicals
  [x == -1/2*sqrt(sqrt(3) + 2), x == 1/2*sqrt(sqrt(3) + 2),
  x == -1/2*sqrt(-sqrt(3) + 2), x == 1/2*sqrt(-sqrt(3) + 2)]

Sage example in ./mpoly.tex, line 1351::

  sage: C.elimination_ideal(y).gens()
  [0]
  sage: H.elimination_ideal(y).gens()
  [0]

Sage example in ./mpoly.tex, line 1431::

  sage: R.<x,y,t> = QQ[]
  sage: Param = R.ideal((1-t^2)-(1+t^2)*x, 2*t-(1+t^2)*y)

Sage example in ./mpoly.tex, line 1437::

  sage: Param.elimination_ideal(t).gens()
  [x^2 + y^2 - 1]

Sage example in ./mpoly.tex, line 1469::

  sage: R.<x,y,t> = QQ[]
  sage: eq = x^2 + (y-t)^2 - 1/2*(t^2+1)
  sage: fig = add((eq(t=k/5)*QQ[x,y]).plot() for k in (-15..15))
  sage: fig.show(aspect_ratio=1, xmin=-2, xmax=2, ymin=-3, ymax=3)

Sage example in ./mpoly.tex, line 1494::

  sage: env = ideal(eq, eq.derivative(t)).elimination_ideal(t)
  sage: env.gens()
  [2*x^2 - 2*y^2 - 1]

Sage example in ./mpoly.tex, line 1502::

  sage: env.change_ring(QQ[x,y]).plot((x,-2,2),(y,-3,3))
  Graphics object consisting of 1 graphics primitive

Sage example in ./mpoly.tex, line 1538::

  sage: R.<x,y,t> = QQ[]
  sage: J = (y-t*x, y-t*(1-x))*R
  sage: (x^2+y^2) - ((1-x)^2+y^2) in J
  False

Sage example in ./mpoly.tex, line 1553::

  sage: R.<x,y,t,u> = QQ[]
  sage: J = (y-t*x, y-t*(1-x), t*u-1)*R
  sage: (x^2+y^2) - ((1-x)^2+y^2) in J
  True

Sage example in ./mpoly.tex, line 1615::

  sage: R.<x,y,t> = QQ[]

Sage example in ./mpoly.tex, line 1618::

  sage: eq.derivative(t).resultant(eq, t)
  x^2 - y^2 - 1/2

Sage example in ./mpoly.tex, line 1636::

  sage: R.<x,y> = QQ[]
  sage: p = y^2 - x; q = y^2 + x
  sage: p.resultant(q, y)
  4*x^2
  sage: ideal(p, q).elimination_ideal(y)
  Ideal (x) of Multivariate Polynomial Ring in x, y over Rational Field

Sage example in ./mpoly.tex, line 1665::

  sage: R.<x,y> = QQ[]
  sage: ((x^2 + y^2)*(x^2 + y^2 + 1)*R).dimension()
  1

Sage example in ./mpoly.tex, line 1690::

  sage: R.<x,y> = QQ[]
  sage: J1 = (x^2 + y^2 - 1, 16*x^2*y^2 - 1)*R
  sage: J1.variety()
  []

Sage example in ./mpoly.tex, line 1705::

  sage: J1.variety(QQbar)[0:2]
  [{y: -0.9659258262890683?, x: -0.2588190451025208?},
  {y: -0.9659258262890683?, x: 0.2588190451025208?}]

Sage example in ./mpoly.tex, line 1759::

  sage: R.<x,y> = PolynomialRing(QQ, order='lex')
  sage: C = ideal(x^2+y^2-1)
  sage: D = ideal((x+y-1)*(x+y+1))
  sage: J = C + D

Sage example in ./mpoly.tex, line 1786::

  sage: J.triangular_decomposition()
  [Ideal (y, x^2 - 1) of Multivariate Polynomial Ring in x, y
  over Rational Field,
  Ideal (y^2 - 1, x) of Multivariate Polynomial Ring in x, y
  over Rational Field]

Sage example in ./mpoly.tex, line 1840::

  sage: D = ideal((x+2*y-1)*(x+2*y+1)); J = C + D
  sage: J.variety()
  [{y: 0, x: 1}, {y: 0, x: -1}, {y: 4/5, x: -3/5}, {y: -4/5, x: 3/5}]
  sage: [T.gens() for T in J.triangular_decomposition()]
  [[y, x^2 - 1], [25*y^2 - 16, 4*x + 3*y]]

Sage example in ./mpoly.tex, line 1855::

  sage: Jy = J.elimination_ideal(x); Jy.gens()
  [25*y^3 - 16*y]

Sage example in ./mpoly.tex, line 1863::

  sage: ys = QQ['y'](Jy.0).roots(); ys
  [(4/5, 1), (0, 1), (-4/5, 1)]
  sage: QQ['x'](J.1(y=ys[0][0])).roots()
  [(-3/5, 1), (-13/5, 1)]

Sage example in ./mpoly.tex, line 1882::

  sage: ys = CDF['y'](Jy.0).roots(); ys # abs tol 2e-15
  [(-0.8, 1), (0.0, 1), (0.8, 1)]
  sage: [CDF['x'](p(y=ys[0][0])).roots() for p in J.gens()] # abs tol 2e-15
  [[(-0.5999999999999999, 1), (0.6000000000000001, 1)], [(0.6000000000000001, 1), (2.600000000000001, 1)]]

Sage example in ./mpoly.tex, line 1911::

  sage: R.<x,y> = QQ[]; J = ideal([ x^7-(100*x-1)^2, y-x^7+1 ])
  sage: J.variety(RealField(51))
  [{y: 396340.890166545, x: -14.1660266425312}]

Sage example in ./mpoly.tex, line 1923::

  sage: J.variety(AA) # py2
  [{x: 0.00999999900000035?, y: -0.999999999999990?},
  {x: 0.01000000100000035?, y: -0.999999999999990?},
  {x: 6.305568998641385?, y: 396340.8901665450?}]
  sage: J.variety(AA) # py3
  [{y: -0.999999999999990?, x: 0.00999999900000035?},
  {y: -0.999999999999990?, x: 0.01000000100000035?},
  {y: 396340.8901665450?, x: 6.305568998641385?}]


Sage example in ./mpoly.tex, line 1983::

  sage: len(J2.variety(QQbar)), J2.vector_space_dimension()
  (4, 8)

Sage example in ./mpoly.tex, line 1993::

  sage: J2.normal_basis()
  [x*y^3, y^3, x*y^2, y^2, x*y, y, x, 1]

Sage example in ./mpoly.tex, line 2187::

  sage: R.<x,y,z,t> = PolynomialRing(QQ, order='lex')

Sage example in ./mpoly.tex, line 2243::

  sage: ((x+y+z)^2).reduce([x-t, y-t^2, z^2-t])
  2*z*t^2 + 2*z*t + t^4 + 2*t^3 + t^2 + t

Sage example in ./mpoly.tex, line 2298::

  sage: R.<x,y> = PolynomialRing(QQ, order='lex')
  sage: (g, h) = (x-y, x-y^2);  p = x*y - x
  sage: p.reduce([g, h])  # two reductions by h
  y^3 - y^2
  sage: p.reduce([h, g])  # two reductions by g
  y^2 - y

Sage example in ./mpoly.tex, line 2311::

  sage: p - y*g + h
  0

Sage example in ./mpoly.tex, line 2575::

  sage: R.<x,y> = PolynomialRing(QQ, order='lex')
  sage: R.ideal(x*y^4, x^2*y^3, x^4*y, x^5).basis_is_groebner()
  True

Sage example in ./mpoly.tex, line 2584::

  sage: R.ideal(x^2+y^2-1, 16*x^2*y^2-1).basis_is_groebner()
  False

Sage example in ./mpoly.tex, line 2610::

  sage: R.ideal(x^2+y^2-1, 16*x^2*y^2-1).groebner_basis()
  [x^2 + y^2 - 1, y^4 - y^2 + 1/16]

Sage example in ./mpoly.tex, line 2618::

  sage: R.ideal(16*x^2*y^2-1).groebner_basis()
  [x^2*y^2 - 1/16]

Sage example in ./mpoly.tex, line 2626::

  sage: R.ideal(x^2+y^2-1, (x+y)^2-1).groebner_basis()
  [x^2 + y^2 - 1, x*y, y^3 - y]

Sage example in ./mpoly.tex, line 2636::

  sage: R_lex.<x,y> = PolynomialRing(QQ, order='lex')
  sage: J_lex = (x*y+x+y^2+1, x^2*y+x*y^2+1)*R_lex; J_lex.gens()
  [x*y + x + y^2 + 1, x^2*y + x*y^2 + 1]
  sage: J_lex.groebner_basis()
  [x - 1/2*y^3 + y^2 + 3/2, y^4 - y^3 - 3*y - 1]

Sage example in ./mpoly.tex, line 2644::

  sage: R_invlex = PolynomialRing(QQ, 'x,y', order='invlex')
  sage: J_invlex = J_lex.change_ring(R_invlex); J_invlex.gens()
  [y^2 + x*y + x + 1, x*y^2 + x^2*y + 1]
  sage: J_invlex.groebner_basis()
  [y^2 + x*y + x + 1, x^2 + x - 1]

Sage example in ./mpoly.tex, line 2651::

  sage: R_drl = PolynomialRing(QQ, 'x,y', order='degrevlex')
  sage: J_drl = J_lex.change_ring(R_drl); J_drl.gens()
  [x*y + y^2 + x + 1, x^2*y + x*y^2 + 1]
  sage: J_drl.groebner_basis()
  [y^3 - 2*y^2 - 2*x - 3, x^2 + x - 1, x*y + y^2 + x + 1]

Sage example in ./mpoly.tex, line 2719::

  sage: p = (x + y)^5
  sage: J_lex.reduce(p)
  17/2*y^3 - 12*y^2 + 4*y - 49/2

Sage example in ./mpoly.tex, line 2726::

  sage: p.reduce(J_lex.groebner_basis())
  17/2*y^3 - 12*y^2 + 4*y - 49/2

Sage example in ./mpoly.tex, line 2732::

  sage: R_lex.quo(J_lex)(p)
  17/2*ybar^3 - 12*ybar^2 + 4*ybar - 49/2

Sage example in ./mpoly.tex, line 2738::

  sage: R_drl.quo(J_drl)(p)
  5*ybar^2 + 17*xbar + 4*ybar + 1

Sage example in ./mpoly.tex, line 2751::

  sage: J_lex.normal_basis()
  [y^3, y^2, y, 1]
  sage: J_invlex.normal_basis()
  [x*y, y, x, 1]
  sage: J_drl.normal_basis()
  [y^2, y, x, 1]

Sage example in ./mpoly.tex, line 2775::

  sage: ideal(16*x^2*y^2-1).dimension()
  1

Sage example in ./mpoly.tex, line 2851::

  sage: R.<t,x,y,z> = PolynomialRing(QQ, order='lex')
  sage: J = R.ideal(t+x+y+z-1, t^2-x^2-y^2-z^2-1, t-x*y)
  sage: [u.polynomial(u.variable(0)) for u in J.groebner_basis()]
  [t + x + y + z - 1,
  (y + 1)*x + y + z - 1,
  (z - 2)*x + y*z - 2*y - 2*z + 1,
  (z - 2)*y^2 + (-2*z + 1)*y - z^2 + z - 1]

Sage example in ./mpoly.tex, line 2970::

  sage: from sage.rings.ideal import Cyclic
  sage: Cyclic(QQ['x,y,z'])
  Ideal (x + y + z, x*y + x*z + y*z, x*y*z - 1) of
  Multivariate Polynomial Ring in x, y, z over Rational Field

Sage example in ./mpoly.tex, line 2980::

  sage: def C(R, n): return Cyclic(PolynomialRing(R, 'x', n))

Sage example in ./mpoly.tex, line 3010::

  sage: p = previous_prime(2^30)
  sage: len(C(GF(p), 6).groebner_basis())
  45

"""

