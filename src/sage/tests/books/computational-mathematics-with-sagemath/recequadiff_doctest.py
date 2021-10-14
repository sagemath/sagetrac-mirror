## -*- encoding: utf-8 -*-
"""
This file (./recequadiff_doctest.sage) was *autogenerated* from ./recequadiff.tex,
with sagetex.sty version 2011/05/27 v2.3.1.
It contains the contents of all the sageexample environments from this file.
You should be able to doctest this file with:
sage -t ./recequadiff_doctest.sage
It is always safe to delete this file; it is not used in typesetting your
document.

Sage example in ./recequadiff.tex, line 182::

  sage: x = var('x')
  sage: y = function('y')(x)

Sage example in ./recequadiff.tex, line 295::

  sage: x = var('x'); y = function('y')(x)

Sage example in ./recequadiff.tex, line 298::

  sage: desolve(diff(y,x) + 3*y == exp(x), y, show_method=True)
  [1/4*(4*_C + e^(4*x))*e^(-3*x), 'linear']

Sage example in ./recequadiff.tex, line 320::

  sage: desolve(y*diff(y,x) == x, y, show_method=True)
  [1/2*y(x)^2 == 1/2*x^2 + _C, 'separable']

Sage example in ./recequadiff.tex, line 332::

  sage: desolve(diff(y,x) == exp(x+y), y, show_method=True)
  [-(e^(x + y(x)) + 1)*e^(-y(x)) == _C, 'exact']

Sage example in ./recequadiff.tex, line 350::

  sage: desolve(diff(y,x)-y == x*y^4, y, show_method=True)
  [e^x/(-1/3*(3*x - 1)*e^(3*x) + _C)^(1/3), 'bernoulli']

Sage example in ./recequadiff.tex, line 371::

  sage: desolve(x^2*diff(y,x) == y^2+x*y+x^2, y, show_method=True)
  [_C*x == e^(arctan(y(x)/x)), 'homogeneous']

Sage example in ./recequadiff.tex, line 399::

  sage: desolve(diff(y,x) == (cos(y)-2*x)/(y+x*sin(y)), y,
  ....:         show_method=True)
  [x^2 - x*cos(y(x)) + 1/2*y(x)^2 == _C, 'exact']

Sage example in ./recequadiff.tex, line 431::

  sage: desolve(diff(y,x) == x*y^2+y/x-1/x^2, y,
  ....:         contrib_ode=True, show_method=True)[1]
  'riccati'

Sage example in ./recequadiff.tex, line 459::

  sage: desolve(y == x*diff(y,x)-diff(y,x)^2, y,
  ....:         contrib_ode=True, show_method=True)
  [[y(x) == -_C^2 + _C*x, y(x) == 1/4*x^2], 'clairault']

Sage example in ./recequadiff.tex, line 487::

  sage: x = var('x'); y = function('y')(x)

Sage example in ./recequadiff.tex, line 491::

  sage: DE = diff(y,x)+2*y == x**2-2*x+3
  sage: desolve(DE, y)
  1/4*((2*x^2 - 2*x + 1)*e^(2*x) - 2*(2*x - 1)*e^(2*x) + 4*_C
  + 6*e^(2*x))*e^(-2*x)

Sage example in ./recequadiff.tex, line 499::

  sage: desolve(DE, y).expand()
  1/2*x^2 + _C*e^(-2*x) - 3/2*x + 9/4

Sage example in ./recequadiff.tex, line 518::

  sage: desolve(DE, y, show_method=True)[1]
  'linear'

Sage example in ./recequadiff.tex, line 526::

  sage: desolve(DE, y, ics=[0,1]).expand()
  1/2*x^2 - 3/2*x - 5/4*e^(-2*x) + 9/4

Sage example in ./recequadiff.tex, line 541::

  sage: x = var('x'); y = function('y')(x)
  sage: desolve(diff(y,x)*log(y) == y*sin(x), y, show_method=True)
  [1/2*log(y(x))^2 == _C - cos(x), 'separable']

Sage example in ./recequadiff.tex, line 552::

  sage: ed = desolve(diff(y,x)*log(y) == y*sin(x), y); ed
  1/2*log(y(x))^2 == _C - cos(x)

Sage example in ./recequadiff.tex, line 565::

  sage: solve(ed, y)
  [y(x) == e^(-sqrt(2*_C - 2*cos(x))), y(x) == e^(sqrt(2*_C - 2*cos(x)))]

Sage example in ./recequadiff.tex, line 581::

  sage: solve(ed, y)[0].substitute(_C==5).rhs()
  Traceback (most recent call last):
  ...
  NameError: name '_C' is not defined

Sage example in ./recequadiff.tex, line 593::

  sage: ed.variables()
  (_C, x)

Sage example in ./recequadiff.tex, line 604::

  sage: c = ed.variables()[0]
  sage: solve(ed, y)[0].substitute(c == 5).rhs()
  e^(-sqrt(-2*cos(x) + 10))

Sage example in ./recequadiff.tex, line 624::

  sage: plot(solve(ed, y)[0].substitute(c == 2).rhs(), x, -3, 3)
  Graphics object consisting of 1 graphics primitive

Sage example in ./recequadiff.tex, line 640::

  sage: P = Graphics()
  sage: for k in range(1,20,2):
  ....:     P += plot(solve(ed, y)[0].substitute(c==1+k/4).rhs(), x, -3, 3)

Sage example in ./recequadiff.tex, line 660::

  sage: P = Graphics()
  sage: for j in [0,1]:
  ....:     for k in range(1,10,2):
  ....:         f = solve(ed,y)[j].substitute(c==2+0.25*k).rhs()
  ....:         P += plot(f, x, -3, 3)
  sage: P
  Graphics object consisting of 10 graphics primitives

Sage example in ./recequadiff.tex, line 731::

  sage: u = function('u')(x)
  sage: y = x*u
  sage: DE = x*diff(y,x) == y + sqrt(x**2 + y**2)

Sage example in ./recequadiff.tex, line 746::

  sage: forget()

Sage example in ./recequadiff.tex, line 750::

  sage: assume(x>0)
  sage: desolve(DE, u)
  x == _C*e^arcsinh(u(x))

Sage example in ./recequadiff.tex, line 773::

  sage: S = desolve(DE,u)._maxima_().ev(logarc=True).sage().solve(u); S
  [u(x) == -(sqrt(u(x)^2 + 1)*_C - x)/_C]

Sage example in ./recequadiff.tex, line 792::

  sage: solu = (x-S[0]*c)^2; solu
  (_C*u(x) - x)^2 == (u(x)^2 + 1)*_C^2
  sage: sol = solu.solve(u); sol
  [u(x) == -1/2*(_C^2 - x^2)/(_C*x)]

Sage example in ./recequadiff.tex, line 801::

  sage: y(x) = x*sol[0].rhs(); y(x)
  -1/2*(_C^2 - x^2)/_C

Sage example in ./recequadiff.tex, line 820::

  sage: c = y(x).variables()[0]
  sage: P = Graphics()
  sage: for k in range(-19,19,2):
  ....:     P += plot(y(x).substitute(c == 1/k), x, 0, 3)
  sage: P
  Graphics object consisting of 19 graphics primitives

Sage example in ./recequadiff.tex, line 870::

  sage: x = var('x'); y = function('y')(x); a, b = var('a, b')
  sage: DE = diff(y,x) - a*y == -b*y**2
  sage: sol = desolve(DE,[y,x]); sol
  -(log(b*y(x) - a) - log(y(x)))/a == _C + x

Sage example in ./recequadiff.tex, line 880::

  sage: Sol = solve(sol, y)[0]; Sol
  log(y(x)) == _C*a + a*x + log(b*y(x) - a)

Sage example in ./recequadiff.tex, line 892::

  sage: Sol(x) = Sol.lhs()-Sol.rhs(); Sol(x)
  -_C*a - a*x - log(b*y(x) - a) + log(y(x))
  sage: Sol = Sol.simplify_log(); Sol(x)
  -_C*a - a*x + log(y(x)/(b*y(x) - a))
  sage: solve(Sol, y)[0].simplify()
  y(x) == a*e^(_C*a + a*x)/(b*e^(_C*a + a*x) - 1)

Sage example in ./recequadiff.tex, line 928::

  sage: x = var('x'); y = function('y')(x)
  sage: DE = diff(y,x,2)+3*y == x^2-7*x+31
  sage: desolve(DE, y).expand()
  1/3*x^2 + _K2*cos(sqrt(3)*x) + _K1*sin(sqrt(3)*x) - 7/3*x + 91/9

Sage example in ./recequadiff.tex, line 938::

  sage: desolve(DE, y, ics=[0,1,2]).expand()
  1/3*x^2 + 13/9*sqrt(3)*sin(sqrt(3)*x) - 7/3*x - 82/9*cos(sqrt(3)*x) + 91/9

Sage example in ./recequadiff.tex, line 948::

  sage: desolve(DE, y, ics=[0,1,-1,0]).expand()
  1/3*x^2 - 7/3*x - 82/9*cos(sqrt(3))*sin(sqrt(3)*x)/sin(sqrt(3))
  + 115/9*sin(sqrt(3)*x)/sin(sqrt(3)) - 82/9*cos(sqrt(3)*x) + 91/9

Sage example in ./recequadiff.tex, line 1022::

  sage: x, t = var('x, t'); f = function('f')(x); g = function('g')(t)
  sage: z = f*g
  sage: eq(x,t) = diff(z,x,2) == diff(z,t); eq(x,t)
  g(t)*diff(f(x), x, x) == f(x)*diff(g(t), t)

Sage example in ./recequadiff.tex, line 1042::

  sage: eqn = eq/z; eqn(x,t)
  diff(f(x), x, x)/f(x) == diff(g(t), t)/g(t)

Sage example in ./recequadiff.tex, line 1061::

  sage: k = var('k')
  sage: eq1(x,t) = eqn(x,t).lhs() == k; eq2(x,t) = eqn(x,t).rhs() == k

Sage example in ./recequadiff.tex, line 1069::

  sage: g(t) = desolve(eq2(x,t),[g,t]); g(t)
  _C*e^(k*t)

Sage example in ./recequadiff.tex, line 1092::

  sage: assume(k>0); desolve(eq1,[f,x])
  _K1*e^(sqrt(k)*x) + _K2*e^(-sqrt(k)*x)

Sage example in ./recequadiff.tex, line 1176::

  sage: x, s = var('x, s'); f = function('f')(x)
  sage: f(x) = sin(x); f.laplace(x,s)
  x |--> 1/(s^2 + 1)

Sage example in ./recequadiff.tex, line 1191::

  sage: X(s) = 1/(s^2-3*s-4)/(s^2+1) + (s-4)/(s^2-3*s-4)
  sage: X(s).inverse_laplace(s, x)
  3/34*cos(x) + 1/85*e^(4*x) + 9/10*e^(-x) - 5/34*sin(x)

Sage example in ./recequadiff.tex, line 1206::

  sage: X(s).partial_fraction()
  1/34*(3*s - 5)/(s^2 + 1) + 9/10/(s + 1) + 1/85/(s - 4)

Sage example in ./recequadiff.tex, line 1221::

  sage: x = var('x'); y = function('y')(x)
  sage: eq = diff(y,x,x) - 3*diff(y,x) - 4*y - sin(x) == 0
  sage: desolve_laplace(eq, y)
  1/85*(17*y(0) + 17*D[0](y)(0) + 1)*e^(4*x) + 1/10*(8*y(0)
  - 2*D[0](y)(0) - 1)*e^(-x) + 3/34*cos(x) - 5/34*sin(x)
  sage: desolve_laplace(eq, y, ics=[0,1,-1])
  3/34*cos(x) + 1/85*e^(4*x) + 9/10*e^(-x) - 5/34*sin(x)

Sage example in ./recequadiff.tex, line 1279::

  sage: x = var('x'); y1 = function('y1')(x)
  sage: y2 = function('y2')(x); y3 = function('y3')(x)
  sage: y = vector([y1, y2, y3])
  sage: A = matrix([[2,-2,0],[-2,0,2],[0,2,2]])
  sage: system = [diff(y[i], x) - (A * y)[i] for i in range(3)]
  sage: desolve_system(system, [y1, y2, y3], ics=[0,2,1,-2])
  [y1(x) == e^(4*x) + e^(-2*x),
   y2(x) == -e^(4*x) + 2*e^(-2*x),
   y3(x) == -e^(4*x) - e^(-2*x)]

Sage example in ./recequadiff.tex, line 1328::

  sage: x = var('x'); y1 = function('y1')(x); y2 = function('y2')(x)
  sage: y = vector([y1,y2])
  sage: A = matrix([[3,-4],[1,3]])
  sage: system = [diff(y[i], x) - (A * y)[i] for i in range(2)]
  sage: desolve_system(system, [y1, y2], ics=[0,2,0])
  [y1(x) == 2*cos(2*x)*e^(3*x), y2(x) == e^(3*x)*sin(2*x)]

Sage example in ./recequadiff.tex, line 1385::

  sage: x = var('x'); u1 = function('u1')(x); u2 = function('u2')(x)
  sage: u3 = function('u3')(x); u4 = function('u4')(x)
  sage: u = vector([u1,u2,u3,u4])
  sage: A = matrix([[0,0,1,0],[0,0,0,1],[2,-6,1,3],[-2,6,1,-1]])
  sage: system = [diff(u[i], x) - (A*u)[i] for i in range(4)]
  sage: sol = desolve_system(system, [u1, u2, u3, u4])

Sage example in ./recequadiff.tex, line 1399::

  sage: sol[0]
  u1(x) == 1/12*(2*u1(0) - 6*u2(0) + 5*u3(0) + 3*u4(0))*e^(2*x)
  + 1/24*(2*u1(0) - 6*u2(0) - u3(0) + 3*u4(0))*e^(-4*x) + 3/4*u1(0)
  + 3/4*u2(0) - 3/8*u3(0) - 3/8*u4(0)
  sage: sol[1]
  u2(x) == -1/12*(2*u1(0) - 6*u2(0) - u3(0) - 3*u4(0))*e^(2*x)
  - 1/24*(2*u1(0) - 6*u2(0) - u3(0) + 3*u4(0))*e^(-4*x) + 1/4*u1(0)
  + 1/4*u2(0) - 1/8*u3(0) - 1/8*u4(0)

Sage example in ./recequadiff.tex, line 1557::

  sage: x = var('x'); f = function('f')(x)
  sage: f(x) = 3.83*x*(1 - x/100000)
  sage: def u(n):
  ....:     if n==0: return(20000)
  ....:     else: return f(u(n-1))

Sage example in ./recequadiff.tex, line 1568::

  sage: def v(n):
  ....:     V = 20000;
  ....:     for k in [1..n]:
  ....:         V = f(V)
  ....:     return V

Sage example in ./recequadiff.tex, line 1582::

  sage: def cloud(u,n):
  ....:     L = [[0,u(0)]];
  ....:     for k in [1..n]:
  ....:         L += [[k,u(k)]]
  ....:     points(L).show()

Sage example in ./recequadiff.tex, line 1598::

  sage: cloud(u,50)

Sage example in ./recequadiff.tex, line 1619::

  sage: def snail(f,x,u0,n,xmin,xmax):
  ....:     u = u0
  ....:     P = plot(x, x, xmin, xmax, color='gray')
  ....:     for i in range(n):
  ....:         P += line([[u,u],[u,f(u)],[f(u),f(u)]], color = 'red')
  ....:         u = f(u)
  ....:     P += f.plot(x, xmin, xmax, color='blue')
  ....:     P.show()

Sage example in ./recequadiff.tex, line 1637::

  sage: f(x) = 3.83*x*(1 - x/100000)
  sage: snail(f,x,20000,100,0,100000)

Sage example in ./recequadiff.tex, line 1687::

  sage: from sympy import Function
  sage: from sympy.abc import n
  sage: u = Function('u')

Sage example in ./recequadiff.tex, line 1700::

  sage: f = u(n+2)-(3/2)*u(n+1)+(1/2)*u(n)

Sage example in ./recequadiff.tex, line 1707::

  sage: from sympy import rsolve
  sage: rsolve(f, u(n), {u(0):-1,u(1):1}) == 3 - 4*2**(-n)
  True

Sage example in ./recequadiff.tex, line 1798::

  sage: from sympy import rsolve_hyper
  sage: from sympy.abc import n
  sage: rsolve_hyper([-2,1],2**(n+2),n)
  2**n*C0 + 2**(n + 2)*(C0 + n/2)

"""

