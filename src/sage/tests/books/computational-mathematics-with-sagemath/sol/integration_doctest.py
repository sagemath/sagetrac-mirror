## -*- encoding: utf-8 -*-
"""
This file (./sol/integration_doctest.sage) was *autogenerated* from ./sol/integration.tex,
with sagetex.sty version 2011/05/27 v2.3.1.
It contains the contents of all the sageexample environments from this file.
You should be able to doctest this file with:
sage -t ./sol/integration_doctest.sage
It is always safe to delete this file; it is not used in typesetting your
document.

Sage example in ./sol/integration.tex, line 22::

  sage: x = var('x')
  sage: def NCRule(n):
  ....:     P = prod([x - j for j in range(n)])
  ....:     return [integrate(P / (x-i), x, 0, n-1) \
  ....:             / (P/(x-i)).subs(x=i) for i in range(n)]

Sage example in ./sol/integration.tex, line 35::

  sage: def QuadNC(f, a, b, n):
  ....:     W = NCRule(n)
  ....:     ret = 0
  ....:     for i in range(n):
  ....:         ret += f(a + (b-a)/(n-1)*i) * W[i]
  ....:     return (b-a)/(n-1)*ret

Sage example in ./sol/integration.tex, line 49::

  sage: QuadNC(lambda u: 1, 0, 1, 12)
  1
  sage: N(QuadNC(sin, 0, pi, 10))
  1.99999989482634

Sage example in ./sol/integration.tex, line 59::

  sage: numerical_integral(x * log(1+x), 0, 1)
  (0.25, 2.7755575615628914e-15)
  sage: N(QuadNC(lambda x: x * log(1+x), 0, 1, 19)) # abs tol 2e-14
  0.250000000000001
  sage: numerical_integral(sqrt(1-x^2), 0, 1) # abs tol 2e-14
  (0.7853981677264822, 9.042725224536535e-07)
  sage: N(pi/4)
  0.785398163397448
  sage: N(QuadNC(lambda x: sqrt(1-x^2), 0, 1, 20))
  0.784586419900198

Sage example in ./sol/integration.tex, line 74::

  sage: [N(QuadNC(lambda x: x * log(1+x), 0, 1, n) - 1/4) # abs tol 3e-16
  ....:  for n in [2, 8, 16]]
  [0.0965735902799726, 1.17408932933522e-7, 2.13449050101566e-13]
  sage: [N(QuadNC(lambda x: sqrt(1-x^2), 0, 1, n) - pi/4)
  ....:  for n in [2, 8, 16]]
  [-0.285398163397448, -0.00524656673640445, -0.00125482109302663]

"""

