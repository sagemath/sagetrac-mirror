## -*- encoding: utf-8 -*-
"""
This file (./sol/graphique_doctest.sage) was *autogenerated* from ./sol/graphique.tex,
with sagetex.sty version 2011/05/27 v2.3.1.
It contains the contents of all the sageexample environments from this file.
You should be able to doctest this file with:
sage -t ./sol/graphique_doctest.sage
It is always safe to delete this file; it is not used in typesetting your
document.

Sage example in ./sol/graphique.tex, line 3::

  sage: reset()

Sage example in ./sol/graphique.tex, line 10::

  sage: t = var('t'); liste = [a + cos(t) for a in srange(0, 2, 0.1)]
  sage: g = polar_plot(liste, (t, 0, 2 * pi)); g.show(aspect_ratio = 1)

Sage example in ./sol/graphique.tex, line 38::

  sage: f = lambda x: abs(x**2 - 1/4)
  sage: def liste_pts(u0, n):
  ....:     u = u0; liste = [[u0,0]]
  ....:     for k in range(n):
  ....:         v, u = u, f(u)
  ....:         liste.extend([[v,u], [u,u]])
  ....:     return(liste)
  sage: g = line(liste_pts(1.1, 8), rgbcolor = (.9,0,0))
  sage: g += line(liste_pts(-.4, 8), rgbcolor = (.01,0,0))
  sage: g += line(liste_pts(1.3, 3), rgbcolor = (.5,0,0))
  sage: g += plot(f, -1, 3, rgbcolor = 'blue')
  sage: g += plot(x, -1, 3, rgbcolor = 'green')
  sage: g.show(aspect_ratio = 1, ymin = -.2, ymax = 3)

Sage example in ./sol/graphique.tex, line 78::

  sage: x = var('x'); y = function('y')
  sage: DE = x^2 * diff(y(x), x) - y(x) == 0
  sage: desolve(DE, y(x))
  _C*e^(-1/x)
  sage: g = plot([c*e^(-1/x) for c in srange(-8, 8, 0.4)], (x, -3, 3))
  sage: y = var('y')
  sage: g += plot_vector_field((x^2, y), (x,-3,3), (y,-5,5))
  sage: g.show()

Sage example in ./sol/graphique.tex, line 124::

  sage: from sage.calculus.desolvers import desolve_system_rk4
  sage: f = lambda x, y: [a*x-b*x*y,-c*y+d*b*x*y]
  sage: x, y, t = var('x, y, t')
  sage: a, b, c, d = 1., 0.1, 1.5, 0.75
  sage: P = desolve_system_rk4(f(x,y), [x,y],\
  ....:        ics=[0,10,5], ivar=t, end_points=15)
  sage: Ql = [[i,j] for i,j,k in P]; p = line(Ql, color='red')
  sage: p += text("Rabbits", (12,37), fontsize=10, color='red')
  sage: Qr = [[i,k] for i,j,k in P]; p += line(Qr, color='blue')
  sage: p += text("Foxes", (12,7), fontsize=10, color='blue')
  sage: p.axes_labels(["time", "population"])
  sage: p.show(gridlines = True)

Sage example in ./sol/graphique.tex, line 154::

  sage: n = 10;  L = srange(6, 18, 12 / n); R = srange(3, 9, 6 / n)
  sage: def g(x,y): v = vector(f(x, y)); return v / v.norm()
  sage: q = plot_vector_field(g(x, y), (x, 0, 60), (y, 0, 36))
  sage: for j in range(n):
  ....:     P = desolve_system_rk4(f(x,y), [x,y],
  ....:            ics=[0,L[j],R[j]], ivar=t, end_points=15)
  ....:     Q = [[j,k] for i,j,k in P]
  ....:     q += line(Q, color=hue(.8-j/(2*n)))
  sage: q.axes_labels(["rabbits", "foxes"]); q.show()

Sage example in ./sol/graphique.tex, line 185::

  sage: from scipy import integrate
  sage: def dX_dt(X, t=0): return [X[1], 0.5*X[1] - X[0] - X[1]^3]
  sage: t = srange(0, 40, 0.01);  x0 = srange(-2, 2, 0.1); y0 = 2.5
  sage: CI = [[i, y0] for i in x0] + [[i, -y0] for i in x0]
  sage: def g(x,y): v = vector(dX_dt([x, y])); return v / v.norm()
  sage: x, y = var('x, y'); n = len(CI)
  sage: q = plot_vector_field(g(x, y), (x, -3, 3), (y, -y0, y0))
  sage: for j in range(n):
  ....:     X = integrate.odeint(dX_dt, CI[j], t)
  ....:     q += line(X, color=(1.7*j/(4*n),1.5*j/(4*n),1-3*j/(8*n)))
  sage: X = integrate.odeint(dX_dt, [0.01,0], t)
  sage: q += line(X, color = 'red'); q.show()

Sage example in ./sol/graphique.tex, line 234::

  sage: from scipy import integrate
  sage: t = srange(0, 40, 0.2)
  sage: n = 35; CI_cart = [[4, .2 * i] for i in range(n)]
  sage: CI = list(map(lambda x: [sqrt(x[0]^2+x[1]^2),
  ....:      pi - arctan(x[1]/x[0])], CI_cart))
  sage: for alpha in [0.1, 0.5, 1, 1.25]:                      # long time
  ....:     dX_dt = lambda X, t=0: [cos(X[1])*(1-1/X[0]^2),
  ....:             -sin(X[1]) * (1/X[0]+1/X[0]^3) + 2*alpha/X[0]^2]
  ....:     q = circle((0, 0), 1, fill=True, rgbcolor='purple')
  ....:     for j in range(n):
  ....:         X = integrate.odeint(dX_dt, CI[j], t)
  ....:         Y = [[u[0]*cos(u[1]), u[0]*sin(u[1])] for u in X]
  ....:         q += line(Y, xmin = -4, xmax = 4, color='blue')
  ....:     q.show(aspect_ratio = 1, axes = False)

"""
