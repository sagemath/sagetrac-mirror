"""
Polygons and triangles in hyperbolic geometry

AUTHORS:

- Hartmut Monien (2011-08)
- Vincent Delecroix (2014-11)
"""
#*****************************************************************************
#       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>,
#                     2014 Vincent Delecroix <20100.delecroix@gmail.com>,
#                     2015 Stefan Kraemer <skraemer@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.plot.bezier_path import BezierPath
from sage.plot.misc import options, rename_keyword
from sage.rings.all import CC
from sage.rings.infinity import infinity
from sage.functions.other import real, imag


class HyperbolicPolygon(BezierPath):
    """
    Primitive class for hyberbolic polygon type.

    See ``hyperbolic_polygon?`` for information about plotting a hyperbolic
    polygon in the complex plane.

    INPUT:

    - ``pts`` -- coordinates of the polygon (as complex numbers)

    - ``options`` -- dict of valid plot options to pass to constructor

    EXAMPLES:

    Note that constructions should use :func:`hyperbolic_polygon` or
    :func:`hyperbolic_triangle`::

         sage: from sage.plot.hyperbolic_polygon import HyperbolicPolygon
         sage: print(HyperbolicPolygon([0, 1/2, I], {}))
         Hyperbolic polygon (0.000000000000000, 0.500000000000000, 1.00000000000000*I)
    """
    def __init__(self, pts, options):
        """
        Initialize HyperbolicPolygon.

        EXAMPLES::

            sage: from sage.plot.hyperbolic_polygon import HyperbolicPolygon
            sage: print(HyperbolicPolygon([0, 1/2, I], {}))
            Hyperbolic polygon (0.000000000000000, 0.500000000000000, 1.00000000000000*I)
        """

        pts = [CC(_) for _ in pts]
        im_list = [imag(_) for _ in pts]
        max_im = max(im_list)
        self.path = []
        if (CC(infinity) in pts):
            idx = pts.index(CC(infinity))
            for i in range(idx, len(pts) - 1):
                self._hyperbolic_arc(pts[i], pts[i + 1], maxY=max_im)
            self._hyperbolic_arc(pts[-1], pts[0], maxY=max_im)
            for i in range(0, idx):
                self._hyperbolic_arc(pts[i], pts[i + 1], maxY=max_im)
        else:
            self._hyperbolic_arc(pts[0], pts[1], first=True)
            for i in range(1, len(pts) - 1):
                self._hyperbolic_arc(pts[i], pts[i + 1])
            self._hyperbolic_arc(pts[-1], pts[0])
        BezierPath.__init__(self, self.path, options)
        self._pts = pts
        self._max_im = max_im

#         pts = [CC(_) for _ in pts]
#         self.path = []
#         self._hyperbolic_arc(pts[0], pts[1], True)
#         for i in range(1, len(pts) - 1):
#             self._hyperbolic_arc(pts[i], pts[i + 1])
#         self._hyperbolic_arc(pts[-1], pts[0])
#         BezierPath.__init__(self, self.path, options)
#         self._pts = pts

    def _max_im_coord(self):
        """
        Function to return the max imaginary coordinate among the vertices.
        """
        return self._max_im

    def _repr_(self):
        """
        String representation of HyperbolicPolygon.

        TESTS::

            sage: from sage.plot.hyperbolic_polygon import HyperbolicPolygon
            sage: HyperbolicPolygon([0, 1/2, I], {})._repr_()
            'Hyperbolic polygon (0.000000000000000, 0.500000000000000, 1.00000000000000*I)'
        """
        return "Hyperbolic polygon ({})".format(", ".join(map(str, self._pts)))

    def _hyperbolic_arc(self, z0, z3, first=False, maxY=0):
        """
        Function to construct Bezier path as an approximation to
        the hyperbolic arc between the complex numbers z0 and z3 in the
        hyperbolic plane.
        """
        z0, z3 = (CC(z0), CC(z3))

        if z0 == CC(infinity):
            z0 = CC(real(z3), maxY + 100)
        else:
            if z3 == CC(infinity):
                z3 = CC(real(z0), maxY + 100)

        p = (abs(z0)*abs(z0)-abs(z3)*abs(z3))/(z0-z3).real()/2
        r = abs(z0-p)

        if abs(z3-z0)/r < 0.1:
            self.path.append([(z0.real(), z0.imag()), (z3.real(), z3.imag())])
            return

        if z0.imag() == 0 and z3.imag() == 0:
            p = (z0.real()+z3.real())/2
            zm = CC(p, r)
            self._hyperbolic_arc(z0, zm, first)
            self._hyperbolic_arc(zm, z3)
            return
        else:
            zm = ((z0+z3)/2-p)/abs((z0+z3)/2-p)*r+p
            t = (8*zm-4*(z0+z3)).imag()/3/(z3-z0).real()
            z1 = z0 + t*CC(z0.imag(), (p-z0.real()))
            z2 = z3 - t*CC(z3.imag(), (p-z3.real()))
        if first:
            self.path.append([(z0.real(), z0.imag()),
                              (z1.real(), z1.imag()),
                              (z2.real(), z2.imag()),
                              (z3.real(), z3.imag())])
            first = False
        else:
            self.path.append([(z1.real(), z1.imag()),
                              (z2.real(), z2.imag()),
                              (z3.real(), z3.imag())])


@rename_keyword(color='rgbcolor')
@options(alpha=1, fill=False, thickness=1, rgbcolor="blue", zorder=2,
         linestyle='solid')
def hyperbolic_polygon(pts, **options):
    r"""
    Return a hyperbolic polygon in the hyperbolic plane with vertices ``pts``.

    Type ``?hyperbolic_polygon`` to see all options.

    INPUT:

    - ``pts`` -- a list or tuple of complex numbers

    OPTIONS:

    - ``alpha`` -- default: 1

    - ``fill`` -- default: ``False``

    - ``thickness`` -- default: 1

    - ``rgbcolor`` -- default: ``'blue'``

    - ``linestyle`` -- (default: ``'solid'``) The style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``,
      or ``'--'``, ``':'``, ``'-'``, ``'-.'``, respectively.

    EXAMPLES:

    Show a hyperbolic polygon with coordinates `-1`, `3i`, `2+2i`, `1+i`::

        sage: hyperbolic_polygon([-1,3*I,2+2*I,1+I])
        Graphics object consisting of 1 graphics primitive

    With more options::

        sage: hyperbolic_polygon([-1,3*I,2+2*I,1+I], fill=True, color='red')
        Graphics object consisting of 1 graphics primitive

    Polygons with ideal vertices are supported::

        sage: P = hyperbolic_polygon([-1,0,+1,oo], fill=True); P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(hyperbolic_polygon([-1,0,+1,oo], fill=True))

    ::

        sage: Q = hyperbolic_triangle(3+2*I,infinity,4+3*I, rgbcolor='orange', fill=True)
        sage: Q
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        Q = hyperbolic_triangle(3+2*I,infinity,4+3*I, rgbcolor='orange', fill=True)
        sphinx_plot(Q)


    When added two hyperbolic polygons with ideal vertices, 'yrange' is
    modified accordingly::

        T = P + Q; T
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_polygon([-1,0,+1,oo], fill=True)
        Q = hyperbolic_triangle(3+2*I,oo,4+3*I, rgbcolor='orange', fill=True)
        T = P + Q
        sphinx_plot(T)

    """
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))
    P = HyperbolicPolygon(pts, options)
    g.add_primitive(P)
    g.set_aspect_ratio(1)
    minmaxData = g.get_axes_range()
    g.set_axes_range(xmin=minmaxData['xmin'],
                     xmax=minmaxData['xmax'],
                     ymin=minmaxData['ymin'],
                     ymax=P._max_im_coord()+1)
    return g


def hyperbolic_triangle(a, b, c, **options):
    """
    Return a hyperbolic triangle in the hyperbolic plane with vertices
    ``(a,b,c)``.

    Type ``?hyperbolic_polygon`` to see all options.

    INPUT:

    - ``a, b, c`` -- complex numbers in the upper half complex plane

    OPTIONS:

    - ``alpha`` -- default: 1

    - ``fill`` -- default: ``False``

    - ``thickness`` -- default: 1

    - ``rgbcolor`` -- default: ``'blue'``

    - ``linestyle`` - (default: ``'solid'``) The style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``,
      or ``'--'``, ``':'``, ``'-'``, ``'-.'``, respectively.

    EXAMPLES:

    Show a hyperbolic triangle with coordinates `0, 1/2+i\sqrt{3}/2` and
    `-1/2+i\sqrt{3}/2`::

         sage: hyperbolic_triangle(0, -1/2+I*sqrt(3)/2, 1/2+I*sqrt(3)/2)
         Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_triangle(0, 0.5*(-1 + I*sqrt(3)), 0.5*(1+I*sqrt(3)))
        sphinx_plot(P)

    A hyperbolic triangle with coordinates `0, 1` and `2+i` and a dashed line::

         sage: hyperbolic_triangle(0, 1, 2+i, fill=true, rgbcolor='red', linestyle='--')
         Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_triangle(0, 1, 2+i, fill=true, rgbcolor='red', linestyle='--')
        sphinx_plot(P)

    A hyperbolic triangle with vertices `I, 1, \infty`::

        sage: P = hyperbolic_triangle(I,1,oo, fill=True);P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_triangle(I,1,oo, fill=True)
        sphinx_plot(P)

    """
    return hyperbolic_polygon((a, b, c), **options)
