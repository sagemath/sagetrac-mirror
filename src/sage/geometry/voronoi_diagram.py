r"""
Voronoi diagram

This module provides the class :class: `voronoi_diagram` for computing the
Voronoi diagram of a finite list of points in \RR^d.
"""

#*****************************************************************************
#       Copyright (C) 2012 Moritz Firsching <moritz@math.fu-berlin.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.all import RDF, QQ, AA
from sage.rings.real_mpfr import RealField_class
from sage.geometry.triangulation.point_configuration import PointConfiguration
from sage.modules.all import vector
from sage.plot.all import line, point, rainbow, plot


class voronoi_diagram(SageObject):
    r"""
    base class for the  Voronoi diagram.
    Computes the Voronoi diagram of a list of points

    INPUT:

    - ``points`` a list of points. Any valid input for the
    :class:`PointConfiguration` will do.

    EXAMPLES:
        Get the Voronoi diagram for some points in \RR^3
        ::
        sage: V=voronoi_diagram([[1,3,.3],[2,-2,1],[-1,2,-.1]]); V
        The Voronoi diagram of 3 points of dimension 3 in the Real Double Field

        sage: voronoi_diagram([])
        The empty Voronoi diagram.

    ALGORITHM:

    We use hyperplanes tangent to the paraboloid one dimension higher to
    get a convex polyhedron and then project back to one dimension lower

    See for example [M2002]

    REFERENCES:

    ..  [M2002]
        Jiri Matousek,
        "Lectures on Discrete Geometry", Springer, Ch.5.7, p.118.

    AUTHORS:

    - Moritz Firsching (2012-09-21)
    """
    def __init__(self, points):
        r"""
        See ``voronoi_diagram`` for full documentation.
        EXAMPLES::
        sage: V=voronoi_diagram([[1,3,3],[2,-2,1],[-1,2,-1]]); V
        The Voronoi diagram of 3 points of dimension 3 in the Rational Field
        """
        self._P={}
        self._points=PointConfiguration(points)
        self._n=self._points.n_points()
        if self._n==0 or self._points.base_ring().is_subring(QQ):
            self._base_ring=QQ
        elif self._points.base_ring() in [RDF,AA]:
            self._base_ring=self._points.base_ring()
        elif isinstance(self._points.base_ring(), RealField_class):
            self._base_ring=RDF
        else:
            raise NotImplementedError('Base ring of the Voronoi diagram must '
                                    +'be one of '+str(QQ)+', ' + str(RDF)+
                                     ', '+str(RDF)+'. Real Field will be converted '
                                     +'to '+str(RDF))


        if self._n > 0:
            self._d=self._points.ambient_dim()
            e=[([sum(vector(i)[k]**2 for k in
            range(self._d))]+[(-2)*vector(i)[l] for l in range(self._d)]+[1])
            for i in self._points]
            e=[[self._base_ring(i) for i in k] for k in e]
            p=Polyhedron(ieqs = e, base_ring=self._base_ring)
        for i in range(self._n):
            equ=p.Hrepresentation(i) #TODO: here we assume that the order of these inequalities is the same as when p was defined two lines above.
            pvert=[[u[k] for k in range(self._d)] for u in equ.incident() if
            u.is_vertex()]
            prays=[[u[k] for k in range(self._d)] for u in equ.incident() if
            u.is_ray()]
            pline=[[u[k] for k in range(self._d)] for u in equ.incident() if
            u.is_line()]
            (self._P)[self._points[i]]=Polyhedron(vertices=pvert, lines=pline, rays=prays,
             base_ring=self._base_ring)

    def points(self):
        r""" Returns the input points (as a PointConfiguration).

        EXAMPLES::

            sage: V=voronoi_diagram([[.5,3],[2,5],[4,5],[4,-1]]); V.points()
            A point configuration in QQ^2 consisting of 4 points. The
            triangulations of this point configuration are assumed to be
            connected, not necessarily
            fine, not necessarily regular.
        """
        return self._points

    def ambient_dim(self):
        r"""
        Returns the ambient dimension of the points.

        EXAMPLES::
            sage: V=voronoi_diagram([[.5,3],[2,5],[4,5],[4,-1]])
            sage: V.ambient_dim()
            2
            sage: V=voronoi_diagram([[1,2,3,4,5,6]]); V.ambient_dim()
            6
        """
        return self._d

    def regions(self):
        r"""
        Returns the Voronoi regions of the Voronoi diagram as a list of
        polyhedra.

        EXAMPLES::
            sage: V=voronoi_diagram([[1,3,.3],[2,-2,1],[-1,2,-.1]]); V.regions()
            {P(1.00000000000000, 3.00000000000000, 0.300000000000000): A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
            P(2.00000000000000, -2.00000000000000, 1.00000000000000): A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 1 vertex, 2 rays, 1 line,
            P(-1.00000000000000, 2.00000000000000, -0.100000000000000): A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 1 vertex, 2 rays, 1 line}
        """
        return self._P

    def base_ring(self):
        r"""
        Returns the base_ring of the regions of the Voronoi diagram.

        EXAMPLES::
            sage: V=voronoi_diagram([[1,3,1],[2,-2,1],[-1,2,1/2]]); V.base_ring()
            Rational Field
            sage: V=voronoi_diagram([[1,3.14],[2,-2/3],[-1,22]]); V.base_ring()
            Real Double Field
            sage: V=voronoi_diagram([[1,3],[2,4]]); V.base_ring()
            Rational Field
        """
        return self._base_ring

    def _repr_(self):
        r"""
        Return a description of the Voronoi diagram.

        EXAMPLES::
            sage: V=voronoi_diagram(polytopes.regular_polygon(3).vertices()); V
            The Voronoi diagram of 3 points of dimension 2 in the Algebraic Real Field
            sage: voronoi_diagram([])
            The empty Voronoi diagram.

        """
        desc = ''
        if self._n:
            desc+= 'The Voronoi diagram of '+str(self._n)
            desc+= ' points of dimension '+str(self.ambient_dim())
            desc+=' in the '+str(self.base_ring())
        else:
            desc+='The empty Voronoi diagram.'

        return desc


    def plot(self, **kwds):
        """
        Return a graphical representation for 2-dimensional Voronoi diagrams.

        INPUT:

        - ``**kwds`` -- optional keyword parameters, passed on as arguments for
        plot().

        OUTPUT:

        A graphics object.

        EXAMPLES::
            sage: P=[[0.671, 0.650],[0.258, 0.767], [0.562, 0.406],\
            [0.254, 0.709], [0.493, 0.879]]

            sage: V=voronoi_diagram(P); S=V.plot()

            sage: show(S,xmin=0,xmax=1,ymin=0,ymax=1, aspect_ratio=1,\
            axes=false)

        Trying to plot a Voronoi diagram of dimension other than 2 gives an
        error::
            sage: voronoi_diagram([[1,2,3],[6,5,4]]).plot()
            Traceback (most recent call last):
            ...
            NotImplementedError: Plotting of 3-dimensional Voronoi diagrams not
                implemented

        """
        if self.ambient_dim()==2:
            S=line([])
            for i,j in enumerate(self._points):
                S+=(self._P[j]).render_solid(color=rainbow(self._n)[i],
                zorder=1)
                S+=point(j, color=rainbow(self._n)[i], pointsize=10,zorder=3)
                S+=point(vector(j), color='black',pointsize=20,zorder=2)
            return plot(S,**kwds)
        raise NotImplementedError('Plotting of '+str(self.ambient_dim())+
                                  '-dimensional Voronoi diagrams not'+
                                  ' implemented')
