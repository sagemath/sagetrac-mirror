r"""
Hyperplane Arrangements in the Space Of Weights

To compute the cohomology of toric line bundles and Klyachko bundles,
one needs to find a subset of weights on which the cohomology is
supported. Note that the weights are the points of the dual lattice
$m\in M$. This amounts to finding the integral points in the compact
regions of a certain hyperplane arrangement in $M_\QQ$.

You should always use the
:meth:`~sage.schemes.toric.sheaf.klyachko.weight_arrangement` method
to construct instances of :class:`WeightArrangement`.

EXAMPLES::

    sage: V = toric_varieties.P1().sheaves.tangent_bundle()
    sage: w = V.weight_arrangement();  w
    Arrangement <-2*t0 + 3 | -2*t0 + 5 | 2*t0 + 3 | 2*t0 + 5>
    sage: type(w)
    <class 'sage.geometry.hyperplane_arrangement.arrangement.WeightArrangement_with_category'>
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of 
#  the License, or (at your option) any later version.  
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.misc.cachefunc import cached_method
from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangementElement


def make_WeightArrangement(hyperplane_arrangement):
    """
    Construct a :class:`WeightArrangement`

    This function is for internal use only.

    INPUT:

    - ``hyperplane_arrangement`` -- a hyperplane arrangement.

    OUTPUT:

    The same hyperplane arrangement, but with its class hierarchy
    extended by :class:`WeightArrangement`.

    EXAMPLES::

        sage: H.<x,y> = HyperplaneArrangements(QQ)
        sage: arrangement = (x | x+1 | x+2 | x+y-3 | x-y+4)
        sage: from sage.schemes.toric.sheaf.weight_arrangement import make_WeightArrangement
        sage: w = make_WeightArrangement(arrangement);  w
        Arrangement of 5 hyperplanes of dimension 2 and rank 2
        sage: type(w)
        <class 'sage.geometry.hyperplane_arrangement.arrangement.WeightArrangement_with_category'>
        sage: w.__class__.__mro__
        (<class '....arrangement.WeightArrangement_with_category'>,
         <class '....arrangement.HyperplaneArrangements_with_category.element_class'>, 
         ...
         <class 'sage.schemes.toric.sheaf.weight_arrangement.WeightArrangement'>, <type 'object'>)
    """
    parent = hyperplane_arrangement.parent()
    if parent.base_ring().characteristic() != 0:
        raise ValueError('characteristic must be zero')
    from sage.structure.dynamic_class import dynamic_class
    bases = [parent.element_class, WeightArrangement]
    cls = dynamic_class('WeightArrangement_with_category', bases, cache=False)
    hyperplane_arrangement.__class__ = cls
    return hyperplane_arrangement

    

class WeightArrangement(object):
    """
    Hyperplane arrangement in the space of weights.

    .. note::

        By construction, there are no integral points on the
        hyperplanes. See
        :meth:`~sage.schemes.toric.sheaf.klyachko.weight_arrangement`
        for details.
    """

    def convex_hull(self):
        """
        Return the convex hull of the vertices
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(self.vertices(exclude_sandwiched=True))

    @cached_method
    def hull_points_by_sign(self):
        """
        Return sets of weights with the same cohomology groups.

        OUTPUT:

        The integral points in the convex hull, 

        EXAMPLES::

            sage: X = toric_varieties.P2()
            sage: L = X.sheaves.line_bundle(-X.K())
            sage: H = L.weight_arrangement()
            sage: pts = H.hull_points_by_sign()
            sage: for sign_vector in sorted(pts.keys()):
            ....:     print sign_vector, ':  ', sorted(pts[sign_vector])
            (-1, 1, -1, 1, 1, 1) :   [(4, -2)]
            (-1, 1, 1, 1, -1, 1) :   [(-2, 4)]
            (-1, 1, 1, 1, 1, 1) :   [(-1, 3), (0, 2), (1, 1), (2, 0), (3, -1)]
            (1, 1, -1, 1, -1, 1) :   [(-2, -2)]
            (1, 1, -1, 1, 1, 1) :   [(-1, -2), (0, -2), (1, -2), (2, -2), (3, -2)]
            (1, 1, 1, 1, -1, 1) :   [(-2, -1), (-2, 0), (-2, 1), (-2, 2), (-2, 3)]
            (1, 1, 1, 1, 1, 1) :   [(-1, -1), (-1, 0), (-1, 1), (-1, 2), (0, -1), 
                                    (0, 0), (0, 1), (1, -1), (1, 0), (2, -1)]
        """
        points_by_sign = dict()
        for p in self.convex_hull().integral_points():
            p_sign = self.sign_vector(p)
            points = points_by_sign.get(p_sign, [])
            points.append(p)
            points_by_sign[p_sign] = points
        return points_by_sign


    @cached_method
    def bounded_regions_points(self):
        """
        Return points in each bounded region.

        OUTPUT:

        A tuple containing tuples of points. Each inner tuple is a
        collection of integral points in one of the bounded regions.

            sage: X = toric_varieties.dP7()
            sage: K = X.sheaves.line_bundle(X.K())
            sage: TX = X.sheaves.tangent_bundle()
            sage: H = (TX * K).weight_arrangement()
            sage: H.plot()

        The following is the correct number of integral points in the
        bounded regions, but it is rather slow for complicated
        hyperplane arrangements::

            sage: sum(len(region.integral_points()) for region in H.bounded_regions()) 
            11

        There are 11 points in the bounded regions, just looking at
        the points in the convex hull would yield too many points::

            sage: len(H.convex_hull().integral_points())
            13

        This method starts with the convex hull and removes the points
        that are in non-compact regions that reach into the convex
        hull::

            sage: sum(map(len, H.bounded_regions_points()))
            11
        """
        points_by_sign = self.hull_points_by_sign()
        return tuple(sorted([tuple(sorted(points)) for points in points_by_sign.values()])) 
