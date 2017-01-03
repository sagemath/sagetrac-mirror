# -*- coding: utf-8 -*-
r"""
Hyperbolic Geodesics

This module implements the abstract base class for geodesics in
hyperbolic space of arbitrary dimension.  It also contains the
implementations for specific models of hyperbolic geometry.

AUTHORS:

- Greg Laun (2013): initial version

EXAMPLES:

We can construct geodesics in the upper half plane model, abbreviated
UHP for convenience::

    sage: g = HyperbolicPlane().UHP().get_geodesic(2, 3)
    sage: g
    Geodesic in UHP from 2 to 3

This geodesic can be plotted using :meth:`plot`, in this example we will show
the axis.

::

    sage: g.plot(axes=True)
    Graphics object consisting of 2 graphics primitives

.. PLOT::

    g = HyperbolicPlane().UHP().get_geodesic(2.0, 3.0)
    sphinx_plot(g.plot(axes=True))

::

    sage: g = HyperbolicPlane().UHP().get_geodesic(I, 3 + I)
    sage: g.length()
    arccosh(11/2)
    sage: g.plot(axes=True)
    Graphics object consisting of 2 graphics primitives

.. PLOT::

    sphinx_plot(HyperbolicPlane().UHP().get_geodesic(I, 3 + I).plot(axes=True))

Geodesics of both types in UHP are supported::

    sage: g = HyperbolicPlane().UHP().get_geodesic(I, 3*I)
    sage: g
    Geodesic in UHP from I to 3*I
    sage: g.plot()
    Graphics object consisting of 2 graphics primitives

.. PLOT::

    sphinx_plot(HyperbolicPlane().UHP().get_geodesic(I, 3*I).plot())

Geodesics are oriented, which means that two geodesics with the same
graph will only be equal if their starting and ending points are
the same::

    sage: g1 = HyperbolicPlane().UHP().get_geodesic(1,2)
    sage: g2 = HyperbolicPlane().UHP().get_geodesic(2,1)
    sage: g1 == g2
    False

.. TODO::

    Implement a parent for all geodesics of the hyperbolic plane?
    Or implement geodesics as a parent in the subobjects category?

"""


#***********************************************************************
#       Copyright (C) 2013 Greg Laun <glaun@math.umd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

from sage.structure.sage_object import SageObject
from sage.symbolic.all import I
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.infinity import infinity
from sage.rings.all import CC, RR
from sage.plot.arc import arc
from sage.plot.line import line
from sage.symbolic.constants import pi
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.functions.other import real, imag, sqrt
from sage.functions.trig import arccos
from sage.functions.log import exp
from sage.functions.hyperbolic import sinh, cosh, arcsinh
from sage.symbolic.ring import SR
from sage.geometry.hyperbolic_space.hyperbolic_constants import EPSILON
from sage.misc.superseded import deprecated_function_alias

from sage.misc.lazy_import import lazy_import
lazy_import('sage.geometry.hyperbolic_space.hyperbolic_isometry',
            'moebius_transform')


class HyperbolicGeodesic(SageObject):
    r"""
    Abstract base class for oriented geodesics that are not necessarily
    complete.

    INPUT:

    - ``start`` -- a HyperbolicPoint or coordinates of a point in
      hyperbolic space representing the start of the geodesic

    - ``end`` -- a HyperbolicPoint or coordinates of a point in
      hyperbolic space representing the end of the geodesic

    EXAMPLES:

    We can construct a hyperbolic geodesic in any model::

        sage: HyperbolicPlane().UHP().get_geodesic(1, 0)
        Geodesic in UHP from 1 to 0
        sage: HyperbolicPlane().PD().get_geodesic(1, 0)
        Geodesic in PD from 1 to 0
        sage: HyperbolicPlane().KM().get_geodesic((0,1/2), (1/2, 0))
        Geodesic in KM from (0, 1/2) to (1/2, 0)
        sage: HyperbolicPlane().HM().get_geodesic((0,0,1), (0,1, sqrt(2)))
        Geodesic in HM from (0, 0, 1) to (0, 1, sqrt(2))

    """

    #####################
    # "Private" Methods #
    #####################

    def __init__(self, model, start, end, **graphics_options):
        r"""
        See :class:`HyperbolicGeodesic` for full documentation.

        EXAMPLES ::

            sage: HyperbolicPlane().UHP().get_geodesic(I, 2 + I)
            Geodesic in UHP from I to I + 2

        """

        self._model = model
        self._start = start
        self._end = end
        self._graphics_options = graphics_options

    @lazy_attribute
    def _cached_geodesic(self):
        r"""
        The representation of the geodesic used for calculations.

        EXAMPLES::

            sage: A = HyperbolicPlane().PD().get_geodesic(0, 1/2)
            sage: A._cached_geodesic
            Geodesic in UHP from I to 3/5*I + 4/5

        """

        M = self._model.realization_of().a_realization()
        return self.to_model(M)

    @lazy_attribute
    def _complete(self):
        r"""
        Return whether the geodesic is complete.  This is used for
        geodesics in non-bounded models.  For these models,
        ``self.complete()`` simply sets ``_complete`` to ``True``.

        EXAMPLES::

            sage: HyperbolicPlane().UHP().get_geodesic(1, -12)._complete
            True
            sage: HyperbolicPlane().UHP().get_geodesic(I, 2 + I)._complete
            False
            sage: HM = HyperbolicPlane().HM()
            sage: g = HM.get_geodesic((0,0,1), (0,1, sqrt(2)))
            sage: g._complete
            False
            sage: g.complete()._complete
            True

        """

        if self._model.is_bounded():
            return (self._start.is_boundary() and self._end.is_boundary())
        return False  # All non-bounded geodesics start life incomplete.

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.get_geodesic(3 + 4*I, I)
            Geodesic in UHP from 4*I + 3 to I

            sage: PD = HyperbolicPlane().PD()
            sage: PD.get_geodesic(1/2 + I/2, 0)
            Geodesic in PD from 1/2*I + 1/2 to 0

            sage: KM = HyperbolicPlane().KM()
            sage: KM.get_geodesic((1/2, 1/2), (0, 0))
            Geodesic in KM from (1/2, 1/2) to (0, 0)

            sage: HM = HyperbolicPlane().HM()
            sage: HM.get_geodesic((0,0,1), (0, 1, sqrt(Integer(2))))
            Geodesic in HM from (0, 0, 1) to (0, 1, sqrt(2))

        """

        msg = "Geodesic in {0} from {1} to {2}"
        return msg.format(self._model.short_name(),
                          self._start.coordinates(),
                          self._end.coordinates())

    def __eq__(self, other):
        r"""
        Return ``True`` if ``self`` is equal to other as an oriented geodesic.

        EXAMPLES::

            sage: g1 = HyperbolicPlane().UHP().get_geodesic(I, 2*I)
            sage: g2 = HyperbolicPlane().UHP().get_geodesic(2*I, I)
            sage: g1 == g2
            False
            sage: g1 == g1
            True

        """

        if not isinstance(other, HyperbolicGeodesic):
            return False
        return (self._model is other._model and
                self._start == other._start and
                self._end == other._end)

    #######################
    # Setters and Getters #
    #######################

    def start(self):
        r"""
        Return the starting point of the geodesic.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 3*I)
            sage: g.start()
            Point in UHP I

        """
        return self._start

    def end(self):
        r"""
        Return the starting point of the geodesic.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 3*I)
            sage: g.end()
            Point in UHP 3*I

        """

        return self._end

    def endpoints(self):
        r"""
        Return a list containing the start and endpoints.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 3*I)
            sage: g.endpoints()
            [Point in UHP I, Point in UHP 3*I]

        """

        return [self._start, self._end]

    def model(self):
        r"""
        Return the model to which the :class:`HyperbolicGeodesic` belongs.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.get_geodesic(I, 2*I).model()
            Hyperbolic plane in the Upper Half Plane Model model

            sage: PD = HyperbolicPlane().PD()
            sage: PD.get_geodesic(0, I/2).model()
            Hyperbolic plane in the Poincare Disk Model model

            sage: KM = HyperbolicPlane().KM()
            sage: KM.get_geodesic((0, 0), (0, 1/2)).model()
            Hyperbolic plane in the Klein Disk Model model

            sage: HM = HyperbolicPlane().HM()
            sage: HM.get_geodesic((0, 0, 1), (0, 1, sqrt(2))).model()
            Hyperbolic plane in the Hyperboloid Model model

        """

        return self._model

    def to_model(self, model):
        r"""
        Convert the current object to image in another model.

        INPUT:

        - ``model`` -- the image model

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: PD = HyperbolicPlane().PD()
            sage: UHP.get_geodesic(I, 2*I).to_model(PD)
            Geodesic in PD from 0 to 1/3*I
            sage: UHP.get_geodesic(I, 2*I).to_model('PD')
            Geodesic in PD from 0 to 1/3*I

        """

        if isinstance(model, str):
            model = getattr(self._model.realization_of(), model)()
        if not model.is_bounded() and self.length() == infinity:
            raise NotImplementedError("cannot convert to an unbounded model")
        start = model(self._start)
        end = model(self._end)
        g = model.get_geodesic(start, end)
        return g

    def graphics_options(self):
        r"""
        Return the graphics options of ``self``.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 2*I, color="red")
            sage: g.graphics_options()
            {'color': 'red'}

        """

        return self._graphics_options

    def update_graphics(self, update=False, **options):
        r"""
        Update the graphics options of ``self``.

        INPUT:

        - ``update`` -- if ``True``, the original option are updated
          rather than overwritten

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(I, 2*I)
            sage: g.graphics_options()
            {}

            sage: g.update_graphics(color = "red"); g.graphics_options()
            {'color': 'red'}

            sage: g.update_graphics(color = "blue"); g.graphics_options()
            {'color': 'blue'}

            sage: g.update_graphics(True, size = 20); g.graphics_options()
            {'color': 'blue', 'size': 20}

        """

        if not update:
            self._graphics_options = {}
        self._graphics_options.update(**options)

    ###################
    # Boolean Methods #
    ###################

    def is_complete(self):
        r"""
        Return ``True`` if ``self`` is a complete geodesic (that is, both
        endpoints are on the ideal boundary) and ``False`` otherwise.

        If we represent complete geodesics using green color and incomplete
        using red colors we have the following graphic:

        .. PLOT::

            UHP = HyperbolicPlane().UHP()
            g = UHP.get_geodesic(1.5*I, 2.5*I)
            h = UHP.get_geodesic(0, I)
            l = UHP.get_geodesic(2, 4)
            m = UHP.get_geodesic(3, infinity)
            G = g.plot(color='red') +\
                text('is_complete()=False',
                     (0, 2),
                     horizontal_alignement='left')
            H = h.plot(color='red') +\
                text('is_complete()=False',
                     (0, 0.5),
                     horizontal_alignement='left')
            L = l.plot(color='green') +\
                text('is_complete()=True',
                     (5, 1.5))
            M = m.plot(color='green') + text('is complete()=True',
                                             (5, 4),
                                             horizontal_alignement='left')
            sphinx_plot(G+H+L+M)

        Notice, that there is no visual indication that the *vertical* geodesic
        is complete

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.get_geodesic(1.5*I, 2.5*I).is_complete()
            False
            sage: UHP.get_geodesic(0, I).is_complete()
            False
            sage: UHP.get_geodesic(3, infinity).is_complete()
            True
            sage: UHP.get_geodesic(2,5).is_complete()
            True

        """

        return self._complete

    def is_asymptotically_parallel(self, other):
        r"""
        Return ``True`` if ``self`` and ``other`` are asymptotically
        parallel and ``False`` otherwise.

        INPUT:

        - ``other`` -- a hyperbolic geodesic

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: h = HyperbolicPlane().UHP().get_geodesic(-2,4)
            sage: g.is_asymptotically_parallel(h)
            True

        .. PLOT::

             g = HyperbolicPlane().UHP().get_geodesic(-2.0,5.0)
             h = HyperbolicPlane().UHP().get_geodesic(-2.0,4.0)
             sphinx_plot(g.plot(color='green')+h.plot(color='green'))

        Ultraparallel geodesics are not asymptotically parallel::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: h = HyperbolicPlane().UHP().get_geodesic(-1,4)
            sage: g.is_asymptotically_parallel(h)
            False

        .. PLOT::

             g = HyperbolicPlane().UHP().get_geodesic(-2.0,5.0)
             h = HyperbolicPlane().UHP().get_geodesic(-1.0,4.0)
             sphinx_plot(g.plot(color='red')+h.plot(color='red'))


        No hyperbolic geodesic is asymptotically parallel to itself::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: g.is_asymptotically_parallel(g)
            False

        """

        p1, p2 = self.complete().endpoints()
        q1, q2 = other.complete().endpoints()
        return ((self != other) and ((p1 in [q1, q2]) or (p2 in [q1, q2])) and
                self.model() is other.model())

    def is_ultra_parallel(self, other):
        r"""
        Return ``True`` if ``self`` and ``other`` are ultra parallel
        and ``False`` otherwise.

        INPUT:

        - ``other`` -- a hyperbolic geodesic

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_geodesic \
            ....:   import *
            sage: g = HyperbolicPlane().UHP().get_geodesic(0,1)
            sage: h = HyperbolicPlane().UHP().get_geodesic(-3,-1)
            sage: g.is_ultra_parallel(h)
            True

        .. PLOT::

             g = HyperbolicPlane().UHP().get_geodesic(0.0,1.1)
             h = HyperbolicPlane().UHP().get_geodesic(-3.0,-1.0)
             sphinx_plot(g.plot(color='green')+h.plot(color='green'))

        ::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: h = HyperbolicPlane().UHP().get_geodesic(2,6)
            sage: g.is_ultra_parallel(h)
            False

        .. PLOT::

             g = HyperbolicPlane().UHP().get_geodesic(-2,5)
             h = HyperbolicPlane().UHP().get_geodesic(2,6)
             sphinx_plot(g.plot(color='red')+h.plot(color='red'))

        ::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: g.is_ultra_parallel(g)
            False

        """

        A = self.reflection_involution()
        B = other.reflection_involution()
        return (A * B).classification() == 'hyperbolic'

    def is_parallel(self, other):
        r"""
        Return ``True`` if the two given hyperbolic geodesics are either
        ultra parallel or asymptotically parallel and``False`` otherwise.

        INPUT:

        - ``other`` -- a hyperbolic geodesic in any model

        OUTPUT:

        ``True`` if the given geodesics are either ultra parallel or
        asymptotically parallel, ``False`` if not.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: h = HyperbolicPlane().UHP().get_geodesic(5,12)
            sage: g.is_parallel(h)
            True

        .. PLOT::

            g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            h = HyperbolicPlane().UHP().get_geodesic(5,12)
            sphinx_plot(g.plot(color='green')+h.plot(color='green'))

        ::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: h = HyperbolicPlane().UHP().get_geodesic(-2,4)
            sage: g.is_parallel(h)
            True

        .. PLOT::

            g = HyperbolicPlane().UHP().get_geodesic(-2.0,5.0)
            h = HyperbolicPlane().UHP().get_geodesic(-2.0,4.0)
            sphinx_plot(g.plot(color='green')+h.plot(color='green'))

        ::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,2)
            sage: h = HyperbolicPlane().UHP().get_geodesic(-1,4)
            sage: g.is_parallel(h)
            False

        .. PLOT::

            g = HyperbolicPlane().UHP().get_geodesic(-2,2)
            h = HyperbolicPlane().UHP().get_geodesic(-1,4)
            sphinx_plot(g.plot(color='red')+h.plot(color='red'))


        No hyperbolic geodesic is either ultra parallel or
        asymptotically parallel to itself::

            sage: g = HyperbolicPlane().UHP().get_geodesic(-2,5)
            sage: g.is_parallel(g)
            False

        """

        A = self.reflection_involution()
        B = other.reflection_involution()
        return (A * B).classification() in ['parabolic', 'hyperbolic']

    def ideal_endpoints(self):
        r"""
        Return the ideal endpoints in bounded models.  Raise a
        ``NotImplementedError`` in models that are not bounded.

        EXAMPLES::

            sage: H = HyperbolicPlane()
            sage: UHP = H.UHP()
            sage: UHP.get_geodesic(1 + I, 1 + 3*I).ideal_endpoints()
            [Boundary point in UHP 1, Boundary point in UHP +Infinity]

            sage: PD = H.PD()
            sage: PD.get_geodesic(0, I/2).ideal_endpoints()
            [Boundary point in PD -I, Boundary point in PD I]

            sage: KM = H.KM()
            sage: KM.get_geodesic((0,0), (0, 1/2)).ideal_endpoints()
            [Boundary point in KM (0, -1), Boundary point in KM (0, 1)]

            sage: HM = H.HM()
            sage: HM.get_geodesic((0,0,1), (1, 0, sqrt(2))).ideal_endpoints()
            Traceback (most recent call last):
            ...
            NotImplementedError: boundary points are not implemented in
             the HM model

        """

        if not self._model.is_bounded():
            errtxt = "boundary points are not implemented in the " + \
                     "{0} model".format(self._model.short_name())
            raise NotImplementedError(errtxt)
        if self.is_complete():
            return self.endpoints()
        return [self._model(k)
                for k in self._cached_geodesic.ideal_endpoints()]

    def complete(self):
        r"""
        Return the geodesic with ideal endpoints in bounded models.  Raise a
        ``NotImplementedError`` in models that are not bounded.
        In the following examples we represent complete geodesics by a dashed
        line.

        EXAMPLES::

            sage: H = HyperbolicPlane()
            sage: UHP = H.UHP()
            sage: UHP.get_geodesic(1 + I, 1 + 3*I).complete()
            Geodesic in UHP from 1 to +Infinity

        .. PLOT::

             g = HyperbolicPlane().UHP().get_geodesic(1 + I, 1 + 3*I)
             h = g.complete()
             sphinx_plot(g.plot()+h.plot(linestyle='dashed'))

        ::

            sage: PD = H.PD()
            sage: PD.get_geodesic(0, I/2).complete()
            Geodesic in PD from -I to I
            sage: PD.get_geodesic(0.25*(-1-I),0.25*(1-I)).complete()
            Geodesic in PD from -0.895806416477617 - 0.444444444444444*I to 0.895806416477617 - 0.444444444444444*I

        .. PLOT::

            PD = HyperbolicPlane().PD()
            g = PD.get_geodesic(0, I/2)
            h = g. complete()
            m = PD.get_geodesic(0.25*(-1-I),0.25*(1-I))
            l = m.complete()
            sphinx_plot(g.plot()+h.plot(linestyle='dashed') +
                        m.plot()+l.plot(linestyle='dashed'))

        ::

            sage: KM = H.KM()
            sage: KM.get_geodesic((0,0), (0, 1/2)).complete()
            Geodesic in KM from (0, -1) to (0, 1)

        .. PLOT::

            g = HyperbolicPlane().KM().get_geodesic((0.0,0.0), (0.0, 0.5))
            h = g.complete()
            sphinx_plot(g.plot()+h.plot(linestyle='dashed'))

        ::

            sage: HM = H.HM()
            sage: HM.get_geodesic((0,0,1), (1, 0, sqrt(2))).complete()
            Geodesic in HM from (0, 0, 1) to (1, 0, sqrt(2))

        .. PLOT::

            g = HyperbolicPlane().HM().get_geodesic((0,0,1), (1, 0, sqrt(2)))
            h = g.complete()
            sphinx_plot(g.plot(color='black')+h.plot(linestyle='dashed',color='black'))

        ::

            sage: g = HM.get_geodesic((0,0,1), (1, 0, sqrt(2))).complete()
            sage: g.is_complete()
            True

        TEST:

        Check that floating points remain floating points through this method::

            sage: H = HyperbolicPlane()
            sage: g = H.UHP().get_geodesic(CC(0,1), CC(2,2))
            sage: gc = g.complete()
            sage: parent(gc.start().coordinates())
            Real Field with 53 bits of precision

        """

        if self._model.is_bounded():
            return self._model.get_geodesic(*self.ideal_endpoints())

        from copy import copy
        g = copy(self)
        g._complete = True
        return g

    def reflection_involution(self):
        r"""
        Return the involution fixing ``self``.

        EXAMPLES::

            sage: H = HyperbolicPlane()
            sage: gU = H.UHP().get_geodesic(2,4)
            sage: RU = gU.reflection_involution(); RU
            Isometry in UHP
            [ 3 -8]
            [ 1 -3]

            sage: RU*gU == gU
            True

            sage: gP = H.PD().get_geodesic(0, I)
            sage: RP = gP.reflection_involution(); RP
            Isometry in PD
            [ 1  0]
            [ 0 -1]

            sage: RP*gP == gP
            True

            sage: gK = H.KM().get_geodesic((0,0), (0,1))
            sage: RK = gK.reflection_involution(); RK
            Isometry in KM
            [-1  0  0]
            [ 0  1  0]
            [ 0  0  1]

            sage: RK*gK == gK
            True

            sage: HM = H.HM()
            sage: g = HM.get_geodesic((0,0,1), (1,0, n(sqrt(2))))
            sage: A = g.reflection_involution()
            sage: B = diagonal_matrix([1, -1, 1])
            sage: bool((B - A.matrix()).norm() < 10**-9)
            True

        The above tests go through the Upper Half Plane.  It remains to
        test that the matrices in the models do what we intend. ::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry \
            ....:   import moebius_transform
            sage: R = H.PD().get_geodesic(-1,1).reflection_involution()
            sage: bool(moebius_transform(R.matrix(), 0) == 0)
            True

        """

        ri = self._cached_geodesic.reflection_involution()
        return ri.to_model(self._model)

    def common_perpendicula(self, other):
        r"""
        Return the unique hyperbolic geodesic perpendicular to two given
        geodesics, if such a geodesic exists.  If none exists, raise a
        ``ValueError``.

        INPUT:

        - ``other`` -- a hyperbolic geodesic in the same model as ``self``

        OUTPUT:

        - a hyperbolic geodesic

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2,3)
            sage: h = HyperbolicPlane().UHP().get_geodesic(4,5)
            sage: g.common_perpendicular(h)
            Geodesic in UHP from 1/2*sqrt(3) + 7/2 to -1/2*sqrt(3) + 7/2

        .. PLOT::

            g = HyperbolicPlane().UHP().get_geodesic(2.0, 3.0)
            h = HyperbolicPlane().UHP().get_geodesic(4.0, 5.0)
            l = g.common_perpendicular(h)
            P = g.plot(color='blue') +\
                h.plot(color='blue') +\
                l.plot(color='orange')
            sphinx_plot(P)

        It is an error to ask for the common perpendicular of two
        intersecting geodesics::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2,4)
            sage: h = HyperbolicPlane().UHP().get_geodesic(3, infinity)
            sage: g.common_perpendicular(h)
            Traceback (most recent call last):
            ...
            ValueError: geodesics intersect; no common perpendicular exists

        """

        if not self.is_parallel(other):
            raise ValueError('geodesics intersect; ' +
                             'no common perpendicular exists')
        cp = self._cached_geodesic.common_perpendicular(other)
        return cp.to_model(self._model)

    def intersection(self, other):
        r"""
        Return the point of intersection of two geodesics (if such a
        point exists).

        INPUT:

        - ``other`` -- a hyperbolic geodesic in the same model as ``self``

        OUTPUT:

        - a hyperbolic point or geodesic

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()

        """

        if self == other:
            return self
        elif self.is_parallel(other):
            raise ValueError("geodesics don't intersect")
        inters = self._cached_geodesic.intersection(other)
        if len(inters) == 2:
            return self
        elif len(inters) == 1:
            return [self._model(inters[0])]
        return []

    def perpendicular_bisector(self):
        r"""
        Return the perpendicular bisector of ``self`` if ``self`` has
        finite length.  Here distance is hyperbolic distance.

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: g = PD.get_geodesic(-0.3+0.4*I,+0.7-0.1*I)
            sage: h = g.perpendicular_bisector()
            sage: P = g.plot(color='blue')+h.plot(color='orange');P
            Graphics object consisting of 4 graphics primitives

        .. PLOT::

            g = HyperbolicPlane().PD().get_geodesic(-0.3+0.4*I,+0.7-0.1*I)
            h = g.perpendicular_bisector()
            sphinx_plot(g.plot(color='blue')+h.plot(color='orange'))

        Complete geodesics cannot be bisected::

            sage: g = HyperbolicPlane().PD().get_geodesic(0, 1)
            sage: g.perpendicular_bisector()
            Traceback (most recent call last):
            ...
            ValueError: the length must be finite

        TEST::

            sage: g = HyperbolicPlane().PD().random_geodesic()
            sage: h = g.perpendicular_bisector()
            sage: bool(h.intersection(g)[0].coordinates() - g.midpoint().coordinates() < 10**-9)
            True

        """

        P = self._cached_geodesic.perpendicular_bisector()
        return P.to_model(self._model)

    def midpoint(self):
        r"""
        Return the (hyperbolic) midpoint of a hyperbolic line segment.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().random_geodesic()
            sage: m = g.midpoint()
            sage: end1, end2 = g.endpoints()
            sage: bool(abs(m.dist(end1) - m.dist(end2)) < 10**-9)
            True

        Complete geodesics have no midpoint::

            sage: HyperbolicPlane().UHP().get_geodesic(0,2).midpoint()
            Traceback (most recent call last):
            ...
            ValueError: the length must be finite

        """

        UHP = self._model.realization_of().a_realization()
        P = self.to_model(UHP).midpoint()
        return self._model(P)

    def dist(self, other):
        r"""
        Return the hyperbolic distance from a given hyperbolic geodesic
        to another geodesic or point.

        INPUT:

        - ``other`` -- a hyperbolic geodesic or hyperbolic point in
          the same model

        OUTPUT:

        - the hyperbolic distance

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2, 4.0)
            sage: h = HyperbolicPlane().UHP().get_geodesic(5, 7.0)
            sage: bool(abs(g.dist(h).n() - 1.92484730023841) < 10**-9)
            True

        If the second object is a geodesic ultraparallel to the first,
        or if it is a point on the boundary that is not one of the
        first object's endpoints, then return +infinity

        ::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2, 2+I)
            sage: p = HyperbolicPlane().UHP().get_point(5)
            sage: g.dist(p)
            +Infinity

        TEST:

        Check that floating points remain floating points in :meth:`dist` ::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.get_geodesic(CC(0,1), CC(2,2))
            sage: UHP.dist(g.start(), g.end())
            1.45057451382258
            sage: parent(_)
            Real Field with 53 bits of precision

        """

        return self._model.dist(self, other)

    def angle(self, other):
        r"""
        Return the angle  between any two given geodesics if they
        intersect.

        INPUT:

        - ``other`` -- a hyperbolic geodesic in the same model as ``self``

        OUTPUT:

        - the angle in radians between the two given geodesics

        EXAMPLES::

            sage: PD = HyperbolicPlane().PD()
            sage: g = PD.get_geodesic(3/5*I + 4/5, 15/17*I + 8/17)
            sage: h = PD.get_geodesic(4/5*I + 3/5, I)
            sage: g.angle(h)
            1/2*pi

        .. PLOT::

            PD = HyperbolicPlane().PD()
            g = PD.get_geodesic(3.0/5.0*I + 4.0/5.0, 15.0/17.0*I + 8.0/17.0)
            h = PD.get_geodesic(4.0/5.0*I + 3.0/5.0, I)
            sphinx_plot(g.plot()+h.plot(color='orange'))

        """

        if self.is_parallel(other):
            raise ValueError("geodesics do not intersect")
        return self._cached_geodesic.angle(other)

    def length(self):
        r"""
        Return the Hyperbolic length of the hyperbolic line segment.

        EXAMPLES::

            sage: g = HyperbolicPlane().UHP().get_geodesic(2 + I, 3 + I/2)
            sage: g.length()
            arccosh(9/4)

        """

        return self._model._dist_points(self._start.coordinates(),
                                        self._end.coordinates())

#***********************************************************************
#                       UHP geodesics
#***********************************************************************


class HyperbolicGeodesicUHP(HyperbolicGeodesic):
    r"""
    Create a geodesic in the upper half plane model.

    The geodesics in this model are represented by circular arcs perpendicular
    to the real axis (half-circles whose origin is on the real axis) and
    straight vertical lines ending on the real axis.

    INPUT:

    - ``start`` -- a :class:`HyperbolicPoint` in hyperbolic space
      representing the start of the geodesic

    - ``end`` -- a :class:`HyperbolicPoint` in hyperbolic space
      representing the end of the geodesic

    EXAMPLES::

        sage: UHP = HyperbolicPlane().UHP()
        sage: g = UHP.get_geodesic(UHP.get_point(I), UHP.get_point(2 + I))
        sage: g = UHP.get_geodesic(I, 2 + I)
        sage: h = UHP.get_geodesic(-1, -1+2*I)

    .. PLOT::

        UHP = HyperbolicPlane().UHP()
        g = UHP.get_geodesic(I, 2 + I)
        h = UHP.get_geodesic(-1, -1+2*I)
        sphinx_plot(g.plot()+h.plot())

    """

    def reflection_involution(self):
        r"""
        Return the isometry of the involution fixing the geodesic ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g1 = UHP.get_geodesic(0, 1)
            sage: g1.reflection_involution()
            Isometry in UHP
            [ 1  0]
            [ 2 -1]
            sage: UHP.get_geodesic(I, 2*I).reflection_involution()
            Isometry in UHP
            [ 1  0]
            [ 0 -1]

        """

        x, y = [real(k.coordinates()) for k in self.ideal_endpoints()]
        if x == infinity:
            M = matrix([[1, -2*y], [0, -1]])
        elif y == infinity:
            M = matrix([[1, -2*x], [0, -1]])
        else:
            M = matrix([[(x+y)/(y-x), -2*x*y/(y-x)], [2/(y-x), -(x+y)/(y-x)]])
        return self._model.get_isometry(M)

    def plot(self, boundary=True, **options):
        r"""
        Plot ``self``.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.get_geodesic(0, 1).plot()
            Graphics object consisting of 2 graphics primitives

        .. PLOT::

            UHP = HyperbolicPlane().UHP()
            g = UHP.get_geodesic(0.0, 1.0).plot()
            sphinx_plot(g)

        ::

            sage: UHP.get_geodesic(I, 3+4*I).plot(linestyle="dashed", color="brown")
            Graphics object consisting of 2 graphics primitives

        .. PLOT::

            UHP = HyperbolicPlane().UHP()
            g = UHP.get_geodesic(I, 3+4*I).plot(linestyle="dashed", color="brown")
            sphinx_plot(g)

        ::

            sage: UHP.get_geodesic(1, infinity).plot(color='orange')
            Graphics object consisting of 2 graphics primitives

        .. PLOT::

            UHP = HyperbolicPlane().UHP()
            g = UHP.get_geodesic(1, infinity).plot(color='orange')
            sphinx_plot(g)

        """

        opts = {'axes': False, 'aspect_ratio': 1}
        opts.update(self.graphics_options())
        opts.update(options)
        end_1, end_2 = [CC(k.coordinates()) for k in self.endpoints()]
        bd_1, bd_2 = [CC(k.coordinates()) for k in self.ideal_endpoints()]
        if (abs(real(end_1) - real(end_2)) < EPSILON) \
                or CC(infinity) in [end_1, end_2]:  # on same vertical line
            # If one of the endpoints is infinity, we replace it with a
            # large finite  point
            if end_1 == CC(infinity):
                end_1 = (real(end_2), (imag(end_2) + 10))
                end_2 = (real(end_2), imag(end_2))
            elif end_2 == CC(infinity):
                end_2 = (real(end_1), (imag(end_1) + 10))
                end_1 = (real(end_1), imag(end_1))
            pic = line((end_1, end_2), **opts)
            if boundary:
                cent = min(bd_1, bd_2)
                bd_dict = {'bd_min': cent - 3, 'bd_max': cent + 3}
                bd_pic = self._model.get_background_graphic(**bd_dict)
                pic = bd_pic + pic
                return pic
        else:
            center = (bd_1 + bd_2) / 2  # Circle center
            radius = abs(bd_1 - bd_2) / 2
            theta1 = CC(end_1 - center).arg()
            theta2 = CC(end_2 - center).arg()
            if abs(theta1 - theta2) < EPSILON:
                theta2 += pi
            pic = arc((real(center), imag(center)), radius,
                      sector=(theta1, theta2), **opts)
            if boundary:
                # We want to draw a segment of the real line.  The
                # computations below compute the projection of the
                # geodesic to the real line, and then draw a little
                # to the left and right of the projection.
                shadow_1, shadow_2 = [real(k) for k in [end_1, end_2]]
                midpoint = (shadow_1 + shadow_2)/2
                length = abs(shadow_1 - shadow_2)
                bd_dict = {'bd_min': midpoint - length, 'bd_max': midpoint +
                           length}
                bd_pic = self._model.get_background_graphic(**bd_dict)
                pic = bd_pic + pic
            return pic

    show = deprecated_function_alias(20530, plot)

    def ideal_endpoints(self):
        r"""
        Determine the ideal (boundary) endpoints of the complete
        hyperbolic geodesic corresponding to ``self``.

        OUTPUT:

        - a list of 2 boundary points

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: UHP.get_geodesic(I, 2*I).ideal_endpoints()
            [Boundary point in UHP 0,
             Boundary point in UHP +Infinity]
            sage: UHP.get_geodesic(1 + I, 2 + 4*I).ideal_endpoints()
            [Boundary point in UHP -sqrt(65) + 9,
             Boundary point in UHP sqrt(65) + 9]

        """

        start = self._start.coordinates()
        end = self._end.coordinates()
        [x1, x2] = [real(k) for k in [start, end]]
        [y1, y2] = [imag(k) for k in [start, end]]
        M = self._model
        # infinity is the first endpoint, so the other ideal endpoint
        # is just the real part of the second coordinate
        if start == infinity:
            return [M.get_point(start), M.get_point(x2)]
        # Same idea as above
        if end == infinity:
            return [M.get_point(x1), M.get_point(end)]
        # We could also have a vertical line with two interior points
        if x1 == x2:
            return [M.get_point(x1), M.get_point(infinity)]
        # Otherwise, we have a semicircular arc in the UHP
        c = ((x1+x2)*(x2-x1) + (y1+y2)*(y2-y1)) / (2*(x2-x1))
        r = sqrt((c - x1)**2 + y1**2)
        return [M.get_point(c - r), M.get_point(c + r)]

    def common_perpendicular(self, other):
        r"""
        Return the unique hyperbolic geodesic perpendicular to ``self``
        and ``other``, if such a geodesic exists; otherwise raise a
        ``ValueError``.

        INPUT:

        - ``other`` -- a hyperbolic geodesic in current model

        OUTPUT:

        - a hyperbolic geodesic

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.get_geodesic(2, 3)
            sage: h = UHP.get_geodesic(4, 5)
            sage: g.common_perpendicular(h)
            Geodesic in UHP from 1/2*sqrt(3) + 7/2 to -1/2*sqrt(3) + 7/2

        .. PLOT::

            UHP = HyperbolicPlane().UHP()
            g = UHP.get_geodesic(2.0, 3.0)
            h = UHP.get_geodesic(4.0, 5.0)
            p = g.common_perpendicular(h)
            sphinx_plot(g.plot(color='blue')+h.plot(color='blue')+p.plot(color='orange'))

        It is an error to ask for the common perpendicular of two
        intersecting geodesics::

            sage: g = UHP.get_geodesic(2, 4)
            sage: h = UHP.get_geodesic(3, infinity)
            sage: g.common_perpendicular(h)
            Traceback (most recent call last):
            ...
            ValueError: geodesics intersect; no common perpendicular exists

        """

        # Make sure both are in the same model
        if other._model is not self._model:
            other = other.to_model(self._model)

        A = self.reflection_involution()
        B = other.reflection_involution()
        C = A * B
        if C.classification() != 'hyperbolic':
            raise ValueError("geodesics intersect; " +
                             "no common perpendicular exists")
        return C.fixed_point_set()

    def intersection(self, other):
        r"""
        Return the point of intersection of ``self`` and ``other``
        (if such a point exists).

        INPUT:

        - ``other`` -- a hyperbolic geodesic in the current model

        OUTPUT:

        - a list of hyperbolic points or a hyperbolic geodesic

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.get_geodesic(3, 5)
            sage: h = UHP.get_geodesic(4, 7)
            sage: g.intersection(h)
            [Point in UHP 2/3*sqrt(-2) + 13/3]

        If the given geodesics do not intersect, the function returns an
        empty list::

            sage: g = UHP.get_geodesic(4, 5)
            sage: h = UHP.get_geodesic(5, 7)
            sage: g.intersection(h)
            []

        If the given geodesics are identical, return that
        geodesic::

            sage: g = UHP.get_geodesic(4 + I, 18*I)
            sage: h = UHP.get_geodesic(4 + I, 18*I)
            sage: g.intersection(h)
            [Boundary point in UHP -1/8*sqrt(114985) - 307/8,
             Boundary point in UHP 1/8*sqrt(114985) - 307/8]

        """

        start_1, end_1 = sorted(self.ideal_endpoints(), key=str)
        start_2, end_2 = sorted(other.ideal_endpoints(), key=str)
        if start_1 == start_2 and end_1 == end_2:  # Unoriented geods are same
            return [start_1, end_1]
        if start_1 == start_2:
            return start_1
        elif end_1 == end_2:
            return end_1
        A = self.reflection_involution()
        B = other.reflection_involution()
        C = A * B
        if C.classification() in ['hyperbolic', 'parabolic']:
            return []
        return C.fixed_point_set()

    def perpendicular_bisector(self):  # UHP
        r"""
        Return the perpendicular bisector of the hyperbolic geodesic ``self``
        if that geodesic has finite length.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.random_geodesic()
            sage: h = g.perpendicular_bisector()
            sage: c = lambda x: x.coordinates()
            sage: bool(c(g.intersection(h)[0]) - c(g.midpoint()) < 10**-9)
            True

        ::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.get_geodesic(1+I,2+0.5*I)
            sage: h = g.perpendicular_bisector()
            sage: show(g.plot(color='blue')+h.plot(color='orange'))

        .. PLOT::

            UHP = HyperbolicPlane().UHP()
            g = UHP.get_geodesic(1+I,2+0.5*I)
            h = g.perpendicular_bisector()
            sphinx_plot(g.plot(color='blue')+h.plot(color='orange'))

        Infinite geodesics cannot be bisected::

            sage: UHP.get_geodesic(0, 1).perpendicular_bisector()
            Traceback (most recent call last):
            ...
            ValueError: the length must be finite

        """

        if self.length() == infinity:
            raise ValueError("the length must be finite")
        start = self._start.coordinates()
        d = self._model._dist_points(start, self._end.coordinates()) / 2
        S = self.complete()._to_std_geod(start)
        T1 = matrix([[exp(d/2), 0], [0, exp(-d/2)]])
        s2 = sqrt(2) * 0.5
        T2 = matrix([[s2, -s2], [s2, s2]])
        isom_mtrx = S.inverse() * (T1 * T2) * S
        # We need to clean this matrix up.
        if (isom_mtrx - isom_mtrx.conjugate()).norm() < 5 * EPSILON:
            # Imaginary part is small.
            isom_mtrx = (isom_mtrx + isom_mtrx.conjugate()) / 2
            # Set it to its real part.
        H = self._model.get_isometry(isom_mtrx)
        return self._model.get_geodesic(H(self._start), H(self._end))

    def midpoint(self):  # UHP
        r"""
        Return the (hyperbolic) midpoint of ``self`` if it exists.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.random_geodesic()
            sage: m = g.midpoint()
            sage: d1 = UHP.dist(m, g.start())
            sage: d2 = UHP.dist(m, g.end())
            sage: bool(abs(d1 - d2) < 10**-9)
            True

        Infinite geodesics have no midpoint::

            sage: UHP.get_geodesic(0, 2).midpoint()
            Traceback (most recent call last):
            ...
            ValueError: the length must be finite

        TESTS:

        This checks :trac:`20330` so that geodesics defined by symbolic
        expressions do not generate runtime errors. ::

            sage: g=HyperbolicPlane().UHP().get_geodesic(-1+I,1+I)
            sage: g.midpoint()
            Point in UHP 1/2*(sqrt(2)*e^(1/2*arccosh(3)) - sqrt(2) + (I - 1)*e^(1/2*arccosh(3)) + I - 1)/((1/4*I - 1/4)*sqrt(2)*e^(1/2*arccosh(3)) - (1/4*I - 1/4)*sqrt(2) + 1/2*e^(1/2*arccosh(3)) + 1/2)

        Check that floating points remain floating points
        in :meth:`midpoint` ::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.get_geodesic(CC(0,1), CC(2,2))
            sage: g.midpoint()
            Point in UHP 0.666666666666667 + 1.69967317119760*I
            sage: parent(g.midpoint().coordinates())
            Complex Field with 53 bits of precision

        """

        from sage.matrix.matrix_symbolic_dense import Matrix_symbolic_dense
        if self.length() == infinity:
            raise ValueError("the length must be finite")

        start = self._start.coordinates()
        end = self._end.coordinates()
        d = self._model._dist_points(start, end) / 2
        S = self.complete()._to_std_geod(start)

        # If the matrix is symbolic then needs to be simplified in order to
        #    make the calculations easier for the symbolic calculus module.
        if isinstance(S, Matrix_symbolic_dense):
            S = S.simplify_full().simplify_full()
        S_1 = S.inverse()
        T = matrix([[exp(d), 0], [0, 1]])
        M = S_1 * T * S
        if ((real(start - end) < EPSILON) or
            (abs(real(start - end)) < EPSILON and
                imag(start - end) < EPSILON)):
            end_p = start
        else:
            end_p = end
        P_3 = moebius_transform(M, end_p)
        return self._model.get_point(P_3)

    def angle(self, other):  # UHP
        r"""
        Return the angle between any two given completed geodesics if
        they intersect.

        INPUT:

        - ``other`` -- a hyperbolic geodesic in the UHP model

        OUTPUT:

        - the angle in radians between the two given geodesics

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.get_geodesic(2, 4)
            sage: h = UHP.get_geodesic(3, 3 + I)
            sage: g.angle(h)
            1/2*pi
            sage: numerical_approx(g.angle(h))
            1.57079632679490

        .. PLOT::

            UHP = HyperbolicPlane().UHP()
            g = UHP.get_geodesic(2, 4)
            h = UHP.get_geodesic(3, 3 + I)
            sphinx_plot(g.plot()+h.plot())

        If the geodesics are identical, return angle 0::

            sage: g.angle(g)
            0

        It is an error to ask for the angle of two geodesics that do not
        intersect::

            sage: g = UHP.get_geodesic(2, 4)
            sage: h = UHP.get_geodesic(5, 7)
            sage: g.angle(h)
            Traceback (most recent call last):
            ...
            ValueError: geodesics do not intersect

        """

        if self.is_parallel(other):
            raise ValueError("geodesics do not intersect")
        # Make sure the segments are complete or intersect
        if (not(self.is_complete() and other.is_complete()) and
                not self.intersection(other)):
            print("Warning: Geodesic segments do not intersect. "
                  "The angle between them is not defined.\n"
                  "Returning the angle between their completions.")

        # Make sure both are in the same model
        if other._model is not self._model:
            other = other.to_model(self._model)

        (p1, p2) = sorted([k.coordinates()
                           for k in self.ideal_endpoints()], key=str)
        (q1, q2) = sorted([k.coordinates()
                           for k in other.ideal_endpoints()], key=str)
        # if the geodesics are equal, the angle between them is 0
        if (abs(p1 - q1) < EPSILON and
                abs(p2 - q2) < EPSILON):
            return 0
        elif p2 != infinity:  # geodesic not a straight line
            # So we send it to the geodesic with endpoints [0, oo]
            T = HyperbolicGeodesicUHP._crossratio_matrix(p1, (p1 + p2) / 2, p2)
        else:
            # geodesic is a straight line, so we send it to the geodesic
            # with endpoints [0,oo]
            T = HyperbolicGeodesicUHP._crossratio_matrix(p1, p1 + 1, p2)
        # b1 and b2 are the endpoints of the image of other
        b1, b2 = [moebius_transform(T, k) for k in [q1, q2]]
        # If other is now a straight line...
        if (b1 == infinity or b2 == infinity):
            # then since they intersect, they are equal
            return 0
        return real(arccos((b1 + b2) / abs(b2 - b1)))

    ##################
    # Helper methods #
    ##################

    @staticmethod
    def _get_B(a):
        r"""
        Helper function to get an appropiate matrix transforming
        (0,1,inf)->(0,I,inf) based on the type of a

        INPUT:

        - ``a`` -- an element to identify the class of the resulting matrix.

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: g = UHP.random_geodesic()
            sage: B = g._get_B(CDF.an_element()); B
            [   1.0    0.0]
            [   0.0 -1.0*I]
            sage: type(B)
            <type 'sage.matrix.matrix_complex_double_dense.Matrix_complex_double_dense'>

      ::

            sage: B = g._get_B(SR(1));  B
            [ 1  0]
            [ 0 -I]
            sage: type(B)
            <type 'sage.matrix.matrix_symbolic_dense.Matrix_symbolic_dense'>

      ::

            sage: B = g._get_B(complex(1));  B
            [   1.0    0.0]
            [   0.0 -1.0*I]
            sage: type(B)
            <type 'sage.matrix.matrix_complex_double_dense.Matrix_complex_double_dense'>

      ::

            sage: B = g._get_B(QQbar(1+I));  B
            [ 1  0]
            [ 0 -I]
            sage: type(B[1,1])
            <class 'sage.rings.qqbar.AlgebraicNumber'>
            sage: type(B)
            <type 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>

        """
        from sage.structure.element import Element
        from sage.symbolic.expression import Expression
        from sage.rings.complex_double import CDF

        if isinstance(a, (int, float, complex)):  # Python number
            a = CDF(a)

        if isinstance(a, Expression):           # symbolic
            P = SR
            zero = SR.zero()
            one = SR.one()
            I = SR("I")
        elif isinstance(a, Element):            # Sage number
            P = a.parent()
            zero = P.zero()
            one = P.one()
            I = P.gen()
            if I.is_one() or (I*I).is_one() or not (-I*I).is_one():
                raise ValueError("invalid number")
        else:
            raise ValueError("not a complex number")

        return matrix(P, 2, [one, zero, zero, -I])

    def _to_std_geod(self, p):
        r"""
        Given the coordinates of a geodesic in hyperbolic space, return the
        hyperbolic isometry that sends that geodesic to the geodesic
        through 0 and infinity that also sends the point ``p`` to `i`.

        INPUT:

        - ``p`` -- the coordinates of the

        EXAMPLES::

            sage: UHP = HyperbolicPlane().UHP()
            sage: (p1, p2) = [UHP.random_point() for k in range(2)]
            sage: g = UHP.get_geodesic(p1, p2)
            sage: m = g.midpoint()
            sage: A = g._to_std_geod(m.coordinates()) # Send midpoint to I.
            sage: A = UHP.get_isometry(A)
            sage: [s, e]= g.complete().endpoints()
            sage: bool(abs(A(s).coordinates()) < 10**-9)
            True
            sage: bool(abs(A(m).coordinates() - I) < 10**-9)
            True
            sage: bool(abs(A(e).coordinates()) > 10**9)
            True

        TEST:

        Check that floating points remain floating points through this method::

            sage: H = HyperbolicPlane()
            sage: g = H.UHP().get_geodesic(CC(0,1), CC(2,2))
            sage: gc = g.complete()
            sage: parent(gc._to_std_geod(g.start().coordinates()))
            Full MatrixSpace of 2 by 2 dense matrices over Complex Field
            with 53 bits of precision

        """

        [s, e] = [k.coordinates() for k in self.complete().endpoints()]
        B = HyperbolicGeodesicUHP._get_B(p)
        # outmat below will be returned after we normalize the determinant.
        outmat = B * HyperbolicGeodesicUHP._crossratio_matrix(s, p, e)
        outmat = outmat / outmat.det().sqrt()
        if (outmat - outmat.conjugate()).norm(1) < 10**-9:
            # Small imaginary part.
            outmat = (outmat + outmat.conjugate()) / 2
            # Set it equal to its real part.
        return outmat

    @staticmethod
    def _crossratio_matrix(p0, p1, p2):  # UHP
        r"""
        Given three points (the list `p`) in `\mathbb{CP}^{1}` in affine
        coordinates, return the linear fractional transformation taking
        the elements of `p` to `0`, `1`, and `\infty`.

        INPUT:

        - a list of three distinct elements
          of `\mathbb{CP}^1` in affine coordinates; that is, each element
          must be a complex number, `\infty`, or symbolic.

        OUTPUT:

        - an element of `\GL(2,\CC)`

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_geodesic \
            ....:   import HyperbolicGeodesicUHP
            sage: from sage.geometry.hyperbolic_space.hyperbolic_isometry \
            ....:   import moebius_transform
            sage: UHP = HyperbolicPlane().UHP()
            sage: (p1, p2, p3) = [UHP.random_point().coordinates()
            ....:                   for k in range(3)]
            sage: A = HyperbolicGeodesicUHP._crossratio_matrix(p1, p2, p3)
            sage: bool(abs(moebius_transform(A, p1)) < 10**-9)
            True
            sage: bool(abs(moebius_transform(A, p2) - 1) < 10**-9)
            True
            sage: bool(moebius_transform(A, p3) == infinity)
            True
            sage: (x,y,z) = var('x,y,z');
            sage: HyperbolicGeodesicUHP._crossratio_matrix(x,y,z)
            [     y - z -x*(y - z)]
            [    -x + y  (x - y)*z]

        """

        if p0 == infinity:
            return matrix([[0, -(p1 - p2)], [-1, p2]])
        elif p1 == infinity:
            return matrix([[1, -p0], [1, -p2]])
        elif p2 == infinity:
            return matrix([[1, -p0], [0, p1 - p0]])
        return matrix([[p1 - p2, (p1 - p2)*(-p0)],
                       [p1 - p0, (p1 - p0)*(-p2)]])

#***********************************************************************
#                       Other geodesics
#***********************************************************************


class HyperbolicGeodesicPD(HyperbolicGeodesic):
    r"""
    A geodesic in the Poincaré disk model.

    Geodesics in this model are represented by segments of circles contained
    within the unit disk that are orthogonal to the boundary of the disk,
    plus all diameters of the disk.

    INPUT:

    - ``start`` -- a :class:`HyperbolicPoint` in hyperbolic space
      representing the start of the geodesic

    - ``end`` -- a :class:`HyperbolicPoint` in hyperbolic space
      representing the end of the geodesic

    EXAMPLES::

        sage: PD = HyperbolicPlane().PD()
        sage: g = PD.get_geodesic(PD.get_point(I), PD.get_point(-I/2))
        sage: g = PD.get_geodesic(I,-I/2)
        sage: h = PD.get_geodesic(-1/2+I/2,1/2+I/2)

    .. PLOT::

         PD = HyperbolicPlane().PD()
         g = PD.get_geodesic(I,-I/2)
         h = PD.get_geodesic(-0.5+I*0.5,0.5+I*0.5)
         sphinx_plot(g.plot()+h.plot(color='green'))

    """

    def plot(self, boundary=True, **options):

        r"""
        Plot ``self``.

        EXAMPLES:

        First some lines::

            sage: PD = HyperbolicPlane().PD()
            sage: PD.get_geodesic(0, 1).plot()
            Graphics object consisting of 2 graphics primitives

        .. PLOT::

            sphinx_plot(HyperbolicPlane().PD().get_geodesic(0, 1).plot())

        ::

            sage: PD.get_geodesic(0, 0.3+0.8*I).plot()
            Graphics object consisting of 2 graphics primitives

        .. PLOT::

            PD = HyperbolicPlane().PD()
            sphinx_plot(PD.get_geodesic(0, 0.3+0.8*I).plot())

        Then some generic geodesics::

            sage: PD.get_geodesic(-0.5, 0.3+0.4*I).plot()
            Graphics object consisting of 2 graphics primitives
            sage: g = PD.get_geodesic(-1, exp(3*I*pi/7))
            sage: G = g.plot(linestyle="dashed",color="red"); G
            Graphics object consisting of 2 graphics primitives
            sage: h = PD.get_geodesic(exp(2*I*pi/11), exp(1*I*pi/11))
            sage: H = h.plot(thickness=6, color="orange"); H
            Graphics object consisting of 2 graphics primitives
            sage: show(G+H)

        .. PLOT::

            PD = HyperbolicPlane().PD()
            PD.get_geodesic(-0.5, 0.3+0.4*I).plot()
            g = PD.get_geodesic(-1, exp(3*I*pi/7))
            G = g.plot(linestyle="dashed",color="red")
            h = PD.get_geodesic(exp(2*I*pi/11), exp(1*I*pi/11))
            H = h.plot(thickness=6, color="orange")
            sphinx_plot(G+H)

        """

        opts = {'axes': False, 'aspect_ratio': 1}
        opts.update(self.graphics_options())
        opts.update(options)
        end_1, end_2 = [CC(k.coordinates()) for k in self.endpoints()]
        bd_1, bd_2 = [CC(k.coordinates()) for k in self.ideal_endpoints()]
        # Check to see if it's a line
        if abs(bd_1 + bd_2) < EPSILON:
            pic = line([end_1, end_2], **opts)
        else:
            # If we are here, we know it's not a line
            # So we compute the center and radius of the circle
            invdet = RR.one() / (real(bd_1)*imag(bd_2) - real(bd_2)*imag(bd_1))
            centerx = (imag(bd_2) - imag(bd_1)) * invdet
            centery = (real(bd_1) - real(bd_2)) * invdet
            center = centerx + I * centery
            radius = RR(abs(bd_1 - center))
            # Now we calculate the angles for the arc
            theta1 = CC(end_1 - center).arg()
            theta2 = CC(end_2 - center).arg()
            theta1, theta2 = sorted([theta1, theta2])
            # Make sure the sector is inside the disk
            if theta2 - theta1 > pi:
                theta1 += 2 * pi
            pic = arc((centerx, centery), radius,
                      sector=(theta1, theta2), **opts)
        if boundary:
            pic += self._model.get_background_graphic()
        return pic

    show = deprecated_function_alias(20530, plot)


class HyperbolicGeodesicKM(HyperbolicGeodesic):
    r"""
    A geodesic in the Klein disk model.

    Geodesics are represented by the chords, straight line segments with ideal
    endpoints on the boundary circle.

    INPUT:

    - ``start`` -- a :class:`HyperbolicPoint` in hyperbolic space
      representing the start of the geodesic

    - ``end`` -- a :class:`HyperbolicPoint` in hyperbolic space
      representing the end of the geodesic

    EXAMPLES::

        sage: KM = HyperbolicPlane().KM()
        sage: g = KM.get_geodesic(KM.get_point((0.1,0.9)), KM.get_point((-0.1,-0.9)))
        sage: g = KM.get_geodesic((0.1,0.9),(-0.1,-0.9))
        sage: h = KM.get_geodesic((-0.707106781,-0.707106781),(0.707106781,-0.707106781))
        sage: P = g.plot(color='orange')+h.plot(); P
        Graphics object consisting of 4 graphics primitives


    .. PLOT::

        KM = HyperbolicPlane().KM()
        g = KM.get_geodesic((0.1,0.9),
                            (-0.1,-0.9))
        h = KM.get_geodesic((-0.707106781,-0.707106781),
                            (0.707106781,-0.707106781))
        sphinx_plot(g.plot(color='orange')+h.plot())

    """

    def plot(self, boundary=True, **options):
        r"""
        Plot ``self``.

        EXAMPLES::

            sage: HyperbolicPlane().KM().get_geodesic((0,0), (1,0)).plot()
            Graphics object consisting of 2 graphics primitives

        .. PLOT::

            KM = HyperbolicPlane().KM()
            sphinx_plot(KM.get_geodesic((0,0), (1,0)).plot())

        """
        opts = {'axes': False, 'aspect_ratio': 1}
        opts.update(self.graphics_options())
        opts.update(options)
        pic = line([k.coordinates() for k in self.endpoints()], **opts)
        if boundary:
            pic += self._model.get_background_graphic()
        return pic

    show = deprecated_function_alias(20530, plot)


class HyperbolicGeodesicHM(HyperbolicGeodesic):
    r"""
    A geodesic in the hyperboloid model.

    Valid points in the hyperboloid model satisfy :math:`x^2 + y^2 - z^2 = -1`

    INPUT:

    - ``start`` -- a :class:`HyperbolicPoint` in hyperbolic space
      representing the start of the geodesic

    - ``end`` -- a :class:`HyperbolicPoint` in hyperbolic space
      representing the end of the geodesic

    EXAMPLES::

        sage: from sage.geometry.hyperbolic_space.hyperbolic_geodesic import *
        sage: HM = HyperbolicPlane().HM()
        sage: p1 = HM.get_point((4, -4, sqrt(33)))
        sage: p2 = HM.get_point((-3,-3,sqrt(19)))
        sage: g = HM.get_geodesic(p1, p2)
        sage: g = HM.get_geodesic((4, -4, sqrt(33)), (-3, -3, sqrt(19)))

    .. PLOT::

        HM = HyperbolicPlane().HM()
        p1 = HM.get_point((4, -4, sqrt(33)))
        p2 = HM.get_point((-3,-3,sqrt(19)))
        g = HM.get_geodesic(p1, p2)
        sphinx_plot(g.plot(color='blue'))

    """

    def plot(self, show_hyperboloid=True, **graphics_options):
        r"""
        Plot ``self``.

        EXAMPLES::

            sage: from sage.geometry.hyperbolic_space.hyperbolic_geodesic \
            ....:    import *
            sage: g = HyperbolicPlane().HM().random_geodesic()
            sage: g.plot()
            Graphics3d Object

        .. PLOT::

            sphinx_plot(HyperbolicPlane().HM().random_geodesic().plot())

        """

        x = SR.var('x')
        opts = self.graphics_options()
        opts.update(graphics_options)
        v1, u2 = [vector(k.coordinates()) for k in self.endpoints()]
        # Lorentzian Gram Shmidt.  The original vectors will be
        # u1, u2 and the orthogonal ones will be v1, v2.  Except
        # v1 = u1, and I don't want to declare another variable,
        # hence the odd naming convention above.
        # We need the Lorentz dot product of v1 and u2.
        v1_ldot_u2 = u2[0]*v1[0] + u2[1]*v1[1] - u2[2]*v1[2]
        v2 = u2 + v1_ldot_u2 * v1
        v2_norm = sqrt(v2[0]**2 + v2[1]**2 - v2[2]**2)
        v2 = v2 / v2_norm
        v2_ldot_u2 = u2[0]*v2[0] + u2[1]*v2[1] - u2[2]*v2[2]
        # Now v1 and v2 are Lorentz orthogonal, and |v1| = -1, |v2|=1
        # That is, v1 is unit timelike and v2 is unit spacelike.
        # This means that cosh(x)*v1 + sinh(x)*v2 is unit timelike.
        hyperbola = cosh(x)*v1 + sinh(x)*v2
        endtime = arcsinh(v2_ldot_u2)
        from sage.plot.plot3d.all import parametric_plot3d
        pic = parametric_plot3d(hyperbola, (x, 0, endtime), **graphics_options)
        if show_hyperboloid:
            pic += self._model.get_background_graphic()
        return pic

    show = deprecated_function_alias(20530, plot)
