r"""
Frames of Okutsu-Montes/Ore-Mac Lane-trees (OM trees)

For a polynomial `\Phi` over a local field, OM trees provide information about
factorization. Each node of such a tree is represented by a :class:`Frame`.
Starting from a single root, whenever the tree branches, new factors of `\Phi`
have been determined.

Each frame contains a polynomial `\phi` which is an approximate factor of
`\Phi`. Getting closer to the leaves of the tree, these approximations are
getting closer to the actual factors of `\Phi`.

AUTHORS:

- Brian Sinclair (2012-02-22): initial version

EXAMPLES:

Consider the following OM tree with respect to `\Z_2`::

    sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
    sage: R.<x> = ZpFM(2, 20, 'terse')[]
    sage: T = OMTree(x^32 + 16)

This OM tree consist of four frames::

    sage: T.leaves()
    [Frame with phi (1 + O(2^20))*x^32 + (0 + O(2^20))*x^31 + (0 + O(2^20))*x^30 + (0 + O(2^20))*x^29 + (0 + O(2^20))*x^28 + (0 + O(2^20))*x^27 + (1048568 + O(2^20))*x^26 + (0 + O(2^20))*x^25 + (1048568 + O(2^20))*x^24 + (0 + O(2^20))*x^23 + (0 + O(2^20))*x^22 + (1048568 + O(2^20))*x^21 + (24 + O(2^20))*x^20 + (0 + O(2^20))*x^19 + (48 + O(2^20))*x^18 + (0 + O(2^20))*x^17 + (24 + O(2^20))*x^16 + (32 + O(2^20))*x^15 + (1048544 + O(2^20))*x^14 + (32 + O(2^20))*x^13 + (1048480 + O(2^20))*x^12 + (0 + O(2^20))*x^11 + (1048496 + O(2^20))*x^10 + (1048544 + O(2^20))*x^9 + (1048544 + O(2^20))*x^8 + (1048512 + O(2^20))*x^7 + (64 + O(2^20))*x^6 + (1048544 + O(2^20))*x^5 + (96 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (96 + O(2^20))*x^2 + (0 + O(2^20))*x + (48 + O(2^20))]
    sage: leaf = T.leaves()[0]
    sage: leaf.prev_frame()
    Frame with phi (1 + O(2^20))*x^16 + (0 + O(2^20))*x^15 + (0 + O(2^20))*x^14 + (0 + O(2^20))*x^13 + (0 + O(2^20))*x^12 + (0 + O(2^20))*x^11 + (1048572 + O(2^20))*x^10 + (0 + O(2^20))*x^9 + (1048572 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (1048572 + O(2^20))*x^5 + (4 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (8 + O(2^20))*x^2 + (0 + O(2^20))*x + (4 + O(2^20))
    sage: leaf.prev_frame().prev_frame()
    Frame with phi (1 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (1048574 + O(2^20))*x^2 + (0 + O(2^20))*x + (1048574 + O(2^20))
    sage: leaf.prev_frame().prev_frame().prev_frame()
    Frame with phi (1 + O(2^20))*x + (0 + O(2^20))

Since there are no branch points in the tree, the polynomial `x^32 + 16` is
irreducible over `\Z_2`.

"""
#*****************************************************************************
#       Copyright (C) 2012-2017 Brian Sinclair <bsinclai@gmail.com>
#                          2017 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.padics.omtree.frameelt import FrameElt, FrameEltTerm
from sage.rings.polynomial.padics.omtree.segment import Segment
from sage.rings.infinity import infinity
from sage.structure.sage_object import SageObject

class Frame(SageObject):
    r"""
    A node in an OM tree.

    .. NOTE::

        Frames need to be seeded by calling :meth:`seed` with an approximation
        `\phi`.

    INPUT:

    - ``Phi`` - the local field polynomial being factored

    - ``prev`` - a :class:`residual_factor.ResidualFactor` (default: ``None``); the residual 
      factor from which this frame was created

    - ``iteration_count`` - Integer (default: ``0``); the number of previous
      iterations of the algorithm run to this point. Essentially this is the
      length of the path to the root of the OM tree but also accounts for
      intermediate frames that have been discarded.

    SEEALSO::

        :meth:`ResidualFactor.next_frame`

    EXAMPLES:

    Creating a frame leaves it unseeded, not having an approximation::

        sage: from sage.rings.polynomial.padics.omtree.frame import *
        sage: Phi = ZpFM(2, 20, 'terse')['x'](x^32 + 16)
        sage: f = Frame(Phi)
        sage: f
        Unseeded Frame regarding (1 + O(2^20))*x^32 + ... + (16 + O(2^20))

    A frame needs to be seeded to compute intermediate values::

        sage: f.seed(Phi.parent().gen())
        sage: f
        Frame with phi (1 + O(2^20))*x + (0 + O(2^20))

    A seeded frame has a newton polygon as a list of :class:`segment.Segment`
    objects::

        sage: f.polygon
        [Segment of length 32 and slope 1/8]

    And additional data associated to it::

        sage: f.polygon[0].psi.polynomial()
        (2 + O(2^20))
        sage: f.polygon[0].factors
        [z + 1]

    New frames in the tree can be created from each residual factor::

        sage: f.polygon[0].factors[0].next_frame()
        Frame with phi (1 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (1048574 + O(2^20))

    Notice that frames created in this way are seeded.

    """
    def __init__(self, Phi, prev=None, iteration_count=0):
        r"""
        .. TODO::

            Fields in this class should be hidden behind methods/properties so
            they get proper docstrings or they should be hidden by prepending a
            ``_`` to their name.

        TESTS::

            sage: from sage.rings.polynomial.padics.omtree.frame import *
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x^32 + 16)
            sage: f = Frame(Phi)
            sage: TestSuite(f).run()

        """
        self._Phi = Phi
        self.prev = prev
        self.O = Phi.base_ring()
        self.Ox = Phi.parent()
        self.R = self.O.residue_class_field()
        self.Rz = PolynomialRing(self.R, names='z')
        self.phi = None
        self.iteration = iteration_count + 1
        if self.is_root(): # that is self.prev is None
            self.E = 1
            self.F = 1
            self.depth = 0
        else:
            self.E = self.prev_frame().E * self.prev.segment.Eplus
            self.F = self.prev_frame().F * self.prev.Fplus
            self.depth = self.prev_frame().depth + 1

    def __eq__(self, other):
        """
        Check whether self is equal to other

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.frame import *
            sage: Phi = ZpFM(2,20,'terse')['x'](x^32+16)
            sage: f = Frame(Phi)
            sage: g = Frame(Phi, iteration_count=4)
            sage: f == g
            False
        """
        return type(self) == type(other) and self.Phi() == other.Phi() and self.phi == other.phi and self.iteration == other.iteration

    def seed(self, phi, length=infinity):
        r"""
        Seed all of the intermediate values of the frame based on the new
        approximation to a factor ``phi``.

        In seeding, the frame will compute a ``phi``-expansion of `\Phi` and
        the newton polygon of higher order. Then the Newton polygon will be
        constructed with its :class:`segment.Segment` and their residual
        polynomials.

        INPUT:

        - ``phi`` - an approximation to a factor of `\Phi`

        - ``length`` - Integer or ``infinity`` (default: ``infinity``); the
          length of the segment that resulted in this frame.  This provides an
          extra stopping point for finding phi-expansions if finite as we only
          need up to ``length`` terms.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.frame import *
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x^32 + 16)
            sage: f = Frame(Phi); f
            Unseeded Frame regarding (1 + O(2^20))*x^32 + ... + (16 + O(2^20))
            sage: f.phi is None
            True
            sage: f.seed(Phi.parent().gen()); f
            Frame with phi (1 + O(2^20))*x + (0 + O(2^20))
            sage: f.phi
            (1 + O(2^20))*x + (0 + O(2^20))

        """
        self.phi = phi

        # Construct the phi expansion of Phi
        self._phi_expansion_as_polynomials = []
        if self.phi.degree() > self.Phi().degree():
            self._phi_expansion_as_polynomials = [self.Phi()]
        q, r = self.Phi().quo_rem(self.phi)
        self._phi_expansion_as_polynomials.append(r)
        while q != 0 and length > len(self._phi_expansion_as_polynomials):
            q, r = q.quo_rem(self.phi)
            self._phi_expansion_as_polynomials.append(r)
        self._phi_expansion_as_elts = [FrameElt(self, a) for a in self._phi_expansion_as_polynomials]
        self.polygon = self._newton_polygon([e.valuation() for e in self._phi_expansion_as_elts]) # list of segments

    def find_psi(self, val):
        r"""
        Find a polynomial (as a :class:`frameelt.FrameElt`) with given valuation.

        INPUT:

        - ``val`` - a rational; the desired valuation. The denominator of
          ``val`` must divide the current level of ramification (``E``).

        OUTPUT:

        a FrameElt with respect to the current frame with valuation ``val``

        EXAMPLES:

            sage: from sage.rings.polynomial.padics.omtree.frame import *
            sage: Phi = ZpFM(2, 20, print_mode='terse')['x'](x^16 + 16)
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen())
            sage: f = f.polygon[0].factors[0].next_frame(); f
            Frame with phi (1 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (1048574 + O(2^20))

        We get a valid :class:`frameelt.FrameElt` with integer exponents as
        long as the denominator of ``val`` divides the current ramification::

            sage: f.E
            4
            sage: f.prev.segment.slope
            1/4
            sage: f.find_psi(7/4)
            [[1*2^1]phi1^3]
            sage: f.find_psi(7/4).polynomial()
            (2 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (0 + O(2^20))

        If the denominator does not divide the ramification, then we cannot
        construct a polynomial of this valuation::

            sage: f.find_psi(3/8)
            Traceback (most recent call last):
            ...
            ValueError: Denominator of given valuation does not divide E

        """

        if self.E % val.denom() != 0:
            raise ValueError("Denominator of given valuation does not divide E")
        psielt = FrameElt(self)
        if self.prev is None:
            psielt.terms = [FrameEltTerm(psielt, self.O(1), val)]
        else:
            vphi = self.prev.segment.slope
            d = self.prev_frame().E
            vprime = val * d
            e = vphi * d
            psimod = e.denom()
            if psimod == 1:
                s = 0
            else:
                s = vprime / e
                s = s % psimod
                val -= s * vphi
            psielt.terms = [FrameEltTerm(psielt, self.prev_frame().find_psi(val), s)]
        return psielt

    def root(self):
        r"""
        Return a root frame of the OM tree.

        .. NOTE::

            The leaves for a polynomial's OM Tree may have different roots.
            Becuase optimal OM Trees remove "improvement frames" (those for
            which neither ramificaiton or inertia is found), an OM Tree may
            technically be a forest.

        .. SEEALSO::

            :meth:`omtree.OMTree.roots`

        EXAMPLES:

        In this example, the root is the default initial frame::

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: from sage.rings.polynomial.padics.omtree.frame import Frame
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x^32 + 16)
            sage: T = OMTree(Phi)
            sage: T.leaves()[0].root()
            Frame with phi (1 + O(2^20))*x + (0 + O(2^20))
            sage: Fr = Frame(Phi)
            sage: Fr.seed(Fr.Ox.gen())
            sage: Fr.root() == Fr
            True
            sage: T.leaves()[0].root() == Fr
            True

        This may not be the case, and we may end up with a forest::

            sage: R.<c> = ZqFM(125, prec = 30, print_mode='terse')
            sage: Rz.<z>=R[]
            sage: g=(z^3 + 2)^5 + 5
            sage: om=OMTree(g)
            sage: om.leaves()[0].root()
            Frame with phi (1 + O(5^30))*z + (931322574615478515623 + O(5^30))
            sage: om.leaves()[1].root()
            Frame with phi (1 + O(5^30))*z + (0 + O(5^30))

        """
        if self.is_root():
            return self
        else:
            return self.prev_frame().root()

    def _newton_polygon(self, a):
        r"""
        Compute the newton polygon of higher order of `\Phi` with respect to
        the valuations of the phi-expansion of `\Phi`.

        This method finds the slopes of each :class:`segment.Segment` as well
        as the points on the line of each segment.

        .. TODO::

            This method should be replaced with the more generic newton polygon
            implementation provided by
            :class:`sage.geometry.newton_polygon.NewtonPolygon`.

        INPUT:

        - ``a`` -- a list; the valuations of the coefficients in the
          phi-expansion of `\Phi`

        OUTPUT:

        a list of segments

        ALGORITHM:

        Monotone chain algorithm,s ee
        http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.frame import *
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x*(x + 1)*(x + 2))
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen())
            sage: v = [e.valuation() for e in f._phi_expansion_as_elts];v
            [20, 1, 0, 0]
            sage: f._newton_polygon(v)
            [Segment of length 1 and slope +Infinity,
             Segment of length 1 and slope 1,
             Segment of length 1 and slope 0]

        """
        def cross(o, a, b):
            # 2D cross product of the vectors oa and ob
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        def find_slope(a, b):
            # the negative of the slope through the points a and b
            # (note that we negate the slope)
            return (a[1]-b[1]) / (b[0]-a[0])

        lower = [(0, a[0])]
        segments = []
        for i in range(1, len(a)):
            p = (i, a[i])
            # We check cross < 0 since we want to retain points on the boundary.
            while len(lower) >= 2 and cross(lower[-2], lower[-1], p) < 0:
                lower.pop()
            lower.append(p)
        if len(lower) <= 1:
            raise ValueError("Not enough vertices")

        # If phi divides Phi then the first segment should have infinite slope.
        if self.phi_divides_Phi():
            for c in xrange(1, len(lower)):
                if lower[c][1] < self.O.precision_cap():
                    break
            else:
                raise ValueError("Entire polygon above precision cap")
            segments.append(Segment(self, infinity, [(0, infinity), lower[c]]))
            lower = lower[c:]
            if len(lower) <= 1:
                return segments

        slope = find_slope(lower[0], lower[1])
        verts = [lower[0], lower[1]]
        for c in xrange(1, len(lower)-1):
            new_slope = find_slope(lower[c], lower[c + 1])
            if new_slope == slope:
                verts.append(lower[c + 1])
            else:
                segments.append(Segment(self, slope, verts))
                slope = new_slope
                verts = [lower[c], lower[c + 1]]
        segments.append(Segment(self, slope, verts))
        return segments

    def prev_frame(self):
        r"""
        Return the parent of this frame in the OM tree.

        Returns ``None`` if this is a root.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.frame import *
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x^32 + 16)
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen())
            sage: f
            Frame with phi (1 + O(2^20))*x + (0 + O(2^20))
            sage: len(f.polygon)
            1
            sage: len(f.polygon[0].factors)
            1
            sage: f = f.polygon[0].factors[0].next_frame()
            sage: f
            Frame with phi (1 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (1048574 + O(2^20))
            sage: f.prev_frame()
            Frame with phi (1 + O(2^20))*x + (0 + O(2^20))

        """
        if self.prev is None:
            return None
        else:
            return self.prev.segment.frame

    def is_root(self):
        r"""
        Return whether this is a root of the OM tree.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.frame import *
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x^32 + 16)
            sage: f = Frame(Phi)
            sage: f.is_root()
            True
            sage: f.seed(Phi.parent().gen())
            sage: f = f.polygon[0].factors[0].next_frame()
            sage: f.is_root()
            False

        """
        return self.prev is None

    def Phi(self):
        r"""
        Return `\Phi`, the local field polynomial being factored.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.frame import *
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x^3 + x)
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen())
            sage: f.Phi()
            (1 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (1 + O(2^20))*x + (0 + O(2^20))

        """
        return self._Phi


    def phi_divides_Phi(self):
        r"""
        Return whether the approximation ``phi`` divides `\Phi`.

        This is done by checking if the constant term of the ``phi``-expansion
        of `\Phi` is zero.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.frame import *
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x^3 + x)
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen())
            sage: f.Phi()
            (1 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (1 + O(2^20))*x + (0 + O(2^20))
            sage: f.phi
            (1 + O(2^20))*x + (0 + O(2^20))
            sage: f.phi_divides_Phi()
            True

        """
        return self._phi_expansion_as_polynomials[0] == 0

    def single_factor_lift(self):
        r"""
        Lift a frame to a factor of the polynomial it approximates.

        INPUT:

        - ``frame`` -- a frame that is the leaf of an OM tree.

        OUTPUT:

        a factor of the polynomial referred to by the input frame

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.frame import Frame
            sage: Kx.<x> = PolynomialRing(ZpFM(3, 20, 'terse'))
            sage: f = (x**3 + x-1)*(x**2 + 1)
            sage: fr = Frame(f)
            sage: fr.seed(x)
            sage: fr = fr.polygon[0].factors[0].next_frame()
            sage: fact = fr.single_factor_lift() ; fact
            (1 + O(3^20))*x + (904752403 + O(3^20))
            sage: f % fact
            0

        REFERENCES:

        - [GNP2012]

        """
        if self.phi_divides_Phi():
            return self.phi

        def _reduce(poly, phi, d):
            # returns poly mod phi and simplifies the denominator of poly
            poly = poly % phi
            if d != 0:
                g = min([d] + [p.valuation() for p in poly])
                if g > 0:
                    poly = poly.parent( [p >> g for p in poly] )
                    d = d - g
            return poly, d

        def _move_elt(a, S):
            # move the element a of R to S which is a changed R
            return S(a.polynomial().list())

        def _move_poly(f, Sy):
            # move the polynomial f over R to a polynomial over S which is a changed R
            return Sy([_move_elt(a, Sy.base()) for a in f.coefficients(sparse=False)])

        prec = self.O.precision_cap()
        LiftRing = self.O.change(prec=2*prec, type="fixed-mod")
        Lifty = PolynomialRing(LiftRing, names='y')
        Phi = _move_poly(self.Phi(), Lifty)
        phi = _move_poly(self.phi, Lifty)

        a0,a1 = self._phi_expansion_as_elts[0:2]

        Psi = self.find_psi(-a1.valuation())
        A0 = Psi * a0
        A1 = Psi * a1

        Psi,Psiden = Psi.polynomial(True)
        Psi = _move_poly(Psi, Lifty)

        C1inv = self.polygon[0].factors[0].lift(1/(A1.residue()))
        C1inv,C1invden = C1inv.polynomial(True)
        C1inv = _move_poly(C1inv, Lifty)
        C1inv,C1invden = _reduce(C1inv, phi, C1invden)

        A0,A0den = A0.polynomial(True)
        A0,A0den = _reduce(_move_poly(A0, Lifty), phi, A0den)

        C,Cden = _reduce(A0*C1inv, phi, A0den + C1invden)
        phi = (phi + C)

        h = 2
        oldphi = None
        while h < prec and phi != oldphi:
            oldphi = phi
            C1,C0 = Phi.quo_rem(phi)

            C0,C0den = _reduce((Psi*C0), phi, Psiden)
            C1,C1den = _reduce((Psi*C1), phi, Psiden)

            x1,x1den = _reduce((LiftRing(2)<<(C1den + C1invden))-C1*C1inv, phi, C1den + C1invden)
            C1inv,C1invden = _reduce(C1inv*x1, phi, C1invden + x1den)

            C,Cden = _reduce((C0*C1inv), phi, C0den + C1invden)

            phi = (phi + C)
            h = 2 * h
        return _move_poly(phi, self.Ox)

    def __repr__(self):
        """
        Return a printable representation of this frame.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.frame import *
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x^32 + 16)
            sage: f = Frame(Phi)
            sage: f.__repr__()
            'Unseeded Frame regarding (1 + O(2^20))*x^32 + ... + (16 + O(2^20))'
            sage: f.seed(Phi.parent().gen())
            sage: f.__repr__()
            'Frame with phi (1 + O(2^20))*x + (0 + O(2^20))'
        """
        if self.phi:
            return 'Frame with phi ' + repr(self.phi)
        else:
            return 'Unseeded Frame regarding ' + repr(self.Phi())
