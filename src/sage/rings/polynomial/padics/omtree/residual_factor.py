r"""
Irreducible factors of residual polynomials of higher order Newton Polygons

Okutsu-Montes/Ore-Mac Lane-trees (OM trees) encode information about the
factorization of a polynomial `\Phi` in the form of a tree. Each node (called a
:class:`frame.Frame`) has a (higher order) Newton polygon associated to it.
To each :class:`segment.Segment` of such a Newton polygon one associates a
polynomial (the :meth:`segment.Segment.residual_polynomial`) over the residue
field (of the previous frame). Each irreducible factor (represented by an
:class:`ResidualFactor`) of these residual polynomials corresponds to a
factor of `\Phi` approximated in this frame. If such an irreducible factor over
the residue field is not linear, it encodes a residue field extension, i.e.,
inertia introduced by the corresponding factor of `\Phi`.

AUTHORS:

- Brian Sinclair and Sebastian Pauli (2012-02-22): initial version

- Julian Rüth (2017-07-19): extended documentation

EXAMPLES:

Consider the following OM tree with respect to `\Z_2`::

    sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
    sage: R.<x> = ZpFM(2, 20, 'terse')[]
    sage: T = OMTree(x^4 + 20*x^3 + 44*x^2 + 80*x + 1040)

This OM tree consist of three frames::

    sage: T.leaves()
    [Frame with phi (1 + O(2^20))*x^4 + (20 + O(2^20))*x^3 + (44 + O(2^20))*x^2 + (80 + O(2^20))*x + (1040 + O(2^20))]
    sage: leaf = T.leaves()[0]
    sage: leaf.prev_frame()
    Frame with phi (1 + O(2^20))*x^2 + (2 + O(2^20))*x + (4 + O(2^20))
    sage: leaf.prev_frame().prev_frame()
    Frame with phi (1 + O(2^20))*x + (0 + O(2^20))

The Newton polygon of the first frame has a residual polynomial that does not
factor into linear factors::

    sage: frame = leaf.prev_frame().prev_frame()
    sage: frame.polygon
    [Segment of length 4 and slope 1]
    sage: segment = frame.polygon[0]
    sage: segment.residual_polynomial().factor()
    (z^2 + z + 1)^2

This is an residual factor of the next frame in the tree which corresponds to
a degree 2 residue field extension::

    sage: frame = leaf.prev_frame()
    sage: frame.prev
    z^2 + z + 1

.. TODO::

    Most of the code in this class only exists to work around the fact that
    there is no proper implementation of towers finite fields (as of 2017-07.)
    Once towers of finite fields are in better shape, this file can probably go
    away and replaced by just an irreducible polynomial over the residue field.
    The only really important method in this file is probably
    :meth:`ResidualFactor.next_frame`.

"""
#*****************************************************************************
#       Copyright (C) 2012-2017 Brian Sinclair <bsinclai@gmail.com>
#                               Sebastian Pauli <s_pauli@uncg.edu>
#                          2017 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.padics.omtree.frameelt import FrameElt
from sage.rings.infinity import infinity

class ResidualFactor:
    r"""
    An irreducible factor of an :meth:`segment.Segment.residual_polynomial` of
    a :class:`segment.Segment` of a higher order Newton polygon.

    INPUT:

    - ``segment`` -- the segment of the Newton polygon of whose residual
      polygon this is a factor

    - ``rho`` -- the irreducible factor over the residue field

    - ``multiplicity`` -- the multiplicity of the factor in the residual polygon

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
        sage: R.<x> = ZpFM(2, 20, 'terse')[]
        sage: t = OMTree(x^4 + 20*x^3 + 44*x^2 + 80*x + 1040)
        sage: factor = t.leaves()[0].prev; factor
        z0^2 + a0*z0 + 1

    TESTS::

        sage: from sage.rings.polynomial.padics.omtree.residual_factor import ResidualFactor
        sage: isinstance(factor, ResidualFactor)
        True

    """
    def __init__(self, segment, rho, multiplicity):
        """
        .. TODO::

            Fields such as ``segments``, ``rho``, ``multiplicity``, ``Fplus``
            could be hidden behind a method/property so they get a proper
            docstring. Most other fields should probably be hidden completely
            by prepending a ``_`` to their names.

        TESTS::

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: R.<x> = ZpFM(2, 20, 'terse')[]
            sage: T = OMTree(x^4 + 20*x^3 + 44*x^2 + 80*x + 1040)
            sage: factor = T.leaves()[0].prev
            sage: TestSuite(factor).run()

        """
        self.segment = segment
        self.rho = rho
        self.multiplicity = multiplicity
        # the degree of inertia contributed by this factor
        self.Fplus = self.rho.degree()

        fr = self.segment.frame
        if self.segment.frame.is_root():
            # In the first frame, so FFbase is the residue class field of O
            self.FFbase = fr.R
        else:
            self.FFbase = fr.prev.FF

        if self.Fplus == 1:
            self.FF = self.FFbase
            self.FFz = PolynomialRing(self.FF, 'z' + str(fr.depth))
            # rho is linear: delta is the root of rho
            self.delta = self.rho.roots()[0][0]
        else:
            # the residue field as an absolute field
            self.FF = GF(self.FFbase.order()**self.Fplus, 'a' + str(fr.depth))
            self.FFz = PolynomialRing(self.FF, 'z' + str(fr.depth))
            self.FFbase_gamma = (self.FFz(self.FFbase.modulus())).roots()[0][0]
            FFrho = self.FFz([self.FFbase_elt_to_FF(a) for a in list(rho)])
            self.gamma = FFrho.roots()[0][0]
            basis = [(self.gamma**j*self.FFbase_gamma**i).polynomial() for j in range(0, self.Fplus) for i in range(0, self.FFbase.degree())]
            self.basis_trans_mat = Matrix([self.FF(b)._vector_() for b in basis])

    def __eq__(self, other):
        r"""
        Return whether this factor equals ``other``.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: R.<x> = ZpFM(2, 20, 'terse')[]
            sage: T = OMTree(x^4 + 20*x^3 + 44*x^2 + 80*x + 1040)
            sage: T.leaves()[0].prev == T.leaves()[0].prev
            True

        """
        return type(self) == type(other) and self.segment == other.segment and self.rho == other.rho and self.multiplicity == other.multiplicity

    def __ne__(self, other):
        r"""
        Return whether this factor does not equal ``other``.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: R.<x> = ZpFM(2, 20, 'terse')[]
            sage: T = OMTree(x^4 + 20*x^3 + 44*x^2 + 80*x + 1040)
            sage: T.leaves()[0].prev != T.leaves()[0].prev
            False

        """
        return not (self == other)

    def __hash__(self):
        r"""
        Return a hash value for this factor.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: k = ZpFM(2, 20, 'terse')
            sage: kx.<x> = k[]
            sage: t = OMTree(x^4 + 20*x^3 + 44*x^2 + 80*x + 1040)
            sage: factor = t.leaves()[0].prev
            sage: hash(factor) == hash(factor)
            True

        """
        return hash((self.segment, self.rho, self.multiplicity))

    def FF_elt_to_FFbase_vector(self, a):
        r"""
        Represent ``a`` as a vector over the ground residue field.
        
        INPUT:

        - ``a`` -- element of our extended residue field

        OUTPUT:

        A list representing ``a`` over the ground field of the latest extension

        .. TODO::

            This methods works around missing functionality for towers of
            finite fields.

        EXAMPLES:

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: R.<x> = ZpFM(2, 20, 'terse')[]
            sage: t = OMTree(x^4 + 20*x^3 + 44*x^2 + 80*x + 1040).leaves()[0].prev_frame()
            sage: t.prev
            z^2 + z + 1
            sage: t.polygon[0].factors[0]
            z0^2 + a0*z0 + 1

        We take elements in the different finite fields and represent them as
        vectors over their base residue field::

            sage: K.<a0> = t.prev.FF; K
            Finite Field in a0 of size 2^2
            sage: t.prev.FF_elt_to_FFbase_vector(a0 + 1)
            [1, 1]
            sage: L.<a1> = t.polygon[0].factors[0].FF;L
            Finite Field in a1 of size 2^4
            sage: t.polygon[0].factors[0].FF_elt_to_FFbase_vector(a1)
            [1, a0 + 1]

        """
        if self.segment.frame.is_root() and self.Fplus == 1:
            return a
        elif self.Fplus == 1:
            return self.segment.frame.prev.FF_elt_to_FFbase_vector(a)
        else:
            basedeg = self.FFbase.degree()
            avec = self.FF(a)._vector_()
            svector = self.basis_trans_mat.solve_left(Matrix(self.FF.prime_subfield(), avec))
            s_list = svector.list()
            s_split = [ s_list[i*basedeg:(i + 1)*basedeg] for i in range(0, self.Fplus)]
            s = [sum([ss[i]*self.FFbase.gen()**i for i in range(0, len(ss))]) for ss in s_split]
            return s

    def FFbase_elt_to_FF(self, b):
        r"""
        Lift an element up from the previous residue field to the current
        extended residue field.

        INPUT:

        - ``b`` -- element in the previous residue field.

        .. TODO::

            This methods works around missing functionality for towers of
            finite fields.

        EXAMPLES:

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: R.<x> = ZpFM(2, 20, 'terse')[]
            sage: t = OMTree(x^4 + 20*x^3 + 44*x^2 + 80*x + 1040).leaves()[0].prev_frame()

        We take elements in the different finite fields and lift them to the
        next residue field upward in the extension tower::

            sage: K.<a0> = t.prev.FF;K
            Finite Field in a0 of size 2^2
            sage: L.<a1> = t.polygon[0].factors[0].FF;L
            Finite Field in a1 of size 2^4
            sage: t.prev.FFbase_elt_to_FF(1)
            1
            sage: t.polygon[0].factors[0].FFbase_elt_to_FF(a0 + 1)
            a1^2 + a1 + 1

        """
        fr = self.segment.frame
        if fr.is_root() and self.Fplus == 1:
            return b
        elif self.Fplus == 1:
            return fr.prev.FFbase_elt_to_FF(b)
        elif fr.F == 1 and self.FFbase.is_prime_field():
            return b * self.FFbase_gamma
        else:
            bvec = b._vector_()
            return sum([ bvec[i]*self.FFbase_gamma**i for i in range(len(bvec))])

    def __repr__(self):
        r"""
        Return a printable representation of this factor.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: k = ZpFM(2, 20, 'terse')
            sage: kx.<x> = k[]
            sage: t = OMTree(x^4 + 20*x^3 + 44*x^2 + 80*x + 1040).leaves()[0].prev_frame()
            sage: t.prev.__repr__()
            'z^2 + z + 1'
            sage: t.polygon[0].factors[0].__repr__()
            'z0^2 + a0*z0 + 1'

        """
        return repr(self.rho)

    def lift(self, delta):
        """
        Return a :class:`sage.rings.polynomial.padics.omtree.frameelt.FrameElt`
        representation of a lift of the residue field element ``delta``.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: k = ZpFM(2, 20, 'terse')
            sage: kx.<x> = k[]
            sage: t = OMTree(x^4 + 20*x^3 + 44*x^2 + 80*x + 1040).leaves()[0].prev_frame()
            sage: K.<a0> = t.prev.FF;K
            Finite Field in a0 of size 2^2
            sage: t.polygon[0].factors[0].lift(a0 + 1)
            [[1*2^0]phi1^0, [1*2^-1]phi1^1]

        """
        fr = self.segment.frame
        if fr.F == 1:
            return FrameElt(fr, fr.Ox(delta))
        elif fr.prev.Fplus == 1:
            return FrameElt(fr, fr.prev.lift(delta), this_exp=0)
        else:
            dvec = fr.prev.FF_elt_to_FFbase_vector(delta)
            return sum([fr.prev.gamma_frameelt**i*FrameElt(fr, fr.prev.lift(dvec[i]), this_exp=0) for i in range(len(dvec)) if dvec[i] != 0])

    def next_frame(self, length=infinity):
        r"""
        Produce the child
        :class:`sage.rings.polynomial.padics.omtree.frame.Frame` corresponding
        to this factor in the OM tree.

        This method generates a new frame with this factor as ``previous`` and
        seeds it with a new approximation with strictly greater valuation than
        the current one.

        INPUT:

        - ``length`` -- an integer or ``infinity`` (default: ``infinity``); the
          length of the segment generating this factor.  This is used to reduce
          the total number of quotient with remainder operations needed in the
          resulting frame.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.frame import Frame
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x^32 + 16)
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen());f
            Frame with phi (1 + O(2^20))*x + (0 + O(2^20))
            sage: f = f.polygon[0].factors[0].next_frame();f
            Frame with phi (1 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (1048574 + O(2^20))
            sage: f = f.polygon[0].factors[0].next_frame();f
            Frame with phi (1 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (1048574 + O(2^20))*x^2 + (0 + O(2^20))*x + (1048574 + O(2^20))
            sage: f = f.polygon[0].factors[0].next_frame();f
            Frame with phi (1 + O(2^20))*x^16 + (0 + O(2^20))*x^15 + (0 + O(2^20))*x^14 + (0 + O(2^20))*x^13 + (0 + O(2^20))*x^12 + (0 + O(2^20))*x^11 + (1048572 + O(2^20))*x^10 + (0 + O(2^20))*x^9 + (1048572 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (1048572 + O(2^20))*x^5 + (4 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (8 + O(2^20))*x^2 + (0 + O(2^20))*x + (4 + O(2^20))

        """
        from frame import Frame
        fr = self.segment.frame
        if self.segment.slope is infinity:
            next = Frame(fr.Phi(), self, fr.iteration)
            self.next = next
            next.seed(fr.phi, length=length)
            return next
        if self.Fplus == 1 and self.segment.Eplus == 1:
            next = Frame(fr.Phi(), fr.prev, fr.iteration)
        else:
            next = Frame(fr.Phi(), self, fr.iteration)
        self.next = next
        self.gamma_frameelt = FrameElt(next, self.segment.psi**-1, self.segment.Eplus)
        if self.Fplus == 1 and fr.F == 1:
            next_phi = fr.phi**self.segment.Eplus - (self.segment.psi.polynomial() * fr.Ox(self.delta))
            self.reduce_elt = FrameElt(next, self.segment.psi * self.lift(self.delta), 0)
            next.seed(next_phi, length=length)
        elif self.Fplus == 1 and self.segment.Eplus == 1:
            delta_elt = self.lift(self.delta)
            next_phi_tail = self.segment.psi * delta_elt.reduce()
            next_phi = fr.phi - next_phi_tail.polynomial()
            self.reduce_elt = FrameElt(next, next_phi_tail, 0)
            next.seed(next_phi, length=length)
        else:
            lifted_rho_coeffs = [self.lift(r) for r in list(self.rho)]
            lifted_rho_coeffs_with_psi = [FrameElt(next, (self.segment.psi**(self.Fplus - i) * lifted_rho_coeffs[i]).reduce(), 0) for i in range(len(lifted_rho_coeffs))]
            phi_elt = FrameElt(next, fr.Ox(1), 1)
            next_phi_tail = sum([phi_elt**(self.segment.Eplus * i) * lifted_rho_coeffs_with_psi[i] for i in range(len(lifted_rho_coeffs_with_psi) - 1)])
            next_phi = (phi_elt**(self.segment.Eplus * self.Fplus) + next_phi_tail).polynomial()
            self.reduce_elt = FrameElt(next) - next_phi_tail
            next.seed(next_phi, length=length)
        return next
