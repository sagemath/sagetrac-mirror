# -*- coding: utf-8 -*-
"""
Interface to the GAP package GBNP

This module is an interface to the `GBNP
<https://gap-packages.github.io/gbnp/>`_ package, which provides algorithms for
computing Groebner bases of noncommutative polynomials with coefficients from a
field implemented in GAP and with respect to the "total degree first then
lexicographical" ordering.

It currently only implements a wrapper for a subset of GBNP. Minimal usage
example::

    sage: from sage.algebras.gbnp import GapFreeAlgebra
    sage: A.<a,b,c> = GapFreeAlgebra(QQ)    # optional - gap_packages
    sage: rels = [a^2-a, b^2-b, c^2-c, a+b+c-1]    # optional - gap_packages
    sage: I = A.ideal(rels)    # optional - gap_packages
    sage: I.is_groebner_basis()    # optional - gap_packages
    False
    sage: I1 = I.groebner_basis()    # optional - gap_packages
    sage: I1    # optional - gap_packages
    Twosided Ideal (-1 + a + b + c, -a + a^2, a*b, b*a, -b + b^2) of Free Algebra on 3 generators (a, b, c) over Rational Field
    sage: I1.is_groebner_basis()    # optional - gap_packages
    True
    sage: I1.reduce(a*b)    # optional - gap_packages
    0
    sage: I1.reduce(c^5)    # optional - gap_packages
    1 - a - b

There are examples where the algorithm to find Groebner basis does not converge::

    sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
    sage: rels = [x * y - y * x - y^2]    # optional - gap_packages
    sage: I = A.ideal(rels)    # optional - gap_packages
    sage: I.is_groebner_basis()    # optional - gap_packages
    False
    sage: I1 = I.groebner_basis(10)    # optional - gap_packages
    sage: I1.is_groebner_basis()    # optional - gap_packages
    False
    sage: I2 = I1.groebner_basis(10)    # optional - gap_packages
    sage: I1.is_groebner_basis()    # optional - gap_packages
    False

However, changing the order of the variables makes the algorithm in this
particular example converge::

    sage: A.<y, x> = GapFreeAlgebra(QQ)    # optional - gap_packages
    sage: rels = [x * y - y * x - y^2]    # optional - gap_packages
    sage: I = A.ideal(rels)    # optional - gap_packages
    sage: I.is_groebner_basis()    # optional - gap_packages
    True

The same example does converge over finite fields::

    sage: A3.<x, y> = GapFreeAlgebra(GF(3))    # optional - gap_packages
    sage: rels = [x * y - y * x - y^2]    # optional - gap_packages
    sage: I = A3.ideal(rels)    # optional - gap_packages
    sage: I.is_groebner_basis()    # optional - gap_packages
    False
    sage: I10 = I.groebner_basis(10)    # optional - gap_packages
    sage: I10.is_groebner_basis()    # optional - gap_packages
    True
    sage: A9.<x, y> = GapFreeAlgebra(GF(9))    # optional - gap_packages
    sage: rels = [x * y - y * x - y^2]    # optional - gap_packages
    sage: I = A9.ideal(rels)    # optional - gap_packages
    sage: I.is_groebner_basis()    # optional - gap_packages
    False
    sage: I10 = I.groebner_basis(10)    # optional - gap_packages
    sage: I10.is_groebner_basis()    # optional - gap_packages
    True

The Cola gene puzzle from
https://gap-packages.github.io/gbnp/doc/chapA.html#X7912E411867E5F8B::

    sage: A1.<A, C, G, T> = GapFreeAlgebra(QQ)    # optional - gap_packages
    sage: rels = [T*C*A*T-T, G*A*G-A*G, C*T*C-T*C, A*G*T*A-A, T*A*T-C*T]    # optional - gap_packages
    sage: I = A1.ideal(rels)    # optional - gap_packages
    sage: I1 = I.groebner_basis(10)    # optional - gap_packages
    sage: I1.is_groebner_basis()    # optional - gap_packages
    True
    sage: milk = T*A*G*C*T*A*G*C*T*A*G*C*T    # optional - gap_packages
    sage: I1.reduce(milk)    # optional - gap_packages
    T
    sage: cola = C*T*G*A*C*T*G*A*C*T    # optional - gap_packages
    sage: I1.reduce(cola)    # optional - gap_packages
    T

Working with quotient algebras::

    sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
    sage: rels = [x^2, y^2]    # optional - gap_packages
    sage: QA = A.quo(rels)    # optional - gap_packages
    sage: QA.is_finite_dimensional()    # optional - gap_packages
    False
    sage: QA.growth()    # optional - gap_packages
    1
    sage: QA.hilbert_series(10)    # optional - gap_packages
    [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

    sage: rels1 = [x^2]    # optional - gap_packages
    sage: QA1 = A.quo(rels1)    # optional - gap_packages
    sage: QA1.is_finite_dimensional()    # optional - gap_packages
    False
    sage: QA1.growth()    # optional - gap_packages
    'exponential growth'
    sage: QA1.hilbert_series(10)    # optional - gap_packages
    [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]

    sage: rels2 = [x^2, y^2, x*y*x*y]    # optional - gap_packages
    sage: QA2 = A.quo(rels2)    # optional - gap_packages
    sage: QA2.is_finite_dimensional()    # optional - gap_packages
    True
    sage: QA2.dim()    # optional - gap_packages
    8
    sage: QA2.hilbert_series(10)    # optional - gap_packages
    [1, 2, 2, 2, 1]
    sage: basis2 = QA2.get_basis()    # optional - gap_packages
    sage: basis2    # optional - gap_packages
    [1, x, y, x*y, y*x, x*y*x, y*x*y, y*x*y*x]
    sage: QA2.get_matrix(0, basis2)    # optional - gap_packages
    [0 1 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0]
    [0 0 0 0 1 0 0 0]
    [0 0 0 0 0 1 0 0]
    [0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 0]
    [0 0 0 0 0 0 0 1]
    [0 0 0 0 0 0 0 0]

Installation
------------

GBNP is one of the GAP packages included in the optional Sage package :ref:`gap_packages <spkg_gap_packages>`. To install it, run the following in the shell::

    $ sage -i gap_packages

It is also possible to manually install GBNP, following the `GAP documentation <https://www.gap-system.org/Manuals/doc/ref/chap76.html>`_. At least on Linux, first download the latest tarball at the time of writing::

    $ wget https://github.com/gap-packages/gbnp/archive/refs/tags/v1.0.5.tar.gz

and verify integrity::

    $ sha256sum v1.0.5.tar.gz
    c9914ad622213bb3fffbdafdd4292b2f4abb7d66f8bad64496a56921b04b1c14  v1.0.5.tar.gz

Unpack the tarball into the ``pkg`` directory of GAP::

    $ sage -sh
    (sage-sh) $ cd $SAGE_ROOT/local/share/gap/pkg/
    (sage-sh) $ tar -xvf /the/path/to/gbnp/v1.0.5.tar.gz
    (sage-sh) $ exit

To check that the package installed correctly, try loading it::

   $ sage
   sage: libgap.load_package("gbnp")    # optional - gap_packages
   true

Links and references:

- https://gap-packages.github.io/gbnp and the `documentation
  <https://gap-packages.github.io/gbnp/doc/chap0.html>`_ page. The old package page
  is `here <https://www.gap-system.org/Packages/gbnp.html>`_.

- `GBNP package repository at GitHub <https://github.com/gap-packages/gbnp>`_

- Arjeh M. Cohen and Dié .A.H. Gijsbers. Noncommutative groebner basis
  computations. Report, 2003, http://www.win.tue.nl/~amc/pub/grobner/gbnp.pdf

- Jan Willem Knopper. GBNP and Vector Enumeration. Internship report, 2004
  http://mathdox.org/gbnp/knopper.pdf
"""
#*****************************************************************************
#    Copyright (C) 2021 Guy Blachar, Tomer Bauer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gap.libgap import libgap
from sage.misc.misc_c import prod
from sage.algebras.free_algebra import FreeAlgebra_generic
from sage.rings.noncommutative_ideals import Ideal_nc
from sage.rings.quotient_ring import QuotientRing_nc
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.all import Matrix
from sage.features.gap import GapPackage


def _sage2gap(elem, gap_alg, s2g):
    """
    Convert an element of a free algebra into a GAP element by generators.

    This is an internal function.

    INPUT:

    - ``elem`` -- an element of a Sage free algebra over some field

    - ``gap_alg`` -- an algebra which is a GAP element

    - ``s2g`` -- a mapping of the *monoid generators* of the free algebra to the
      generators of ``gap_alg``

    OUTPUT:

    A ``GapElement`` which is a non-commutative polynomial representing
    ``elem``.

    EXAMPLES::

        sage: from sage.algebras.gbnp import _sage2gap, _gap2sage
        sage: A.<sa, sb> = FreeAlgebra(QQ)
        sage: msa, msb = A.monoid().gens()
        sage: elem = 3*sa^2 + 2/5*sb*sa*sb
        sage: GapA = libgap.FreeAssociativeAlgebraWithOne(libgap(QQ), ["a", "b"])
        sage: ga, gb = GapA.GeneratorsOfAlgebraWithOne()
        sage: gap_elem = _sage2gap(elem, GapA, {msa: ga, msb: gb}); gap_elem
        (3)*a^2+(2/5)*b*a*b
        sage: _gap2sage(libgap.GP2NP(gap_elem), A) == elem    # optional - gap_packages
        True
    """
    gap_elem = gap_alg.ZeroImmutable()
    if not elem.is_zero():
        free_coeff = 0
        if elem.trailing_monomial().is_one():
            free_coeff = elem.trailing_coefficient()
            elem -= elem.trailing_term()
        gap_elem += free_coeff * gap_alg.OneImmutable()
        for mon, coeff in elem:
            gap_elem += coeff * prod([s2g[x]**exponent for x, exponent in mon])
    return gap_elem


def _gap2sage(gap_l, sage_alg):
    """
    Convert a GBNP non-commutative polynomial to an element of a Sage algebra.

    This is an internal function.

    INPUT:

    - ``gap_l`` -- a ``GapElement`` which is in GBNP's non-commutative
      polynomial (NP) format. It is a list of two lists of the same length. The
      first of the two lists encodes the monomials, and the second list holds
      the coefficients of the base field. See the `GBNP
      documentation<https://gap-packages.github.io/gbnp/doc/chap2.html#X7FDF3E5E7F33D3A2>`_
      for the full description of the NP format.

    - ``sage_alg`` -- a Sage algebra, such as the free algebra. It needs to be
      generated by at least the number of generators of the algebra of
      ``elem``

    EXAMPLES::

        sage: from sage.algebras.gbnp import _sage2gap, _gap2sage
        sage: F = GF(5)
        sage: gap_elem = [[[2, 1, 2], [1, 1]], [libgap(F(3)), libgap(F(2))]]
        sage: A.<a, b> = FreeAlgebra(GF(5))
        sage: elem = _gap2sage(gap_elem, A); elem
        2*a^2 + 3*b*a*b

    TESTS::

        sage: GapA = libgap.FreeAssociativeAlgebraWithOne(libgap(F), ["a", "b"])
        sage: ga, gb = GapA.GeneratorsOfAlgebraWithOne()
        sage: ma, mb = A.monoid().gens()
        sage: f1 = _sage2gap(elem, GapA, {ma: ga, mb: gb})
        sage: f2 = libgap.NP2GP(gap_elem, GapA)   # optional - gap_packages
        sage: _gap2sage(libgap.GP2NP(f1 - f2), A)   # optional - gap_packages
        0
    """
    sage_gens = sage_alg.gens()
    sage_elem = sage_alg.zero()
    for vars_l, coeff in zip(*gap_l):
        sage_elem += coeff.sage() * prod([sage_gens[i - 1] for i in vars_l], sage_alg.one())
    return sage_elem


class GapIdeal(Ideal_nc):
    def __init__(self, *args, **kwds):
        """
        A non-commutative two-sided ideal implemented in GAP.

        INPUT:

        - ``ring`` -- the ring of the ideal. Should be a free algebra over some field

        - ``gens`` -- the generators of the ideal

        EXAMPLES::

            sage: from sage.algebras.gbnp import GapFreeAlgebra, GapIdeal
            sage: A.<x, y> = FreeAlgebra(QQ)    # optional - gap_packages
            sage: I = GapIdeal(A, [x*y - y*x])    # optional - gap_packages
            sage: I    # optional - gap_packages
            Twosided Ideal (x*y - y*x) of Free Algebra on 2 generators (x, y) over Rational Field

            sage: A4.<x, y> = FreeAlgebra(GF(4))    # optional - gap_packages
            sage: I = GapIdeal(A4, [x^2, y^3])    # optional - gap_packages
            sage: I    # optional - gap_packages
            Twosided Ideal (x^2, y^3) of Free Algebra on 2 generators (x, y) over Finite Field in z2 of size 2^2

        One can also construct GapIdeals using GapFreeAlgebra::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x])    # optional - gap_packages
            sage: I    # optional - gap_packages
            Twosided Ideal (x*y - y*x) of Free Algebra on 2 generators (x, y) over Rational Field
        """
        GapPackage("GBNP", spkg="gap_packages").require()
        libgap.load_package("GBNP")
        libgap.SetInfoLevel(libgap.InfoGBNP, 0)
        libgap.SetInfoLevel(libgap.InfoGBNPTime, 0)

        Ideal_nc.__init__(self, *args, **kwds)

        sage_alg = self.ring()
        sage_alg_gens = sage_alg.monoid().gens()
        self._gap_algebra = libgap.FreeAssociativeAlgebraWithOne(libgap(sage_alg.base_ring()),
                                        *[str(t) for t in sage_alg_gens])

        gap_gens = libgap.GeneratorsOfAlgebraWithOne(self._gap_algebra)
        self._s2g = dict(zip(sage_alg_gens, gap_gens))
        self._gap_rels = libgap.GP2NPList([_sage2gap(x, self._gap_algebra, self._s2g) for x in self.gens()])

    def groebner_basis(self, max_iters=10, strong=True):
        """
        Computes a Groebner basis for the ideal.

        INPUT:

        - ``max_iters`` (default: 10) -- the number of iterations for the
          Buchberger's Algorithm. If 0, the calculations will continue until it
          terminates (but might not terminate at all)

        - ``strong`` (default: ``True``) -- whether to compute a strong Groebner basis

        OUTPUT:

        A new ideal, with the Groebner basis of the given ideal as its generators.

        EXAMPLES::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x])    # optional - gap_packages
            sage: GI = I.groebner_basis()    # optional - gap_packages
            sage: GI    # optional - gap_packages
            Twosided Ideal (-x*y + y*x) of Free Algebra on 2 generators (x, y) over Rational Field

        Another example, in which the new ideal contains more generators::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y-1, x*y^2-x^2*y])    # optional - gap_packages
            sage: GI = I.groebner_basis()    # optional - gap_packages
            sage: GI    # optional - gap_packages
            Twosided Ideal (-x + y, -1 + x^2) of Free Algebra on 2 generators (x, y) over Rational Field

        Note that for some examples, the function may not converge unless you add the optional parameter ``max_iters``::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x - y^2])    # optional - gap_packages
            sage: GI = I.groebner_basis(10)    # optional - gap_packages
            sage: GI.is_groebner_basis()    # optional - gap_packages
            False
        """
        if strong:
            if max_iters == 0:
                GB = libgap.SGrobner(self._gap_rels)
            else:
                GB = libgap.SGrobner(self._gap_rels, max_iters)['G']
        else:
            if max_iters == 0:
                GB = libgap.Grobner(self._gap_rels)
            else:
                GB = libgap.Grobner(self._gap_rels, max_iters)['G']
        return GapIdeal(self.ring(), [_gap2sage(x, self.ring()) for x in GB])

    def is_groebner_basis(self, strong=True):
        """
        Return ``True`` if the generators of the given ideal form a Groebner
        basis, else ``False``.

        INPUT:

        - ``strong`` (default: ``True``) -- whether to check for a strong Groebner basis

        OUTPUT:

        ``True`` if the generators of the ideal form a Groebner basis, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x])    # optional - gap_packages
            sage: I.is_groebner_basis()    # optional - gap_packages
            True

            sage: A3.<x, y> = GapFreeAlgebra(GF(3))    # optional - gap_packages
            sage: I = A3.ideal([x*y - y*x - y^2])    # optional - gap_packages
            sage: I.is_groebner_basis()    # optional - gap_packages
            False
            sage: GI = I.groebner_basis()    # optional - gap_packages
            sage: GI.is_groebner_basis()    # optional - gap_packages
            True
        """
        if strong:
            return bool(libgap.IsStrongGrobnerBasis(self._gap_rels))
        return bool(libgap.IsGrobnerBasis(self._gap_rels))

    def reduce(self, elem, check=True):
        """
        Given an element of the free algebra, reduces it to a normal form using the given generators.

        INPUT:

        - ``elem`` -- an element of the free algebra
        - ``check`` (default: ``True``) -- if ``True``, checks whether the given
          generators form a Groebner basis, and if not compute a new Groebner
          basis (might not terminate!)

        OUTPUT:

        The reduced form of ``elem`` with respect to the given generators of
        the ideal.  Note that if the generators do not form a Groebner basis,
        this is not a normal form.

        EXAMPLES:

        We begin with the commutative polynomial ring, in which reducing simply
        means adding the powers of each element::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x])    # optional - gap_packages
            sage: I.reduce(y^3*x*y^2*x^2)    # optional - gap_packages
            x^3*y^5

        Another example, in which the original basis is not a Groebner basis::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x^2*y - x*y, x*y^2 + x^3])    # optional - gap_packages
            sage: I.is_groebner_basis()    # optional - gap_packages
            False
            sage: I.reduce(x^4 - x^3, check=False)    # optional - gap_packages
            -x^3 + x^4
            sage: I.reduce(x^4 - x^3, check=True)    # optional - gap_packages
            0
        """
        GB = self
        if check:
            GB = GB.groebner_basis(max_iters=0)

        gap_rels = GB._gap_rels

        gap_elem = libgap.GP2NP(_sage2gap(elem, GB._gap_algebra, GB._s2g))
        return _gap2sage(libgap.StrongNormalFormNP(gap_elem, gap_rels), self.ring())


class GapQuotientRing(QuotientRing_nc):
    def __init__(self, R, I, names=None, category=None):
        """
        A quotient ring of non-commutative rings implemented in GAP.

        INPUT:

        - ``R`` -- the cover ring

        - ``I`` -- the defining ideal of the quotient

        EXAMPLES:

        The commutative polynomial ring is a quotient of the free algebra::

            sage: from sage.algebras.gbnp import GapFreeAlgebra, GapQuotientRing
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x])    # optional - gap_packages
            sage: QA = GapQuotientRing(A, I)    # optional - gap_packages
            sage: QA    # optional - gap_packages
            Quotient of Free Algebra on 2 generators (x, y) over Rational Field by the ideal (x*y - y*x)

        It can also be defined using the functions of a GapFreeAlgebra::

            sage: A.quo(I)    # optional - gap_packages
            Quotient of Free Algebra on 2 generators (x, y) over Rational Field by the ideal (x*y - y*x)
        """
        GapPackage("GBNP", spkg="gap_packages").require()
        libgap.load_package("GBNP")
        libgap.SetInfoLevel(libgap.InfoGBNP, 0)
        libgap.SetInfoLevel(libgap.InfoGBNPTime, 0)

        QuotientRing_nc.__init__(self, R, I, names, category=category)

        sage_alg_gens = R.monoid().gens()
        self._gap_algebra = libgap.FreeAssociativeAlgebraWithOne(libgap(self.base_ring()),
                                *[str(t) for t in sage_alg_gens])

        gap_gens = libgap.GeneratorsOfAlgebraWithOne(self._gap_algebra)
        self._s2g = dict(zip(sage_alg_gens, gap_gens))

        sage_ideal = self.defining_ideal()
        self._gap_rels = libgap.GP2NPList([_sage2gap(x, self._gap_algebra, self._s2g) for x in sage_ideal.gens()])
        self._gap_ideal = GapIdeal(R, I.gens())

    def get_basis(self, maxno=0):
        r"""
        Return a basis for the quotient algebra.

        INPUT:

        - ``maxno`` (default: 0) -- if nonzero, computes a basis until it has
          at least this number of elements

        OUTPUT:

        A basis (or a partial set from the basis, if ``maxno`` is given) for
        the quotient algebra.

        EXAMPLES:

        The commutative algebra with two generators `x, y` such that
        `x^3 = y^3 = 0` is finite dimensional, with basis `\{x^i y^j | 0 \le i, j < 3\}`::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x, x^3, y^3])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.get_basis()    # optional - gap_packages
            [1, x, y, x^2, x*y, y^2, x^2*y, x*y^2, x^2*y^2]

        However, the commutative polynomial algebra with two generators `x, y`
        such that `x^3 = 0` is infinite dimensional. We can compute a partial
        basis::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x, x^3])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.get_basis(10)    # optional - gap_packages
            [1, x, y, x^2, x*y, y^2, x^2*y, x*y^2, y^3, x^2*y^2, x*y^3, y^4]
        """
        res = libgap.BaseQA(self._gap_rels, self.cover_ring().ngens(), maxno)
        return [_gap2sage(x, self.cover_ring()) for x in res]

    def dim(self):
        """
        Return the dimension of the quotient algebra.

        EXAMPLES:

        The commutative algebra with two generators `x, y` such that `x^3 = y^3 = 0`
        has dimension 9, and indeed::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x, x^3, y^3])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.dim()    # optional - gap_packages
            9

        However, the commutative polynomial algebra with two generators is infinite dimensional::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.dim()    # optional - gap_packages
            +Infinity
        """
        if self.is_finite_dimensional():
            return Integer(libgap.DimQA(self._gap_rels, self.cover_ring().ngens()))
        return infinity

    def get_matrix(self, gen_index, basis):
        """
        Given a basis for the quotient algebra, returns the matrix for the
        multiplication of a specific generator.

        INPUT:

        - ``gen_index`` -- the index of the generator

        - ``basis`` -- a basis for the quotient algebra

        OUTPUT:

        The matrix which corresponds to multiplication by the generator at
        index ``gen_index``, with respect to the basis ``basis``.  Each row
        correpsonds to the multiplication of the generator with a basis
        element.

        EXAMPLES:

        The commutative algebra with two generators `x, y` such that 
        `x^3 = y^3 = 0` has the basis

        .. MATH::

            [1, x, y, x^2, x y, y^2, x^2 y, x y^2, x^2 y^2].

        Multiplying the basis by `x`, we get the elements

        .. MATH::

            [x, x^2, x y, 0, x^2 y, x y^2, 0, x^2 y^2, 0].

        Indeed::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x, x^3, y^3])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: B = QA.get_basis()    # optional - gap_packages
            sage: QA.get_matrix(0, B)    # optional - gap_packages
            [0 1 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 1 0]
            [0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 1]
            [0 0 0 0 0 0 0 0 0]
        """
        gap_basis = libgap.GP2NPList([_sage2gap(x, self._gap_algebra, self._s2g) for x in basis])
        return Matrix(self.base_ring(), libgap.MatrixQA(gen_index+1, gap_basis, self._gap_rels))

    def get_matrices(self, basis):
        """
        Return the matrices corresponding to multiplication by each generator
        of the algebra.

        INPUT:

        - ``basis`` -- a basis for the quotient algebra

        OUTPUT:

        The list of matrices which correspond to multiplication of the basis
        ``basis`` by each generator with respect to the basis ``basis``.  Each
        row correpsonds to the multiplication of the generator with a basis
        element.

        EXAMPLES:

        The commutative algebra with two generators `x, y` such that
        `x^3 = y^3 = 0` has the basis

        .. MATH::

            [1, x, y, x^2, x y, y^2, x^2 y, x y^2, x^2 y^2].

        Multiplying the basis by `x`, we get the elements

        .. MATH::

            [x, x^2, x y, 0, x^2 y, x y^2, 0, x^2 y^2, 0].

        Multiplying the basis by y, we get the elements

        .. MATH::

            [y, x y, y^2, x^2 y, x y^2, 0, x^2 y^2, 0, 0].

        Indeed::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x, x^3, y^3])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: B = QA.get_basis()    # optional - gap_packages
            sage: QA.get_matrices(B)    # optional - gap_packages
            [
            [0 1 0 0 0 0 0 0 0]  [0 0 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0]  [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]  [0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 1 0 0]  [0 0 0 0 0 0 0 1 0]
            [0 0 0 0 0 0 0 1 0]  [0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 1]
            [0 0 0 0 0 0 0 0 1]  [0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0], [0 0 0 0 0 0 0 0 0]
            ]
        """
        gap_basis = libgap.GP2NPList([_sage2gap(x, self._gap_algebra, self._s2g) for x in basis])
        return [Matrix(self.base_ring(), M) for M in libgap.MatricesQA(self.ngens(), gap_basis, self._gap_rels)]

    def reduce(self, elem, check=True):
        """
        Given an element of the free algebra, reduces it to a normal form using
        the generators of the defining ideal.

        INPUT:

        - ``elem`` -- an element of the free algebra

        - ``check`` (default: ``True``) -- if ``True``, checks whether the
          generators form a Groebner basis, and if not compute a new Groebner
          basis (might not terminate!)

        OUTPUT:

        The reduced form of ``elem`` with respect to the given generators of
        the ideal.  Note that if the generators do not form a Groebner basis,
        this is not a normal form.

        EXAMPLES:

        We begin with the commutative polynomial ring, in which reducing simply
        means adding the powers of each element::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.reduce(y^3*x*y^2*x^2)    # optional - gap_packages
            x^3*y^5

        Another example, in which the original basis is not a Groebner basis::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x^2*y - x*y, x*y^2 + x^3])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.reduce(x^4 - x^3, check=False)    # optional - gap_packages
            -x^3 + x^4
            sage: QA.reduce(x^4 - x^3, check=True)    # optional - gap_packages
            0
        """
        return self._gap_ideal.reduce(elem, check=check)

    def get_leading_monomials(self, gap_obj=False):
        """
        Return the leading monomials of the generators of the defining ideal.

        INPUT:

        - ``gap_obj`` (default: ``False``) -- if set to ``True``, returns the monomials
          as GBNP objects

        OUTPUT:

        A list of the leading monomials of the generators of the defining ideal.

        EXAMPLES::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x, 2*x^2*y - y^2])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.get_leading_monomials()    # optional - gap_packages
            [y*x, x^2*y]
        """
        if len(self._gap_rels) == 0:
            return []
        lms = libgap.LMonsNP(self._gap_rels)
        if gap_obj:
            return lms
        return [_gap2sage([[x], [libgap(1)]], self.cover_ring()) for x in lms]

    def growth(self, exact=True):
        """
        Determines the growth of the quotient algebra.

        INPUT:

        - ``exact`` (default: ``True``) -- if set to ``False`` and the growth is
          polynomial, returns a list of the possible degrees

        OUTPUT:

        - If the growth is polynomial and exact is ``True``, the degree of the growth.
        - If the growth is polynomial and exact is ``False``, a list of possible degrees.
        - If the growth is exponential, the string "exponential growth".

        Note that if the generators of the defining ideal do not form a Groebner basis, the result may be wrong.

        EXAMPLES:

        Any commutative algebra has polynomial growth. For example, the commutative polynomial algebra with two generators has quadratic growth::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.growth()    # optional - gap_packages
            2

        Also the Weyl algebra has quadratic growth::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x - 1])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.growth()    # optional - gap_packages
            2

        However, taking the free algebra modulo the relation x^2==0 has exponential growth::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x^2])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.growth()    # optional - gap_packages
            'exponential growth'

        With parameter ``exact`` set to ``False``, the result may be a list::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y^2, y*x^2])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.growth(False)    # optional - gap_packages
            [2, 3]
            sage: QA.growth(True)    # optional - gap_packages
            2

        If the generators of the ideal do not form a Groebner basis, the result may be wrong, as the following example shows::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x - y^2])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.growth()    # optional - gap_packages
            'exponential growth'
            sage: I2 = A.ideal([x*y - y*x - x^2])    # optional - gap_packages
            sage: QA2 = A.quo(I2)    # optional - gap_packages
            sage: QA2.growth()    # optional - gap_packages
            2
        """
        res = libgap.DetermineGrowthQA(self.get_leading_monomials(gap_obj=True), self.ngens(), exact)
        if libgap.IsString(res):
            return str(res)
        elif libgap.IsList(res):
            return [Integer(c) for c in res]
        return Integer(res)

    def is_finite_dimensional(self):
        """
        Return ``True`` if the quotient algebra is finite dimensional, otherwise
        returns ``False``.

        EXAMPLES:

        The commutative algebra with two generators `x, y` such that
        `x^3 = y^3 = 0` is finite dimensional::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x, x^3, y^3])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.is_finite_dimensional()    # optional - gap_packages
            True

        However, the noncommutative algebra with two generators `x, y` such
        that `x^3 = y^3 = 0` is not finite dimensional::

            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x^3, y^3])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.is_finite_dimensional()    # optional - gap_packages
            False
        """
        return bool(libgap.FinCheckQA(self.get_leading_monomials(gap_obj=True), self.ngens()))

    def hilbert_series(self, deg):
        r"""
        Return the first ``deg`` coefficients of the Hilbert series for the given
        quotient algebra.

        INPUT:

        - ``deg`` -- The maximal degree to compute in the Hilbert series

        OUTPUT:

        The first ``deg`` coefficients of the Hilbert series for the given
        quotient algebra.

        EXAMPLES:

        The commutative polynomial algebra with two generators has Hilbert
        series `1+2t+3t^2+4t^3+\cdots`::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA.hilbert_series(10)    # optional - gap_packages
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        """
        return [Integer(c) for c in libgap.HilbertSeriesQA(self.get_leading_monomials(gap_obj=True), self.ngens(), deg)]


class GapFreeAlgebra(FreeAlgebra_generic):
    def __init__(self, R, n=None, names=None):
        """
        A non-commutative free algebra implemented in GAP. The current implementation only supports the Rationals field and finite fields.

        INPUT:

        - ``R`` -- the base ring of the algebra. Should be a field

        - ``n`` (default: ``None``) -- the number of generators of the algebra.
          If None, ``names`` must be given.

        - ``names`` (default: ``None``) -- the names for the generators of the
          algebra. If None, ``n`` must be given.

        EXAMPLES:

        One can initialize the free algebra in any of the following ways::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x0, x1> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: A    # optional - gap_packages
            Free Algebra on 2 generators (x0, x1) over Rational Field
            sage: A = GapFreeAlgebra(QQ, 2)    # optional - gap_packages
            sage: A    # optional - gap_packages
            Free Algebra on 2 generators (x0, x1) over Rational Field
            sage: A = GapFreeAlgebra(QQ, 2, 'x')    # optional - gap_packages
            sage: A    # optional - gap_packages
            Free Algebra on 2 generators (x0, x1) over Rational Field
        """
        GapPackage("GBNP", spkg="gap_packages").require()
        libgap.load_package("GBNP")
        libgap.SetInfoLevel(libgap.InfoGBNP, 0)
        libgap.SetInfoLevel(libgap.InfoGBNPTime, 0)

        assert n is not None or names is not None, "A GBNP free algebra must be provided with either number of generators or with names for the generators"
        if n is None:
            n = len(names)
        elif names is None:
            names = ['x{}'.format(i) for i in range(n)]
        FreeAlgebra_generic.__init__(self, R, n, names)
        self._gap_algebra = libgap.FreeAssociativeAlgebraWithOne(libgap(self.base_ring()),
                                        self.variable_names())

    def ideal(self, *gens, **kwds):
        """
        Return the ideal generated by the elements in ``gens``.

        INPUT:

        - ``gens`` -- list or tuple of generators (or several input arguments)

        OUTPUT:

        The ideal generated by ``gens``, implemented in GAP.

        EXAMPLES::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x, x^2])    # optional - gap_packages
            sage: I    # optional - gap_packages
            Twosided Ideal (x*y - y*x, x^2) of Free Algebra on 2 generators (x, y) over Rational Field
        """
        I = super(FreeAlgebra_generic, self).ideal(*gens, **kwds)
        return GapIdeal(self, I.gens())

    def quotient(self, rels):
        """
        Return a quotient algebra.

        INPUT:

        - ``rels`` -- a list of the generators for the defining ideal of the quotient

        OUTPUT:

        The quotient algebra of the free algebra modulo ``rels``, implemented in GAP.

        EXAMPLES::

            sage: from sage.algebras.gbnp import GapFreeAlgebra
            sage: A.<x, y> = GapFreeAlgebra(QQ)    # optional - gap_packages
            sage: I = A.ideal([x*y - y*x, x^2])    # optional - gap_packages
            sage: QA = A.quo(I)    # optional - gap_packages
            sage: QA    # optional - gap_packages
            Quotient of Free Algebra on 2 generators (x, y) over Rational Field by the ideal (x*y - y*x, x^2)
        """
        if type(rels) == GapIdeal:
            return GapQuotientRing(self, rels)
        return GapQuotientRing(self, GapIdeal(self, rels))

    quo = quotient
