# -*- coding: utf-8 -*-
"""
`p`-adic valuations on number fields and their subrings and completions.

AUTHORS:

- Julian Rüth (2013-03-16): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013-2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from valuation import DiscreteValuation
from value_group import DiscreteValueSemigroup
from mapped_valuation import FiniteExtensionFromLimitValuation
from sage.structure.factory import UniqueFactory
from sage.misc.cachefunc import cached_method
from sage.misc.fast_methods import WithEqualityById

from sage.rings.all import infinity

class PadicValuationFactory(UniqueFactory):
    """
    Create a ``prime``-adic valuation on ``R``.

    INPUT:

    - ``R`` -- a subring of a number field or a subring of a local field in
      characteristic zero.

    - ``prime`` -- a prime that does not split, a discrete (pseudo-)valuation,
      a fractional ideal, or ``None`` (default: ``None``)

    EXAMPLES:

    For integers and rational numbers, ``prime`` is just a prime of the
    integers::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: pAdicValuation(ZZ, 3)
        3-adic valuation

        sage: pAdicValuation(QQ, 3)
        3-adic valuation

    ``prime`` may be ``None`` for local rings::

        sage: pAdicValuation(Qp(2))
        2-adic valuation

        sage: pAdicValuation(Zp(2))
        2-adic valuation

    But it must be specified in all other cases::

        sage: pAdicValuation(ZZ)
        Traceback (most recent call last):
        ...
        ValueError: prime must be specified for this ring

    For number fields, ``prime`` can be an integer that is completely ramified
    in ``R``::

        sage: pAdicValuation(GaussianIntegers().fraction_field(), 2)
        2-adic valuation

    For number fields, ``prime`` can be an integer that is unramified in ``R``:

        sage: pAdicValuation(GaussianIntegers().fraction_field(), 3)
        3-adic valuation

    The same applies if ``R`` is a subring of a number field::
    
        sage: pAdicValuation(GaussianIntegers(), 3)
        3-adic valuation

    However, this is only supported if ``prime`` does not factor into
    pairwise distinct factors::

        sage: pAdicValuation(GaussianIntegers(), 5)
        Traceback (most recent call last):
        ...
        ValueError: The valuation Gauss valuation induced by 5-adic valuation does not approximate a unique extension of 5-adic valuation with respect to x^2 + 1

    When ``R`` is an absolute or relative number field, or a subring thereof,
    ``prime`` can also be specified by providing a valuation on the base ring
    that has a unique extension::

        sage: pAdicValuation(CyclotomicField(5), pAdicValuation(ZZ, 5))
        5-adic valuation

    When the extension is not unique, this does not work::

        sage: pAdicValuation(GaussianIntegers(), pAdicValuation(ZZ, 5))
        Traceback (most recent call last):
        ...
        ValueError: The valuation Gauss valuation induced by 5-adic valuation does not approximate a unique extension of 5-adic valuation with respect to x^2 + 1

    For a number field which is of the form `K[x]/(G)`, you can specify a
    valuation by providing a discrete pseudo-valuation on `K[x]` which sends
    `G` to `\infty`. This lets us specify which extension of the 5-adic
    valuation we care about in the above example::

        sage: R.<x> = QQ[]
        sage: v = pAdicValuation(GaussianIntegers(), GaussValuation(R, pAdicValuation(QQ, 5)).augmentation(x + 2, infinity))
        sage: w = pAdicValuation(GaussianIntegers(), GaussValuation(R, pAdicValuation(QQ, 5)).augmentation(x + 1/2, infinity))
        sage: v == w
        False

    Note that you get the same valuation, even if you write down the
    pseudo-valuation differently::

        sage: ww = pAdicValuation(GaussianIntegers(), GaussValuation(R, pAdicValuation(QQ, 5)).augmentation(x + 3, infinity))
        sage: w is ww
        True

    The valuation ``prime`` does not need to send the defining polynomial `G`
    to `\infty`. It is sufficient if it singles out one of the valuations on
    the number field.  This is important if the prime only factors over the
    completion, i.e., if it is not possible to write down one of the factors
    within the number field::

        sage: v = GaussValuation(R, pAdicValuation(QQ, 5)).augmentation(x + 3, 1)
        sage: pAdicValuation(GaussianIntegers().fraction_field(), v)
        [ 5-adic valuation, v(x + 3) = 1 ]-adic valuation

    Finally, ``prime`` can also be a fractional ideal of a number field if it
    singles out an extension of a `p`-adic valuation of the base field::

        sage: R = GaussianIntegers()
        sage: I = R.fraction_field().gen()
        sage: pAdicValuation(R, R.fractional_ideal(I + 1))
        2-adic valuation

    It can sometimes be beneficial to define a number field extension as a
    quotient of a polynomial ring (since number field extensions always compute
    an absolute polynomial defining the extension which can be very costly)::

        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^2 + 1)
        sage: R.<x> = K[]
        sage: L.<b> = R.quo(x^2 + a)
        sage: pAdicValuation(L, 2)
        2-adic valuation

    """
    def create_key_and_extra_args(self, R, prime=None, approximants=None):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(QQ, 2) # indirect doctest
            2-adic valuation

        """
        from sage.rings.all import ZZ, QQ
        from sage.rings.padics.padic_generic import pAdicGeneric
        from sage.rings.number_field.number_field import is_NumberField
        from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing

        if R.characteristic() != 0:
            # We do not support equal characteristic yet
            raise ValueError("R must be a ring of characteristic zero.")

        if R is ZZ or R is QQ:
            return self.create_key_for_integers(R, prime), {}
        elif isinstance(R, pAdicGeneric):
            return self.create_key_for_local_ring(R, prime), {}
        elif is_NumberField(R.fraction_field()) or is_PolynomialQuotientRing(R):
            return self.create_key_and_extra_args_for_number_field(R, prime, approximants=approximants)
        else:
            raise NotImplementedError("p-adic valuations not implemented for %r"%(R,))

    def create_key_for_integers(self, R, prime):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(QQ, 2) # indirect doctest
            2-adic valuation

        """
        from sage.rings.all import ZZ
        if prime is None:
            raise ValueError("prime must be specified for this ring")
        from valuation import DiscretePseudoValuation
        if isinstance(prime, DiscretePseudoValuation):
            prime = prime.uniformizer()
        if prime not in ZZ or not ZZ(prime).is_prime():
            raise ValueError("prime must be a prime in the integers but %s is not"%(prime,))
        return R, prime

    def create_key_for_local_ring(self, R, prime):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(Qp(2)) # indirect doctest
            2-adic valuation

        """
        # We do not care much about the value of prime since there is only one
        # reasonable p-adic valuation here
        if prime is not None:
            if prime in R:
                if R(prime).valuation() <= 0:
                    raise ValueError("prime must be an element of positive valuation")
            elif prime(R.prime()) <= 0:
                raise ValueError("prime must be an element of positive valuation")

        return (R,)

    def create_key_and_extra_args_for_number_field(self, R, prime, approximants):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(GaussianIntegers(), 2) # indirect doctest
            2-adic valuation

        """
        K, L, G = self._normalize_number_field_data(R)

        from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
        from valuation import DiscretePseudoValuation
        if isinstance(prime, DiscretePseudoValuation):
            return self.create_key_and_extra_args_for_number_field_from_valuation(R, prime, prime, approximants=approximants)
        elif prime in K:
            return self.create_key_and_extra_args_for_number_field_from_valuation(R, pAdicValuation(K, prime), prime, approximants=approximants)
        elif prime in L or isinstance(prime, NumberFieldFractionalIdeal):
            return self.create_key_and_extra_args_for_number_field_from_ideal(R, L.fractional_ideal(prime), prime)
        else:
            raise ValueError("prime must be a discrete pseudo-valuation, a prime in the base ring, or a fractional ideal")

    def create_key_and_extra_args_for_number_field_from_valuation(self, R, v, prime, approximants):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``v``.

        .. NOTE::

            ``prime``, the original parameter that was passed to
            :meth:`create_key_and_extra_args``, is only used to provide more
            meaningful error messages

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(GaussianIntegers(), pAdicValuation(ZZ, 2)) # indirect doctest
            2-adic valuation

        TESTS:

        We can extend to the field of fractions of a quotient ring::

            sage: R.<x> = ZZ[]
            sage: S = R.quo(x^2 + 1)
            sage: v = pAdicValuation(S, 2)
            sage: R.<x> = QQ[]
            sage: S = R.quo(x^2 + 1)
            sage: v = pAdicValuation(S, v)

        """
        K, L, G = self._normalize_number_field_data(R)

        if v.domain().is_subring(G.parent()):
            # v is defined on a subring of K[x].
            # We try to lift v to a pseudo-valuation on K[x].
            if _fraction_field(v.domain()) is not _fraction_field(G.parent()):
                # First, we lift valuations defined on subrings of K to
                # valuations on K[x].
                if v.domain().is_subring(K):
                    if v.domain() is not K:
                        v = pAdicValuation(K, v)
                    from gauss_valuation import GaussValuation
                    v = GaussValuation(G.parent(), v)
            if v.domain() != G.parent():
                # Then, we lift valuations defined on polynmial rings which are
                # subrings of K[x] to K[x]
                v = v.extension(G.parent())
        elif _fraction_field(v.domain()) == L:
            # v is defined on a ring whose field of fractions is L
            v = v._base_valuation._initial_approximation.change_domain(G.parent())
        else:
            raise NotImplementedError("can not rewrite %r which is defined on %r as a pseudo-valuation on %r"%(v, v.domain(), G.parent()))
            

        assert(v.domain() is G.parent())

        # To obtain uniqueness of p-adic valuations, we need a canonical
        # description of v. We consider all extensions of vK to L and select
        # the one approximated by v.
        vK = v.restriction(v.domain().base_ring()).extension(K)
        if approximants is None:
            approximants = vK.mac_lane_approximants(G)
        approximants = [approximant.extension(v.domain()) for approximant in approximants]
        approximant = vK.mac_lane_approximant(G, v, approximants=tuple(approximants))

        return (R, approximant, L.construction()), {'approximants': approximants}

    def create_key_and_extra_args_for_number_field_from_ideal(self, R, I, prime):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``I``.

        .. NOTE::

            ``prime``, the original parameter that was passed to
            :meth:`create_key_and_extra_args``, is only used to provide more
            meaningful error messages

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(GaussianIntegers(), GaussianIntegers().ideal(2)) # indirect doctest
            2-adic valuation

        """
        K, L, G = self._normalize_number_field_data(R)

        # To obtain uniqueness of p-adic valuations, we need a canonical
        # description of v. We consider all extensions of vK to L and select
        # the one approximated by v.
        # Of course, this only works if I comes from a single prime downstairs.
        p = I.relative_norm()
        F = p.factor()
        if len(F) != 1:
            raise ValueError("%r does not lie over a single prime of %r"%(I, K))
        vK = pAdicValuation(K, F[0][0])
        candidates = vK.mac_lane_approximants(G)

        candidates_for_I = [c for c in candidates if all(c(g.polynomial()) > 0 for g in I.gens())]
        assert(len(candidates_for_I) > 0) # This should not be possible, unless I contains a unit
        if len(candidates_for_I) > 1:
            raise ValueError("%s does not single out a unique extension of %s to %s"%(prime, vK, L))
        else:
            # equality of number fields has it quirks since it says that two
            # fields are == even if they are distinguishable (because they come
            # from different constructions.)
            # Including structure() into the key seems to be a way to distinguish such cases properly.
            # This used to be an issue but seems to be fixed, namely, the
            # absolute_field of a number field was deemed equivalent to the
            # directly created absolute field, even though the absolute_field
            # carried the information where it came from
            return (R, candidates_for_I[0], L.construction()), {'approximants': candidates}

    def _normalize_number_field_data(self, R):
        r"""
        Helper method which returns the defining data of the number field
        ``R``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: K = R.quo(x^2 + 1)
            sage: pAdicValuation._normalize_number_field_data(K)
            (Rational Field,
             Univariate Quotient Polynomial Ring in xbar over Rational Field with modulus x^2 + 1,
             x^2 + 1)

        """
        from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing
        from sage.rings.number_field.number_field import is_NumberField
        from sage.rings.fraction_field import is_FractionField
        if is_NumberField(R.fraction_field()):
            L = R.fraction_field()
            G = L.relative_polynomial()
            K = L.base_ring()
        elif is_PolynomialQuotientRing(R):
            from sage.categories.all import NumberFields
            if R.base_ring().fraction_field() not in NumberFields():
                raise NotImplementedError("can not normalize quotients over %r"%(R.base_ring(),))
            L = R.fraction_field()
            K = R.base_ring().fraction_field()
            G = R.modulus().change_ring(K)
        else:
            raise NotImplementedError("can not normalize %r"%(R,))

        return K, L, G


    def create_object(self, version, key, **extra_args):
        r"""
        Create a `p`-adic valuation from ``key``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(ZZ, 5) # indirect doctest
            5-adic valuation

        """
        from sage.rings.all import ZZ, QQ
        from sage.rings.padics.padic_generic import pAdicGeneric
        from valuation_space import DiscretePseudoValuationSpace
        from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing
        from sage.rings.number_field.number_field import is_NumberField
        R = key[0]
        K = R.fraction_field()
        parent = DiscretePseudoValuationSpace(R)
        if isinstance(R, pAdicGeneric):
            assert(len(key)==1)
            return parent.__make_element_class__(pAdicValuation_padic)(parent)
        elif R is ZZ or R is QQ:
            prime = key[1]
            assert(len(key)==2)
            return parent.__make_element_class__(pAdicValuation_int)(parent, prime)
        else:
            v = key[1]
            _ = key[2] # ignored
            approximants = extra_args['approximants']
            parent = DiscretePseudoValuationSpace(R)
            if is_NumberField(K):
                G = K.relative_polynomial()
            elif is_PolynomialQuotientRing(R):
                G = R.modulus()
            else:
                raise NotImplementedError
            return parent.__make_element_class__(pAdicFromLimitValuation)(parent, v, G.change_ring(R.base_ring()), approximants)

pAdicValuation = PadicValuationFactory("pAdicValuation")

class pAdicValuation_base(DiscreteValuation):
    """
    Abstract base class for `p`-adic valuations.

    INPUT:

    - ``ring`` -- an integral domain

    - ``p`` -- a rational prime over which this valuation lies, not
      necessarily a uniformizer for the valuation

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: pAdicValuation(ZZ, 3)
        3-adic valuation

        sage: pAdicValuation(QQ, 5)
        5-adic valuation

     For `p`-adic rings, ``p`` has to match the `p` of the ring.

        sage: v = pAdicValuation(Zp(3), 2); v
        Traceback (most recent call last):
        ...
        ValueError: prime must be an element of positive valuation

    TESTS::

        sage: TestSuite(pAdicValuation(ZZ, 3)).run() # long time
        sage: TestSuite(pAdicValuation(QQ, 5)).run() # long time
        sage: TestSuite(pAdicValuation(Zp(5), 5)).run() # long time

    """
    def __init__(self, parent, p):
        """
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: isinstance(pAdicValuation(ZZ, 3), pAdicValuation_base)
            True

        """
        DiscreteValuation.__init__(self, parent)

        from sage.rings.all import ZZ
        self._p = ZZ(p)

    def p(self):
        """
        Return the `p` of this `p`-adic valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(GaussianIntegers(), 2).p()
            2

        """
        return self._p

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        OUTPUT:

        An element of the :meth:`residue_field`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 3)
            sage: v.reduce(4)
            1

        """
        x = self.domain().coerce(x)

        if self(x) < 0:
            raise ValueError("reduction is only defined for elements of non-negative valuation")

        return self.residue_field()(x)

    def lift(self, x):
        """
        Lift ``x`` from the residue field to the domain of this valuation.

        INPUT:

        - ``x`` -- an element of the :meth:`residue_field`

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 3)
            sage: xbar = v.reduce(4)
            sage: v.lift(xbar)
            1

        """
        x = self.residue_field().coerce(x)

        return self.domain()(x)

    def is_unramified(self, G, include_steps=False, assume_squarefree=False):
        """
        Return whether ``G`` defines a single unramified extension of the
        completion of the domain of this valuation.

        INPUT:

        - ``G`` -- a monic squarefree polynomial over the domain of this valuation

        - ``include_steps`` -- a boolean (default: ``False``); whether to
          include the approximate valuations that were used to determine the
          result in the return value.

        - ``assume_squarefree`` -- a boolean (default: ``False``); whether to
          assume that ``G`` is square-free over the completion of the domain of
          this valuation. Setting this to ``True`` can significantly improve
          the performance.

        EXAMPLES:

        We consider an extension as unramified if its ramification index is 1.
        Hence, a trivial extension is unramified::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: v = pAdicValuation(QQ, 2)
            sage: v.is_unramified(x)
            True

        If ``G`` remains irreducible in reduction, then it defines an
        unramified extension::

            sage: v.is_unramified(x^2 + x + 1)
            True

        However, even if ``G`` factors, it might define an unramified
        extension::

            sage: v.is_unramified(x^2 + 2*x + 4)
            True

        """
        R = G.parent()

        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if not is_PolynomialRing(R) or R.base_ring() is not self.domain() or not len(R.gens()) == 1 or not G.is_monic():
            raise ValueError("G must be a monic univariate polynomial over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from gauss_valuation import GaussValuation

        steps = [ GaussValuation(R, self) ]
        while True:
            v = steps[-1]
            if v.E() > 1:
                ret = False
                break
            if v.F() == G.degree():
                ret = True
                break

            assert v(G) is not infinity
            if v.is_key(G):
                ret = True
                break

            next = v.mac_lane_step(G, assume_squarefree=True)
            if len(next)>1:
                ret = False
                break
            steps.append(next[0])

        if include_steps:
            return ret, steps
        else:
            return ret

    def is_totally_ramified(self, G, include_steps=False, assume_squarefree=False):
        """
        Return whether ``G`` defines a single totally ramified extension of the
        completion of the domain of this valuation.

        INPUT:

        - ``G`` -- a monic squarefree polynomial over the domain of this valuation

        - ``include_steps`` -- a boolean (default: ``False``); where to include
          the valuations produced during the process of checking whether ``G``
          is totally ramified in the return value

        - ``assume_squarefree`` -- a boolean (default: ``False``); whether to
          assume that ``G`` is square-free over the completion of the domain of
          this valuation. Setting this to ``True`` can significantly improve
          the performance.

        ALGORITHM:

        This is a simplified version of :meth:`mac_lane_approximants`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: k=Qp(5,4)
            sage: v = pAdicValuation(k)
            sage: R.<x>=k[]
            sage: G = x^2 + 1
            sage: v.is_totally_ramified(G)
            False
            sage: G = x + 1
            sage: v.is_totally_ramified(G)
            True
            sage: G = x^2 + 2
            sage: v.is_totally_ramified(G)
            False
            sage: G = x^2 + 5
            sage: v.is_totally_ramified(G)
            True
            sage: v.is_totally_ramified(G, include_steps=True)
            (True, [Gauss valuation induced by 5-adic valuation, [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x) = 1/2 ]])

        """
        R = G.parent()

        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if not is_PolynomialRing(R) or R.base_ring() is not self.domain() or not len(R.gens()) == 1 or not G.is_monic():
            raise ValueError("G must be a monic univariate polynomial over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from gauss_valuation import GaussValuation

        steps = [ GaussValuation(R, self) ]
        while True:
            v = steps[-1]
            if v.F() > 1:
                ret = False
                break
            if v.E() == G.degree():
                ret = True
                break

            assert v(G) is not infinity
            if v.is_key(G):
                ret = False
                break

            next = v.mac_lane_step(G, assume_squarefree=True)
            if len(next)>1:
                ret = False
                break
            steps.append(next[0])

        if include_steps:
            return ret, steps
        else:
            return ret

    def change_domain(self, ring):
        r"""
        Change the domain of this valuation to ``ring`` if possible.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: v.change_domain(QQ).domain()
            Rational Field

        """
        return pAdicValuation(ring, self.p())

    def _extensions_to_quotient(self, ring, approximants=None):
        r"""
        Return the extensions of this valuation to an integral quotient over
        the domain of this valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R.<x> = QQ[]
            sage: pAdicValuation(QQ, 2)._extensions_to_quotient(R.quo(x^2 + x + 1))
            [2-adic valuation]

        """
        from valuation_space import DiscretePseudoValuationSpace
        parent = DiscretePseudoValuationSpace(ring)
        approximants = approximants or self.mac_lane_approximants(ring.modulus().change_ring(self.domain()), assume_squarefree=True)
        return [pAdicValuation(ring, approximant, approximants) for approximant in approximants]

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: v.extensions(GaussianIntegers())
            [2-adic valuation]

        TESTS::

            sage: R.<a> = QQ[]
            sage: L.<a> = QQ.extension(x^3 - 2)
            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b^2 + 2*b + a)
            sage: pAdicValuation(M, 2)
            2-adic valuation

        Check that we can extend to a field written as a quotient::

            sage: R.<x> = QQ[]
            sage: K.<a> = QQ.extension(x^2 + 1)
            sage: R.<y> = K[]
            sage: L.<b> = R.quo(x^2 + a)
            sage: pAdicValuation(QQ, 2).extensions(L)
            [2-adic valuation]

        """
        if self.domain() is ring:
            return [self]
        domain_fraction_field = _fraction_field(self.domain())
        if domain_fraction_field is not self.domain():
            if domain_fraction_field.is_subring(ring):
                return pAdicValuation(domain_fraction_field, self).extensions(ring)
        if self.domain().is_subring(ring):
            from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing
            if is_PolynomialQuotientRing(ring):
                if is_PolynomialQuotientRing(self.domain()):
                    if self.domain().modulus() == ring.modulus():
                        base_extensions = self._base_valuation.extensions(self._base_valuation.domain().change_ring(self._base_valuation.domain().base_ring().fraction_field()))
                        return [pAdicValuation(ring, base._initial_approximation) for base in base_extensions]
                if ring.base_ring() is self.domain():
                    from sage.categories.all import IntegralDomains
                    if ring in IntegralDomains():
                        return self._extensions_to_quotient(ring)
                else:
                    return sum([w.extensions(ring) for w in self.extensions(ring.base_ring())], [])
            from sage.rings.number_field.number_field import is_NumberField
            if is_NumberField(ring.fraction_field()):
                if ring.base_ring().fraction_field() is self.domain().fraction_field():
                    from valuation_space import DiscretePseudoValuationSpace
                    parent = DiscretePseudoValuationSpace(ring)
                    approximants = self.mac_lane_approximants(ring.fraction_field().relative_polynomial().change_ring(self.domain()), assume_squarefree=True)
                    return [pAdicValuation(ring, approximant, approximants) for approximant in approximants]
                if ring.base_ring() is not ring and self.domain().is_subring(ring.base_ring()):
                    return sum([w.extensions(ring) for w in self.extensions(ring.base_ring())], [])
        return super(pAdicValuation_base, self).extensions(ring)

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(GaussianIntegers(), 2)
            sage: v.restriction(ZZ)
            2-adic valuation

        """
        if ring is self.domain():
            return self

        if not ring.is_subring(self.domain()):
            raise ValueError("ring must be a subring of the domain of this valuation but %r is not a subring of %r"%(ring, self.domain()))

        return pAdicValuation(ring, self.p())

    @cached_method
    def value_semigroup(self):
        r"""
        Return the value semigroup of this valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(GaussianIntegers(), 2)
            sage: v.value_semigroup()
            Additive Abelian Semigroup generated by 1/2

        """
        from sage.categories.all import Fields
        v = self(self.uniformizer())
        if self.domain() in Fields():
            return DiscreteValueSemigroup([-v,v])
        else:
            return DiscreteValueSemigroup([v])


class pAdicValuation_padic(pAdicValuation_base):
    """
    The `p`-adic valuation of a complete `p`-adic ring.

    INPUT:

    - ``R`` -- a `p`-adic ring

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: v = pAdicValuation(Qp(2)); v #indirect doctest
        2-adic valuation

    TESTS::

        sage: TestSuite(v).run() # optional: integrated, long time

    """
    def __init__(self, parent):
        """
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: from sage.rings.padics.padic_valuation import padicValuation_padic # optional: integrated
            sage: isinstance(pAdicValuation(Qp(2)), pAdicValuation_padic)
            True

        """
        pAdicValuation_base.__init__(self, parent, parent.domain().prime())

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element of the domain of this valuation

        OUTPUT:

        An element of the :meth:`residue_field`.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R = Zp(3)
            sage: pAdicValuation(Zp(3)).reduce(R(4))
            1

        """
        x = self.domain().coerce(x)
        return self.residue_field()(x.residue())

    def lift(self, x):
        """
        Lift ``x`` from the :meth:`residue_field` to the :meth:`domain` of this
        valuation.

        INPUT:

        - ``x`` -- an element of the :meth:`residue_field`

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R = Zp(3)
            sage: v = pAdicValuation(R)
            sage: xbar = v.reduce(R(4))
            sage: v.lift(xbar)
            1 + O(3^20)

        """
        x = self.residue_field().coerce(x)
        return self.domain()(x).lift_to_precision()

    def uniformizer(self):
        """
        Return a uniformizer of this valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(Zp(3))
            sage: v.uniformizer()
            3 + O(3^21)

        """
        return self.domain().uniformizer()

    def element_with_valuation(self, v):
        """
        Return an element of valuation ``v``.

        INPUT:

        - ``v`` -- an element of the :meth:`value_semigroup` of this valuation

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R = Zp(3)
            sage: v = pAdicValuation(Zp(3))
            sage: v.element_with_valuation(3)
            3^3 + O(3^23)

        """
        from sage.rings.all import QQ, ZZ
        v = QQ(v)
        if v not in self.value_semigroup():
            raise ValueError("%r is not in the value semigroup of %r"%(v, self))
        v = ZZ(v * self.domain().ramification_index())
        return self.domain().one() << v

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(ZZ, 3)._repr_()
            '3-adic valuation'

        """
        return "%s-adic valuation"%(self.p())

    def _call_(self, x):
        r"""
        Evaluate this valuation at ``x``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: K = Qp(3)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - 3)
            sage: pAdicValuation(L, 3)(3)
            1

        """
        return x.ordp()

    def residue_ring(self):
        r"""
        Return the residue field of this valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(Qq(9, names='a'), 3).residue_ring()
            Finite Field in a0 of size 3^2

        """
        return self.domain().residue_field()

    def shift(self, x, s):
        r"""
        Shift ``x`` in its expansion with respect to :meth:`uniformizer` by
        ``s`` "digits".

        For non-negative ``s``, this just returns ``x`` multiplied by a
        power of the uniformizer `\pi`.

        For negative ``s``, it does the same but when not over a field, it
        drops coefficients in the `\pi`-adic expension which have negative
        valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R = ZpCA(2)
            sage: v = pAdicValuation(R)
            sage: v.shift(R.one(), 1)
            2 + O(2^20)
            sage: v.shift(R.one(), -1)
            O(2^19)

        """
        from sage.rings.all import ZZ
        x = self.domain().coerce(x)
        s = self.value_group()(s)
        v = ZZ(s / self.domain().ramification_index())
        return x << v

    def simplify(self, x, error=None, force=False):
        r"""
        Return a simplified version of ``x``.

        Produce an element which differs from ``x`` by an element of
        valuation strictly greater than the valuation of ``x`` (or strictly
        greater than ``error`` if set.)

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        - ``error`` -- a rational, infinity, or ``None`` (default: ``None``),
          the error allowed to introduce through the simplification

        - ``force`` -- ignored

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: R = Zp(2)
            sage: v = pAdicValuation(R, 2)
            sage: v.simplify(6)
            2 + O(2^21)
            sage: v.simplify(6, error=0)
            0

        """
        x = self.domain().coerce(x)

        if error is None:
            error = self(x)
        from sage.rings.all import infinity
        if error is infinity:
            return x
        normalized_error = (error / self.value_group().gen()).ceil()
        return x.add_bigoh(normalized_error + 1).lift_to_precision()


class pAdicValuation_int(pAdicValuation_base):
    r"""
    A `p`-adic valuation on the integers or the rationals.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: v = pAdicValuation(ZZ, 3); v
        3-adic valuation

    TESTS::

        sage: TestSuite(v).run() # long time

    """
    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(ZZ, 3)._repr_()
            '3-adic valuation'

        """
        return "%s-adic valuation"%(self.p())

    def _call_(self, x):
        """
        Evaluate this valuation at ``x``.

        INPUT::

        - ``x`` --  an element in the domain of this valuation

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: pAdicValuation(ZZ, 3)(9)
            2

        """
        if x.is_zero():
            # x.valuation() is a factor 10 slower when computing the valuation
            # of a rational zero than when computing the valuation of another
            # small rational. Special casing this is a factor 100 faster.
            return infinity
        return x.valuation(self._p)

    def uniformizer(self):
        """
        Return a uniformizer of this `p`-adic valuation, i.e., `p` as an
        element of the domain.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 3)
            sage: v.uniformizer()
            3

        """
        return self.domain()(self.p())

    def residue_ring(self):
        """
        Return the residue field of this valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 3)
            sage: v.residue_ring()
            Finite Field of size 3

        """
        from sage.rings.all import GF
        return GF(self.p())

    def _ge_(self, other):
        r"""
        Return whether this valuation is greater than or equal than ``other``
        everywhere.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: w = TrivialValuation(ZZ)
            sage: v >= w
            True

        """
        if other.is_trivial():
            return other.is_discrete_valuation()
        if isinstance(other, pAdicValuation_int):
            return self.p() == other.p()
        return super(pAdicValuation_base, self)._ge_(other)

    def _relative_size(self, x):
        r"""
        Return an estimate on the coefficient size of ``x``.

        The number returned is an estimate on the factor between the number of
        bits used by ``x`` and the minimal number of bits used by an element
        congruent to ``x``.

        This is used by :meth:`simplify` to decide whether simplification of
        coefficients is going to lead to a significant shrinking of the
        coefficients of ``x``.

        EXAMPLES:: 

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: v._relative_size(2)
            2
            sage: v._relative_size(2**20)
            21

        """
        x = self.domain().coerce(x)
        return x.numerator().nbits() + x.denominator().nbits() - 1

    def simplify(self, x, error=None, force=False, size_heuristic_bound=32):
        r"""
        Return a simplified version of ``x``.

        Produce an element which differs from ``x`` by an element of
        valuation strictly greater than the valuation of ``x`` (or strictly
        greater than ``error`` if set.)

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        - ``error`` -- a rational, infinity, or ``None`` (default: ``None``),
          the error allowed to introduce through the simplification

        - ``force`` -- ignored

        - ``size_heuristic_bound` -- when ``force`` is not set, the expected
          factor by which the ``x`` need to shrink to perform an actual
          simplification (default: 32)

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(ZZ, 2)
            sage: v.simplify(6, force=True)
            2
            sage: v.simplify(6, error=0, force=True)
            0

        """
        if not force and self._relative_size(x) <= size_heuristic_bound:
            return x

        x = self.domain().coerce(x)

        v = self(x)
        if error is None:
            error = v
        from sage.rings.all import infinity
        if error is infinity:
            return x
        if error < v:
            return self.domain().zero()
        from sage.rings.all import QQ
        error = QQ(error).ceil()
        
        from sage.rings.all import Qp
        precision_ring = Qp(self.p(), error + 1 - v)
        reduced = precision_ring(x)
        if error - v >= 5:
            # If there is not much relative precision left, it is better to
            # just go with the integer/rational lift. The rational
            # reconstruction is likely not smaller.
            try:
                reconstruction = reduced.rational_reconstruction()
                if reconstruction in self.domain():
                    return self.domain()(reconstruction)
            except ArithmeticError:pass
        
        return self.domain()(reduced.lift())


class pAdicFromLimitValuation(FiniteExtensionFromLimitValuation, pAdicValuation_base):
    r"""
    A `p`-adic valuation on a number field or a subring thereof, i.e., a
    valuation that extends the `p`-adic valuation on the integers.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
        sage: v = pAdicValuation(GaussianIntegers(), 3); v
        3-adic valuation

    TESTS::

        sage: TestSuite(v).run() # long time

    """
    def __init__(self, parent, approximant, G, approximants):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(GaussianIntegers(), 3)
            sage: isinstance(v, pAdicFromLimitValuation)
            True

        """
        FiniteExtensionFromLimitValuation.__init__(self, parent, approximant, G, approximants)
        pAdicValuation_base.__init__(self, parent, approximant.restriction(approximant.domain().base_ring()).p())

    def _to_base_domain(self, f):
        r"""
        Return ``f``, an element of the domain of this valuation, as an element
        of the domain of the underlying limit valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(GaussianIntegers(), 3)
            sage: I = GaussianIntegers().fraction_field().gen()
            sage: v._to_base_domain(I)
            x

        """
        polynomial = f.polynomial() if hasattr(f,'polynomial') else f.lift()
        return polynomial(self._base_valuation.domain().gen())

    def _from_base_domain(self, f):
        r"""
        Return ``f``, an element of the domain of this valuation, as an element
        of the domain of the underlying limit valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(GaussianIntegers(), 3)
            sage: v._from_base_domain(v._base_valuation.domain().gen())
            I

        """
        return self.domain()(f)

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from mac_lane import * # optional: standalone
            sage: v = pAdicValuation(GaussianIntegers(), 3)
            sage: v.extensions(v.domain().fraction_field())
            [3-adic valuation]

        """
        if ring is self.domain().fraction_field():
            if self.domain() is not self.domain().fraction_field():
                base_ring = self.domain().base_ring()
                base_valuation = self.restriction(base_ring).extension(base_ring.fraction_field())
                G = ring.relative_polynomial()
                approximant = self._base_valuation.change_domain(G.parent())._initial_approximation
                return [pAdicValuation(ring, approximant)]
        return super(pAdicFromLimitValuation, self).extensions(ring)

def _fraction_field(ring):
    r"""
    Return a fraction field of ``ring``.

    EXAMPLES:

    This works around some annoyances with ``ring.fraction_field()``::

        sage: R.<x> = ZZ[]
        sage: S = R.quo(x^2 + 1)
        sage: S.fraction_field()
        Fraction Field of Univariate Quotient Polynomial Ring in xbar over Integer Ring with modulus x^2 + 1

        sage: from mac_lane.padic_valuation import _fraction_field
        sage: _fraction_field(S)
        Univariate Quotient Polynomial Ring in xbar over Rational Field with modulus x^2 + 1

    """
    from sage.categories.all import Fields
    if ring in Fields():
        return ring

    from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing
    if is_PolynomialQuotientRing(ring):
        from sage.categories.all import IntegralDomains
        if ring in IntegralDomains():
            return ring.base().change_ring(ring.base_ring().fraction_field()).quo(ring.modulus())
    return ring.fraction_field()
