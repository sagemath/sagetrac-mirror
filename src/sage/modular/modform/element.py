# -*- coding: utf-8 -*-
"""
Elements of modular forms spaces

Class hierarchy:

- :class:`ModularForm_abstract`

  - :class:`Newform`

  - :class:`ModularFormElement`

    - :class:`ModularFormElement_elliptic_curve`

    - :class:`EisensteinSeries`

"""
# ****************************************************************************
#       Copyright (C) 2004-2008 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import absolute_import, division
from six.moves import range

import sage.modular.hecke.element as element

from sage.arith.all import lcm, divisors, moebius, sigma, factor, crt
from sage.arith.srange import xsrange
from sage.matrix.constructor import matrix
from sage.misc.all import prod
from sage.misc.cachefunc import cached_method
from sage.misc.misc import verbose
from sage.misc.superseded import deprecated_function_alias
from sage.modular.dirichlet import DirichletGroup
from sage.modular.modsym.modsym import ModularSymbols
from sage.modular.modsym.p1list import lift_to_sl2z
from sage.modular.modsym.space import is_ModularSymbolsSpace
from sage.modules.free_module_element import vector
from sage.rings.all import ZZ, QQ, Integer, RealField, ComplexField
from sage.rings.fast_arith import prime_range
from sage.rings.morphism import RingHomomorphism
from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
from sage.structure.element import coercion_model, ModuleElement


def is_ModularFormElement(x):
    """
    Return ``True`` if x is a modular form.

    EXAMPLES::

        sage: from sage.modular.modform.element import is_ModularFormElement
        sage: is_ModularFormElement(5)
        False
        sage: is_ModularFormElement(ModularForms(11).0)
        True
    """
    return isinstance(x, ModularFormElement)


def delta_lseries(prec=53, max_imaginary_part=0,
                  max_asymp_coeffs=40, algorithm=None):
    r"""
    Return the L-series of the modular form `\Delta`.

    If algorithm is "gp", this returns an interface to Tim
    Dokchitser's program for computing with the L-series of the
    modular form `\Delta`.

    If algorithm is "pari", this returns instead an interface to Pari's
    own general implementation of L-functions.

    INPUT:

    - ``prec`` -- integer (bits precision)

    - ``max_imaginary_part`` -- real number

    - ``max_asymp_coeffs`` -- integer

    - ``algorithm`` -- optional string: 'gp' (default), 'pari'

    OUTPUT:

    The L-series of `\Delta`.

    EXAMPLES::

        sage: L = delta_lseries()
        sage: L(1)
        0.0374412812685155

        sage: L = delta_lseries(algorithm='pari')
        sage: L(1)
        0.0374412812685155
    """
    if algorithm is None:
        algorithm = 'pari'

    if algorithm == 'gp':
        from sage.lfunctions.all import Dokchitser
        L = Dokchitser(conductor=1, gammaV=[0, 1], weight=12, eps=1,
                       prec=prec)
        s = 'tau(n) = (5*sigma(n,3)+7*sigma(n,5))*n/12-35*sum(k=1,n-1,(6*k-4*(n-k))*sigma(k,3)*sigma(n-k,5));'
        L.init_coeffs('tau(k)', pari_precode=s,
                      max_imaginary_part=max_imaginary_part,
                      max_asymp_coeffs=max_asymp_coeffs)
        L.set_coeff_growth('2*n^(11/2)')
        L.rename('L-series associated to the modular form Delta')
        return L
    elif algorithm == 'pari':
        from sage.lfunctions.pari import LFunction, lfun_delta
        return LFunction(lfun_delta(), prec=prec)

    raise ValueError('algorithm must be "gp" or "pari"')



class ModularForm_abstract(ModuleElement):
    """
    Constructor for generic class of a modular form. This
    should never be called directly; instead one should
    instantiate one of the derived classes of this
    class.
    """
    def group(self):
        """
        Return the group for which ``self`` is a modular form.

        EXAMPLES::

            sage: ModularForms(Gamma1(11), 2).gen(0).group()
            Congruence Subgroup Gamma1(11)
        """
        return self.parent().group()

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: (ModularForms(Gamma1(9),2).6).weight()
            2
        """
        return self.parent().weight()

    def level(self):
        """
        Return the level of self.

        EXAMPLES::

            sage: ModularForms(25,4).0.level()
            25
        """
        return self.parent().level()

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: ModularForms(25,4).0._repr_()
            'q + O(q^6)'

            sage: ModularForms(25,4).3._repr_()
            'q^4 + O(q^6)'
        """
        return str(self.q_expansion())

    def __call__(self, x, prec=None):
        """
        Evaluate the q-expansion of this modular form at x.

        EXAMPLES::

            sage: f = ModularForms(DirichletGroup(17).0^2,2).2

            sage: q = f.q_expansion().parent().gen()
            sage: f(q^2 + O(q^7))
            q^2 + (-zeta8^2 + 2)*q^4 + (zeta8 + 3)*q^6 + O(q^7)

            sage: f(0)
            0
        """
        return self.q_expansion(prec)(x)

    @cached_method
    def valuation(self):
        """
        Return the valuation of self (i.e. as an element of the power
        series ring in q).

        EXAMPLES::

            sage: ModularForms(11,2).0.valuation()
            1
            sage: ModularForms(11,2).1.valuation()
            0
            sage: ModularForms(25,6).1.valuation()
            2
            sage: ModularForms(25,6).6.valuation()
            7
        """
        v = self.qexp().valuation()
        if v != self.qexp().prec():
            return v
        return self.qexp(self.parent().sturm_bound()).valuation()

    def qexp(self, prec=None):
        """
        Same as ``self.q_expansion(prec)``.

        .. SEEALSO:: :meth:`q_expansion`

        EXAMPLES::

            sage: CuspForms(1,12).0.qexp()
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
        """
        return self.q_expansion(prec)

    def __eq__(self, other):
        """
        Compare self to other.

        EXAMPLES::

            sage: f = ModularForms(6,4).0
            sage: g = ModularForms(23,2).0
            sage: f == g ## indirect doctest
            False
            sage: f == f
            True
            sage: f == loads(dumps(f))
            True
        """
        if not isinstance(other, ModularFormElement) or \
           self.ambient_module() != other.ambient_module():
            return False
        else:
            return self.element() == other.element()

    def __ne__(self, other):
        """
        Return True if ``self != other``.

        EXAMPLES::

            sage: f = Newforms(Gamma1(30), 2, names='a')[1]
            sage: g = ModularForms(23, 2).0
            sage: f != g
            True
            sage: f != f
            False

        TESTS:

        The following used to fail (see :trac:`18068`)::

            sage: f != loads(dumps(f))
            False
        """
        return not (self == other)

    def _compute(self, X):
        """
        Compute the coefficients of `q^n` of the power series of self,
        for `n` in the list `X`.  The results are not cached.  (Use
        coefficients for cached results).

        EXAMPLES::

            sage: f = ModularForms(18,2).1
            sage: f.q_expansion(20)
            q + 8*q^7 + 4*q^10 + 14*q^13 - 4*q^16 + 20*q^19 + O(q^20)
            sage: f._compute([10,17])
            [4, 0]
            sage: f._compute([])
            []
        """
        if not isinstance(X, list) or not X:
            return []
        bound = max(X)
        q_exp = self.q_expansion(bound + 1)
        return [q_exp[i] for i in X]

    def coefficients(self, X):
        """
        The coefficients a_n of self, for integers n>=0 in the list
        X. If X is an Integer, return coefficients for indices from 1
        to X.

        This function caches the results of the compute function.

        TESTS::

            sage: e = DirichletGroup(11).gen()
            sage: f = EisensteinForms(e, 3).eisenstein_series()[0]
            sage: f.coefficients([0,1])
            [15/11*zeta10^3 - 9/11*zeta10^2 - 26/11*zeta10 - 10/11,
            1]
            sage: f.coefficients([0,1,2,3])
            [15/11*zeta10^3 - 9/11*zeta10^2 - 26/11*zeta10 - 10/11,
            1,
            4*zeta10 + 1,
            -9*zeta10^3 + 1]
            sage: f.coefficients([2,3])
            [4*zeta10 + 1,
            -9*zeta10^3 + 1]

        Running this twice once revealed a bug, so we test it::

            sage: f.coefficients([0,1,2,3])
            [15/11*zeta10^3 - 9/11*zeta10^2 - 26/11*zeta10 - 10/11,
            1,
            4*zeta10 + 1,
            -9*zeta10^3 + 1]
        """
        try:
            self.__coefficients
        except AttributeError:
            self.__coefficients = {}
        if isinstance(X, Integer):
            X = list(range(1, X + 1))
        Y = [n for n in X  if n not in self.__coefficients]
        v = self._compute(Y)
        for i in range(len(v)):
            self.__coefficients[Y[i]] = v[i]
        return [ self.__coefficients[x] for x in X ]

    def __getitem__(self, n):
        """
        Returns the `q^n` coefficient of the `q`-expansion of self or
        returns a list containing the `q^i` coefficients of self
        where `i` is in slice `n`.

        EXAMPLES::

            sage: f = ModularForms(DirichletGroup(17).0^2,2).2
            sage: f.__getitem__(10)
            zeta8^3 - 5*zeta8^2 - 2*zeta8 + 10
            sage: f[30]
            -2*zeta8^3 - 17*zeta8^2 + 4*zeta8 + 29
            sage: f[10:15]
            [zeta8^3 - 5*zeta8^2 - 2*zeta8 + 10,
            -zeta8^3 + 11,
            -2*zeta8^3 - 6*zeta8^2 + 3*zeta8 + 9,
            12,
            2*zeta8^3 - 7*zeta8^2 + zeta8 + 14]
        """
        if isinstance(n, slice):
            if n.stop is None:
                return self.q_expansion().list()[n]
            else:
                return self.q_expansion(n.stop+1).list()[n]
        else:
            return self.q_expansion(n+1)[int(n)]

    def padded_list(self, n):
        """
        Return a list of length n whose entries are the first n
        coefficients of the q-expansion of self.

        EXAMPLES::

            sage: CuspForms(1,12).0.padded_list(20)
            [0, 1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920, 534612, -370944, -577738, 401856, 1217160, 987136, -6905934, 2727432, 10661420]
        """
        return self.q_expansion(n).padded_list(n)


    def _latex_(self):
        """
        Return the LaTeX expression of self.

        EXAMPLES::

            sage: ModularForms(25,4).0._latex_()
            'q + O(q^{6})'

            sage: ModularForms(25,4).4._latex_()
            'q^{5} + O(q^{6})'
        """
        return self.q_expansion()._latex_()

    def base_ring(self):
        """
        Return the base_ring of self.

        EXAMPLES::

            sage: (ModularForms(117, 2).13).base_ring()
            Rational Field
            sage: (ModularForms(119, 2, base_ring=GF(7)).12).base_ring()
            Finite Field of size 7
        """
        return self.parent().base_ring()

    def character(self, compute=True):
        """
        Return the character of self. If ``compute=False``, then this will
        return None unless the form was explicitly created as an element of a
        space of forms with character, skipping the (potentially expensive)
        computation of the matrices of the diamond operators.

        EXAMPLES::

            sage: ModularForms(DirichletGroup(17).0^2,2).2.character()
            Dirichlet character modulo 17 of conductor 17 mapping 3 |--> zeta8

            sage: CuspForms(Gamma1(7), 3).gen(0).character()
            Dirichlet character modulo 7 of conductor 7 mapping 3 |--> -1
            sage: CuspForms(Gamma1(7), 3).gen(0).character(compute = False) is None
            True
            sage: M = CuspForms(Gamma1(7), 5).gen(0).character()
            Traceback (most recent call last):
            ...
            ValueError: Form is not an eigenvector for <3>
        """
        chi = self.parent().character()
        if (chi is not None) or (not compute):
            return chi
        else: # do the expensive computation
            G = DirichletGroup(self.parent().level(), base_ring = self.parent().base_ring())
            gens = G.unit_gens()
            i = self.valuation()
            vals = []
            for g in gens:
                df = self.parent().diamond_bracket_operator(g)(self)
                if df != (df[i] / self[i]) * self:
                    raise ValueError("Form is not an eigenvector for <%s>" % g)
                vals.append(df[i] / self[i])
            return G(vals)

    def __bool__(self):
        """
        Return ``True`` if ``self`` is nonzero, and ``False`` if not.

        EXAMPLES::

            sage: bool(ModularForms(25,6).6)
            True
        """
        return not self.element().is_zero()

    __nonzero__ = __bool__

    def prec(self):
        """
        Return the precision to which self.q_expansion() is
        currently known. Note that this may be 0.

        EXAMPLES::

            sage: M = ModularForms(2,14)
            sage: f = M.0
            sage: f.prec()
            0

            sage: M.prec(20)
            20
            sage: f.prec()
            0
            sage: x = f.q_expansion() ; f.prec()
            20
        """
        try:
            return self.__q_expansion[0]
        except AttributeError:
            return 0

    def q_expansion(self, prec=None):
        r"""
        The `q`-expansion of the modular form to precision `O(q^\text{prec})`.
        This function takes one argument, which is the integer prec.

        EXAMPLES:

        We compute the cusp form `\Delta`::

            sage: delta = CuspForms(1,12).0
            sage: delta.q_expansion()
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)

        We compute the `q`-expansion of one of the cusp forms of level 23::

            sage: f = CuspForms(23,2).0
            sage: f.q_expansion()
            q - q^3 - q^4 + O(q^6)
            sage: f.q_expansion(10)
            q - q^3 - q^4 - 2*q^6 + 2*q^7 - q^8 + 2*q^9 + O(q^10)
            sage: f.q_expansion(2)
            q + O(q^2)
            sage: f.q_expansion(1)
            O(q^1)
            sage: f.q_expansion(0)
            O(q^0)
            sage: f.q_expansion(-1)
            Traceback (most recent call last):
            ...
            ValueError: prec (= -1) must be non-negative
        """
        if prec is None:
            prec = self.parent().prec()
        prec = Integer(prec)
        try:
            current_prec, f = self.__q_expansion
        except AttributeError:
            current_prec = 0
            f = self.parent()._q_expansion_ring()(0, 0)

        if current_prec == prec:
            return f
        elif current_prec > prec:
            return f.add_bigoh(prec)
        else:
            f = self._compute_q_expansion(prec)
            self.__q_expansion = (prec, f)
            return f

    def atkin_lehner_eigenvalue(self, d=None, embedding=None):
        """
        Return the eigenvalue of the Atkin-Lehner operator `W_d`
        acting on ``self``.

        INPUT:

        - ``d`` -- a positive integer exactly dividing the level `N`
          of ``self``, i.e. `d` divides `N` and is coprime to `N/d`
          (default: `d = N`)

        - ``embedding`` -- (optional) embedding of the base ring of
          ``self`` into another ring

        OUTPUT:

        The Atkin-Lehner eigenvalue of `W_d` on ``self``.  This is
        returned as an element of the codomain of ``embedding`` if
        specified, and in (a suitable extension of) the base field of
        ``self`` otherwise.

        If ``self`` is not an eigenform for `W_d`, a ``ValueError`` is
        raised.

        .. SEEALSO::

            :meth:`sage.modular.hecke.module.HeckeModule_free_module.atkin_lehner_operator`
            (especially for the conventions used to define the operator `W_d`).

        EXAMPLES::

            sage: CuspForms(1, 12).0.atkin_lehner_eigenvalue()
            1
            sage: CuspForms(2, 8).0.atkin_lehner_eigenvalue()
            Traceback (most recent call last):
            ...
            NotImplementedError: Don't know how to compute Atkin-Lehner matrix acting on this space (try using a newform constructor instead)
        """
        raise NotImplementedError("Don't know how to compute Atkin-Lehner matrix acting on this space" \
            + " (try using a newform constructor instead)")

    # The methods period() and lseries() below currently live
    # in ModularForm_abstract so they are inherited by Newform (which
    # does *not* derive from ModularFormElement).

    def period(self, M, prec=53):
        r"""
        Return the period of ``self`` with respect to `M`.

        INPUT:

        - ``self`` -- a cusp form `f` of weight 2 for `Gamma_0(N)`

        - ``M`` -- an element of `\Gamma_0(N)`

        - ``prec`` -- (default: 53) the working precision in bits.  If
          `f` is a normalised eigenform, then the output is correct to
          approximately this number of bits.

        OUTPUT:

        A numerical approximation of the period `P_f(M)`.  This period
        is defined by the following integral over the complex upper
        half-plane, for any `\alpha` in `\Bold{P}^1(\QQ)`:

        .. MATH::

            P_f(M) = 2 \pi i \int_\alpha^{M(\alpha)} f(z) dz.

        This is independent of the choice of `\alpha`.

        EXAMPLES::

            sage: C = Newforms(11, 2)[0]
            sage: m = C.group()(matrix([[-4, -3], [11, 8]]))
            sage: C.period(m)
            -0.634604652139776 - 1.45881661693850*I

            sage: f = Newforms(15, 2)[0]
            sage: g = Gamma0(15)(matrix([[-4, -3], [15, 11]]))
            sage: f.period(g)  # abs tol 1e-15
            2.17298044293747e-16 - 1.59624222213178*I

        If `E` is an elliptic curve over `\QQ` and `f` is the newform
        associated to `E`, then the periods of `f` are in the period
        lattice of `E` up to an integer multiple::

            sage: E = EllipticCurve('11a3')
            sage: f = E.newform()
            sage: g = Gamma0(11)([3, 1, 11, 4])
            sage: f.period(g)
            0.634604652139777 + 1.45881661693850*I
            sage: omega1, omega2 = E.period_lattice().basis()
            sage: -2/5*omega1 + omega2
            0.634604652139777 + 1.45881661693850*I

        The integer multiple is 5 in this case, which is explained by
        the fact that there is a 5-isogeny between the elliptic curves
        `J_0(5)` and `E`.

        The elliptic curve `E` has a pair of modular symbols attached
        to it, which can be computed using the method
        :meth:`sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.modular_symbol`.
        These can be used to express the periods of `f` as exact
        linear combinations of the real and the imaginary period of `E`::

            sage: s = E.modular_symbol(sign=+1)
            sage: t = E.modular_symbol(sign=-1, implementation="sage")
            sage: s(3/11), t(3/11)
            (1/10, 1/2)
            sage: s(3/11)*omega1 + t(3/11)*2*omega2.imag()*I
            0.634604652139777 + 1.45881661693850*I

        ALGORITHM:

        We use the series expression from [Cre1997]_, Chapter II,
        Proposition 2.10.3.  The algorithm sums the first `T` terms of
        this series, where `T` is chosen in such a way that the result
        would approximate `P_f(M)` with an absolute error of at most
        `2^{-\text{prec}}` if all computations were done exactly.

        Since the actual precision is finite, the output is currently
        *not* guaranteed to be correct to ``prec`` bits of precision.

        TESTS::

            sage: C = Newforms(11, 2)[0]
            sage: g = Gamma0(15)(matrix([[-4, -3], [15, 11]]))
            sage: C.period(g)
            Traceback (most recent call last):
            ...
            TypeError: matrix [-4 -3]
                              [15 11]
            is not an element of Congruence Subgroup Gamma0(11)

            sage: f = Newforms(Gamma0(15), 4)[0]
            sage: f.period(g)
            Traceback (most recent call last):
            ...
            ValueError: period pairing only defined for cusp forms of weight 2

            sage: S = Newforms(Gamma1(17), 2, names='a')
            sage: f = S[1]
            sage: g = Gamma1(17)([18, 1, 17, 1])
            sage: f.period(g)
            Traceback (most recent call last):
            ...
            NotImplementedError: period pairing only implemented for cusp forms of trivial character

            sage: E = ModularForms(Gamma0(4), 2).eisenstein_series()[0]
            sage: gamma = Gamma0(4)([1, 0, 4, 1])
            sage: E.period(gamma)
            Traceback (most recent call last):
            ...
            NotImplementedError: Don't know how to compute Atkin-Lehner matrix acting on this space (try using a newform constructor instead)

            sage: E = EllipticCurve('19a1')
            sage: M = Gamma0(19)([10, 1, 19, 2])
            sage: E.newform().period(M)  # abs tol 1e-14
            -1.35975973348831 + 1.09365931898146e-16*I

        """
        R = RealField(prec)

        N = self.level()
        if not self.character().is_trivial():
            raise NotImplementedError('period pairing only implemented for cusp forms of trivial character')
        if self.weight() != 2:
            raise ValueError('period pairing only defined for cusp forms of weight 2')
        if not isinstance(self, (Newform, ModularFormElement_elliptic_curve)):
            print('Warning: not a newform, precision not guaranteed')

        M = self.group()(M)
        # coefficients of the matrix M
        (b, c, d) = (M.b(), M.c() / N, M.d())
        if d == 0:
            return R.zero()
        if d < 0:
            (b, c, d) = (-b, -c, -d)

        twopi = 2 * R.pi()
        I = R.complex_field().gen()
        rootN = R(N).sqrt()

        eps = self.atkin_lehner_eigenvalue()
        mu_N = (-twopi / rootN).exp()
        mu_dN = (-twopi / d / rootN).exp()
        mu_d = (twopi * I / d).exp()

        # We bound the tail of the series by means of the triangle
        # inequality and the following bounds (tau(n) = #divisors(n)):
        #       mu_N <= mu_dN
        #   |a_n(f)| <= tau(n)*sqrt(n)  (holds if f is a newform)
        #     tau(n) <= sqrt(3)*sqrt(n) for all n >= 1
        # This gives a correct but somewhat coarse lower bound on the
        # number of terms needed.  We ignore rounding errors.
        numterms = (((1 - mu_dN) * R(2)**(-prec)
                     / ((abs(eps - 1) + 2) * R(3).sqrt())).log()
                    / mu_dN.log()).ceil()
        coeff = self.coefficients(numterms)

        return sum((coeff[n - 1] / n)
                   *((eps - 1) * mu_N ** n
                     + mu_dN ** n * (mu_d ** (n * b) - eps * mu_d ** (n * c)))
                   for n in range(1, numterms + 1))

    def lseries(self, embedding=0, prec=53,
                         max_imaginary_part=0,
                         max_asymp_coeffs=40,
                         conjugate=None):
        r"""
        Return the L-series of the weight k cusp form
        `f` on `\Gamma_0(N)`.

        This actually returns an interface to Tim Dokchitser's program for
        computing with the L-series of the cusp form.

        INPUT:

        - ``embedding`` - either an embedding of the coefficient field of self
          into `\CC`, or an integer `i` between 0 and D-1 where D is the degree
          of the coefficient field (meaning to pick the `i`-th embedding).
          (Default: 0)

        - ``prec`` - integer (bits precision). Default: 53.

        - ``max_imaginary_part`` - real number. Default: 0.

        - ``max_asymp_coeffs`` - integer. Default: 40.

        - ``conjugate`` -- deprecated synonym for ``embedding``.

        For more information on the significance of the last three arguments,
        see :mod:`~sage.lfunctions.dokchitser`.

        .. note::

            If an explicit embedding is given, but this embedding is specified
            to smaller precision than ``prec``, it will be automatically
            refined to precision ``prec``.

        OUTPUT:

        The L-series of the cusp form, as a
        :class:`sage.lfunctions.dokchitser.Dokchitser` object.

        EXAMPLES::

            sage: f = CuspForms(2,8).newforms()[0]
            sage: L = f.lseries()
            sage: L
            L-series associated to the cusp form q - 8*q^2 + 12*q^3 + 64*q^4 - 210*q^5 + O(q^6)
            sage: L(1)
            0.0884317737041015
            sage: L(0.5)
            0.0296568512531983

        As a consistency check, we verify that the functional equation holds::

            sage: abs(L.check_functional_equation()) < 1.0e-20
            True

        For non-rational newforms we can specify an embedding of the coefficient field::

            sage: f = Newforms(43, names='a')[1]
            sage: K = f.hecke_eigenvalue_field()
            sage: phi1, phi2 = K.embeddings(CC)
            sage: L = f.lseries(embedding=phi1)
            sage: L
            L-series associated to the cusp form q + a1*q^2 - a1*q^3 + (-a1 + 2)*q^5 + O(q^6), a1=-1.41421356237310
            sage: L(1)
            0.620539857407845
            sage: L = f.lseries(embedding=1)
            sage: L(1)
            0.921328017272472

        For backward-compatibility, ``conjugate`` is accepted as a synonym for ``embedding``::

            sage: f.lseries(conjugate=1)
            doctest:...: DeprecationWarning: The argument 'conjugate' for 'lseries' is deprecated -- use the synonym 'embedding'
            See http://trac.sagemath.org/19668 for details.
            L-series associated to the cusp form q + a1*q^2 - a1*q^3 + (-a1 + 2)*q^5 + O(q^6), a1=1.41421356237310

        An example with a non-real coefficient field (`\QQ(\zeta_3)`
        in this case)::

            sage: f = Newforms(Gamma1(13), 2, names='a')[0]
            sage: f.lseries(embedding=0)(1)
            0.298115272465799 - 0.0402203326076734*I
            sage: f.lseries(embedding=1)(1)
            0.298115272465799 + 0.0402203326076732*I

        We compute with the L-series of the Eisenstein series `E_4`::

            sage: f = ModularForms(1,4).0
            sage: L = f.lseries()
            sage: L(1)
            -0.0304484570583933
            sage: L = eisenstein_series_lseries(4)
            sage: L(1)
            -0.0304484570583933

        Consistency check with delta_lseries (which computes coefficients in pari)::

            sage: delta = CuspForms(1,12).0
            sage: L = delta.lseries()
            sage: L(1)
            0.0374412812685155
            sage: L = delta_lseries()
            sage: L(1)
            0.0374412812685155

        We check that :trac:`5262` is fixed::

            sage: E = EllipticCurve('37b2')
            sage: h = Newforms(37)[1]
            sage: Lh = h.lseries()
            sage: LE = E.lseries()
            sage: Lh(1), LE(1)
            (0.725681061936153, 0.725681061936153)
            sage: CuspForms(1, 30).0.lseries().eps
            -1.00000000000000

        We check that :trac:`25369` is fixed::

            sage: f5 = Newforms(Gamma1(4), 5, names='a')[0]; f5
            q - 4*q^2 + 16*q^4 - 14*q^5 + O(q^6)
            sage: L5 = f5.lseries()
            sage: abs(L5.check_functional_equation()) < 1e-15
            True
            sage: abs(L5(4) - (gamma(1/4)^8/(3840*pi^2)).n()) < 1e-15
            True

        We can change the precision (in bits)::

            sage: f = Newforms(389, names='a')[0]
            sage: L = f.lseries(prec=30)
            sage: abs(L(1)) < 2^-30
            True
            sage: L = f.lseries(prec=53)
            sage: abs(L(1)) < 2^-53
            True
            sage: L = f.lseries(prec=100)
            sage: abs(L(1)) < 2^-100
            True

            sage: f = Newforms(27, names='a')[0]
            sage: L = f.lseries()
            sage: L(1)
            0.588879583428483
        """
        from sage.lfunctions.all import Dokchitser

        # compute the requested embedding
        C = ComplexField(prec)
        if conjugate is not None:
            from sage.misc.superseded import deprecation
            deprecation(19668, "The argument 'conjugate' for 'lseries' is deprecated -- use the synonym 'embedding'")
            embedding=conjugate
        K = self.base_ring()
        if isinstance(embedding, RingHomomorphism):
            # Target of embedding might have precision less than desired, so
            # need to refine
            emb = NumberFieldEmbedding(K, C, embedding(K.gen()))
        else:
            emb = self.base_ring().embeddings(C)[embedding]

        # key = (prec, max_imaginary_part, max_asymp_coeffs)
        l = self.weight()
        N = self.level()

        # get global root number
        w = self.atkin_lehner_eigenvalue(N, embedding=emb)
        e = ~C.gen()**l * w

        if self.is_cuspidal():
            poles = []  # cuspidal
        else:
            poles = [l] # non-cuspidal

        L = Dokchitser(conductor = N,
                       gammaV = [0, 1],
                       weight = l,
                       eps = e,
                       poles = poles,
                       prec = prec)
        # Find out how many coefficients of the Dirichlet series are needed
        # in order to compute to the required precision
        num_coeffs = L.num_coeffs()
        coeffs = self.q_expansion(num_coeffs+1).padded_list()[1:]

        # renormalize so that coefficient of q is 1
        b = coeffs[0]
        if b != 1:
            invb = 1/b
            coeffs = (invb*c for c in coeffs)

        v = [emb(c) for c in coeffs]
        w = [c.conjugate() for c in v]
        L.init_coeffs(v=v, w=w,
                      max_imaginary_part=max_imaginary_part,
                      max_asymp_coeffs=max_asymp_coeffs)
        L.check_functional_equation()
        if K == QQ:
            L.rename('L-series associated to the cusp form %s'%self)
        else:
            L.rename('L-series associated to the cusp form %s, %s=%s' \
                % (self, K.variable_name(), emb(K.gen())))
        return L

    cuspform_lseries = deprecated_function_alias(16917, lseries)

    def symsquare_lseries(self, chi=None, embedding=0, prec=53):
        r"""
        Compute the symmetric square L-series of this modular form, twisted by
        the character `\chi`.

        INPUT:

        - ``chi`` -- Dirichlet character to twist by, or None (default None,
          interpreted as the trivial character).
        - ``embedding`` -- embedding of the coefficient field into `\RR` or
          `\CC`, or an integer `i` (in which case take the `i`-th embedding)
        - ``prec`` -- The desired precision in bits (default 53).

        OUTPUT: The symmetric square L-series of the cusp form, as a
        :class:`sage.lfunctions.dokchitser.Dokchitser` object.

        EXAMPLES::

            sage: CuspForms(1, 12).0.symsquare_lseries()(22)
            0.999645711124771

        An example twisted by a nontrivial character::

            sage: psi = DirichletGroup(7).0^2
            sage: L = CuspForms(1, 16).0.symsquare_lseries(psi)
            sage: L(22)
            0.998407750967420 - 0.00295712911510708*I

        An example with coefficients not in `\QQ`::

            sage: F = Newforms(1, 24, names='a')[0]
            sage: K = F.hecke_eigenvalue_field()
            sage: phi = K.embeddings(RR)[0]
            sage: L = F.symsquare_lseries(embedding=phi)
            sage: L(5)
            verbose -1 (...: dokchitser.py, __call__) Warning: Loss of 8 decimal digits due to cancellation
            -3.57698266793901e19

        TESTS::

            sage: CuspForms(1,16).0.symsquare_lseries(prec=200).check_functional_equation().abs() < 1.0e-80
            True
            sage: CuspForms(1, 12).0.symsquare_lseries(prec=1000)(22) # long time (20s)
            0.999645711124771397835729622033153189549796658647254961493709341358991830134499267117001769570658192128781135161587571716303826382489492569725002840546129937149159065273765309218543427544527498868033604310899372849565046516553245752253255585377793879866297612679545029546953895098375829822346290125161

        Check that :trac:`23247` is fixed::

            sage: F = Newforms(1,12)[0]
            sage: chi = DirichletGroup(7).0
            sage: abs(F.symsquare_lseries(chi).check_functional_equation()) < 1e-5
            True

        AUTHORS:

        - Martin Raum (2011) -- original code posted to sage-nt
        - David Loeffler (2015) -- added support for twists, integrated into
          Sage library
        """
        from sage.lfunctions.all import Dokchitser
        weight = self.weight()
        C = ComplexField(prec)
        if self.level() != 1:
            raise NotImplementedError("Symmetric square L-functions only implemented for level 1")

        # compute the requested embedding
        if isinstance(embedding, RingHomomorphism):
            # Target of embedding might have precision less than desired, so
            # need to refine
            K = self.base_ring()
            emb = NumberFieldEmbedding(K, ComplexField(prec), embedding(K.gen()))
        else:
            emb = self.base_ring().embeddings(ComplexField(prec))[embedding]

        if chi is None:
            eps = 1
            N = 1
        else:
            assert chi.is_primitive()
            chi = chi.change_ring(C)
            eps = chi.gauss_sum()**3 / chi.base_ring()(chi.conductor())**QQ( (3, 2) )
            N = chi.conductor()**3

        if (chi is None) or chi.is_even():
            L = Dokchitser(N, [0, 1, -weight + 2], 2 * weight - 1,
                           eps, prec=prec)
        else:
            L = Dokchitser(N, [0, 1, -weight + 1], 2 * weight - 1,
                           eps * C((0, 1)), prec=prec)
        lcoeffs_prec = L.num_coeffs()

        t = verbose("Computing %s coefficients of F" % lcoeffs_prec, level=1)
        F_series = [u**2 for u in self.qexp(lcoeffs_prec + 1).list()[1:]]
        verbose("done", t, level=1)

        # utility function for Dirichlet convolution of series
        def dirichlet_convolution(A, B):
            return [sum(A[d-1] * B[n/d - 1] for d in divisors(n))
                for n in range(1, 1 + min(len(A), len(B)))]

        # The Dirichlet series for \zeta(2 s - 2 k + 2)
        riemann_series = [ n**(weight - 1) if n.is_square() else 0
                       for n in xsrange(1, lcoeffs_prec + 1) ]
        # The Dirichlet series for 1 / \zeta(s - k + 1)
        mu_series = [ moebius(n) * n**(weight - 1) for n in xsrange(1, lcoeffs_prec + 1) ]
        conv_series = dirichlet_convolution(mu_series, riemann_series)

        dirichlet_series = dirichlet_convolution(conv_series, F_series)

        # If the base ring is QQ we pass the coefficients to GP/PARI as exact
        # rationals. Otherwise, need to use the embedding.
        if self.base_ring() != QQ:
            dirichlet_series = [emb(cf) for cf in dirichlet_series]

        if chi is not None:
            pari_precode_chi = str(chi.values()) + "[n%" + str(chi.conductor()) + "+1]; "
        else:
            pari_precode_chi = "1"

        pari_precode = "hhh(n) = " + str(dirichlet_series) + "[n] * " + pari_precode_chi

        L.init_coeffs( "hhh(k)", w="conj(hhh(k))",
            pari_precode=pari_precode)

        return L

    def petersson_norm(self, embedding=0, prec=53):
        r"""
        Compute the Petersson scalar product of f with itself:

        .. MATH::

            \langle f, f \rangle = \int_{\Gamma_0(N) \backslash \mathbb{H}} |f(x + iy)|^2 y^k\, \mathrm{d}x\, \mathrm{d}y.

        Only implemented for N = 1 at present. It is assumed that `f` has real
        coefficients. The norm is computed as a special value of the symmetric
        square L-function, using the identity

        .. MATH::

            \langle f, f \rangle = \frac{(k-1)! L(\mathrm{Sym}^2 f, k)}{2^{2k-1} \pi^{k+1}}

        INPUT:

        - ``embedding``: embedding of the coefficient field into `\RR` or
          `\CC`, or an integer `i` (interpreted as the `i`-th embedding)
          (default: 0)
        - ``prec`` (integer, default 53): precision in bits

        EXAMPLES::

            sage: CuspForms(1, 16).0.petersson_norm()
            verbose -1 (...: dokchitser.py, __call__) Warning: Loss of 2 decimal digits due to cancellation
            2.16906134759063e-6

        The Petersson norm depends on a choice of embedding::

            sage: set_verbose(-2, "dokchitser.py") # disable precision-loss warnings
            sage: F = Newforms(1, 24, names='a')[0]
            sage: F.petersson_norm(embedding=0)
            0.000107836545077234
            sage: F.petersson_norm(embedding=1)
            0.000128992800758160

        TESTS:

        Verify that the Petersson norm is a quadratic form::

            sage: F, G = CuspForms(1, 24).basis()
            sage: X = lambda u: u.petersson_norm(prec=100)
            sage: (X(F + G) + X(F - G) - 2*X(F) - 2*X(G)).abs() < 1e-25
            True
        """
        pi = RealField(prec).pi()
        L = self.symsquare_lseries(prec=prec, embedding=embedding)
        k = self.weight()
        return (ZZ(k - 1).factorial() / 2**(2*k - 1) / pi**(k+1)) * L(k).real_part()

    def _q_expansion_bound(self, eps):
        r"""
        This function takes as input a modular form, ``self`` and a
        Dirichlet character ``eps`` and returns an integer bound such
        that if ``self`` and its twist by ``eps`` have the same
        q-expansion up to this bound, then they are equal.

        The bound is taken from [Mu1997]_. See also [Shi1971]_, Proposition
        3.64.

        INPUT:

        - ``eps`` -- a Dirichlet character

        OUTPUT:

        A positive integer.

        EXAMPLES:

        Here is an example that can easily be checked by hand. ::

            sage: M = ModularForms(Gamma0(11), 2)
            sage: C = M.cuspidal_submodule()
            sage: f = C.gens()[0]
            sage: F = CyclotomicField(5)
            sage: D = DirichletGroup(11, F)
            sage: eps = D.gens()[0]
            sage: f._q_expansion_bound(eps)
            22

        The level of ``self`` does not have to be related to the conductor
        of eps. ::

            sage: M = ModularForms(Gamma0(1), 12)
            sage: C = M.cuspidal_submodule()
            sage: Delta = C.gens()[0]
            sage: F = CyclotomicField(12)
            sage: D = DirichletGroup(13, F)
            sage: eps = D.gens()[0]
            sage: Delta._q_expansion_bound(eps)
            182
        """
        chi = self.character()
        M = lcm([self.level(), eps.conductor()**2,
                 chi.conductor() * eps.conductor()])
        y = QQ(self.weight()) / QQ(12) * M
        for p in M.prime_divisors():
            y *= (1 + 1/QQ(p))
        return y.ceil()

    @cached_method
    def has_cm(self):
        r"""
        Return whether the modular form ``self`` has complex multiplication.

        OUTPUT:

        Boolean

        .. SEEALSO::

            - :meth:`cm_discriminant` (to return the CM field)
            - :meth:`sage.schemes.elliptic_curves.ell_rational_field.has_cm`

        EXAMPLES::

            sage: G = DirichletGroup(21); eps = G.0 * G.1
            sage: Newforms(eps, 2)[0].has_cm()
            True

        This example illustrates what happens when
        candidate_characters(self) is the empty list. ::

            sage: M = ModularForms(Gamma0(1), 12)
            sage: C = M.cuspidal_submodule()
            sage: Delta = C.gens()[0]
            sage: Delta.has_cm()
            False

        We now compare the function has_cm between elliptic curves and
        their associated modular forms. ::

            sage: E = EllipticCurve([-1, 0])
            sage: f = E.modular_form()
            sage: f.has_cm()
            True
            sage: E.has_cm() == f.has_cm()
            True

        Here is a non-cm example coming from elliptic curves. ::

            sage: E = EllipticCurve('11a')
            sage: f = E.modular_form()
            sage: f.has_cm()
            False
            sage: E.has_cm() == f.has_cm()
            True
        """
        N = self.level()
        M = self.character().conductor()

        for p in N.prime_factors():
            if M % p and N.valuation(p) == 1:
                verbose("Form is Steinberg at %s, cannot be CM" % p, level=1)
                return False
        cand_chars = [(x, self._q_expansion_bound(x)) for x in DirichletGroup(N, QQ) if x.is_odd()]

        verbose("Conductors of candidate characters: %s" % (", ".join(str(x[0].conductor()) for x in cand_chars)), level=1)
        verbose("Qexp bounds: %s" % (", ".join(str(x[1]) for x in cand_chars)), level=1)
        # If there are no candidate characters, then self cannot have CM.
        if not cand_chars:
            return False

        # Test each prime and discard characters for which eps(p) != 1 when f[p] != 0.
        p = ZZ(2)
        while p <= min(B for (eps, B) in cand_chars):
            verbose("Checking p = %s (%s candidate characters left))" % (p, len(cand_chars)), level=1)
            # We only have to test the CM condition at primes that do not
            # divide the level of self.
            if not self.level() % p:
                p = p.next_prime()
                continue

            # Evaluating characters is cheap, while computing f[p] is
            # expensive, so if eps(p) = 1 for all p, then we don't bother to
            # compute f[p].
            cand_chars = [(eps, B) for (eps, B) in cand_chars if (eps(p) == 1 or self[p] == 0)]

            if len(cand_chars) == 0:
                # f doesn't have CM
                return False

            # go on to next prime
            p = p.next_prime()

        B0 = min(B for (eps, B) in cand_chars)
        C = [eps for (eps, B) in cand_chars if B == B0]
        if len(C) > 1:
            # can't happen (except in weight 1, which isn't implemented yet
            # anyway)
            raise ArithmeticError("Got multiple characters in has_cm")
        self.__cm_char = C[0].primitive_character()
        return True

    def cm_discriminant(self):
        r"""
        Return the discriminant of the CM field associated to this form. An
        error will be raised if the form isn't of CM type.

        EXAMPLES::

            sage: Newforms(49, 2)[0].cm_discriminant()
            -7
            sage: CuspForms(1, 12).gen(0).cm_discriminant()
            Traceback (most recent call last):
            ...
            ValueError: Not a CM form
        """
        if not self.has_cm():
            raise ValueError("Not a CM form")
        return -self.__cm_char.conductor()

class Newform(ModularForm_abstract):
    def __init__(self, parent, component, names, check=True):
        r"""
        Initialize a Newform object.

        INPUT:

        - ``parent`` - An ambient cuspidal space of modular forms for
          which self is a newform.

        - ``component`` - A simple component of a cuspidal modular
          symbols space of any sign corresponding to this newform.

        - ``check`` - If check is ``True``, check that parent and
          component have the same weight, level, and character, that
          component has sign 1 and is simple, and that the types are
          correct on all inputs.

        EXAMPLES::

            sage: sage.modular.modform.element.Newform(CuspForms(11,2), ModularSymbols(11,2,sign=1).cuspidal_subspace(), 'a')
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)

            sage: f = Newforms(DirichletGroup(5).0, 7,names='a')[0]; f[2].trace(f.base_ring().base_field())
            -5*zeta4 - 5
        """
        from .space import is_ModularFormsSpace
        if check:
            if not is_ModularFormsSpace(parent):
                raise TypeError("parent must be a space of modular forms")
            if not is_ModularSymbolsSpace(component):
                raise TypeError("component must be a space of modular symbols")
            if parent.group() != component.group():
                raise ValueError("parent and component must be defined by the same congruence subgroup")
            if parent.weight() != component.weight():
                raise ValueError("parent and component must have the same weight")
            if not component.is_cuspidal():
                raise ValueError("component must be cuspidal")
            if not component.is_simple():
                raise ValueError("component must be simple")
        extension_field = component.eigenvalue(1,name=names).parent()
        if extension_field != parent.base_ring(): # .degree() != 1 and rings.is_NumberField(extension_field):
            assert extension_field.base_field() == parent.base_ring()
            extension_field = parent.base_ring().extension(extension_field.relative_polynomial(), names=names)
        self.__name = names
        ModuleElement.__init__(self, parent.base_extend(extension_field))
        self.__modsym_space = component
        self.__hecke_eigenvalue_field = extension_field

    def _name(self):
        """
        Return the name of the generator of the Hecke eigenvalue field
        of self. Note that a name exists even when this field is QQ.

        EXAMPLES::

            sage: [ f._name() for f in Newforms(38,4,names='a') ]
            ['a0', 'a1', 'a2']
        """
        return self.__name

    def _compute_q_expansion(self, prec):
        """
        Return the q-expansion of self to precision prec.

        EXAMPLES::

            sage: forms = Newforms(31, 6, names='a')
            sage: forms[0]._compute_q_expansion(10)
            q + a0*q^2 + (5/704*a0^4 + 43/704*a0^3 - 61/88*a0^2 - 197/44*a0 + 717/88)*q^3 + (a0^2 - 32)*q^4 + (-31/352*a0^4 - 249/352*a0^3 + 111/22*a0^2 + 218/11*a0 - 2879/44)*q^5 + (-1/352*a0^4 - 79/352*a0^3 - 67/44*a0^2 + 13/22*a0 - 425/44)*q^6 + (17/88*a0^4 + 133/88*a0^3 - 405/44*a0^2 - 1005/22*a0 - 35/11)*q^7 + (a0^3 - 64*a0)*q^8 + (39/352*a0^4 + 441/352*a0^3 - 93/44*a0^2 - 441/22*a0 - 5293/44)*q^9 + O(q^10)
            sage: forms[0]._compute_q_expansion(15)
            q + a0*q^2 + (5/704*a0^4 + 43/704*a0^3 - 61/88*a0^2 - 197/44*a0 + 717/88)*q^3 + (a0^2 - 32)*q^4 + (-31/352*a0^4 - 249/352*a0^3 + 111/22*a0^2 + 218/11*a0 - 2879/44)*q^5 + (-1/352*a0^4 - 79/352*a0^3 - 67/44*a0^2 + 13/22*a0 - 425/44)*q^6 + (17/88*a0^4 + 133/88*a0^3 - 405/44*a0^2 - 1005/22*a0 - 35/11)*q^7 + (a0^3 - 64*a0)*q^8 + (39/352*a0^4 + 441/352*a0^3 - 93/44*a0^2 - 441/22*a0 - 5293/44)*q^9 + (15/176*a0^4 - 135/176*a0^3 - 185/11*a0^2 + 311/11*a0 + 2635/22)*q^10 + (-291/704*a0^4 - 3629/704*a0^3 + 1139/88*a0^2 + 10295/44*a0 - 21067/88)*q^11 + (-75/176*a0^4 - 645/176*a0^3 + 475/22*a0^2 + 1503/11*a0 - 5651/22)*q^12 + (207/704*a0^4 + 2977/704*a0^3 + 581/88*a0^2 - 3307/44*a0 - 35753/88)*q^13 + (-5/22*a0^4 + 39/11*a0^3 + 763/22*a0^2 - 2296/11*a0 - 2890/11)*q^14 + O(q^15)
        """
        return self.modular_symbols(1).q_eigenform(prec, names=self._name())

    def __eq__(self, other):
        """
        Return True if self equals other, and False otherwise.

        EXAMPLES::

            sage: f1, f2 = Newforms(17,4,names='a')
            sage: f1.__eq__(f1)
            True
            sage: f1.__eq__(f2)
            False

        We test comparison of equal newforms with different parents
        (see :trac:`18478`)::

            sage: f = Newforms(Gamma1(11), 2)[0]; f
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)
            sage: g = Newforms(Gamma0(11), 2)[0]; g
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)
            sage: f == g
            True

            sage: f = Newforms(DirichletGroup(4)[1], 5)[0]; f
            q - 4*q^2 + 16*q^4 - 14*q^5 + O(q^6)
            sage: g = Newforms(Gamma1(4), 5)[0]; g
            q - 4*q^2 + 16*q^4 - 14*q^5 + O(q^6)
            sage: f == g
            True

        """
        if (not isinstance(other, ModularForm_abstract)
            or self.weight() != other.weight()):
            return False
        if isinstance(other, Newform):
            if (self.level() != other.level() or
                self.character() != other.character()):
                return False
            # The two parents may have different Sturm bounds in case
            # one of them is a space of cusp forms with character
            # (possibly trivial, i.e. for the group Gamma0(n)) and the
            # other is a space of cusp forms for Gamma1(n).  It is
            # safe to take the smaller bound because we have checked
            # that the characters agree.
            bound = min(self.parent().sturm_bound(),
                        other.parent().sturm_bound())
            return self.q_expansion(bound) == other.q_expansion(bound)
        # other is a ModularFormElement
        return self.element() == other

    @cached_method
    def abelian_variety(self):
        """
        Return the abelian variety associated to self.

        EXAMPLES::

            sage: Newforms(14,2)[0]
            q - q^2 - 2*q^3 + q^4 + O(q^6)
            sage: Newforms(14,2)[0].abelian_variety()
            Newform abelian subvariety 14a of dimension 1 of J0(14)
            sage: Newforms(1, 12)[0].abelian_variety()
            Traceback (most recent call last):
            ...
            TypeError: f must have weight 2

        """
        from sage.modular.abvar.abvar_newform import ModularAbelianVariety_newform
        return ModularAbelianVariety_newform(self)

    def hecke_eigenvalue_field(self):
        r"""
        Return the field generated over the rationals by the
        coefficients of this newform.

        EXAMPLES::

            sage: ls = Newforms(35, 2, names='a') ; ls
            [q + q^3 - 2*q^4 - q^5 + O(q^6),
            q + a1*q^2 + (-a1 - 1)*q^3 + (-a1 + 2)*q^4 + q^5 + O(q^6)]
            sage: ls[0].hecke_eigenvalue_field()
            Rational Field
            sage: ls[1].hecke_eigenvalue_field()
            Number Field in a1 with defining polynomial x^2 + x - 4
        """
        return self.__hecke_eigenvalue_field

    def coefficient(self, n):
        """
        Return the coefficient of `q^n` in the power series of self.

        INPUT:

        - ``n`` - a positive integer

        OUTPUT:

        - the coefficient of `q^n` in the power series of self.

        EXAMPLES::

            sage: f = Newforms(11)[0]; f
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)
            sage: f.coefficient(100)
            -8

            sage: g = Newforms(23, names='a')[0]; g
            q + a0*q^2 + (-2*a0 - 1)*q^3 + (-a0 - 1)*q^4 + 2*a0*q^5 + O(q^6)
            sage: g.coefficient(3)
            -2*a0 - 1
        """
        return self.modular_symbols(1).eigenvalue(n, self._name())

    def _compute(self, X):
        """
        Compute the coefficients of `q^n` of the power series of self,
        for `n` in the list `X`.  The results are not cached.  (Use
        coefficients for cached results).

        EXAMPLES::

            sage: f = Newforms(39,4,names='a')[1] ; f
            q + a1*q^2 - 3*q^3 + (2*a1 + 5)*q^4 + (-2*a1 + 14)*q^5 + O(q^6)
            sage: f._compute([2,3,7])
            [a1, -3, -2*a1 + 2]
            sage: f._compute([])
            []

        Check that :trac:`20793` is fixed::

            sage: f = Newforms(83, 2, names='a')[1]; f
            q + a1*q^2 + (1/2*a1^4 - 1/2*a1^3 - 7/2*a1^2 + 3/2*a1 + 4)*q^3 + (a1^2 - 2)*q^4 + (-1/2*a1^5 - 1/2*a1^4 + 9/2*a1^3 + 7/2*a1^2 - 8*a1 - 2)*q^5 + O(q^6)
            sage: K = f.hecke_eigenvalue_field(); K
            Number Field in a1 with defining polynomial x^6 - x^5 - 9*x^4 + 7*x^3 + 20*x^2 - 12*x - 8
            sage: l = f.coefficients(20); l[-1]
            -a1^4 + 5*a1^2 - 4
            sage: l[-1].parent() is K
            True
        """
        M = self.modular_symbols(1)
        return [M.eigenvalue(x, name=self._name()) for x in X]

    def element(self):
        """
        Find an element of the ambient space of modular forms which
        represents this newform.

        .. note::

           This can be quite expensive. Also, the polynomial defining
           the field of Hecke eigenvalues should be considered random,
           since it is generated by a random sum of Hecke
           operators. (The field itself is not random, of course.)

        EXAMPLES::

            sage: ls = Newforms(38,4,names='a')
            sage: ls[0]
            q - 2*q^2 - 2*q^3 + 4*q^4 - 9*q^5 + O(q^6)
            sage: ls # random
            [q - 2*q^2 - 2*q^3 + 4*q^4 - 9*q^5 + O(q^6),
            q - 2*q^2 + (-a1 - 2)*q^3 + 4*q^4 + (2*a1 + 10)*q^5 + O(q^6),
            q + 2*q^2 + (1/2*a2 - 1)*q^3 + 4*q^4 + (-3/2*a2 + 12)*q^5 + O(q^6)]
            sage: type(ls[0])
            <class 'sage.modular.modform.element.Newform'>
            sage: ls[2][3].minpoly()
            x^2 - 9*x + 2
            sage: ls2 = [ x.element() for x in ls ]
            sage: ls2 # random
            [q - 2*q^2 - 2*q^3 + 4*q^4 - 9*q^5 + O(q^6),
            q - 2*q^2 + (-a1 - 2)*q^3 + 4*q^4 + (2*a1 + 10)*q^5 + O(q^6),
            q + 2*q^2 + (1/2*a2 - 1)*q^3 + 4*q^4 + (-3/2*a2 + 12)*q^5 + O(q^6)]
            sage: type(ls2[0])
            <class 'sage.modular.modform.cuspidal_submodule.CuspidalSubmodule_g0_Q_with_category.element_class'>
            sage: ls2[2][3].minpoly()
            x^2 - 9*x + 2
        """
        S = self.parent()
        return S(self.q_expansion(S.sturm_bound()))

    def is_cuspidal(self):
        """
        Return True. For compatibility with elements of modular forms spaces.

        EXAMPLES::

            sage: Newforms(11, 2)[0].is_cuspidal()
            True
        """
        return True

    def modular_symbols(self, sign=0):
        """
        Return the subspace with the specified sign of the space of
        modular symbols corresponding to this newform.

        EXAMPLES::

            sage: f = Newforms(18,4)[0]
            sage: f.modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 18 for Gamma_0(18) of weight 4 with sign 0 over Rational Field
            sage: f.modular_symbols(1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 11 for Gamma_0(18) of weight 4 with sign 1 over Rational Field
        """
        return self.__modsym_space.modular_symbols_of_sign(sign)

    @cached_method
    def modsym_eigenspace(self, sign=0):
        """
        Return a submodule of dimension 1 or 2 of the ambient space of
        the sign 0 modular symbols space associated to ``self``,
        base-extended to the Hecke eigenvalue field, which is an
        eigenspace for the Hecke operators with the same eigenvalues
        as this newform, *and* is an eigenspace for the star
        involution of the appropriate sign if the sign is not 0.

        EXAMPLES::

            sage: N = Newform("37a")
            sage: N.modular_symbols(0)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
            sage: M = N.modular_symbols(0)
            sage: V = N.modsym_eigenspace(1); V
            Vector space of degree 5 and dimension 1 over Rational Field
            Basis matrix:
            [ 0  1 -1  1  0]
            sage: V.0 in M.free_module()
            True
            sage: V=N.modsym_eigenspace(-1); V
            Vector space of degree 5 and dimension 1 over Rational Field
            Basis matrix:
            [   0    0    0    1 -1/2]
            sage: V.0 in M.free_module()
            True
        """
        M = self.modular_symbols(sign=0)
        if sign != 0:
            Ms = M.sign_submodule(sign)
            r = 1
        else:
            Ms = M
            r = 2
        # silly thing: can't do Ms.eigenvector(), even when Ms is simple,
        # because it can't be relied on to choose the coefficient fields
        # consistently
        A = M.ambient()
        X = Ms.free_module().base_extend(self.hecke_eigenvalue_field())
        p = ZZ(2)
        while X.rank() > r:
            assert p <= M.sturm_bound()
            X = (A.hecke_matrix(p).base_extend(self.hecke_eigenvalue_field()) - self[p]).kernel_on(X)
            p = p.next_prime()

        # should really return a modular symbol submodule object, but these are
        # not implemented over non-minimal base rings
        return X

    def _defining_modular_symbols(self):
        """
        Return the modular symbols space corresponding to self.

        EXAMPLES::

            sage: Newforms(43,2,names='a')
            [q - 2*q^2 - 2*q^3 + 2*q^4 - 4*q^5 + O(q^6),
            q + a1*q^2 - a1*q^3 + (-a1 + 2)*q^5 + O(q^6)]
            sage: [ x._defining_modular_symbols() for x in Newforms(43,2,names='a') ]
            [Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field]
            sage: ModularSymbols(43,2,sign=1).cuspidal_subspace().new_subspace().decomposition()
            [
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field
            ]
        """
        return self.__modsym_space

    def number(self):
        """
        Return the index of this space in the list of simple, new,
        cuspidal subspaces of the full space of modular symbols for
        this weight and level.

        EXAMPLES::

            sage: Newforms(43, 2, names='a')[1].number()
            1
        """
        return self._defining_modular_symbols().ambient().cuspidal_subspace().new_subspace().decomposition().index(self._defining_modular_symbols())

    def __bool__(self):
        """
        Return ``True``, as newforms are never zero.

        EXAMPLES::

            sage: bool(Newforms(14,2)[0])
            True
        """
        return True

    __nonzero__ = __bool__

    def character(self):
        r"""
        The nebentypus character of this newform (as a Dirichlet character with
        values in the field of Hecke eigenvalues of the form).

        EXAMPLES::

            sage: Newforms(Gamma1(7), 4,names='a')[1].character()
            Dirichlet character modulo 7 of conductor 7 mapping 3 |--> 1/2*a1
            sage: chi = DirichletGroup(3).0; Newforms(chi, 7)[0].character() == chi
            True
        """
        return self._defining_modular_symbols().q_eigenform_character(self._name())

    ###########################
    # Atkin--Lehner operators #
    ###########################

    def _atkin_lehner_eigenvalue_from_qexp(self, Q):
        """
        Return the arithmetically-normalized `W_Q`-pseudoeigenvalue of
        ``self``, using a formula based on `q`-expansions (Theorem 2.1 of
        [AL1978]_).

        INPUT:

        - ``self`` -- a newform `f`

        - ``Q`` -- an integer exactly dividing the level of ``self``

        .. NOTE::

            This method assumes that the `Q`-th coefficient in the
            `q`-expansion of ``self`` is non-zero.

        TESTS::

            sage: f = Newforms(Gamma0(18), 4)[0]; f
            q + 2*q^2 + 4*q^4 - 6*q^5 + O(q^6)
            sage: f._atkin_lehner_eigenvalue_from_qexp(2)
            -2
            sage: f._atkin_lehner_eigenvalue_from_qexp(9)
            Traceback (most recent call last):
            ...
            ValueError: a_Q must be nonzero

        An example with odd weight::

            sage: f = Newforms(Gamma1(15), 3, names='a')[2]; f
            q + a2*q^2 + (-a2 - 2)*q^3 - q^4 - a2*q^5 + O(q^6)
            sage: f._atkin_lehner_eigenvalue_from_qexp(5)
            a2
        """
        if Q == 1:
            return ZZ(1)
        a_Q = self[Q]
        if not a_Q:
            raise ValueError("a_Q must be nonzero")

        l = ZZ(1)
        M = self.character().conductor()
        for p, e in Q.factor():
            if p.divides(M): # principal series at p
                l *= (p**(self.weight() - 2) / self[p])**e
            else: # special at p
                l *= -self[p]
        return l

    def _atkin_lehner_eigenvalue_from_modsym(self, Q):
        """
        Return the arithmetically-normalized `W_Q`-pseudoeigenvalue of
        ``self``, using the action of `W_Q` on modular symbols.

        INPUT:

        - ``self`` -- a newform `f`

        - ``Q`` -- a positive integer exactly dividing the level of ``self``

        .. NOTE::

            This algorithm is only valid if, for every prime `p` dividing `Q`,
            the `p`-part of the level of `f` is strictly larger than the
            `p`-part of the conductor of its character.

        EXAMPLES::

            sage: F = Newforms(Gamma1(15), 3, names='a')[2]
            sage: F._atkin_lehner_eigenvalue_from_modsym(5)
            a2
            sage: _ == F._atkin_lehner_eigenvalue_from_qexp(5)
            True
        """
        if Q == 1:
            return ZZ(1)

        S = self._defining_modular_symbols()
        A = S.ambient()
        N = self.level()
        M = self.character().conductor()
        Q0 = M // M.prime_to_m_part(Q)

        W = self.group().atkin_lehner_matrix(Q)
        if Q0 == 1:
            L = [ W ]
        else:
            L = []
            for a in xsrange(Q0):
                if a.gcd(Q0) > 1:
                    continue
                aa = crt(a, 1, Q, N.prime_to_m_part(Q))
                diam = matrix(ZZ, 2, lift_to_sl2z(0,aa,N) )
                L.append( (W * diam * matrix(QQ, 2, [1,a/Q0,0,1]) ).change_ring(ZZ) )

        W = A._matrix_of_operator_on_modular_symbols(A, [x.list() for x in L])
        e = S.dual_eigenvector(names=self._name())
        i = e.nonzero_positions()[0]
        w = (W*e)[i]/e[i]
        if W * e != w * e:
            raise ArithmeticError("Bug in Atkin--Lehner computation: eigenspace not invariant")
        sign = prod([eps(-1) for eps in self.character().decomposition() if eps.conductor().divides(Q)])
        return w / Q0 / sign * self.character()(crt(1, Q//Q0, Q, N//Q))

    def atkin_lehner_action(self, d=None, normalization='analytic', embedding=None):
        r"""
        Return the result of the Atkin-Lehner operator `W_d` on this form `f`,
        in the form of a constant `\lambda_d(f)` and a normalized newform `f'`
        such that

        .. math::

            f \mid W_d = \lambda_d(f) f'.

        See :meth:`atkin_lehner_eigenvalue` for further details.

        EXAMPLES::

            sage: f = Newforms(DirichletGroup(30).1^2, 2, names='a')[0]
            sage: emb = f.base_ring().complex_embeddings()[0]
            sage: for d in divisors(30):
            ....:     print(f.atkin_lehner_action(d, embedding=emb))
            (1.00000000000000, q + a0*q^2 - a0*q^3 - q^4 + (a0 - 2)*q^5 + O(q^6))
            (-1.00000000000000*I, q + a0*q^2 - a0*q^3 - q^4 + (a0 - 2)*q^5 + O(q^6))
            (1.00000000000000*I, q + a0*q^2 - a0*q^3 - q^4 + (a0 - 2)*q^5 + O(q^6))
            (-0.894427190999916 + 0.447213595499958*I, q - a0*q^2 + a0*q^3 - q^4 + (-a0 - 2)*q^5 + O(q^6))
            (1.00000000000000, q + a0*q^2 - a0*q^3 - q^4 + (a0 - 2)*q^5 + O(q^6))
            (-0.447213595499958 - 0.894427190999916*I, q - a0*q^2 + a0*q^3 - q^4 + (-a0 - 2)*q^5 + O(q^6))
            (0.447213595499958 + 0.894427190999916*I, q - a0*q^2 + a0*q^3 - q^4 + (-a0 - 2)*q^5 + O(q^6))
            (-0.894427190999916 + 0.447213595499958*I, q - a0*q^2 + a0*q^3 - q^4 + (-a0 - 2)*q^5 + O(q^6))

        The above computation can also be done exactly::

            sage: K.<z> = CyclotomicField(20)
            sage: f = Newforms(DirichletGroup(30).1^2, 2, names='a')[0]
            sage: emb = f.base_ring().embeddings(CyclotomicField(20, 'z'))[0]
            sage: for d in divisors(30):
            ....:     print(f.atkin_lehner_action(d, embedding=emb))
            (1, q + a0*q^2 - a0*q^3 - q^4 + (a0 - 2)*q^5 + O(q^6))
            (z^5, q + a0*q^2 - a0*q^3 - q^4 + (a0 - 2)*q^5 + O(q^6))
            (-z^5, q + a0*q^2 - a0*q^3 - q^4 + (a0 - 2)*q^5 + O(q^6))
            (-2/5*z^7 + 4/5*z^6 + 1/5*z^5 - 4/5*z^4 - 2/5*z^3 - 2/5, q - a0*q^2 + a0*q^3 - q^4 + (-a0 - 2)*q^5 + O(q^6))
            (1, q + a0*q^2 - a0*q^3 - q^4 + (a0 - 2)*q^5 + O(q^6))
            (4/5*z^7 + 2/5*z^6 - 2/5*z^5 - 2/5*z^4 + 4/5*z^3 - 1/5, q - a0*q^2 + a0*q^3 - q^4 + (-a0 - 2)*q^5 + O(q^6))
            (-4/5*z^7 - 2/5*z^6 + 2/5*z^5 + 2/5*z^4 - 4/5*z^3 + 1/5, q - a0*q^2 + a0*q^3 - q^4 + (-a0 - 2)*q^5 + O(q^6))
            (-2/5*z^7 + 4/5*z^6 + 1/5*z^5 - 4/5*z^4 - 2/5*z^3 - 2/5, q - a0*q^2 + a0*q^3 - q^4 + (-a0 - 2)*q^5 + O(q^6))

        We can compute the eigenvalue of `W_{p^e}` in certain cases
        where the `p`-th coefficient of `f` is zero::

            sage: f = Newforms(169, names='a')[0]; f
            q + a0*q^2 + 2*q^3 + q^4 - a0*q^5 + O(q^6)
            sage: f[13]
            0
            sage: f.atkin_lehner_eigenvalue(169)
            -1

        An example showing the non-multiplicativity of the pseudo-eigenvalues::

            sage: chi = DirichletGroup(18).0^4
            sage: f = Newforms(chi, 2)[0]
            sage: w2, _ = f.atkin_lehner_action(2); w2
            zeta6
            sage: w9, _ = f.atkin_lehner_action(9); w9
            -zeta18^4
            sage: w18,_ = f.atkin_lehner_action(18); w18
            -zeta18
            sage: w18 == w2 * w9 * chi( crt(2, 9, 9, 2) )
            True

        TESTS::

            sage: K.<a> = QuadraticField(1129)
            sage: f = Newforms(Gamma0(20), 8, base_ring=K)[2]; f
            q + (2*a - 10)*q^3 + 125*q^5 + O(q^6)
            sage: f.atkin_lehner_action(2)
            (-1, q + (2*a - 10)*q^3 + 125*q^5 + O(q^6))

            sage: f = Newforms(Gamma1(11), 2)[0]
            sage: f.atkin_lehner_action(2)
            Traceback (most recent call last):
            ...
            ValueError: d (= 2) does not divide the level (= 11)
        """
        # normalise d
        N = self.level()
        if d is None:
            d = N
        else:
            d = ZZ(d)
            if not d.divides(N):
                raise ValueError('d (= {}) does not divide the level (= {})'.format(d, N))
            d = N // N.prime_to_m_part(d)

        if d == 1:
            w = self.base_ring().one()
            if embedding is not None:
                w = embedding(w)
            return w, self

        eps_d = prod([eps.extend(d) for eps in self.character().decomposition() if eps.modulus().divides(d)])
        if eps_d.is_trivial():
            g = self
        else:
            g = self.twist(~eps_d, level=N)
        l = self.atkin_lehner_eigenvalue(d, normalization, embedding)
        return l, g

    def atkin_lehner_eigenvalue(self, d=None, normalization='analytic', embedding=None):
        r"""
        Return the pseudo-eigenvalue of the Atkin-Lehner operator `W_d`
        acting on this form `f`.

        INPUT:

        - ``d`` -- a positive integer exactly dividing the level `N` of `f`,
          i.e. `d` divides `N` and is coprime to `N/d`. The default is `d = N`.

          If `d` does not divide `N` exactly, then it will be replaced with a
          multiple `D` of `d` such that `D` exactly divides `N` and `D` has the
          same prime factors as `d`. An error will be raised if `d` does not
          divide `N`.

        - ``normalization`` -- either ``'analytic'`` (the default) or
          ``'arithmetic'``; see below.

        - ``embedding`` -- (optional) embedding of the coefficient field of `f`
          into another ring. Ignored if `'normalization='arithmetic'``.

        OUTPUT:

        The Atkin-Lehner pseudo-eigenvalue of `W_d` on `f`, as an element of
        the coefficient field of `f`, or the codomain of ``embedding`` if
        specified.

        As defined in [AL1978]_, the pseudo-eigenvalue is the constant
        `\lambda_d(f)` such that

        ..math::

            f \mid W_d = \lambda_d(f) f'

        where `f'` is some normalised newform (not necessarily equal to `f`).

        If ``normalisation='analytic'`` (the default), this routine will
        compute `\lambda_d`, using the conventions of [AL1978]_ for the weight
        `k` action, which imply that `\lambda_d` has complex absolute value 1.
        However, with these conventions `\lambda_d` is not in the Hecke
        eigenvalue field of `f` in general, so it is often necessary to specify
        an embedding of the eigenvalue field into a larger ring (which needs to
        contain roots of unity of sufficiently large order, and a square root
        of `d` if `k` is odd).

        If ``normalisation='arithmetic'`` we compute instead the quotient

        ..math::

            d^{k/2-1} \lambda_d(f) \varepsilon_{N/d}(d / d_0) / G(\varepsilon_d),

        where `G(\varepsilon_d)` is the Gauss sum of the `d`-primary part of
        the nebentype of `f` (more precisely, of its associated primitive
        character), and `d_0` its conductor. This ratio is always in the Hecke
        eigenvalue field of `f` (and can be computed using only arithmetic in
        this field), so specifying an embedding is not needed, although we
        still allow it for consistency.

        (Note that if `k = 2` and `\varepsilon` is trivial, both
        normalisations coincide.)

        .. SEEALSO::

            - :meth:`sage.modular.hecke.module.atkin_lehner_operator`
              (especially for the conventions used to define the operator
              `W_d`)

            - :meth:`atkin_lehner_action`, which returns both the
              pseudo-eigenvalue and the newform `f'`.

        EXAMPLES::

            sage: [x.atkin_lehner_eigenvalue() for x in ModularForms(53).newforms('a')]
            [1, -1]

            sage: f = Newforms(Gamma1(15), 3, names='a')[2]; f
            q + a2*q^2 + (-a2 - 2)*q^3 - q^4 - a2*q^5 + O(q^6)
            sage: f.atkin_lehner_eigenvalue(5)
            Traceback (most recent call last):
            ...
            ValueError: Unable to compute square root. Try specifying an embedding into a larger ring
            sage: L = f.hecke_eigenvalue_field(); x = polygen(QQ); M.<sqrt5> = L.extension(x^2 - 5)
            sage: f.atkin_lehner_eigenvalue(5, embedding=M.coerce_map_from(L))
            1/5*a2*sqrt5
            sage: f.atkin_lehner_eigenvalue(5, normalization='arithmetic')
            a2

            sage: Newforms(DirichletGroup(5).0^2, 6, names='a')[0].atkin_lehner_eigenvalue()
            Traceback (most recent call last):
            ...
            ValueError: Unable to compute Gauss sum. Try specifying an embedding into a larger ring

        TESTS:

        Check that the bug reported at :trac:`18061` is fixed::

            sage: K.<i> = CyclotomicField(4)
            sage: f = Newforms(DirichletGroup(30, QQ).1, 2, K)[0]
            sage: f.atkin_lehner_eigenvalue(embedding=K.embeddings(QQbar)[0])
            -0.8944271909999159? - 0.4472135954999580?*I

        Check that :trac:`24086` is fixed::

            sage: f = Newforms(24, 4)[0]
            sage: f.atkin_lehner_eigenvalue(8)
            -1
            sage: f.atkin_lehner_eigenvalue(3)
            -1

        A case where the eigenvalue isn't in the coefficient field of `f`::

            sage: chi = DirichletGroup(7, QQ).0
            sage: f = Newforms(chi, 3)[0]
            sage: f.atkin_lehner_eigenvalue()
            Traceback (most recent call last):
            ...
            ValueError: Unable to compute square root. Try specifying an embedding into a larger ring
            sage: emb = f.hecke_eigenvalue_field().embeddings(QQbar)[0]
            sage: f.atkin_lehner_eigenvalue(embedding=emb)
            0.?e-18 - 1.000000000000000?*I

        A case where the embeddings really matter::

            sage: chi2 = chi.extend(63)
            sage: g = Newforms(chi2, 3, names='a')[2]
            sage: g.atkin_lehner_eigenvalue(7)
            Traceback (most recent call last):
            ...
            ValueError: Unable to compute Gauss sum. Try specifying an embedding into a larger ring
            sage: g.atkin_lehner_eigenvalue(7, embedding=g.hecke_eigenvalue_field().embeddings(QQbar)[0])
            0.?e-18 + 1.000000000000000?*I
        """
        if normalization not in ['arithmetic', 'analytic']:
            raise ValueError("unknown normalization %s" % normalization)
        N = self.level()
        if d is None:
            d = N
        d = ZZ(d)
        if N % d:
            raise ValueError("d should divide N")
        d = N // N.prime_to_m_part(d)

        d1 = d2 = d3 = 1
        for (p, e) in d.factor():
            if self[p] == 0:
                d1 *= p**e
            elif self.character().conductor().valuation(p) == e:
                d2 *= p**e
            else:
                d3 *= p**e

        verbose("computing W_%s using modsym, W_%s using qexp, W_%s using both" % (d1, d2, d3), level=2)
        w1 = self._atkin_lehner_eigenvalue_from_modsym(d1)
        w2 = self._atkin_lehner_eigenvalue_from_qexp(d2)
        w3a = self._atkin_lehner_eigenvalue_from_modsym(d3)
        w3b = self._atkin_lehner_eigenvalue_from_qexp(d3)
        assert w3a == w3b
        w = w1*w2*w3a

        if embedding is None:
            R = self.hecke_eigenvalue_field()
            embedding = R.Hom(R).identity()
        else:
            R = embedding.codomain()

        if normalization == 'arithmetic':
            return embedding(w)
        else:
            # get rid of the normalisation factors
            from sage.modular.dirichlet import trivial_character

            epsd = prod([eps.extend(d) for eps in self.character().decomposition() if eps.modulus().divides(d)], trivial_character(d))
            epsd = epsd.primitive_character()
            d0 = epsd.modulus()
            epsdR = epsd.change_ring(embedding)
            if d0 > 1:
                try:
                    G = epsdR.gauss_sum()
                except NotImplementedError:
                    raise ValueError("Unable to compute Gauss sum. Try specifying an embedding into a larger ring")
            else:
                G = R(1)
            if not R(d**(self.weight()-2)).is_square():
                raise ValueError("Unable to compute square root. Try specifying an embedding into a larger ring")
            ratio = R(d**(self.weight()-2)).sqrt() * embedding(self.character()( crt(1, d//d0, d, N//d) )) / G
            return embedding(w) / ratio

    def twist(self, chi, level=None, check=True):
        r"""
        Return the twist of the newform ``self`` by the Dirichlet
        character ``chi``.

        If ``self`` is a newform `f` with character `\epsilon` and
        `q`-expansion

        .. MATH::

            f(q) = \sum_{n=1}^\infty a_n q^n,

        then the twist by `\chi` is the unique newform `f\otimes\chi`
        with character `\epsilon\chi^2` and `q`-expansion

        .. MATH::

            (f\otimes\chi)(q) = \sum_{n=1}^\infty b_n q^n

        satisfying `b_n = \chi(n) a_n` for all but finitely many `n`.

        INPUT:

        - ``chi`` -- a Dirichlet character. Note that Sage must be able to
          determine a common base field into which both the Hecke eigenvalue
          field of self, and the field of values of ``chi``, can be embedded.

        - ``level`` -- (optional) the level `N` of the twisted form.
          By default, the algorithm tries to compute `N` using
          [AL1978]_, Theorem 3.1.

        - ``check`` -- (optional) boolean; if ``True`` (default), ensure that
          the space of modular symbols that is computed is genuinely simple and
          new. This makes it less likely that a wrong result is returned if an
          incorrect ``level`` is specified.

        OUTPUT:

        The form `f\otimes\chi` as an element of the set of newforms
        for `\Gamma_1(N)` with character `\epsilon\chi^2`.

        EXAMPLES::

            sage: G = DirichletGroup(3, base_ring=QQ)
            sage: Delta = Newforms(SL2Z, 12)[0]; Delta
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
            sage: Delta.twist(G[0]) == Delta
            True
            sage: Delta.twist(G[1])  # long time (about 5 s)
            q + 24*q^2 - 1472*q^4 - 4830*q^5 + O(q^6)

            sage: M = CuspForms(Gamma1(13), 2)
            sage: f = M.newforms('a')[0]; f
            q + a0*q^2 + (-2*a0 - 4)*q^3 + (-a0 - 1)*q^4 + (2*a0 + 3)*q^5 + O(q^6)
            sage: f.twist(G[1])
            q - a0*q^2 + (-a0 - 1)*q^4 + (-2*a0 - 3)*q^5 + O(q^6)

            sage: f = Newforms(Gamma1(30), 2, names='a')[1]; f
            q + a1*q^2 - a1*q^3 - q^4 + (a1 - 2)*q^5 + O(q^6)
            sage: f.twist(f.character())
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot calculate 5-primary part of the level of the twist of q + a1*q^2 - a1*q^3 - q^4 + (a1 - 2)*q^5 + O(q^6) by Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -1
            sage: f.twist(f.character(), level=30)
            q - a1*q^2 + a1*q^3 - q^4 + (-a1 - 2)*q^5 + O(q^6)

        TESTS:

        We test that feeding inappropriate values of the ``level`` parameter is handled gracefully::

            sage: chi = DirichletGroup(1)[0]
            sage: Delta.twist(chi, level=3)
            Traceback (most recent call last):
            ...
            ValueError: twist of q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6) by Dirichlet character modulo 1 of conductor 1 is not a newform of level 3

        Twisting and twisting back works::

            sage: f = Newforms(11)[0]
            sage: chi = DirichletGroup(5).0
            sage: f.twist(chi).twist(~chi, level=11) == f
            True

        AUTHORS:

        - Peter Bruin (April 2015)

        """
        from sage.modular.all import CuspForms
        R = coercion_model.common_parent(self.base_ring(), chi.base_ring())
        N = self.level()
        epsilon = self.character()
        chi = chi.primitive_character()
        if level is None:
            N_epsilon = epsilon.conductor()
            N_chi = chi.conductor()
            G = DirichletGroup(N_epsilon.lcm(N_chi), base_ring=R)
            epsilon_chi = G(epsilon) * G(chi)
            N_epsilon_chi = epsilon_chi.conductor()
            for q in N_chi.prime_divisors():
                # See [AL1978], Theorem 3.1.
                alpha = N_epsilon.valuation(q)
                beta = N_chi.valuation(q)
                gamma = N.valuation(q)
                delta = max(alpha + beta, 2*beta, gamma)
                if delta == gamma and max(alpha + beta, 2*beta) < gamma:
                    continue
                if delta > gamma and N_epsilon_chi.valuation(q) == max(alpha, beta):
                    continue
                raise NotImplementedError('cannot calculate %s-primary part of the level of the twist of %s by %s'
                                          % (q, self, chi))
            level = lcm([N, N_epsilon * N_chi, N_chi**2])

        # determine the character of the twisted form
        G = DirichletGroup(lcm([N, chi.modulus(), level]), base_ring=R)
        eps_new = (G(epsilon) * G(chi)**2).restrict(level)

        # create an ambient space
        D = ModularSymbols(eps_new, self.weight(), base_ring=R, sign=1).new_submodule()
        S = CuspForms(eps_new, self.weight(), base_ring=R)

        # pull out the eigenspace
        for p in prime_range(500):
            if p.divides(N) or p.divides(chi.level()):
                continue
            D = (D.hecke_operator(p) - self[p]*chi(p)).kernel()
            if D.rank() == 1:
                break
            if D.is_zero():
                raise ValueError('twist of %s by %s is not a newform of level %s' % (self, chi, level))
        else:
            raise RuntimeError('unable to identify modular symbols for twist of %s by %s' % (self, chi))
        return Newform(S, D, names='_', check=check)

class ModularFormElement(ModularForm_abstract, element.HeckeModuleElement):
    def __init__(self, parent, x, check=True):
        r"""
        An element of a space of modular forms.

        INPUT:

        - ``parent`` - ModularForms (an ambient space of modular forms)

        - ``x`` - a vector on the basis for parent

        - ``check`` - if check is ``True``, check the types of the
          inputs.

        OUTPUT:

        - ``ModularFormElement`` - a modular form

        EXAMPLES::

            sage: M = ModularForms(Gamma0(11),2)
            sage: f = M.0
            sage: f.parent()
            Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
        """
        from .space import ModularFormsSpace
        if not isinstance(parent, ModularFormsSpace):
            raise TypeError("First argument must be an ambient space of modular forms.")
        element.HeckeModuleElement.__init__(self, parent, x)

    def _compute_q_expansion(self, prec):
        """
        Computes the q-expansion of self to precision prec.

        EXAMPLES::

            sage: f = EllipticCurve('37a').modular_form()
            sage: f.q_expansion() ## indirect doctest
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)

            sage: f._compute_q_expansion(10)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + 6*q^9 + O(q^10)
        """
        return self.parent()._q_expansion(element = self.element(), prec=prec)

    def _add_(self, other):
        """
        Add self to other.

        EXAMPLES::

            sage: f = ModularForms(DirichletGroup(17).0^2,2).2
            sage: g = ModularForms(DirichletGroup(17).0^2,2).1
            sage: f
            q + (-zeta8^2 + 2)*q^2 + (zeta8 + 3)*q^3 + (-2*zeta8^2 + 3)*q^4 + (-zeta8 + 5)*q^5 + O(q^6)

            sage: g
            1 + (-14/73*zeta8^3 + 57/73*zeta8^2 + 13/73*zeta8 - 6/73)*q^2 + (-90/73*zeta8^3 + 64/73*zeta8^2 - 52/73*zeta8 + 24/73)*q^3 + (-81/73*zeta8^3 + 189/73*zeta8^2 - 3/73*zeta8 + 153/73)*q^4 + (72/73*zeta8^3 + 124/73*zeta8^2 + 100/73*zeta8 + 156/73)*q^5 + O(q^6)

            sage: f+g ## indirect doctest
            1 + q + (-14/73*zeta8^3 - 16/73*zeta8^2 + 13/73*zeta8 + 140/73)*q^2 + (-90/73*zeta8^3 + 64/73*zeta8^2 + 21/73*zeta8 + 243/73)*q^3 + (-81/73*zeta8^3 + 43/73*zeta8^2 - 3/73*zeta8 + 372/73)*q^4 + (72/73*zeta8^3 + 124/73*zeta8^2 + 27/73*zeta8 + 521/73)*q^5 + O(q^6)
        """
        return ModularFormElement(self.parent(), self.element() + other.element())

    def __mul__(self, other):
        r"""
        Calculate the product self * other.

        This tries to determine the
        characters of self and other, in order to avoid having to compute a
        (potentially very large) Gamma1 space. Note that this might lead to
        a modular form that is defined with respect to a larger subgroup than
        the factors are.

        An example with character::

            sage: f = ModularForms(DirichletGroup(3).0, 3).0
            sage: f * f
            1 + 108*q^2 + 144*q^3 + 2916*q^4 + 8640*q^5 + O(q^6)
            sage: (f*f).parent()
            Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(3) of weight 6 over Rational Field
            sage: (f*f*f).parent()
            Modular Forms space of dimension 4, character [-1] and weight 9 over Rational Field

        An example where the character is computed on-the-fly::

            sage: f = ModularForms(Gamma1(3), 5).0
            sage: f*f
            1 - 180*q^2 - 480*q^3 + 8100*q^4 + 35712*q^5 + O(q^6)
            sage: (f*f).parent()
            Modular Forms space of dimension 4 for Congruence Subgroup Gamma0(3) of weight 10 over Rational Field

            sage: f = ModularForms(Gamma1(3), 7).0
            sage: f*f
            q^2 - 54*q^4 + 128*q^5 + O(q^6)
            sage: (f*f).parent()
            Modular Forms space of dimension 5 for Congruence Subgroup Gamma0(3) of weight 14 over Rational Field

        An example with no character::

            sage: f = ModularForms(Gamma1(5), 2).0
            sage: f*f
            1 + 120*q^3 - 240*q^4 + 480*q^5 + O(q^6)
            sage: (f*f).parent()
            Modular Forms space of dimension 5 for Congruence Subgroup Gamma1(5) of weight 4 over Rational Field

        TESTS:

        This shows that the issue at :trac:`7548` is fixed::

            sage: M = CuspForms(Gamma0(5*3^2), 2)
            sage: f = M.basis()[0]
            sage: 2*f
            2*q - 2*q^4 + O(q^6)
            sage: f*2
            2*q - 2*q^4 + O(q^6)
        """

        # boring case: scalar multiplication
        if not isinstance(other, ModularFormElement):
            return element.HeckeModuleElement.__mul__(self, other)

        # first ensure the levels are equal
        if self.level() != other.level():
            raise NotImplementedError("Cannot multiply forms of different levels")

        # find out about characters
        try:
            eps1 = self.character()
            verbose("character of left is %s" % eps1)
            eps2 = other.character()
            verbose("character of right is %s" % eps2)
            newchar = eps1 * eps2
            verbose("character of product is %s" % newchar)
        except (NotImplementedError, ValueError):
            newchar = None
            verbose("character of product not determined")

        # now do the math
        from .constructor import ModularForms
        if newchar is not None:
            verbose("creating a parent with char")
            newparent = ModularForms(newchar, self.weight() + other.weight(),
                                     base_ring=newchar.base_ring())
            verbose("parent is %s" % newparent)
        else:
            newparent = ModularForms(self.group(),
                                     self.weight() + other.weight(),
                                     base_ring=ZZ)
        m = newparent.sturm_bound()
        newqexp = self.qexp(m) * other.qexp(m)

        return newparent.base_extend(newqexp.base_ring())(newqexp)

    modform_lseries = deprecated_function_alias(16917,
            ModularForm_abstract.lseries)

    def atkin_lehner_eigenvalue(self, d=None, embedding=None):
        """
        Return the result of the Atkin-Lehner operator `W_d` on
        ``self``.

        INPUT:

        - ``d`` -- a positive integer exactly dividing the level `N` of
          ``self``, i.e. `d` divides `N` and is coprime to `N/d`. (Default: `d
          = N`)

        - ``embedding`` -- ignored (but accepted for compatibility with
          :meth:`Newform.atkin_lehner_eigenvalue`)

        OUTPUT:

        The Atkin-Lehner eigenvalue of `W_d` on ``self``. If ``self`` is not an
        eigenform for `W_d`, a ``ValueError`` is raised.

        .. SEEALSO::

            For the conventions used to define the operator `W_d`, see
            :meth:`sage.modular.hecke.module.HeckeModule_free_module.atkin_lehner_operator`.

        EXAMPLES::

            sage: CuspForms(1, 30).0.atkin_lehner_eigenvalue()
            1
            sage: CuspForms(2, 8).0.atkin_lehner_eigenvalue()
            Traceback (most recent call last):
            ...
            NotImplementedError: Don't know how to compute Atkin-Lehner matrix acting on this space (try using a newform constructor instead)
        """
        if d is None:
            d = self.level()
        try:
            f = self.parent().atkin_lehner_operator(d)(self)
        except NotImplementedError:
            raise NotImplementedError("Don't know how to compute Atkin-Lehner matrix acting on this space" \
                + " (try using a newform constructor instead)")
        w = self.element().nonzero_positions()[0]
        t = f.element()[w] / self.element()[w]
        if f.element() == self.element() * t:
            return t
        else:
            raise ValueError("%r is not an eigenform for W_%r" % (self, d))

    def twist(self, chi, level=None):
        r"""
        Return the twist of the modular form ``self`` by the Dirichlet
        character ``chi``.

        If ``self`` is a modular form `f` with character `\epsilon`
        and `q`-expansion

        .. MATH::

            f(q) = \sum_{n=0}^\infty a_n q^n,

        then the twist by `\chi` is a modular form `f_\chi` with
        character `\epsilon\chi^2` and `q`-expansion

        .. MATH::

            f_\chi(q) = \sum_{n=0}^\infty \chi(n) a_n q^n.

        INPUT:

        - ``chi`` -- a Dirichlet character

        - ``level`` -- (optional) the level `N` of the twisted form.
          By default, the algorithm chooses some not necessarily
          minimal value for `N` using [AL1978]_, Proposition 3.1,
          (See also [Kob1993]_, Proposition III.3.17, for a simpler
          but slightly weaker bound.)

        OUTPUT:

        The form `f_\chi` as an element of the space of modular forms
        for `\Gamma_1(N)` with character `\epsilon\chi^2`.

        EXAMPLES::

            sage: f = CuspForms(11, 2).0
            sage: f.parent()
            Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
            sage: f.q_expansion(6)
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)
            sage: eps = DirichletGroup(3).0
            sage: eps.parent()
            Group of Dirichlet characters modulo 3 with values in Cyclotomic Field of order 2 and degree 1
            sage: f_eps = f.twist(eps)
            sage: f_eps.parent()
            Cuspidal subspace of dimension 9 of Modular Forms space of dimension 16 for Congruence Subgroup Gamma0(99) of weight 2 over Cyclotomic Field of order 2 and degree 1
            sage: f_eps.q_expansion(6)
            q + 2*q^2 + 2*q^4 - q^5 + O(q^6)

        Modular forms without character are supported::

            sage: M = ModularForms(Gamma1(5), 2)
            sage: f = M.gen(0); f
            1 + 60*q^3 - 120*q^4 + 240*q^5 + O(q^6)
            sage: chi = DirichletGroup(2)[0]
            sage: f.twist(chi)
            60*q^3 + 240*q^5 + O(q^6)

        The base field of the twisted form is extended if necessary::

            sage: E4 = ModularForms(1, 4).gen(0)
            sage: E4.parent()
            Modular Forms space of dimension 1 for Modular Group SL(2,Z) of weight 4 over Rational Field
            sage: chi = DirichletGroup(5)[1]
            sage: chi.base_ring()
            Cyclotomic Field of order 4 and degree 2
            sage: E4_chi = E4.twist(chi)
            sage: E4_chi.parent()
            Modular Forms space of dimension 10, character [-1] and weight 4 over Cyclotomic Field of order 4 and degree 2

        REFERENCES:

        - [AL1978]_

        - [Kob1993]_

        AUTHORS:

        - \L. J. P. Kilford (2009-08-28)

        - Peter Bruin (2015-03-30)

        """
        from sage.modular.all import CuspForms, ModularForms
        from sage.rings.all import PowerSeriesRing
        R = coercion_model.common_parent(self.base_ring(), chi.base_ring())
        N = self.level()
        Q = chi.modulus()
        try:
            epsilon = self.character()
        except ValueError:
            epsilon = None
        constructor = CuspForms if self.is_cuspidal() else ModularForms
        if epsilon is not None:
            if level is None:
                # See [AL1978], Proposition 3.1.
                level = lcm([N, epsilon.conductor() * Q, Q**2])
            G = DirichletGroup(level, base_ring=R)
            M = constructor(G(epsilon) * G(chi)**2, self.weight(), base_ring=R)
        else:
            from sage.modular.arithgroup.all import Gamma1
            if level is None:
                # See [AL1978], Proposition 3.1.
                level = lcm([N, Q]) * Q
            M = constructor(Gamma1(level), self.weight(), base_ring=R)
        bound = M.sturm_bound() + 1
        S = PowerSeriesRing(R, 'q')
        f_twist = S([self[i] * chi(i) for i in range(bound)], prec=bound)
        return M(f_twist)


class ModularFormElement_elliptic_curve(ModularFormElement):
    r"""
    A modular form attached to an elliptic curve over `\QQ`.
    """
    def __init__(self, parent, E):
        """
        Modular form attached to an elliptic curve as an element
        of a space of modular forms.

        EXAMPLES::

            sage: E = EllipticCurve('5077a')
            sage: f = E.modular_form()
            sage: f
            q - 2*q^2 - 3*q^3 + 2*q^4 - 4*q^5 + O(q^6)
            sage: f.q_expansion(10)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 4*q^5 + 6*q^6 - 4*q^7 + 6*q^9 + O(q^10)
            sage: f.parent()
            Modular Forms space of dimension 423 for Congruence Subgroup Gamma0(5077) of weight 2 over Rational Field

            sage: E = EllipticCurve('37a')
            sage: f = E.modular_form() ; f
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)
            sage: f == loads(dumps(f))
            True
        """
        ModularFormElement.__init__(self, parent, None)
##                                    parent.find_in_space( E.q_expansion(parent.hecke_bound()) ))
        self.__E = E


    def elliptic_curve(self):
        """
        Return elliptic curve associated to self.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: f = E.modular_form()
            sage: f.elliptic_curve()
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: f.elliptic_curve() is E
            True
        """
        return self.__E

    def _compute_element(self):
        """
        Compute ``self`` as a linear combination of the basis elements
        of parent.

        EXAMPLES::

            sage: EllipticCurve('11a1').modular_form()._compute_element()
            (1, 0)
            sage: EllipticCurve('389a1').modular_form()._compute_element()
            (1, -2, -2, 2, -3, 4, -5, 0, 1, 6, -4, -4, -3, 10, 6, -4, -6, -2, 5, -6, 10, 8, -4, 0, 4, 6, 4, -10, -6, -12, 4, 8, 0)
        """
        M = self.parent()
        S = M.cuspidal_subspace()
##        return S.find_in_space( self.__E.q_expansion( S.q_expansion_basis()[0].prec() ) ) + [0] * ( M.dimension() - S.dimension() )
        return vector(S.find_in_space(self.__E.q_expansion(S.sturm_bound())) + [0] * (M.dimension() - S.dimension()))

    def _compute_q_expansion(self, prec):
        r"""
        The `q`-expansion of the modular form to precision `O(q^\text{prec})`.

        This function takes one argument, which is the integer ``prec``.

        EXAMPLES::

            sage: E = EllipticCurve('11a') ; f = E.modular_form()
            sage: f._compute_q_expansion(10)
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10)

            sage: f._compute_q_expansion(30)
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 - 2*q^10 + q^11 - 2*q^12 + 4*q^13 + 4*q^14 - q^15 - 4*q^16 - 2*q^17 + 4*q^18 + 2*q^20 + 2*q^21 - 2*q^22 - q^23 - 4*q^25 - 8*q^26 + 5*q^27 - 4*q^28 + O(q^30)

            sage: f._compute_q_expansion(10)
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10)
        """
        return self.__E.q_expansion(prec)

    def atkin_lehner_eigenvalue(self, d=None, embedding=None):
        """
        Return the result of the Atkin-Lehner operator `W_d` on
        ``self``.

        INPUT:

        - ``d`` -- a positive integer exactly dividing the level `N` of
          ``self``, i.e. `d` divides `N` and is coprime to `N/d`. (Defaults to
          `d = N` if not given.)

        - ``embedding`` -- ignored (but accepted for compatibility with
          :meth:`Newform.atkin_lehner_action`)

        OUTPUT:

        The Atkin-Lehner eigenvalue of `W_d` on ``self``. This is either `1` or
        `-1`.

        EXAMPLES::

            sage: EllipticCurve('57a1').newform().atkin_lehner_eigenvalue()
            1
            sage: EllipticCurve('57b1').newform().atkin_lehner_eigenvalue()
            -1
            sage: EllipticCurve('57b1').newform().atkin_lehner_eigenvalue(19)
            1
        """
        if d == self.level():
            w = -self.__E.root_number()
        else:
            # The space of modular symbols attached to E is
            # one-dimensional.
            w = self.__E.modular_symbol_space().atkin_lehner_operator(d).matrix()[0,0]
        return w


class EisensteinSeries(ModularFormElement):
    """
    An Eisenstein series.

    EXAMPLES::

        sage: E = EisensteinForms(1,12)
        sage: E.eisenstein_series()
        [
        691/65520 + q + 2049*q^2 + 177148*q^3 + 4196353*q^4 + 48828126*q^5 + O(q^6)
        ]
        sage: E = EisensteinForms(11,2)
        sage: E.eisenstein_series()
        [
        5/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
        ]
        sage: E = EisensteinForms(Gamma1(7),2)
        sage: E.set_precision(4)
        sage: E.eisenstein_series()
        [
        1/4 + q + 3*q^2 + 4*q^3 + O(q^4),
        1/7*zeta6 - 3/7 + q + (-2*zeta6 + 1)*q^2 + (3*zeta6 - 2)*q^3 + O(q^4),
        q + (-zeta6 + 2)*q^2 + (zeta6 + 2)*q^3 + O(q^4),
        -1/7*zeta6 - 2/7 + q + (2*zeta6 - 1)*q^2 + (-3*zeta6 + 1)*q^3 + O(q^4),
        q + (zeta6 + 1)*q^2 + (-zeta6 + 3)*q^3 + O(q^4)
        ]
    """
    def __init__(self, parent, vector, t, chi, psi):
        """
        An Eisenstein series.

        EXAMPLES::

            sage: E = EisensteinForms(1,12) ## indirect doctest
            sage: E.eisenstein_series()
            [
            691/65520 + q + 2049*q^2 + 177148*q^3 + 4196353*q^4 + 48828126*q^5 + O(q^6)
            ]
            sage: E = EisensteinForms(11,2)
            sage: E.eisenstein_series()
            [
            5/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6)
            ]
            sage: E = EisensteinForms(Gamma1(7),2)
            sage: E.set_precision(4)
            sage: E.eisenstein_series()
            [
            1/4 + q + 3*q^2 + 4*q^3 + O(q^4),
            1/7*zeta6 - 3/7 + q + (-2*zeta6 + 1)*q^2 + (3*zeta6 - 2)*q^3 + O(q^4),
            q + (-zeta6 + 2)*q^2 + (zeta6 + 2)*q^3 + O(q^4),
            -1/7*zeta6 - 2/7 + q + (2*zeta6 - 1)*q^2 + (-3*zeta6 + 1)*q^3 + O(q^4),
            q + (zeta6 + 1)*q^2 + (-zeta6 + 3)*q^3 + O(q^4)
            ]
        """
        N = parent.level()
        K = parent.base_ring()
        if chi.parent().modulus() != N or psi.parent().modulus() != N:
            raise ArithmeticError("Incompatible moduli")
        if chi.parent().base_ring() != K or psi.parent().base_ring() != K:
            raise ArithmeticError("Incompatible base rings")
        t = int(t)
        #if not isinstance(t, int): raise TypeError, "weight must be an int"
        if parent.weight() == 2 and chi.is_trivial() \
               and psi.is_trivial() and t==1:
            raise ArithmeticError("If chi and psi are trivial and k=2, then t must be >1.")
        ModularFormElement.__init__(self, parent, vector)
        self.__chi = chi
        self.__psi = psi
        self.__t   = t

    def _compute_q_expansion(self, prec=None):
        """
        Compute the q-expansion of self to precision prec.

        EXAMPLES::

            sage: EisensteinForms(11,2).eisenstein_series()[0]._compute_q_expansion(10)
            5/12 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + 12*q^6 + 8*q^7 + 15*q^8 + 13*q^9 + O(q^10)
        """
        if prec is None:
            prec = self.parent().prec()
        F = self._compute(list(range(prec)))
        R = self.parent()._q_expansion_ring()
        return R(F, prec)

    def _compute(self, X):
        r"""
        Compute the coefficients of `q^n` of the power series of self,
        for `n` in the list `X`.  The results are not cached.  (Use
        coefficients for cached results).

        EXAMPLES::

            sage: e = DirichletGroup(11).gen()
            sage: f = EisensteinForms(e, 3).eisenstein_series()[0]
            sage: f._compute([3,4,5])
            [-9*zeta10^3 + 1,
             16*zeta10^2 + 4*zeta10 + 1,
             25*zeta10^3 - 25*zeta10^2 + 25*zeta10 - 24]

        """
        if self.weight() == 2 and (self.__chi.is_trivial() and self.__psi.is_trivial()):
            return self.__compute_weight2_trivial_character(X)
        else: # general case
            return self.__compute_general_case(X)

    def __compute_weight2_trivial_character(self, X):
        r"""
        Compute coefficients for self an Eisenstein series of the form
        `E_2 - t*E_2(q^t)`. Computes `a_n` for each `n \in X`.

        EXAMPLES::

            sage: EisensteinForms(14,2).eisenstein_series()[0]._EisensteinSeries__compute_weight2_trivial_character([0])
            [1/24]
            sage: EisensteinForms(14,2).eisenstein_series()[0]._EisensteinSeries__compute_weight2_trivial_character([0,4,11,38])
            [1/24, 1, 12, 20]
        """
        F = self.base_ring()
        v = []
        t = self.__t
        for n in X:
            if n < 0:
                pass
            elif n == 0:
                v.append(F(t - 1) / F(24))
            else:
                an = sigma(n, 1)
                if n % t == 0:
                    an -= t * sigma(n // t, 1)
                v.append(an)
        return v

    def __compute_general_case(self, X):
        r"""
        Return the list coefficients of `q^n` of the power series of self,
        for `n` in the list `X`.  The results are not cached.  (Use
        coefficients for cached results).

        General case (except weight 2, trivial character, where this
        is wrong!)  `\chi` is a primitive character of conductor `L`
        `\psi` is a primitive character of conductor `M` We have
        `MLt \mid N`, and

        .. MATH::

          E_k(chi,psi,t) =
           c_0 + sum_{m \geq 1}[sum_{n|m} psi(n) * chi(m/n) * n^(k-1)] q^{mt},

        with `c_0=0` if `L>1`, and `c_0=L(1-k,psi)/2` if `L=1` (that
        second `L` is an `L`-function `L`).

        EXAMPLES::

            sage: e = DirichletGroup(11).gen()
            sage: f = EisensteinForms(e, 3).eisenstein_series()[0]
            sage: f._EisensteinSeries__compute_general_case([1])
            [1]
            sage: f._EisensteinSeries__compute_general_case([2])
            [4*zeta10 + 1]
            sage: f._EisensteinSeries__compute_general_case([0,1,2])
            [15/11*zeta10^3 - 9/11*zeta10^2 - 26/11*zeta10 - 10/11, 1, 4*zeta10 + 1]
        """
        c0, chi, psi, K, n, t, L, M = self.__defining_parameters()
        zero = K.zero()
        k = self.weight()
        v = []
        for i in X:
            if i == 0:
                v.append(c0)
                continue
            if i % t != 0:
                v.append(zero)
            else:
                m = i // t
                v.append(sum([psi(d) * chi(m / d) * d ** (k - 1)
                              for d in divisors(m)]))
        return v

    @cached_method
    def __defining_parameters(self):
        r"""
        Return defining parameters for ``self``.

        EXAMPLES::

            sage: EisensteinForms(11,2).eisenstein_series()[0]._EisensteinSeries__defining_parameters()
            (-1/24, Dirichlet character modulo 1 of conductor 1, Dirichlet character modulo 1 of conductor 1, Rational Field, 2, 11, 1, 1)
        """
        chi = self.__chi.primitive_character()
        psi = self.__psi.primitive_character()
        k = self.weight()
        t = self.__t
        L = chi.conductor()
        M = psi.conductor()
        K = chi.base_ring()
        n = K.zeta_order()
        if L == 1:
            c0 = K(-psi.bernoulli(k)) / K(2 * k)
        else:
            c0 = K.zero()
        return (c0, chi, psi, K, n, t, L, M)

    def chi(self):
        """
        Return the parameter chi associated to self.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].chi()
            Dirichlet character modulo 17 of conductor 17 mapping 3 |--> zeta16
        """
        return self.__chi

    def psi(self):
        """
        Return the parameter psi associated to self.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].psi()
             Dirichlet character modulo 17 of conductor 1 mapping 3 |--> 1
        """
        return self.__psi

    def t(self):
        """
        Return the parameter t associated to self.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].t()
            1
        """
        return self.__t

    def parameters(self):
        """
        Return chi, psi, and t, which are the defining parameters of self.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].parameters()
            (Dirichlet character modulo 17 of conductor 17 mapping 3 |--> zeta16, Dirichlet character modulo 17 of conductor 1 mapping 3 |--> 1, 1)
        """
        return self.__chi, self.__psi, self.__t

    def L(self):
        """
        Return the conductor of self.chi().

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].L()
            17
        """
        return self.__chi.conductor()

    def M(self):
        """
        Return the conductor of self.psi().

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].M()
            1
        """
        return self.__psi.conductor()

    @cached_method
    def character(self):
        """
        Return the character associated to self.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].character()
            Dirichlet character modulo 17 of conductor 17 mapping 3 |--> zeta16

            sage: chi = DirichletGroup(7)[4]
            sage: E = EisensteinForms(chi).eisenstein_series() ; E
            [
            -1/7*zeta6 - 2/7 + q + (2*zeta6 - 1)*q^2 + (-3*zeta6 + 1)*q^3 + (-2*zeta6 - 1)*q^4 + (5*zeta6 - 4)*q^5 + O(q^6),
            q + (zeta6 + 1)*q^2 + (-zeta6 + 3)*q^3 + (zeta6 + 2)*q^4 + (zeta6 + 4)*q^5 + O(q^6)
            ]
            sage: E[0].character() == chi
            True
            sage: E[1].character() == chi
            True

        TESTS::

            sage: [ [ f.character() == chi for f in EisensteinForms(chi).eisenstein_series() ] for chi in DirichletGroup(17) ]
            [[True], [], [True, True], [], [True, True], [], [True, True], [], [True, True], [], [True, True], [], [True, True], [], [True, True], []]

            sage: [ [ f.character() == chi for f in EisensteinForms(chi).eisenstein_series() ] for chi in DirichletGroup(16) ]
            [[True, True, True, True, True], [], [True, True], [], [True, True, True, True], [], [True, True], []]
        """
        return self.__chi * self.__psi

    def new_level(self):
        """
        Return level at which self is new.

        EXAMPLES::

            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].level()
            17
            sage: EisensteinForms(DirichletGroup(17).0,99).eisenstein_series()[1].new_level()
            17
            sage: [ [x.level(), x.new_level()] for x in EisensteinForms(DirichletGroup(60).0^2,2).eisenstein_series() ]
            [[60, 2], [60, 3], [60, 2], [60, 5], [60, 2], [60, 2], [60, 2], [60, 3], [60, 2], [60, 2], [60, 2]]
        """
        if self.__chi.is_trivial() and self.__psi.is_trivial() and self.weight() == 2:
            return factor(self.__t)[0][0]
        return self.L() * self.M()
