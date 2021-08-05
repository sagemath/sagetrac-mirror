r"""
Graded quasimodular forms ring

Let `E_2` be the weight 2 Eisenstein series defined by

.. MATH::

    E_2(z) = 1 - \frac{2k}{B_k} \sum_{n=1}^{\infty} \sigma(n) q^n

where `\sigma` is the sum of divisors function and `q = \mathrm{exp}(2\pi i z)`
is the classical parameter at infinity, with `\mathrm{im}(z)>0`. This weight 2
Eisenstein series is not a modular forms as it does not statisfies the
modularity condition:

.. MATH::

    z^2 E_2(-1/z) = E_2(z) + \frac{2k}{4\pi i B_k z}.

`E_2` is a quasimodular form of weight 2. General quasimodular forms of given
weight can also be defined. We denotes by `QM` the graded ring of quasimodular
forms for the full modular group `\mathrm{SL}_2(\ZZ)`.

The SageMath implementation of the graded ring of quasimodular forms uses the
following isomorphism:

.. MATH::

    QM \cong M_* [E_2]

where `M_* \cong \CC[E_4, E_6]` is the graded ring of modular forms for
`\mathrm{SL}_2(\ZZ)`. (see :meth:`sage.modular.modform.ring.ModularFormRing`).

EXAMPLES::

    sage: QM = QuasiModularForms(1); QM
    Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field
    sage: QM.gens()
    [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
    1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
    1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
    sage: E2 = QM.0; E4 = QM.1; E6 = QM.2
    sage: E2 * E4 + E6
    2 - 288*q - 20304*q^2 - 185472*q^3 - 855216*q^4 - 2697408*q^5 + O(q^6)
    sage: E2.parent()
    Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field

.. NOTE:

    Only the ring of quasimodular forms for the full modular group have been
    implemented.

AUTHORS:

- DAVID AYOTTE (2021-03-18): initial version

"""

# ****************************************************************************
#       Copyright (C) 2021 DAVID AYOTTE
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modular.arithgroup.all import Gamma0, is_CongruenceSubgroup
from sage.modular.modform.element import GradedModularFormElement, ModularFormElement
from sage.modular.modform.space import ModularFormsSpace
from sage.modular.modform.ring import ModularFormsRing

from sage.rings.all import Integer, QQ, ZZ
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.graded_algebras import GradedAlgebras


from .element import QuasiModularFormsElement

class QuasiModularForms(Parent, UniqueRepresentation):
    r"""
    The graded ring of quasimodular forms for the full modular group
    `{\rm SL}_2(\ZZ)`, with coefficients in a ring.

    EXAMPLES::

        sage: QM = QuasiModularForms(1); QM
        Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field
        sage: QM.gens()
        [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
        1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
        1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]

    It is possible to access the weight 2 Eisenstein series::

        sage: QM.weigt_2_eisenstein_series()
        1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
    """
    Element = QuasiModularFormsElement
    def __init__(self, group=1, base_ring=QQ, name='E2'):
        r"""
        INPUT:

        - ``group`` (default: `{\rm SL}_2(\ZZ)`) -- a congruence subgroup of
          `{\rm SL}_2(\ZZ)`, or a positive integer `N` (interpreted as
          `\Gamma_0(N)`).

        - ``base_ring`` (ring, default: `\QQ`) -- a base ring, which should be
          `\QQ`, `\ZZ`, or the integers mod `p` for some prime `p`.

        - ``name`` (str, default: ``'E2'``) -- a variable name corresponding to
          the weight 2 Eisenstein series.

        TESTS:

            sage: M = QuasiModularForms(1)
            sage: M.group()
            Modular Group SL(2,Z)
            sage: M.base_ring()
            Rational Field
            sage: M = QuasiModularForms(1, Integers(5))
            sage: M.base_ring()
            Ring of integers modulo 5
            sage: QuasiModularForms(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: space of quasimodular forms are only implemented for the full modular group
            sage: QuasiModularForms(Integers(5))
            Traceback (most recent call last):
            ...
            ValueError: Group (=Ring of integers modulo 5) should be a congruence subgroup
            sage: M2 = QuasiModularForms(1, GF(7))
            sage: M == M2
            False

        ::

            sage: TestSuite(QuasiModularForms(1)).run()
        """
        if not isinstance(name, str):
            raise TypeError("`name` must be a string")
        #check if the group is SL2(Z)
        if isinstance(group, (int, Integer)):
            if group>1:
                raise NotImplementedError("space of quasimodular forms are only implemented for the full modular group")
            group = Gamma0(1)
        elif not is_CongruenceSubgroup(group):
            raise ValueError("Group (=%s) should be a congruence subgroup" % group)
        elif group is not Gamma0(1):
            raise NotImplementedError("space of quasimodular forms are implemented for the full modular group")

        #Check if the base ring is a field
        #For some reasons, there is a problem when computing a basis of ModularForms
        if not base_ring.is_field():
            raise ValueError("The base ring must be a field")

        self.__group = group
        self.__base_ring = base_ring
        self.__modular_forms_subring = ModularFormsRing(group, base_ring)
        self.__polynomial_subring = self.__modular_forms_subring[name]
        Parent.__init__(self, base=base_ring, category=GradedAlgebras(base_ring))

    def group(self):
        r"""
        Return the congruence subgroup for which this is the ring of
        quasimodular forms.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.group() is SL2Z
            True
            sage: QM = QuasiModularForms(Gamma0(1)); QM
            Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field

        Higher level congruence subgroups are not yet implemented::

            sage: QuasiModularForms(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: space of quasimodular forms are only implemented for the full modular group
        """
        return self.__group

    def base_ring(self):
        r"""
        Return the coefficient ring of this quasimodular forms ring.

        EXAMPLES::

            sage: QuasiModularForms(1).base_ring()
            Rational Field
            sage: QuasiModularForms(1, base_ring=Integers(5)).base_ring()
            Ring of integers modulo 5
        """
        return self.__base_ring

    def modular_forms_subring(self):
        r"""
        Return the subring of modular forms of this ring of quasimodular forms.

        EXAMPLES::

            sage: QuasiModularForms(1).modular_forms_subring()
            Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field
        """
        return self.__modular_forms_subring

    def modular_forms_of_weight(self, weight):
        r"""
        Return the space of modular forms on this group of the given weight.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.modular_forms_of_weight(12)
            Modular Forms space of dimension 2 for Modular Group SL(2,Z) of weight 12 over Rational Field
        """
        return self.__modular_forms_subring.modular_forms_of_weight(weight)

    def quasimodular_forms_of_weight(self, weight):
        r"""
        Return the space of quasimodular forms on this group of the given weight.

        INPUT:

        - ``weight`` (int, Integer)

        OUTPUT: A quasimodular forms space of the given weight.

        EXAMPLES::

            sage: QuasiModularForms(1).quasimodular_forms_of_weight(4)
            Traceback (most recent call last):
            ...
            NotImplementedError: spaces of quasimodular forms of fixed weight not yet implemented

        """
        raise NotImplementedError("spaces of quasimodular forms of fixed weight not yet implemented")

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: QuasiModularForms(1)._repr_()
            'Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field'
            sage: QuasiModularForms(1, base_ring=Integers(13))._repr_()
            'Ring of Quasimodular Forms for Modular Group SL(2,Z) over Ring of integers modulo 13'
        """
        return "Ring of Quasimodular Forms for %s over %s" % (self.group(), self.base_ring())

    def _coerce_map_from_(self, M):
        r"""
        Code to make QuasiModularForms work well with coercion framework.

        TESTS::

            sage: E2 = QuasiModularForms(1).0
            sage: M = ModularFormsRing(1)
            sage: E2 + M.0
            2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
            sage: M.0 + E2
            2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
            sage: 1 + E2
            2 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: E2 + 1
            2 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: f = ModularForms(1, 12).0
            sage: E2 + f
            1 - 23*q - 96*q^2 + 156*q^3 - 1640*q^4 + 4686*q^5 + O(q^6)
            sage: f + E2
            1 - 23*q - 96*q^2 + 156*q^3 - 1640*q^4 + 4686*q^5 + O(q^6)
        """
        if isinstance(M, (ModularFormsRing, ModularFormsSpace)):
            if M.group() == self.group() and self.has_coerce_map_from(M.base_ring()):
                return True
        if self.base_ring().has_coerce_map_from(M):
            return True
        return False

    def _element_constructor_(self, datum):
        r"""
        The call method of self.

        INPUT:

        - ``datum`` - list, GradedModularFormElement, ModularFormElement,
          Polynomial, base ring element

        OUTPUT: QuasiModularFormElement

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: M = QM.modular_forms_subring()
            sage: m12 = QM.modular_forms_of_weight(12)
            sage: QM([M.0, M.1])
            2 - 288*q - 2448*q^2 + 319104*q^3 + 3681936*q^4 + 21775680*q^5 + O(q^6)
            sage: QM([m12.0, m12.1])
            1 + 49627/691*q + 132611664/691*q^2 + 8380115796/691*q^3 - 13290096200/691*q^4 - 4248043226454/691*q^5 + O(q^6)
            sage: QM([])
            Traceback (most recent call last):
            ...
            ValueError: the given list should be non-empty
            sage: QM(M.0)
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: QM(m12.0)
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
            sage: y = polygen(QQ)
            sage: QM(y)
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: QM(1 + y + y^2)
            3 - 72*q + 360*q^2 + 3168*q^3 + 9288*q^4 + 21456*q^5 + O(q^6)
            sage: QM(1)
            1
            sage: QM(1/2)
            1/2
            sage: QM('E2')
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from <class 'str'> to Univariate Polynomial Ring in E2 over Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field
        """
        if isinstance(datum, list):
            if len(datum) == 0:
                raise ValueError("the given list should be non-empty")
            for idx, f in enumerate(datum):
                if not isinstance(f, (GradedModularFormElement, ModularFormElement)):
                    raise ValueError("one list element is not a modular form")
                datum[idx] = self.__modular_forms_subring(f) #to ensure that every forms is a GradedModularFormElement
            datum = self.__polynomial_subring(datum)
        elif isinstance(datum, (GradedModularFormElement, ModularFormElement)):
            datum = self.__modular_forms_subring(datum) # GradedModularFormElement
            datum = self.__polynomial_subring(datum)
        elif isinstance(datum, Polynomial):
            datum = self.__polynomial_subring(datum.coefficients(sparse=False))
        else:
            datum = self.__polynomial_subring.coerce(datum)
        return self.element_class(self, datum)

    def weigt_2_eisenstein_series(self):
        r"""
        Return the weight 2 Eisenstein series.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: E2 = QM.weigt_2_eisenstein_series(); E2
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: E2.parent()
            Ring of Quasimodular Forms for Modular Group SL(2,Z) over Rational Field
        """
        return self(self.__polynomial_subring.gen())

    def gens(self):
        r"""
        Return a list of generators of the quasimodular forms ring. Note that
        the generators of the modular forms subring is given are the one given
        by the method
        :meth: `sage.modular.modform.ring.ModularFormsRing.gen_forms`

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.gens()
            [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
            sage: QM.modular_forms_subring().gen_forms()
            [1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
        """
        gen_list = [self.weigt_2_eisenstein_series()]
        for f in self.__modular_forms_subring.gen_forms():
            gen_list.append(self(f))
        return gen_list

    generators = gens # alias

    def ngens(self):
        r"""
        Return the number of generators of the given graded quasimodular forms
        ring.

        EXAMPLES::

            sage: QuasiModularForms(1).ngens()
            3
        """
        return len(self.gens())

    def gen(self, n):
        r"""
        Return the `n`-th generator of the quasimodular forms ring.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.0
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: QM.1
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: QM.2
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)
        """
        return self.gens()[n]

    def zero(self):
        r"""
        Return the zero element of this ring.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.zero()
            0
            sage: QM.zero().is_zero()
            True
        """
        return self.element_class(self, self.__polynomial_subring.zero())

    def one(self):
        r"""
        Return the one element of this ring.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.one()
            1
            sage: QM.one().is_one()
            True
        """
        return self.element_class(self, self.__polynomial_subring.one())

    def some_elements(self):
        r"""
        Return a list of generators of ``self``.

        EXAMPLES::

            sage: QuasiModularForms(1).some_elements()
            [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
        """
        return self.gens()

    def polygen(self):
        r"""
        Return the generator of this quasimodular form space as a polynomial
        ring over the modular form subring. Note that this generator correspond
        to the weight-2 Eisenstein series. The default name of this generator is
        'E2'.

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: QM.polygen()
            E2
            sage: QuasiModularForms(1, name='X').polygen()
            X
            sage: QM.polygen().parent()
            Univariate Polynomial Ring in E2 over Ring of Modular Forms for Modular Group SL(2,Z) over Rational Field
        """
        return self.__polynomial_subring.gen()

    def polynomial_ring(self, names='E2, E4, E6'):
        r"""
        Return a multivariate polynomial ring isomorphic to the given graded
        quasimodular forms ring. In the case of the full modular group, this
        ring is `R[E_2, E_4, E_6]` where `E_2`, `E_4` and `E_6` have degrees 2,
        4 and 6 respectively.

        INPUT:

        - ``names`` (str, default: ``'E2, E4, E6'``) -- a list or tuple of names
        (strings), or a comma separated string. Correspond to the names of the
        variables;

        OUTPUT: A multivariate polynomial ring in the variables ``names``

        EXAMPLES:

            sage: QM = QuasiModularForms(1)
            sage: P.<E2, E4, E6> = QM.polynomial_ring(); P
            Multivariate Polynomial Ring in E2, E4, E6 over Rational Field
            sage: E2.degree()
            2
            sage: E4.degree()
            4
            sage: E6.degree()
            6
            sage: P.<x, y, z, w> = QQ[]
            sage: QM.from_polynomial(x+y+z+w)
            Traceback (most recent call last):
            ...
            ValueError: the number of variables (4) of the given polynomial cannot exceed the number of generators (3) of the quasimodular forms ring
        """
        return PolynomialRing(self.base_ring(), 3, names, order=TermOrder('wdeglex', [ZZ(2), ZZ(4), ZZ(6)]))

    def from_polynomial(self, polynomial):
        r"""
        Convert the given polynomial `P(X, Y, Z)` to the graded quasiform
        `P(E_2, E_4, E_6)` where `E_2`, `E_4` and `E_6` are the generators given
        by :meth:`~sage.modular.quasimodform.ring.QuasiModularForms.gens`.

        INPUT:

        - ``plynomial`` -- A multivariate polynomial

        OUTPUT: the graded quasimodular forms `P(E_2, E_4, E_6)`

        EXAMPLES::

            sage: QM = QuasiModularForms(1)
            sage: P.<x, y, z> = QQ[]
            sage: QM.from_polynomial(x)
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: QM.from_polynomial(x) == QM.0
            True
            sage: QM.from_polynomial(y) == QM.1
            True
            sage: QM.from_polynomial(z) == QM.2
            True
            sage: QM.from_polynomial(x^2 + y + x*z + 1)
            4 - 336*q - 2016*q^2 + 322368*q^3 + 3691392*q^4 + 21797280*q^5 + O(q^6)

        TESTS::

            sage: QuasiModularForms(1).from_polynomial('x')
            Traceback (most recent call last):
            ...
            TypeError: the input must be a polynomial
        """
        if not isinstance(polynomial, (MPolynomial, Polynomial)):
            raise TypeError('the input must be a polynomial')
        poly_parent = polynomial.parent()
        nb_var = poly_parent.ngens()
        if nb_var > self.ngens():
            raise ValueError("the number of variables (%s) of the given polynomial cannot exceed the number of generators (%s) of the quasimodular forms ring" % (nb_var, self.ngens()))
        gens_dict = {poly_parent.gen(i):self.gen(i) for i in range(0, nb_var)}
        return self(polynomial.subs(gens_dict))

    def differentiation_operator(self, f):
        r"""
        Compute the formal derivative `q\frac{d}{dq}` of the q-expansion of a
        quasimodular form `f`

        INPUT:

        - ``f`` -- a power serie in corresponding to the q-expansion of a
          quasimodular form.

        OUTPUT:

        The power series `q\frac{d}{dq}(f)`

        EXAMPLES::

            sage: M = QuasiModularForms()
            sage: D = M.differentiation_operator
            sage: B = M.gens()
            sage: P = B[0]; Q = B[1]; R = B[2]
            sage: D(P)
            -24*q - 144*q^2 - 288*q^3 - 672*q^4 - 720*q^5 + O(q^6)
            sage: (P.q_expansion()^2 - Q.q_expansion())/12
            -24*q - 144*q^2 - 288*q^3 - 672*q^4 - 720*q^5 + O(q^6)
            sage: D(Q)
            240*q + 4320*q^2 + 20160*q^3 + 70080*q^4 + 151200*q^5 + O(q^6)
            sage: (P.q_expansion()*Q.q_expansion() - R.q_expansion())/3
            240*q + 4320*q^2 + 20160*q^3 + 70080*q^4 + 151200*q^5 + O(q^6)
            sage: D(R)
            -504*q - 33264*q^2 - 368928*q^3 - 2130912*q^4 - 7877520*q^5 + O(q^6)
            sage: (P.q_expansion()*R.q_expansion() - Q.q_expansion()^2)/2
            -504*q - 33264*q^2 - 368928*q^3 - 2130912*q^4 - 7877520*q^5 + O(q^6)

        TODO:: This method need some work. It should return a QuasiModularFormsElement (not a power series in q)
        """
        q = f.q_expansion().parent().gen()
        return q*f.q_expansion().derivative()
