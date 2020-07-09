r"""
Univariate Skew Polynomials

This module provides the
:class:`~sage.rings.polynomial.skew_polynomial_element.SkewPolynomial`.
In the class hierarchy in Sage, the locution *Skew Polynomial* is used 
for a Ore polynomial without twisting derivation.

.. WARNING::

    The current semantics of
    :meth:`~sage.rings.polynomial.skew_polynomial_element.SkewPolynomial.__call__`
    are experimental, so a warning is thrown when a skew polynomial is evaluated
    for the first time in a session. See the method documentation for details.

    TESTS::

        sage: R.<t> = QQ[]
        sage: sigma = R.hom([t+1])
        sage: S.<x> = R['x',sigma]
        sage: a = 2*(t + x) + 1
        sage: a(t^2)
        doctest:...: FutureWarning: This class/method/function is marked as experimental. 
        It, its functionality or its interface might change without a formal deprecation.
        See http://trac.sagemath.org/13215 for details.
        2*t^3 + 3*t^2 + 4*t + 2
        sage: a(t)
        2*t^2 + 3*t + 2

AUTHORS:

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests and
  refactored classes and methods

- Johan Rosenkilde (2016-08-03): changes for bug fixes, docstring and
  doctest errors

"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

import re
from cysignals.signals cimport sig_check

from sage.rings.infinity import infinity
from sage.structure.factorization import Factorization
from sage.structure.element cimport Element, RingElement, AlgebraElement, ModuleElement
from sage.structure.parent cimport Parent
from sage.structure.parent_gens cimport ParentWithGens
from sage.misc.abstract_method import abstract_method
from sage.categories.homset import Hom
from sage.categories.fields import Fields
from sage.rings.integer cimport Integer
from cpython.object cimport PyObject_RichCompare
from sage.categories.map cimport Map
from sage.rings.morphism cimport Morphism, RingHomomorphism
from sage.rings.polynomial.polynomial_element cimport _dict_to_list
from sage.structure.element import coerce_binop
from sage.misc.superseded import experimental

from sage.rings.polynomial.ore_polynomial_element cimport OrePolynomial
from sage.rings.polynomial.ore_polynomial_element cimport OrePolynomial_generic_dense


cdef class SkewPolynomial_generic_dense(OrePolynomial_generic_dense):
    r"""
    Generic implementation of dense skew polynomial supporting any valid base
    ring and twisting morphism.
    """
    cpdef left_power_mod(self, exp, modulus):
        r"""
        Return the remainder of ``self**exp`` in the left euclidean division
        by ``modulus``.

        INPUT:

        - ``exp`` -- an Integer

        - ``modulus`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        Remainder of ``self**exp`` in the left euclidean division
        by ``modulus``.

        REMARK:

        The quotient of the underlying skew polynomial ring by the
        principal ideal generated by ``modulus`` is in general *not*
        a ring.

        As a consequence, Sage first computes exactly ``self**exp``
        and then reduce it modulo ``modulus``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: a.left_power_mod(100,modulus)
            (4*t^2 + t + 1)*x^2 + (t^2 + 4*t + 1)*x + 3*t^2 + 3*t
        """
        cdef SkewPolynomial_generic_dense r
        if not isinstance(exp, Integer):
            try:
                exp = Integer(exp)
            except TypeError:
                raise TypeError("non-integral exponents not supported")

        if len(self._coeffs) <= 1:
            return self.parent()(self._coeffs[0]**exp)
        if exp == 0:
            return self.parent().one()
        if exp < 0:
            return (~self).left_power_mod(-exp, modulus)

        if self == self.parent().gen():
            P = self.parent()
            R = P.base_ring()
            v = [R.zero()]*exp + [R.one()]
            r = <OrePolynomial_generic_dense>self._parent(v)
        else:
            r = <OrePolynomial_generic_dense>self._new_c(list(self._coeffs), self._parent)
            r._inplace_pow(exp)

        if modulus:
            _, r = r._left_quo_rem(modulus)
        return r

    cpdef right_power_mod(self, exp, modulus):
        r"""
        Return the remainder of ``self**exp`` in the right euclidean division
        by ``modulus``.

        INPUT:

        - ``exp`` -- an Integer

        - ``modulus`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        Remainder of ``self**exp`` in the right euclidean division
        by ``modulus``.

        REMARK:

        The quotient of the underlying skew polynomial ring by the
        principal ideal generated by ``modulus`` is in general *not*
        a ring.

        As a consequence, Sage first computes exactly ``self**exp``
        and then reduce it modulo ``modulus``.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: b = a^10  # short form for ``a._pow_(10)``
            sage: b == a*a*a*a*a*a*a*a*a*a
            True
            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: br = a.right_power_mod(10,modulus); br
            (t^2 + t)*x^2 + (3*t^2 + 1)*x + t^2 + t
            sage: rq, rr = b.right_quo_rem(modulus)
            sage: br == rr
            True
            sage: a.right_power_mod(100,modulus)
            (2*t^2 + 3)*x^2 + (t^2 + 4*t + 2)*x + t^2 + 2*t + 1
        """
        cdef SkewPolynomial_generic_dense r
        if not isinstance(exp, Integer):
            try:
                exp = Integer(exp)
            except TypeError:
                raise TypeError("non-integral exponents not supported")

        if len(self._coeffs) <= 1:
            return self.parent()(self._coeffs[0]**exp)
        if exp == 0:
            return self.parent().one()
        if exp < 0:
            return (~self).right_power_mod(-exp, modulus)

        if self == self.parent().gen():
            P = self.parent()
            R = P.base_ring()
            v = [R.zero()]*exp + [R.one()]
            r = <OrePolynomial_generic_dense>self._parent(v)
        else:
            r = <OrePolynomial_generic_dense>self._new_c(list(self._coeffs), self._parent)
            r._inplace_pow(exp)

        if modulus:
            _, r = r._right_quo_rem(modulus)
        return r

    def __pow__(self, exp, modulus):
        r"""
        Return the remainder of ``self**exp`` in the left euclidean
        division by ``modulus``.

        INPUT:

        - ``exp`` -- an Integer

        - ``modulus`` -- a skew polynomial in the same ring as ``self``

        OUTPUT:

        Remainder of ``self**exp`` in the right euclidean division
        by ``modulus``.

        REMARK:

        The quotient of the underlying skew polynomial ring by the
        principal ideal generated by ``modulus`` is in general *not*
        a ring.

        As a consequence, Sage first computes exactly ``self**exp``
        and then reduce it modulo ``modulus``.

        .. SEEALSO::

            :meth:`~sage.rings.polynomial.skew_polynomial_element._pow_`

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: b = a^10
            sage: b == a*a*a*a*a*a*a*a*a*a
            True
            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: bmod = a.right_power_mod(10,modulus); bmod
            (t^2 + t)*x^2 + (3*t^2 + 1)*x + t^2 + t
            sage: rq, rr = b.right_quo_rem(modulus)
            sage: bmod == rr
            True
        """
        return self.right_power_mod(exp, modulus)

    def __call__(self, eval_pt):
        r"""
        Evaluate ``self`` at ``eval_pt`` using operator evaluation.

        Given a skew polynomial `p(x) = \sum_{i=0}^d a_i * x^i`, we define
        the evaluation `p(r)` to be `\sum_{i=0}^d a_i * \sigma^i(r)`, where
        `\sigma` is the twisting morphism of the skew polynomial ring.

        INPUT:

        - ``eval_pt`` -- element of the base ring of ``self``

        OUTPUT:

        The operator evaluation of ``self`` at ``eval_pt``.

        .. TODO::

            Currently, only "operator evaluation" of skew polynomials is
            implemented (see :meth:`.operator_eval`).
            There are two other notions of evaluation of a skew polynomial
            `p(x)` at some element `a` of the base ring. First, the value
            of the polynomial can be defined as the remainder of the right
            division of `p(x)` by `x-a`. Second, the value can be given by
            the formula, `p(a) = \sum_{i=0}^{m-1} B_{i} * p(\beta_{i})`
            where `m` is the degree of the base ring (`F_{q^m}`) of the skew
            polynomial ring, `B_{i}` is the `i`-th element in the vector
            representation of `a` in `F_{q}` and`\beta_{i}` is the `i`-th
            element of the corresponding basis of `F_{q^m}` over `F_{q}`.
            
            The current calling convention might change in the future to
            accommodate these. Therefore, the current method has been
            marked as experimental.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x + 1
            sage: a(t^2)
            t^3 + 3*t^2 + t
            sage: b = x^2 + t*x^3 + t^2*x + 1
            sage: b(2*t + 3)
            2*t^3 + 7*t^2 + 13*t + 10
        """
        return self._call(eval_pt)

    @experimental(trac_number=13215)
    def _call(self, eval_pt):
        r"""
        Helper function for the :meth:`__call__` method to accommodate
        the ``@experimental`` decorator.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<x> = k['x',Frob]
            sage: a = 3*t^2*x^2 + (t + 1)*x + 2
            sage: a(t) #indirect test
            2*t^2 + 2*t + 3
        """
        return self.operator_eval(eval_pt)

    cpdef operator_eval(self, eval_pt):
        r"""
        Evaluate ``self`` at ``eval_pt`` by the operator evaluation
        method.

        INPUT:

        - ``eval_pt`` -- element of the base ring of ``self``

        OUTPUT:

        The value of the polynomial at the point specified by the argument.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<x> = k['x',Frob]
            sage: a = 3*t^2*x^2 + (t + 1)*x + 2
            sage: a(t) #indirect test
            2*t^2 + 2*t + 3
            sage: a.operator_eval(t)
            2*t^2 + 2*t + 3

        Evaluation points outside the base ring is usually not possible due to the twisting morphism::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t*x + 1
            sage: a.operator_eval(1/t)
            Traceback (most recent call last):
            ...
            TypeError: 1/t fails to convert into the map's domain Univariate Polynomial Ring in t over Rational Field, but a `pushforward` method is not properly implemented
        """
        cdef RingHomomorphism sigma = self._parent.twisting_morphism()
        cdef list coefficients = self.list()
        cdef RingElement ret = self.base_ring().zero()
        cdef RingElement a = eval_pt
        for c in coefficients:
            ret += c * a
            a = sigma(a)
        return ret

    def conjugate(self, n):
        r"""
        Return ``self`` conjugated by `x^n`, where `x` is the
        variable of ``self``.

        The conjugate is obtained from ``self`` by applying the `n`-th iterate
        of the twisting morphism to each of its coefficients.

        INPUT:

        - `n` -- an integer, the power of conjugation

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: K = R.fraction_field()
            sage: sigma = K.hom([1 + 1/t])
            sage: S.<x> = K['x',sigma]
            sage: a = t*x^3 + (t^2 + 1)*x^2 + 2*t
            sage: b = a.conjugate(2); b
            ((2*t + 1)/(t + 1))*x^3 + ((5*t^2 + 6*t + 2)/(t^2 + 2*t + 1))*x^2 + (4*t + 2)/(t + 1)
            sage: x^2*a == b*x^2
            True

        In principle, negative values for `n` are allowed, but Sage needs to be
        able to invert the twisting morphism::

            sage: b = a.conjugate(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: inverse not implemented for morphisms of Fraction Field of Univariate Polynomial Ring in t over Rational Field

        Here is a working example::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: T.<y> = k['y',Frob]
            sage: u = T.random_element(); u
            (2*t^2 + 3)*y^2 + (4*t^2 + t + 4)*y + 2*t^2 + 2
            sage: v = u.conjugate(-1); v
            (3*t^2 + t)*y^2 + (4*t^2 + 2*t + 4)*y + 3*t^2 + t + 4
            sage: u*y == y*v
            True
        """
        r = self._new_c([self._parent.twisting_morphism(n)(x) for x in self.list()],
                        self._parent, 0)
        return r

    def multi_point_evaluation(self, eval_pts):
        """
        Evaluate ``self`` at list of evaluation points.

        INPUT:

        - ``eval_pts`` -- list of points at which ``self`` is to be evaluated

        OUTPUT:

        List of values of ``self`` at the ``eval_pts``.

        .. TODO::

            This method currently trivially calls the evaluation function
            repeatedly. If fast skew polynomial multiplication is available, an
            asymptotically faster method is possible using standard divide and
            conquer techniques and
            :meth:`sage.rings.polynomial.skew_polynomial_ring.SkewPolynomialRing_general.minimal_vanishing_polynomial`.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: eval_pts = [1, t, t^2]
            sage: c = a.multi_point_evaluation(eval_pts); c
            [t + 1, 3*t^2 + 4*t + 4, 4*t]
            sage: c == [ a(e) for e in eval_pts ]
            True
        """
        return [ self(e) for e in eval_pts ]

    cpdef ModuleElement _lmul_(self, Element right):
        r"""
        Return the product ``self * right``.

        INPUT:

        - ``right`` -- an element of the base ring

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x + t
            sage: b = t
            sage: a * b
            (t + 1)*x + t^2
            sage: a * b == b * a
            False
        """
        if right == 0:
            return self._parent.zero()
        cdef list x = (<SkewPolynomial_generic_dense>self)._coeffs
        cdef Py_ssize_t i
        twisting_morphism = self._parent._morphism
        r = self._new_c([ (twisting_morphism**i)(right)*x[i] for i from 0 <= i < len(x) ],
                        self._parent, 0)
        return r

    cpdef ModuleElement _rmul_(self, Element left):
        r"""
        Return the product ``left * self``.

        INPUT:

        - ``left`` -- an element of the base ring

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = t
            sage: b = x + t
            sage: a * b
            t*x + t^2
            sage: a * b == b * a
            False
        """
        if left == 0:
            return self.parent().zero()
        cdef list x = (<SkewPolynomial_generic_dense>self)._coeffs
        cdef Py_ssize_t i
        r = self._new_c([ left*x[i] for i from 0 <= i < len(x) ], self._parent, 0)
        return r

    cpdef _mul_(self, right):
        r"""
        Return the product ``self * right``.

        INPUT:

        - ``right`` -- a Ore polynomial in the same ring as ``self``

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: a = x^2 + t; a
            x^2 + t
            sage: b = x^2 + (t + 1)*x; b
            x^2 + (t + 1)*x
            sage: a * b
            x^4 + (t + 3)*x^3 + t*x^2 + (t^2 + t)*x
            sage: a * b == b * a
            False

        TESTS::

            sage: S(0)*a, (S(0)*a).list()
            (0, [])
        """
        cdef list x = (<SkewPolynomial_generic_dense>self)._coeffs
        cdef list y = (<SkewPolynomial_generic_dense>right)._coeffs
        cdef Py_ssize_t i, k, start, end
        cdef Py_ssize_t dx = len(x)-1, dy = len(y)-1
        parent = self._parent
        if dx == -1:
            return self # = zero
        elif dy == -1:
            return right # = zero
        elif dx == 0:
            c = x[0]
            r = self._new_c([c*a for a in y], parent, 0)
            return r
        cdef list coeffs = []
        for k from 0 <= k <= dx+dy:
            start = 0 if k <= dy else k-dy
            end = k if k <= dx else dx
            sum = x[start] * parent.twisting_morphism(start)(y[k-start])
            for i from start < i <= end:
                sum += x[i] * parent.twisting_morphism(i)(y[k-i])
            coeffs.append(sum)
        r = self._new_c(coeffs, parent, 0)
        return r

    cdef void _inplace_rmul(self, SkewPolynomial_generic_dense right):
        r"""
        Replace ``self`` by ``self*right`` (only for internal use).

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: a.left_power_mod(100,modulus)  # indirect doctest
            (4*t^2 + t + 1)*x^2 + (t^2 + 4*t + 1)*x + 3*t^2 + 3*t
        """
        cdef list x = self._coeffs
        cdef list y = right._coeffs
        cdef Py_ssize_t i, k, start, end
        cdef Py_ssize_t d1 = len(x)-1, d2 = len(y)-1
        parent = self._parent
        if d2 == -1:
            self._coeffs = [ ]
        elif d1 >= 0:
            for k from d1 < k <= d1+d2:
                start = 0 if k <= d2 else k-d2
                sum = x[start] * parent.twisting_morphism(start)(y[k-start])
                for i from start < i <= d1:
                    sum += x[i] * parent.twisting_morphism(i)(y[k-i])
                x.append(sum)
            for k from d1 >= k >= 0:
                start = 0 if k <= d2 else k-d2
                end = k if k <= d1 else d1
                sum = x[start] * parent.twisting_morphism(start)(y[k-start])
                for i from start < i <= end:
                    sum += x[i] * parent.twisting_morphism(i)(y[k-i])
                x[k] = sum

    cdef void _inplace_pow(self, Py_ssize_t n):
        r"""
        Replace ``self`` by ``self**n`` (only for internal use).

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x + t
            sage: modulus = x^3 + t*x^2 + (t+3)*x - 2
            sage: a.left_power_mod(100,modulus)  # indirect doctest
            (4*t^2 + t + 1)*x^2 + (t^2 + 4*t + 1)*x + 3*t^2 + 3*t
        """
        while n & 1 == 0:
            self._inplace_rmul(self)
            n = n >> 1
        cdef SkewPolynomial_generic_dense selfpow = <SkewPolynomial_generic_dense>self._new_c(list(self._coeffs), self._parent)
        n = n >> 1
        while n != 0:
            selfpow._inplace_rmul(selfpow)
            if n&1 == 1:
                self._inplace_rmul(selfpow)
            n = n >> 1

    cdef _left_quo_rem(self, OrePolynomial other):
        r"""
        Return the quotient and remainder of the left euclidean
        division of ``self`` by ``other`` (C implementation).
        """
        sig_check()
        cdef list a = list(self._coeffs)
        cdef list b = (<SkewPolynomial_generic_dense>other)._coeffs
        cdef Py_ssize_t i, j
        cdef Py_ssize_t da = self.degree(), db = other.degree()
        if da < db:
            return (self._new_c([], self._parent), self)
        try:
            inv = self.base_ring()(~b[db])
        except (ZeroDivisionError, TypeError):
            raise NotImplementedError("the leading coefficient of the divisor is not invertible")
        cdef list q = [ ]
        parent = self._parent
        for i from da-db >= i >= 0:
            try:
                c = parent.twisting_morphism(-db)(inv*a[i+db])
                for j from 0 <= j < db:
                    a[i+j] -= b[j] * parent.twisting_morphism(j)(c)
            except Exception:
                raise NotImplementedError("inversion of the twisting morphism %s" % parent.twisting_morphism())
            q.append(c)
        q.reverse()
        return (self._new_c(q, parent), self._new_c(a[:db], parent, 1))

    cdef _right_quo_rem(self, OrePolynomial other):
        r"""
        Return the quotient and remainder of the right euclidean
        division of ``self`` by ``other`` (C implementation).
        """
        sig_check()
        cdef list a = list(self._coeffs)
        cdef list b = (<SkewPolynomial_generic_dense>other)._coeffs
        cdef Py_ssize_t i, j
        cdef Py_ssize_t da = self.degree(), db = other.degree()
        parent = self._parent
        if da < db:
            return (self._new_c([],parent), self)
        try:
            inv = self.base_ring()(~b[db])
        except (ZeroDivisionError, TypeError):
            raise NotImplementedError("the leading coefficient of the divisor"
                                      " is not invertible")
        cdef list q = [ ]
        parent = self._parent
        for i from da-db >= i >= 0:
            c = parent.twisting_morphism(i)(inv) * a[i+db]
            for j from 0 <= j < db:
                a[i+j] -= c * parent.twisting_morphism(i)(b[j])
            q.append(c)
        q.reverse()
        return (self._new_c(q, parent), self._new_c(a[:db], parent, 1))
