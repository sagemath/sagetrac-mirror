# -*- coding: utf-8 -*-
r"""
An element in an indexed free module.

AUTHORS:

- Travis Scrimshaw (03-2017): Moved code from :mod:`sage.combinat.free_module`.
"""

#*****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function

from sage.structure.element cimport have_same_parent
from sage.structure.sage_object cimport richcmp, richcmp_not_equal, rich_to_bool
from cpython.object cimport Py_NE, Py_EQ

from sage.misc.misc import repr_lincomb
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.typeset.ascii_art import AsciiArt, empty_ascii_art
from sage.typeset.unicode_art import UnicodeArt, empty_unicode_art
from sage.categories.all import Category, Sets, ModulesWithBasis
from sage.data_structures.blas_dict cimport add, negate, scal, axpy

# TODO: move the content of this class to CombinatorialFreeModule.Element and ModulesWithBasis.Element
cdef class IndexedFreeModuleElement(Element):
    def __init__(self, M, x):
        """
        Create a combinatorial module element. This should never be
        called directly, but only through the parent combinatorial
        free module's :meth:`__call__` method.

        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']; f
            B['a'] + 3*B['c']
            sage: f == loads(dumps(f))
            True
        """
        Element.__init__(self, M)
        self._monomial_coefficients = x
        self._hash_set = False

    def __iter__(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: [i for i in sorted(f)]
            [('a', 1), ('c', 3)]

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1]) + s([3])
            sage: [i for i in sorted(a)]
            [([2, 1], 1), ([3], 1)]
        """
        return self._monomial_coefficients.iteritems()

    def __contains__(self, x):
        """
        Returns whether or not a combinatorial object x indexing a basis
        element is in the support of self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: 'a' in f
            True
            sage: 'b' in f
            False

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1]) + s([3])
            sage: Partition([2,1]) in a
            True
            sage: Partition([1,1,1]) in a
            False
        """
        return x in self._monomial_coefficients and self._monomial_coefficients[x] != 0

    def __hash__(self):
        """
        Return the hash value for ``self``.

        The result is cached.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: hash(f)
            6429418278783588506           # 64-bit
            726440090                     # 32-bit

            sage: F = RootSystem(['A',2]).ambient_space()
            sage: f = F.simple_root(0)
            sage: hash(f)
            6920829894162680369           # 64-bit
            -528971215                    # 32-bit

        This uses the recipe that was proposed for frozendicts in `PEP
        0416 <http://legacy.python.org/dev/peps/pep-0416/>`_ (and adds
        the hash of the parent). This recipe relies on the hash
        function for frozensets which uses tricks to mix the hash
        values of the items in case they are similar.

        .. TODO::

            It would be desirable to make the hash value depend on the
            hash value of the parent. See :trac:`15959`.
        """
        if not self._hash_set:
            self._hash = hash(frozenset(self._monomial_coefficients.items()))
            self._hash_set = True
        return self._hash

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: loads(dumps(F.an_element())) == F.an_element()
            True
        """
        return (_unpickle_element, (self.parent(), self._monomial_coefficients))

    cpdef dict monomial_coefficients(self, bint copy=True):
        """
        Return the internal dictionary which has the combinatorial objects
        indexing the basis as keys and their corresponding coefficients as
        values.

        INPUT:

        - ``copy`` -- (default: ``True``) if ``self`` is internally
          represented by a dictionary ``d``, then make a copy of ``d``;
          if ``False``, then this can cause undesired behavior by
          mutating ``d``

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: d = f.monomial_coefficients()
            sage: d['a']
            1
            sage: d['c']
            3

        To run through the monomials of an element, it is better to
        use the idiom::

            sage: for (t,c) in f:
            ....:     print("{} {}".format(t,c))
            a 1
            c 3

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1])+2*s([3,2])
            sage: d = a.monomial_coefficients()
            sage: type(d)
            <... 'dict'>
            sage: d[ Partition([2,1]) ]
            1
            sage: d[ Partition([3,2]) ]
            2
        """
        if copy:
            return dict(self._monomial_coefficients)
        return self._monomial_coefficients

    def _sorted_items_for_printing(self):
        """
        Returns the items (i.e terms) of ``self``, sorted for printing

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 2*B['c'] + 3 * B['b']
            sage: f._sorted_items_for_printing()
            [('a', 1), ('b', 3), ('c', 2)]
            sage: F.print_options(sorting_reverse=True)
            sage: f._sorted_items_for_printing()
            [('c', 2), ('b', 3), ('a', 1)]
            sage: F.print_options(sorting_reverse=False) #reset to original state

        .. SEEALSO:: :meth:`_repr_`, :meth:`_latex_`, :meth:`print_options`
        """
        print_options = self.parent().print_options()
        v = self._monomial_coefficients.items()
        try:
            v.sort(key=lambda monomial_coeff:
                        print_options['sorting_key'](monomial_coeff[0]),
                   reverse=print_options['sorting_reverse'])
        except Exception: # Sorting the output is a plus, but if we can't, no big deal
            pass
        return v

    def _repr_(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='F')
            sage: e = F.basis()
            sage: e['a'] + 2*e['b'] # indirect doctest
            F['a'] + 2*F['b']
            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='')
            sage: e = F.basis()
            sage: e['a'] + 2*e['b'] # indirect doctest
            ['a'] + 2*['b']
            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='', scalar_mult=' ', bracket=False)
            sage: e = F.basis()
            sage: e['a'] + 2*e['b'] # indirect doctest
            'a' + 2 'b'

        Controling the order of terms by providing a comparison
        function on elements of the support::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'],
            ....:                             sorting_reverse=True)
            sage: e = F.basis()
            sage: e['a'] + 3*e['b'] + 2*e['c']
            2*B['c'] + 3*B['b'] + B['a']

            sage: F = CombinatorialFreeModule(QQ, ['ac', 'ba', 'cb'],
            ....:                             sorting_key=lambda x: x[1])
            sage: e = F.basis()
            sage: e['ac'] + 3*e['ba'] + 2*e['cb']
            3*B['ba'] + 2*B['cb'] + B['ac']
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult=self.parent()._print_options['scalar_mult'],
                            repr_monomial = self.parent()._repr_term,
                            strip_one = True)

    def _ascii_art_(self):
        """
        TESTS::

            sage: M = QuasiSymmetricFunctions(QQ).M()
            sage: ascii_art(M[1,3]**2)  # indirect doctest
            4*M      + 2*M       + 2*M      + 2*M       + 2*M       + M
                 ***      ******        ***         ***         ***     ******
               ***        *             *        ****         ***      **
               *          *           ***        *           **
               *                      *
            sage: ascii_art(M.zero())
            0
        """
        from sage.misc.misc import coeff_repr
        terms = self._sorted_items_for_printing()
        scalar_mult = self.parent()._print_options['scalar_mult']
        repr_monomial = self.parent()._ascii_art_term
        strip_one = True

        if repr_monomial is None:
            repr_monomial = str

        s = empty_ascii_art # ""
        first = True

        if scalar_mult is None:
            scalar_mult = "*"

        for (monomial,c) in terms:
            b = repr_monomial(monomial) # PCR
            if c != 0:
                break_points = []
                coeff = coeff_repr(c, False)
                if coeff != "0":
                    if coeff == "1":
                        coeff = ""
                    elif coeff == "-1":
                        coeff = "-"
                    elif b._l > 0:
                        if len(coeff) > 0 and monomial == 1 and strip_one:
                            b = empty_ascii_art # ""
                        else:
                            b = AsciiArt([scalar_mult]) + b
                    if not first:
                        if len(coeff) > 0 and coeff[0] == "-":
                            coeff = " - %s"%coeff[1:]
                        else:
                            coeff = " + %s"%coeff
                        break_points = [2]
                    else:
                        coeff = "%s"%coeff
                s += AsciiArt([coeff], break_points) + b
                first = False
        if first:
            return AsciiArt(["0"])
        elif s == empty_ascii_art:
            return AsciiArt(["1"])
        else:
            return s

    def _unicode_art_(self):
        """
        TESTS::

            sage: M = QuasiSymmetricFunctions(QQ).M()
            sage: unicode_art(M[1,1]**2)  # indirect doctest
            6*M   + 2*M    + 2*M    + 2*M    + M
               ┌┐      ┌┬┐       ┌┐       ┌┐     ┌┬┐
               ├┤      ├┼┘      ┌┼┤       ├┤    ┌┼┼┘
               ├┤      ├┤       ├┼┘      ┌┼┤    └┴┘
               ├┤      └┘       └┘       └┴┘
               └┘
        """
        from sage.misc.misc import coeff_repr
        terms = self._sorted_items_for_printing()
        scalar_mult = self.parent()._print_options['scalar_mult']
        repr_monomial = self.parent()._unicode_art_term
        strip_one = True

        if repr_monomial is None:
            repr_monomial = str

        s = empty_unicode_art  # ""
        first = True

        if scalar_mult is None:
            scalar_mult = "*"

        for (monomial, c) in terms:
            b = repr_monomial(monomial)  # PCR
            if c != 0:
                break_points = []
                coeff = coeff_repr(c, False)
                if coeff != "0":
                    if coeff == "1":
                        coeff = ""
                    elif coeff == "-1":
                        coeff = "-"
                    elif b._l > 0:
                        if len(coeff) > 0 and monomial == 1 and strip_one:
                            b = empty_unicode_art  # ""
                        else:
                            b = UnicodeArt([scalar_mult]) + b
                    if not first:
                        if len(coeff) > 0 and coeff[0] == "-":
                            coeff = " - %s" % coeff[1:]
                        else:
                            coeff = " + %s" % coeff
                        break_points = [2]
                    else:
                        coeff = "%s" % coeff
                s += UnicodeArt([coeff], break_points) + b
                first = False
        if first:
            return "0"
        elif s == empty_unicode_art:
            return UnicodeArt(["1"])
        else:
            return s

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: latex(f)
            B_{a} + 3B_{c}

        ::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: a = 2 + QS3([2,1,3])
            sage: latex(a) #indirect doctest
            2[1, 2, 3] + [2, 1, 3]

       ::

            sage: F = CombinatorialFreeModule(QQ, ['a','b'], prefix='beta', latex_prefix='\\beta')
            sage: x = F.an_element()
            sage: x
            2*beta['a'] + 2*beta['b']
            sage: latex(x)
            2\beta_{a} + 2\beta_{b}

        Controling the order of terms by providing a comparison
        function on elements of the support::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'],
            ....:                             sorting_reverse=True)
            sage: e = F.basis()
            sage: latex(e['a'] + 3*e['b'] + 2*e['c'])
            2B_{c} + 3B_{b} + B_{a}

            sage: F = CombinatorialFreeModule(QQ, ['ac', 'ba', 'cb'],
            ....:                             sorting_key=lambda x: x[1])
            sage: e = F.basis()
            sage: latex(e['ac'] + 3*e['ba'] + 2*e['cb'])
            3B_{ba} + 2B_{cb} + B_{ac}
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult       = self.parent()._print_options['scalar_mult'],
                            latex_scalar_mult = self.parent()._print_options['latex_scalar_mult'],
                            repr_monomial = self.parent()._latex_term,
                            is_latex=True, strip_one=True)

    cpdef _richcmp_(self, other, int op):
        """
        Rich comparison for equal parents.

        EXAMPLES::

            sage: F1 = CombinatorialFreeModule(QQ, [1, 2, 3])
            sage: F2 = CombinatorialFreeModule(QQ, [1, 2, 3], prefix = "g")
            sage: F1.zero() == F1.zero()
            True
            sage: F1.zero() == F1.an_element()
            False
            sage: F1.an_element() == F1.an_element()
            True
            sage: F1.an_element() is None
            False

        ::

            sage: F3 = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: F3.an_element() != F3.an_element()
            False
            sage: F3.an_element() != F3.zero()
            True

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1])
            sage: b = s([1,1,1])
            sage: cmp(a,b) #indirect doctest
            1

        .. TODO::

            Currently, if ``self`` and ``other`` do not have the same parent,
            seemingly equal elements do not evaluate equal, since conversions
            between different modules have not been established.

        ::

            sage: F1.zero() == 0
            True
            sage: F1(0)
            0

        ::

            sage: F1.zero() == F2.zero()
            False
            sage: F1(F2.zero())
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= 0) an element of self (=Free module generated by {1, 2, 3} over Rational Field)
            sage: F = AlgebrasWithBasis(QQ).example()
            sage: F.one() == 1
            True
            sage: 1 == F.one()
            True
            sage: 2 * F.one() == int(2)
            True
            sage: int(2) == 2 * F.one()
            True

            sage: S = SymmetricFunctions(QQ); s = S.s(); p = S.p()
            sage: p[2] == s[2] - s[1, 1]
            True
            sage: p[2] == s[2]
            False

        This feature is disputable, in particular since it can make
        equality testing costly. It may be removed at some point.

        Equality testing can be a bit tricky when the order of terms
        can vary because their indices are incomparable with
        ``cmp``. The following test did fail before :trac:`12489` ::

            sage: F = CombinatorialFreeModule(QQ, Subsets([1,2,3]))
            sage: x = F.an_element()
            sage: (x+F.zero()).terms()  # random
            [2*B[{1}], 3*B[{2}], B[{}]]
            sage: x.terms()             # random
            [2*B[{1}], B[{}], 3*B[{2}]]
            sage: x+F.zero() == x
            True

        TESTS::

            sage: TestSuite(F1).run()
            sage: TestSuite(F).run()
        """
        if op in [Py_EQ, Py_NE]:
            return richcmp(self._monomial_coefficients, other._monomial_coefficients, op)

        if self._monomial_coefficients == other._monomial_coefficients:
            return rich_to_bool(op, 0)

        v = sorted(self._monomial_coefficients.items())
        w = sorted(other._monomial_coefficients.items())
        return richcmp_not_equal(v, w, op)

    cpdef _add_(self, other):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a'] + 3*B['c']
            B['a'] + 3*B['c']

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: s([2,1]) + s([5,4]) # indirect doctest
            s[2, 1] + s[5, 4]
            sage: a = s([2,1]) + 0
            sage: len(a.monomial_coefficients())
            1
        """
        F = self.parent()
        return F._from_dict(add(self._monomial_coefficients,
                                other._monomial_coefficients),
                            remove_zeros=False)

    cpdef _neg_(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: -f
            -B['a'] - 3*B['c']

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: -s([2,1]) # indirect doctest
            -s[2, 1]
        """
        F = self.parent()
        return F._from_dict(negate(self._monomial_coefficients),
                            remove_zeros=False)

    cpdef _sub_(self, other):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a'] - 3*B['c']
            B['a'] - 3*B['c']

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: s([2,1]) - s([5,4]) # indirect doctest
            s[2, 1] - s[5, 4]
        """
        F = self.parent()
        return F._from_dict(axpy(-1,
                                 other._monomial_coefficients,
                                 self._monomial_coefficients),
                            remove_zeros=False)

    cpdef _coefficient_fast(self, m):
        """
        Return the coefficient of ``m`` in ``self``, where ``m`` is key in
        ``self._monomial_coefficients``.

        EXAMPLES::

            sage: p = Partition([2,1])
            sage: q = Partition([1,1,1])
            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s(p)
            sage: a._coefficient_fast([2,1])
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'list'

        ::

            sage: a._coefficient_fast(p)
            1
            sage: a._coefficient_fast(q)
            0
        """
        return self._monomial_coefficients.get(m, self.base_ring().zero())

    def __getitem__(self, m):
        """
        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: p = Partition([2,1])
            sage: q = Partition([1,1,1])
            sage: a = s(p)
            sage: a[p]
            1
            sage: a[q]
            0
        """
        return self._coefficient_fast(m)

    def _vector_(self, new_base_ring=None):
        """
        Returns ``self`` as a dense vector

        INPUT:

        - ``new_base_ring`` -- a ring (default: ``None``)

        OUTPUT: a dense :func:`FreeModule` vector

        .. WARNING:: This will crash/run forever if ``self`` is infinite dimensional!

        .. SEEALSO::

            - :func:`vector`
            - :meth:`CombinatorialFreeModule.get_order`
            - :meth:`CombinatorialFreeModule.from_vector`
            - :meth:`CombinatorialFreeModule._dense_free_module`

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f._vector_()
            (1, 0, -3)

        One can use equivalently::

            sage: f.to_vector()
            (1, 0, -3)
            sage: vector(f)
            (1, 0, -3)

        More examples::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = 2*QS3([1,2,3]) + 4*QS3([3,2,1])
            sage: a._vector_()
            (2, 0, 0, 0, 0, 4)
            sage: a.to_vector()
            (2, 0, 0, 0, 0, 4)
            sage: vector(a)
            (2, 0, 0, 0, 0, 4)
            sage: a == QS3.from_vector(a.to_vector())
            True

        If ``new_base_ring`` is specified, then a vector over
        ``new_base_ring`` is returned::

            sage: a._vector_(RDF)
            (2.0, 0.0, 0.0, 0.0, 0.0, 4.0)

        .. NOTE::

            :trac:`13406`: the current implementation has been optimized, at
            the price of breaking the encapsulation for FreeModule
            elements creation, with the following use case as metric,
            on a 2008' Macbook Pro::

                sage: F = CombinatorialFreeModule(QQ, range(10))
                sage: f = F.an_element()
                sage: %timeit f._vector_()   # not tested
                625 loops, best of 3: 17.5 micros per loop

             Other use cases may call for different or further
             optimizations.
        """
        parent = self.parent()
        dense_free_module = parent._dense_free_module(new_base_ring)
        d = self._monomial_coefficients
        zero = dense_free_module.base_ring().zero()
        return dense_free_module.element_class(dense_free_module,
                                               [d.get(m, zero) for m in parent.get_order()],
                                               coerce=True, copy=False)

    to_vector = _vector_

    cpdef _acted_upon_(self, scalar, bint self_on_left):
        """
        Return the action of ``scalar`` (an element of the base ring) on
        ``self``.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a']*(1/2)  # indirect doctest
            1/2*B['a']
            sage: B['a']/2
            1/2*B['a']
            sage: B['a']*2      # indirect doctest
            2*B['a']
            sage: B['a']*int(2) # indirect doctest
            2*B['a']

            sage: 1/2*B['a']
            1/2*B['a']
            sage: 2*B['a']      # indirect doctest
            2*B['a']
            sage: int(2)*B['a'] # indirect doctest
            2*B['a']

        TESTS::

            sage: F.get_action(QQ, operator.mul, True)
            Right action by Rational Field on Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: F.get_action(QQ, operator.mul, False)
            Left action by Rational Field on Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: F.get_action(ZZ, operator.mul, True)
            Right action by Integer Ring on Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: F.get_action(F, operator.mul, True)
            sage: F.get_action(F, operator.mul, False)

        This also works when a coercion of the coefficient is needed, for
        example with polynomials or fraction fields (:trac:`8832`)::

            sage: P.<q> = QQ['q']
            sage: V = CombinatorialFreeModule(P, Permutations())
            sage: el = V(Permutation([3,1,2]))
            sage: (3/2)*el
            3/2*B[[3, 1, 2]]

            sage: P.<q> = QQ['q']
            sage: F = FractionField(P)
            sage: V = CombinatorialFreeModule(F, Words())
            sage: w = Words()('abc')
            sage: (1+q)*V(w)
            (q+1)*B[word: abc]
            sage: ((1+q)/q)*V(w)
            ((q+1)/q)*B[word: abc]

        .. TODO::

            Add non commutative tests.
        """
        # With the current design, the coercion model does not have
        # enough information to detect a priori that this method only
        # accepts scalars; so it tries on some elements(), and we need
        # to make sure to report an error.
        if isinstance(scalar, Element) and scalar.parent() is not self.base_ring():
            # Temporary needed by coercion (see Polynomial/FractionField tests).
            if self.base_ring().has_coerce_map_from(scalar.parent()):
                scalar = self.base_ring()( scalar )
            else:
                return None

        F = self.parent()
        D = self._monomial_coefficients
        return F._from_dict(scal(scalar, D, factor_on_left=not self_on_left),
                            remove_zeros=False)

    def _lmul_(self, other):
        """
        For backward compatibility.

        EXAMPLES::

            sage: C = CombinatorialFreeModule(QQ, [1,2,3])
            sage: C.an_element()._lmul_(2)
            4*B[1] + 4*B[2] + 6*B[3]
        """
        return self._acted_upon_(other, True)

    def _rmul_(self, other):
        """
        For backward compatibility.

        EXAMPLES::

            sage: C = CombinatorialFreeModule(QQ, [1,2,3])
            sage: C.an_element()._rmul_(2)
            4*B[1] + 4*B[2] + 6*B[3]
        """
        return self._acted_upon_(other, False)

    def __truediv__(self, x):
        """
        Division by coefficients.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, [1,2,3])
            sage: x = F._from_dict({1:2, 2:3})
            sage: x/2
            B[1] + 3/2*B[2]

        ::

            sage: F = CombinatorialFreeModule(QQ, [1,2,3])
            sage: B = F.basis()
            sage: f = 2*B[2] + 4*B[3]
            sage: f/2
            B[2] + 2*B[3]
        """
        if not self.base_ring().is_field():
            return self.map_coefficients(lambda c: _divide_if_possible(c, x))

        F = self.parent()
        x = self.base_ring()( x )
        x_inv = x ** -1
        D = self._monomial_coefficients
        return F._from_dict(scal(x_inv, D),
                            remove_zeros=False)

    __div__ = __truediv__


def _divide_if_possible(x, y):
    """
    EXAMPLES::

        sage: from sage.modules.with_basis.indexed_free_module_element import _divide_if_possible
        sage: _divide_if_possible(4, 2)
        2
        sage: _.parent()
        Integer Ring

    ::

        sage: _divide_if_possible(4, 3)
        Traceback (most recent call last):
        ...
        ValueError: 4 is not divisible by 3
    """
    q, r = x.quo_rem(y)
    if r != 0:
        raise ValueError("%s is not divisible by %s"%(x, y))
    else:
        return q

def _unpickle_element(C, d):
    """
    Unpickle an element in ``C`` given by ``d``.

    EXAMPLES::

        sage: from sage.modules.with_basis.indexed_free_module_element import _unpickle_element
        sage: C = CombinatorialFreeModule(QQ, [1,2,3])
        sage: _unpickle_element(C, {1: -2, 3: -12})
        -2*B[1] - 12*B[3]
    """
    return C._from_dict(d, coerce=False, remove_zeros=False)

