"""
Free modules
"""
#*****************************************************************************
#       Copyright (C) 2007      Mike Hansen <mhansen@gmail.com>,
#                     2007-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#                     2010      Christian Stump <christian.stump@univie.ac.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element, have_same_parent
from sage.structure.parent import Parent
from sage.structure.indexed_generators import IndexedGenerators
from sage.modules.module import Module
from sage.misc.misc import repr_lincomb
from sage.misc.cachefunc import cached_method
from sage.misc.all import lazy_attribute
from sage.rings.all import Integer
import sage.structure.element
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.categories.poor_man_map import PoorManMap
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.category import Category, JoinCategory
from sage.categories.tensor import tensor, TensorProductsCategory
from sage.combinat.cartesian_product import CartesianProduct
from sage.categories.all import Sets
from sage.combinat.dict_addition import dict_addition, dict_linear_combination
from sage.sets.family import Family, FiniteFamily
from sage.typeset.ascii_art import AsciiArt, empty_ascii_art

# TODO: move the content of this class to CombinatorialFreeModule.Element and ModulesWithBasis.Element

class CombinatorialFreeModuleElement(Element):
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

    @cached_method
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
        return hash(frozenset(self._monomial_coefficients.items()))

    def monomial_coefficients(self):
        """
        Return the internal dictionary which has the combinatorial objects
        indexing the basis as keys and their corresponding coefficients as
        values.

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
            ...       print t,c
            a 1
            c 3

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1])+2*s([3,2])
            sage: d = a.monomial_coefficients()
            sage: type(d)
            <type 'dict'>
            sage: d[ Partition([2,1]) ]
            1
            sage: d[ Partition([3,2]) ]
            2
        """
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
            sage: F.print_options(generator_cmp = lambda x,y: -cmp(x,y))
            sage: f._sorted_items_for_printing()
            [('c', 2), ('b', 3), ('a', 1)]
            sage: F.print_options(generator_cmp=cmp) #reset to original state

        .. seealso:: :meth:`_repr_`, :meth:`_latex_`, :meth:`print_options`
        """
        print_options = self.parent().print_options()
        v = self._monomial_coefficients.items()
        try:
            v.sort(cmp = print_options['generator_cmp'],
                   key = lambda monomial_coeff: monomial_coeff[0])
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
            ...                               generator_cmp = lambda x,y: cmp(y,x))
            sage: e = F.basis()
            sage: e['a'] + 3*e['b'] + 2*e['c']
            2*B['c'] + 3*B['b'] + B['a']

            sage: F = CombinatorialFreeModule(QQ, ['ac', 'ba', 'cb'],
            ...                               generator_cmp = lambda x,y: cmp(x[1],y[1]))
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
        i = 0

        if scalar_mult is None:
            scalar_mult = "*"

        all_atomic = True
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
            return "0"
        elif s == empty_ascii_art:
            return AsciiArt(["1"])
        else:
            return s

    def _latex_(self):
        """
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
            ...                               generator_cmp = lambda x,y: cmp(y,x))
            sage: e = F.basis()
            sage: latex(e['a'] + 3*e['b'] + 2*e['c'])
            2B_{c} + 3B_{b} + B_{a}

            sage: F = CombinatorialFreeModule(QQ, ['ac', 'ba', 'cb'],
            ...                               generator_cmp = lambda x,y: cmp(x[1],y[1]))
            sage: e = F.basis()
            sage: latex(e['ac'] + 3*e['ba'] + 2*e['cb'])
            3B_{ba} + 2B_{cb} + B_{ac}
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult       = self.parent()._print_options['scalar_mult'],
                            latex_scalar_mult = self.parent()._print_options['latex_scalar_mult'],
                            repr_monomial = self.parent()._latex_term,
                            is_latex=True, strip_one = True)

    def __eq__(self, other):
        """
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
        if have_same_parent(self, other):
            return self._monomial_coefficients == other._monomial_coefficients
        from sage.structure.element import get_coercion_model
        import operator
        try:
            return get_coercion_model().bin_op(self, other, operator.eq)
        except TypeError:
            return False

    def __ne__(left, right):
        """
        EXAMPLES::

            sage: F1 = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: F1.an_element() != F1.an_element()
            False
            sage: F1.an_element() != F1.zero()
            True
        """
        return not left == right

    def __cmp__(left, right):
        """
        The ordering is the one on the underlying sorted list of
        (monomial,coefficients) pairs.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1])
            sage: b = s([1,1,1])
            sage: cmp(a,b) #indirect doctest
            1
        """
        if have_same_parent(left, right) and left._monomial_coefficients == right._monomial_coefficients:
            return 0
        v = sorted(mc for mc in left._monomial_coefficients.items() if mc[1] != 0)
        w = sorted(mc for mc in right._monomial_coefficients.items() if mc[1] != 0)
        return cmp(v, w)

    def _add_(self, other):
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
        assert have_same_parent(self, other)

        F = self.parent()
        return F._from_dict( dict_addition( [ self._monomial_coefficients, other._monomial_coefficients ] ), remove_zeros=False )

    def _neg_(self):
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
        return F._from_dict( dict_linear_combination( [ ( self._monomial_coefficients, -1 ) ] ), remove_zeros=False )

    def _sub_(self, other):
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
        assert have_same_parent(self, other)
        F = self.parent()
        return F._from_dict( dict_linear_combination( [ ( self._monomial_coefficients, 1 ), (other._monomial_coefficients, -1 ) ] ), remove_zeros=False )

    def _coefficient_fast(self, m, default=None):
        """
        Returns the coefficient of m in self, where m is key in
        self._monomial_coefficients.

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
            sage: a._coefficient_fast(p, 2)
            1
            sage: a._coefficient_fast(q)
            0
            sage: a._coefficient_fast(q, 2)
            2
            sage: a[p]
            1
            sage: a[q]
            0
        """
        if default is None:
            default = self.base_ring()(0)
        return self._monomial_coefficients.get(m, default)

    __getitem__ = _coefficient_fast

    def coefficient(self, m):
        """
        EXAMPLES::

            sage: s = CombinatorialFreeModule(QQ, Partitions())
            sage: z = s([4]) - 2*s([2,1]) + s([1,1,1]) + s([1])
            sage: z.coefficient([4])
            1
            sage: z.coefficient([2,1])
            -2
            sage: z.coefficient(Partition([2,1]))
            -2
            sage: z.coefficient([1,2])
            Traceback (most recent call last):
            ...
            AssertionError: [1, 2] should be an element of Partitions
            sage: z.coefficient(Composition([2,1]))
            Traceback (most recent call last):
            ...
            AssertionError: [2, 1] should be an element of Partitions

        Test that coefficient also works for those parents that do not yet have an element_class::

            sage: G = DihedralGroup(3)
            sage: F = CombinatorialFreeModule(QQ, G)
            sage: hasattr(G, "element_class")
            False
            sage: g = G.an_element()
            sage: (2*F.monomial(g)).coefficient(g)
            2
        """
        # NT: coefficient_fast should be the default, just with appropriate assertions
        # that can be turned on or off
        C = self.parent()._indices
        assert m in C, "%s should be an element of %s"%(m, C)
        if hasattr(C, "element_class") and not isinstance(m, C.element_class):
            m = C(m)
        return self._coefficient_fast(m)


    def is_zero(self):
        """
        Returns True if and only self == 0.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.is_zero()
            False
            sage: F.zero().is_zero()
            True

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: s([2,1]).is_zero()
            False
            sage: s(0).is_zero()
            True
            sage: (s([2,1]) - s([2,1])).is_zero()
            True
        """
        BR = self.parent().base_ring()
        zero = BR( 0 )
        return all( v == zero for v in self._monomial_coefficients.values() )

    def __len__(self):
        """
        Returns the number of basis elements of self with nonzero
        coefficients.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: len(f)
            2

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: len(z)
            4
        """
        return self.length()

    def length(self):
        """
        Returns the number of basis elements of self with nonzero
        coefficients.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.length()
            2

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.length()
            4
        """
        BR = self.parent().base_ring()
        zero = BR( 0 )
        return len( [ key for key, coeff in self._monomial_coefficients.iteritems() if coeff != zero ] )

    def support(self):
        """
        Returns a list of the combinatorial objects indexing the basis
        elements of self which non-zero coefficients.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.support()
            ['a', 'c']

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.support()
            [[1], [1, 1, 1], [2, 1], [4]]
        """
        BR = self.parent().base_ring()
        zero = BR( 0 )

        supp = sorted([ key for key, coeff in self._monomial_coefficients.iteritems() if coeff != zero ])

        return supp

    def monomials(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 2*B['c']
            sage: f.monomials()
            [B['a'], B['c']]

            sage: (F.zero()).monomials()
            []
        """
        P = self.parent()
        BR = P.base_ring()
        zero = BR( 0 )
        one = BR( 1 )

        supp = sorted([ key for key, coeff in self._monomial_coefficients.iteritems() if coeff != zero ])

        return [ P._from_dict( { key : one }, remove_zeros=False ) for key in supp ]

    def terms(self):
        """
        Returns a list of the terms of ``self``

        .. seealso:: :meth:`monomials`

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 2*B['c']
            sage: f.terms()
            [B['a'], 2*B['c']]
        """
        BR = self.parent().base_ring()
        zero = BR( 0 )
        v = sorted([ ( key, value ) for key, value in self._monomial_coefficients.iteritems() if value != zero ])
        from_dict = self.parent()._from_dict
        return [ from_dict( { key : value } ) for key,value in v ]

    def coefficients(self):
        """
        Returns a list of the coefficients appearing on the basis elements in
        self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f.coefficients()
            [1, -3]

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: z = s([4]) + s([2,1]) + s([1,1,1]) + s([1])
            sage: z.coefficients()
            [1, 1, 1, 1]
        """
        BR = self.parent().base_ring()
        zero = BR( 0 )
        v = sorted([ ( key, value ) for key, value in self._monomial_coefficients.iteritems() if value != zero ])
        return [ value for key,value in v ]

    def _vector_(self, new_base_ring=None):
        """
        Returns ``self`` as a dense vector

        INPUT:

        - ``new_base_ring`` -- a ring (default: None)

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
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a._vector_()
            (2, 0, 0, 0, 0, 4)
            sage: a.to_vector()
            (2, 0, 0, 0, 0, 4)
            sage: vector(a)
            (2, 0, 0, 0, 0, 4)
            sage: a == QS3.from_vector(a.to_vector())
            True

        If ''new_base_ring'' is specified, then a vector over
        ''new_base_ring'' is returned::

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
        return dense_free_module.element_class(dense_free_module,
                                               [d.get(m, 0) for m in parent.get_order()],
                                               coerce=True, copy=False)

    to_vector = _vector_

    def _acted_upon_(self, scalar, self_on_left = False):
        """
        Returns the action of a scalar on self

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

        TODO:
         - add non commutative tests
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
        if self_on_left:
            D = dict_linear_combination( [ ( D, scalar ) ], factor_on_left = False )
        else:
            D = dict_linear_combination( [ ( D, scalar ) ] )

        return F._from_dict( D, remove_zeros=False )

    # For backward compatibility
    _lmul_ = _acted_upon_
    _rmul_ = _acted_upon_

    def __div__(self, x, self_on_left=False ):
        """
        Division by coefficients

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
        if self.base_ring().is_field():
            F = self.parent()
            x = self.base_ring()( x )
            x_inv = x**-1
            D = self._monomial_coefficients
            if self_on_left:
                D = dict_linear_combination( [ ( D, x_inv ) ], factor_on_left=False )
            else:
                D = dict_linear_combination( [ ( D, x_inv ) ] )

            return F._from_dict( D, remove_zeros=False )
        else:
            return self.map_coefficients(lambda c: _divide_if_possible(c, x))

def _divide_if_possible(x, y):
    """
    EXAMPLES::

        sage: from sage.combinat.free_module import _divide_if_possible
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

class CombinatorialFreeModule(UniqueRepresentation, Module, IndexedGenerators):
    r"""
    Class for free modules with a named basis

    INPUT:

    - ``R`` - base ring

    - ``basis_keys`` - list, tuple, family, set, etc. defining the
      indexing set for the basis of this module

    - ``element_class`` - the class of which elements of this module
      should be instances (optional, default None, in which case the
      elements are instances of
      :class:`CombinatorialFreeModuleElement`)

    - ``category`` - the category in which this module lies (optional,
      default None, in which case use the "category of modules with
      basis" over the base ring ``R``); this should be a subcategory
      of :class:`ModulesWithBasis`

    For the options controlling the printing of elements, see
    :class:`~sage.structure.indexed_generators.IndexedGenerators`.

    .. NOTE::

        These print options may also be accessed and modified using the
        :meth:`print_options` method, after the module has been defined.

    EXAMPLES:

    We construct a free module whose basis is indexed by the letters a, b, c::

        sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
        sage: F
        Free module generated by {'a', 'b', 'c'} over Rational Field

    Its basis is a family, indexed by a, b, c::

        sage: e = F.basis()
        sage: e
        Finite family {'a': B['a'], 'c': B['c'], 'b': B['b']}

    ::

        sage: [x for x in e]
        [B['a'], B['b'], B['c']]
        sage: [k for k in e.keys()]
        ['a', 'b', 'c']

    Let us construct some elements, and compute with them::

        sage: e['a']
        B['a']
        sage: 2*e['a']
        2*B['a']
        sage: e['a'] + 3*e['b']
        B['a'] + 3*B['b']

    Some uses of
    :meth:`sage.categories.commutative_additive_semigroups.CommutativeAdditiveSemigroups.ParentMethods.summation`
    and :meth:`.sum`::

        sage: F = CombinatorialFreeModule(QQ, [1,2,3,4])
        sage: F.summation(F.monomial(1), F.monomial(3))
        B[1] + B[3]

        sage: F = CombinatorialFreeModule(QQ, [1,2,3,4])
        sage: F.sum(F.monomial(i) for i in [1,2,3])
        B[1] + B[2] + B[3]

    Note that free modules with a given basis and parameters are unique::

        sage: F1 = CombinatorialFreeModule(QQ, (1,2,3,4))
        sage: F1 is F
        True

    The identity of the constructed free module depends on the order of the
    basis and on the other parameters, like the prefix. Note that :class:`CombinatorialFreeModule` is
    a :class:`~sage.structure.unique_representation.UniqueRepresentation`. Hence,
    two combinatorial free modules evaluate equal if and only if they are
    identical::

        sage: F1 = CombinatorialFreeModule(QQ, (1,2,3,4))
        sage: F1 is F
        True
        sage: F1 = CombinatorialFreeModule(QQ, [4,3,2,1])
        sage: F1 == F
        False
        sage: F2 = CombinatorialFreeModule(QQ, [1,2,3,4], prefix='F')
        sage: F2 == F
        False

    Because of this, if you create a free module with certain parameters and
    then modify its prefix or other print options, this affects all modules
    which were defined using the same parameters.
    ::

        sage: F2.print_options(prefix='x')
        sage: F2.prefix()
        'x'
        sage: F3 = CombinatorialFreeModule(QQ, [1,2,3,4], prefix='F')
        sage: F3 is F2   # F3 was defined just like F2
        True
        sage: F3.prefix()
        'x'
        sage: F4 = CombinatorialFreeModule(QQ, [1,2,3,4], prefix='F', bracket=True)
        sage: F4 == F2   # F4 was NOT defined just like F2
        False
        sage: F4.prefix()
        'F'

        sage: F2.print_options(prefix='F') #reset for following doctests

    The constructed module is in the category of modules with basis
    over the base ring::

        sage: CombinatorialFreeModule(QQ, Partitions()).category()
        Category of vector spaces with basis over Rational Field

    If furthermore the index set is finite (i.e. in the category
    ``Sets().Finite()``), then the module is declared as being finite
    dimensional::

        sage: CombinatorialFreeModule(QQ, [1,2,3,4]).category()
        Category of finite dimensional vector spaces with basis over Rational Field
        sage: CombinatorialFreeModule(QQ, Partitions(3),
        ....:                         category=Algebras(QQ).WithBasis()).category()
        Category of finite dimensional algebras with basis over Rational Field

    See :mod:`sage.categories.examples.algebras_with_basis` and
    :mod:`sage.categories.examples.hopf_algebras_with_basis` for
    illustrations of the use of the ``category`` keyword, and see
    :class:`sage.combinat.root_system.weight_space.WeightSpace` for an
    example of the use of ``element_class``.

    Customizing print and LaTeX representations of elements::

        sage: F = CombinatorialFreeModule(QQ, ['a','b'], prefix='x')
        sage: original_print_options = F.print_options()
        sage: sorted(original_print_options.items())
        [('bracket', None), ('generator_cmp', <built-in function cmp>),
         ('latex_bracket', False), ('latex_prefix', None),
         ('latex_scalar_mult', None), ('prefix', 'x'),
         ('scalar_mult', '*'), ('tensor_symbol', None)]

        sage: e = F.basis()
        sage: e['a'] - 3 * e['b']
        x['a'] - 3*x['b']

        sage: F.print_options(prefix='x', scalar_mult=' ', bracket='{')
        sage: e['a'] - 3 * e['b']
        x{'a'} - 3 x{'b'}
        sage: latex(e['a'] - 3 * e['b'])
        x_{a} - 3 x_{b}

        sage: F.print_options(latex_prefix='y')
        sage: latex(e['a'] - 3 * e['b'])
        y_{a} - 3  y_{b}

        sage: F.print_options(generator_cmp = lambda x,y: -cmp(x,y))
        sage: e['a'] - 3 * e['b']
        -3 x{'b'} + x{'a'}
        sage: F.print_options(**original_print_options) # reset print options

        sage: F = CombinatorialFreeModule(QQ, [(1,2), (3,4)])
        sage: e = F.basis()
        sage: e[(1,2)] - 3 * e[(3,4)]
        B[(1, 2)] - 3*B[(3, 4)]

        sage: F.print_options(bracket=['_{', '}'])
        sage: e[(1,2)] - 3 * e[(3,4)]
        B_{(1, 2)} - 3*B_{(3, 4)}

        sage: F.print_options(prefix='', bracket=False)
        sage: e[(1,2)] - 3 * e[(3,4)]
        (1, 2) - 3*(3, 4)

    TESTS:

    Before :trac:`14054`, combinatorial free modules violated the unique
    parent condition. That caused a problem. The tensor product construction
    involves maps, but maps check that their domain and the parent of a
    to-be-mapped element are identical (not just equal). However, the tensor
    product was cached by a :class:`~sage.misc.cachefunc.cached_method`, which
    involves comparison by equality (not identity). Hence, the last line of
    the following example used to fail with an assertion error::

        sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix="F")
        sage: G = CombinatorialFreeModule(ZZ, [1,2,3,4], prefix="G")
        sage: f =   F.monomial(1) + 2 * F.monomial(2)
        sage: g = 2*G.monomial(3) +     G.monomial(4)
        sage: tensor([f, g])
        2*F[1] # G[3] + F[1] # G[4] + 4*F[2] # G[3] + 2*F[2] # G[4]
        sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix='x')
        sage: G = CombinatorialFreeModule(ZZ, [1,2,3,4], prefix='y')
        sage: f =   F.monomial(1) + 2 * F.monomial(2)
        sage: g = 2*G.monomial(3) +     G.monomial(4)
        sage: tensor([f, g])
        2*x[1] # y[3] + x[1] # y[4] + 4*x[2] # y[3] + 2*x[2] # y[4]
    """

    @staticmethod
    def __classcall_private__(cls, base_ring, basis_keys, category = None, prefix="B", **keywords):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: G = CombinatorialFreeModule(QQ, ('a','b','c'))
            sage: F is G
            True

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'], latex_bracket=['LEFT', 'RIGHT'])
            sage: F.print_options()['latex_bracket']
            ('LEFT', 'RIGHT')

            sage: F is G
            False

        We check that the category is properly straightened::

            sage: F  = CombinatorialFreeModule(QQ, ['a','b'])
            sage: F1 = CombinatorialFreeModule(QQ, ['a','b'], category = ModulesWithBasis(QQ))
            sage: F2 = CombinatorialFreeModule(QQ, ['a','b'], category = [ModulesWithBasis(QQ)])
            sage: F3 = CombinatorialFreeModule(QQ, ['a','b'], category = (ModulesWithBasis(QQ),))
            sage: F4 = CombinatorialFreeModule(QQ, ['a','b'], category = (ModulesWithBasis(QQ),CommutativeAdditiveSemigroups()))
            sage: F5 = CombinatorialFreeModule(QQ, ['a','b'], category = (ModulesWithBasis(QQ),Category.join((LeftModules(QQ), RightModules(QQ)))))
            sage: F1 is F, F2 is F, F3 is F, F4 is F, F5 is F
            (True, True, True, True, True)

            sage: G  = CombinatorialFreeModule(QQ, ['a','b'], category = AlgebrasWithBasis(QQ))
            sage: F is G
            False
        """
        if isinstance(basis_keys, (list, tuple)):
            basis_keys = FiniteEnumeratedSet(basis_keys)
        category = ModulesWithBasis(base_ring).or_subcategory(category)
        # bracket or latex_bracket might be lists, so convert
        # them to tuples so that they're hashable.
        bracket = keywords.get('bracket', None)
        if isinstance(bracket, list):
            keywords['bracket'] = tuple(bracket)
        latex_bracket = keywords.get('latex_bracket', None)
        if isinstance(latex_bracket, list):
            keywords['latex_bracket'] = tuple(latex_bracket)
        return super(CombinatorialFreeModule, cls).__classcall__(cls, base_ring, basis_keys, category = category, prefix=prefix, **keywords)

    Element = CombinatorialFreeModuleElement

    def __init__(self, R, basis_keys, element_class = None, category = None, prefix="B", **kwds):
        r"""
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])

            sage: F.category()
            Category of finite dimensional vector spaces with basis over Rational Field

        One may specify the category this module belongs to::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'], category=AlgebrasWithBasis(QQ))
            sage: F.category()
            Category of finite dimensional algebras with basis over Rational Field

            sage: F = CombinatorialFreeModule(GF(3), ['a','b','c'],
            ....:                             category=(Modules(GF(3)).WithBasis(), Semigroups()))
            sage: F.category()
            Join of Category of finite semigroups
                 and Category of finite dimensional modules with basis over Finite Field of size 3
                 and Category of vector spaces with basis over Finite Field of size 3

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'], category = FiniteDimensionalModulesWithBasis(QQ))
            sage: F.basis()
            Finite family {'a': B['a'], 'c': B['c'], 'b': B['b']}
            sage: F.category()
            Category of finite dimensional vector spaces with basis over Rational Field

            sage: TestSuite(F).run()

        TESTS:

        Regression test for :trac:10127`: ``self._indices`` needs to be
        set early enough, in case the initialization of the categories
        use ``self.basis().keys()``. This occured on several occasions
        in non trivial constructions. In the following example,
        :class:`AlgebrasWithBasis` constructs ``Homset(self,self)`` to
        extend by bilinearity method ``product_on_basis``, which in
        turn triggers ``self._repr_()``::

            sage: class MyAlgebra(CombinatorialFreeModule):
            ...       def _repr_(self):
            ...           return "MyAlgebra on %s"%(self.basis().keys())
            ...       def product_on_basis(self,i,j):
            ...           pass
            sage: MyAlgebra(ZZ, ZZ, category = AlgebrasWithBasis(QQ))
            MyAlgebra on Integer Ring

        A simpler example would be welcome!

        We check that unknown options are caught::

            sage: CombinatorialFreeModule(ZZ, [1,2,3], keyy=2)
            Traceback (most recent call last):
            ...
            ValueError: keyy is not a valid print option.
        """
        #Make sure R is a ring with unit element
        from sage.categories.all import Rings
        if R not in Rings():
            raise TypeError("Argument R must be a ring.")

        if element_class is not None:
            self.Element = element_class

        # The following is needed by e.g. root systems that don't call
        # the classcall and passes lists as basis_keys
        if isinstance(basis_keys, (list, tuple)):
            basis_keys = FiniteEnumeratedSet(basis_keys)

        # ignore the optional 'key' since it only affects CachedRepresentation
        kwds.pop('key', None)
        # This needs to be first as per #10127
        if 'monomial_cmp' in kwds:
            kwds['generator_cmp'] = kwds['monomial_cmp']
            del kwds['monomial_cmp']
        IndexedGenerators.__init__(self, basis_keys, prefix, **kwds)

        if category is None:
            category = ModulesWithBasis(R)
        elif isinstance(category, tuple):
            category = Category.join(category)
        if basis_keys in Sets().Finite():
            category = category.FiniteDimensional()

        Parent.__init__(self, base = R, category = category,
                        # Could we get rid of this?
                        element_constructor = self._element_constructor_)

        self._order = None

    # For backwards compatibility
    _repr_term = IndexedGenerators._repr_generator
    _latex_term = IndexedGenerators._latex_generator

    def _ascii_art_term(self, m):
        r"""
        Return an ascii art representation of the term indexed by ``m``.

        TESTS::

            sage: R = NonCommutativeSymmetricFunctions(QQ).R()
            sage: ascii_art(R.one())  # indirect doctest
            1
        """
        from sage.typeset.ascii_art import AsciiArt
        try:
            if m == self.one_basis():
                return AsciiArt(["1"])
        except StandardError:
            pass
        return IndexedGenerators._ascii_art_generator(self, m)

    # mostly for backward compatibility
    @lazy_attribute
    def _element_class(self):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: F._element_class
            <class 'sage.combinat.free_module.CombinatorialFreeModule_with_category.element_class'>
        """
        return self.element_class

    def _an_element_(self):
        """
        EXAMPLES::

            sage: CombinatorialFreeModule(QQ, ("a", "b", "c")).an_element()
            2*B['a'] + 2*B['b'] + 3*B['c']
            sage: CombinatorialFreeModule(QQ, ("a", "b", "c"))._an_element_()
            2*B['a'] + 2*B['b'] + 3*B['c']
            sage: CombinatorialFreeModule(QQ, ()).an_element()
            0
            sage: CombinatorialFreeModule(QQ, ZZ).an_element()
            3*B[-1] + B[0] + 3*B[1]
            sage: CombinatorialFreeModule(QQ, RR).an_element()
            B[1.00000000000000]
        """
        # Try a couple heuristics to build a not completely trivial
        # element, while handling cases where R.an_element is not
        # implemented, or R has no iterator, or R has few elements.
        x = self.zero()
        I = self.basis().keys()
        R = self.base_ring()
        try:
            x = x + self.monomial(I.an_element())
        except Exception:
            pass
        try:
            g = iter(self.basis().keys())
            for c in range(1,4):
                x = x + self.term(next(g), R(c))
        except Exception:
            pass
        return x

    # What semantic do we want for containment?
    # Accepting anything that can be coerced is not reasonnable, especially
    # if we allow coercion from the enumerated set.
    # Accepting anything that can be converted is an option, but that would
    # be expensive. So far, x in self if x.parent() == self

    def __contains__(self, x):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ,["a", "b"])
            sage: G = CombinatorialFreeModule(ZZ,["a", "b"])
            sage: F.monomial("a") in F
            True
            sage: G.monomial("a") in F
            False
            sage: "a" in F
            False
            sage: 5/3 in F
            False
        """
        return sage.structure.element.parent(x) == self # is self?

    def _element_constructor_(self, x):
        """
        Coerce x into self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ,["a", "b"])
            sage: F(F.monomial("a"))   # indirect doctest
            B['a']

        Do not rely on the following feature which may be removed in the future::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3([2,3,1])     # indirect doctest
            [2, 3, 1]

        instead, use::

            sage: P = QS3.basis().keys()
            sage: QS3.monomial(P([2,3,1]))   # indirect doctest
            [2, 3, 1]

        or:
            sage: B = QS3.basis()
            sage: B[P([2,3,1])]
            [2, 3, 1]

        TODO: The symmetric group algebra (and in general,
        combinatorial free modules on word-like object could instead
        provide an appropriate short-hand syntax QS3[2,3,1]).

        Rationale: this conversion is ambiguous in situations like::

            sage: F = CombinatorialFreeModule(QQ,[0,1])

        Is ``0`` the zero of the base ring, or the index of a basis
        element?  I.e. should the result be ``0`` or ``B[0]``?

            sage: F = CombinatorialFreeModule(QQ,[0,1])
            sage: F(0) # this feature may eventually disappear
            0
            sage: F(1)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= 1) an element of Free module generated by ... over Rational Field

        It is preferable not to rely either on the above, and instead, use::

            sage: F.zero()
            0

        Note that, on the other hand, conversions from the ground ring
        are systematically defined (and mathematically meaningful) for
        algebras.

        Conversions between distinct free modules are not allowed any
        more::

            sage: F = CombinatorialFreeModule(ZZ, ["a", "b"]);      F.rename("F")
            sage: G = CombinatorialFreeModule(QQ, ["a", "b"]);      G.rename("G")
            sage: H = CombinatorialFreeModule(ZZ, ["a", "b", "c"]); H.rename("H")
            sage: G(F.monomial("a"))
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= B['a']) an element of self (=G)
            sage: H(F.monomial("a"))
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= B['a']) an element of self (=H)

        Here is a real life example illustrating that this yielded
        mathematically wrong results::

            sage: S = SymmetricFunctions(QQ)
            sage: s = S.s(); p = S.p()
            sage: ss = tensor([s,s]); pp = tensor([p,p])
            sage: a = tensor((s[2],s[2]))

        The following originally used to yield ``p[[2]] # p[[2]]``, and if
        there was no natural coercion between ``s`` and ``p``, this would
        raise a ``NotImplementedError``. Since :trac:`15305`, this takes the
        coercion between ``s`` and ``p`` and lifts it to the tensor product. ::

            sage: pp(a)
            1/4*p[1, 1] # p[1, 1] + 1/4*p[1, 1] # p[2] + 1/4*p[2] # p[1, 1] + 1/4*p[2] # p[2]

        Extensions of the ground ring should probably be reintroduced
        at some point, but via coercions, and with stronger sanity
        checks (ensuring that the codomain is really obtained by
        extending the scalar of the domain; checking that they share
        the same class is not sufficient).

        TESTS:

        Conversion from the ground ring is implemented for algebras::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3(2)
            2*[1, 2, 3]
        """
        R = self.base_ring()
        eclass = self.element_class

        #Coerce ints to Integers
        if isinstance(x, int):
            x = Integer(x)

        if x in R:
            if x == 0:
                return self.zero()
            else:
                raise TypeError("do not know how to make x (= %s) an element of %s"%(x, self))
        #x is an element of the basis enumerated set;
        # This is a very ugly way of testing this
        elif ((hasattr(self._indices, 'element_class') and
               isinstance(self._indices.element_class, type) and
               isinstance(x, self._indices.element_class))
              or (sage.structure.element.parent(x) == self._indices)):
            return self.monomial(x)
        elif x in self._indices:
            return self.monomial(self._indices(x))
        else:
            if hasattr(self, '_coerce_end'):
                try:
                    return self._coerce_end(x)
                except TypeError:
                    pass
            raise TypeError("do not know how to make x (= %s) an element of self (=%s)"%(x,self))

    def _an_element_impl(self):
        """
        Returns an element of self, namely the zero element.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F._an_element_impl()
            0
            sage: _.parent() is F
            True
        """
        return self.element_class(self, {})

    def dimension(self):
        """
        Returns the dimension of the combinatorial algebra (which is given
        by the number of elements in the basis).

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.dimension()
            3
            sage: F.basis().cardinality()
            3
            sage: F.basis().keys().cardinality()
            3

        ::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: s.dimension()
            +Infinity
        """
        return self._indices.cardinality()

    def gens(self):
        """
        Return a tuple consisting of the basis elements of ``self``.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, ['a', 'b', 'c'])
            sage: F.gens()
            (B['a'], B['b'], B['c'])
        """
        return tuple(self.basis().values())

    def set_order(self, order):
        """
        Sets the order of the elements of the basis.

        If :meth:`set_order` has not been called, then the ordering is
        the one used in the generation of the elements of self's
        associated enumerated set.

        .. WARNING:: Many cached methods depend on this order, in
           particular for constructing subspaces and quotients.
           Changing the order after some computations have been
           cached does not invalidate the cache, and is likely to
           introduce inconsistencies.

        EXAMPLES::

            sage: QS2 = SymmetricGroupAlgebra(QQ,2)
            sage: b = list(QS2.basis().keys())
            sage: b.reverse()
            sage: QS2.set_order(b)
            sage: QS2.get_order()
            [[2, 1], [1, 2]]
        """
        self._order = order
        from sage.combinat.ranker import rank_from_list
        self._rank_basis = rank_from_list(self._order)

    @cached_method
    def get_order(self):
        """
        Returns the order of the elements in the basis.

        EXAMPLES::

            sage: QS2 = SymmetricGroupAlgebra(QQ,2)
            sage: QS2.get_order()                    # note: order changed on 2009-03-13
            [[2, 1], [1, 2]]
        """
        if self._order is None:
            self.set_order(self.basis().keys().list())
        return self._order

    def get_order_cmp(self):
        """
        Return a comparison function on the basis indices that is
        compatible with the current term order.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
            sage: Acmp = A.get_order_cmp()
            sage: sorted(A.basis().keys(), Acmp)
            ['x', 'y', 'a', 'b']
            sage: A.set_order(list(reversed(A.basis().keys())))
            sage: Acmp = A.get_order_cmp()
            sage: sorted(A.basis().keys(), Acmp)
            ['b', 'a', 'y', 'x']
        """
        self.get_order()
        return self._order_cmp

    def _order_cmp(self, x, y):
        """
        Compare `x` and `y` w.r.t. the term order.

        INPUT:

        - ``x``, ``y`` -- indices of the basis of ``self``

        OUTPUT:

        `-1`, `0`, or `1` depending on whether `x<y`, `x==y`, or `x>y`
        w.r.t. the term order.

        EXAMPLES::

            sage: A = CombinatorialFreeModule(QQ, ['x','y','a','b'])
            sage: A.set_order(['x', 'y', 'a', 'b'])
            sage: A._order_cmp('x', 'y')
            -1
            sage: A._order_cmp('y', 'y')
            0
            sage: A._order_cmp('a', 'y')
            1
        """
        return cmp( self._rank_basis(x), self._rank_basis(y) )


    @cached_method
    def _dense_free_module(self, base_ring=None):
        """
        Returns a dense free module of the same dimension

        INPUT:

        - ``base_ring`` -- a ring or ``None``

        If ``base_ring`` is None, then the base ring of ``self`` is used.

        This method is mostly used by
        :meth:`CombinatorialFreeModule.Element._vector_`

        .. SEEALSO:: :meth:`from_vector`

        EXAMPLES::

            sage: C = CombinatorialFreeModule(QQ['x'], ['a','b','c']); C
            Free module generated by {'a', 'b', 'c'} over Univariate Polynomial Ring in x over Rational Field
            sage: C._dense_free_module()
            Ambient free module of rank 3 over the principal ideal domain Univariate Polynomial Ring in x over Rational Field
            sage: C._dense_free_module(QQ['x,y'])
            Ambient free module of rank 3 over the integral domain Multivariate Polynomial Ring in x, y over Rational Field
        """
        if base_ring is None:
            base_ring = self.base_ring()
        from sage.modules.free_module import FreeModule
        return FreeModule(base_ring, self.dimension())

    def from_vector(self, vector):
        """
        Build an element of ``self`` from a (sparse) vector

        .. SEEALSO:: :meth:`get_order`, :meth:`CombinatorialFreeModule.Element._vector_`

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: b = QS3.from_vector(vector((2, 0, 0, 0, 0, 4))); b
            2*[1, 2, 3] + 4*[3, 2, 1]
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a == b
            True
        """
        cc = self.get_order()
        return self._from_dict(dict( (cc[index], coeff) for (index,coeff) in vector.iteritems()))

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: XQ = SchubertPolynomialRing(QQ)
            sage: XZ = SchubertPolynomialRing(ZZ)
            sage: XQ == XZ #indirect doctest
            False
            sage: XQ == XQ
            True
        """
        if not isinstance(other, self.__class__):
            return -1
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        return 0

    def _apply_module_morphism( self, x, on_basis, codomain=False ):
        """
        Returns the image of ``x`` under the module morphism defined by
        extending :func:`on_basis` by linearity.

        INPUT:


        - ``x`` : a element of self

        - ``on_basis`` - a function that takes in a combinatorial
          object indexing a basis element and returns an element of the
          codomain

        - ``codomain`` (optional) - the codomain of the morphism, otherwise it is computed
          using :func:`on_basis`

        If ``codomain`` is not specified, then the function tries to compute the codomain
        of the module morphism by finding the image of one of the elements in the
        support, hence :func:`on_basis` should return an element whose parent is the
        codomain.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([3]) + s([2,1]) + s([1,1,1])
            sage: b = 2*a
            sage: f = lambda part: Integer( len(part) )
            sage: s._apply_module_morphism(a, f) #1+2+3
            6
            sage: s._apply_module_morphism(b, f) #2*(1+2+3)
            12
            sage: s._apply_module_morphism(s(0), f)
            0
            sage: s._apply_module_morphism(s(1), f)
            0
            sage: s._apply_module_morphism(s(1), lambda part: len(part), ZZ)
            0
            sage: s._apply_module_morphism(s(1), lambda part: len(part))
            Traceback (most recent call last):
            ...
            ValueError: Codomain could not be determined
        """

        if x == self.zero():
            if not codomain:
                B = Family(self.basis())
                try:
                    z = B.first()
                except StopIteration:
                    raise ValueError('Codomain could not be determined')
                codomain = on_basis(z).parent()
            return codomain.zero()
        else:
            if not codomain:
                keys = x.support()
                key = keys[0]
                try:
                    codomain = on_basis(key).parent()
                except Exception:
                    raise ValueError('Codomain could not be determined')

            if hasattr( codomain, 'linear_combination' ):
                return codomain.linear_combination( ( on_basis( key ), coeff ) for key, coeff in x._monomial_coefficients.iteritems() )
            else:
                return_sum = codomain.zero()
                for key, coeff in x._monomial_coefficients.iteritems():
                    return_sum += coeff * on_basis( key )
                return return_sum

    def _apply_module_endomorphism(self, x, on_basis):
        """
        This takes in a function from the basis elements to the elements of
        self and applies it linearly to a. Note that
        _apply_module_endomorphism does not require multiplication on
        self to be defined.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).schur()
            sage: f = lambda part: 2*s(part.conjugate())
            sage: s._apply_module_endomorphism( s([2,1]) + s([1,1,1]), f)
            2*s[2, 1] + 2*s[3]
        """
        return self.linear_combination( ( on_basis( key ), coeff ) for key, coeff in x._monomial_coefficients.iteritems() )

    def sum(self, iter_of_elements):
        """
        overrides method inherited from commutative additive monoid as it is much faster on dicts directly

        INPUT:

        - ``iter_of_elements``: iterator of elements of ``self``

        Returns the sum of all elements in ``iter_of_elements``

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ,[1,2])
            sage: f = F.an_element(); f
            2*B[1] + 2*B[2]
            sage: F.sum( f for _ in range(5) )
            10*B[1] + 10*B[2]
        """

        D = dict_addition( element._monomial_coefficients for element in iter_of_elements )
        return self._from_dict( D, remove_zeros=False )

    def linear_combination( self, iter_of_elements_coeff, factor_on_left=True ):
        """
        INPUT:

        - ``iter_of_elements_coeff`` -- iterator of pairs ``(element, coeff)``
          with ``element`` in ``self`` and ``coeff`` in ``self.base_ring()``

        - ``factor_on_left`` (optional) -- if ``True``, the coefficients are
          multiplied from the left if ``False``, the coefficients are
          multiplied from the right

        Returns the linear combination `\lambda_1 v_1 + ... + \lambda_k v_k`
        (resp.  the linear combination `v_1 \lambda_1 + ... + v_k \lambda_k`)
        where ``iter_of_elements_coeff`` iterates through the sequence
        `(\lambda_1, v_1) ... (\lambda_k, v_k)`.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ,[1,2])
            sage: f = F.an_element(); f
            2*B[1] + 2*B[2]
            sage: F.linear_combination( (f,i) for i in range(5) )
            20*B[1] + 20*B[2]
        """
        return self._from_dict( dict_linear_combination( ( ( element._monomial_coefficients, coeff ) for element, coeff in iter_of_elements_coeff ), factor_on_left=factor_on_left ), remove_zeros=False )

    def term(self, index, coeff=None):
        """
        Constructs a term in ``self``

        INPUT:

        - ``index`` -- the index of a basis element
        - ``coeff`` -- an element of the coefficient ring (default: one)

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.term('a',3)
            3*B['a']
            sage: F.term('a')
            B['a']

        Design: should this do coercion on the coefficient ring?
        """
        if coeff is None:
            coeff = self.base_ring().one()
        return self._from_dict( {index: coeff} )

    def _monomial(self, index):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F._monomial('a')
            B['a']
        """
        return self._from_dict( {index: self.base_ring().one()}, remove_zeros = False )

    # This is generic, and should be lifted into modules_with_basis
    @lazy_attribute
    def monomial(self):
        """
        INPUT:

         - ''i''

        Returns the basis element indexed by i

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.monomial('a')
            B['a']

        F.monomial is in fact (almost) a map::

            sage: F.monomial
            Term map from {'a', 'b', 'c'} to Free module generated by {'a', 'b', 'c'} over Rational Field
        """
        # Should use a real Map, as soon as combinatorial_classes are enumerated sets, and therefore parents
        return PoorManMap(self._monomial, domain=self._indices, codomain=self, name="Term map")

    def _sum_of_monomials(self, indices):
        """
        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F._sum_of_monomials(['a', 'b'])
            B['a'] + B['b']
        """
        # TODO: optimize by calling directly _from_dict if we
        # know that all indices are distinct as sum_of_terms;
        # otherwise, maybe call dict_addition directly
        return self.sum(self.monomial(index) for index in indices)

    @lazy_attribute
    def sum_of_monomials(self):
        """
        INPUT:

         - ''indices'' -- an list (or iterable) of indices of basis elements

        Returns the sum of the corresponding basis elements

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.sum_of_monomials(['a', 'b'])
            B['a'] + B['b']

            sage: F.sum_of_monomials(['a', 'b', 'a'])
            2*B['a'] + B['b']

        F.sum_of_monomials is in fact (almost) a map::

            sage: F.sum_of_monomials
            A map to Free module generated by {'a', 'b', 'c'} over Rational Field
        """
        # domain = sets of self.combinatorial_class(),
        return PoorManMap(self._sum_of_monomials, codomain = self)

    def sum_of_terms(self, terms, distinct=False):
        """
        Constructs a sum of terms of ``self``

        INPUT:

        - ``terms`` -- a list (or iterable) of pairs (index, coeff)
        - ``distinct`` -- whether the indices are guaranteed to be distinct (default: ``False``)

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.sum_of_terms([('a',2), ('c',3)])
            2*B['a'] + 3*B['c']

        If ``distinct`` is True, then the construction is optimized::

            sage: F.sum_of_terms([('a',2), ('c',3)], distinct = True)
            2*B['a'] + 3*B['c']

        .. warning::

            Use ``distinct=True`` only if you are sure that the
            indices are indeed distinct::

                sage: F.sum_of_terms([('a',2), ('a',3)], distinct = True)
                3*B['a']

        Extreme case::

            sage: F.sum_of_terms([])
            0
        """
        if distinct:
            return self._from_dict(dict(terms))
        return self.sum(self.term(index, coeff) for (index, coeff) in terms)

    def monomial_or_zero_if_none(self, i):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.monomial_or_zero_if_none('a')
            B['a']
            sage: F.monomial_or_zero_if_none(None)
            0
        """
        if i is None:
            return self.zero()
        return self.monomial(i)

    @cached_method
    def zero(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.zero()
            0
        """
        return self._from_dict({}, remove_zeros=False)

    def _from_dict( self, d, coerce=False, remove_zeros=True ):
        """
        Construct an element of ``self`` from an `{index: coefficient}` dictionary

        INPUT:

        - ``d`` -- a dictionary ``{index: coeff}`` where each ``index`` is the
          index of a basis element and each ``coeff`` belongs to the
          coefficient ring ``self.base_ring()``

        - ``coerce`` -- a boolean (default: ``False``), whether to coerce the
          ``coeff``s to the coefficient ring

        - ``remove_zeros`` -- a boolean (default: ``True``), if some
          ``coeff``s may be zero and should therefore be removed

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ).elementary()
            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = e([2,1]) + e([1,1,1]); a
            e[1, 1, 1] + e[2, 1]
            sage: s._from_dict(a.monomial_coefficients())
            s[1, 1, 1] + s[2, 1]

        If the optional argument ``coerce`` is ``True``, then the
        coefficients are coerced into the base ring of ``self``::

            sage: part = Partition([2,1])
            sage: d = {part:1}
            sage: a = s._from_dict(d,coerce=True); a
            s[2, 1]
            sage: a.coefficient(part).parent()
            Rational Field

        With ``remove_zeros=True``, zero coefficients are removed::

            sage: s._from_dict({part:0})
            0

        .. warning::

            With ``remove_zeros=True``, it is assumed that no
            coefficient of the dictionary is zero. Otherwise, this may
            lead to illegal results::

                sage: list(s._from_dict({part:0}, remove_zeros=False))
                [([2, 1], 0)]
        """
        assert isinstance(d, dict)
        if coerce:
            R = self.base_ring()
            d = dict( (key, R(coeff)) for key,coeff in d.iteritems())
        if remove_zeros:
            d = dict( (key, coeff) for key, coeff in d.iteritems() if coeff)
        return self.element_class( self, d )

    @cached_method
    def tensor_unit(self, category = None, **keywords):
        if not category:
            category = self.category()
        return TensorUnit(category=category.TensorProducts(), **keywords)

class CombinatorialFreeModule_Tensor(CombinatorialFreeModule):
        """
        Tensor Product of Free Modules

        This class should not be constructed directly.
        Please use `tensor`.

        EXAMPLES:

        We construct two free modules, assign them short names, and construct their tensor product::

            sage: F = CombinatorialFreeModule(ZZ, [1,2]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [3,4]); G.__custom_name = "G"
            sage: T = tensor([F, G]); T
            F # G

            sage: T.category()
            Category of tensor products of modules with basis over Integer Ring
            sage: T.construction() # todo: not implemented
            [tensor, ]

        T is a free module, with same base ring as F and G::

            sage: T.base_ring()
            Integer Ring

        The basis of T is indexed by tuples of basis indices of F and G::

            sage: T.basis().keys()
            Image of Cartesian product of {1, 2}, {3, 4} by <type 'tuple'>
            sage: T.basis().keys().list()
            [(1, 3), (1, 4), (2, 3), (2, 4)]

        FIXME: Should elements of a CartesianProduct be tuples (making them hashable)? M.S. votes: Yes!

        Here are the basis elements themselves::

            sage: T.basis().cardinality()
            4
            sage: list(T.basis())
            [B[1] # B[3], B[1] # B[4], B[2] # B[3], B[2] # B[4]]

        The tensor product is associative and flattens sub tensor products::

            sage: H = CombinatorialFreeModule(ZZ, [5,6]); H.rename("H")
            sage: tensor([F, tensor([G, H])])
            F # G # H
            sage: tensor([tensor([F, G]), H])
            F # G # H
            sage: tensor([F, G, H])
            F # G # H

        We now compute the tensor product of elements of free modules::

            sage: f =   F.monomial(1) + 2 * F.monomial(2)
            sage: g = 2*G.monomial(3) +     G.monomial(4)
            sage: h =   H.monomial(5) +     H.monomial(6)
            sage: tensor([f, g])
            2*B[1] # B[3] + B[1] # B[4] + 4*B[2] # B[3] + 2*B[2] # B[4]

        Again, the tensor product is associative on elements::

            sage: tensor([f, tensor([g, h])]) == tensor([f, g, h])
            True
            sage: tensor([tensor([f, g]), h]) == tensor([f, g, h])
            True

        Note further that the tensor product spaces need not preexist::

            sage: t = tensor([f, g, h])
            sage: t.parent()
            F # G # H


        TESTS::

            sage: tensor([tensor([F, G]), H]) == tensor([F, G, H])
            True
            sage: tensor([F, tensor([G, H])]) == tensor([F, G, H])
            True

        Implementation notes:

        Tensor products are implemented for free modules with distinguished basis over a commutative ring `R`.
        They are implemented as modules with basis with index set given by tuples
        of the index sets of the tensor factors. They are automatically flattened so that there is no
        nesting of tensor products.

        Although the one-fold tensor product::

            sage: M = CombinatorialFreeModule(ZZ, [1,2])
            sage: T1M = tensor([M])
            sage: T1M == M
            False

        of a module is isomorphic to the module, they are different to the code: the category of the one-fold
        tensor is a subclass of ModulesWithBasis(R).TensorProducts() since it was created by the tensor construction,
        whereas the original module has category subclassing ModulesWithBasis(R).

        Mathematically there is a zero-fold tensor product which is the unit object in the tensor category.
        It is the base ring `R` viewed as a free module over itself. It is implemented as a special free
        module of rank 1::

            sage: MR = M.tensor_unit()
            sage: T1M == tensor([M, MR])
            True

        """
        @staticmethod
        def __classcall_private__(cls, modules, category, **options):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2])
                sage: G = CombinatorialFreeModule(ZZ, [3,4])
                sage: H = CombinatorialFreeModule(ZZ, [4])
                sage: tensor([tensor([F, G]), H]) == tensor([F, G, H])
                True
                sage: tensor([F, tensor([G, H])]) == tensor([F, G, H])
                True
                sage: U = F.tensor_unit()
                sage: tensor([U, G, U, F]) == tensor([G,F])
                True
                sage: U == tensor([U,U])
                True

            """
            # category should be a subcategory of a tensor category
            # this code finds the first tensor category in a join and returns its base category
            base_category = None
            if isinstance(category, Category_realization_of_parent):
                category = category.base()
            if isinstance(category, JoinCategory):
                for supercategory in category.super_categories():
                    if isinstance(supercategory, TensorProductsCategory):
                        base_category = supercategory.base_category()
                        break
            elif isinstance(category, TensorProductsCategory):
                base_category = category.base_category()
            if base_category is None:
                raise TypeError, "Category (%s) should be a subcategory of tensor products of modules"%category
            R = base_category.base_ring()
            if len(modules) > 0:
                tensor_unit = modules[0].tensor_unit(category=base_category)
                # TODO: need more accurate check for tensor units; this won't match if category disagrees
                modules = [module for module in modules if module != tensor_unit]
                if len(modules) == 0:
                    return tensor_unit
                # flatten the tensor product
                modules = sum([module._sets if isinstance(module, CombinatorialFreeModule_Tensor) else (module,) for module in modules],())
            return super(CombinatorialFreeModule_Tensor, cls).__classcall__(cls, modules, category, **options)

        def __init__(self, modules, category, **options):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2]); F
                F
                sage: tensor([F,F])
                F # F

            """
            from sage.categories.tensor import tensor
            # save the flattened tensor factors
            self._sets = modules
            base_category = None
            if isinstance(category, JoinCategory):
                for supercategory in category.super_categories():
                    if isinstance(supercategory, TensorProductsCategory):
                        if base_category is None:
                            base_category = supercategory.base_category()
                        else:
                            base_category = Category.join((base_category, supercategory.base_category()))
            elif isinstance(category, TensorProductsCategory):
                base_category = category.base_category()
            if base_category is None:
                raise TypeError, "Category (%s) should be a subcategory of tensor products of modules"%category
            R = base_category.base_ring()
            category = Category.join((category,ModulesWithBasis(R)))
            CombinatorialFreeModule.__init__(self, R, CartesianProduct(*[module.basis().keys() for module in modules]).map(tuple), category=category, **options)

            # the following is not the best option, but it's better than nothing.
            self._print_options['tensor_symbol'] = options.get('tensor_symbol', tensor.symbol)

        def _repr_(self):
            """
            This is customizable by setting
            ``self.print_options('tensor_symbol'=...)``.

            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3])
                sage: G = CombinatorialFreeModule(ZZ, [1,2,3,8])
                sage: F.rename("F")
                sage: tensor([F])
                F
                sage: G.rename("G")
                sage: T = tensor([F, G])
                sage: T # indirect doctest
                F # G
                sage: T.print_options(tensor_symbol= ' @ ')  # note the spaces
                sage: T # indirect doctest
                F @ G

            To avoid a side\--effect on another doctest, we revert the change::

                sage: T.print_options(tensor_symbol= ' # ')
            """
            from sage.categories.tensor import tensor
            if hasattr(self, "_print_options"):
                symb = self._print_options['tensor_symbol']
                if symb is None:
                    symb = tensor.symbol
            else:
                symb = tensor.symbol
            return symb.join(["%s"%module for module in self._sets])
            # TODO: make this overridable by setting _name

        def _ascii_art_(self, term):
            """
            TESTS::

                sage: R = NonCommutativeSymmetricFunctions(QQ).R()
                sage: Partitions.global_options(diagram_str="#", convention="french")
                sage: ascii_art(tensor((R[1,2], R[3,1,2])))
                R   # R
                 #     ###
                 ##      #
                         ##
            """
            from sage.categories.tensor import tensor
            if hasattr(self, "_print_options"):
                symb = self._print_options['tensor_symbol']
                if symb is None:
                    symb = tensor.symbol
            else:
                symb = tensor.symbol
            it = iter(zip(self._sets, term))
            module, t = next(it)
            rpr = module._ascii_art_term(t)
            for (module,t) in it:
                rpr += AsciiArt([symb], [len(symb)])
                rpr += module._ascii_art_term(t)
            return rpr

        _ascii_art_term = _ascii_art_

        def _latex_(self):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3])
                sage: G = CombinatorialFreeModule(ZZ, [1,2,3,8])
                sage: F.rename("F")
                sage: G.rename("G")
                sage: latex(tensor([F, F, G])) # indirect doctest
                \text{\texttt{F}} \otimes \text{\texttt{F}} \otimes \text{\texttt{G}}
                sage: F._latex_ = lambda : "F"
                sage: G._latex_ = lambda : "G"
                sage: latex(tensor([F, F, G])) # indirect doctest
                F \otimes F \otimes G
            """
            from sage.misc.latex import latex
            symb = " \\otimes "
            return symb.join(["%s"%latex(module) for module in self._sets])

        def _repr_term(self, term):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix="F")
                sage: G = CombinatorialFreeModule(ZZ, [1,2,3,4], prefix="G")
                sage: f =   F.monomial(1) + 2 * F.monomial(2)
                sage: g = 2*G.monomial(3) +     G.monomial(4)
                sage: tensor([f, g]) # indirect doctest
                2*F[1] # G[3] + F[1] # G[4] + 4*F[2] # G[3] + 2*F[2] # G[4]
            """
            from sage.categories.tensor import tensor
            symb = self._print_options['tensor_symbol']
            if symb is None:
                symb = tensor.symbol
            return symb.join(module._repr_term(t) for (module, t) in zip(self._sets, term))

        def _latex_term(self, term):
            """
            TESTS::

                sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix='x')
                sage: G = CombinatorialFreeModule(ZZ, [1,2,3,4], prefix='y')
                sage: f =   F.monomial(1) + 2 * F.monomial(2)
                sage: g = 2*G.monomial(3) +     G.monomial(4)
                sage: latex(tensor([f, g])) # indirect doctest
                2x_{1} \otimes y_{3} + x_{1} \otimes y_{4} + 4x_{2} \otimes y_{3} + 2x_{2} \otimes y_{4}
            """
            symb = " \\otimes "
            return symb.join(module._latex_term(t) for (module, t) in zip(self._sets, term))

        @cached_method
        def tensor_constructor(self, modules, **keywords):
            r"""
            INPUT:

             - ``modules`` -- a tuple `(F_1,\dots,F_n)` of
               free modules whose tensor product is self

            Returns the canonical multilinear morphism from
            `F_1 \times \dots \times F_n` to `F_1 \otimes \dots \otimes F_n`

            EXAMPLES::

                sage: F = CombinatorialFreeModule(ZZ, [1,2]); F.__custom_name = "F"
                sage: G = CombinatorialFreeModule(ZZ, [3,4]); G.__custom_name = "G"
                sage: H = CombinatorialFreeModule(ZZ, [5,6]); H.rename("H")

                sage: f =   F.monomial(1) + 2 * F.monomial(2)
                sage: g = 2*G.monomial(3) +     G.monomial(4)
                sage: h =   H.monomial(5) +     H.monomial(6)
                sage: FG  = tensor([F, G   ])
                sage: phi_fg = FG.tensor_constructor((F, G))
                sage: phi_fg(f,g)
                2*B[1] # B[3] + B[1] # B[4] + 4*B[2] # B[3] + 2*B[2] # B[4]

                sage: FGH = tensor([F, G, H])
                sage: phi_fgh = FGH.tensor_constructor((F, G, H))
                sage: phi_fgh(f, g, h)
                2*B[1] # B[3] # B[5] + 2*B[1] # B[3] # B[6] + B[1] # B[4] # B[5] + B[1] # B[4] # B[6] + 4*B[2] # B[3] # B[5] + 4*B[2] # B[3] # B[6] + 2*B[2] # B[4] # B[5] + 2*B[2] # B[4] # B[6]

                sage: phi_fg_h = FGH.tensor_constructor((FG, H))
                sage: phi_fg_h(phi_fg(f, g), h)
                2*B[1] # B[3] # B[5] + 2*B[1] # B[3] # B[6] + B[1] # B[4] # B[5] + B[1] # B[4] # B[6] + 4*B[2] # B[3] # B[5] + 4*B[2] # B[3] # B[6] + 2*B[2] # B[4] # B[5] + 2*B[2] # B[4] # B[6]
            """
            # a list l such that l[i] is True if modules[i] is readily a tensor product
            is_tensor = [isinstance(module, CombinatorialFreeModule_Tensor) for module in modules]
            # the tensor_constructor, on basis elements
            result = self.monomial * CartesianProductWithFlattening(is_tensor) #.
            # TODO: make this into an element of Hom( A x B, C ) when those will exist
            for i in range(0, len(modules)):
                result = modules[i]._module_morphism(result, position = i, codomain = self)
            return result

        def _tensor_of_elements(self, elements, **keywords):
            """
            Returns the tensor product of the specified elements.
            The result should be in self.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(ZZ, [1,2]); F.__custom_name = "F"
                sage: G = CombinatorialFreeModule(ZZ, [3,4]); G.__custom_name = "G"
                sage: H = CombinatorialFreeModule(ZZ, [5,6]); H.rename("H")

                sage: f =   F.monomial(1) + 2 * F.monomial(2)
                sage: g = 2*G.monomial(3) +     G.monomial(4)
                sage: h =   H.monomial(5) +     H.monomial(6)

                sage: GH  = tensor([G, H])
                sage: gh = GH._tensor_of_elements([g, h]); gh
                2*B[3] # B[5] + 2*B[3] # B[6] + B[4] # B[5] + B[4] # B[6]

                sage: FGH = tensor([F, G, H])
                sage: FGH._tensor_of_elements([f, g, h])
                2*B[1] # B[3] # B[5] + 2*B[1] # B[3] # B[6] + B[1] # B[4] # B[5] + B[1] # B[4] # B[6] + 4*B[2] # B[3] # B[5] + 4*B[2] # B[3] # B[6] + 2*B[2] # B[4] # B[5] + 2*B[2] # B[4] # B[6]

                sage: FGH._tensor_of_elements([f, gh])
                2*B[1] # B[3] # B[5] + 2*B[1] # B[3] # B[6] + B[1] # B[4] # B[5] + B[1] # B[4] # B[6] + 4*B[2] # B[3] # B[5] + 4*B[2] # B[3] # B[6] + 2*B[2] # B[4] # B[5] + 2*B[2] # B[4] # B[6]

                sage: tensor([F])._tensor_of_elements([f])
                B[1] + 2*B[2]

            """
            return self.tensor_constructor(tuple(element.parent() for element in elements), **keywords)(*elements)

        def _a_basis_key(self):
            r"""
            A basis key.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(ZZ,[1,2])
                sage: G = CombinatorialFreeModule(ZZ,[3,4,5])
                sage: FG = tensor([F,G])
                sage: FG._a_basis_key()
                (1, 3)
            """
            return tuple([x.an_element() for x in self.basis().keys().cc.iters])

        def a_monomial(self):
            r"""
            A monomial of ``self``.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(ZZ,[1,2])
                sage: G = CombinatorialFreeModule(ZZ,[3,4,5])
                sage: FG = tensor([F,G])
                sage: FG.a_monomial()
                B[1] # B[3]
            """
            return self.monomial(self._a_basis_key())

        def an_element(self):
            r"""
            An element of ``self``.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(ZZ,[1,2])
                sage: G = CombinatorialFreeModule(ZZ,[3])
                sage: FG = tensor([F,G])
                sage: FG.an_element()
                4*B[1] # B[3] + 4*B[2] # B[3]

            """
            return self(tensor([module.an_element() for module in self._sets]))

        def _coerce_map_from_(self, R):
            """
            Return ``True`` if there is a coercion from ``R`` into ``self`` and
            ``False`` otherwise.  The things that coerce into ``self`` are:

            - Anything with a coercion into ``self.base_ring()``.

            - A tensor algebra whose factors have a coercion into the
              corresponding factors of ``self``.

            TESTS::

                sage: C = CombinatorialFreeModule(ZZ, ZZ)
                sage: C2 = CombinatorialFreeModule(ZZ, NN)
                sage: M = C.module_morphism(lambda x: C2.monomial(abs(x)), codomain=C2)
                sage: M.register_as_coercion()
                sage: C2(C.basis()[3])
                B[3]
                sage: C2(C.basis()[3] + C.basis()[-3])
                2*B[3]
                sage: S = C.tensor(C)
                sage: S2 = C2.tensor(C2)
                sage: S2.has_coerce_map_from(S)
                True
                sage: S.has_coerce_map_from(S2)
                False
                sage: S.an_element()
                9*B[-1] # B[-1] + 3*B[-1] # B[0] + 9*B[-1] # B[1] + 3*B[0] # B[-1] + B[0] # B[0] + 3*B[0] # B[1] + 9*B[1] # B[-1] + 3*B[1] # B[0] + 9*B[1] # B[1]
                sage: S2(S.an_element())
                B[0] # B[0] + 6*B[0] # B[1] + 6*B[1] # B[0] + 36*B[1] # B[1]

            ::

                sage: C = CombinatorialFreeModule(ZZ, Set([1,2]))
                sage: D = CombinatorialFreeModule(ZZ, Set([2,4]))
                sage: f = C.module_morphism(on_basis=lambda x: D.monomial(2*x), codomain=D)
                sage: f.register_as_coercion()
                sage: T = tensor((C,C))
                sage: p = D.an_element()
                sage: T(tensor((p,p)))
                Traceback (most recent call last):
                ...
                NotImplementedError
                sage: T = tensor((D,D))
                sage: p = C.an_element()
                sage: T(tensor((p,p)))
                4*B[2] # B[2] + 4*B[2] # B[4] + 4*B[4] # B[2] + 4*B[4] # B[4]
            """
            if R in ModulesWithBasis(self.base_ring()).TensorProducts() \
                    and isinstance(R, CombinatorialFreeModule_Tensor) \
                    and len(R._sets) == len(self._sets) \
                    and all(self._sets[i].has_coerce_map_from(M)
                            for i,M in enumerate(R._sets)):
                modules = R._sets
                vector_map = [self._sets[i]._internal_coerce_map_from(M)
                              for i,M in enumerate(modules)]
                return R.module_morphism(lambda x: self._tensor_of_elements(
                        [vector_map[i](M.monomial(x[i]))
                         for i,M in enumerate(modules)]),
                                         codomain=self)

            return super(CombinatorialFreeModule_Tensor, self)._coerce_map_from_(R)

        def _tensor_of_maps(self, maps, **keywords):
            r"""
            Tensor product of maps.

            EXAMPLES::

                sage: from sage.categories.morphism import SetMorphism
                sage: Sym = SymmetricFunctions(QQ)
                sage: s = Sym.s()
                sage: id_s = s._identity_map()
                sage: antipode_map = SetMorphism(Hom(s,s),lambda x: s.antipode(x))
                sage: iS = tensor([id_s, antipode_map])
                sage: ss = tensor([s,s])
                sage: coproduct_map = s.module_morphism(on_basis=s.coproduct_on_basis, codomain=ss)
                sage: product_map = ss.module_morphism(on_basis=lambda x: s[x[0]]*s[x[1]], codomain=s)
                sage: U = s.tensor_unit()
                sage: counit_map = SetMorphism(Hom(s,U),lambda x: s.counit(x)*U.one())
                sage: unit_map = U.module_morphism(on_basis=lambda x: s.one(), codomain=s)
                sage: b = s.an_element(); b
                2*s[] + 2*s[1] + 3*s[2]
                sage: (unit_map * counit_map)(b)
                2*s[]
                sage: (product_map*iS*coproduct_map)(b)
                2*s[]

            """
            tensor_unit = self.tensor_unit()
            # delete any map from the tensor unit to itself
            maps = [map for map in maps if map.domain() != tensor_unit or map.codomain() != tensor_unit]
            if len(maps) == 0:
                return tensor_unit._identity_map()
            category = self.category()
            base_category = None
            if isinstance(category, JoinCategory):
                for supercategory in category.super_categories():
                    if isinstance(supercategory,TensorProductsCategory):
                        if base_category is None:
                            base_category = supercategory.base_category()
                        else:
                            base_category = Category.join((base_category, supercategory.base_category()))
            elif isinstance(category, TensorProductsCategory):
                base_category = category.base_category()
            if base_category is None:
                raise TypeError, "Should be a tensor category of modules"
            R = base_category.base_ring()
            module_tensor_category = ModulesWithBasis(R).TensorProducts()
            codomains = [map.codomain() for map in maps]
            # flag which codomains are tensors
            codomain_is_tensor = tuple(codomain in module_tensor_category for codomain in codomains)
            # find the sizes of the tensors. Size zero means the tensor factor consists of scalars
            codomain_n_tensor = tuple(0 if codomains[i] == tensor_unit else len(codomains[i]._sets) if codomain_is_tensor[i] else 1 for i in range(len(maps)))

            domains = [map.domain() for map in maps]
            domain_is_tensor = tuple(domain in module_tensor_category for domain in domains)
            domain_n_tensor = tuple(0 if domains[i] == tensor_unit else len(domains[i]._sets) if domain_is_tensor[i] else 1 for i in range(len(maps)))
            # the function that splits one tuple into tuples that plug into the maps
            key_splitter = CartesianProductWithUnflattening(domain_is_tensor, domain_n_tensor)
            # make a tuple whose items are (i,f) where i is the index into the tuple and f is the scalar-valued map
            scalar_part = tuple([(i,maps[i]) for i in range(len(maps)) if codomain_n_tensor[i] == 0])
            # make a tuple whose items are (i, b, f) where i is the index, b is boolean
            # where False means that f will be a constant vector and True means that f will be a vector-valued function
            vector_part = tuple((i, domain_n_tensor[i], maps[i] if domain_n_tensor[i] != 0 else maps[i](maps[i].domain().one())) \
              for i in range(len(maps)) if codomain_n_tensor[i] != 0)

            R = self.base_ring()
            def on_basis(long_key):
                keys = key_splitter(long_key)
                monomials = [domains[i].monomial(keys[i]) if domain_n_tensor[i] != 0 else () for i in range(len(maps))]
                return R.prod([f(monomials[i]) for (i,f) in scalar_part]) * tensor([f(monomials[i]) if d!=0 else f for (i,d,f) in vector_part])
            return self.module_morphism(on_basis=on_basis, codomain=tensor(codomains, category=base_category.TensorProducts()))

class CartesianProductWithFlattening(object):
    """
    A class for cartesian product constructor, with partial flattening
    """
    def __init__(self, flatten):
        """
        INPUT:

         - ``flatten`` -- a tuple of booleans

        This constructs a callable which accepts ``len(flatten)``
        arguments, and builds a tuple out them. When ``flatten[i]``,
        the i-th argument itself should be a tuple which is flattened
        in the result.

            sage: from sage.combinat.free_module import CartesianProductWithFlattening
            sage: CartesianProductWithFlattening([True, False, True, True])
            <sage.combinat.free_module.CartesianProductWithFlattening object at ...>

        """
        self._flatten = flatten

    def __call__(self, *indices):
        """
        EXAMPLES::

            sage: from sage.combinat.free_module import CartesianProductWithFlattening
            sage: cp = CartesianProductWithFlattening([True, False, True, True])
            sage: cp((1,2), (3,4), (5,6), (7,8))
            (1, 2, (3, 4), 5, 6, 7, 8)
            sage: cp((1,2,3), 4, (5,6), (7,8))
            (1, 2, 3, 4, 5, 6, 7, 8)

        """
        return sum( (i if flatten else (i,) for (i,flatten) in zip(indices, self._flatten) ), ())

class CartesianProductWithUnflattening(object):
    """
    A class for cartesian product constructor.
    """
    def __init__(self, unflatten, sizes):
        r"""
        This class computes the inverse of flattening,
        splitting a tuple of length equal to `sum(sizes)`
        into a tuple of tuples or single elements,
        such that the length of the `i`-th tuple is `sizes[i]`.

        INPUT:

         - ``unflatten`` -- a tuple of booleans
         - ``sizes`` -- a tuple of positive integers

        """
        assert len(unflatten) == len(sizes)
        # unflatten[i] must be true except when size[i] == 1.
        assert all(unflatten_this or size == 1 for (unflatten_this, size) in zip(unflatten,sizes))
        def partial_sums(sz):
            if len(sz) == 1:
                return sz
            else:
                ps = partial_sums(sz[:-1])
                return ps + tuple([ps[-1]+sz[-1]])
        if len(sizes) == 0:
            self._bounds = tuple([])
        else:
            upper_bounds = partial_sums(sizes)
            lower_bounds = tuple([0])+upper_bounds[:-1]
            self._bounds = tuple(zip(unflatten, lower_bounds, upper_bounds))

    def __call__(self, tup):
        """
        EXAMPLES::

            sage: from sage.combinat.free_module import CartesianProductWithUnflattening
            sage: un = CartesianProductWithUnflattening((True, False, True, True), (3,1,2,2))
            sage: un((1,2,3,4,5,6,7,8))
            ((1, 2, 3), 4, (5, 6), (7, 8))
            sage: un = CartesianProductWithUnflattening((True, False, True, True, True), (3,1,0,1,3))
            sage: un((1,2,3,4,5,6,7,8))
            ((1, 2, 3), 4, (), (5,), (6, 7, 8))

        """
        return tuple([tup[low:high] if unflatten else tup[low] for (unflatten,low,high) in self._bounds])

# TODO: find a way to avoid this hack to allow for cross references
CombinatorialFreeModule.Tensor = CombinatorialFreeModule_Tensor

class CombinatorialFreeModule_TensorGrouped(CombinatorialFreeModule_Tensor):
    r"""
    Concrete class for tensor products that remembers grouping.

    INPUT:

    - modules -- the modules to be tensored
    - keywords -- optional keyword parameters

    Implementation: This class remembers the tensor factors used in its creation.
    The actual object is the flattened tensor product provided by
    :class:`CombinatorialFreeModule_Tensor`.

    The methods :meth:`.indices_to_index` and :meth:`.index_to_indices` provide bijections between
    the grouped tuples that naturally index the grouped tensor product, and the flattened indices used in the
    implementation.

    """

    def __init__(self, modules, category, **keywords):
        # save the tensor factors
        self._factors = modules
        self._is_tensor = tuple(isinstance(module, CombinatorialFreeModule_Tensor) for module in modules)
        self._n_factors = tuple(len(module.factors()) if isinstance(module, CombinatorialFreeModule_TensorGrouped) else len(module._sets) if isinstance(module, CombinatorialFreeModule_Tensor) else 1 for module in modules)
        self._unflattening_function = CartesianProductWithUnflattening(self._is_tensor, self._n_factors)
        self._flattening_function = CartesianProductWithFlattening(self._is_tensor)
        flattened_modules = sum([module._sets if isinstance(module, CombinatorialFreeModule_Tensor) else (module,) for module in modules],())
        CombinatorialFreeModule_Tensor.__init__(self, flattened_modules, category, **keywords)

    def _repr_(self):
        """
        TESTS::

            sage: W=WeylGroup("A2",prefix="s")
            sage: A=W.algebra(ZZ); A.rename("A")
            sage: A2=tensor([A,A]); A2
            A # A
            sage: tensor([A2,A2])
            (A # A) # (A # A)
            sage: tensor([A, tensor([A,A])])
            A # (A # A)
            sage: tensor([A2,A2], category=ModulesWithBasis(ZZ))
            A # A # A # A

        """
        from sage.categories.tensor import tensor
        if hasattr(self, "_print_options"):
            symb = self._print_options['tensor_symbol']
            if symb is None:
                symb = tensor.symbol
        else:
            symb = tensor.symbol
        return symb.join(["(%s)"%self.factor(i) if self._n_factors[i] > 1 else "%s"%self.factor(i) for i in range(len(self.factors()))])

    def factors(self):
        r"""
        The list of tensor factors of `self`.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: A2 = tensor([A,A]) # calls this class 
            sage: A4 = tensor([A2,A2])
            sage: A4.factors()
            (A # A, A # A)

        """
        return self._factors

    def factor(self, i):
        r"""
        The `i`-th factor of ``self``.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: A2 = tensor([A,A])
            sage: A4 = tensor([A2,A2])
            sage: A4.factor(1)
            A # A
        """
        return self._factors[i]

    def indices_to_index(self):
        r"""
        The bijection from the tuple of indices of bases of `self.factors()` to the
        index set of `self.basis()`.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: A2 = tensor([A,A])
            sage: A4 = tensor([A2,A2])
            sage: a = (r([1]),r([2]))
            sage: b = (r([2,1]),r([1,2]))
            sage: A4.indices_to_index()(a,b)
            (s1, s2, s2*s1, s1*s2)

        """
        return self._flattening_function

    def index_to_indices(self):
        r"""
        The bijection from the index set of `self.basis()` to the tuples of indices
        of bases of `self.factors()`.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: A2 = tensor([A,A])
            sage: A4 = tensor([A2,A2])
            sage: A4.index_to_indices()((r([1]),r([2]),r([2,1]),r([1,2])))
            ((s1, s2), (s2*s1, s1*s2))
        """
        return self._unflattening_function

    def monomial_grouped(self, tup):
        r"""
        The monomial indexed by a grouped tuple.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: A2 = tensor([A,A])
            sage: A4 = tensor([A2,A2])
            sage: m = A4.monomial_grouped(((r([1]),r([2])),(r([1,2]),r([2,1])))); m
            B[s1] # B[s2] # B[s1*s2] # B[s2*s1]

        """
        return self.monomial(self.indices_to_index()(*tup))

    def tensor_module_morphism(self, on_basis, codomain, category=None, **keywords):
        r"""
        Returns a module morphism whose domain is a tensor product.

        This is analogous to :meth:`.module_morphism` except that only the ``on_basis`` option is
        implemented and the inputs for the ``on_basis`` function are grouped tuples instead of the
        flattened tuples secretly used in the implementation of tensor products.

        EXAMPLES::

            sage: W = WeylGroup(['A',2],prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: A2 = tensor([A,A])
            sage: A3 = tensor([A2,A])
            sage: R = lambda (a,b): (b,a)
            sage: RI = A3.tensor_module_morphism(on_basis = lambda (ab, c): A3.monomial_grouped((R(ab),c)), codomain=A3)
            sage: pair = (r([1,2]), r([2]))
            sage: RI(A3.monomial_grouped((pair,r([1]))))
            B[s2] # B[s1*s2] # B[s1]

        """
        return self.module_morphism(on_basis = lambda x: on_basis(self.index_to_indices()(x)), category=category, codomain=codomain, **keywords)

    @cached_method
    def _flat_tensor(self):
        r"""
        The flattened tensor product of ``self``.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: A2 = tensor([A,A])
            sage: A4 = tensor([A2,A2])
            sage: AAAA = A4._flat_tensor(); AAAA
            A # A # A # A
            sage: isinstance(AAAA, CombinatorialFreeModule.TensorGrouped)
            False

        """
        return tensor(self.factors(), category=ModulesWithBasis(self.base_ring()).TensorProducts())

    def from_direct_product(self, tup):
        r"""
        Given a tuple of elements of the factors of ``self``, return the corresponding element of ``self``.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: A2 = tensor([A,A])
            sage: A4 = tensor([A2,A2])
            sage: a = A2.monomial((r([1]),r([2])))
            sage: b = A2.monomial((r([2,1]),r([1,2])))
            sage: A4.from_direct_product((a,b))
            B[s1] # B[s2] # B[s2*s1] # B[s1*s2]
        """
        if len(tup) != len(self.factors()):
            raise ValueError, "Number of components of tuple is incorrect"
        try:
            return self(tensor([the_factor(element) for the_factor, element in zip(self.factors(), tup)], category=ModulesWithBasis(self.base_ring())))
        except:
            raise ValueError, "Cannot coerce tuple (%s) into tensor product %s"%(tup, self)

    def an_element(self):
        r"""
        An element of ``self``.

        EXAMPLES::

            sage: W = WeylGroup("A1",prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: A2 = tensor([A,A])
            sage: A4 = tensor([A2,A2])
            sage: A4.an_element()
            81*B[s1] # B[s1] # B[s1] # B[s1] + 27*B[s1] # B[s1] # B[s1] # B[1] + 27*B[s1] # B[s1] # B[1] # B[s1] + 9*B[s1] # B[s1] # B[1] # B[1] + 27*B[s1] # B[1] # B[s1] # B[s1] + 9*B[s1] # B[1] # B[s1] # B[1] + 9*B[s1] # B[1] # B[1] # B[s1] + 3*B[s1] # B[1] # B[1] # B[1] + 27*B[1] # B[s1] # B[s1] # B[s1] + 9*B[1] # B[s1] # B[s1] # B[1] + 9*B[1] # B[s1] # B[1] # B[s1] + 3*B[1] # B[s1] # B[1] # B[1] + 9*B[1] # B[1] # B[s1] # B[s1] + 3*B[1] # B[1] # B[s1] # B[1] + 3*B[1] # B[1] # B[1] # B[s1] + B[1] # B[1] # B[1] # B[1]

        """
        return self.from_direct_product(tuple([the_factor.an_element() for the_factor in self.factors()]))

CombinatorialFreeModule.TensorGrouped = CombinatorialFreeModule_TensorGrouped

class CombinatorialFreeModule_CartesianProduct(CombinatorialFreeModule):
    """
    An implementation of cartesian products of modules with basis

    EXAMPLES:

    We construct two free modules, assign them short names, and construct their cartesian product::

        sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
        sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
        sage: H = CombinatorialFreeModule(ZZ, [4,7]); H.__custom_name = "H"
        sage: S = cartesian_product([F, G])
        sage: S
        F (+) G
        sage: S.basis()
        Lazy family (Term map from Disjoint union of Family ({4, 5}, {4, 6}) to F (+) G(i))_{i in Disjoint union of Family ({4, 5}, {4, 6})}

    Note that the indices of the basis elements of F and G intersect non
    trivially. This is handled by forcing the union to be disjoint::

        sage: list(S.basis())
        [B[(0, 4)], B[(0, 5)], B[(1, 4)], B[(1, 6)]]

    We now compute the cartesian product of elements of free modules::

        sage: f =   F.monomial(4) + 2 * F.monomial(5)
        sage: g = 2*G.monomial(4) +     G.monomial(6)
        sage: h =   H.monomial(4) +     H.monomial(7)
        sage: cartesian_product([f,g])
        B[(0, 4)] + 2*B[(0, 5)] + 2*B[(1, 4)] + B[(1, 6)]
        sage: cartesian_product([f,g,h])
        B[(0, 4)] + 2*B[(0, 5)] + 2*B[(1, 4)] + B[(1, 6)] + B[(2, 4)] + B[(2, 7)]
        sage: cartesian_product([f,g,h]).parent()
        F (+) G (+) H

    TODO: choose an appropriate semantic for cartesian products of cartesian products (associativity?)::

        sage: S = cartesian_product([cartesian_product([F, G]), H]) # todo: not implemented
        F (+) G (+) H
    """

    def __init__(self, modules, **options):
        r"""
        TESTS::

            sage: F = CombinatorialFreeModule(ZZ, [2,4,5])
            sage: G = CombinatorialFreeModule(ZZ, [2,4,7])
            sage: C = cartesian_product([F, G]) ; C
            Free module generated by {2, 4, 5} over Integer Ring (+) Free module generated by {2, 4, 7} over Integer Ring
            sage: TestSuite(C).run()
        """
        assert(len(modules) > 0) # TODO: generalize to a family or tuple
        R = modules[0].base_ring()
        assert(all(module in ModulesWithBasis(R)) for module in modules)
        # should check the base ring
        self._sets = modules
        CombinatorialFreeModule.__init__(self, R,
            DisjointUnionEnumeratedSets(
                [module.basis().keys() for module in modules], keepkey=True),
            **options)

    def _sets_keys(self):
        """
        In waiting for self._sets.keys()

        TESTS::

            sage: F = CombinatorialFreeModule(ZZ, [2,4,5])
            sage: G = CombinatorialFreeModule(ZZ, [2,4,7])
            sage: CP = cartesian_product([F, G])
            sage: CP._sets_keys()
            [0, 1]
        """
        return range(len(self._sets))

    def _repr_(self):
        """
        TESTS::

        sage: F = CombinatorialFreeModule(ZZ, [2,4,5])
        sage: CP = cartesian_product([F, F]); CP  # indirect doctest
        Free module generated by {2, 4, 5} over Integer Ring (+) Free module generated by {2, 4, 5} over Integer Ring
        sage: F.__custom_name = "F"; CP
        F (+) F
        """
        from sage.categories.cartesian_product import cartesian_product
        return cartesian_product.symbol.join(["%s"%module for module in self._sets])
        # TODO: make this overridable by setting _name

    @cached_method
    def cartesian_embedding(self, i):
        """
        Return the natural embedding morphism of the ``i``-th
        cartesian factor (summand) of ``self`` into ``self``.

        INPUTS:

         - ``i`` -- an integer

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: phi = S.cartesian_embedding(0)
            sage: phi(F.monomial(4) + 2 * F.monomial(5))
            B[(0, 4)] + 2*B[(0, 5)]
            sage: phi(F.monomial(4) + 2 * F.monomial(6)).parent() == S
            True

        TESTS::

            sage: phi(G.monomial(4))
            Traceback (most recent call last):
            ...
            AssertionError
        """
        assert i in self._sets_keys()
        return self._sets[i]._module_morphism(lambda t: self.monomial((i, t)),
                                              codomain = self)

    summand_embedding = cartesian_embedding

    @cached_method
    def cartesian_projection(self, i):
        """
        Return the natural projection onto the `i`-th cartesian factor
        (summand) of ``self``.

        INPUTS:

         - ``i`` -- an integer

        EXAMPLE::

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: x = S.monomial((0,4)) + 2 * S.monomial((0,5)) + 3 * S.monomial((1,6))
            sage: S.cartesian_projection(0)(x)
            B[4] + 2*B[5]
            sage: S.cartesian_projection(1)(x)
            3*B[6]
            sage: S.cartesian_projection(0)(x).parent() == F
            True
            sage: S.cartesian_projection(1)(x).parent() == G
            True
        """
        assert i in self._sets_keys()
        module = self._sets[i]
        return self._module_morphism(lambda j_t: module.monomial(j_t[1]) if i == j_t[0] else module.zero(), codomain = module)

    summand_projection = cartesian_projection

    def _cartesian_product_of_elements(self, elements):
        """
        Return the cartesian product of the elements.

        INPUT:

        - ``elements`` -- a tuple with one element of each cartesian
          factor of ``self``

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: f =   F.monomial(4) + 2 * F.monomial(5)
            sage: g = 2*G.monomial(4) +     G.monomial(6)
            sage: S._cartesian_product_of_elements([f, g])
            B[(0, 4)] + 2*B[(0, 5)] + 2*B[(1, 4)] + B[(1, 6)]
            sage: S._cartesian_product_of_elements([f, g]).parent() == S
            True

        """
        return self.sum(self.summand_embedding(i)(elements[i]) for i in self._sets_keys())

    def cartesian_factors(self):
        """
        Return the factors of the cartesian product.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [4,5]); F.__custom_name = "F"
            sage: G = CombinatorialFreeModule(ZZ, [4,6]); G.__custom_name = "G"
            sage: S = cartesian_product([F, G])
            sage: S.cartesian_factors()
            (F, G)

        """
        return self._sets

    class Element(CombinatorialFreeModule.Element): # TODO: get rid of this inheritance
        pass

CombinatorialFreeModule.CartesianProduct = CombinatorialFreeModule_CartesianProduct

class TensorUnitElement(CombinatorialFreeModuleElement):
    def coefficient(self, m):
        r"""
        Return the coefficient of `m` in ``self``.

        EXAMPLES::

            sage: M=CombinatorialFreeModule(ZZ,[1,2])
            sage: U=M.tensor_unit()
            sage: u=U.an_element(); u
            2*B[()]
            sage: u.coefficient(())
            2

        """
        #assert m == tuple([])
        return self[m]

    def the_coefficient(self):
        r"""
        Return the sole coefficient of ``self``.

        Note that ``self`` is a free module with only one basis element whose index is the 0-tuple.
        """
        return self[()]

class TensorUnit(CombinatorialFreeModule.Tensor):
    r"""
    The unit in the tensor category of ModulesWithBasis over a base ring, realized by
    a zerofold `class`:CombinatorialFreeModule_Tensor`.

    This is just the base ring itself, constructed as a free module of rank 1.
    Optionally, by specifying a category, it can be viewed as an algebra, coalgebra, hopf algebra, etc.
    All the structure maps are trivial.

    This class is implemented by a direct call to CombinatorialFreeModule_Tensor.

    """

    def __init__(self, category=None, **keywords):
        """
        INPUTS:

        - category -- tensor category of result
        - keywords -- optional keyword arguments

        EXAMPLES::

            sage: M = CombinatorialFreeModule(ZZ,[1,2])
            sage: M.tensor_unit() # indirect doctest
            The unit object in Category of tensor products of modules with basis over Integer Ring

        """
        self._one_key = tuple([])
        CombinatorialFreeModule_Tensor.__init__(self, tuple([]), category, **keywords)
        assert self._one_key == self.basis().keys().an_element()

    def _repr_(self):
        """
        EXAMPLES::

            sage: M = CombinatorialFreeModule(ZZ,[1,2])
            sage: M.tensor_unit()
            The unit object in Category of tensor products of modules with basis over Integer Ring

        Should this also describe the implementation class?
        """
        return "The unit object in %s"%self.category()

    def _repr_term(self, term):
        return CombinatorialFreeModule._repr_term(self, term)

    def one_basis(self):
        """
        Returns the one of the group, which index the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(ZZ,[1])
            sage: A = M.tensor_unit(category=AlgebrasWithBasis(ZZ))
            sage: A.one_basis()
            ()
            sage: A.one()
            B[()]

        """
        return self._one_key

    def one(self):
        """
        Returns the unit element in ``self`` viewed as an algebra.

            sage: M = CombinatorialFreeModule(ZZ,[1])
            sage: A = M.tensor_unit(category=AlgebrasWithBasis(ZZ))
            sage: A.one()
            B[()]

        """
        return self.monomial(self._one_key)

    def an_element(self):
        r"""
        An element of ``self``.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(ZZ,[1])
            sage: A = M.tensor_unit()
            sage: A.an_element()
            2*B[()]

        """
        return 2 * self.one()

    def product_on_basis(self, g1, g2):
        r"""
        Product, on basis elements, as per :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis`.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(ZZ,[1])
            sage: A = M.tensor_unit(category=AlgebrasWithBasis(ZZ))
            sage: A.one_basis()
            ()
            sage: a = A.one_basis()
            sage: A.product_on_basis(a, a)
            B[()]
            sage: ma = A.monomial(a); ma
            B[()]
            sage: ma*ma
            B[()]

        """
        return self.one()

    def algebra_generators(self):
        r"""
        The generators of this algebra, as per :meth:`Algebras.ParentMethods.algebra_generators`.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(ZZ,[1])
            sage: A = M.tensor_unit(category=AlgebrasWithBasis(ZZ))
            sage: A.algebra_generators()
            Finite family {(): B[()]}

        """
        return Family(tuple([self._one_key]), self.monomial)

    def coproduct_on_basis(self, g):
        r"""
        Coproduct, on basis elements, as per :meth:`HopfAlgebrasWithBasis.ParentMethods.coproduct_on_basis`.

        Since we are taking the coproduct of the base ring and compressing tensor products, the
        codomain is automatically isomorphed to a single copy of the base ring and the
        coproduct is the identity map.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(ZZ,[1])
            sage: A = M.tensor_unit(category=HopfAlgebrasWithBasis(ZZ))
            sage: A.coproduct_on_basis(A.one_basis()) == A.one()
            True

        """
        return self.one()

    def counit_on_basis(self, g):
        r"""
        Counit, on basis elements, as per :meth:`HopfAlgebrasWithBasis.ParentMethods.counit_on_basis`.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(ZZ,[1])
            sage: A = M.tensor_unit(category=HopfAlgebrasWithBasis(ZZ))
            sage: A.counit_on_basis(A.one_basis()) == A.base_ring().one()
            True

        """
        return self.base_ring().one()

    def antipode_on_basis(self, g):
        r"""
        Antipode, on basis elements, as per :meth:`HopfAlgebrasWithBasis.ParentMethods.antipode_on_basis`.

        EXAMPLES::

            sage: M = CombinatorialFreeModule(ZZ,[1])
            sage: A = M.tensor_unit(category=HopfAlgebrasWithBasis(ZZ))
            sage: A.antipode_on_basis(A.one_basis()) == A.one()
            True

        """
        return self.one()

    Element = TensorUnitElement
