r"""
The WQSym Operad

AUTHORS:

- Florent Hivert (2011)
"""

#*****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@lri.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.element import RingElement
from sage.categories.all import SymmetricOperads, OperadsWithBasis
from sage.symbolic.ring import SR
from sage.misc.misc_c import prod
from sage.misc.cachefunc import cached_method
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words
from sage.combinat.words.shuffle_product import ShuffleProduct_overlapping



# TODO: generalize to non standard case:
# SP = DisjointUnionEnumeratedSets(Family(NN, OrderedSetPartitions))

class WQSymNewOperad(CombinatorialFreeModule):
    r"""
    The WQSym operads.

    The basis is indexed by words-partitions.
    """

    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = WQSymNewOperad(QQ); A
            The WQSym operad over Rational Field
            sage: TestSuite(A).run()

            sage: A.composition(A[{1,2}, {3}],
            ...                 A[{2}, {1}, {3}], 1)
            A[12,13,3] + A[12,3,1,3] + 2*A[12,1,3,3] + A[12,1,3]
        """
        CombinatorialFreeModule.__init__(self, R, Words(),
                                         category = OperadsWithBasis(R))

    def _repr_(self):
        """
        EXAMPLES::

            sage: WQSymNewOperad(QQ)        # indirect doctest
            The WQSym operad over Rational Field
        """
        return "The WQSym operad over %s"%(self.base_ring())

    def species(self):
        """
        The species of ordered partitions

        EXAMPLES::

            sage: f=WQSymNewOperad(QQ).species()
            sage: f.generating_series().coefficients(5)
            [0, 1, 3/2, 13/6, 25/8]
        """
        from sage.combinat.species.library import *
        E=SetSpecies().restricted(min=1)
        L=LinearOrderSpecies().restricted(min=1)
        return L(E)

    def __getitem__(self, args):
        """
        Allows for the shortcut ``A[{1,2}, {3}]``

        EXAMPLES::

            sage: A = WQSymNewOperad(QQ)
            sage: A[{1,2}, {3}]
            A[12,3]
        """
        return self._from_key(args)

    def _repr_term(self, t):
        """
        """
        res = ",".join("".join(str(x) for x in part) for part in t)
        return "A["+res+"]"

    def _fix_key(self, k):
        """
        EXAMPLES::

            sage: WQSymNewOperad(QQ)._fix_key([[1,2],[3]])
            word: frozenset([1, 2]),frozenset([3])
        """
        return self.basis().keys()(tuple(map(frozenset, k)))

    def _from_key(self, k):
        """
        EXAMPLES::

            sage: WQSymNewOperad(QQ)._from_key([[1,2],[3]])
            A[12,3]
            sage: WQSymNewOperad(QQ)._from_key("abc")
            A[a,b,c]

        """
        k = self._fix_key(k)
        return self._element_constructor(k)

    def _an_element_(self):
        """
        EXAMPLES::

            sage: A = WQSymNewOperad(QQ)
            sage: A.an_element()
            2*A[1] + A[12,3]
        """
        return ( self._from_key(({1,2},{3})) +
                 self._from_key(({1},)) + self._from_key(({1},)) )

    @cached_method
    def one_basis(self, letter=1):
        """
        Returns the word of length one, which index the one of this operad,
        as per :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = WQSymNewOperad(QQ)
            sage: A.one_basis("a")
            A[a]
        """
        return self._from_key(({letter}))


    # def operad_generators(self):
    #     """
    #     EXAMPLES::

    #     sage: WQSymNewOperad(QQ).operad_generators()
    #     Finite family {'zinbiel_product': B[word: 12]}
    #     """
    #     from sage.sets.family import Family
    #     return Family(dict([("zinbiel_product",
    #                         self._from_key([1,2]))]))

    def composition_on_basis_list(self,x,y,i):
        """
        Returns the composition of two words x o_i y as a list of
        words. i must be a label of x.

        EXAMPLES::

            sage: A = WQSymNewOperad(QQ)
            sage: A.composition_on_basis_list(
            ...       A._fix_key([[1,2],[3],[4]]),
            ...       A._fix_key([[5,6], [7]]), 3)
            [word: frozenset([1, 2]),frozenset([3, 5, 6]),frozenset([4]),frozenset([7]), word: frozenset([1, 2]),frozenset([3, 5, 6]),frozenset([7]),frozenset([4]), word: frozenset([1, 2]),frozenset([3, 5, 6]),frozenset([4, 7])]

        TESTS::

            sage: A.composition_on_basis_list(
            ...       A._fix_key([[1,2],[3],[4]]),
            ...       A._fix_key([[5,6]]), 5)
            Traceback (most recent call last):
            ...
            ValueError: The composition index is not present.
        """
        for pos, st in enumerate(x):
            if i in st: break
        else:
            raise ValueError, "The composition index is not present."
        start = x[:pos]+Words()((x[pos]|y[0] - {i},))
        end   = Words()(x[pos+1:])
        return map(lambda z: start+z,
              ShuffleProduct_overlapping(end, Words()(y[1:]),
                                       sum_func=lambda x, y: x|y,
                                       neutral = frozenset()).list())

    def composition_on_basis(self,x, y,i):
        """
        Composition of basis elements, as per :meth:`OperadsWithBasis.ParentMethods.composition_on_basis`.

        EXAMPLES::

            sage: A = WQSymNewOperad(QQ)
            sage: A.composition_on_basis(
            ...       A._fix_key([[1,2],[3],[4]]),
            ...       A._fix_key([[5,6],[7]]), 3)
            A[12,356,7,4] + A[12,356,47] + A[12,356,4,7]
            sage: A.composition(A[{1,2},{3},{4}], A[{5,6},{7}], 3)
            A[12,356,47] + A[12,356,7,4] + A[12,356,4,7]
        """
        return self.sum_of_monomials(
            t for t in self.composition_on_basis_list(x,y,i))


class WQSymFractionOperadElement(RingElement):
    """
    """

    def __init__(self, parent, frac, degree):
        """
        TESTS::

            sage: A = WQSymFractionOperad()
            sage: A.element_class(A, A.var(2), 3)
            Z2
            sage: A = WQSymFractionOperad()
            sage: el = A.from_set_partition(Set([(1,2)]))
            sage: el.orbit()
            [1/(Z1*Z2 - 1)]
        """
        self._frac = frac
        self._degree = degree
        RingElement.__init__(self, parent)

    def _repr_(self):
        """
        EXAMPLES::

            sage: A = WQSymFractionOperad()
            sage: A.from_set_partition(Set([(1,2),(3,)])).__repr__()
            '1/((Z1*Z2 - 1)*(Z1*Z2*Z3 - 1))'
        """
        return repr(self._frac)

    def _latex_(self):
        """
        """
        return self.nice_frac()._latex_()

    def degree(self):
        """

        EXAMPLES::

            sage: A = WQSymFractionOperad()
            sage: A.from_set_partition(Set([(1,2),(3,)])).degree()
            3
        """
        return self._degree

    def __cmp__(self, other):
        """
            sage: A = WQSymFractionOperad()
            sage: el1 = A.from_set_partition(Set([(1,2),(3,)]))
            sage: el2 = A.from_set_partition(Set([(1,3),(2,)]))
            sage: el1 == el1
            True
            sage: el1 == el2
            False
        """
        return (cmp(self.__class__, other.__class__) or
                cmp(self._degree, other._degree) or
                cmp(self._frac, other._frac))

    def _add_(self, other):
        """
            sage: A = WQSymFractionOperad()
            sage: el1 = A.from_set_partition(Set([(1,2),(3,)]))
            sage: el2 = A.from_set_partition(Set([(1,3),(2,)]))
            sage: el1+el2
            (Z1*Z2 + Z2 - 2)/((Z2 - 1)*(Z1*Z2 - 1)*(Z1*Z2*Z3 - 1))
        """
        parent = self.parent()
        # assert(self._degree == other._degree)
        return self.__class__(parent,
                              (self._frac + other._frac).factor(),
                              max(self._degree, other._degree))

    def compose(self, g, i):
        """
        EXAMPLES::

            sage: A = WQSymFractionOperad()
            sage: el1 = A.from_set_partition(Set([(1,2),(3,)])); el1
            1/((Z1*Z2 - 1)*(Z1*Z2*Z3 - 1))
            sage: el1.compose(A.one(), 1)
            1/((Z1*Z2 - 1)*(Z1*Z2*Z3 - 1))

        TESTS::

            sage: el1.compose(A.one(), 2) == el1
            True
            sage: el1.compose(A.one(), 3) == el1
            True
        """
        parent = self.parent()
        var = parent.var
        n = self._degree
        m = g._degree
        prd = prod(var(k) for k in range(i,i+m))
        ssf = dict()
        ssf[var(i)] = prd
        for k in range(1, n-i+1):
            ssf[var(i+k)] = var(i+k+m-1)
        ssg = dict()
        for k in range(1, m+1):
            ssg[var(k)] = var(k+i-1)
        return self.__class__(parent,
                   (g._frac.subs(ssg)*self._frac.subs(ssf)*(prd-1)).factor(),
                              m+n-1)

    def permute(self, perm):
        """
        EXAMPLES::

            sage: A = WQSymFractionOperad()
            sage: el1 = A.from_set_partition(Set([(1,2),(3,)])); el1
            1/((Z1*Z2 - 1)*(Z1*Z2*Z3 - 1))
            sage: el1.permute(Permutation([3,1,2]))
            1/((Z1*Z3 - 1)*(Z1*Z2*Z3 - 1))
        """
        parent = self.parent()
        var = parent.var
        # assert(len(perm)==self._degree)
        return self.__class__(parent,
                   self._frac.subs(dict((var(i), var(perm(i)))
                                        for i in range(1,self._degree+1))),
                              self._degree)

    def nice_frac(self):
        n = self.degree()
        return self._frac.subs(self.parent().nice_frac_subs(n))

    def num_den_subsets(self):
        num, den = self.nice_frac().numerator_denominator()
        num = num.operands()
        den = den.operands()
        dct = self.parent().num_den_subsets(self.degree())
        return tuple(dct[i] for i in num), tuple(dct[i] for i in den)


from sage.structure.unique_representation import UniqueRepresentation
class WQSymFractionOperad(UniqueRepresentation, Parent):
    def __init__(self, var_name = "Z", R=SR):
        """
        EXAMPLES::

            sage: A = WQSymFractionOperad(); A
            The WQSym operad over Symbolic Ring
            sage: TestSuite(A).run()
        """
        self._var_name = var_name
        self._base = R
        Parent.__init__(self, category = SymmetricOperads(SR))

    def _repr_(self):
        """
        EXAMPLES::

            sage: WQSymFractionOperad()
            The WQSym operad over Symbolic Ring
        """
        return "The WQSym operad over %s"%(self.base_ring())

    def _an_element_(self):
        """
        EXAMPLES::

            sage: WQSymFractionOperad().an_element()
            (Z1 - 1)/(Z1*Z2 - 1)
        """
        return self._element_constructor_(
            (self.var(1)-1)/(self.var(1)*self.var(2)-1), 2)


    def zero(self):
        """
        EXAMPLES::

            sage: WQSymFractionOperad().zero()
            0
        """
        return self._element_constructor_(0, 0)

    def composition(self, a, b, i):
        """
        """
        return a.compose(b, i)

    def var(self, i):
        r"""
        EXAMPLES::

            sage: A = WQSymFractionOperad(); A.var(2)
            Z2
        """
        from sage.calculus.var import var
        return var("%s%s"%(self._var_name, i))

    @cached_method
    def nice_frac_subs(self, n):
        from sage.combinat.subset import Subsets
        return dict((prod(self.var(i) for i in st)-1,
                     SR.symbol("["+"".join(map(str, st))+"]"))
                    for st in Subsets(n) if st)

    symbdict = dict()

    @cached_method
    def num_den_subsets(self, n):
        from sage.combinat.subset import Subsets
        for st in Subsets(n):
            self.symbdict[SR.symbol("["+"".join(map(str, st))+"]")] = st
        return self.symbdict

    @cached_method
    def one(self, *args):
        """
        EXAMPLES::

            sage: A = WQSymFractionOperad()
            sage: A.one()
            1/(Z1 - 1)
        """
        return self._element_constructor(1/(self.var(1) - 1), 1)

    def from_set_partition(self, sp):
        """
        EXAMPLES::

            sage: A = WQSymFractionOperad()
            sage: A.from_set_partition(Set([(1,2),(3,)]))
            1/((Z1*Z2 - 1)*(Z1*Z2*Z3 - 1))
            sage: [A.from_set_partition(p) for p in OrderedSetPartitions(3)]
            [1/((Z1 - 1)*(Z1*Z2 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z1 - 1)*(Z1*Z3 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z2 - 1)*(Z1*Z2 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z3 - 1)*(Z1*Z3 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z2 - 1)*(Z2*Z3 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z3 - 1)*(Z2*Z3 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z1 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z2 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z3 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z1*Z2 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z1*Z3 - 1)*(Z1*Z2*Z3 - 1)), 1/((Z2*Z3 - 1)*(Z1*Z2*Z3 - 1)), 1/(Z1*Z2*Z3 - 1)]
        """
        den = 1
        res = 1
        deg = 0
        for p in sp:
            deg += len(p)
            den = den*prod(self.var(i) for i in p)
            res = res/(den-1)
        return self._element_constructor_(res, deg)

    def _element_constructor_(self, *args, **opts):
        return self.element_class(self, *args, **opts)

    Element = WQSymFractionOperadElement

# Back compat for pickling:
# WQSymOperad = WQSymFractionOperad
# WQSymOperadElement = WQSymFractionOperad.Element
