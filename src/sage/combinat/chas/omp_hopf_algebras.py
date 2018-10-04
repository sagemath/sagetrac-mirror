"""
Hopf Algebras on Ordered Multiset Partitions into Sets

AUTHORS:

- Aaron Lauve (03-30-2018): initial implementation, following ``wqsym.py``.
"""
#*****************************************************************************
#       Copyright (C) 2018 Aaron Lauve <lauve at math.luc.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************
from six.moves import range
from functools import reduce

from sage.misc.bindable_class import BindableClass
from sage.misc.cachefunc import cached_method
#from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.misc_c import prod
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.realizations import Category_realization_of_parent
#from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.rings import Rings
from sage.categories.fields import Fields

from sage.matrix.matrix_space import MatrixSpace
from sage.sets.set import Set
from sage.rings.all import ZZ

from sage.functions.other import factorial

from sage.combinat.permutation import Permutations_mset
from sage.combinat.free_module import CombinatorialFreeModule
#from .bases import OMPBases, MultiplicativeOMPBases, OMPBasis_abstract
from sage.combinat.multiset_partition_into_sets_ordered import OrderedMultisetPartitionsIntoSets #OrderedMultisetPartitionsIntoSets, OrderedMultisetPartitionsIntoSets_n, OrderedMultisetPartitionsIntoSets_alph_d
from sage.sets.set import Set, Set_object
from sage.combinat.posets.posets import Poset
from sage.combinat.sf.sf import SymmetricFunctions

class OMPBases(Category_realization_of_parent):
    r"""
    The category of bases of `OMPSym` and `OMPQSym`.
    """
    def __init__(self, base, graded):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base`` -- an instance of `OMPSym` or `OMPQSym`
        - ``graded`` -- boolean; if the basis is graded or filtered

        TESTS::

            sage: from sage.combinat.chas.omp_hopf_algebras import OMPBases
            sage: OMPSym = OMPNonCommutativeSymmetricFunctions(ZZ)
            sage: bases = OMPBases(OMPSym, True)
            sage: OMPSym.H() in bases
            True
        """
        # TODO: do I need the next two lines?
        self._A = base._A
        self._order_grading = base._order_grading
        self._graded = graded
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.omp_hopf_algebras import OMPBases
            sage: OMPSym = OMPNonCommutativeSymmetricFunctions(ZZ)
            sage: OMPBases(OMPSym, True)
            Category of graded bases of Free Hopf Algebra on Finite Sets over Integer Ring
            sage: OMPQSym = OMPQuasiSymmetricFunctions(QQ, alphabet=[2,3,4], order_grading=False)
            sage: OMPBases(OMPQSym, False)
            Category of filtered bases of Dual of Free Hopf Algebra on Finite Sets over Rational Field
        """
        if self._graded:
            type_str = "graded"
        else:
            type_str = "filtered"
        return "Category of {} bases of {}".format(type_str, self.base())

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.omp_hopf_algebras import OMPBases
            sage: bases = OMPBases(OMPNonCommutativeSymmetricFunctions(ZZ), True)
            sage: bases.super_categories()
            [Category of realizations of Free Hopf Algebra on Finite Sets
                 over Integer Ring,
             Join of Category of realizations of hopf algebras over Integer Ring
                 and Category of graded algebras over Integer Ring,
             Category of graded connected hopf algebras with basis over Integer Ring]

            sage: bases = OMPBases(OMPQuasiSymmetricFunctions(QQ), False)
            sage: bases.super_categories()
            [Category of realizations of Dual of Free Hopf Algebra on Finite Sets
                 over Rational Field,
             Join of Category of realizations of hopf algebras over Rational Field
                 and Category of filtered algebras over Rational Field,
             Join of Category of filtered connected hopf algebras with basis over Rational Field
                 and Category of filtered algebras over Rational Field]
        """
        R = self.base().base_ring()
        cat = HopfAlgebras(R).Graded().WithBasis()
        #cat = GradedHopfAlgebras(R).WithBasis()
        if self._graded:
            cat = cat.Graded()
        else:
            cat = cat.Filtered()
        return [self.base().Realizations(),
                HopfAlgebras(R).Graded().Realizations(),
                cat.Connected()]
        #return [self.base().Realizations(),
        #        GradedHopfAlgebras(R).Realizations(),
        #        cat.Connected()]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of `OMPSym` or `OMPQSym`.

            EXAMPLES::

                sage: OMPSym = OMPNonCommutativeSymmetricFunctions(ZZ)
                sage: OMPSym.H()
                Free Hopf Algebra on Finite Sets over Integer Ring in the Homogeneous basis
                sage: OMPQSymA = OMPQuasiSymmetricFunctions(QQ, alphabet=[2,3,4])
                sage: OMPQSymA.M()
                Dual of Free Hopf Algebra on Finite Sets on alphabet {1, 4} over Rational Field with order grading in the Monomial basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def __getitem__(self, A):
            """
            Return the basis element indexed by ``A``.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            .. TODO::

                - Fix the output below. It should be ``H[{3,4,5}]``:
                    sage: H[[[3,4,5]]]
                    H[{3}, {4}, {5}]
                - Allow for input such as `H[3,4,0,5]` and `H[3,4,5]`?

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
                sage: H[[1, 3, 2]]
                H[{1,2,3}]
                sage: H[[1,3],[1]]
                H[{1,3}, {1}]
                sage: M = OMPQuasiSymmetricFunctions(QQ).M()
                sage: M[[[1, 3],[1]]]
                M[{1,3}, {1}]
            """
            try:
                return self.monomial(self._indices(A))
            except TypeError:
                return self.monomial(self._indices([A]))

        def is_field(self, proof=True):
            """
            Return whether ``self`` is a field.

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
                sage: H.is_field()
                False
            """
            return self._A == []

        def is_commutative(self):
            """
            Return whether ``self`` is commutative.

            EXAMPLES::

                sage: OMPNonCommutativeSymmetricFunctions(ZZ).H().is_commutative()
                False
                sage: OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2]).H().is_commutative()
                True
                sage: OMPQuasiSymmetricFunctions(ZZ).M().is_commutative()
                True
            """
            # TODO: refactor code.
            # if self is OMPQSym, return True

            if self.base_ring().is_zero():
                return True
            if self._A is not None and len(self._A) == 1:
                return True
            return False

        def one_basis(self):
            """
            Return the index of the unit.

            EXAMPLES::

                sage: A = OMPNonCommutativeSymmetricFunctions(QQ).H()
                sage: A.one_basis()
                []
            """
            OMP = self.basis().keys()
            return OMP([])

        def degree_on_basis(self, A):
            """
            Return the degree of an ordered multiset partition into sets
            ``A`` in Hopf Algebra of Ordered Multiset Partitions or its dual.

            If order grading is used, then the degree is ``A.order()``.
            Otherwise, the degree is ``A.size()``.

            EXAMPLES::

                sage: HA = OMPNonCommutativeSymmetricFunctions(QQ).H()
                sage: co = OrderedMultisetPartitionIntoSets([[2,1],[1,4]])
                sage: co, co.size(), co.order()
                ([{1,2}, {1,4}], 8, 4)
                sage: HA.degree_on_basis(co)
                8
                sage: HB = OMPNonCommutativeSymmetricFunctions(QQ, alphabet=[1,2,4]).H()
                sage: HB.degree_on_basis(co)
                4
                sage: HC = OMPNonCommutativeSymmetricFunctions(QQ, alphabet=[1,2,4], order_grading=False).H()
                sage: HC.degree_on_basis(co)
                8
            """
            if self._order_grading:
                return A.order()
            else:
                return A.size()

    class ElementMethods:
        def duality_pairing(self, x):
            r"""
            Compute the pairing between an element of ``self`` and an
            element of the dual.

            INPUT:

            - ``self`` -- an element of `OMPSym` or `OMPQSym`.
            - ``x`` -- an element of ``self.parent().dual_basis()``.

            OUTPUT:

            - an element of the base ring of ``self.parent()``

            EXAMPLES::

                sage: M = OMPQuasiSymmetricFunctions(QQ).M()
                sage: H = M.dual_basis()
                sage: matrix([[M(A).duality_pairing(H(B)) for A in OrderedMultisetPartitionsIntoSets(3)] for B in OrderedMultisetPartitionsIntoSets(3)])
                [1 0 0 0 0]
                [0 1 0 0 0]
                [0 0 1 0 0]
                [0 0 0 1 0]
                [0 0 0 0 1]
                sage: (M[[1,2]]*M[[3]]).duality_pairing(2*H[[1,2],[3]] - H[[1,2,3]] + 2*H[[3],[1,2]])
                3
                sage: P = OMPNonCommutativeSymmetricFunctions(QQ).P()
                sage: matrix([[M(A).duality_pairing(P(B)) for A in OrderedMultisetPartitionsIntoSets(3)] for B in OrderedMultisetPartitionsIntoSets(3)])
                [_ _ _ _ _]
                [_ _ _ _ _]
                [_ _ _ _ _]
                [_ _ _ _ _]
                [_ _ _ _ _]
            """
            x = self.parent().dual_basis()(x)
            return sum(coeff * x[I] for (I, coeff) in self)

class OMPBasis_abstract(CombinatorialFreeModule, BindableClass):
    """
    Abstract base class for bases of `OMPSym` and `OMPQSym`.

    This must define two attributes:

    - ``_prefix`` -- the basis prefix
    - ``_basis_name`` -- the name of the basis (must match one
      of the names that the basis can be constructed from OMPSym or OMPQSym)
    """
    def __init__(self, alg, graded=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: OMPNonCommutativeSymmetricFunctions(ZZ).H()
            Free Hopf Algebra on Finite Sets over Integer Ring in the Homogeneous basis

            sage: OMPQuasiSymmetricFunctions(ZZ, alphabet=[2,3,4]).M()
            Dual of Free Hopf Algebra on Finite Sets on alphabet {2, 3, 4} over Integer Ring with order grading in the Monomial basis
            sage: OMPQuasiSymmetricFunctions(ZZ, alphabet=[2,3,4], order_grading=False).M()
            Dual of Free Hopf Algebra on Finite Sets on alphabet {2, 3, 4} over Integer Ring in the Monomial basis

        TESTS::

            sage: H = OMPNonCommutativeSymmetricFunctions(ZZ).H()
            sage: TestSuite(H).run()  # long time
            sage: H = OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2]).H()
            sage: TestSuite(H).run()  # long time
            sage: H = OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2,3,4], order_grading=False).H()
            sage: TestSuite(H).run()  # long time

            sage: M = OMPQuasiSymmetricFunctions(QQ).M()
            sage: TestSuite(M).run()  # long time
            sage: M = OMPQuasiSymmetricFunctions(QQ, alphabet=[2,3,4]).M()
            sage: TestSuite(M).run()  # long time
            sage: M = OMPQuasiSymmetricFunctions(QQ, alphabet=[2,3,4], order_grading=False).M()
            sage: TestSuite(M).run()  # long time
        """
        self._A = alg._A
        self._order_grading = alg._order_grading
        if self._A:
            _OMPs = OrderedMultisetPartitionsIntoSets(alphabet=self._A)
            _OMPs._repr_ = lambda : "Ordered Multiset Partitions into Sets over alphabet %s"%self._A
        else:
            _OMPs = OrderedMultisetPartitionsIntoSets()
        CombinatorialFreeModule.__init__(self, alg.base_ring(),
                                         _OMPs,
                                         category=OMPBases(alg, graded),
                                         bracket="", prefix=self._prefix)

    def _repr_term(self, x):
        """
        Return a string representation of an element of `OMPSym` or `OMPQSym`
        in the basis ``self``.

        TESTS::

            sage: M = OMPQuasiSymmetricFunctions(QQ).M()
            sage: elt = M[[1,2]]*M[[1]]; elt
            ???
        """
        return self._prefix + repr(x).replace(", ", ",")  #.replace("}, {", "}{")

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are elements of a parent
        of similar type (`OMPSym` or `OMPQSym`) over a base with a
        coercion map into ``self.base_ring()``.

        TODO:

            - consider being more accommodating with coercions:
              + suppose x is an element of H1 and H1._A is not a subset of H2._A.
              + if each omp in support of x is over an alphabet contained in H2._A,
              + then allow H2(x) instead of throwing a TypeError.

            - add similar method in :class:`rCompBasis_OMPSym` that also absorbs NSYM.
            - add similar method in :class:`rCompBasis_OMPQSym` that also absorbs NCQSYM??

        EXAMPLES::

            sage: H = OMPNonCommutativeSymmetricFunctions(GF(7)).H(); H
            Free Hopf Algebra on Finite Sets over Finite Field of size 7 in the Homogeneous basis

        Elements of the Homogeneous basis of OMPSym canonically coerce in::

            sage: a, b = H([[1]]), H([[2,3]])
            sage: H.coerce(a+b) == a+b
            True

        Elements of the Homogeneous basis of OMPSym over `\ZZ` coerce in,
        since `\ZZ` coerces to `\GF{7}`::

            sage: HZ = OMPNonCommutativeSymmetricFunctions(ZZ).H(); HZ
            Free Hopf Algebra on Finite Sets over Integer Ring in the Homogeneous basis
            sage: aZ, bZ = HZ([[1]]), HZ([[2,3]])
            sage: c = H.coerce(aZ+bZ); c
            H[{2,3}] + H[{1}]
            sage: c.parent() is H
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so `OMPSym` over `\GF{7}`
        does not coerce to the same algebra over `\ZZ`::

            sage: HZ.coerce(b)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Free Hopf Algebra on Finite Sets
             over Finite Field of size 7 in the Homogeneous basis to
             Free Hopf Algebra on Finite Sets over Integer Ring in
             the Homogeneous basis

        Elements of the Homogeneous basis of OMPSym over the same base ring
        but with different alphabet coerce when one is contained in the other.
        This is independent of the grading chosen::

            sage: H1 = OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2,4]).H(); H1
            Free Hopf Algebra on Finite Sets on alphabet {2, 4} over the
             Integer Ring with order grading in the Homogeneous basis
            sage: H2 = OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2,3,4], order_grading=False).H()
            sage: a1, b1 = H1([[2]]), H1([[2,4]])
            sage: a2, b2 = H2([[2]]), H2([[2,3]])
            sage: c = H.coerce(a1+b1+a2+b2); c
            2*H[{2}] + H[{2,4}] + H[{2,3}]
            sage: c.parent() is H
            True
            sage: b
            sage: H1(a)  # support not built over parent's alphabet
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= H[{1}]) an element of self
             (=Free Hopf Algebra on Finite Sets on alphabet {2, 4} over the Integer Ring
             with order grading in the Homogeneous basis)
            sage: H1(a2)  # could coerce, but parent's alphabet is incompatible
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= H[{2}]) an element of self
             (=Free Hopf Algebra on Finite Sets on alphabet {2, 4} over the Integer Ring
             with order grading in the Homogeneous basis)
            sage: H2(a1+b1)  # parent's alphabet is a subset of that of H2
            H[{2}] + H[{2,4}]

        TESTS::

            sage: HZ = OMPNonCommutativeSymmetricFunctions(ZZ).H()
            sage: HQ = OMPNonCommutativeSymmetricFunctions(QQ).H()
            sage: H1 = OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2,4]).H()
            sage: H2 = OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2,3,4], order_grading=False).H()
            sage: H3 = OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2,3,4]).H()
            sage: HZ.has_coerce_map_from(HQ)
            False
            sage: HQ.has_coerce_map_from(HZ)
            True
            sage: HZ.has_coerce_map_from(QQ)
            False
            sage: HQ.has_coerce_map_from(QQ) and HQ.has_coerce_map_from(ZZ)
            True
            sage: HZ.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
            sage: all(HZ.has_coerce_map_from(HA) for HA in (H1, H2, H3))
            True
            sage: any(HA.has_coerce_map_from(HZ) for HA in (H1, H2, H3))
            False
            sage: all(HA.has_coerce_map_from(HB) for HA in (H2, H3) for HB in (H1, H2, H3))
            True
        """
        # Hopf algebras on Ordered Multiset Partitions into Sets
        # in any base should coerce in when alphabets are compatible.
        if isinstance(R, OMPBasis_abstract):
            if R.realization_of() == self.realization_of():
                return True

            elif self.base_ring().has_coerce_map_from(R.base_ring()):
                selfA, otherA = self._A, R._A
                if selfA is None:
                    cmp = True
                elif otherA is None:
                    return False
                else:
                    cmp = all(a in selfA for a in otherA)

            else:
                return False

            if cmp:
                if self._basis_name == R._basis_name: # The same basis
                    def coerce_base_ring(self, x):
                        return self._from_dict(x.monomial_coefficients())
                    return coerce_base_ring
                else: # lift that basis up, then coerce over
                    target = getattr(self.realization_of(), R._basis_name)()
                    return self._coerce_map_via([target], R)
        return super(OMPBasis_abstract, self)._coerce_map_from_(R)

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
            sage: H.an_element()
            H[{4}] + 2*H[{1}, {1,2}]
        """
        if self._A:
            A = self._A
            return self([[A[0]], A]) + 2*self([[A[0]], A[:len(A)//2], A[len(A)//2:]])
        else:
            return self([[4]]) + 2*self([[1], [1, 2]])

    def some_elements(self):
        """
        Return some elements of (dual of) Free Hopf Algebra on Finite Sets.

        EXAMPLES::

            sage: M = OMPQuasiSymmetricFunctions(QQ).M()
            sage: M.some_elements()
            [M[], M[{1}], 1/2*M[] + M[{4}] + 2*M[{1}, {1,2}],
             ???]
        """
        u = self.one()
        if self._A:
            o = self([[self._A[0]]])
        else:
            o = self([[1]])
        x = o.leading_support()
        s = self.base_ring().an_element()
        a = self.an_element()
        _a = dict(a).iteritems()
        fake_commutator = self.sum_of_terms([(x+y, c) for (y,c) in _a]) \
                    - self.sum_of_terms([(y+x, -c) for (y,c) in _a])
        return [u, o, s + a, o*o, fake_commutator]



#################
class OMPBasis_OMPSym(OMPBasis_abstract):
    """
    Add methods for `OMPSym` beyond those appearing in ``OMPBasis_abstract``
    """
    def _Coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        .. TODO::

            - describe code added here:
                if R is NSYM, bring things in, else punt to the default method.
            - add "See :meth:`OMPBasis_abstract._coerce_map_from_` for additional coercions."
        """
        # see fqsym.py: :meth:`FQSymBasis_abstract._coerce_map_from_` for hints
        from sage.combinat.ncsf_qsym.ncsf import NonCommutativeSymmetricFunctions
        S = NonCommutativeSymmetricFunctions(self.base_ring()).Complete()
        if S.has_coerce_map_from(R):
            print("do something")
            return True
        else:
            return OMPBasis_abstract._coerce_map_from_(self, R)


    class Element(OMPBasis_abstract.Element):
        def to_symmetric_function(self):
            r"""
            The projection of ``self`` to the ring of symmetric functions.

            .. TODO::

                check and add to docstring: vsp map?, alg map?, coalg map? Hopf map?
            """
            return self.to_noncommutative_symmetric_function().to_symmetric_function()

        def to_noncommutative_symmetric_function(self):
            r"""
            The projection of ``self`` to the ring `NSym` of
            noncommutative symmetric functions.

            .. TODO::

                - check and add to docstring: vsp map?, alg map?, coalg map? Hopf map?
                - Fix what follows (docstring and code).

            There is a canonical projection `\pi : OMPSym \to NSym`
            that sends ...
            This `\pi` is a Hopf algebra homomorphism.

            OUTPUT:

            - an element of the Hopf algebra of noncommutative symmetric functions
              in the Complete basis

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
                sage: (H[[1,3],[1]] - 2*H[[1],[1,3]]).to_noncommutative_symmetric_function()
                ???
                sage: X, Y = H[[1,3]], H[[1]]
                sage: X.to_noncommutative_symmetric_function() * Y.to_noncommutative_symmetric_function() == (X*Y).to_noncommutative_symmetric_function()
                True

                sage: X = H[[1,3], [1], [1,4]]; CX = X.to_noncommutative_symmetric_function()
                sage: C2 = tensor([CX.parent(), CX.parent()]); c2 = C2.zero()
                sage: XX = X.coproduct()
                sage: for ((A,B), c) in XX:
                ....:     c2 += c * H(A).to_noncommutative_symmetric_function().tensor(H(B).to_noncommutative_symmetric_function()
                sage: c2 == CX.coproduct()
                True
            """
            # TODO: fix the map, which is just a placeholder for the true definition.
            from sage.combinat.ncsf_qsym.nsym import NonCommutativeSymmetricFunctions
            S = NonCommutativeSymmetricFunctions(self.parent().base_ring()).Complete()
            H = self.parent().realization_of().H()
            return S.sum_of_terms((A.shape_from_cardinality(), coeff/prod(factorial(len(a)) for a in A))
                                  for (A, coeff) in H(self))


class OMPNonCommutativeSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    Free Hopf Algebra on Finite Sets (OMPSym).

    The Hopf algebra OMPSym is a free algebra built on finite subsets
    of positive integers. Its coproduct is defined on subsets `K` via

    .. MATH::

        \Delta(K) = \sum_{I\sqcup J = K} I \otimes J,

    then extended multiplicatively and linearly. See [LM2018]_.

    The product and coproduct are analogous to the complete/homogeneous
    basis for the Hopf algebras :class:`NonCommutativeSymmetricFunctions`
    and :class:`SymmetricFunctions`, so we use `H` as prefix for this
    natural basis. E.g., for two subsets `K` and `L`, we write

    .. MATH::

        H_{[K]} \cdot H_{[L]} = H_{[K, L]}.

    This Hopf algebra is implemented in the Homogeneous basis as a
    :class:`CombinatorialFreeModule` whose basis elements are indexed
    by *ordered multiset partitions*, which are simply lists of subsets
    of positive integers.

    There are two natural gradings on OMP: grading by the *size* or by
    the *order* of ordered multiset partitions. Only the first of these
    yields finite dimensional slices if a finite alphabet is not also
    provided. Hence, we allow an alphabet to be passed as an optional
    argument. (The user may also elect to use the size grading, even if
    a finite alphabet is passed.)

    .. SEEALSO::

        :class:`OMPQuasiSymmetricFunctions`
        :class:`OrderedMultisetPartitionsIntoSets`

    REFERENCES:

    - [LM2018]_

    EXAMPLES:

    Standard use of Free Hopf Algebra on Finite Sets:

    We begin by first creating the ring of `OMPSym` and the bases that are
    analogues of the usual non-commutative symmetric functions::

        sage: OMPSym = OMPNonCommutativeSymmetricFunctions(QQ)
        sage: H = OMPSym.H()
        sage: P = OMPSym.P() # ???
        sage: H
        Free Hopf Algebra on Finite Sets over Rational Field
         in the Homogeneous basis

    The basis is indexed by ordered multiset partitions, so we create an
    element and convert it between these bases::

        sage: elt = H(OrderedMultisetPartitionIntoSets([[2], [1,3]])) - 2*H(OrderedMultisetPartitionIntoSets([[1,2,3]])); elt
        -2*H[{1,2,3}] + H[{2}, {1,3}]
        sage: P(elt)
        ???

    There is also a shorthand for creating elements. We note that we must
    use ``P[[]]`` to create the empty ordered multiset partition due to
    python's syntax. ::

        sage: elth = H[[2], [1,3]] - 2*H[[1,2,3]]; elth
        -2*H[{1,2,3}] + H[{2}, {1,3}]
        sage: eltp = P[[]] + elth; eltp
        ???

    Restricting the alphabet and alternate grading:

    Instead of working over ordered multiset partitions of positive integers,
    we may restrict each block's elements to come from a finite alphabet `A`::

        sage: OMPSymA = OMPNonCommutativeSymmetricFunctions(QQ, alphabet=[3,4,5])
        sage: HA = OMPSymA.HA()
        sage: HA.an_element()
        ???
        sage: H.an_element()
        ???
        sage: H.an_element() in OMPSymA
        False???
        sage: HA.an_element() in OMPSym
        True??? # I'm happy with either answer, I suppose

    The homogeneous component of degree `d` of the basis is indexed by
    ordered multiset partitions over `A` of order `d`::

        sage: HA.basis(3).keys().list()
        ???

    If the keyword argument ``order_grading`` is set to ``False`` (or, if the
    ``alphabet`` keyword is not used), then homogeneous degree is determined
    by size of ordered multiset partitions::

        sage: OMPSymB = OMPNonCommutativeSymmetricFunctions(QQ, alphabet=[3,4,5], order_grading=False)
        sage: HB = OMPSymB.HB()
        sage: HB.basis(12).keys().list()
        ???

    TESTS::

        sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
        sage: a = H[OrderedMultisetPartitionIntoSets([[1]])]
        sage: b = H[OrderedMultisetPartitionsIntoSets(1)([[1]])]
        sage: c = H[[1]]
        sage: a == b == c
        True
        sage: H = OMPNonCommutativeSymmetricFunctions(QQ, alphabet=[2]).H()
        sage: a = H[OrderedMultisetPartitionIntoSets([[2]])]
        sage: b = H[OrderedMultisetPartitionsIntoSets(2)([[2]])]
        sage: c = H[[2]]
        sage: d = H[OrderedMultisetPartitionsIntoSets([1,2], 1)([[2]])]
        sage: a == b == c == d
        True

    .. TODO::

        - Add/correct the powersum basis.
        - Add/correct maps to NCSym, NSym, Sym.
        - Implement dual basis of `OMPSym` (which I call `OMPQSym` below).
    """
    @staticmethod
    def __classcall_private__(cls, R, alphabet=None, order_grading=None):
        """
        Normalize passed arguments.
        """
        # ensure R is a commutative ring
        if R not in CommutativeRings():
            raise ValueError("argument R (=%s) must be a commutative ring" % repr(R))

        # make alphabet hashable
        if alphabet is not None:
            if alphabet in ZZ:
                _A = Set(range(1,ZZ(alphabet)+1))
            else:
                _A = Set(ZZ(a) for a in alphabet)
            if _A == Set() or any(a not in ZZ for a in _A):
                raise ValueError("keyword alphabet was converted to %s, which must be a nonempty set of positive integers"%(self._A))
        else:
            # treat as ``_A`` is PositiveIntegers
            _A = None

        # pick an appropriate value for grading
        if order_grading is None:
            if _A:
                _order_grading = True
            else:
                _order_grading = False
        else:
            _order_grading = order_grading

        return super(OMPNonCommutativeSymmetricFunctions, cls).__classcall__(cls, R, _A, _order_grading)

    def __init__(self, R, alphabet, order_grading):
        """
        Initialize ``self``.

        INPUT:

        - ``R`` -- a commutative ring
        - ``alphabet`` -- a finite set of positive integers (value ``None``
          is used to represent all positive integers)
        - ``order_grading`` -- boolean: ``False`` indicates grading by size,
          while ``True`` indicates grading by order

        TESTS::

            sage: A = OMPNonCommutativeSymmetricFunctions(QQ)
            sage: TestSuite(A).run()  # long time
        """
        category = GradedHopfAlgebras(R).Connected()  # TODO: add Cocommutative
        Parent.__init__(self, base=R, category=category.WithRealizations())
        self._A = alphabet
        self._order_grading = order_grading

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: OMPNonCommutativeSymmetricFunctions(ZZ)
            Free Hopf Algebra on Finite Sets over Integer Ring
            sage: OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2,3])
            Free Hopf Algebra on Finite Sets on alphabet {2, 3} over Integer Ring with order grading
            sage: OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2,3], order_grading=False)
            Free Hopf Algebra on Finite Sets on alphabet {2, 3} over Integer Ring
        """
        OMP = "Free Hopf Algebra on Finite Sets"
        if self._A:
            OMP += " on alphabet %s"%self._A
        OMP += " over the %s"%self.base_ring()
        if self._order_grading:
            OMP += " with order grading"
        return OMP

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the `H` basis).

        EXAMPLES::

            sage: OMPNonCommutativeSymmetricFunctions(QQ).a_realization()
            Free Hopf Algebra on Finite Sets over Rational Field in the Homogeneous basis
        """
        return self.H()

    _shorthands = tuple(['H', 'P'])

    def dual(self):
        r"""
        Return the Hopf algebra that is the graded dual of
        Free Hopf Algebra on Finite Sets.

        EXAMPLES::

            sage: OMPNonCommutativeSymmetricFunctions(QQ).dual()
            Dual of Free Hopf Algebra on Finite Sets over the Rational Field
        """
        return OMPQuasiSymmetricFunctions(self.base_ring(), self._A, self._order_grading)

    class H(OMPBasis_OMPSym):
        r"""
        The Homogeneous basis of Free Hopf Algebra on Finite Sets.

        The family `(H_\mu)`, as `\mu` ranges over all ordered multiset
        partitions into sets, is identified with the natural "free algebra
        on finite sets" basis of Free Hopf Algebra on Finite Sets. That is,
        product is given by concatenation of ordered multiset partitions.
        E.g., for two subsets `K` and `L`, we write

        .. MATH::

            H_{[K]} \cdot H_{[L]} = H_{[K, L]}.

        The coproduct formula for ordered multiset partitions with one block `K`
        is given by

        .. MATH::

            \Delta(H_{[K]}) = \sum_{I\sqcup J = K} H_{[I]} \otimes H_{[J]},

        which is then extended multiplicatively to all of `(H_\mu)`.

        EXAMPLES::

            sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
            sage: H[[2], [1,3]] - 2*H[[1,2,4]]
            H[{2},{1,3}] - 2*H[{1,2,4}]

            sage: H[[2], [1,3]] * H[[1,2],[4]]
            H[{2},{1,3},{1,2},{4}]

            sage: H[[4]].coproduct()
            H[] # H[{4}] + H[{4}] # H[]
            sage: H[[1,2]].coproduct()
            H[] # H[{1,2}] + H[{1}] # H[{2}] + H[{2}] # H[{1}] + H[{1,2}] # H[]
            sage: H[[4],[1,2]].coproduct()
            H[] # H[{4},{1,2}] + H[{1}] # H[{4},{2}] + H[{2}] # H[{4},{1}]
             + H[{1,2}] # H[{4}] + H[{4}] # H[{1,2}] + H[{4},{1}] # H[{2}]
             + H[{4},{2}] # H[{1}] + H[{4},{1,2}] # H[]

        TESTS::

            sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
            sage: TestSuite(H).run()  # long time

        .. TODO: Do we need these tests, or are they checked within TestSuite?::

            sage: tests = []
            sage: for ha in H.basis(3):
            ....:     tests.append(H.zero() == sum(H[a1].antipode() * H[a2] * coeff \
            ....:           for ((a1,a2),coeff) in (ha.coproduct())))
            ....:     tests.append(H.zero() == sum(H(a1) * H(a2).antipode() * coeff \
            ....:           for ((a1,a2),coeff) in (ha.coproduct())))
            ....:     for hb in [H[[]]] + list(H.basis(1)) + list(H.basis(2)):
            ....:         tests.append((ha * hb).coproduct() == ha.coproduct() * hb.coproduct())
            sage: all(tests)
        """
        _prefix = "H"
        _basis_name = "Homogeneous"

        def product_on_basis(self, A, B):
            r"""
            Return the (associative) `*` product of the basis elements
            of ``self`` indexed by the ordered multiset partitions
            `A` and `B`.

            INPUT:

            - ``A, B`` -- two ordered multiset partitions

            OUTPUT:

            - The basis element of ``self`` indexed by the concatenation
              `A + B`. (This product is non-commutative.)

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
                sage: H[{2,3}] * H[{2,3}, {2}]
                H[{2,3}, {2,3}, {2}]

            TESTS::

                sage: one = OrderedMultisetPartitionIntoSets([])
                sage: all(H.product_on_basis(one, z) == H(z) == H.basis()[z] for z in OrderedMultisetPartitionsIntoSets(3))
                True
                sage: all(H.product_on_basis(z, one) == H(z) == H.basis()[z] for z in OrderedMultisetPartitionsIntoSets(3))
                True
            """
            return self.monomial(A + B)

        def coproduct_on_basis(self, A):
            r"""
            Return the coproduct, defined on subsets `K` via
            .. MATH::

                \sum_{I \sqcup J} = K I \otimes J,

            then extended multiplicatively and linearly. This coproduct is cocommutative.
            It is analogous to the Complete basis ring in the Hopf algebra of
            :class:`non-commutative symmetric functions<NonCommutativeSymmetricFunctions>`.

            INPUT:

            - ``A`` -- a set partition

            OUTPUT:

            - The coproduct applied to the element of OMP indexed by ``A`` expressed in the
              Homogeneous basis.

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
                sage: H[{2,3}].coproduct()
                H[] # H[{2,3}] + H[{2}] # H[{3}] + H[{3}] # H[{2}] + H[{2,3}] # H[]
                sage: H[{2,3}, {2}].coproduct()
                H[] # H[{2,3}, {2}] + H[{2}] # H[{3}, {2}] + H[{3}] # H[{2}, {2}] + H[{2,3}] # H[{2}] + H[{2}] # H[{2,3}] + H[{2}, {2}] # H[{3}] + H[{3}, {2}] # H[{2}] + H[{2,3}, {2}] # H[]
            """
            return self.tensor_square()._from_dict( A.split_blocks(2) )

        def _Antipode_on_basis(self, A):
            r"""
            Return the antipode applied to a Homogeneous basis element.

            If `A` is an ordered multiset partition with a single block then
            the antipode is given by

            .. MATH::

                S(A) = \sum_{B \prec A} (-1)^{\ell(B)} B,

            where we sum over all ordered multiset partitions that refine `A`.

            If `A` has more than one block, then we extend using the fact that
            the antipode is an algebra antimorphism.

            INPUT:

            - ``A`` -- an ordered multiset partition

            .. TODO::

                - Check that this method is faster than the default
                  implementation coming from the fact that ``self`` is
                  an instance of a connected Hopf algebra.
                - Check that it is correct!

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
                sage: H[{2,3}].antipode()
                -H[{2,3}] + H[{2}, {3}] + H[{3}, {2}]
                sage: H[{2,3}, {2}].antipode()
                H[{2}, {2,3}] - H[{2}, {2}, {3}] - H[{2}, {3}, {2}]
                sage: HH = H[[1,3],[5],[1,4]].coproduct()
                sage: HH.apply_multilinear_morphism(lambda x,y: x.antipode()*y)
                0
            """
            out = []
            P = self.basis().keys()
            AA = [P([a]).finer() for a in A.reversal()]
            for refinement in cartesian_product(AA):
                C = reduce(lambda a,b: a + b, refinement, P([]))
                ell = len(C)
                out.append((C, (-1)**ell))
            return self.sum_of_terms(out)

        def dual_basis(self):
            r"""
            Return the dual basis to the `\mathbf{H}` basis.

            The dual basis to the `\mathbf{H}` basis is the Monomial basis
            of Dual of Free Hopf Algebra on Finite Sets.

            OUTPUT:

            - the Monomial basis of Dual of Free Hopf Algebra on Finite Sets

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
                sage: H.dual_basis()
                Dual of Free Hopf Algebra on Finite Sets over the Rational Field in the Monomial basis
            """
            return self.realization_of().dual().M()

    Homogeneous = H

    class P(OMPBasis_OMPSym):
        r"""
        The Powersum basis of Free Hopf Algebra on Finite Sets.

        The family `(P_\mu)`, as `\mu` ranges over all ordered multiset
        partitions of nonnegative integers, is called the *Powersum basis*
        here. It is defined via a unitriangular change of basis from the
        Homogeneous basis of `OMPSym` via *some partial order* on
        ordered multiset partitions.

        The principle feature of the Powersum basis is that for subsets `K`,
        `P_{[K]}` is primitive.

        EXAMPLES::

            sage: P = OMPNonCommutativeSymmetricFunctions(QQ).Powersum()
            sage: P[[2], [1,3]] - 2*P[[1,2,3]]
            P[{2}, {1,3}] - 2*P[{1,2,3}]
            sage: P[[2,3]].coproduct()
            ???
            sage: P[[2,3],[1,2]].coproduct()
            ???
            sage: P[[2,3]] * P[[1,2]]
            ???
            sage: p = P.an_element(); p
            ???
            sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
            sage: H(p)
            ???

        TESTS::

            sage: TestSuite(P).run()  # long time
            sage: all(P(H(p)) == p for p in P.some_elements())
            True
        """
        _prefix = "P"
        _basis_name = "Powersum"

        def __init__(self, alg):
            """
            Initialize ``self``.
            """
            OMPBasis_OMPSym.__init__(self, alg)

            # Register coercions
            # TODO: check (uni)triangular claim(s)
            H = self.realization_of().H()
            phi = self.module_morphism(self._P_to_H_on_basis, codomain=H, unitriangular="upper")
            phi.register_as_coercion()
            (~phi).register_as_coercion()

        def _P_to_H_on_basis(self, A):
            """
            Return `P_A` in terms of the Homogeneous basis.

            INPUT:

            - ``A`` -- an ordered multiset partition

            OUTPUT:

            - An element of the Homogeneous basis

            .. TODO::

                - Fix this method! As a placeholder, I am
                  currently just taking a definition from `NCSym`.
                - What is "remove_zeros" meant to do? (should I use it?)

            EXAMPLES::

                sage: P = OMPNonCommutativeSymmetricFunctions(QQ).P()
                sage: A = OrderedMultisetPartitionIntoSets([[1,2,3]])
                sage: P._P_to_H_on_basis(A)
                ???
                sage: B = OrderedMultisetPartitionIntoSets([[1,2],[3]])
                sage: P._P_to_H_on_basis(B)
                ???
            """
            H = self.realization_of().H()
            P_refine = Poset((A.finer(), lambda a,b: a.is_finer(b)))  #A.parent().is_finer
            #print A, A.finer(), P_refine
            c = abs(prod((-1)**(i-1) * factorial(i-1) for i in A.shape_from_size()))
            R = self.base_ring()
            return H._from_dict({B: R(P_refine.moebius_function(B, A) / R(c))
                                 for B in list(P_refine)}) #, remove_zeros=False)

        def primitive(self, A, i=None):
            r"""
            Return the primitive associated to ``A`` in ``self``.

            Suppose `A` is an ordered multiset partition with a single block.
            Return the primitive

            .. MATH::

                p(A) = \sum_{B} (-1)^{\ell(B)-1} B,

            where the sum is over all ordered set partitions `B` that are finer
            than `A` and that have `i \in A` in their first block.

            .. TODO::

                If `A` has more than one block, then create some sort of
                Lie polynomial from its constituent blocks.

            INPUT:

            - ``A`` -- an ordered multiset partition
            - ``i`` -- (default: A[0][0]) index in the base multiset
                       for ``A`` specifying which family of primitives output
                       should belong to

            OUTPUT:

            - an element in the basis ``self``

            EXAMPLES::

                sage: P = OMPNonCommutativeSymmetricFunctions(QQ).Powersum()
                sage: elt = P.primitive(OrderedMultisetPartitionIntoSets([[1,3], [1,2]])); elt
                -p{{1, 2}, {3}} + p{{1, 3}, {2}}
                sage: elt.coproduct()
                -p{} # p{{1, 2}, {3}} + p{} # p{{1, 3}, {2}} - p{{1, 2}, {3}} # p{} + p{{1, 3}, {2}} # p{}
                sage: p.primitive(SetPartition([[1], [2,3]]))
                0
                sage: p.primitive(SetPartition([]))
                p{}
            """
            if len(A) == 0:
                return self.zero()
            elif len(A) > 1:
                # TODO: finish the job!
                return None
            else:
                # TODO:
                # Is ``.finer()`` the best/only set to sum over?
                # What about setting strong=True?
                P = OrderedMultisetPartitionsIntoSets()
                if i not in A[0]:
                    i = A[0][0]
                return self.sum_of_terms((finer, (-1)**len(finer))
                                    for finer in A.finer() if i in finer[0])

    Powersum = P


#################
class OMPBasis_OMPQSym(OMPBasis_abstract):
    """
    Add methods for `OMPQSym` beyond those appearing in ``OMPBasis_abstract``
    """
    def _Coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        .. TODO::

            - poke around for interesting (known) subalgebras
            - describe code added here:
                if R is ???, bring things in, else punt to the default method.
            - add "See :meth:`OMPBasis_abstract._coerce_map_from_` for additional coercions."
        """
        # see fqsym.py: :meth:`FQSymBasis_abstract._coerce_map_from_` for hints
        from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
        M = QuasiSymmetricFunctions(self.base_ring()).Monomial()
        if M.has_coerce_map_from(R):
            print("do something")
            return True
        else:
            return OMPBasis_abstract._coerce_map_from_(self, R)


    class Element(OMPBasis_abstract.Element):
        def is_symmetric(self):
            r"""
            Meant to model the inclusion of SYM inside QSYM, though there is
            no notion of quasisymmetric polynomials related to
            Dual to Hopf Algebra on Ordered Multiset Partitions
            as of yet.

            Determine if a `OMPQSym` function, expressed in the
            `\mathbf{M}` basis, is symmetric.

            A function `f` in the `\mathbf{M}` basis shall be deemed *symmetric*
            if for each ordered multiset partition `A` in its support, all derangements of `A`
            are also in the support of `f` with the same coefficient.

            .. MATH::

                f = \sum_{\lambda} c_{\lambda}
                        \sum_{\lambda(A) = \lambda} \mathbf{M}_A,

            where the first sum is over unordered multiset partitions into sets,
            and the second sum is over all ordered multiset partitions `A` with
            a common multiset of blocks (equal to `\lambda`).

            OUTPUT:

            - ``True`` if whenever ``sorted(A) == sorted(B)`` for `A,B` in
              the support of ``M(self)``, the corresponding coefficients are equal;
              ``False`` otherwise

            EXAMPLES::

                sage: M = OMPQuasiSymmetricFunctions(QQ).M()
                sage: elt = M.sum_of_derangements([[2],[2,3],[1]])
                sage: elt.is_symmetric()
                True
                sage: elt -= 3*M[[1],[1]]
                sage: elt.is_symmetric()
                True
                sage: elt += M([[2],[2,3]])
                sage: elt.is_symmetric()
                False
                sage: elt += M([[2,3],[2]])
                sage: elt.is_symmetric()
                True

            TESTS::

                sage: elt = M[[]]
                sage: elt.is_symmetric()
                True
            """
            OMP = self.parent()._indices
            M = self.parent().realization_of().M()
            d = {}
            for A, coeff in M(self):
                la = OMP(sorted(A))
                if la not in d:
                    d[la] = [coeff, 1]
                else:
                    if d[la][0] != coeff:
                        return False
                    d[la][1] += 1
            # Make sure we've seen each ordered multiset partition in the derangement class
            return all(d[la][1] == Permutations_mset(la).cardinality() for la in d)

        def to_quasisymmetric_function(self):
            r"""
            The projection of ``self`` to the ring of quasisymmetric functions.

            .. TODO::

                check and add to docstring: vsp map?, alg map?, coalg map? Hopf map?

            OUTPUT:

            - an element of the quasisymmetric functions in the monomial basis
            """
            M = QuasiSymmetricFunctions(self.parent().base_ring()).Monomial()
            qM = self.parent().realization_of().M()
            terms = []
            for (omp, coeff) in qM(self):
                if self._order_grading:
                    I = omp.shape_from_cardinality()
                    coeff = coeff * zee(I)
                else:
                    I = omp.shape_from_size()
                    coeff = coeff * zee(I)
                terms.append((I, coeff))
            return M.sum_of_terms(terms)


############################
class OMPQuasiSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The Hopf algebra that is graded dual to Free Hopf Algebra on Finite Sets (OMPSym).

    .. TODO: modify/correct documentation (simply copied from OMPNonCommutativeSymmetricFunctions)

    The Hopf algebra OMPSym is a free algebra built on finite subsets
    of positive integers. Its coproduct is defined on subsets `K` via

    .. MATH::

        \Delta(K) = \sum_{I\sqcup J = K} I \otimes J,

    then extended multiplicatively and linearly. See [LM]_.

    The product and coproduct are analogous to the complete/homogeneous
    basis for the Hopf algebras :class:`NonCommutativeSymmetricFunctions`
    and :class:`SymmetricFunctions`, so we use `H` as prefix for this
    natural basis. E.g., for two subsets `K` and `L`, we write

    .. MATH::

        H_{[K]} \cdot H_{[L]} = H_{[K, L]}.

    This Hopf algebra is implemented in the Homogeneous basis as a
    :class:`CombinatorialFreeModule` whose basis elements are indexed
    by *ordered multiset partitions*, which are simply lists of subsets
    of positive integers.

    There are two natural gradings on OMP: grading by the *size* or by
    the *order* of ordered multiset partitions. Only the first of these
    yields finite dimensional slices if a finite alphabet is not also
    provided. Hence, we allow an alphabet to be passed as an optional
    argument. (The user may also elect to use the size grading, even if
    a finite alphabet is passed.)

    .. SEEALSO::

        :class:`OMPNonCommutativeSymmetricFunctions`
        :class:`OrderedMultisetPartitionsIntoSets`

    EXAMPLES:

    Standard use of Free Hopf Algebra on Finite Sets:

    We begin by first creating the ring of `OMPSym` and the bases that are
    analogues of the usual non-commutative symmetric functions::

        sage: OMPSym = OMPNonCommutativeSymmetricFunctions(QQ)
        sage: H = OMPSym.H()
        sage: P = OMPSym.P() # ???
        sage: H
        Free Hopf Algebra on Finite Sets over Rational Field
         in the Homogeneous basis

    The basis is indexed by ordered multiset partitions, so we create an
    element and convert it between these bases::

        sage: elt = H(OrderedMultisetPartitionIntoSets([[2], [1,3]])) - 2*H(OrderedMultisetPartitionIntoSets([[1,2,3]])); elt
        -2*H[{1,2,3}] + H[{2}, {1,3}]
        sage: P(elt)
        ???

    There is also a shorthand for creating elements. We note that we must
    use ``P[[]]`` to create the empty ordered multiset partition due to
    python's syntax. ::

        sage: elth = H[[2], [1,3]] - 2*H[[1,2,3]]; elth
        -2*H[{1,2,3}] + H[{2}, {1,3}]
        sage: eltp = P[[]] + elth; eltp
        ???

    Restricting the alphabet and alternate grading:

    Instead of working over ordered multiset partitions of positive integers,
    we may restrict each block's elements to come from a finite alphabet `A`::

        sage: OMPSymA = OMPNonCommutativeSymmetricFunctions(QQ, alphabet=[3,4,5])
        sage: HA = OMPSymA.HA()
        sage: HA.an_element()
        ???
        sage: H.an_element()
        ???
        sage: H.an_element() in OMPSymA
        False???
        sage: HA.an_element() in OMPSym
        True??? # I'm happy with either answer, I suppose

    The homogeneous component of degree `d` of the basis is indexed by
    ordered multiset partitions over `A` of order `d`::

        sage: HA.basis(3).keys().list()
        ???

    If the keyword argument ``order_grading`` is set to ``False`` (or, if the
    ``alphabet`` keyword is not used), then homogeneous degree is determined
    by size of ordered multiset partitions::

        sage: OMPSymB = OMPNonCommutativeSymmetricFunctions(QQ, alphabet=[3,4,5], order_grading=False)
        sage: HB = OMPSymB.HB()
        sage: HB.basis(12).keys().list()
        ???

    TESTS::

        sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
        sage: a = H[OrderedMultisetPartitionIntoSets([[1]])]
        sage: b = H[OrderedMultisetPartitionsIntoSets(1)([[1]])]
        sage: c = H[[1]]
        sage: a == b == c
        True
        sage: H = OMPNonCommutativeSymmetricFunctions(QQ, alphabet=[2]).H()
        sage: a = H[OrderedMultisetPartitionIntoSets([[2]])]
        sage: b = H[OrderedMultisetPartitionsIntoSets(2)([[2]])]
        sage: c = H[[2]]
        sage: d = H[OrderedMultisetPartitionsIntoSets([1,2], 1)([[2]])]
        sage: a == b == c == d
        True

    .. TODO::

        - Add dual to the powersum basis of OMPSym.
        - Add/correct maps to/from NCSym, NCQSym, Sym, QSym.
    """
    @staticmethod
    def __classcall_private__(cls, R, alphabet=None, order_grading=None):
        """
        Normalize passed arguments.
        """
        # ensure R is a commutative ring
        if R not in CommutativeRings():
            raise ValueError("argument R (=%s) must be a commutative ring" % repr(R))

        # make alphabet hashable
        if alphabet is not None:
            if alphabet in ZZ:
                _A = Set(range(1,ZZ(alphabet)+1))
            else:
                _A = Set(ZZ(a) for a in alphabet)
            if _A == Set() or any(a not in ZZ for a in _A):
                raise ValueError("keyword alphabet was converted to %s, which must be a nonempty set of positive integers"%(self._A))
        else:
            # treat as ``_A`` is PositiveIntegers
            _A = None

        # pick an appropriate value for grading
        if order_grading is None:
            if _A:
                _order_grading = True
            else:
                _order_grading = False
        else:
            _order_grading = order_grading

        return super(OMPQuasiSymmetricFunctions, cls).__classcall__(cls, R, _A, _order_grading)

    def __init__(self, R, alphabet, order_grading):
        """
        Initialize ``self``.

        INPUT:

        - ``R`` -- a commutative ring
        - ``alphabet`` -- a finite set of positive integers (value ``None``
          is used to represent all positive integers)
        - ``order_grading`` -- boolean: ``False`` indicates grading by size,
          while ``True`` indicates grading by order

        TESTS::

            sage: A = OMPQuasiSymmetricFunctions(QQ)
            sage: TestSuite(A).run()  # long time
        """
        category = GradedHopfAlgebras(R).Connected().Commutative()
        Parent.__init__(self, base=R, category=category.WithRealizations())
        self._A = alphabet
        self._order_grading = order_grading


    def _repr_(self):
        r"""
        EXAMPLES::

            sage: OMPQuasiSymmetricFunctions(ZZ)
            Dual of Free Hopf Algebra on Finite Sets over Integer Ring
        """
        return "Dual to " + repr(self.dual())

    def a_realization(self):
        r"""
        Return the realization of the `\mathbf{M}` basis of ``self``.

        EXAMPLES::

            sage: OMPQuasiSymmetricFunctions(QQ).a_realization()
            Dual of Free Hopf Algebra on Finite Sets over the Rational Field in the Monomial basis
        """
        return self.M()

    _shorthands = tuple(['M'])

    def dual(self):
        r"""
        Return the Free Hopf Algebra on Finite Sets.

        EXAMPLES::

            sage: OMPD = OMPQuasiSymmetricFunctions(QQ)
            sage: OMPD.dual()
            Free Hopf Algebra on Finite Sets over the Rational Field
        """
        from . import OMPNonCommutativeSymmetricFunctions
        return OMPNonCommutativeSymmetricFunctions(self.base_ring(), self._A, self._order_grading)

    class M(OMPBasis_OMPQSym):
        r"""
        The Monomial basis of Dual to Hopf Algebra on Ordered Multiset Partitions.


        The family `(M_\mu)`, as `\mu` ranges over all ordered multiset
        partitions into sets, is the basis that is graded dual to the
        Homogeneous basis of Hopf Algebra on Ordered Multiset Partitions.
        By analogy with :class:`QuasiSymmetricFunctions`, it is called the
        *Monomial basis* here.

        .. TODO: expand the docstring (mention quasi shuffle and deconcatenation)

        EXAMPLES::

            sage: M = OMPQuasiSymmetricFunctions(QQ).M()
            sage: M[[2], [1,3]] - 2*M[[1,2,4]]
            M[{2}, {1,3}] - 2*M[{1,2,4}]

            sage: M[[2], [1,3]] * M[[2,4]]
            M[{2},{1,3},{2,4}] + M[{2},{1,2,3,4}] + M[{2},{2,4},{1,3}]
             + M[{2,4},{2},{1,3}]
            sage: M[[5], [1,3]] * M[[2,4]]
            M[{5},{1,3},{2,4}] + M[{5},{1,2,3,4}] + M[{5},{2,4},{1,3}]
             + M[{2,4,5},{1,3}] + M[{2,4},{5},{1,3}]

            sage: M[[2], [1,3], [1,2,4]].coproduct()
            M[] # M[{2},{1,3},{1,2,4}] + M[{2}] # M[{1,3},{1,2,4}]
             + M[{2},{1,3}] # M[{1,2,4}] + M[{2},{1,3},{1,2,4}] # M[]

        TESTS::

            sage: M = OMPQuasiSymmetricFunctions(QQ).M()
            sage: TestSuite(M).run()  # long time

        .. TODO: Do we need these tests, or are they checked within TestSuite?::

            sage: tests = []
            sage: for ma in M.basis(3):
            ....:     tests.append(M.zero() == sum(M[a1].antipode() * M[a2] * coeff \
            ....:           for ((a1,a2),coeff) in (ma.coproduct())))
            ....:     tests.append(M.zero() == sum(M(a1) * M(a2).antipode() * coeff \
            ....:           for ((a1,a2),coeff) in (ma.coproduct())))
            ....:     for mb in [M[[]]] + list(M.basis(1)) + list(M.basis(2)):
            ....:         tests.append((ma * mb).coproduct() == ma.coproduct() * mb.coproduct())
            sage: all(tests)
        """
        _prefix = "M"
        _basis_name = "Monomial"

        def product_on_basis(self, A, B):
            r"""
            The product on `\mathbf{M}` basis elements.

            The product on the `\mathbf{M}` is the dual to the coproduct on the
            `\mathbf{H}` basis.  On the basis `\mathbf{M}` it is defined as

            .. MATH::

                ???

            where the sum is over all possible ... ???
            This product is commutative.

            INPUT:

            - ``A``, ``B`` -- ordered multiset partitions into sets

            OUTPUT:

            - an element of the `\mathbf{M}` basis

            EXAMPLES::

                sage: M = HopfAlgebraOnOrderedMultisetPartitions(QQ).dual().M()
                sage: A = OrderedMultisetComposition(from_zero_list=[2,1,3,0,1,2]); A
                [{1,2,3}, {1,2}]
                sage: B = OrderedMultisetComposition([[3,4]]); B
                [{3,4}]
                sage: M.product_on_basis(A, B)
                M[{1,2,3}, {1,2}, {3,4}] + M[{1,2,3}, {3,4}, {1,2}] + M[{3,4}, {1,2,3}, {1,2}] + M[{1,2,3}, {1,2,3,4}]
                sage: C = OrderedMultisetComposition([[4,5]]); C
                [{4,5}]
                sage: M.product_on_basis(A, C)
                M[{1,2,3}, {1,2}, {4,5}] + M[{1,2,3}, {4,5}, {1,2}] + M[{4,5}, {1,2,3}, {1,2}] + M[{1,2,3}, {1,2,4,5}] + M[{1,2,3,4,5}, {1,2}]
                sage: M.product_on_basis(A, OrderedMultisetComposition([]))
                M[{1,2,3}, {1,2}]
            """
            terms = A.shuffle_product(B, overlap=True)
            return self.sum_of_terms([(s, 1) for s in terms])

        def coproduct_on_basis(self, A):
            r"""
            Return the coproduct of a `\mathbf{M}` basis element.

            The coproduct on the basis element `\mathbf{M}_A` is the sum over
            tensor product terms `\mathbf{M}_B \otimes \mathbf{M}_C` where
            `A=B+C`. (Here `+` is concatenation of ordered multiset partitions.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - The coproduct applied to the `\mathbf{M}` basis element indexed by ``A``
              expressed in the `\mathbf{M}` basis.

            EXAMPLES::

                sage: M = HopfAlgebraOnOrderedMultisetPartitions(QQ).dual().M()
                sage: M[[2], [2,3]].coproduct()
                M[] # M[{2}, {2,3}] + M[{2}] # M[{2,3}] + M[{2}, {2,3}] # M[]
                sage: M.coproduct_on_basis(OrderedMultisetComposition([]))
                M[] # M[]
            """
            return self.tensor_square().sum_of_terms(A.deconcatenate(2), distinct=True)

        def _Antipode_on_basis(self, A):
            r"""
            Return the antipode applied to the basis element indexed by ``A``.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - an element in the basis ``self``
            """
            pass

        def dual_basis(self):
            r"""
            Return the dual basis to the `\mathbf{M}` basis.

            The dual basis to the `\mathbf{M}` basis is the Homogeneous basis
            of Free Hopf Algebra on Finite Sets.

            OUTPUT:

            - the Homogeneous basis of Free Hopf Algebra on Finite Sets

            EXAMPLES::

                sage: M = OMPQuasiSymmetricFunctions(QQ).M()
                sage: M.dual_basis()
                Free Hopf Algebra on Finite Sets over the Rational Field in the Homogeneous basis
            """
            return self.realization_of().dual().H()

        def sum_of_derangements(self, A):
            """
            Return the sum, in the Monomial basis, over all ordered multiset
            partitions whose blocks coincide with those of ``A``.

            .. TODO:

                - investigate if this behaves like `Sym` inside `QSym`.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - an element of ``self``

            EXAMPLES::

                sage: M = OMPQuasiSymmetricFunctions(QQ).M()
                sage: M.sum_of_derangements([[2,1],[1]])
                M[{1},{1,2}] + M[{1,2},{1}]
            """
            OMP = self._indices
            return self.sum_of_terms([(OMP(pi),1) \
                        for pi in Permutations_mset(list(A))], distinct=True)


        class Element(CombinatorialFreeModule.Element):
            r"""
            An element in the `\mathbf{M}` basis of OMP^*.
            """

####### Auxillary Functions #################
def zee(x):
    if not x in Partitions():
        x = Partitions()(sorted(x, reverse=True))
    return x.centralizer_size()


##############
OMPSym = OMPNonCommutativeSymmetricFunctions(QQ)
H = OMPSym.H()
H