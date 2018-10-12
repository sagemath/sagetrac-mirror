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
from sage.misc.misc_c import prod
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.cartesian_product import cartesian_product
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.fields import Fields

from sage.matrix.constructor import matrix
from sage.sets.set import Set
from sage.rings.all import ZZ

from sage.functions.other import factorial

from sage.combinat.posets.posets import Poset
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.sf.sfa import zee
from sage.combinat.permutation import Permutations_mset
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.set_partition import SetPartitions_set
from sage.combinat.multiset_partition_into_sets_ordered import OrderedMultisetPartitionsIntoSets

class OMPBases(Category_realization_of_parent):
    r"""
    The category of bases of `OMPNSym` and `OMPQSym`.
    """
    def __init__(self, base, graded):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base`` -- an instance of `OMPNSym` or `OMPQSym`
        - ``graded`` -- boolean; if the basis is graded or filtered

        TESTS::

            sage: from sage.combinat.chas.omp_hopf_algebras import OMPBases
            sage: OMPNSym = OMPNonCommutativeSymmetricFunctions(ZZ)
            sage: bases = OMPBases(OMPNSym, True)
            sage: OMPNSym.H() in bases
            True
        """
        self._graded = graded
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.omp_hopf_algebras import OMPBases
            sage: OMPNSym = OMPNonCommutativeSymmetricFunctions(ZZ)
            sage: OMPBases(OMPNSym, True)
            Category of graded bases of Free Hopf Algebra on Finite Sets over the Integer Ring
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
        if self._graded:
            cat = cat.Graded()
        else:
            cat = cat.Filtered()
        return [self.base().Realizations(),
                HopfAlgebras(R).Graded().Realizations(),
                cat.Connected()]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of `OMPNSym` or `OMPQSym`.

            EXAMPLES::

                sage: OMPNSym = OMPNonCommutativeSymmetricFunctions(ZZ)
                sage: OMPNSym.H()
                Free Hopf Algebra on Finite Sets over the Integer Ring in the Homogeneous basis
                sage: OMPQSymA = OMPQuasiSymmetricFunctions(QQ, alphabet=[2,3,4])
                sage: OMPQSymA.M()
                Dual of Free Hopf Algebra on Finite Sets on alphabet {2, 3, 4} over the Rational Field with order grading in the Monomial basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def __getitem__(self, A):
            """
            Return the basis element indexed by ``A``.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

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

            - ``self`` -- an element of `OMPNSym` or `OMPQSym`.
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
                [ 1  0  0  0  0]
                [ 0  1  0 -1  0]
                [ 0  0  1  0  0]
                [ 0  0  0  1  0]
                [ 0  0  0  0  1]
            """
            x = self.parent().dual_basis()(x)
            return sum(coeff * x[I] for (I, coeff) in self)

###### Common Basis Methods for OMPNSym & OMPQSym ############
class OMPBasis_abstract(CombinatorialFreeModule, BindableClass):
    """
    Abstract base class for bases of `OMPNSym` and `OMPQSym`.

    This must define two attributes:

    - ``_prefix`` -- the basis prefix
    - ``_basis_name`` -- the name of the basis (must match one
      of the names that the basis can be constructed from OMPNSym or OMPQSym)
    """
    def __init__(self, alg, graded=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: OMPNonCommutativeSymmetricFunctions(ZZ).H()
            Free Hopf Algebra on Finite Sets over the Integer Ring in the Homogeneous basis

            sage: OMPQuasiSymmetricFunctions(ZZ, alphabet=[2,3,4]).M()
            Dual of Free Hopf Algebra on Finite Sets on alphabet {2, 3, 4} over the Integer Ring with order grading in the Monomial basis
            sage: OMPQuasiSymmetricFunctions(ZZ, alphabet=[2,3,4], order_grading=False).M()
            Dual of Free Hopf Algebra on Finite Sets on alphabet {2, 3, 4} over the Integer Ring in the Monomial basis

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
        Return a string representation of an element of `OMPNSym` or `OMPQSym`
        in the basis ``self``.

        TESTS::

            sage: M = OMPQuasiSymmetricFunctions(QQ).M()
            sage: elt = M[[1,2]]*M[[1]]; elt
            M[{1}, {1,2}] + M[{1,2}, {1}]
        """
        return self._prefix + repr(x) #.replace(", ", ",")

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are elements of a parent
        of similar type (`OMPNSym` or `OMPQSym`) over a base with a
        coercion map into ``self.base_ring()``.

        TODO:

            - consider being more accommodating with coercions:
              + suppose x is an element of H1 and H1._A is not a subset of H2._A.
              + if each omp in support of x is over an alphabet contained in H2._A,
              + then allow H2(x) instead of throwing a TypeError.

            - add similar method in :class:`rCompBasis_OMPNSym` that also absorbs NSYM.
            - add similar method in :class:`rCompBasis_OMPQSym` that also absorbs NCQSYM??

        EXAMPLES::

            sage: H = OMPNonCommutativeSymmetricFunctions(GF(7)).H(); H
            Free Hopf Algebra on Finite Sets over the Finite Field
             of size 7 in the Homogeneous basis

        Elements of the Homogeneous basis of OMPNSym canonically coerce in::

            sage: a, b = H([[1]]), H([[2,3]])
            sage: H.coerce(a+b) == a+b
            True

        Elements of the Homogeneous basis of OMPNSym over `\ZZ` coerce in,
        since `\ZZ` coerces to `\GF{7}`::

            sage: HZ = OMPNonCommutativeSymmetricFunctions(ZZ).H(); HZ
            Free Hopf Algebra on Finite Sets over the Integer Ring in the Homogeneous basis
            sage: aZ, bZ = HZ([[1]]), HZ([[2,3]])
            sage: c = H.coerce(aZ+bZ); c
            H[{2,3}] + H[{1}]
            sage: c.parent() is H
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so `OMPNSym` over `\GF{7}`
        does not coerce to the same algebra over `\ZZ`::

            sage: HZ.coerce(b)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Free Hopf Algebra on Finite Sets
             over the Finite Field of size 7 in the Homogeneous basis to
             Free Hopf Algebra on Finite Sets over the Integer Ring in
             the Homogeneous basis

        Elements of the Homogeneous basis of OMPNSym over the same base ring
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
            sage: a
            H[{1}]
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
            2*H[{1}, {1,2}] + H[{4}]
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
            [M[], M[{1}], 1/2*M[] + 2*M[{1}, {1,2}] + M[{4}], 2*M[{1}, {1}],
             M[{1}, {4}] + 2*M[{1}, {1}, {1,2}] + 2*M[{1}, {1,2}, {1}] + M[{4}, {1}]]
        """
        u = self.one()
        if self._A:
            o = self([[self._A[0]]])
        else:
            o = self([[1]])
        x = o.leading_support()
        s = self.base_ring().an_element()
        a = self.an_element()
        fake_commutator = self.sum_of_terms([(x+y, c) for (y,c) in a]) \
                    - self.sum_of_terms([(y+x, -c) for (y,c) in a])
        return [u, o, s + a, o*o, fake_commutator]



###### Basis Methods for OMPNSym ############
class OMPBasis_OMPNSym(OMPBasis_abstract):
    """
    Add methods for `OMPNSym` beyond those appearing in ``OMPBasis_abstract``

    .. TODO::

        - Consider adding maps from `NSym` to `OMPNSym`. Any powersum basis
          there could be mapped to certain Powersum elements in `OMPNSym`.
        - If so, then consider updating the inherited method ``_coerce_map_from_``.
        - Consider adding maps to/from NCSym using powersums ... not entirely natural.
    """
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
        if self.base_ring().is_zero():
            return True
        if self._A is not None and len(self._A) == 1:
            return True
        return False

    def primitive(self, K, i=None, within_self=True):
        r"""
        Return a primitive associated to ``K`` in ``self``.

        Suppose `K` is a nonempty set with distinguished element `i`.
        Return the primitive

        .. MATH::

            p_i(K) = \sum_{B} (-1)^{\ell(B)-1} H_B,

        where the sum is over all ordered set partitions `B` that are finer
        than `[K]` and that have `i` in their first block. If ``within_self``
        is ``True``, return this primitive as an element of ``self``.

        INPUT:

        - ``K`` -- a nonempty set
        - ``i`` -- (default: min(K)) distinguished element in ``K``
          around which the primitive is built.

        OUTPUT:

        - an element in the basis ``self``

        EXAMPLES::

            sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
            sage: P = OMPNonCommutativeSymmetricFunctions(QQ).Powersum()
            sage: H.primitive({2,3,4})
            H[{2}, {4}, {3}] + H[{2}, {3}, {4}] - H[{2}, {3,4}]
             - H[{2,3}, {4}] - H[{2,4}, {3}] + H[{2,3,4}]
            sage: _ == H(P[[2,3,4]])
            True
            sage: elt = H.primitive({2,3,4}, i=3); elt
            H[{3}, {4}, {2}] + H[{3}, {2}, {4}] - H[{3}, {2,4}]
             - H[{2,3}, {4}] - H[{3,4}, {2}] + H[{2,3,4}]
            sage: elt.coproduct() - elt.tensor(H.one()) - H.one().tensor(elt)
            0

            sage: P.primitive({2,3,4}, within_self=False)
            H[{2}, {4}, {3}] + H[{2}, {3}, {4}] - H[{2}, {3,4}]
             - H[{2,3}, {4}] - H[{2,4}, {3}] + H[{2,3,4}]
            sage: P(_) == P.primitive({2,3,4}) == P[[2,3,4]]
            True
            sage: P.primitive({2,3,4}).coproduct()
            P[{2,3,4}] # P[] + P[] # P[{2,3,4}]
        """
        if [K] not in self._indices:
            raise ValueError("{0} must be a nonempty set compatible with {1}".format(K, self._indices))
        if i is None:
            i = min(K)

        # build a primitive in H
        H = self.realization_of().H()
        A = self._indices([K])
        primitive = H.sum_of_terms([(B, (-1)**(1+len(B))) for B in A.finer() if i in B[0]])

        # return in the given basis
        if within_self:
            primitive = self(primitive)
        return primitive

    class Element(OMPBasis_abstract.Element):
        def to_symmetric_function(self):
            r"""
            Return the image of ``self`` in the homogeneous basis of the
            ring of symmetric functions.

            There is a natural Hopf algebra map from the `(H_\mu)` basis of
            Free Hopf Algebra on Finite Sets to the `(h_\lambda)` basis of the
            ring of symmetric functions. It suffices to describe this map on
            free algebra generators (then extend multiplicatively).

            Given an ordered multiset partition `\mu` with *one block* of
            cardinality `k`, the map is given by

            .. MATH::

                H_{\mu} \mapsto k! h_k.

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
                sage: h = SymmetricFunctions(QQ).homogeneous()
                sage: (H[[2,3]] * H[[3,4,5]]).to_symmetric_function()
                12*h[3, 2]
                sage: H[[2,3]].to_symmetric_function() * H[[3,4,5]].to_symmetric_function()
                12*h[3, 2]

                sage: D_Phi = H[[2,3],[3,4]].to_symmetric_function().coproduct()
                sage: HH = H[[2,3],[3,4]].coproduct()
                sage: Phi_D = HH.apply_multilinear_morphism(lambda x,y: x.to_symmetric_function().tensor(y.to_symmetric_function()))
                sage: D_Phi == Phi_D
                True
            """
            return self.to_noncommutative_symmetric_function().to_symmetric_function()

        def to_noncommutative_symmetric_function(self):
            r"""
            Return the image of ``self`` in the Complete basis of the
            ring of noncommutative symmetric functions.

            There is a natural Hopf algebra map from the `(H_\mu)` basis of
            Free Hopf Algebra on Finite Sets to the `(S_\alpha)` basis of
            the ring of noncommutative symmetric functions. It suffices to
            describe this map on free algebra generators (then extend
            multiplicatively).

            Given an ordered multiset partition `\mu` with *one block* of
            cardinality `k`, the map is given by

            .. MATH::

                H_{\mu} \mapsto k! S_k.

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
                sage: S = NonCommutativeSymmetricFunctions(QQ).Complete()
                sage: (H[[2,3,5]] * H[[3,4]]).to_noncommutative_symmetric_function()
                12*S[3, 2]
                sage: H[[2,3,5]].to_noncommutative_symmetric_function() * H[[3,4]].to_noncommutative_symmetric_function()
                12*S[3, 2]

                sage: D_Phi = H[[2,3],[3,4]].to_noncommutative_symmetric_function().coproduct()
                sage: HH = H[[2,3],[3,4]].coproduct()
                sage: Phi_D = HH.apply_multilinear_morphism(lambda x,y: x.to_noncommutative_symmetric_function().tensor(y.to_noncommutative_symmetric_function()))
                sage: D_Phi == Phi_D
                True
            """
            from sage.combinat.ncsf_qsym.ncsf import NonCommutativeSymmetricFunctions
            S = NonCommutativeSymmetricFunctions(self.parent().base_ring()).Complete()
            H = self.parent().realization_of().H()
            return S.sum_of_terms((A.shape_from_cardinality(),
                                    coeff*prod(factorial(len(a)) for a in A))
                                        for (A, coeff) in H(self))

###### The Hopf algebra OMPNSym ############
class OMPNonCommutativeSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    Free Hopf Algebra on Finite Sets (OMPNSym).

    INPUT:

    - ``R`` -- a base ring for the algebra

    OPTIONAL KEYWORDS:

    - ``alphabet`` -- a list of positive integers or an integer `n`.
      In the latter case, the alphabet is taken to be `\{1,\ldots,\}`.
      If alphabet is passed, then basis keys are restricted to those
      ordered multiset partitions into sets supported by ``alphabet``.

    - ``order_grading`` -- a boolean that defaults to ``True`` if an
      alphabet is passed, and to ``False`` otherwise. Determines the
      grading function to be used for `OMPQSym`. Let `A` be an ordered
      multiset partition into sets. If this keyword is ``True``, then
      use ``A.order()`` as a degree function; else use ``A.size()``.

    The Hopf algebra `OMPNSym` is a free algebra built on finite subsets
    of positive integers. Its coproduct is defined on subsets `K` via

    .. MATH::

        \Delta(K) = \sum_{I\sqcup J = K} I \otimes J,

    then extended multiplicatively and linearly. See [LM2018]_.

    The product and coproduct are analogous to the complete (homogeneous)
    basis for the Hopf algebras :class:`NonCommutativeSymmetricFunctions`
    and :class:`SymmetricFunctions`, so we use `H` as prefix for this
    natural basis. E.g., for two subsets `K` and `L`, we write

    .. MATH::

        H_{[K]} \cdot H_{[L]} = H_{[K, L]}.

    This Hopf algebra is implemented in the Homogeneous basis as a
    :class:`CombinatorialFreeModule` whose basis elements are indexed
    by *ordered multiset partitions into sets*, which are simply lists
    of subsets of positive integers.

    There are two natural gradings on `OMPNSym`: grading by the *size* or
    by the *order* of ordered multiset partitions. Only the first of these
    yields finite dimensional slices if a finite alphabet is not also
    provided. Hence, we allow an alphabet to be passed as an optional
    keyword argument. (The user may also elect to use the size grading,
    even if a finite alphabet is passed.)

    In the case that graded pieces are finite dimensional, then one may
    speak of the graded dual Hopf algebra `OMPQSym`.

    .. SEEALSO::

        :class:`OMPQuasiSymmetricFunctions`
        :class:`OrderedMultisetPartitionsIntoSets`

    REFERENCES:

    - [LM2018]_


    .. RUBRIC:: Standard use of Free Hopf Algebra on Finite Sets

    We begin by first creating the ring of `OMPNSym` and two bases that
    are analogues of the usual non-commutative symmetric functions::

        sage: OMPNSym = OMPNonCommutativeSymmetricFunctions(QQ)
        sage: H = OMPNSym.H()
        sage: P = OMPNSym.P()
        sage: H
        Free Hopf Algebra on Finite Sets over the Rational Field
         in the Homogeneous basis

    The homogeneous component of degree `d` of any basis is indexed by
    ordered multiset partitions into sets of size `d`::

        sage: H.basis(3).keys()
        Ordered Multiset Partitions into Sets of integer 3

    So we create an element and convert it between these bases::

        sage: OMP = OrderedMultisetPartitionIntoSets
        sage: elt = H(OMP([[2], [1,3]])) - 2*H(OMP([[1,2,3]])); elt
        H[{2}, {1,3}] - 2*H[{1,2,3}]
        sage: P(elt)
        -2*P[{1}, {2,3}] + P[{2}, {1}, {3}] + P[{2}, {1,3}] - 2*P[{1}, {2}, {3}]
         - 2*P[{1,2}, {3}] - 2*P[{1,3}, {2}] - 2*P[{1,2,3}]

    Given an element of `OMPNSym` expressed in two different bases, the two
    expressions are treated as equal::

        sage: P(elt) == elt
        True
        sage: dict(P(elt)) == dict(elt)
        False

    There is also a shorthand for creating elements. We note that we must
    use ``H[[]]`` to create the empty ordered multiset partition due to
    python's syntax::

        sage: elth = H[[2], [1,3]] - 2*H[[1,2,3]] + H[[]]; elth
        H[] + H[{2}, {1,3}] - 2*H[{1,2,3}]

    Finally, expressions are returned in only one basis, even if given
    in mixed bases::

        sage: eltp = 2*P[[1,2,3]] + elth; eltp
        -2*P[{1}, {2}, {3}] - 2*P[{1}, {2,3}] - 2*P[{1,2}, {3}]
         - 2*P[{1,3}, {2}] + P[] + P[{2}, {1}, {3}] + P[{2}, {1,3}]

        sage: elth * P[[1],[2],[3]]
        H[{1}, {2}, {3}] + H[{2}, {1,3}, {1}, {2}, {3}]
         - 2*H[{1,2,3}, {1}, {2}, {3}]

    .. RUBRIC:: Restricting the alphabet and using the order grading

    When the keyword ``alphabet`` is used to select an alphabet `A`, then
    the homogeneous component of degree `d` of the basis is indexed by
    ordered multiset partitions over `A` of order `d`::

        sage: OMPNSymA = OMPNonCommutativeSymmetricFunctions(QQ, alphabet=[3,4,5])
        sage: HA = OMPNSymA.H()
        sage: HA.basis(2).keys()
        Ordered Multiset Partitions into Sets of order 2 over alphabet {3, 4, 5}

    The keyword argument ``order_grading`` may be used to override the default
    grading::

        sage: OMPNSym_order = OMPNonCommutativeSymmetricFunctions(QQ, order_grading=True)
        sage: H_order = OMPNSym_order.H()
        sage: H_order.basis(2).keys()
        #??? wrong!!! Ordered Multiset Partitions into Sets of integer 2

        sage: OMPNSymA_size = OMPNonCommutativeSymmetricFunctions(QQ, alphabet=[3,4,5], order_grading=False)
        sage: H_size = OMPNSymA_size.H()
        sage: H_size.basis(2).keys()
        #??? wrong!!! Ordered Multiset Partitions into Sets of order 2 over alphabet {3, 4, 5}

    *Warning:* Realizations of `OMPNSym` such as ``OMPNSym_order`` above will
    have infinite dimensional spaces for their graded slices `OMPNSym_d`.

    .. RUBRIC:: Moving between realizations of `OMPNSym`

    If parents' alphabets are compatible, coerce elements between different
    realizations of `OMPNSym`::

        sage: HA.an_element()
        2*H[{3}, {3}, {4,5}] + H[{3}, {3,4,5}]
        sage: H.an_element()
        2*H[{1}, {1,2}] + H[{4}]
        sage: H(HA.an_element()) in H
        True
        sage: HA(H.an_element()) in HA
        Traceback (most recent call last):
        ...
        TypeError: do not know how to make x (= 2*H[{1}, {1,2}] + H[{4}])
         an element of self (=Free Hopf Algebra on Finite Sets on alphabet
         {3, 4, 5} over the Rational Field with order grading in the Homogeneous basis)

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

            - sort out the degree issues above. (search for "wrong!!!")
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
                raise ValueError("keyword alphabet (=%s) must be a nonempty set of positive integers"%(_A))
        else:
            # treat ``_A`` as PositiveIntegers
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
        category = HopfAlgebras(R).Graded().Connected()  # TODO: add Cocommutative
        Parent.__init__(self, base=R, category=category.WithRealizations())
        self._A = alphabet
        self._order_grading = order_grading

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: OMPNonCommutativeSymmetricFunctions(ZZ)
            Free Hopf Algebra on Finite Sets over the Integer Ring
            sage: OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2,3])
            Free Hopf Algebra on Finite Sets on alphabet {2, 3} over
             the Integer Ring with order grading
            sage: OMPNonCommutativeSymmetricFunctions(ZZ, alphabet=[2,3], order_grading=False)
            Free Hopf Algebra on Finite Sets on alphabet {2, 3} over the Integer Ring
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
            Free Hopf Algebra on Finite Sets over the Rational Field in the Homogeneous basis
        """
        return self.H()

    _shorthands = tuple(['H', 'P', 'R'])

    def dual(self):
        r"""
        Return the Hopf algebra that is the graded dual of
        Free Hopf Algebra on Finite Sets.

        EXAMPLES::

            sage: OMPNonCommutativeSymmetricFunctions(QQ).dual()
            Dual of Free Hopf Algebra on Finite Sets over the Rational Field
        """
        return OMPQuasiSymmetricFunctions(self.base_ring(), self._A, self._order_grading)

    class H(OMPBasis_OMPNSym):
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
            H[{2}, {1,3}] - 2*H[{1,2,4}]

            sage: H[[2], [1,3]] * H[[1,2],[4]]
            H[{2}, {1,3}, {1,2}, {4}]

            sage: H[[4]].coproduct()
            H[] # H[{4}] + H[{4}] # H[]
            sage: H[[1,2]].coproduct()
            H[] # H[{1,2}] + H[{1}] # H[{2}] + H[{2}] # H[{1}] + H[{1,2}] # H[]
            sage: H[[4],[1,2]].coproduct()
            H[] # H[{4}, {1,2}] + H[{4}] # H[{1,2}]
             + H[{4}, {1}] # H[{2}] + H[{4}, {2}] # H[{1}]
             + H[{4}, {1,2}] # H[] + H[{2}] # H[{4}, {1}]
             + H[{1}] # H[{4}, {2}] + H[{1,2}] # H[{4}]

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
            True
        """
        _prefix = "H"
        _basis_name = "Homogeneous"

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

        def product_on_basis(self, A, B):
            r"""
            Return the (associative) `*` product of the basis elements
            of ``self`` indexed by the ordered multiset partitions
            `A` and `B`.

            INPUT:

            - ``A, B`` -- two ordered multiset partitions into sets

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

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - The coproduct applied to the element of OMP indexed by ``A`` expressed in the
              Homogeneous basis.

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
                sage: H[{2,3}].coproduct()
                H[] # H[{2,3}] + H[{2}] # H[{3}] + H[{3}] # H[{2}] + H[{2,3}] # H[]
                sage: H[{2,3}, {2}].coproduct()
                H[] # H[{2,3}, {2}] + H[{2}, {2}] # H[{3}]
                 + H[{3}] # H[{2}, {2}] + H[{2}] # H[{3}, {2}]
                 + H[{3}, {2}] # H[{2}] + H[{2}] # H[{2,3}]
                 + H[{2,3}] # H[{2}] + H[{2,3}, {2}] # H[]
            """
            return self.tensor_square()._from_dict( A.split_blocks(2) )

        def antipode_on_basis(self, A):
            r"""
            Return the antipode applied to a Homogeneous basis element.

            If `A` is an ordered multiset partition with a single block then
            the antipode is given by

            .. MATH::

                S(H_A) = \sum_{B \prec A} (-1)^{\ell(B)} H_B,

            where the sum is over all refinements of `A`.

            If `A` has more than one block, then we extend using the fact that
            the antipode is an algebra antimorphism.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            EXAMPLES::

                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
                sage: H[{2,3}].antipode()
                H[{3}, {2}] + H[{2}, {3}] - H[{2,3}]
                sage: H[{2,3}, {2}].antipode()
                -H[{2}, {2}, {3}] - H[{2}, {3}, {3}] + H[{2}, {2,3}]

            TESTS::

                sage: H[[]].antipode()
                H[]
                sage: H[[1]].antipode()
                -H[{1}]
                sage: HH = H[[1,2,3],[5],[1,4]].coproduct()
                sage: HH.apply_multilinear_morphism(lambda x,y: x.antipode()*y)
                0
            """
            if not A:
                return self.monomial(A)

            out = []
            P = self.basis().keys()
            AA = [P([a]).finer() for a in A.reversal()]
            for refinement in cartesian_product(AA):
                C = reduce(lambda a,b: a + b, refinement, P([]))
                ell = len(C)
                out.append((C, (-1)**ell))
            return self.sum_of_terms(out)

    Homogeneous = H

    class P(OMPBasis_OMPNSym):
        r"""
        The Powersum basis of Free Hopf Algebra on Finite Sets.

        The family `(P_\mu)`, as `\mu` ranges over all ordered multiset
        partitions into sets, is a multiplicative basis of `OMPNSym` called
        the *Powersum basis* here. It is defined via a unitriangular change
        of basis from the Homogeneous basis as described below.

        The principle feature of the Powersum basis is that for subsets `K`,
        `P_{[K]}` is primitive.

        Given a nonempty subset `K`, we identify a special element `k_0\in K` and take

        .. MATH::

            P_{[K]} = \sum_{\pi in F^*(K)} H_{\pi},

        where the sum is over all refinements of the ordered set partition `[K]`
        that contain `k_0` in the first block. (Here

        EXAMPLES::

            sage: P = OMPNonCommutativeSymmetricFunctions(QQ).Powersum()
            sage: P[[2], [1,3]] - 2*P[[1,2,3]]
            P[{2}, {1,3}] - 2*P[{1,2,3}]
            sage: P[[2,3]].coproduct()
            P[] # P[{2,3}] + P[{2,3}] # P[]
            sage: P[[2,3],[1,2]].coproduct()
            P[{2,3}, {1,2}] # P[] + P[{1,2}] # P[{2,3}]
             + P[{2,3}] # P[{1,2}] + P[] # P[{2,3}, {1,2}]
            sage: P[[2,3]] * P[[1,2]]
            P[{2,3}, {1,2}]
            sage: p = P.an_element(); p
            2*P[{1}, {1,2}] + P[{4}]
            sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
            sage: H(p)
            H[{4}] - 2*H[{1}, {1}, {2}] + 2*H[{1}, {1,2}]

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
            OMPBasis_OMPNSym.__init__(self, alg)

            # Register coercions
            H = self.realization_of().H()
            phi = self.module_morphism(self._P_to_H, codomain=H, unitriangular="upper")
            phi.register_as_coercion()
            phi_inv = H.module_morphism(self._H_to_P, codomain=self, unitriangular="upper")
            phi_inv.register_as_coercion()

        def _P_to_H(self, A):
            """
            Return `P_A` in terms of the Homogeneous basis.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - An element of the Homogeneous basis

            EXAMPLES::

                sage: P = OMPNonCommutativeSymmetricFunctions(QQ).P()
                sage: A = OrderedMultisetPartitionIntoSets([])
                sage: P._P_to_H(A)
                H[]
                sage: A = OrderedMultisetPartitionIntoSets([[1], [3], [3]])
                sage: P._P_to_H(A)
                H[{1}, {3}, {3}]

                sage: B = OrderedMultisetPartitionIntoSets([[1,3], [2]])
                sage: P._P_to_H(B)
                -H[{1}, {3}, {2}] + H[{1,3}, {2}]
                sage: B = OrderedMultisetPartitionIntoSets([[1,3], [1,4]])
                sage: P._P_to_H(B)
                H[{1}, {3}, {1}, {4}] - H[{1}, {3}, {1,4}]
                 - H[{1,3}, {1}, {4}] + H[{1,3}, {1,4}]

                sage: C = OrderedMultisetPartitionIntoSets([[1,3,4]])
                sage: P._P_to_H(C)
                H[{1}, {3}, {4}] + H[{1}, {4}, {3}] - H[{1}, {3,4}]
                 - H[{1,3}, {4}] - H[{1,4}, {3}] + H[{1,3,4}]
            """
            H = self.realization_of().H()
            if not A:
                return H.one()

            # If `A` has one part, make `P_A` primitive.
            if len(A) == 1:
                return H.primitive(A[0], min(A[0]), within_self=False)

            # Else, extend multiplicatively to a basis
            OMP = self._indices
            return prod(H(self([a])) for a in A)

        def _H_to_P(self, A):
            """
            Return `H_A` in terms of the Powersum basis.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - An element of the Powersum basis

            TODO:: improve the examples!

            EXAMPLES::

                sage: P = OMPNonCommutativeSymmetricFunctions(QQ).P()
                sage: A = OrderedMultisetPartitionIntoSets([])
                sage: P._H_to_P(A)
                P[]
                sage: A = OrderedMultisetPartitionIntoSets([[1], [3], [3]])
                sage: P._H_to_P(A)
                P[{1}, {3}, {3}]

                sage: B = OrderedMultisetPartitionIntoSets([[1,3], [2]])
                sage: P._H_to_P(B)
                P[{1}, {3}, {2}] + P[{1,3}, {2}]
                sage: B = OrderedMultisetPartitionIntoSets([[1,3], [1,4]])
                sage: P._H_to_P(B)
                P[{1}, {3}, {1}, {4}] + P[{1}, {3}, {1,4}]
                 + P[{1,3}, {1}, {4}] + P[{1,3}, {1,4}]

                sage: C = OrderedMultisetPartitionIntoSets([[1,3,4]])
                sage: P._H_to_P(C)
                P[{1}, {3}, {4}] + P[{1}, {3,4}] + P[{1,4}, {3}]
                 + P[{1,3}, {4}] + P[{1,3,4}]
            """
            # base cases
            if A._order == A.length():
                return self.monomial(A)

            # sum over the set partitions of each block of A
            OMP = self._indices
            terms = []
            for tup in cartesian_product([SetPartitions_set(block) for block in A]):
                terms.append(sum(map(OMP, tup)))
            return self.sum_of_monomials(terms)

        def dual_basis(self):
            r"""
            Return the dual basis to the `\mathbf{P}` basis.

            The dual basis to the `\mathbf{P}` basis is the dual Powersum
            basis of Dual of Free Hopf Algebra on Finite Sets.

            OUTPUT:

            - the dual Powersum basis of Dual of Free Hopf Algebra on Finite Sets

            EXAMPLES::

                sage: P = OMPNonCommutativeSymmetricFunctions(QQ).P()
                sage: P.dual_basis()
                Dual of Free Hopf Algebra on Finite Sets over the Rational Field in the dual Powersum basis
            """
            return self.realization_of().dual().PowersumDual()

        def product_on_basis(self, A, B):
            r"""
            Return the (associative) `*` product of the basis elements
            of ``self`` indexed by the ordered multiset partitions
            `A` and `B`.

            INPUT:

            - ``A, B`` -- two ordered multiset partitions into sets

            OUTPUT:

            - The basis element of ``self`` indexed by the concatenation
              `A + B`. (This product is non-commutative.)

            EXAMPLES::

                sage: P = OMPNonCommutativeSymmetricFunctions(QQ).Powersum()
                sage: P[{2,3}] * P[{1}, {4,5}]
                P[{2,3}, {1}, {4,5}]

            TESTS::

                sage: OMP = OrderedMultisetPartitionsIntoSets
                sage: one = OMP()([])
                sage: all(P.product_on_basis(one, z) == P(z) == P.basis()[z] for z in [one] + list(OMP(3)))
                True
                sage: all(P.product_on_basis(z, one) == P(z) == P.basis()[z] for z in OMP(3))
                True
                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
                sage: all(P[A] * P[B] == P( H(P[A])*H(P[B]) ) for A in OMP(3) for B in OMP(3))  # indirect doctest
                True
            """
            return self.monomial(A + B)

        def coproduct_on_basis(self, A):
            r"""
            Return the coproduct in the basis `(P_\mu)`.

            Partitions `A` with one block yield primitives, i.e.,

            .. MATH::

                \Delta(P_A) = P_A \otimes P_{[]} + P_{[]} \otimes P_A.

            This formula is extended multiplicatively to all of `(P_\mu)`.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - The coproduct applied to the element of `OMPNSym` indexed by ``A``
              expressed in the Powersum basis.

            EXAMPLES::

                sage: P = OMPNonCommutativeSymmetricFunctions(QQ).Powersum()
                sage: P[{2,3,4}].coproduct()
                P[] # P[{2,3,4}] + P[{2,3,4}] # P[]
                sage: P[{2,3,4}, {1,2}].coproduct()
                P[{2,3,4}, {1,2}] # P[] + P[{1,2}] # P[{2,3,4}]
                 + P[{2,3,4}] # P[{1,2}] + P[] # P[{2,3,4}, {1,2}]

            TESTS::

                sage: OMP = OrderedMultisetPartitionsIntoSets
                sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
                sage: all(P.coproduct_on_basis(A) == P.tensor_square()( H(P[A]).coproduct() ) for A in OMP(3))
                True
            """
            z = self.one().leading_support()
            if not A:
                return self.tensor_square().sum_of_terms([((z,z), 1)])
            if len(A) == 1:
                return self.tensor_square().sum_of_terms([((A,z), 1), ((z,A), 1)])

            else:
                # `P` is a multiplicative basis of primitives...
                # use the rule `\Delta(xy) = \Delta(x) * \Delta(y)`
                X = self._indices
                return reduce(lambda x,y: x*y, [self.monomial(X([a])).coproduct() for a in A], self(z).tensor(self(z)))

        def antipode_on_basis(self, A):
            r"""
            Return the antipode applied to a Powersum basis element.

            If `A` is an ordered multiset partition with a single block then
            the antipode is given by `S(P_A) = -P_A`.

            If `A` has more than one block, then we extend using the fact that
            the antipode is an algebra antimorphism.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            EXAMPLES::

                sage: P = OMPNonCommutativeSymmetricFunctions(QQ).Powersum()
                sage: P[{2,3}].antipode()
                -P[{2,3}]
                sage: P[{2,3}, {2}].antipode()
                P[{2}, {2,3}]

            TESTS::

                sage: P[[]].antipode()
                P[]
                sage: PP = P[[1,2,3],[5],[1,4]].coproduct()
                sage: PP.apply_multilinear_morphism(lambda x,y: x.antipode()*y)
                0
            """
            return (-1)**A.length() * self.monomial(A.reversal())

    Powersum = P

    class R(OMPBasis_OMPNSym):
        r"""
        The Ribbon basis of Free Hopf Algebra on Finite Sets.

        The family `(R_\mu)`, as `\mu` ranges over all ordered multiset
        partitions into sets, is called the *Ribbon basis* here. It
        is defined in relation to the Monomial basis of `OMPNSym` as "sum
        above `\mu` in the strong refinement order."

        EXAMPLES::

            sage: R = OMPNonCommutativeSymmetricFunctions(QQ).Ribbon()
            sage: H = OMPNonCommutativeSymmetricFunctions(QQ).Homogeneous()
            sage: x = R[[1,2],[4,5],[3]]
            sage: H(x)
            H[{1,2}, {4,5}, {3}] - H[{1,2,4,5}, {3}]
            sage: x = R[[1,2],[3,4],[5]]
            sage: H(x)
            -H[{1,2}, {3,4,5}] - H[{1,2,3,4}, {5}]
             + H[{1,2,3,4,5}] + H[{1,2}, {3,4}, {5}]

            sage: R[[2,3],[1],[2,4]] * R[[2,3],[5]]
            R[{2,3}, {1}, {2,4}, {2,3}, {5}]
            sage: R[[2,3],[1],[2,4]] * R[[3,6],[5]]
            R[{2,3}, {1}, {2,4}, {3,6}, {5}]
            sage: R[[2,3],[1],[2,4]] * R[[5,6],[5]]
            R[{2,3}, {1}, {2,4}, {5,6}, {5}] + R[{2,3}, {1}, {2,4,5,6}, {5}]

            sage: R[[2,3]].coproduct()
            R[] # R[{2,3}] + R[{2}] # R[{3}] + R[{3}] # R[{2}] + R[{2,3}] # R[]
            sage: R[[2,3],[2]].coproduct()
            R[] # R[{2,3}, {2}] + R[{2}, {2}] # R[{3}]
             + R[{3}] # R[{2}, {2}] + R[{2}] # R[{3}, {2}]
             + R[{3}, {2}] # R[{2}] + R[{2}] # R[{2,3}]
             + R[{2,3}] # R[{2}] + R[{2,3}, {2}] # R[]

        TESTS::

            sage: TestSuite(R).run()  # long time
            sage: all(R(H(p)) == p for p in R.some_elements())
            True
        """
        _prefix = "R"
        _basis_name = "Ribbon"

        def __init__(self, alg):
            """
            Initialize ``self``.
            """
            OMPBasis_OMPNSym.__init__(self, alg)

            # Register coercions
            H = self.realization_of().H()
            phi = self.module_morphism(self._R_to_H, codomain=H, triangular="lower")
            phi.register_as_coercion()
            #(~phi).register_as_coercion()
            phi_inv = H.module_morphism(self._H_to_R, codomain=self, triangular="lower")
            phi_inv.register_as_coercion()

        def _R_to_H(self, A):
            """
            Return `R_A` in terms of the Homogeneous basis.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - An element of the Homogeneous basis

            EXAMPLES::

                sage: R = OMPNonCommutativeSymmetricFunctions(QQ).Ribbon()
                sage: A = OrderedMultisetPartitionIntoSets([])
                sage: R._R_to_H(A)
                H[]
                sage: A = OrderedMultisetPartitionIntoSets([[4], [3], [3]])
                sage: R._R_to_H(A)
                H[{4}, {3}, {3}]

                sage: B = OrderedMultisetPartitionIntoSets([[1,3], [2]])
                sage: R._R_to_H(B)
                H[{1,3}, {2}]
                sage: B = OrderedMultisetPartitionIntoSets([[1,3], [4]])
                sage: R._R_to_H(B)
                H[{1,3}, {4}] - H[{1,3,4}]

                sage: C = OrderedMultisetPartitionIntoSets([[1,3], [4,5], [6]])
                sage: R._R_to_H(C)
                H[{1,3}, {4,5}, {6}] - H[{1,3}, {4,5,6}]
                 - H[{1,3,4,5}, {6}] + H[{1,3,4,5,6}]
            """
            H = self.realization_of().H()
            PS = Poset({B:B.finer_pred(strong=True) for B in A.fatter(strong=True)})
            return H.sum_of_terms([(B, PS.moebius_function(A, B)) for B in PS])

        def _H_to_R(self, A):
            """
            Return `H_A` in terms of the Ribbon basis.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - An element of the Ribbon basis

            EXAMPLES::

                sage: R = OMPNonCommutativeSymmetricFunctions(QQ).R()
                sage: A = OrderedMultisetPartitionIntoSets([])
                sage: R._H_to_R(A)
                R[]
                sage: A = OrderedMultisetPartitionIntoSets([[1,3], [3]])
                sage: R._H_to_R(A)
                R[{1,3}, {3}]
                sage: A = OrderedMultisetPartitionIntoSets([[1,3], [2,4]])
                sage: R._H_to_R(A)
                R[{1,3}, {2,4}]

                sage: B = OrderedMultisetPartitionIntoSets([[1,2], [3]])
                sage: R._H_to_R(B)
                R[{1,2}, {3}] + R[{1,2,3}]

                sage: C = OrderedMultisetPartitionIntoSets([[1,2], [3,4], [5]])
                sage: R._H_to_R(C)
                R[{1,2}, {3,4}, {5}] + R[{1,2}, {3,4,5}]
                 + R[{1,2,3,4}, {5}] + R[{1,2,3,4,5}]
            """
            # base cases
            if A.length() <= 1:
                return self.monomial(A)

            # sum over strongly fatter elements
            terms = A.fatter(strong=True)
            return self.sum_of_terms([(B,1) for B in terms])

        def dual_basis(self):
            r"""
            Return the dual basis to the `\mathbf{R}` basis.

            The dual basis to the `\mathbf{R}` basis is the Fundamental
            basis of Dual of Free Hopf Algebra on Finite Sets.

            OUTPUT:

            - the Fundamental basis of Dual of Free Hopf Algebra on Finite Sets

            EXAMPLES::

                sage: R = OMPNonCommutativeSymmetricFunctions(QQ).Ribbon()
                sage: R.dual_basis()
                Dual of Free Hopf Algebra on Finite Sets over the Rational Field in the Fundamental basis
            """
            return self.realization_of().dual().Fundamental()

    Ribbon = R


###### Basis Methods for OMPQSym ############
class OMPBasis_OMPQSym(OMPBasis_abstract):
    """
    Add methods for `OMPQSym` beyond those appearing in ``OMPBasis_abstract``

    .. TODO::

        - Are there any maps to/from `WQSym`?

        - Are there any natural maps from `QSym`?:
          The transpose of the natural map from `OMPNSym` to `NSym` is not a map
          here (e.g., the image of `M_2` is `sum_{i<j} M_[{i,j}]`. And while
          restricting the alphabet makes the sum finite, it also breaks the
          product axiom `phi(a*b) = phi(a)*phi(b)`.

        - There is at least one character on `OMPQSym`, defined on the dual Powersum
          basis by sending `Pd_{A}` to `1` of `A` has length at most `1`, and to `0`
          otherwise. Consider using this to create a ``.to_quasisymmetric_function``
          method according to Theorem 4.1 of [AguBerSot2006].

    REFERENCES:

    - [AguBerSot2006] \M. Aguiar, N. Bergeron, F. Sottile. *Combinatorial
                      Hopf algebras and generalized Dehn–Sommerville relations*.
                      Compositio Math. 142 (2006) 1–30.
    """
    def is_commutative(self):
        """
        Return whether ``self`` is commutative.

        EXAMPLES::

            sage: OMPQuasiSymmetricFunctions(ZZ).M().is_commutative()
            True
        """
        return True

    class Element(OMPBasis_abstract.Element):
        def is_symmetric(self):
            r"""
            Determine if a `OMPQSym` function, expressed in the
            `\mathbf{M}` basis, is symmetric.

            A function `f` in the `\mathbf{M}` basis is deemed *symmetric*
            if for each ordered multiset partition `A` in its support, all
            derangements of `A` are also in the support of `f` with the
            same coefficient. That is,

            .. MATH::

                f = \sum_{\lambda} c_{\lambda}
                        \sum_{sort(A) = \lambda} \mathbf{M}_A,

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
                sage: elt += 2*M([[2],[2,3]])
                sage: elt.is_symmetric()
                False
                sage: elt += 2*M([[2,3],[2]])
                sage: elt.is_symmetric()
                True

            The symmetric functions form a subalgebra of `OMPQSym`, as
            these small examples illustrate::

                sage: M = OMPQuasiSymmetricFunctions(QQ).M()
                sage: x = M[[2,4]]; x.is_symmetric()
                True
                sage: y = M.sum_of_derangements([[1], [2]]); y.is_symmetric()
                True
                sage: x * y
                M[{2}, {1}, {2,4}] + M[{2}, {2,4}, {1}] + M[{2}, {1,2,4}]
                 + M[{2,4}, {1}, {2}] + M[{2,4}, {2}, {1}] + M[{1}, {2}, {2,4}]
                 + M[{1}, {2,4}, {2}] + M[{1,2,4}, {2}]
                sage: (x * y).is_symmetric()
                True
                sage: x = M.sum_of_derangements([[2], [3,4], [5]])
                sage: y = M.sum_of_derangements([[1], [2]])
                sage: x.is_symmetric(), y.is_symmetric()
                (True, True)
                sage: (x * y).is_symmetric()
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
                la = tuple(sorted(map(lambda k: tuple(sorted(k)), A)))
                if la not in d:
                    d[la] = [coeff, 1]
                else:
                    if d[la][0] != coeff:
                        return False
                    d[la][1] += 1
            # Make sure we've seen each ordered multiset partition in the derangement class
            return all(d[la][1] == Permutations_mset(la).cardinality() for la in d)

###### The Hopf algebra OMPQSym ############
class OMPQuasiSymmetricFunctions(UniqueRepresentation, Parent):
    r"""
    The Hopf algebra (`OMPQSym`) that is the graded dual to
    Free Hopf Algebra on Finite Sets (OMPNSym).

    INPUT:

    - ``R`` -- a base ring for the algebra

    OPTIONAL KEYWORDS:

    - ``alphabet`` -- a list of positive integers or an integer `n`.
      In the latter case, the alphabet is taken to be `\{1,\ldots,\}`.
      If alphabet is passed, then basis keys are restricted to those
      ordered multiset partitions into sets supported by ``alphabet``.

    - ``order_grading`` -- a boolean that defaults to ``True`` if an
      alphabet is passed, and to ``False`` otherwise. Determines the
      grading function to be used for `OMPQSym`. Let `A` be an ordered
      multiset partition into sets. If this keyword is ``True``, then
      use ``A.order()`` as a degree function; else use ``A.size()``.

    In addition to being the dual of the Hopf algebra `OMPNSym`, this
    Hopf algebra has a realization as quasi-symmetric functions in
    variables `x_{i,K}` satisfying the following relations:

    - `x_{i, K}` and `x_{j, L}` commute if `i \neq j`;
    - `x_{i, K} * x_{i, L}` equals zero when `K \cap L` is nonempty,
      and equals `x_{i, K \cup L}` otherwise.

    The quasi-symmetrizing action is on the first subscript, so, e.g.,

    .. MATH::

        (23) *\ast x_{1,\{2,5\}}x_{3, \{1,4,5\}} = x_{1,\{2,5\}}x_{2, \{1,4,5\}}

    The resulting quasi-symmetrized polynomial is denoted here by
    `M_{[\{2,5\}, \{1,4,5\}]}`. As the number of variables `n` tends to
    infinity, the product rule in this "monomial basis" becomes that
    of the basis `H^*`, where `H` is the Homogeneous basis of `OMPNSym`:

    .. MATH::

        M_A \cdot M_B = \sum_{C} M_C,

    where the sum is over "quasi-shuffles" of the ordered multiset
    partitions `A` and `B`. (That is, shuffle the blocks `a` of `A`
    with the blocks `b` of `B`, as well as merging `a\cup b` when
    `a\cap b` is empty.)

    A natural coproduct may be defined for these quasi-symmetric functions
    using an alphabet doubling procedure. The resulting formula is
    precisely the one holding for `H^*`:

    .. MATH::

        \Delta(M_C) = \sum_{A+B = C} M_A \otimes M_B,

    where the sum is taken over deconcatenations of the ordered multiset
    partition `C`.

    This Hopf algebra is implemented in the Monomial basis as a
    :class:`CombinatorialFreeModule` whose basis elements are indexed
    by *ordered multiset partitions into sets*, which are simply lists
    of subsets of positive integers.

    There are two natural gradings on `OMPQSym`: grading by the *size* or
    by the *order* of ordered multiset partitions. Only the first of these
    yields finite dimensional slices if a finite alphabet is not also
    provided. Hence, we allow an alphabet to be passed as an optional
    keyword argument. (The user may also elect to use the size grading,
    even if a finite alphabet is passed.)

    In the case that graded pieces are finite dimensional, then one may
    speak of the graded dual Hopf algebra `OMPNSym`.

    .. SEEALSO::

        :class:`OMPNonCommutativeSymmetricFunctions`
        :class:`OrderedMultisetPartitionsIntoSets`

    .. RUBRIC:: Standard use of Dual of Free Hopf Algebra on Finite Sets

    We begin by first creating the ring of `OMPQSym` and two bases that
    are analogues of the usual non-commutative symmetric functions::

        sage: OMPQSym = OMPQuasiSymmetricFunctions(QQ)
        sage: M = OMPQSym.M()
        sage: Pd = OMPQSym.Pd()
        sage: M
        Dual of Free Hopf Algebra on Finite Sets over the Rational Field
         in the Monomial basis

    The homogeneous component of degree `d` of any basis is indexed by
    ordered multiset partitions into sets of size `d`::

        sage: M.basis(3).keys()
        Ordered Multiset Partitions into Sets of integer 3

    So we create an element and convert it between these bases::

        sage: OMP = OrderedMultisetPartitionIntoSets
        sage: elt = M(OMP([[1,3], [2]])) + 9*M(OMP([[2,3], [1]])); elt
        9*M[{2,3}, {1}] + M[{1,3}, {2}]
        sage: Pd(elt)
        9*Pd[{2,3}, {1}] + Pd[{1,3}, {2}] - Pd[{1,2,3}]

    Given an element of `OMPNSym` expressed in two different bases, the two
    expressions are treated as equal::

        sage: Pd(elt) == elt
        True
        sage: dict(Pd(elt)) == dict(elt)
        False

    There is also a shorthand for creating elements. We note that we must
    use ``M[[]]`` to create the empty ordered multiset partition due to
    python's syntax::

        sage: eltm = M[[2], [1,3]] - 2*M[[1,2],[3]] + M[[]]; eltm
        M[] + M[{2}, {1,3}] - 2*M[{1,2}, {3}]

    Finally, expressions are returned in only one basis, even if given
    in mixed bases::

        sage: eltp = Pd[[1,2,3]] + eltm; eltp
        Pd[] + Pd[{2}, {1,3}] - 2*Pd[{1,2}, {3}] + 3*Pd[{1,2,3}]

        sage: eltm * Pd[[1,2,3]]
        M[{2}, {1,3}, {1,2,3}] + M[{2}, {1,2,3}, {1,3}]
         - 2*M[{1,2}, {3}, {1,2,3}] - 2*M[{1,2}, {1,2,3}, {3}]
         + M[{1,2,3}, {2}, {1,3}] - 2*M[{1,2,3}, {1,2}, {3}] + M[{1,2,3}]

    .. RUBRIC:: Restricting the alphabet and using the order grading

    When the keyword ``alphabet`` is used to select an alphabet `A`, then
    the homogeneous component of degree `d` of the basis is indexed by
    ordered multiset partitions over `A` of order `d`::

        sage: OMPQSymA = OMPQuasiSymmetricFunctions(QQ, alphabet=[3,4,5])
        sage: MA = OMPQSymA.M()
        sage: MA.basis(3).keys()
        Ordered Multiset Partitions into Sets of order 3 over alphabet {3, 4, 5}

    The keyword argument ``order_grading`` may be used to override the default
    grading::

        sage: OMPQSym_order = OMPQuasiSymmetricFunctions(QQ, order_grading=True)
        sage: M_order = OMPQSym_order.M()
        sage: M_order.basis(2).keys()
        #??? wrong!!! Ordered Multiset Partitions into Sets of integer 2

        sage: OMPQSymA_size = OMPQuasiSymmetricFunctions(QQ, alphabet=[3,4,5], order_grading=False)
        sage: M_size = OMPQSymA_size.M()
        sage: M_size.basis(2).keys()
        #??? wrong!!! Ordered Multiset Partitions into Sets of order 2 over alphabet {3, 4, 5}

    *Warning:* Realizations of `OMPQSym` such as ``OMPQSym_order`` above will
    have infinite dimensional spaces for their graded slices `OMPQSym_d`.

    .. RUBRIC:: Moving between realizations of `OMPQSym`

    If parents' alphabets are compatible, coerce elements between different
    realizations of `OMPQSym`::

        sage: MA.an_element()
        2*M[{3}, {3}, {4,5}] + M[{3}, {3,4,5}]
        sage: M.an_element()
        2*M[{1}, {1,2}] + M[{4}]
        sage: M(MA.an_element()) in M
        True
        sage: MA(M.an_element()) in MA
        Traceback (most recent call last):
        ...
        TypeError: do not know how to make x (= 2*M[{1}, {1,2}] + M[{4}])
         an element of self (=Dual of Free Hopf Algebra on Finite Sets on alphabet
         {3, 4, 5} over the Rational Field with order grading in the Monomial basis)

    TESTS::

        sage: M = OMPQuasiSymmetricFunctions(QQ).M()
        sage: a = M[OrderedMultisetPartitionIntoSets([[1]])]
        sage: b = M[OrderedMultisetPartitionsIntoSets(1)([[1]])]
        sage: c = M[[1]]
        sage: a == b == c
        True
        sage: M = OMPQuasiSymmetricFunctions(QQ, alphabet=[2]).M()
        sage: a = M[OrderedMultisetPartitionIntoSets([[2]])]
        sage: b = M[OrderedMultisetPartitionsIntoSets(2)([[2]])]
        sage: c = M[[2]]
        sage: d = M[OrderedMultisetPartitionsIntoSets([1,2], 1)([[2]])]
        sage: a == b == c == d
        True
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
                raise ValueError("keyword alphabet (=%s) must be a nonempty set of positive integers"%(_A))
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
        category = HopfAlgebras(R).Graded().Connected().Commutative()
        Parent.__init__(self, base=R, category=category.WithRealizations())
        self._A = alphabet
        self._order_grading = order_grading


    def _repr_(self):
        r"""
        EXAMPLES::

            sage: OMPQuasiSymmetricFunctions(ZZ)
            Dual of Free Hopf Algebra on Finite Sets over the Integer Ring
        """
        return "Dual of " + repr(self.dual())

    def a_realization(self):
        r"""
        Return the realization of the `\mathbf{M}` basis of ``self``.

        EXAMPLES::

            sage: OMPQuasiSymmetricFunctions(QQ).a_realization()
            Dual of Free Hopf Algebra on Finite Sets over the Rational Field in the Monomial basis
        """
        return self.M()

    _shorthands = tuple(['M', 'Pd', 'F'])

    def dual(self):
        r"""
        Return the Free Hopf Algebra on Finite Sets.

        EXAMPLES::

            sage: OMPD = OMPQuasiSymmetricFunctions(QQ)
            sage: OMPD.dual()
            Free Hopf Algebra on Finite Sets over the Rational Field
        """
        return OMPNonCommutativeSymmetricFunctions(self.base_ring(), self._A, self._order_grading)

    class M(OMPBasis_OMPQSym):
        r"""
        The Monomial basis of Dual of Hopf Algebra on Ordered Multiset Partitions.

        The family `(M_\mu)`, as `\mu` ranges over all ordered multiset
        partitions into sets, is the basis that is graded dual to the
        Homogeneous basis of Hopf Algebra on Ordered Multiset Partitions.
        By analogy with :class:`QuasiSymmetricFunctions`, it is called the
        *Monomial basis* here.

        Analogous to the Monomial basis of :class:`QuasiSymmetricFunctions`,
        product and coproduct rules for `(M_\mu)` are given by "quasi-shuffle"
        and "deconcatenate", respectively.

        EXAMPLES::

            sage: M = OMPQuasiSymmetricFunctions(QQ).M()
            sage: M[[2], [1,3]] - 2*M[[1,2,4]]
            M[{2}, {1,3}] - 2*M[{1,2,4}]

            sage: M[[2], [1,3]] * M[[2,4]]
            M[{2}, {2,4}, {1,3}] + M[{2}, {1,3}, {2,4}]
             + M[{2}, {1,2,3,4}] + M[{2,4}, {2}, {1,3}]
            sage: M[[5], [1,3]] * M[[2,4]]
            M[{5}, {2,4}, {1,3}] + M[{5}, {1,3}, {2,4}] + M[{5}, {1,2,3,4}]
             + M[{2,4}, {5}, {1,3}] + M[{2,4,5}, {1,3}]

            sage: M[[2], [1,3], [1,2,4]].coproduct()
            M[] # M[{2}, {1,3}, {1,2,4}] + M[{2}] # M[{1,3}, {1,2,4}]
             + M[{2}, {1,3}] # M[{1,2,4}] + M[{2}, {1,3}, {1,2,4}] # M[]

        TESTS::

            sage: M = OMPQuasiSymmetricFunctions(QQ).M()
            sage: TestSuite(M).run()  # long time

            sage: M = OMPQuasiSymmetricFunctions(QQ).M()
            sage: H = OMPNonCommutativeSymmetricFunctions(QQ).H()
            sage: A = OrderedMultisetPartitionIntoSets([[1,2],[3]])
            sage: B = OrderedMultisetPartitionIntoSets([[3]])
            sage: C = OrderedMultisetPartitionIntoSets([[1,2,3],[3]])
            sage: D = OrderedMultisetPartitionIntoSets([[1,2],[3],[3]])

            sage: M[A] * M[B]
            2*M[{1,2}, {3}, {3}] + M[{3}, {1,2}, {3}] + M[{1,2,3}, {3}]
            sage: (M[A] * M[B]).duality_pairing(H[C])
            1
            sage: (M[A] * M[B]).duality_pairing(H[D])
            2

            sage: M[C].coproduct()
            M[] # M[{1,2,3}, {3}] + M[{1,2,3}] # M[{3}] + M[{1,2,3}, {3}] # M[]
            sage: [M[pi].duality_pairing(H[A]) * M[mu].duality_pairing(H[B]) \
                   for ((pi,mu), c) in M[C].coproduct()]
            [0, 0, 0]
            sage: M[D].coproduct()
            M[] # M[{1,2}, {3}, {3}] + M[{1,2}] # M[{3}, {3}]
             + M[{1,2}, {3}] # M[{3}] + M[{1,2}, {3}, {3}] # M[]
            sage: [M[pi].duality_pairing(H[A]) * M[mu].duality_pairing(H[B]) \
                   for ((pi,mu), c) in M[D].coproduct()]
            [0, 0, 0, 1]

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
            True
        """
        _prefix = "M"
        _basis_name = "Monomial"

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

        def product_on_basis(self, A, B):
            r"""
            The product on `\mathbf{M}` basis elements.

            This commutative product on the `\mathbf{M}` basis is the dual to
            the coproduct on the `\mathbf{H}` basis. On the basis `\mathbf{M}`
            it is defined as

            .. MATH::

                M_A \cdot M_B = \sum_{C} M_C,

            where the sum is over all overlapping shuffles of the blocks of
            `A` with the blocks of `B`. This is (ultimately) implemented via
            :class:`~sage.combinat.shuffle.ShuffleProduct_overlapping` with
            the ``+`` operation for the overlapping defined as:
            the union of sets, if the summands are disjoint; and ``None``
            otherwise.

            INPUT:

            - ``A``, ``B`` -- ordered multiset partitions into sets

            OUTPUT:

            - an element of the `\mathbf{M}` basis

            EXAMPLES::

                sage: M = OMPNonCommutativeSymmetricFunctions(QQ).dual().M()
                sage: A = OrderedMultisetPartitionsIntoSets()._from_list([2,1,3,0,1,2]); A
                [{1,2,3}, {1,2}]
                sage: B = OrderedMultisetPartitionIntoSets([[3,4]]); B
                [{3,4}]
                sage: M.product_on_basis(A, B)
                M[{3,4}, {1,2,3}, {1,2}] + M[{1,2,3}, {3,4}, {1,2}]
                 + M[{1,2,3}, {1,2}, {3,4}] + M[{1,2,3}, {1,2,3,4}]
                sage: C = OrderedMultisetPartitionIntoSets([[4,5]]); C
                [{4,5}]
                sage: M.product_on_basis(A, C)
                M[{4,5}, {1,2,3}, {1,2}] + M[{1,2,3}, {4,5}, {1,2}] + M[{1,2,3}, {1,2}, {4,5}] + M[{1,2,3}, {1,2,4,5}] + M[{1,2,3,4,5}, {1,2}]
                sage: M.product_on_basis(A, OrderedMultisetPartitionIntoSets([]))
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

                sage: M = OMPNonCommutativeSymmetricFunctions(QQ).dual().M()
                sage: M[[2], [2,3]].coproduct()
                M[] # M[{2}, {2,3}] + M[{2}] # M[{2,3}] + M[{2}, {2,3}] # M[]
                sage: M.coproduct_on_basis(OrderedMultisetPartitionIntoSets([]))
                M[] # M[]
            """
            return self.tensor_square().sum_of_monomials(A.deconcatenate(2))

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
                M[{1}, {1,2}] + M[{1,2}, {1}]
            """
            OMP = self._indices
            A = OMP(A)
            return self.sum_of_terms([(OMP(pi),1) \
                        for pi in Permutations_mset(list(A))], distinct=True)

    Monomial = M

    class Pd(OMPBasis_OMPQSym):
        r"""
        The dual Powersum basis of Dual of Free Hopf Algebra on Finite Sets.

        The family `(P_\mu)`, as `\mu` ranges over all ordered multiset
        partitions into sets, is called the *dual Powersum basis* here. It
        is defined in relation to the Monomial basis of `OMPQSym` as the
        transpose of the relationship between the bases `(H_\mu)` and
        `(P_\nu)` of `OMPNSym`.

        .. TODO:: give a more explicit description.

        EXAMPLES::

            sage: Pd = OMPQuasiSymmetricFunctions(QQ).PowersumDual()
            sage: p = Pd[[2,4], [3]] - 2*Pd[[1,2,3]]; p
            -2*Pd[{1,2,3}] + Pd[{2,4}, {3}]

            sage: M = OMPQuasiSymmetricFunctions(QQ).Monomial()
            sage: M(p)
            M[{2,4}, {3}] - 2*M[{1,2,3}] + M[{2,3,4}]

            sage: Pd[[2,3]].coproduct()
            Pd[] # Pd[{2,3}] + Pd[{2,3}] # Pd[]
            sage: Pd[[2,3],[1,2]].coproduct()
            Pd[] # Pd[{2,3}, {1,2}] + Pd[{2,3}] # Pd[{1,2}] + Pd[{2,3}, {1,2}] # Pd[]

            sage: Pd[[2,3]] * Pd[[1,2]]
            Pd[{2,3}, {1,2}] + Pd[{1,2}, {2,3}]

        TESTS::

            sage: TestSuite(Pd).run()  # long time
            sage: all(Pd(M(p)) == p for p in Pd.some_elements())
            True
        """
        _prefix = "Pd"
        _basis_name = "dual Powersum"

        def __init__(self, alg):
            """
            Initialize ``self``.
            """
            OMPBasis_OMPQSym.__init__(self, alg)

            # Register coercions
            M = self.realization_of().M()
            phi = self.module_morphism(self._Pd_to_M, codomain=M, triangular="lower")
            phi.register_as_coercion()
            #(~phi).register_as_coercion()
            phi_inv = M.module_morphism(self._M_to_Pd, codomain=self, triangular="lower")
            phi_inv.register_as_coercion()

        def _Pd_to_M(self, A):
            """
            Return `M_A` in terms of the Monomial basis.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - An element of the dual Powersum basis

            EXAMPLES::

                sage: Pd = OMPQuasiSymmetricFunctions(QQ).Pd()
                sage: A = OrderedMultisetPartitionIntoSets([[1,2,3]])
                sage: Pd._Pd_to_M(A)
                M[{1,2,3}]
                sage: A = OrderedMultisetPartitionIntoSets([[2,3,4],[1]])
                sage: Pd._Pd_to_M(A)
                M[{2,3,4}, {1}]

                sage: B = OrderedMultisetPartitionIntoSets([[1,3,4],[2]])
                sage: Pd._Pd_to_M(B)
                M[{1,3,4}, {2}] + M[{1,2,3,4}]
                sage: B = OrderedMultisetPartitionIntoSets([[1,3,6],[2,4],[5]])
                sage: Pd._Pd_to_M(B)
                M[{1,3,6}, {2,4}, {5}] + M[{1,3,6}, {2,4,5}]
                 + M[{1,2,3,4,6}, {5}] + M[{1,2,3,4,5,6}]
            """
            M = self.realization_of().M()
            # base cases
            if A.length() <= 1:
                return M.monomial(A)

            # Sum over certain elements in ``A.fatter()``.
            # A merging of adjacent blocks is allowed if their min elements are in order.
            # This condition always reaches back to `A`, so one cannot just
            # iterate through elements of `A`.fatter().
            def glueings(tupl_of_sets):
                terms = [tupl_of_sets]
                for A in terms:
                    for i in range(len(A)-1):
                        A1,A2 = A[i:i+2]
                        if min(A1) < min(A2) and len(A1)+len(A2) == len(A1.union(A2)):
                            B = A[:i] + (A1.union(A2),) + A[i+2:]
                            terms.append(tuple([b for b in B if b]))
                # kill undesired repeats within added_terms with ``distinct=True``
                return set(terms)
            descents = [i+1 for i in range(len(A)-1) if min(A[i]) >= min(A[i+1])]
            descents = [0] + descents + [len(A)]
            factored_A = [tuple(A[descents[i]:descents[i+1]]) for i in range(len(descents)-1)]
            out = M.zero()
            for AA in cartesian_product([glueings(ll) for ll in factored_A]):
                out += M.monomial(sum(map(self._indices, AA)))
            return out


        def _M_to_Pd(self, A):
            """
            Return `M_A` in terms of the dual Powersum basis.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - An element of the dual Powersum basis

            EXAMPLES::

                sage: Pd = OMPQuasiSymmetricFunctions(QQ).Pd()
                sage: A = OrderedMultisetPartitionIntoSets([[1,2,3]])
                sage: Pd._M_to_Pd(A)
                Pd[{1,2,3}]
                sage: A = OrderedMultisetPartitionIntoSets([[2,3,4],[1]])
                sage: Pd._M_to_Pd(A)
                Pd[{2,3,4}, {1}]

                sage: B = OrderedMultisetPartitionIntoSets([[1,3,4],[2]])
                sage: Pd._M_to_Pd(B)
                -Pd[{1,2,3,4}] + Pd[{1,3,4}, {2}]
                sage: B = OrderedMultisetPartitionIntoSets([[1,3,6],[2,4],[5]])
                sage: Pd._M_to_Pd(B)
                -Pd[{1,3,6}, {2,4,5}] - Pd[{1,2,3,4,6}, {5}]
                 + Pd[{1,2,3,4,5,6}] + Pd[{1,3,6}, {2,4}, {5}]
            """
            # base cases
            if A.length() <= 1:
                return self.monomial(A)

            # signed sum over certain elements in ``A.fatter()``
            # a merging of adjacent blocks is allowed if their min elements are in order
            terms = [(A,1)]
            for (A,c) in terms:
                for i in range(A.length()-1):
                    A1,A2 = A[i:i+2]
                    if min(A1) < min(A2) and len(A1)+len(A2) == len(A1.union(A2)):
                        B = A[:i] + [A1.union(A2)] + A[i+2:]
                        terms.append((self._indices([b for b in B if b]),-c))
            # kill undesired repeats within added_terms with ``distinct=True``
            return self.sum_of_terms(terms, distinct=True)

        def dual_basis(self):
            r"""
            Return the dual basis to the `\mathbf{Pd}` basis.

            The dual basis to the `\mathbf{Pd}` basis is the Powersum
            basis of Free Hopf Algebra on Finite Sets.

            OUTPUT:

            - the Powersum basis of Free Hopf Algebra on Finite Sets

            EXAMPLES::

                sage: Pd = OMPQuasiSymmetricFunctions(QQ).PowersumDual()
                sage: Pd.dual_basis()
                Free Hopf Algebra on Finite Sets over the Rational Field in the Powersum basis
            """
            return self.realization_of().dual().Powersum()

        def product_on_basis(self, A, B):
            r"""
            Return the (associative, commutative) `*` product of the
            basis elements of ``self`` indexed by the ordered multiset
            partitions `A` and `B`.

            INPUT:

            - ``A, B`` -- two ordered multiset partitions into sets

            OUTPUT:

            - The shuffle of the blocks of `A` with the blocks of `B`.

            EXAMPLES::

                sage: Pd = OMPQuasiSymmetricFunctions(QQ, alphabet=[2,3,4,5]).PowersumDual()
                sage: Pd[{2,3}] * Pd[{3,5}, {4}]
                Pd[{2,3}, {3,5}, {4}] + Pd[{3,5}, {2,3}, {4}] + Pd[{3,5}, {4}, {2,3}]

            TESTS::

                sage: OMP = OrderedMultisetPartitionsIntoSets(alphabet=[2,3,4])
                sage: M = Pd.realization_of().Monomial()
                sage: one = OMP([])
                sage: all(Pd.product_on_basis(one, z) == Pd(z) == Pd.basis()[z] for z in [one] + list(OMP.subset(3)))
                True
                sage: all(Pd.product_on_basis(z, one) == Pd(z) == Pd.basis()[z] for z in OMP.subset(3))
                True
                sage: all(Pd[A] * Pd[B] == Pd( M(Pd[A])*M(Pd[B]) ) for A in OMP.subset(3) for B in OMP.subset(2))  # indirect doctest
                True
            """
            terms = A.shuffle_product(B, overlap=False)
            return self.sum_of_terms([(s, 1) for s in terms])

        def coproduct_on_basis(self, A):
            r"""
            Return the coproduct in the basis `(Pd_\mu)`.

            As for the Monomial basis, the rule for coproducts in the
            dual Powersum basis is "deconcatenate".

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            EXAMPLES::

                sage: Pd = OMPQuasiSymmetricFunctions(QQ).PowersumDual()
                sage: Pd[{2,3,4}].coproduct()
                Pd[] # Pd[{2,3,4}] + Pd[{2,3,4}] # Pd[]
                sage: Pd[{1}, {2,3,4}, {1,2}].coproduct()
                Pd[] # Pd[{1}, {2,3,4}, {1,2}] + Pd[{1}] # Pd[{2,3,4}, {1,2}]
                 + Pd[{1}, {2,3,4}] # Pd[{1,2}] + Pd[{1}, {2,3,4}, {1,2}] # Pd[]

            TESTS::

                sage: OMP = OrderedMultisetPartitionsIntoSets
                sage: M = Pd.realization_of().Monomial()
                sage: all(Pd.coproduct_on_basis(A) == Pd.tensor_square()( M(Pd[A]).coproduct() ) for A in OMP(4))
                True
            """
            return self.tensor_square().sum_of_monomials(A.deconcatenate(2))

        def antipode_on_basis(self, A):
            r"""
            Return the antipode applied to a dual Powersum basis element.

            If `A` is an ordered multiset partition with a single block then
            the antipode is given by `S(Pd_A) = -Pd_A`.

            If `A` has more than one block, then we extend using the fact that
            the antipode is an algebra antimorphism.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            EXAMPLES::

                sage: Pd = OMPQuasiSymmetricFunctions(QQ).PowersumDual()
                sage: Pd[{2,3}].antipode()
                -Pd[{2,3}]
                sage: Pd[{2,3}, {2}].antipode()
                Pd[{2}, {2,3}]

            TESTS::

                sage: Pd[[]].antipode()
                Pd[]
                sage: PP = Pd[[1,2,3],[5],[1,4]].coproduct()
                sage: PP.apply_multilinear_morphism(lambda x,y: x.antipode()*y)
                0
            """
            return (-1)**A.length() * self.monomial(A.reversal())

    PowersumDual = Pd

    class F(OMPBasis_OMPQSym):
        r"""
        The Fundamental basis of Dual of Free Hopf Algebra on Finite Sets.

        The family `(F_\mu)`, as `\mu` ranges over all ordered multiset
        partitions into sets, is called the *Fundamental basis* here. It
        is defined in relation to the Monomial basis of `OMPQSym` as "sum
        above `\mu` in the strong refinement order."

        .. TODO::

            - decide if this should be called the "fundamental".
              Pros: defined via refinement; and product and coproduct are
              reminiscent of those for Gessel's basis for quasisymmetric
              functions.
              Con: the antipode does not act as a signed permutation matrix
              in this basis, as it does for Gessel's basis.

        EXAMPLES::

            sage: F = OMPQuasiSymmetricFunctions(QQ).Fundamental()
            sage: F[[2,3]].coproduct()
            F[] # F[{2,3}] + F[{2}] # F[{3}] + F[{2,3}] # F[]
            sage: F[[2,3], [1,2]].coproduct()
            F[] # F[{2,3}, {1,2}] + F[{2}] # F[{3}, {1,2}]
             + F[{2,3}] # F[{1,2}] + F[{2,3}, {1}] # F[{2}]
             + F[{2,3}, {1,2}] # F[]

            sage: F[[2,3]] * F[[1,2]]
            F[{2}, {1,3}, {2}] + F[{2}, {1,2,3}] + 2*F[{1,2}, {2,3}]
            + F[{2,3}, {1,2}] + F[{1,2,3}, {2}]

            sage: M = OMPQuasiSymmetricFunctions(QQ).Monomial()
            sage: p = F[[1,2],[3],[4,5]]
            sage: M(p)
            M[{1}, {2}, {3}, {4}, {5}] + M[{1}, {2}, {3}, {4,5}]
            + M[{1,2}, {3}, {4}, {5}] + M[{1,2}, {3}, {4,5}]
            sage: q = F[[1,2],[3],[1,4,5]]
            sage: M(q)
            M[{1}, {2}, {3}, {1}, {4}, {5}] + M[{1}, {2}, {3}, {1}, {4,5}]
            + M[{1}, {2}, {3}, {1,4}, {5}] + M[{1}, {2}, {3}, {1,4,5}]
            + M[{1,2}, {3}, {1}, {4}, {5}] + M[{1,2}, {3}, {1}, {4,5}]
            + M[{1,2}, {3}, {1,4}, {5}] + M[{1,2}, {3}, {1,4,5}]

        TESTS::

            sage: TestSuite(F).run()  # long time
            sage: all(F(M(p)) == p for p in F.some_elements())
            True
        """
        _prefix = "F"
        _basis_name = "Fundamental"

        def __init__(self, alg):
            """
            Initialize ``self``.
            """
            OMPBasis_OMPQSym.__init__(self, alg)

            # Register coercions
            M = self.realization_of().M()
            phi = self.module_morphism(self._F_to_M, codomain=M, triangular="lower")
            phi.register_as_coercion()
            #(~phi).register_as_coercion()
            phi_inv = M.module_morphism(self._M_to_F, codomain=self, triangular="lower")
            phi_inv.register_as_coercion()

        def _F_to_M(self, A):
            """
            Return `F_A` in terms of the Monomial basis.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - An element of the Monomial basis

            EXAMPLES::

                sage: F = OMPQuasiSymmetricFunctions(QQ).F()
                sage: A = OrderedMultisetPartitionIntoSets([])
                sage: F._F_to_M(A)
                M[]
                sage: A = OrderedMultisetPartitionIntoSets([[1], [3], [3]])
                sage: F._F_to_M(A)
                M[{1}, {3}, {3}]

                sage: B = OrderedMultisetPartitionIntoSets([[1,3], [2]])
                sage: F._F_to_M(B)
                M[{1}, {3}, {2}] + M[{1,3}, {2}]
                sage: B = OrderedMultisetPartitionIntoSets([[1,3], [1,4]])
                sage: F._F_to_M(B)
                M[{1}, {3}, {1}, {4}] + M[{1}, {3}, {1,4}]
                 + M[{1,3}, {1}, {4}] + M[{1,3}, {1,4}]

                sage: C = OrderedMultisetPartitionIntoSets([[1,3,4]])
                sage: F._F_to_M(C)
                M[{1}, {3}, {4}] + M[{1}, {3,4}] + M[{1,3}, {4}] + M[{1,3,4}]
            """
            M = self.realization_of().M()
            # base cases
            if A.length() == A.order():
                return M.monomial(A)

            # sum over strongly finer elements
            terms = A.finer(strong=True)
            return M.sum_of_terms([(B,1) for B in terms])

        def _M_to_F(self, A):
            """
            Return `M_A` in terms of the Fundamental basis.

            INPUT:

            - ``A`` -- an ordered multiset partition into sets

            OUTPUT:

            - An element of the Fundamental basis

            EXAMPLES::

                sage: F = OMPQuasiSymmetricFunctions(QQ).F()
                sage: A = OrderedMultisetPartitionIntoSets([])
                sage: F._M_to_F(A)
                F[]
                sage: A = OrderedMultisetPartitionIntoSets([[1], [3], [3]])
                sage: F._M_to_F(A)
                F[{1}, {3}, {3}]

                sage: B = OrderedMultisetPartitionIntoSets([[1,3], [2]])
                sage: F._M_to_F(B)
                -F[{1}, {3}, {2}] + F[{1,3}, {2}]
                sage: B = OrderedMultisetPartitionIntoSets([[1,3], [1,4]])
                sage: F._M_to_F(B)
                F[{1}, {3}, {1}, {4}] - F[{1}, {3}, {1,4}]
                 - F[{1,3}, {1}, {4}] + F[{1,3}, {1,4}]

                sage: C = OrderedMultisetPartitionIntoSets([[1,3,4]])
                sage: F._M_to_F(C)
                F[{1}, {3}, {4}] - F[{1}, {3,4}] - F[{1,3}, {4}] + F[{1,3,4}]
            """
            PS = Poset({B:B.finer_succ(strong=True) for B in A.finer(strong=True)})
            return self.sum_of_terms([(B, PS.moebius_function(A, B)) for B in PS])

        def dual_basis(self):
            r"""
            Return the dual basis to the `\mathbf{F}` basis.

            The dual basis to the `\mathbf{F}` basis is the Ribbon
            basis of Free Hopf Algebra on Finite Sets.

            OUTPUT:

            - the Ribbon basis of Free Hopf Algebra on Finite Sets

            EXAMPLES::

                sage: F = OMPQuasiSymmetricFunctions(QQ).Fundamental()
                sage: F.dual_basis()
                Free Hopf Algebra on Finite Sets over the Rational Field in the Ribbon basis
            """
            return self.realization_of().dual().Ribbon()

    Fundamental = F
