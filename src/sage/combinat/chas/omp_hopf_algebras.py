"""
Hopf Algebra on Ordered Multiset Partitions

AUTHORS:

- Aaron Lauve (2018): initial implementation, following ``wqsym.py``.
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
#from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.cartesian_product import cartesian_product

from sage.matrix.matrix_space import MatrixSpace
from sage.sets.set import Set
from sage.rings.all import Integer, ZZ, QQ

from sage.functions.other import factorial

from sage.combinat.free_module import CombinatorialFreeModule
#from .bases import OMPBases, MultiplicativeOMPBases, OMPBasis_abstract
from sage.combinat.multiset_partition_ordered import OrderedMultisetPartition, OrderedMultisetPartitions, OrderedMultisetPartitions_n, OrderedMultisetPartitions_A
from sage.sets.set import Set, Set_object
from sage.combinat.posets.posets import Poset
from sage.combinat.sf.sf import SymmetricFunctions

class OMPBasis_abstract(CombinatorialFreeModule, BindableClass):
    """
    Abstract base class for bases of `OMPSym` and `OMPQSym`.

    This must define two attributes:

    - ``_prefix`` -- the basis prefix
    - ``_basis_name`` -- the name of the basis (must match one
      of the names that the basis can be constructed from `OMPSym`)
    """
    def __init__(self, alg, graded=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: HopfAlgebraOnOrderedMultisetPartitions(ZZ).H()
            Hopf algebra of ordered multiset partitions over Integer Ring in the Homogeneous basis
            sage: HopfAlgebraOnOrderedMultisetPartitions(ZZ, alphabet=[2,3,4]).H()
            Hopf algebra of ordered multiset partitions on alphabet {2, 3, 4} over Integer Ring with order grading in the Homogeneous basis

        TESTS::

            sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
            sage: TestSuite(H).run()  # long time
            sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ, alphabet=[2,3,4]).H()
            sage: TestSuite(H).run()  # long time
        """
        self._A = alg._A
        self._order_grading = alg._order_grading
        if self._A:
            _OMPs = OrderedMultisetPartitions(alphabet=self._A)
            _OMPs._repr_ = lambda : "Ordered Multiset Partitions over alphabet %s"%Set(self._A)
        else:
            _OMPs = OrderedMultisetPartitions()
        #def key_func_omp(A):
        #    # not sure we need to do anything special here!
        #    return A
        CombinatorialFreeModule.__init__(self, alg.base_ring(),
                                         _OMPs,
                                         category=OMPBases(alg, graded),
                                         #sorting_key=key_func_omp,
                                         bracket="", prefix=self._prefix)

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - elements of `OMPSym` over a base with
          a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: H = HopfAlgebraOnOrderedMultisetPartitions(GF(7)).H(); H
            Hopf Algebra on Ordered Multiset Partitions over Finite Field of size 7 in the Homogeneous basis

        Elements of the Homogeneous basis of OMPSym canonically coerce in::

            sage: x, y = H([[1]]), H([[2,1]])
            sage: H.coerce(x+y) == x+y
            True

        Elements of the Homogeneous basis of OMPSym over `\ZZ` coerces in,
        since `\ZZ` coerces to `\GF{7}`::

            sage: HZ = HopfAlgebraOnOrderedMultisetPartitions(ZZ).H(); HZ
            sage: xZ, yZ = HZ([[1]]), HZ([[2,1]])
            sage: z = H.coerce(xZ+yZ); z
            H[{1}] + H[{1,2}]
            sage: z.parent() is H
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so `OMPSym` over `\GF{7}`
        does not coerce to the same algebra over `\ZZ`::

            sage: HZ.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Hopf Algebra on Ordered Multiset
             Partitions over Finite Field of size 7 in the Homogeneous basis to
             Hopf Algebra on Ordered Multiset Partitions over Integer Ring in
             the Monomial basis

        TESTS::

            sage: HZ = HopfAlgebraOnOrderedMultisetPartitions(ZZ).H()
            sage: HQ = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
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
        """
        # Hopf algebra on ordered multiset partitions
        # over any base that coerces in:
        if isinstance(R, OMPBasis_abstract):
            if R.realization_of() == self.realization_of():
                return True
            if not self.base_ring().has_coerce_map_from(R.base_ring()):
                return False
            if self._basis_name == R._basis_name: # The same basis
                def coerce_base_ring(self, x):
                    return self._from_dict(x.monomial_coefficients())
                return coerce_base_ring
            # Otherwise lift that basis up and then coerce over
            target = getattr(self.realization_of(), R._basis_name)()
            return self._coerce_map_via([target], R)
        return super(OMPBasis_abstract, self)._coerce_map_from_(R)

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
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
        Return some elements of the Hopf Algebra on Ordered Multiset Partitions.

        EXAMPLES::

            sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
            sage: H.some_elements()
            [H[], H[{1}], 1/2*H[] + H[{4}] + 2*H[{1}, {1,2}],
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
        oa_ao = self.sum_of_terms([(x+y, c) for (y,c) in _a]) \
                    - self.sum_of_terms([(y+x, -c) for (y,c) in _a])
        return [u, o, s + a, oa_ao]

class HopfAlgebraOnOrderedMultisetPartitions(UniqueRepresentation, Parent):
    r"""
    Hopf Algebra on Ordered Multiset Partitions (OMP).

    The Hopf algebra OMP is a free algebra built on finite subsets
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

        :class:`OrderedMultisetPartitions`.

    REFERENCES:

    - [LM2018]_

    EXAMPLES:

    Standard use of Hopf Algebra on Ordered Multiset Partitions:

    We begin by first creating the ring of `OMPSym` and the bases that are
    analogues of the usual non-commutative symmetric functions::

        sage: OMPSym = HopfAlgebraOnOrderedMultisetPartitions(QQ)
        sage: H = OMPSym.H()
        sage: P = OMPSym.P() # ???
        sage: H
        Hopf Algebra on Ordered Multiset Partitions over Rational Field
         in the Homogeneous basis

    The basis is indexed by ordered multiset partitions, so we create an
    element and convert it between these bases::

        sage: elt = H(OrderedMultisetPartition([[2], [1,3]])) - 2*H(OrderedMultisetPartition([[1,2,3]])); elt
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

        sage: OMPSymA = HopfAlgebraOnOrderedMultisetPartitions(QQ, alphabet=[3,4,5])
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

        sage: OMPSymB = HopfAlgebraOnOrderedMultisetPartitions(QQ, alphabet=[3,4,5], order_grading=False)
        sage: HB = OMPSymB.HB()
        sage: HB.basis(12).keys().list()
        ???

    TESTS::

        sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
        sage: a = H[OrderedMultisetPartition([[1]])]
        sage: b = H[OrderedMultisetPartitions(1)([[1]])]
        sage: c = H[[1]]
        sage: a == b == c
        True
        sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ, alphabet=[2]).H()
        sage: a = H[OrderedMultisetPartition([[2]])]
        sage: b = H[OrderedMultisetPartitions(2)([[2]])]
        sage: c = H[[2]]
        sage: d = H[OrderedMultisetPartitions([1,2], 1)([[2]])]
        sage: a == b == c == d
        True

    .. TODO::

        - Add/correct the powersum basis.
        - Add/correct maps to NCSym, NSym, Sym.
        - Implement dual basis of `OMPSym` (which I call `OMPQSym` below).
        - Implement the alphabet keyword and order grading.
           (I believe the issue stems from the fact that ``OrderedMultisetPartitions_A``
            takes two arguments, while ``CombinatorialFreeModule`` expects its basis
            keys to take one argument, as ``OrderedMultisetPartitions_n`` does.)
    """
    @staticmethod
    def __classcall_private__(cls, R, alphabet=None, order_grading=True):
        """
        Normalize the input to ensure a unique representation.

        TESTS::

            sage: OMP1 = HopfAlgebraOnOrderedMultisetPartitions(QQ, alphabet=[2,3,5])
            sage: OMP2 = HopfAlgebraOnOrderedMultisetPartitions(QQ, alphabet=Set([2,3,5]))
            sage: OMP1 is OMP2
            True
        """
        if alphabet is not None:
            if alphabet in ZZ:
                alphabet = Set(range(1,alphabet+1))
            else:
                alphabet = Set(alphabet)
            if alphabet == Set() or any(a not in ZZ for a in alphabet):
                raise ValueError("keyword alphabet was converted to %s, which must be a nonempty set of positive integers"%(Set(alphabet)))
        else:
            alphabet = None
        return super(HopfAlgebraOnOrderedMultisetPartitions, cls).__classcall__(cls, R, alphabet, order_grading)

    def __init__(self, R, alphabet=None, order_grading=True):
        """
        Initialize ``self``.

        TESTS::

            sage: A = HopfAlgebraOnOrderedMultisetPartitions(QQ)
            sage: TestSuite(A).run()  # long time
        """
        if alphabet is not None:
            #from warnings import warn
            #warn("alphabet keyword is not yet implemented", RuntimeWarning)
            if alphabet in ZZ:
                self._A = Set(range(1,alphabet+1))
            else:
                self._A = Set(alphabet)
            if self._A == Set() or any(a not in ZZ for a in self._A):
                raise ValueError("keyword alphabet was converted to %s, which must be a nonempty set of positive integers"%(Set(self._A)))
        else:
            self._A = None
        #category = HopfAlgebras(R).Graded().Connected()  # TODO: add Cocommutative
        category = GradedHopfAlgebras(R).Connected()  # TODO: add Cocommutative
        Parent.__init__(self, base=R, category=category.WithRealizations())
        if self._A:
            self._order_grading = order_grading
        else:
            self._order_grading = False

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: HopfAlgebraOnOrderedMultisetPartitions(ZZ)
            Hopf Algebra on Ordered Multiset Partitions over Integer Ring
            sage: HopfAlgebraOnOrderedMultisetPartitions(ZZ, alphabet=[2,3])
            Hopf Algebra on Ordered Multiset Partitions on alphabet {2, 3} over Integer Ring with order grading
            sage: HopfAlgebraOnOrderedMultisetPartitions(ZZ, alphabet=[2,3], order_grading=False)
            Hopf Algebra on Ordered Multiset Partitions on alphabet {2, 3} over Integer Ring
        """
        OMP = "Hopf Algebra on Ordered Multiset Partitions"
        if self._A:
            OMP += " on alphabet %s"%Set(self._A)
        OMP += " over the %s"%self.base_ring()
        if self._order_grading:
            OMP += " with order grading"
        return OMP

    def a_realization(self):
        r"""
        Return a particular realization of ``self`` (the Homogeneous basis).

        EXAMPLES::

            sage: HopfAlgebraOnOrderedMultisetPartitions(QQ).a_realization()
            Hopf Algebra on Ordered Multiset Partitions over Rational Field in the Homogeneous basis
        """
        return self.H()

    _shorthands = tuple(['H', 'P'])

    def dual(self):
        r"""
        Return the dual Hopf algebra to Hopf Algebra on Ordered Multiset Partitions.

        EXAMPLES::

            sage: HopfAlgebraOnOrderedMultisetPartitions(QQ).dual()
            Dual to Hopf Algebra on Ordered Multiset Partitions over the Rational Field
        """
        return HopfAlgebraOnOrderedMultisetPartitionsDual(self.base_ring(), self._A, self._order_grading)

    class H(OMPBasis_abstract):
        r"""
        The Homogeneous basis of Hopf Algebra on Ordered Multiset Partitions.

        The family `(H_\mu)`, as `\mu` ranges over all ordered multiset
        partitions of nonnegative integers, is called the *Homogeneous basis*
        here. It is the natural "free algebra on finite sets" basis for
        Hopf Algebra on Ordered Multiset Partitions.

        EXAMPLES::

            sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
            sage: H[[2], [1,3]] - 2*H[[1,2,3]]
            H[{2}, {1,3}] - 2*H[{1,2,3}]
            sage: H[[2,3]].coproduct()
            H[] # H[{2,3}] + H[{2}] # H[{3}] + H[{3}] # H[{2}] + H[{2,3}] # H[]

        TESTS::

            sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
            sage: TestSuite(H).run()  # long time
        """
        _prefix = "H"
        _basis_name = "Homogeneous"

        def dual_basis(self):
            """
            Return the dual basis, Monomial basis of Dual to
            Hopf Algebra on Ordered Multiset Partitions.
            """
            return self.realization_of().dual().M()

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

                sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).Homogeneous()
                sage: H[{2,3}] * H[{2,3}, {2}]
                H[{2,3}, {2,3}, {2}]

            TESTS::

                sage: one = OrderedMultisetPartition([])
                sage: all(H.product_on_basis(one, z) == H(z) == H.basis()[z] for z in OrderedMultisetPartitions(3))
                True
                sage: all(H.product_on_basis(z, one) == H(z) == H.basis()[z] for z in OrderedMultisetPartitions(3))
                True
            """
            # TODO: if alphabet is not None, then one should check that B has valid blocks.
            return self.monomial(A + B)

        def coproduct_on_basis(self, A):
            r"""
            Return the coproduct, defined on subsets `K` via
            .. MATH::

                \sum_{I\sqcup J = K I \otimes J,

            then extended multiplicatively and linearly. This coproduct is cocommutative.
            It is analogous to the Complete basis ring in the Hopf algebra of
            :class:`non-commutative symmetric functions<NonCommutativeSymmetricFunctions>`.

            INPUT:

            - ``A`` -- a set partition

            OUTPUT:

            - The coproduct applied to the element of OMP indexed by ``A`` expressed in the
              Homogeneous basis.

            EXAMPLES::

                sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).Homogeneous()
                sage: H[{2,3}].coproduct()
                H[] # H[{2,3}] + H[{2}] # H[{3}] + H[{3}] # H[{2}] + H[{2,3}] # H[]
                sage: H[{2,3}, {2}].coproduct()
                H[] # H[{2,3}, {2}] + H[{2}] # H[{3}, {2}] + H[{3}] # H[{2}, {2}] + H[{2,3}] # H[{2}] + H[{2}] # H[{2,3}] + H[{2}, {2}] # H[{3}] + H[{3}, {2}] # H[{2}] + H[{2,3}, {2}] # H[]
            """
            return self.tensor_square()._from_dict( A.split(2) )

        def antipode_on_basis(self, A):
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

            EXAMPLES::

                sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).Homogeneous()
                sage: H[{2,3}].antipode()
                H[{2}, {3}] + H[{3}, {2}] - H[{2,3}]
                sage: H[{2},{3}].antipode()
                H[{3},{2}]
                sage: H[{2,3,4}].antipode()
                H[{2,4}, {3}] + H[{3}, {2,4}] - H[{2}, {4}, {3}] - H[{2}, {3}, {4}]
                 - H[{3}, {2}, {4}] - H[{4}, {3}, {2}] + H[{3,4}, {2}]
                 + H[{2}, {3,4}] - H[{4}, {2}, {3}] + H[{4}, {2,3}]
                 - H[{3}, {4}, {2}] + H[{2,3}, {4}] - H[{2,3,4}]
                sage: HH = H[[1,2,3,5],[1,3]].coproduct()
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

        class Element(CombinatorialFreeModule.Element):
            r"""
            An element in the `\mathbf{M}` basis of OMP^*.
            """
            def to_symmetric_function(self):
                r"""
                The projection of ``self`` to the ring of symmetric functions.
                """
                return self.to_noncommutative_symmetric_function().to_symmetric_function()

            def to_noncommutative_symmetric_function(self):
                r"""
                The projection of ``self`` to the ring `NSym` of
                noncommutative symmetric functions.

                .. TODO::

                    Fix what follows (docstring and code).

                There is a canonical projection `\pi : OMPSym \to NSym`
                that sends ...
                This `\pi` is a Hopf algebra homomorphism.

                OUTPUT:

                - an element of the Hopf algebra of noncommutative symmetric functions
                  in the Complete basis

                EXAMPLES::

                    sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
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
                from sage.combinat.ncsf_qsym.ncsf import NonCommutativeSymmetricFunctions
                S = NonCommutativeSymmetricFunctions(self.parent().base_ring()).Complete()
                H = self.parent().realization_of().H()
                return S.sum_of_terms((A.shape_from_cardinality(), coeff/prod(factorial(len(a)) for a in A))
                                      for (A, coeff) in H(self))

    Homogeneous = H

    class P(OMPBasis_abstract):
        r"""
        The Powersum basis of Hopf Algebra on Ordered Multiset Partitions.

        The family `(P_\mu)`, as `\mu` ranges over all ordered multiset
        partitions of nonnegative integers, is called the *Powersum basis*
        here. It is defined via a unitriangular change of basis from the
        Homogeneous basis of `OMPSym` via *some partial order* on
        ordered multiset partitions.

        The principle feature of the Powersum basis is that for subsets `K`,
        `P_{[K]}` is primitive.

        EXAMPLES::

            sage: P = HopfAlgebraOnOrderedMultisetPartitions(QQ).Powersum()
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
            sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).Homogeneous()
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
            OMPBasis_abstract.__init__(self, alg, True)

            # Register coercions
            H = self.realization_of().H()
            phi = self.module_morphism(self._P_to_H_on_basis, codomain=H, unitriangular="upper")
            # TODO: check the above unitriangular claim
            phi.register_as_coercion()
            (~phi).register_as_coercion() # TODO: figure out what tilde-phi signifies

        #@cached_method
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
                - Should we cache this method?

            EXAMPLES::

                sage: P = HopfAlgebraOnOrderedMultisetPartitions(QQ).P()
                sage: A = OrderedMultisetPartition([[1,2,3]])
                sage: P._P_to_H_on_basis(A)
                ???
                sage: B = OrderedMultisetPartition([[1,2],[3]])
                sage: P._P_to_H_on_basis(B)
                ???
            """
            OMP = self.basis().keys()
            def lt(x, y):
                if x == y:
                    return False
                return OMP(x).is_finer(OMP(y))

            H = self.realization_of().H()
            P_refine = Poset((A.finer(), lt))
            print P_refine
            c = abs(prod((-1)**(i-1) * factorial(i-1) for i in A.shape_from_size()))
            R = self.base_ring()
            print {B: R(P_refine.moebius_function(B, A) / ZZ(c))
                                 for B in P_refine}
            return H._from_dict({OMP(B): R(P_refine.moebius_function(B, A) / ZZ(c))
                                 for B in P_refine}) #, remove_zeros=False)

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

                sage: P = HopfAlgebraOnOrderedMultisetPartitions(QQ).Powersum()
                sage: elt = P.primitive(OrderedMultisetPartition([[1,3], [1,2]])); elt
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
                P = OrderedMultisetPartitions()
                if i not in A[0]:
                    i = A[0][0]
                return self.sum_of_terms((finer, (-1)**len(finer))
                                    for finer in A.finer() if i in finer[0])

    Powersum = P

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

            sage: from sage.combinat.chas.omp import OMPBases
            sage: OMPSym = HopfAlgebraOnOrderedMultisetPartitions(ZZ)
            sage: bases = OMPBases(OMPSym, True)
            sage: OMPSym.H() in bases
            True
        """
        # TODO: do I need the next two lines?
        #self._A = base._A
        #self._order_grading = base._order_grading
        self._graded = graded
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.chas.omp import OMPBases
            sage: OMPSym = HopfAlgebraOnOrderedMultisetPartitions(ZZ)
            sage: OMPBases(OMPSym, True)
            Category of graded bases of Hopf Algebra on Ordered Multiset Partitions over Integer Ring
            sage: WQSymBases(OMPSym, False)
            Category of filtered bases of Hopf Algebra on Ordered Multiset Partitions over Integer Ring
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

            sage: from sage.combinat.chas.omp import OMPBases
            sage: OMPSym = HopfAlgebraOnOrderedMultisetPartitions(ZZ)
            sage: bases = OMPBases(OMPSym, True)
            sage: bases.super_categories()
            [Category of realizations of Hopf Algebra on Ordered Multiset Partitions
                 over Integer Ring,
             Join of Category of realizations of hopf algebras over Integer Ring
                 and Category of graded algebras over Integer Ring,
             Category of graded connected hopf algebras with basis over Integer Ring]

            sage: bases = OMPBases(OMPSym, True)
            sage: bases.super_categories()
            [Category of realizations of Hopf Algebra on Ordered Multiset Partitions
                 over Integer Ring,
             Join of Category of realizations of hopf algebras over Integer Ring
                 and Category of graded algebras over Integer Ring,
             Join of Category of filtered connected hopf algebras with basis over Integer Ring
                 and Category of graded algebras over Integer Ring]
        """
        R = self.base().base_ring()
        #cat = HopfAlgebras(R).Graded().WithBasis()
        cat = GradedHopfAlgebras(R).WithBasis()
        if self._graded:
            cat = cat.Graded()
        else:
            cat = cat.Filtered()
        return [self.base().Realizations(),
                #HopfAlgebras(R).Graded().Realizations(),
                GradedHopfAlgebras(R).Realizations(),
                cat.Connected()]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of `OMPSym` or `OMPQSym`.

            EXAMPLES::

                sage: OMPSym = HopfAlgebraOnOrderedMultisetPartitions(ZZ)
                sage: OMPSym.H()
                Hopf Algebra on Ordered Multiset Partitions over Integer Ring in the Homogeneous basis
                sage: OMPSymA = HopfAlgebraOnOrderedMultisetPartitions(ZZ, alphabet=[1,4])
                sage: OMPSymA.H()
                Hopf Algebra on Ordered Multiset Partitions on alphabet {1, 4} over Integer Ring with order grading in the Homogeneous basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def __getitem__(self, A):
            """
            Return the basis element indexed by ``A``.

            INPUT:

            - ``A`` -- an ordered multiset partition

            .. TODO::

                - Fix the output below. It should be ``H[{3,4,5}]``:
                    sage: H[[[3,4,5]]]
                    H[{3}, {4}, {5}]
                - Allow for input such as `H[3,4,0,5]` and `H[3,4,5]`?

            EXAMPLES::

                sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
                sage: H[[1, 3, 2]]
                H[{1,2,3}]
                sage: H[[1,3],[1]]
                H[{1,3}, {1}]
                sage: H[OrderedMultisetPartition([2,0,1,4,0,3,5])]
                H[{2}, {1,4}, {3,5}]
            """
            try:
                return self.monomial(self._indices(A))
            except TypeError:
                return self.monomial(self._indices([A]))

        def is_field(self, proof=True):
            """
            Return whether ``self`` is a field.

            EXAMPLES::

                sage: H = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
                sage: H.is_field()
                False
            """
            return False

        def is_commutative(self):
            """
            Return whether ``self`` is commutative.

            EXAMPLES::

                sage: H = HopfAlgebraOnOrderedMultisetPartitions(ZZ).H()
                sage: H.is_commutative()
                False
            """
            return self.base_ring().is_zero()

        def one_basis(self):
            """
            Return the index of the unit.

            EXAMPLES::

                sage: A = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
                sage: A.one_basis()
                []
            """
            P = self.basis().keys()
            return P([])

        def degree_on_basis(self, A):
            """
            Return the degree of an ordered multiset partition ``A`` in
            Hopf Algebra of Ordered Multiset Partitions or its dual.

            If order grading is used, then the degree is ``A.order()``.
            Else, the degree is ``A.size()``.

            EXAMPLES::

                sage: HA = HopfAlgebraOnOrderedMultisetPartitions(QQ).H()
                sage: co = OrderedMultisetPartition([[2,1],[1,4]])
                sage: co, co.size(), co.order()
                ([{1,2}, {1,4}], 8, 4)
                sage: HA.degree_on_basis(co)
                8
                sage: HB = HopfAlgebraOnOrderedMultisetPartitions(QQ, alphabet=[1,2,4]).H()
                sage: HB.degree_on_basis(co)
                4
                sage: HC = HopfAlgebraOnOrderedMultisetPartitions(QQ, alphabet=[1,2,4], order_grading=False).H()
                sage: HC.degree_on_basis(co)
                8
            """
            # TODO: should we keep next line?
            #A = self.basis().keys()(A)
            if self._order_grading:
                return A.order()
            else:
                return A.size()

    class ElementMethods:
        def duality_pairing(self, x):
            r"""
            Compute the pairing between element ``self`` and an
            element of the dual.

            INPUT:

            - ``x`` -- an element of ``self.parent().dual_basis()``.

            OUTPUT:

            - an element of the base ring of ``self.parent()``

            EXAMPLES::

                sage: DOMP = HopfAlgebraOnOrderedMultisetPartitionsDual(QQ)
                sage: M = DOMP.M()
                sage: h = M.dual_basis()
                sage: matrix([[M(A).duality_pairing(H(B)) for A in OrderedMultisetPartitions(3)] for B in OrderedMultisetPartitions(3)])
                [1 0 0 0 0]
                [0 1 0 0 0]
                [0 0 1 0 0]
                [0 0 0 1 0]
                [0 0 0 0 1]
                sage: (M[[1,2],[3]] + 3*M[[1,3],[2]]).duality_pairing(2*H[[1,3],[2]] + H[[1,2,3]] + 2*H[[1,2],[3]])
                8
                sage: P = HopfAlgebraOnOrderedMultisetPartitions(QQ).P()
                sage: matrix([[M(A).duality_pairing(P(B)) for A in OrderedMultisetPartitions(3)] for B in OrderedMultisetPartitions(3)])
                [_ _ _ _ _]
                [_ _ _ _ _]
                [_ _ _ _ _]
                [_ _ _ _ _]
                [_ _ _ _ _]
                sage: (2*M[[2,3],[3]] + 3*M[[3],[2],[3]] + M[[3],[2,3]]).duality_pairing(P[[2,3],[3]] + 3*P[[3] [2,3]])
                ???
            """
            P = self.parent()
            x = P.dual_basis()(x)
            return sum(coeff * x[A] for (A, coeff) in self)

############################
class HopfAlgebraOnOrderedMultisetPartitionsDual(UniqueRepresentation, Parent):
    r"""
    The Hopf dual to Hopf Algebra on Ordered Multiset Partitions.

    See Example 14 of [LM2018]_ for its dual.
    """
    @staticmethod
    def __classcall_private__(cls, R, alphabet=None, order_grading=True):
        """
        Normalize the input to ensure a unique representation.

        TESTS::

            sage: OMPQ1 = HopfAlgebraOnOrderedMultisetPartitionsDual(QQ, alphabet=[2,3,5])
            sage: OMPQ2 = HopfAlgebraOnOrderedMultisetPartitionsDual(QQ, alphabet=Set([2,3,5]))
            sage: OMPQ1 is OMPQ2
            True
        """
        if alphabet is not None:
            if alphabet in ZZ:
                alphabet = Set(range(1,alphabet+1))
            else:
                alphabet = Set(alphabet)
            if alphabet == Set() or any(a not in ZZ for a in alphabet):
                raise ValueError("keyword alphabet was converted to %s, which must be a nonempty set of positive integers"%(Set(alphabet)))
        else:
            alphabet = None
        return super(HopfAlgebraOnOrderedMultisetPartitionsDual, cls).__classcall__(cls, R, alphabet, order_grading)

    def __init__(self, R, alphabet=None, order_grading=True):
        """
        Initialize ``self``.

        TESTS::

            sage: A = HopfAlgebraOnOrderedMultisetPartitionsDual(QQ)
            sage: TestSuite(A).run()  # long time
        """
        if alphabet is not None:
            #from warnings import warn
            #warn("alphabet keyword is not yet implemented", RuntimeWarning)
            if alphabet in ZZ:
                self._A = Set(range(1,alphabet+1))
            else:
                self._A = Set(alphabet)
            if self._A == Set() or any(a not in ZZ for a in self._A):
                raise ValueError("keyword alphabet was converted to %s, which must be a nonempty set of positive integers"%(Set(self._A)))
        else:
            self._A = None
        category = GradedHopfAlgebras(R).Commutative().Connected()
        Parent.__init__(self, base=R, category=category.WithRealizations())
        if self._A:
            self._order_grading = order_grading
        else:
            self._order_grading = False


    def _repr_(self):
        r"""
        EXAMPLES::

            sage: HopfAlgebraOnOrderedMultisetPartitions(ZZ).dual()
            Dual to Hopf Algebra on Ordered Multiset Partitions over the Integer Ring
        """
        return "Dual to Hopf Algebra on Ordered Multiset Partitions over the %s"%self.base_ring()

    def a_realization(self):
        r"""
        Return the realization of the `\mathbf{M}` basis of ``self``.

        EXAMPLES::

            sage: HopfAlgebraOnOrderedMultisetPartitions(QQ).dual().a_realization()
            Dual to Hopf Algebra on Ordered Multiset Partitions over the Rational Field in the m basis
        """
        return self.M()

    _shorthands = tuple(['M'])

    def dual(self):
        r"""
        Return the dual Hopf algebra to Hopf Algebra on Ordered Multiset Partitions.

        EXAMPLES::

            sage: OMPD = HopfAlgebraOnOrderedMultisetPartitions(QQ).dual()
            sage: OMPD.dual()
            Hopf Algebra on Ordered Multiset Partitions over the Rational Field
        """
        return HopfAlgebraOnOrderedMultisetPartitions(self.base_ring())

    class M(OMPBasis_abstract):
        r"""
        The Dual to Hopf Algebra on Ordered Multiset Partitions in the `\mathbf{M}` basis.

        EXAMPLES::

            sage: OMPD = HopfAlgebraOnOrderedMultisetPartitions(QQ).dual()
            sage: M = OMPD.M()

        """
        _prefix = "M"
        _basis_name = "Monomial"

        def dual_basis(self):
            """
            Return the dual basis, Homogeneous basis of
            Hopf Algebra on Ordered Multiset Partitions.
            """
            return self.realization_of().dual().H()

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

            - ``A``, ``B`` -- ordered multiset partitions

            OUTPUT:

            - an element of the `\mathbf{M}` basis

            EXAMPLES::

                sage: M = HopfAlgebraOnOrderedMultisetPartitions(QQ).dual().M()
                sage: A = OrderedMultisetPartition(from_zero_list=[2,1,3,0,1,2]); A
                [{1,2,3}, {1,2}]
                sage: B = OrderedMultisetPartition([[3,4]]); B
                [{3,4}]
                sage: M.product_on_basis(A, B)
                M[{1,2,3}, {1,2}, {3,4}] + M[{1,2,3}, {3,4}, {1,2}] + M[{3,4}, {1,2,3}, {1,2}] + M[{1,2,3}, {1,2,3,4}]
                sage: C = OrderedMultisetPartition([[4,5]]); C
                [{4,5}]
                sage: M.product_on_basis(A, C)
                M[{1,2,3}, {1,2}, {4,5}] + M[{1,2,3}, {4,5}, {1,2}] + M[{4,5}, {1,2,3}, {1,2}] + M[{1,2,3}, {1,2,4,5}] + M[{1,2,3,4,5}, {1,2}]
                sage: M.product_on_basis(A, OrderedMultisetPartition([]))
                M[{1,2,3}, {1,2}]
            """
            # TODO: if alphabet is not None, then one should check that B has valid blocks.
            terms = A.shuffle_product(B, overlap=True)
            return self.sum_of_terms([(s, 1) for s in terms])

        def coproduct_on_basis(self, A):
            r"""
            Return the coproduct of a `\mathbf{M}` basis element.

            The coproduct on the basis element `\mathbf{M}_A` is the sum over
            tensor product terms `\mathbf{M}_B \otimes \mathbf{M}_C` where
            `A=B+C`. (Here `+` is concatenation of ordered multiset partitions.

            INPUT:

            - ``A`` -- a ordered multiset partition

            OUTPUT:

            - The coproduct applied to the `\mathbf{M}` basis element indexed by ``A``
              expressed in the `\mathbf{M}` basis.

            EXAMPLES::

                sage: M = HopfAlgebraOnOrderedMultisetPartitions(QQ).dual().M()
                sage: M[[2], [2,3]].coproduct()
                M[] # M[{2}, {2,3}] + M[{2}] # M[{2,3}] + M[{2}, {2,3}] # M[]
                sage: M.coproduct_on_basis(OrderedMultisetPartition([]))
                M[] # M[]
            """
            return self.tensor_square().sum_of_monomials(A.deconcatenate(2), distinct=True)

        def _antipode_on_basis(self, A):
            r"""
            Return the antipode applied to the basis element indexed by ``A``.

            .. TODO::

                Develop a formula for the antipode that is faster than
                the default implementation.

            INPUT:

            - ``A`` -- an ordered multiset partition

            OUTPUT:

            - an element in the basis ``self``
            """
            pass

        def sum_of_derangements(self, A):
            """
            Return the sum over all ordered multiset partitions whose multiset of blocks
            coincide with those of `A`.

            INPUT:

            - ``A`` -- a ordered multiset partition

            OUTPUT:

            - an element of ``self``

            EXAMPLES::

                sage: M = HopfAlgebraOnOrderedMultisetPartitions(QQ).dual().M()
                sage: M.sum_of_partitions([[2,1],[1]])
                M[{1}, {1,2}] + M[{1,2}, {1}]
            """
            A = OrderedMultisetPartitions(A)
            P = OrderedMultisetPartitions
            return self.sum_of_terms([(P(omp), 1) for omp in Permutations(list(A))], distinct=True)

        def sum_of_shape_from_cardinality(self, A, deg=None):
            """
            Return the sum over all ordered multiset partitions of degree ``deg`` whose underlying composition shape (from cardinality) is ``A``.

            INPUT:

            - ``deg`` -- a nonnegative integer
            - ``A`` -- a composition

            OUTPUT:

            - an element of ``self``

            EXAMPLES::

                sage: M = HopfAlgebraOnOrderedMultisetPartitions(QQ).dual().M()
                sage: M.sum_of_shape_from_cardinality(5,[2,1])
                M[{1,3}, {1}] + M[{1,2}, {2}]
                sage: M.sum_of_shape_from_cardinality(6,[2,1])
                M[{2,3}, {1}] + M[{1,3}, {2}] + M[{1,4}, {1}] + M[{1,2}, {3}]
            """
            # NOTE: This function is *really* meant to be used with ``deg`` parameter specified.
            if deg is None:
                deg = sum([binomial(k,2) for k in A])
            A = Composition(A)
            P = OrderedMultisetPartitions(deg)
            return self.sum_of_terms([(omp, 1) for omp in P if P.shape_from_cardinality() == A], distinct=True)

        def sum_of_shape_from_size(self, A):
            """
            Return the sum over all ordered multiset partitions of |A| whose underlying composition shape (from size) is ``A``.

            INPUT:

            - ``A`` -- a composition

            OUTPUT:

            - an element of ``self``

            EXAMPLES::

                sage: M = HopfAlgebraOnOrderedMultisetPartitions(QQ).dual().M()
                sage: M.sum_of_shape_from_cardinality([3,1])
                M[{1,2}, {1}] + M[{3}, {1}]
                sage: M.sum_of_shape_from_cardinality([6,2])
                M[{1,2,3}, {2}] + M[{1,5}, {2}] + M[{2,4}, {2}] + M[{6}, {2}]
            """
            A = Composition(A)
            P = OrderedMultisetPartitions(sum(A))
            return self.sum_of_terms([(omp, 1) for omp in P if P.shape_from_size() == A], distinct=True)

        class Element(CombinatorialFreeModule.Element):
            r"""
            An element in the `\mathbf{M}` basis of OMP^*.
            """
            def is_symmetric(self):
                r"""
                Meant to model the inclusion of SYM inside QSYM, though there is no notion of
                quasisymmetric polynomials related to Dual to Hopf Algebra on Ordered Multiset Partitions
                as of yet.

                Determine if a `OMP^*` function, expressed in the
                `\mathbf{M}` basis, is symmetric.

                A function `f` in the `\mathbf{M}` basis shall be deemed symmetric
                if for each ordered multiset partition `A` in its support, all derangements of `A`
                are also in the support of `f` with the same coefficient.

                .. MATH::

                    f = \sum_{\lambda} c_{\lambda} \prod_i m_i(\lambda)!
                    \sum_{\lambda(A) = \lambda} \mathbf{M}_A

                where the second sum is over all ordered multiset partitions `A` whose
                shape `\lambda(A)` is equal to `\lambda` and `m_i(\mu)` is
                the multiplicity of `i` in the partition `\mu`.

                OUTPUT:

                - ``True`` if `sorted(A)=sorted(B)` implies the coefficients of
                  `\mathbf{M}_A` and `\mathbf{M}_B` are equal, ``False`` otherwise

                EXAMPLES::

                    sage: M = HopfAlgebraOnOrderedMultisetPartitions(QQ).dual().M()
                    sage: elt = M.sum_of_derangements([[2],[2,3],[1]])
                    sage: elt.is_symmetric()
                    True
                    sage: elt -= 3*M.sum_of_partitions([[1],[1]])
                    sage: elt.is_symmetric()
                    True
                    sage: elt += M([[2],[1],[2,3]])
                    sage: elt.is_symmetric()
                    False
                    sage: elt = M[[1,3],[2]]
                    sage: elt.is_symmetric()
                    False
                    sage: elt = M[[1,3],[2]] + M[[2],[1,3]] + 2* M[[2,3],[2,3]]
                    sage: elt.is_symmetric()
                    True
                """
                P = OrderedMultisetPartitions()
                d = {}
                R = self.base_ring()
                for A, coeff in self:
                    la = P(sorted(A))
                    if la not in d:
                        d[la] = [coeff, 1]
                    else:
                        if d[la][0] != coeff:
                            return False
                        d[la][1] += 1
                # Make sure we've seen each ordered multiset partition in the derangement class
                return all(d[la][1] == Permutations(la).cardinality() for la in d)

    Monomial = M

############################
OMPSym = HopfAlgebraOnOrderedMultisetPartitions(QQ)
H = OMPSym.H()
hx = H[[2,3],[2,1]]
OP = OrderedMultisetPartitions()
print H
print "hx =", hx
print
print "hx^2 = ", hx * hx
print "hx.coproduct() = ", hx.coproduct()
print
print "hx.antipode()", hx.antipode()
print
# OMPSymA = HopfAlgebraOnOrderedMultisetPartitions(QQ, alphabet=[2,3,5])
# HA = OMPSymA.H()
# OPA = OrderedMultisetPartitions(alphabet=[2,3,5])
# print HA
# hy = HA[[2,3],[2,5]]
# print "hy =", hy
# print "hy^2 = ", hy * hy
# print "hy.coproduct() = ", hy.coproduct()
# print
# print "hy.antipode()", hy.antipode()
# print "hy._antipode()", hy._antipode()
# print
# DOMPSym = HopfAlgebraOnOrderedMultisetPartitionsDual(QQ)
# M = DOMPSym.M()
# mx = M[[2,3],[2,1]]
# OP = OrderedMultisetPartitions()
# print M
# print "mx =", mx
# print "mx^2 = ", mx * mx
# print "mx.coproduct() = ", mx.coproduct()
# #print
# #print "mx.antipode()", mx.antipode()
# #print "mx._antipode()", mx._antipode()
# print
# DOMPSymA = HopfAlgebraOnOrderedMultisetPartitionsDual(QQ, alphabet=[2,3,5])
# MA = DOMPSymA.M()
# my = MA[[2,3],[2,5]]
# OPA = OrderedMultisetPartitions(alphabet=[2,3,5])
# print MA
# print "my =", my
# print "my^2 = ", my * my
# print "my.coproduct() = ", my.coproduct()
# #print
# #print "my.antipode()", my.antipode()
# #print "my._antipode()", my._antipode()
# print
