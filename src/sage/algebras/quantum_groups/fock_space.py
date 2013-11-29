r"""
Fock Space

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.categories.category_o import CategoryOInt # Should be quantum version
from sage.categories.realizations import Realizations, Category_realization_of_parent

from sage.rings.all import ZZ, QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.combinat.partition import Partition, Partitions, _Partitions
from sage.combinat.partition_tuple import PartitionTuples
from sage.combinat.root_system.root_system import RootSystem
from sage.algebras.quantum_groups.quantum_group import QuantumGroup
from sage.algebras.quantum_groups.q_numbers import q_factorial

class FockSpace(CombinatorialFreeModule):
    r"""
    The Fock space of `U_q(\widehat{\mathfrak{sl}}_n)` of level `1` of
    weight `\Lambda_r`.
    """
    @staticmethod
    def __classcall_private__(cls, n, r=[0], q=None, base_ring=None):
        """
        Standardize input to ensure a unique representation.
        """
        if q is None:
            base_ring = PolynomialRing(QQ, 'q')
            q = base_ring.gen(0)
        if base_ring is None:
            base_ring = q.parent()
        base_ring = FractionField(base_ring)
        q = base_ring(q)
        M = IntegerModRing(n)
        if r in ZZ:
            r = [r]
        r = map(M, r)
        return super(FockSpace, cls).__classcall__(cls, n, tuple(r), q, base_ring)

    def __init__(self, n, r, q, base_ring):
        """
        Initialize ``self``.
        """
        self._n = n
        self._q = q
        self._r = r
        # If the cell x is above the cell y
        if len(r) == 1: # For partitions
            self._above = lambda x,y: x[1] > y[1]
        else: # For partition tuples
            self._above = lambda x,y: x[0] > y[0] or (x[0] == y[0] and x[1] > y[1])
        self._q_group = QuantumGroup(['A', n-1, 1], q, R=base_ring)
        self._addable = lambda la,i: filter(lambda x: la.content(*x, multicharge=self._r) == i,
                                            la.outside_corners())
        self._removable = lambda la,i: filter(lambda x: la.content(*x, multicharge=self._r) == i,
                                              la.corners())
        CombinatorialFreeModule.__init__(self, base_ring, PartitionTuples(len(r)),
                                         prefix='', bracket=['|', '>'],
                                         latex_bracket=['\\lvert', '\\rangle'],
                                         monomial_cmp=lambda x,y: -cmp(x,y),
                                         category=CategoryOInt(self._q_group))

    def __getitem__(self, i):
        """
        Return the basis element indexed by ``i``.

        INPUT:

        - ``i`` -- a partition

        EXAMPLES::

            sage: F = FockSpace(2)
            sage: F[[]]
            |[]>
            sage: F[1]
            |[1]>
            sage: F[2,2,1]
            |[2, 2, 1]>
        """
        if i in ZZ:
            i = [i]
        return self.monomial(self._indices(i))

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Fock space of rank {} of type {} over {}".format(self._n, self._r, self.base_ring())

    def basic_representation(self):
        """
        Return the basic representation `B(\Lambda_r)` in ``self``.
        """
        if len(self._r) != 1:
            raise NotImplementedError("only implemented for level 1")
        return BasicRepresentation(self)

    def module_generator(self):
        """
        Return the module generator of ``self``.
        """
        if len(self._r) == 1:
            return self.monomial(_Partitions([]))
        return self.monomial( self._indices([[]]*len(self._r)) )

    class Element(CombinatorialFreeModuleElement):
        def e(self, i):
            """
            Apply the action of `e_i` on ``self``.
            """
            P = self.parent()
            N_left = lambda la, x, i: \
                sum(1 for y in P._addable(la, i) if P._above(x, y)) \
                - sum(1 for y in P._removable(la, i) if P._above(x, y))
            q = P._q
            return P.sum_of_terms( [( la.remove_cell(*x), c * q**(-N_left(la, x, i)) )
                                    for la,c in self for x in P._removable(la, i)] )

        def f(self, i):
            """
            Apply the action of `f_i` on ``self``.
            """
            P = self.parent()
            N_right = lambda la, x, i: \
                sum(1 for y in P._addable(la, i) if P._above(y, x)) \
                - sum(1 for y in P._removable(la, i) if P._above(y, x))
            q = P._q
            return P.sum_of_terms( [(la.add_cell(*x), c * q**N_right(la, x, i))
                                for la,c in self for x in P._addable(la, i)] )

        def h(self, i):
            """
            Apply the action of `h_i` on ``self``.
            """
            P = self.parent()
            N_i = lambda la, i: len(P._addable(la, i)) - len(P._removable(la, i))
            q = P._q
            return P.sum_of_terms([(la, c * q**N_i(la, i) * P.monomial(la))
                                   for la,c in self])

        def d(self):
            """
            Apply the action of `d` on ``self``.
            """
            q = self.parent()._q
            return self.parent().sum_of_terms(
                [(la, c * q**sum(1 for x in la.cells() if la.content(*x, multicharge=self._r) == 0))
                 for la,c in self])

        def _l_action_(self, g):
            """
            Return the left action of ``g`` on ``self``.
            """
            if g not in self._q_group:
                raise TypeError("{0} must be an element of {1}".format(g, self._q_group))

            s = self.zero()
            for cg, monomial in g:
                c = self
                for x in monomial.variables():
                    if x in self._q_group._e:
                        c = c.e(self._q_group._e.index(x))
                    elif x in self._q_group._f:
                        c = c.f(self._q_group._f.index(x))
                    elif x in self._q_group._h:
                        c = c.e(self._q_group._h.index(x))
                    else: # d
                        c = c.d()
                s += cg * c
            return s

class BasicRepresentation(Parent, UniqueRepresentation):
    """
    Basic representation `B(\Lambda_0)` of `U_q(\widehat{\mathfrak{sl}}_n)`.
    """
    @staticmethod
    def __classcall_private__(cls, n, r=0, q=None, base_ring=None):
        """
        Standardize input to ensure a unique representation.
        """
        if isinstance(n, FockSpace):
            return super(BasicRepresentation, cls).__classcall__(cls, n)
        return super(BasicRepresentation, cls).__classcall__(cls, FockSpace(n, [r], q, base_ring))

    def __init__(self, fock_space):
        """
        Initialize ``self``.
        """
        self._n = fock_space._n
        self._q = fock_space._q
        self._fock = fock_space
        Parent.__init__(self, base=self._fock.base_ring(),
                        category=CategoryOInt(self._fock._q_group).WithRealizations())

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        ct = self._fock._q_group.cartan_type()
        P = ct.root_system().weight_lattice()
        wt = P.sum_of_monomials(self._fock._r)
        return "Basic representation of {} of weight {}".format(ct, wt)

    def a_realization(self):
        """
        Return a realization of ``self``.
        """
        return self.G()

    class A(CombinatorialFreeModule, BindableClass):
        """
        The `A` basis of the Fock space which is the approximation of the
        lower global crystal basis.
        """
        def __init__(self, basic):
            """
            Initialize ``self``.
            """
            self._basis_name = "approximation"
            self._n = basic._n
            self._q = basic._q
            CombinatorialFreeModule.__init__(self, self._q.parent(), _Partitions,
                                             prefix='A', bracket=False,
                                             monomial_cmp=lambda x,y: -cmp(x,y),
                                             category=BasicRepresentationBases(basic))
            self.module_morphism(self._A_to_fock_basis,
                                 triangular='upper', unitriangular=True,
                                 codomain=basic._fock).register_as_coercion()

        @cached_method
        def _A_to_fock_basis(self, la):
            """
            Return the `A` basis indexed by ``la`` in the natural basis.
            """
            fock = self.realization_of()._fock
            cur = fock.module_generator()
            corners = la.corners()
            cells = set(la.cells())
            q = self._q
            k = self._n - 1 # This is sl_{k+1}
            r = fock._r[0]
            b = ZZ.zero()
            while any(c[0]*k + c[1] >= b for c in corners): # While there is some cell left to count
                power = 0
                i = -b + r # This will be converted to a mod n number
                for x in range(0, b // k + 1):
                    if (x, b-x*k) in cells:
                        power += 1
                        cur = cur.f(i)
                cur /= q_factorial(power, q)
                b += 1
            return cur

    approximation = A

    class G(CombinatorialFreeModule, BindableClass):
        """
        The lower global crystal basis living inside of Fock space.
        """
        def __init__(self, basic):
            """
            Initialize ``self``.
            """
            self._basis_name = "lower global crystal"
            self._n = basic._n
            self._q = basic._q
            CombinatorialFreeModule.__init__(self, self._q.parent(), _Partitions,
                                             prefix='G', bracket=False,
                                             monomial_cmp=lambda x,y: -cmp(x,y),
                                             category=BasicRepresentationBases(basic))
            self.module_morphism(self._G_to_fock_basis,
                                 triangular='upper', unitriangular=True,
                                 codomain=basic._fock).register_as_coercion()

        @cached_method
        def _G_to_fock_basis(self, la):
            """
            Return the `G` basis indexed by ``la`` in the natural basis.
            """
            # Special case for the empty partition
            if len(la) == 0:
                return self.realization_of()._fock.one()
            cur = self.realization_of().A()._A_to_fock_basis(la)
            mu = la.next()
            while mu:
                if cur[mu] != 0:
                    d = cur[mu].denominator()
                    k = d.degree()
                    n = cur[mu].numerator()
                    if k != 0 or n.degree() == 0:
                        gamma = sum(n[i] * (self._q**(i-k) + self._q**(k-i))
                                    for i in range(min(n.degree(), k)))
                        gamma += n[k]
                        cur -= gamma * self._G_to_fock_basis(mu)
                mu = mu.next()
            return cur

        def _repr_(self):
            """
            Return a string representation of ``self``.
            """
            return "Lower global crystal basis of {}".format(self.realization_of()._fock)

    lower_global_crystal_basis = G

class BasicRepresentationBases(Category_realization_of_parent):
    r"""
    The category of bases of a basic representation.
    """
    def __init__(self, base):
        r"""
        Initialize the bases of a basic representation.

        INPUT:

        - ``base`` -- a basic representation

        TESTS::

            sage: from sage.algebras.quantum_groups.fock_space import BasicRepresentationBases
            sage: F = FockSpace(2)
            sage: bases = BasicRepresentationBases(F)
            sage: TestSuite(bases).run()
        """
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        r"""
        Returns the representation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.fock_space import BasicRepresentationBases
            sage: F = FockSpace(2)
            sage: BasicRepresentationBases(F)
            Category of bases of Fock space of type B2 in 1,-1 over Integer Ring
        """
        return "Category of bases of %s" % self.base()

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_groups.fock_space import BasicRepresentationBases
            sage: F = FockSpace(2)
            sage: bases = BasicRepresentationBases(F)
            sage: bases.super_categories()
            [Category of finite dimensional algebras with basis over Integer Ring,
             Category of realizations of Fock space of type B2 in 1,-1 over Integer Ring]
        """
        return [CategoryOInt(self.base()._fock._q_group), Realizations(self.base())]

    class ParentMethods:
        def _repr_(self):
            """
            Text representation of this basis of Fock space.

            EXAMPLES::

                sage: F = FockSpace(2)
                sage: F.natural_basis()
                Fock space of type B2 in 1,-1 over Integer Ring in the standard basis
                sage: F.A()
                Fock space of type B2 in 1,-1 over Integer Ring in the C basis
            """
            return "%s in the %s basis"%(self.realization_of(), self._basis_name)

        def highest_weight_vector(self):
            """
            Return the highest weight vector of ``self``.
            """
            return self.monomial(Partition([]))

        def __getitem__(self, i):
            """
            Return the basis element indexed by ``i``.

            INPUT:

            - ``i`` -- a partition

            EXAMPLES::

                sage: F = FockSpace(2)
                sage: A = F.A()
                sage: A[[]]
                A[]
                sage: A[3]
                A[3]
                sage: A[2,2,1]
                A[2, 2, 1]
                sage: G = F.G()
                sage: G[[]]
                G[]
                sage: G[3]
                G[3]
                sage: G[2,2,1]
                G[2, 2, 1]
            """
            if i in ZZ:
                i = [i]
            i = Partition(i)
            if max(i.to_exp()) >= self._n:
                raise ValueError("{} is not a {}-regular partition".format(i, self._n))
            return self.monomial(i)

