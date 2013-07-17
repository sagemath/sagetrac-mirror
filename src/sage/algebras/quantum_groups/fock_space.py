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
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.combinat.partition import Partition, Partitions
from sage.algebras.quantum_groups.quantum_group import QuantumGroup
from sage.algebras.quantum_groups.q_numbers import q_factorial

class FockSpace(CombinatorialFreeModule):
    """
    The Fock space of `U_q(\widehat{\mathfrak{sl}}_n)`.
    """
    @staticmethod
    def __classcall_private__(cls, n, q=None, base_ring=None):
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
        return super(FockSpace, cls).__classcall__(cls, n, q, base_ring)

    def __init__(self, n, q, base_ring):
        """
        Initialize ``self``.
        """
        self._n = n
        self._q = q
        self._q_group = QuantumGroup(['A', n-1, 1], q, R=base_ring)
        CombinatorialFreeModule.__init__(self, self._q_group.base_ring(), Partitions(),
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
            sage: N = F.natural_basis()
            sage: N[[]]
            |[]>
            sage: N[1]
            |[1]>
            sage: N[2,2,1]
            |[2, 2, 1]>
        """
        if i in ZZ:
            i = [i]
        return self.monomial(Partition(i))

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Fock space of rank {} over {}".format(self._n, self._q_group.base_ring())

    def basic_representation(self):
        """
        Return the basic representation `B(\Lambda_0)` in ``self``.
        """
        return BasicRepresentation(self._n, self._q, self.base_ring())

    class Element(CombinatorialFreeModuleElement):
        def e(self, i):
            """
            Apply the action of `e_i` on ``self``.
            """
            n = self.parent()._n
            N_left = lambda la, x, i: \
                len( filter(lambda y: y[1] < x[1], la.addable_cells_residue(i, n)) ) \
                - len( filter(lambda y: y[1] < x[1], la.removable_cells_residue(i, n)) )
            q = self.parent()._q
            return self.parent().sum_of_terms( [( la.remove_cell(*x), c * q**(-N_left(la, x, i)) )
                                for la,c in self for x in la.removable_cells_residue(i, n)] )

        def f(self, i):
            """
            Apply the action of `f_i` on ``self``.
            """
            n = self.parent()._n
            N_right = lambda la, x, i: \
                len( filter(lambda y: y[1] > x[1], la.addable_cells_residue(i, n)) ) \
                - len( filter(lambda y: y[1] > x[1], la.removable_cells_residue(i, n)) )
            q = self.parent()._q
            return self.parent().sum_of_terms( [(la.add_cell(*x), c * q**N_right(la, x, i))
                                for la,c in self for x in la.addable_cells_residue(i, n)] )

        def h(self, i):
            """
            Apply the action of `h_i` on ``self``.
            """
            N_i = lambda la, i: \
                len(la.addable_cells_residue(i, self.parent()._n)) \
                - len(la.removable_cells_residue(i, self.parent()._n))
            q = self.parent()._q
            return self.parent().sum_of_terms([(la, c * q**N_i(la, i) * self.parent().monomial(la)) for la,c in self])

        def d(self):
            """
            Apply the action of `d` on ``self``.
            """
            q = self.parent()._q
            return self.parent().sum_of_terms(
                [(la, c * q**len( filter(lambda x: la.content(*x) % self.parent()._n == 0, la.cells()) ))
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
    def __classcall_private__(cls, n, q=None, base_ring=None):
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
        return super(BasicRepresentation, cls).__classcall__(cls, n, q, base_ring)

    def __init__(self, n, q, base_ring):
        """
        Initialize ``self``.
        """
        self._n = n
        self._q = q
        self._fock = FockSpace(n, q, base_ring)
        Parent.__init__(self, base=self._fock._q_group.base_ring(),
                        category=CategoryOInt(self._fock._q_group).WithRealizations())

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Basic representation of {}".format(self._fock._q_group.cartan_type())

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
            self._fock = basic._fock
            self._n = basic._n
            self._q = basic._q
            CombinatorialFreeModule.__init__(self, self._q.parent(), Partitions(),
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
            cur = self._fock[[]]
            corners = la.corners()
            cells = set(la.cells())
            q = self._q
            k = self._n - 1 # This is sl_{k+1}
            b = ZZ.zero()
            while any(c[0]*k + c[1] >= b for c in corners): # While there is some cell left to count
                power = 0
                i = -b % (k+1)
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
            self._fock = basic._fock
            self._n = basic._n
            self._q = basic._q
            CombinatorialFreeModule.__init__(self, self._q.parent(), Partitions(),
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
                return self._fock.one()
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
            return "Lower global crystal basis of {}".format(self._fock)

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

