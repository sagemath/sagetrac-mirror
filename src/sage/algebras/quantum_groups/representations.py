r"""
Representations of Quantum Groups

There are the following representations for quantum groups:

- `U_q(\mathfrak{g})`-representations from highest weight crystals,

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

#from sage.structure.element import parent
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass

from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.finite_dimensional_modules_with_basis import FiniteDimensionalModulesWithBasis
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.realizations import Realizations, Category_realization_of_parent

from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.combinat.crystals.tensor_product import TensorProductOfCrystals
from sage.algebras.quantum_groups.q_numbers import q_number, q_factorial

#####################################################################
## Base classes

class BasisBaseClass(CombinatorialFreeModule):
    r"""
    Base class for bases of a highest weight
    `U_q(\mathfrak{g})`-representation from a crystal.
    """
    def __init__(self, V, prefix, name):
        """
        Initialize ``self``.
        """
        self._basis_name = name
        cat = HighestWeightRepresentationBasesCategory(V)
        CombinatorialFreeModule.__init__(self, V.base_ring(), V._crystal,
                                         category=cat,
                                         prefix=prefix, bracket=False)

    @cached_method
    def highest_weight_vector(self, b, index_set=None):
        r"""
        Return the highest weight vector corresponding to the highest weight
        crystal element ``b``.
        """
        if index_set is None:
            index_set = self.basis().keys().index_set()
        if not b.is_highest_weight(index_set):
            raise ValueError("{} is not a {} highest weight element".format(b, index_set))
        cur = self.monomial(b)
        index_set = self.basis().keys().index_set()
        for i in index_set:
            next = cur.e(i)
            while next:
                m,c = next.leading_item()
                elt = self.monomial(m).f(i)
                print "cur:", cur, '\n'
                print "next:", next, '\n'
                print "elt:", elt
                print "\n----------\n"
                S = elt.support()
                for s in cur.support():
                    if s in S:
                        elt = self.sum(elt.terms()[:S.index(s)])
                        break
                print "elt red:", elt, '\n'
                m = elt.leading_monomial()
                print "leading monomial:", m
                print "\n=====================\n"
                cur -= c / m.e(i).leading_coefficient() * m
                next = cur.e(i)
        return cur

#####################################################################
## Natural basis

class NaturalBasisFromCrystal(BasisBaseClass):
    r"""
    A highest weight `U_q(\mathfrak{g})`-representation from a crystal in
    the natrual basis.

    EXAMPLES::

        sage: C = CrystalOfTableaux(['A',1], shape=[2])
        sage: V = HighestWeightRepresentationFromCrystal(C)
        sage: N = V.natural()

    We check the coercion going to the divided power basis::

        sage: D = V.divided_power()
        sage: all(N(D(N[c])) == N[c] for c in C)
        True
    """
    def __init__(self, V, prefix='v'):
        r"""
        Initialize ``self``.

        TESTS::

            sage: C = CrystalOfTableaux(['A',1], shape=[2])
            sage: V = HighestWeightRepresentationFromCrystal(C)
            sage: TestSuite(V.natural()).run()
        """
        BasisBaseClass.__init__(self, V, prefix, "natural")
        phi = self.module_morphism(diagonal=self._divided_power_scalar,
                                   codomain=V.divided_power())
        phi.register_as_coercion()
        (~phi).register_as_coercion()

    def _divided_power_scalar(self, b):
        r"""
        Return the scalar to send ``b`` to the divided power basis.

        EXAMPLES::

            sage: C = CrystalOfTableaux(['A',1], shape=[2])
            sage: V = HighestWeightRepresentationFromCrystal(C)
            sage: N = V.natural()
            sage: map(N._divided_power_scalar, C)
        """
        q = self.realization_of()._q
        coeff = self.base_ring().one()
        path = b.to_highest_weight()[1]

        if not path:
            return coeff

        last = path[-1]
        p = 0
        for i in path:
            if i != last:
                coeff *= q_factorial(p, q) #q_{last}?
                p = 1
                last = i
            else:
                p += 1
        return coeff * q_factorial(p, q) #q_{last}?

    def e_on_basis(self, i, b):
        r"""
        Return the action of `e_i` on the basis element indexed by ``b``.
        """
        next = b.e(i)
        if next is None:
            return self.zero()

        q = self.realization_of()._q # q_i?
        h = self.basis().keys().weight_lattice_realization().simple_coroots()
        b = next
        path = b.to_highest_weight()[1]
        p = h[i].scalar(b.weight())
        coeff = (q**p - q**-p) / (q - q**-1)
        for j in path:
            next = next.e(j)
            if i == j:
                p = h[i].scalar(next.weight())
                coeff += (q**p - q**-p) / (q - q**-1)
        return self.term(b, coeff)

    def f_on_basis(self, i, b):
        r"""
        Return the action of `f_i` on the basis element indexed by ``b``.
        """
        next = b.f(i)
        if next is None:
            return self.zero()
        return self.monomial(next)

    class Element(CombinatorialFreeModuleElement):
        r"""
        An element in the divided power basis.
        """
        def e_kashiwara(self, i):
            r"""
            Apply the Kashiwara operator `\tilde{e}_i` to ``self``.
            """
            P = self.parent()
            q = P.realization_of()._q # q_i?
            terms = {}
            for b,c in self:
                next = b.e(i)
                if next is None:
                    continue
                terms[next] = c * q_number(b.epsilon(i), q)
            return P._from_dict(terms)

        def f_kashiwara(self, i):
            r"""
            Apply the Kashiwara operator `\tilde{f}_i` to ``self``.
            """
            P = self.parent()
            q = P.realization_of()._q # q_i?
            terms = {}
            for b,c in self:
                b = b.f(i)
                if b is None:
                    continue
                terms[b] = c / q_number(b.epsilon(i), q)
            return P._from_dict(terms)

class NaturalBasisFromTensorProductOfCrystals(NaturalBasisFromCrystal):
    r"""
    Natural basis for a tensor product of `U_q(\mathfrak{g})`-crystals.
    """
    def e_on_basis(self, i, b):
        r"""
        Return the action of `e_i` on the basis element indexed by ``b``.
        """
        TP = self.basis().keys()
        h = TP.weight_lattice_realization().simple_coroots()

        ret = []
        k = len(TP.crystals)
        coeff = self.base_ring().one()
        q = self.realization_of()._q # q_i?
        for n in range(k): # We are using the anti-Kashiwara convention
            next = NaturalBasisFromCrystal.e_on_basis(self, i, b[n])
            for m,c in next:
                ret.append((TP(*(b[:n] + [m] + b[n+1:])), c * coeff))
            coeff *= q**-b[n].weight().scalar(h[i])
        return self.sum_of_terms(ret)

    def f_on_basis(self, i, b):
        r"""
        Return the action of `f_i` on the basis element indexed by ``b``.
        """
        TP = self.basis().keys()
        h = TP.weight_lattice_realization().simple_coroots()

        ret = []
        k = len(TP.crystals)
        coeff = self.base_ring().one()
        q = self.realization_of()._q # q_i?
        for n in reversed(range(k)): # We are using the anti-Kashiwara convention
            next = NaturalBasisFromCrystal.f_on_basis(self, i, b[n])
            for m,c in next:
                ret.append((TP(*(b[:n] + [m] + b[n+1:])), c * coeff))
            coeff *= q**b[n].weight().scalar(h[i])
        return self.sum_of_terms(ret)

    Element = CombinatorialFreeModuleElement

#####################################################################
## Divided power basis

class DividedPowerBasisFromCrystal(BasisBaseClass):
    r"""
    A `U_q(\mathfrak{g})`-representation given by a crystal in the
    divided power basis.

    Recall that a divided power operator is given by:

    .. MATH::

        X_i^{(r)} = \frac{X_i^r}{[r]_{d_i}!}

    where

    .. MATH::

        [r]_d = \frac{q^{dr} - v^{-rd}}{v^d - v^{-d}}.


    EXAMPLES::

        sage: C = CrystalOfTableaux(['A',1], shape=[2])
        sage: V = HighestWeightRepresentationFromCrystal(C)
        sage: D = V.divided_power()

    We check the coercion going to the natural basis::

        sage: N = V.natural()
        sage: all(D(N(D[c])) == D[c] for c in C)
        True
    """
    def __init__(self, V, prefix='d'):
        r"""
        Initialize ``self``.

        TESTS::

            sage: C = CrystalOfTableaux(['A',1], shape=[2])
            sage: V = HighestWeightRepresentationFromCrystal(C)
            sage: TestSuite(V.divided_power()).run()
        """
        BasisBaseClass.__init__(self, V, prefix, "divided power")

    def e_on_basis(self, i, b):
        r"""
        Return the action of `e_i` on the basis element indexed by ``b``.
        """
        next = b.e(i)
        if next is None:
            return self.zero()

        q = self.realization_of()._q # q_i?
        h = self.basis().keys().weight_lattice_realization().simple_coroots()
        b = next
        path = b.to_highest_weight()[1]
        p = h[i].scalar(b.weight())
        coeff = (q**p - q**-p) / (q - q**-1)
        for j in path:
            next = next.e(j)
            if i == j:
                p = h[i].scalar(next.weight())
                coeff += (q**p - q**-p) / (q - q**-1)
        return self.term(b, coeff / q_number(b.epsilon(i)+1, q))

    def f_on_basis(self, i, b):
        r"""
        Return the action of `f_i` on the basis element indexed by ``b``.
        """
        next = b.f(i)
        if next is None:
            return self.zero()

        q = self.realization_of()._q # q_i?
        return self.term( next, q_number(next.epsilon(i), q) )

    class Element(CombinatorialFreeModuleElement):
        r"""
        An element in the divided power basis.
        """
        def e_kashiwara(self, i):
            r"""
            Apply the Kashiwara operator `\tilde{e}_i` to ``self``.
            """
            P = self.parent()
            terms = {}
            for b,c in self:
                next = b.e(i)
                if next is None:
                    continue
                terms[next] = c
            return P._from_dict(terms)

        def f_kashiwara(self, i):
            r"""
            Apply the Kashiwara operator `\tilde{f}_i` to ``self``.
            """
            P = self.parent()
            terms = {}
            for b,c in self:
                next = b.f(i)
                if next is None:
                    continue
                terms[next] = c
            return P._from_dict(terms)

class DividedPowerBasisFromTensorProductOfCrystals(DividedPowerBasisFromCrystal):
    r"""
    Divided power basis for a tensor product of `U_q(\mathfrak{g})`-crystals.
    """
    def e_on_basis(self, i, b):
        r"""
        Return the action of `e_i` on the basis element indexed by ``b``.
        """
        TP = self.basis().keys()
        h = TP.weight_lattice_realization().simple_coroots()

        ret = []
        k = len(TP.crystals)
        coeff = self.base_ring().one()
        q = self.realization_of()._q # q_i?
        for n in range(k): # We are using the anti-Kashiwara convention
            next = DividedPowerBasisFromCrystal.e_on_basis(self, i, b[n])
            for m,c in next:
                ret.append((TP(*(b[:n] + [m] + b[n+1:])), c * coeff))
            coeff *= q**-b[n].weight().scalar(h[i])
        return self.sum_of_terms(ret)

    def f_on_basis(self, i, b):
        r"""
        Return the action of `f_i` on the basis element indexed by ``b``.
        """
        TP = self.basis().keys()
        h = TP.weight_lattice_realization().simple_coroots()

        ret = []
        k = len(TP.crystals)
        coeff = self.base_ring().one()
        q = self.realization_of()._q # q_i?
        for n in reversed(range(k)): # We are using the anti-Kashiwara convention
            next = DividedPowerBasisFromCrystal.f_on_basis(self, i, b[n])
            for m,c in next:
                ret.append((TP(*(b[:n] + [m] + b[n+1:])), c * coeff))
            coeff *= q**b[n].weight().scalar(h[i])
        return self.sum_of_terms(ret)

    Element = CombinatorialFreeModuleElement

#####################################################################
## Highest weight representation class

class HighestWeightRepresentationFromCrystal(Parent, UniqueRepresentation):
    r"""
    A highest weight `U_q(\mathfrak{g})`-representation given from a highest
    weight crystal.

    INPUT:

    - ``crystal`` -- a crystal
    - ``q`` -- (optional) the variable `q`; default is the `q` of `QQ(q)`
    - ``base`` -- (optional) the base ring; default is the parent of ``q``

    EXAMPLES::

        sage: C = CrystalOfTableaux(['A',1], shape=[2])
        sage: V = HighestWeightRepresentationFromCrystal(C)
    """
    @staticmethod
    def __classcall_private__(cls, crystal, q=None, base=None):
        """
        Normalize input to ensure a unique representation.
        """
        if crystal not in HighestWeightCrystals():
            raise ValueError("only for highest weight crystals")

        if base is None:
            if q is None:
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                from sage.rings.all import QQ
                q = PolynomialRing(QQ, 'q').gen(0)
            base = q.parent()
        else:
            q = base('q')
        # Make sure it's in the fraction field
        base = base.fraction_field()
        q = base(q)

        return super(HighestWeightRepresentationFromCrystal, cls).__classcall__(cls, crystal, q)

    def __init__(self, crystal, q):
        """
        Initialize ``self``.
        """
        R = q.parent()
        self._q = q
        self._crystal = crystal

        if crystal in FiniteCrystals():
            self._category = FiniteDimensionalModulesWithBasis(R)
        else:
            self._category = ModulesWithBasis(R)
        Parent.__init__(self, base=R, category=self._category.WithRealizations())

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Representation of {} over {}".format(self._crystal, self.base_ring())

    def natural(self, prefix='v'):
        """
        The natural basis.
        """
        if isinstance(self._crystal, TensorProductOfCrystals):
            return NaturalBasisFromTensorProductOfCrystals(self, prefix)
        return NaturalBasisFromCrystal(self, prefix)

    def divided_power(self, prefix='d'):
        """
        The divided power basis.
        """
        if isinstance(self._crystal, TensorProductOfCrystals):
            return DividedPowerBasisFromTensorProductOfCrystals(self, prefix)
        return DividedPowerBasisFromCrystal(self, prefix)

    class canonical(BasisBaseClass, BindableClass):
        """
        The canonical (or global) basis of a highest weight module
        in `\mathcal{O}^q_{int}`.
        """
        def __init__(self, V, prefix='G'):
            """
            Initialize ``self``.
            """
            BasisBaseClass.__init__(self, V, prefix, "canonical")

            # Coercions to/from the natural basis
            phi = self.module_morphism(self._to_natural_basis,
                                       codomain=V.natural(), triangular='lower')
            phi.register_as_coercion()
            (~phi).register_as_coercion()
            # Coercions to/from the divided power basis
            phi = self.module_morphism(self._to_divided_power_basis,
                                       codomain=V.divided_power(), triangular='lower')
            phi.register_as_coercion()
            (~phi).register_as_coercion()

        def _latex_(self):
            r"""
            Return a `\LaTeX` representation of ``self``.
            """
            from sage.misc.latex import latex
            return "G\\left( {} \\right)".format(latex(self._V)) # FIXME

        @cached_method
        def highest_weight_vector(self, b):
            """
            Return the highest weight vector corresponding to the highest weight
            crystal element ``b``.
            """
            if not b.is_highest_weight():
                raise ValueError("{} is not a highest weight element".format(b))
            return self.monomial(b)

        @cached_method
        def _to_natural_basis(self, m):
            """
            Return the basis element indexed by ``m`` in the natural basis
            of the tensor product.

            EXAMPLES:

            We construct the examples in [HK02]_ Section 6.2::

                from sage.algebras.quantum_groups.representations import GlobalCrystalBasis
                sage: Q = QuantumGroup(['A',1])
                sage: B = CrystalOfTableaux(['A',1],shape=[2])
                sage: T = TensorProductOfCrystals(B,B)
                sage: G = GlobalCrystalBasis(Q, T)
                sage: map(lambda x: G.to_natural_basis(x), T)
            """
            B = self.basis().keys()
            N = self.realization_of().natural()
            # Special case for a single factor
            if not isinstance(B, TensorProductOfCrystals):
                return N.term(m, N._divided_power_scalar(m))

            # Construct the maximal vector
            hw, path = m.to_highest_weight()
            cur = N.highest_weight_vector(hw)

            # Now build down the crystal until we reach the desired node
            if not path:
                return cur

            last = path[0]
            p = 0
            coeff = self.base_ring().one()
            q = self.realization_of()._q
            for i in reversed(path):
                if i != last:
                    coeff *= q_factorial(p, q) #q_{last}?
                    p = 1
                    last = i
                else:
                    p += 1
                cur = cur.f(i)
            return cur / (coeff * q_factorial(p, q)) #q_{last}?

        @cached_method
        def _to_divided_power_basis(self, m):
            """
            Return the basis element indexed by ``m`` in the natural basis
            of the tensor product.

            EXAMPLES:

            We construct the examples in [HK02]_ Section 6.2::

                sage: C = CrystalOfTableaux(['A',1], shape=[2])
                sage: T = TensorProductOfCrystals(C, C)
                sage: V = HighestWeightRepresentationFromCrystal(T)
                sage: G = V.canonical()
                sage: map(lambda x: G._to_divided_power_basis(x), T)
            """
            B = self.basis().keys()
            D = self.realization_of().divided_power()
            # Special case for a single factor
            if not isinstance(B, TensorProductOfCrystals):
                return D.monomial(m)

            # Construct the maximal vector
            hw, path = m.to_highest_weight()
            cur = D.highest_weight_vector(hw)

            # Now build down the crystal until we reach the desired node
            if not path:
                return cur

            for i in reversed(path):
                cur = cur.f(i)
            return cur

        def e_on_basis(self, i, b):
            r"""
            Return the action of `e_i` on the basis element indexed by ``b``.
            """
            return self(self._to_natural_basis(b).e(i))

        def f_on_basis(self, i, b):
            r"""
            Return the action of `f_i` on the basis element indexed by ``b``.
            """
            return self(self._to_natural_basis(b).f(i))

#####################################################################
## The category for the bases

class HighestWeightRepresentationBasesCategory(Category_realization_of_parent):
    r"""
    The category of bases of a highest weight
    `U_q(\mathfrak{g})`-representation.
    """
    def __init__(self, base):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base`` -- a highest weight `U_q(\mathfrak{g})`-representation

        TESTS::
        """
        Category_realization_of_parent.__init__(self, base)

    def super_categories(self):
        r"""
        The super categories of ``self``.

        EXAMPLES::

            sage: bases = H._BasesCategory()
            sage: bases.super_categories()
            [Category of realizations of Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring,
             Category of finite dimensional algebras with basis over Integer Ring]
        """
        return [Realizations(self.base()), self.base()._category]

    def _repr_(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: H._BasesCategory()
            Category of bases of Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring
        """
        return "Category of bases of {}".format(self.base())

    class ParentMethods:
        r"""
        This class collects code common to all the various bases. In most
        cases, these are just default implementations that will get
        specialized in a basis.
        """
        def _repr_(self):
            """
            Text representation of this basis of Iwahori-Hecke algebra.

            EXAMPLES::

                sage: H.T()
                Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring in the T-basis
                sage: H.C()
                Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring in the C-basis
                sage: H.Cp()
                Iwahori-Hecke algebra of type B2 in 1,-1 over Integer Ring in the Cp-basis
            """
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

        def __getitem__(self, i):
            """
            Return the basis element indexed by ``i``.

            INPUT:

            - ``i`` -- the index

            OUTPUT:

            If ``i`` is in the crystal ``B``, then return ``i`` the basis
            element indexed by ``i``; otherwise return ``B[i]``.

            EXAMPLES::
            """
            B = self.realization_of()._crystal
            if i in B:
                return self.monomial(i)
            return self.monomial(self._crystal[i])

    class ElementMethods:
        def e(self, i):
            """
            Return the action of `e_i` on ``self``.
            """
            P = self.parent()
            return P.linear_combination( (P.e_on_basis(i, b), c) for b,c in self )

        def f(self, i):
            """
            Return the action of `f_i` on ``self``.
            """
            P = self.parent()
            return P.linear_combination( (P.f_on_basis(i, b), c) for b,c in self )

        def h(self, i):
            """
            Return the action of `h_i` on ``self``.
            """
            P = self.parent()
            q = P.realization_of()._q
            h = P.basis().keys().weight_lattice_realization().simple_coroots()
            return P.sum_of_terms([(m, c * q**h[i].scalar(m.weight())) for m,c in self])

        def e_string(self, alpha):
            r"""
            Apply `e_{i_r} \cdots e_{i_1}` to ``self`` where ``alpha``
            is `(i_1, \ldots, i_r)`
            """
            cur = self
            for i in alpha:
                cur = cur.e(i)
                if not cur:
                    break
            return cur

        def f_string(self, alpha):
            r"""
            Apply `f_{i_r} \cdots f_{i_1}` to ``self`` where ``alpha``
            is `(i_1, \ldots, i_r)`
            """
            cur = self
            for i in alpha:
                cur = cur.f(i)
                if not cur:
                    break
            return cur

        def e_divided_power(self, i, k):
            r"""
            Apply the divided power operator `e_i^{(k)}` to ``self``.

            The divided power operator is defined as

            .. MATH::

                e_i^{(k)} := \frac{e_i^k}{[k]_{q_i}!}.
            """
            cur = self
            for dummy in range(k):
                cur = cur.e(i)
                if not cur:
                    break
            q = self.parent().realization_of()._q # q_i?
            return cur / q_factorial(k, q)

        def f_divided_power(self, i, k):
            r"""
            Apply the divided power operator `f_i^{(k)}` to ``self``.

            The divided power operator is defined as

            .. MATH::

                f_i^{(k)} := \frac{f_i^k}{[k]_{q_i}!}.
            """
            cur = self
            for dummy in range(k):
                cur = cur.f(i)
                if not cur:
                    break
            q = self.parent().realization_of()._q # q_i?
            return cur / q_factorial(k, q)

