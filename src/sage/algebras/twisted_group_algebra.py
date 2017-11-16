"""
Twisted group algebras

Author: Mark Shimozono
"""
#*****************************************************************************
#  Copyright (C) 2015 Mark Shimozono <mshimo at math.vt.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.category import Category
from sage.categories.rings import Rings
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.sets_cat import Sets

class TwistedGroupAlgebraElement(CombinatorialFreeModule.Element):
    def acted_on(self, g):
        r"""
        Left action of group element `g` on ``self``.
        """
        par = self.parent()
        return par.sum_of_terms([(g*w,par.action()[g](c)) for (w,c) in self],distinct=True)

    def _mul_(self, other):
        par = self.parent()
        return par.sum([c*other.acted_on(w) for (w,c) in self])

class TwistedGroupAlgebra(CombinatorialFreeModule):
    r"""
    Twisted group algebra.

    INPUT:

        - `R` -- a not-necessarily commutative ring
        - `G` -- a group with multiplicative notation
        - `action` -- a family from  `G` to Aut(R).

    The action must define a group homomorphism from `G` to
    the automorphism group of `R`. Syntactically,
    `action(g)` need only be a function from `R` to `R`.

    As a left `R`-module the twisted group algebra is just the
    group algebra `R[G]`. However the product is given by

    ..MATH::

        (r_1 g_1)*(r_2 g_2) = (r_1 g_1(r_2), g_1 g_2)

    A ring map coercion is supplied from `R` to ``self``.
    For elements of `G` one may use :meth:`monomial`.

    ..TODO::

        Handle additive groups.
    """
    def __init__(self, R, G, action, prefix=None):
        """
        TESTS::

            sage: from sage.algebras.twisted_group_algebra import TwistedGroupAlgebra
            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: S = EquivariantCohomologyPoint(['A',2])
            sage: W = S.weyl_group()
            sage: action = Family(W, lambda w: S.weyl_automorphism(w))
            sage: A = TwistedGroupAlgebra(S, W, action)
            sage: TestSuite(A).run()
        """
        assert R in Rings()
        self._R = R
        self._G = G
        self._action = action
        if prefix is None:
            prefix = "X"
        CombinatorialFreeModule.__init__(self, R, G, prefix=prefix, category=ModulesWithBasis(R)&Rings())
        self._base_ring_morphism = SetMorphism(Hom(R, self, Rings()), lambda x: self.term(G.one(), x))
        self._base_ring_morphism.register_as_coercion()
        self._G_to_self = SetMorphism(Hom(self._G, self, Sets()), lambda w: self.monomial(w))
        self._G_to_self.register_as_coercion()

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.algebras.twisted_group_algebra import TwistedGroupAlgebra
            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: S = EquivariantCohomologyPoint(['A',2])
            sage: W = S.weyl_group()
            sage: action = Family(W, lambda w: S.weyl_automorphism(w))
            sage: TwistedGroupAlgebra(S, W, action, "W")
            Twisted group ring of Equivariant cohomology of a point for the maximal torus of type ['A', 2] and Weyl Group of type ['A', 2] (as a matrix group acting on the weight space)
        """
        return "Twisted group ring of %s and %s"%(self._R, self._G)

    def action(self):
        r"""
        The function for the action of the group on the coefficient ring.

        EXAMPLES::

            sage: from sage.algebras.twisted_group_algebra import TwistedGroupAlgebra
            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: S = EquivariantCohomologyPoint(['A',2],'ambient')
            sage: W = S.weyl_group()
            sage: action = Family(W, lambda w: S.weyl_automorphism(w))
            sage: A = TwistedGroupAlgebra(S, W, action, "W")
            sage: w = A.group().an_element(); w
            s1*s2
            sage: r = A.coefficient_ring().an_element(); r
            x0
            sage: action[w](r*r+1)
            x1^2 + 1
        """
        return self._action

    def coefficient_ring(self):
        r"""
        The coefficient ring of ``self``.

        EXAMPLES::

            sage: from sage.algebras.twisted_group_algebra import TwistedGroupAlgebra
            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: S = EquivariantCohomologyPoint(['A',2],'ambient')
            sage: W = S.weyl_group()
            sage: action = Family(W, lambda w: S.weyl_automorphism(w))
            sage: A = TwistedGroupAlgebra(S, W, action, "W")
            sage: A.coefficient_ring()
            Equivariant cohomology of a point for the maximal torus of type ['A', 2]
        """
        return self._R

    def group(self):
        r"""
        The group of ``self``.

        EXAMPLES::

            sage: from sage.algebras.twisted_group_algebra import TwistedGroupAlgebra
            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: S = EquivariantCohomologyPoint(['A',2],'ambient')
            sage: W = S.weyl_group()
            sage: action = Family(W, lambda w: S.weyl_automorphism(w))
            sage: A = TwistedGroupAlgebra(S, W, action, "W")
            sage: A.group()
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
        """
        return self._G

    @cached_method
    def one_basis(self):
        r"""
        The index of identity basis element.

        EXAMPLES::

            sage: from sage.algebras.twisted_group_algebra import TwistedGroupAlgebra
            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: S = EquivariantCohomologyPoint(['A',2],'ambient')
            sage: W = S.weyl_group()
            sage: action = Family(W, lambda w: S.weyl_automorphism(w))
            sage: TwistedGroupAlgebra(S, W, action, "W").one_basis()
            1
        """
        return self.group().one()

    def an_element(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.algebras.twisted_group_algebra import TwistedGroupAlgebra
            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: S = EquivariantCohomologyPoint(['A',2],'ambient')
            sage: W = S.weyl_group()
            sage: action = Family(W, lambda w: S.weyl_automorphism(w))
            sage: TwistedGroupAlgebra(S, W, action, "W").an_element()
            W[s1*s2] + x0*W[s1] + W[s2] + W[1]
        """
        group = self.group()
        unit_coef = self.coefficient_ring().one()
        group_elements = [x for x in group.some_elements()]
        coef = self.coefficient_ring().an_element()
        terms = [(group_elements[0],coef)]+[(x, unit_coef) for x in group_elements[1:]]
        return self.sum_of_terms(terms)

    Element = TwistedGroupAlgebraElement
