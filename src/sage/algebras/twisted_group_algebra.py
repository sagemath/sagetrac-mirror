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

class TwistedGroupAlgebra(CombinatorialFreeModule):
    r"""
    Twisted group algebra.

    INPUT:

        - `R` -- a not-necessarily commutative ring
        - `G` -- a group with multiplicative notation
        - `action` -- a function from `G \times R` to `R`.

    The action must define a group homomorphism from `G` to
    the automorphism group of `R`.

    As a left `R`-module the twisted group algebra is just the
    group algebra `R[G]`. However the product is given by

    ..MATH::

        (r_1 g_1)*(r_2 g_2) = (r_1 g_1(r_2), g_1 g_2)

    ..TODO::

        Handle additive groups.

    ..WARNING::

    We must avoid writing expressions of the form `g*r` since sage changes
    it to `r*g` rather than the correct `g(r)*g`.
    """

    def __init__(self, R, G, action):
        """
        TESTS::

            sage: from sage.categories.examples.twisted_group_algebra import MyTwistedGroupAlgebra
            sage: A = MyTwistedGroupAlgebra("A2")
            sage: TestSuite(A).run()
        """
        assert R in Rings()
        self._R = R
        self._G = G
        self._action = action
        CombinatorialFreeModule.__init__(self, R, G, prefix="X",category=Category.join([ModulesWithBasis(R),Rings()]))
        self._base_ring_morphism = SetMorphism(Hom(R, self, Rings()), lambda x: x * self.one())
        self._base_ring_morphism.register_as_coercion()

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.examples.twisted_group_algebra import MyTwistedGroupAlgebra
            sage: MyTwistedGroupAlgebra("A2")
            Twisted group ring of Fraction Field of Group algebra of the Ambient space of the Root system of type ['A', 2] over Integer Ring and Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
        """
        return "Twisted group ring of %s and %s"%(self._R, self._G)

    def action(self):
        r"""
        The function for the action of the group on the coefficient ring.

        EXAMPLES::

            sage: from sage.categories.examples.twisted_group_algebra import MyTwistedGroupAlgebra
            sage: A = MyTwistedGroupAlgebra("A2")
            sage: action = A.action()
            sage: w = A.group().an_element(); w
            s1*s2
            sage: r = A.coefficient_ring().an_element(); r
            x0
            sage: action(w,1/r)
            1/x1
        """
        return self._action

    def coefficient_ring(self):
        r"""
        The coefficient ring of ``self``.

        EXAMPLES::

            sage: from sage.categories.examples.twisted_group_algebra import MyTwistedGroupAlgebra
            sage: MyTwistedGroupAlgebra("A2").coefficient_ring()
            Fraction Field of Group algebra of the Ambient space of the Root system of type ['A', 2] over Integer Ring
        """
        return self._R

    def group(self):
        r"""
        The group of ``self``.

        EXAMPLES::

            sage: from sage.categories.examples.twisted_group_algebra import MyTwistedGroupAlgebra
            sage: MyTwistedGroupAlgebra("A2").group()
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
        """
        return self._G

    @cached_method
    def one(self):
        r"""
        The unit of the twisted group algebra.

        EXAMPLES::

            sage: from sage.categories.examples.twisted_group_algebra import MyTwistedGroupAlgebra
            sage: MyTwistedGroupAlgebra("A2").one()
            X[1]
        """
        return self.monomial(self._G.one())

    def mul_terms(self, (g,r),(h,s)):
        r"""
        Multiply `(rg)*(sh)`.

        Answer is `r g(s) (gh)`.

        EXAMPLES::

            sage: from sage.categories.examples.twisted_group_algebra import MyTwistedGroupAlgebra
            sage: A = MyTwistedGroupAlgebra("A2")
            sage: r = A.coefficient_ring().an_element(); r
            x0
            sage: g = A.group().an_element(); g
            s1*s2
            sage: A.mul_terms((g,r),(g,r))
            x0*x1*X[s2*s1]
        """
        return self.term(g*h,r*self.action()(g,s))

    def product(self, a, b):
        r"""
        The product in ``self``.

        EXAMPLES::

            sage: from sage.categories.examples.twisted_group_algebra import MyTwistedGroupAlgebra
            sage: A = MyTwistedGroupAlgebra("A2")
            sage: a = A.an_element(); a
            X[s1*s2] + x0*X[s1] + X[s2] + X[1]
            sage: g = A.group().an_element(); g
            s1*s2
            sage: a * A.monomial(g)
            X[s1*s2*s1] + X[s1*s2] + X[s2*s1] + x0*X[s2]
            sage: A.monomial(g) * a
            x1*X[s1*s2*s1] + X[s1*s2] + X[s2*s1] + X[s1]
        """
        return self.sum([self.mul_terms(t1,t2) for t1 in a for t2 in b])

    def an_element(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.categories.examples.twisted_group_algebra import MyTwistedGroupAlgebra
            sage: MyTwistedGroupAlgebra("A2").an_element()
            X[s1*s2] + x0*X[s1] + X[s2] + X[1]
        """
        group = self.group()
        unit_coef = self.coefficient_ring().one()
        group_elements = [x for x in group.some_elements()]
        coef = self.coefficient_ring().an_element()
        terms = [(group_elements[0],coef)]+[(x, unit_coef) for x in group_elements[1:]]
        return self.sum_of_terms(terms)
