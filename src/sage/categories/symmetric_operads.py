r"""
Symmetric Operads
"""
#*****************************************************************************
#  Copyright (C) 2011 Floent Hivert (CNRS) <Florent.Hivert@lri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import Operads, SymmetricSetOperads
from sage.categories.cartesian_product import cartesian_product
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute

class SymmetricOperads(Category_over_base_ring):
    """
    The category of symmetric operads

    EXAMPLES::

      sage: SymmetricOperads(ZZ)
      Category of symmetric operads over Integer Ring
      sage: SymmetricOperads(ZZ).super_categories()
      [Category of operads over Integer Ring, Category of symmetric set operads]

    TESTS::

        sage: C = SymmetricOperads(ZZ)
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: SymmetricOperads(QQ).super_categories()
            [Category of operads over Rational Field, Category of symmetric set operads]
        """
        R = self.base_ring()
        return [Operads(R), SymmetricSetOperads()]

    class ParentMethods:

        # Should be in symmetric_operad_with_basis
        @abstract_method(optional = True)
        def symmetric_group_action_on_basis(self, basis_elem, perm):
            """
            returns the action of ``perm`` over ``elem``
            """

        @lazy_attribute
        def symmetric_group_action(self):
            """
            returns the action of ``perm`` over ``elem``
            """
            if self.composition_on_basis is not NotImplemented:
                return self_module_morphism(
                    self.symmetric_group_action_on_basis,
                    position = 0, codomain = self)
            else:
                return NotImplemented


