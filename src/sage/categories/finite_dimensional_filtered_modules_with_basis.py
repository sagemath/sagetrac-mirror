r"""
Finite Dimensional Filtered Modules With Basis
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.filtered_modules import FilteredModulesCategory

class FiniteDimensionalFilteredModulesWithBasis(FilteredModulesCategory):
    class ParentMethods:
        def homogeneous_component_basis(self, d):
            """
            Return a basis for the ``d``-th graded component of ``self``.

            EXAMPLES::

                sage: cat = ModulesWithBasis(ZZ).Filtered().FiniteDimensional()
                sage: C = CombinatorialFreeModule(ZZ, ['a', 'b'], category=cat)
                sage: C.degree_on_basis = lambda x: 1 if x == 'a' else 2
                sage: C.homogeneous_component_basis(1)
                Finite family {'a': B['a']}
                sage: C.homogeneous_component_basis(2)
                Finite family {'b': B['b']}
            """
            from sage.sets.family import Family
            try:
                S = self._indices.subset(size=d)
            except (AttributeError, ValueError, TypeError):
                S = [i for i in self._indices if self.degree_on_basis(i) == d]
            return Family(S, self.monomial)

