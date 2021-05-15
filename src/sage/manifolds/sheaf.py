r"""
Presheaves and Sheaves over Manifolds

AUTHORS:

- Michael Jung (2021): initial version

"""

#******************************************************************************
#       Copyright (C) 2021 Michael Jung <m.jung at vu.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.abstract_method import abstract_method

class Presheaf(SageObject):
    r"""

    """
    def __init__(self, manifold, section_set_constructor):
        r"""

        EXAMPLES::

            sage: from sage.manifolds.sheaf import Presheaf
            sage: M = Manifold(2, 'M')
            sage: F = Presheaf(M, lambda d: d.scalar_field_algebra()); F
            ...
            sage: U = M.open_subset()
            sage: F(U)
            ...

        """
        if manifold is not self._manifold:
            raise ValueError
        self._manifold = manifold
        self._section_set_constructor = section_set_constructor

    def __call__(self, open_subset):
        r"""
        Return the set of sections of ``self`` w.r.t. ``open_subset``.

        """
        if not open_subset.is_subset(self._manifold):
            raise ValueError
        return self._section_set_constructor(open_subset)

    def section_set(self, open_subset):
        r"""
        Return the set of sections of ``self`` w.r.t. ``open_subset``.

        """
        return self(open_subset)

    def restriction_morphism(self, from_open_subset, to_open_subset):
        r"""
        Return the restriction morphism from ``from_open_subset`` to
        ``to_open_subset``.

        """
        if not from_open_subset.is_subset(self._manifold) or not to_open_subset.is_subset(self._manifold):
            raise ValueError
        from sage.categories.morphism import SetMorphism
        from sage.categories.homset import Hom
        from sage.categories.sets_cat import Sets
        sec_dom = self.section_set(from_open_subset)
        sec_codom = self.section_set(to_open_subset)
        return SetMorphism(Hom(sec_dom, sec_codom, Sets()),
                           lambda x: x.restrict(to_open_subset))

class PresheafSection(SageObject):
    r"""

    """
    def __init__(self, domain):
        r"""

        """
        self._domain = domain

    def restrict(self, open_subset):
        r"""
        Return the restriction of ``self`` to ``open_subset``.

        """
        # try to get restriction from restriction graph
        # ...
        # if that fails, create restriction from scratch
        return self._construct_restriction_(open_subset)

    @abstract_method
    def _construct_restriction_(self, open_subset):
        r"""
        Construct a restriction of ``self`` to ``open_subset`` from scratch.

        """