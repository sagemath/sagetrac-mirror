r"""
Convex Sets
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.abstract_method import abstract_method

class ConvexSet_base(SageObject):

    """
    Abstract base class for convex sets.
    """

    @abstract_method
    def is_empty(self):
        r"""
        Test whether ``self`` is the empty set

        OUTPUT:

        Boolean.
        """

    @abstract_method
    def is_universe(self):
        r"""
        Test whether ``self`` is the whole ambient space

        OUTPUT:

        Boolean.
        """

    @abstract_method
    def is_full_dimensional(self):
        r"""
        Return whether ``self`` is full dimensional.

        OUTPUT:

        Boolean. Whether the polyhedron is not contained in any strict
        affine subspace.

        """

    @abstract_method
    def is_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        """

    def is_relatively_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        """
        if self.is_open():
            return True
        raise NotImplementedError

    @abstract_method
    def is_closed(self):
        r"""
        Return whether ``self`` is closed.

        OUTPUT:

        Boolean.

        """

    def closure(self):
        r"""
        Return the topological closure of ``self``.
        """
        if self.is_closed():
            return self
        raise NotImplementedError

    def interior(self):
        r"""
        Return the topological interior of ``self``.
        """
        if self.is_closed():
            return self
        raise NotImplementedError

    @abstract_method(optional=True)
    def affine_hull(self):
        r"""
        Return the affine hull of ``self``.
        """


class ConvexSet_closed(ConvexSet_base):

    r"""
    Abstract base class for closed convex sets.
    """

    def is_closed(self):
        r"""
        Return whether ``self`` is closed.

        OUTPUT:

        Boolean.
        """
        return True

    def is_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        """
        return self.is_empty() or self.is_universe()


class ConvexSet_relatively_open(ConvexSet_base):

    r"""
    Abstract base class for relatively open sets.
    """

    def is_relatively_open(self):
        r"""
        Return whether ``self`` is open.

        OUTPUT:

        Boolean.

        """
        return True
