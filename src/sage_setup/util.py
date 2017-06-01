r"""
Utility functions for building Sage
"""

#*****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def stable_uniq(L):
    """
    Given an iterable L, remove duplicate items from L by keeping only
    the last occurance of any item.

    The items must be hashable.

    EXAMPLES::

        sage: from sage_setup.util import stable_uniq
        sage: stable_uniq( (1, 2, 3, 4, 5, 6, 3, 7, 5, 1, 5, 9) )
        [2, 4, 6, 3, 7, 1, 5, 9]
    """
    D = {}
    for pos, item in enumerate(L):
        D[item] = pos  # Store the last position where an item appears
    return sorted(D, key=lambda item: D[item])
