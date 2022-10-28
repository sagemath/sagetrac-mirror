"""
Mutability -- Pyrex Implementation
"""

##########################################################################
#
#   Sage: Open Source Mathematical Software
#
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
##########################################################################

cdef class Mutability:
    cdef public bint _is_immutable
    cdef _require_mutable(self)
    cdef _require_immutable(self)
    cdef bint is_immutable(self)
    cdef bint is_mutable(self)