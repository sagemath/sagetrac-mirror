#*****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef fibers(f, domain)

from sage.structure.parent cimport Parent
from sage.structure.list_clone cimport ClonableIntArray

cdef class FiniteSetMap_MN(ClonableIntArray):
    cdef _setimage(self, int i, int j)
    cdef _getimage(self, int i)
    cdef setimage(self, i, j)
    cdef getimage(self, i)
    cdef domain(self)
    cdef codomain(self)
    cdef image_set(self)
    cdef fibers(self)
    cdef items(self)
    cdef FiniteSetMap_MN _compose_internal_(self, FiniteSetMap_MN other,
                                             Parent resParent)
    cdef check(self)

cdef class FiniteSetMap_Set(FiniteSetMap_MN): pass

cdef FiniteSetMap_Set FiniteSetMap_Set_from_list(cls, parent, lst)
cdef FiniteSetMap_Set FiniteSetMap_Set_from_dict(cls, parent, d)

cdef class FiniteSetEndoMap_N(FiniteSetMap_MN): pass
cdef class FiniteSetEndoMap_Set(FiniteSetMap_Set): pass
