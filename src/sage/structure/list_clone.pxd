#*****************************************************************************
#  Copyright (C) 2009-2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element cimport Element


# Cython-0.17.2 disallows inline cdef in non-final classes
# This restriction will be lifted at one point, then we can set
# some of the methods to be inline again, that is,
# revert the patch form http://trac.sagemath.org/13740

cdef class ClonableElement(Element):
    cdef bint _is_immutable
    cdef bint _needs_check
    cdef long int  _hash

    cdef bint _require_mutable(self) except -2
    cdef bint is_mutable(self)
    cdef bint is_immutable(self)
    cdef set_immutable(self)

    cdef _set_mutable(self)

    cdef ClonableElement clone(self, bint check=?)

cdef class ClonableArray(ClonableElement):
    cdef list _list

    cdef list _get_list(self)
    cdef _set_list(self, list lst)
    cdef ClonableArray __copy__(self)
    cdef check(self)
    cdef object _getitem(self, int key)
    cdef _setitem(self, int key, value)
    cdef int index(self, key, start=*, stop=*) except -1
    cdef int count(self, key) except -1
    cdef long int _hash_(self) except? -1

cdef class ClonableList(ClonableArray):
    cdef append(self, el)
    cdef extend(self, it)
    cdef insert(self, int index, el)
    cdef pop(self, int index=*)
    cdef remove(self, el)

cdef class NormalizedClonableList(ClonableList):
    cdef normalize(self)

cdef class ClonableIntArray(ClonableElement):
    cdef int _len
    cdef int* _list

    cdef _alloc_(self, int size)
    cdef ClonableIntArray __copy__(self)
    cdef check(self)
    cdef object _getitem(self, int key)
    cdef _setitem(self, int item, value)
    cdef int index(self, int item) except -1
    cdef long int _hash_(self) except? -1
    cdef list list(self)
