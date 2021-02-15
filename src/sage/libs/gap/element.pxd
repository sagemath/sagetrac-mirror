#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from gappy.core cimport Gap
from gappy.gapobj cimport GapObj
from gappy.gap_includes cimport Obj

from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.structure.sage_object cimport SageObject

cdef GapElement make_any_gap_element(parent, GapObj obj)
cdef GapElement make_GapElement(parent, GapObj obj)
cdef GapElement_List make_GapElement_List(parent, GapObj obj)
cdef GapElement_Record make_GapElement_Record(parent, GapObj obj)
cdef GapElement_Integer make_GapElement_Integer(parent, GapObj obj)
cdef GapElement_Rational make_GapElement_Rational(parent, GapObj obj)
cdef GapElement_String make_GapElement_String(parent, GapObj obj)
cdef GapElement_Boolean make_GapElement_Boolean(parent, GapObj obj)
cdef GapElement_Function make_GapElement_Function(parent, GapObj obj)
cdef GapElement_Permutation make_GapElement_Permutation(parent, GapObj obj)


cdef class GapElement(RingElement):

    # the pointer to the gappy GapObj
    cdef GapObj obj

    # comparison
    cdef bint _compare_by_id
    cpdef _set_compare_by_id(self)
    cpdef _assert_compare_by_id(self)

    cdef _initialize(self, parent, GapObj obj)
    cpdef is_bool(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _mod_(self, other)
    cpdef _pow_(self, other)

    cpdef GapElement deepcopy(self, bint mut)

cdef class GapElement_Integer(GapElement):
    pass

cdef class GapElement_Rational(GapElement):
    pass

cdef class GapElement_IntegerMod(GapElement):
    cpdef GapElement_Integer lift(self)

cdef class GapElement_FiniteField(GapElement):
    cpdef GapElement_Integer lift(self)

cdef class GapElement_Cyclotomic(GapElement):
    pass

cdef class GapElement_Ring(GapElement):
    pass

cdef class GapElement_String(GapElement):
    pass

cdef class GapElement_Boolean(GapElement):
    pass

cdef class GapElement_Function(GapElement):
    pass

cdef class GapElement_MethodProxy(GapElement_Function):
    cdef GapElement first_argument

cdef class GapElement_Record(GapElement):
    pass

cdef class GapElement_List(GapElement):
    pass

cdef class GapElement_Permutation(GapElement):
    pass
