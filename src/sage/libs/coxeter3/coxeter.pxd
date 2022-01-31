#*****************************************************************************
#       Copyright (C) 2009-2013 Mike Hansen <mhansen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object cimport SageObject
from .decl cimport *

cdef class String:
    cdef c_String x

cdef class Type:
    cdef c_Type x

cdef class CoxGroup(SageObject):
    cdef c_CoxGroup* x
    cdef object cartan_type
    cdef dict in_ordering
    cdef dict out_ordering
    cpdef object full_context(self)

cdef class CoxGroupElement:
    cdef c_CoxWord word
    cdef c_CoxGroup* group
    cdef CoxGroup _parent_group
    cdef CoxGroupElement _new(self)
    cpdef CoxGroup parent_group(self)

cdef class CoxGraph:
    cdef c_CoxGraph x
