#*****************************************************************************
#       Copyright (C) 2014 Sage
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#
###############################################################################

cimport sage.structure.parent

cdef class RecursivelyEnumeratedSet(sage.structure.parent.Parent):
    cdef readonly _seeds
    cdef public successors
    cdef readonly str _enumeration
    cdef readonly _max_depth
    cdef readonly _graded_component
    cdef readonly _graded_component_it

    cpdef seeds(self)
    cpdef graded_component(self, depth)

cdef class RecursivelyEnumeratedSet_symmetric(RecursivelyEnumeratedSet):
    cdef set _get_next_graded_component(self, set A, set B)

cdef class RecursivelyEnumeratedSet_graded(RecursivelyEnumeratedSet):
    cdef set _get_next_graded_component(self, set B)

