#*****************************************************************************
#  Copyright (C) 2012 Florent Hivert <Florent.Hivert at lri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.type cimport type

from sage.misc.nested_class cimport NestedClassMetaclass

cdef class ClasscallMetaclass(NestedClassMetaclass):
    cdef object classcall
    cdef object classget
    cdef object classcontains
