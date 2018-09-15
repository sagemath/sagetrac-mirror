#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .types cimport libGAP_Obj


############################################################################
### Hooking into the GAP memory management #################################
############################################################################

cdef class ObjWrapper(object):
    cdef libGAP_Obj value

cdef ObjWrapper wrap_obj(libGAP_Obj obj)

# returns the refcount dictionary for debugging purposes
cpdef get_owned_objects()

# Reference count GAP objects that you want to prevent from being
# garbage collected
cdef void reference_obj(libGAP_Obj obj)
cdef void dereference_obj(libGAP_Obj obj)

# callback from the GAP memory manager so we can mark all_gap_elements.values()
cdef void gasman_callback()


############################################################################
### Initialization of libGAP ###############################################
############################################################################

cdef initialize()


############################################################################
### Evaluate string in GAP #################################################
############################################################################

# Evaluate a string
cdef libGAP_Obj gap_eval(str gap_string) except? NULL


############################################################################
### Debug functions ########################################################
############################################################################

# Return details of the GAP memory pool
cpdef memory_usage()
