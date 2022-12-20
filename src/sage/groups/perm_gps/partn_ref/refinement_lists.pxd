# sage_setup: distribution = sagemath-gap
#*****************************************************************************
#       Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#       Copyright (C) 2009 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .data_structures cimport *


# name of the three functions to customize
cdef int refine_list(PartitionStack *, void *, int *, int)
cdef int compare_lists(int *, int *, void *, void *, int)
cdef bint all_list_children_are_equivalent(PartitionStack *, void *)
