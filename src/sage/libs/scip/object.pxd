"""
The base type for all objects in the SCIP wrapper.

AUTHOR: 

- Martin Albrecht (2010-11, initial version)
"""

##############################################################################
#       Copyright (C) 2010 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.structure.sage_object cimport SageObject
from scip cimport SCIP

cdef class SCIPObject(SageObject):
    cdef SCIP _scip

