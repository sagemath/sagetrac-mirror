"""
Variable objects for SCIP.

AUTHOR: 

- Martin Albrecht (2010-11, initial version)
"""

##############################################################################
#       Copyright (C) 2010 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################


from decl cimport *
from scip cimport SCIP

cdef class Variable(SCIPObject):
    def __cinit__(self, SCIP scip, name, lb, ub, obj, long vartype):
        """
        EXAMPLE::

            sage: from sage.libs.scip import *                 # optional - SCIP
            sage: from sage.libs.scip.variable import Variable # optional - SCIP
            sage: scip = SCIP(name="99 problems")              # optional - SCIP
            sage: X = Variable(scip, "X", None, None, None, vartype=scip.VARTYPE_INTEGER); X # optional - SCIP
            X
            sage: X.is_integer()                               # optional - SCIP
            True
        """
        self._scip = scip

        cdef SCIP_RETCODE _status_
        cdef double _lb, _ub, _obj

        if lb is None:
            if vartype == SCIP_VARTYPE_BINARY:
                _lb = 0.0
            else:
                _lb = -SCIPinfinity(self._scip._scip)
        else:
            _lb = float(lb)

        if ub is None:
            if vartype == SCIP_VARTYPE_BINARY:
                _ub = 1.0
            else:
                _ub = SCIPinfinity(self._scip._scip)
        else:
            _ub = float(ub)

        if obj is None:
            _obj = 0.0
        else:
            _obj = obj

        if name is None:
            name = "x"

        _status_ = SCIPcreateVar(self._scip._scip, &self._var, name, _lb, _ub, _obj, vartype, 1, 0, NULL, NULL, NULL, NULL, NULL)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d initialising SCIP variable."%_status_)

        _status_ = SCIPaddVar(self._scip._scip, self._var)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d adding variable to SCIP instance."%_status_)        

    def __dealloc__(self):
        """
        TESTS::

            sage: from sage.libs.scip.scip import SCIP          # optional - SCIP
            sage: from sage.libs.scip.variable import Variable  # optional - SCIP
            sage: s = SCIP()                                    # optional - SCIP
            sage: for i in range(100):                          # optional - SCIP
            ...       v = Variable(s, str(i), float(i), float(i+1), None, vartype=s.VARTYPE_INTEGER) # optional - SCIP
            ...       del v                                     # optional - SCIP
        """
        if not self._scip or not self._var:
            return
        cdef SCIP_RETCODE _status_ 
        _status_ = SCIPreleaseVar(self._scip._scip, &self._var)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d deleting SCIP variable."%_status_)
                
    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import *                 # optional - SCIP
            sage: from sage.libs.scip.variable import Variable # optional - SCIP
            sage: scip = SCIP(name="99 problems")              # optional - SCIP
            sage: X = Variable(scip, "X", None, None, None, vartype=scip.VARTYPE_INTEGER) # optional - SCIP
            sage: repr(X) # indirect doctest                   # optional - SCIP
            'X'
        """
        return self._var.name

    def is_binary(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import *                 # optional - SCIP
            sage: from sage.libs.scip.variable import Variable # optional - SCIP
            sage: scip = SCIP("99 problems")                   # optional - SCIP
            sage: X = Variable(scip, "X", None, None, None, vartype=scip.VARTYPE_BINARY); X # optional - SCIP
            X
            sage: X.is_binary()                                # optional - SCIP
            True
        """
        return SCIPvarGetType(self._var) == SCIP_VARTYPE_BINARY

    def is_continuous(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import *                 # optional - SCIP
            sage: from sage.libs.scip.variable import Variable # optional - SCIP
            sage: scip = SCIP("99 problems")                   # optional - SCIP
            sage: X = Variable(scip, "X", None, None, None, vartype=scip.VARTYPE_CONTINUOUS); X # optional - SCIP
            X
            sage: X.is_continuous()                            # optional - SCIP
            True
        """
        return SCIPvarGetType(self._var) == SCIP_VARTYPE_CONTINUOUS

    def is_integer(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import *                 # optional - SCIP
            sage: from sage.libs.scip.variable import Variable # optional - SCIP
            sage: scip = SCIP("99 problems")                   # optional - SCIP
            sage: X = Variable(scip, "X", None, None, None, vartype=scip.VARTYPE_INTEGER); X # optional - SCIP
            X
            sage: X.is_integer()                               # optional - SCIP
            True
        """
        return SCIPvarGetType(self._var) == SCIP_VARTYPE_INTEGER

    def bounds(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import *                 # optional - SCIP
            sage: from sage.libs.scip.variable import Variable # optional - SCIP
            sage: scip = SCIP("99 problems")                   # optional - SCIP
            sage: X = Variable(scip, "X", None, None, None, vartype=scip.VARTYPE_INTEGER) # optional - SCIP
            sage: X.bounds()                                   # optional - SCIP
            (None, None)

            sage: X = Variable(scip, "X", -1, 3, None, vartype=scip.VARTYPE_INTEGER) # optional - SCIP
            sage: X.bounds()                                   # optional - SCIP
            (-1.0, 3.0)

            sage: X = Variable(scip, "X", 2^15, 2^16, None, vartype=scip.VARTYPE_INTEGER) # optional - SCIP
            sage: X.bounds()                                   # optional - SCIP
            (32768.0, 65536.0)
        """
        if SCIPvarIsTransformed(self._var):
            lb = SCIPvarGetLbGlobal(self._var)
            ub = SCIPvarGetUbGlobal(self._var)
        else:
            lb = SCIPvarGetLbOriginal(self._var)
            ub = SCIPvarGetUbOriginal(self._var)
        if ub == SCIPinfinity(self._scip._scip):
            ub = None
        if lb == -SCIPinfinity(self._scip._scip):
            lb = None
        return lb,ub

    def __hash__(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import *                 # optional - SCIP
            sage: from sage.libs.scip.variable import Variable # optional - SCIP
            sage: scip = SCIP("99 problems")                   # optional - SCIP
            sage: x = Variable(scip, "X", None, None, None, vartype=scip.VARTYPE_CONTINUOUS) # optional - SCIP
            sage: y = Variable(scip, "X", None,  1.0, None, vartype=scip.VARTYPE_CONTINUOUS) # optional - SCIP
            sage: z = Variable(scip, "X", None, None,  1.0, vartype=scip.VARTYPE_CONTINUOUS) # optional - SCIP

            sage: hash(x) == hash(y) # optional - SCIP
            False
            sage: hash(y) == hash(z) # optional - SCIP
            False
            sage: hash(z) == hash(x) # optional - SCIP
            True
        """
        cdef double lb = SCIPvarGetLbOriginal(self._var)
        cdef double ub = SCIPvarGetUbOriginal(self._var)

        return hash(self._var.name) ^ SCIPvarGetType(self._var)<<10 ^ hash(lb + 13.37 * ub)

    def index(self):
        """
        Return the index of this variable.

        EXAMPLE::

            sage: from sage.libs.scip import *                 # optional - SCIP
            sage: from sage.libs.scip.variable import Variable # optional - SCIP
            sage: scip = SCIP("99 problems")                   # optional - SCIP
            sage: x = scip.add_variable(); x                   # optional - SCIP
            0
            sage: scip.variable(x).index()                     # optional - SCIP
            0
        """
        return self._var.index
