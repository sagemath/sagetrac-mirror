"""
Constraint objects for SCIP.

AUTHOR: 

- Martin Albrecht (2010-11, initial version)
"""

##############################################################################
#  Copyright (C) 2010,2011 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

include "sage/ext/stdsage.pxi"

from decl cimport *

from scip cimport SCIP
from variable cimport Variable

cdef class Constraint(SCIPObject):
    """
    A superclass for all constraints available in SCIP.
    """
    def __dealloc__(self):
        if self._cons:
            SCIPreleaseCons(self._scip._scip, &self._cons)
    
cdef class LinearConstraint(Constraint):
    """
    SCIP's linear constraint.
    """
    def __cinit__(self, SCIP scip, constraints, name=None, lower_bound=None, upper_bound=None,
                  initial=True, 
                  separate=True, 
                  enforce=True, 
                  check=True,
                  propagate=True,
                  local=False,
                  modifiable=False,
                  dynamic=False,
                  removable=False,
                  stickingatnode=False):
        """
        Create a new linear constraint.

        INPUT:

        - ``scip`` - the SCIP instance this constraint will live
          in. The constraint is *not* added to the SCIP problem.

        - ``constraints`` - list of (variable,coefficients) tuples

        - ``name`` - some name (default: ``None``)

        - ``lower_bound`` - the lower bound for this constraint

        - ``upper_bound`` - the upper bound for this constraint

        - ``initial`` - should the LP relaxation of constraint be in
          the initial LP? Usually set to ``True``. Set to `False`` for
          'lazy constraints'. (default: ``True``, SCIP specific)

        - ``separate`` - should the constraint be separated during LP
          processing? (default: ``True``, SCIP specific)

        - ``enforce`` - should the constraint be enforced during node
          processing? ``True`` for model constraints, ``False`` for
          additional, redundant constraints.  (default: ``True``, SCIP
          specific)

        - ``check`` - should the constraint be checked for
          feasibility? ``True`` for model constraints, ``False`` for
          additional, redundant constraints.  (default: ``True``, SCIP
          specific)

        - ``propagate`` - should the constraint be propagated during
          node processing? (default: `True``, SCIP specific)

        - ``local`` - is constraint only valid locally? Usually set to
          ``False``. Has to be set to ``True``, e.g., for branching
          constraints.  (default: ``False``, SCIP specific)

        - ``modifiable`` - is constraint modifiable (subject to column
          generation)? Usually set to ``False``. In column generation
          applications, set to ``True`` if pricing adds coefficients
          to this constraint.  (default: ``False``, SCIP specific)

        - ``dynamic`` - Is constraint subject to aging? Usually set to
          ``False``. Set to ``True`` for own cuts which are seperated
          as constraints.  (default: ``False``, SCIP specific)

        - ``removable`` - should the relaxation be removed from the LP
          due to aging or cleanup? Usually set to ``False``. Set to
          ``True`` for 'lazy constraints' and 'user cuts'. (default:
          ``False``, SCIP specific)

        - ``stickingatnode`` - should the constraint always be kept at
          the node where it was added, even if it may be moved to a
          more global node? Usually set to ``False``. Set to ``True``
          to for constraints that represent node data. (default:
          ``False``, SCIP specific)

        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import LinearConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: LinearConstraint(scip, [ (x,1.0), (y,2.0) ], lower_bound=-1.0, upper_bound=1.0) # optional - SCIP
            -1.000000 <= 1.000000*x0 + 2.000000*x1 <= 1.000000
        """
        self._scip = scip

        if name is None:
            name = "<some constraint>"
        self._name = name

        cdef double _lhs, _rhs
        cdef nlinterms = len(constraints)
        cdef SCIP_VAR **linvars
        cdef double *lincoefs


        _lhs = -SCIPinfinity(self._scip._scip) if lower_bound is None else float(lower_bound)
        _rhs =  SCIPinfinity(self._scip._scip) if upper_bound is None else float(upper_bound)

        for v,c in constraints:
            if not (PY_TYPE_CHECK(v, Variable) and (<Variable>v)._scip is self._scip):
                raise TypeError("Variable must be of type SCIP Variable of this SCIP instance")

        if nlinterms > 0:
            linvars = <SCIP_VAR**>sage_malloc(sizeof(SCIP_VAR*)*nlinterms)
            lincoefs = <double*>sage_malloc(sizeof(double)*nlinterms)

            for i,(v,c) in enumerate(constraints):
                lincoefs[i] = float(c)
                linvars[i] = (<Variable>v)._var
        else:
            linvars = NULL
            lincoefs = NULL

        cdef SCIP_RETCODE _retstat_
        _retstat_ = SCIPcreateConsLinear(self._scip._scip, &self._cons, name, 
                                         nlinterms, linvars, lincoefs,
                                         _lhs, _rhs, initial, separate, enforce, check,
                                         propagate, local, modifiable, dynamic, removable, 
                                         stickingatnode)

        if nlinterms:
            sage_free(linvars)
            sage_free(lincoefs)

        if _retstat_ != SCIP_OKAY:
            raise RuntimeError("Error %d initialising SCIP linear constraint."%_retstat_)

    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import LinearConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable(integer=True)) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: LinearConstraint(scip, [ (x,1.0), (y,2.0) ], lower_bound=None, upper_bound=1.0) # optional - SCIP, indirect doctest
            -inf <= 1.000000*x0 + 2.000000*x1 <= 1.000000
        """
        cdef double lhs = SCIPgetLhsLinear(self._scip._scip, self._cons)
        cdef double rhs = SCIPgetRhsLinear(self._scip._scip, self._cons)

        cdef int i
        cdef int n = SCIPgetNVarsLinear(self._scip._scip, self._cons)
        cdef SCIP_VAR **v = SCIPgetVarsLinear(self._scip._scip, self._cons)
        cdef double *c = SCIPgetValsLinear(self._scip._scip, self._cons)

        if lhs == -SCIPinfinity(self._scip._scip):
            lhs_s = "-inf <= "
        else:
            lhs_s = "%f <= "%(lhs,)

        if rhs == SCIPinfinity(self._scip._scip):
            rhs_s = " <= inf"
        else:
            rhs_s = " <= %f"%(rhs,)

        c_s = []
        for i in range(n):
            c_s.append("%f*%s"%(c[i],v[i].name))
        c_s = " + ".join(c_s)

        return "%s%s%s"%(lhs_s, c_s, rhs_s)

    def coefficients(self):
        """
        Return the list of coefficients.

        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import LinearConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable(integer=True)) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: l = LinearConstraint(scip, [ (x,1.0), (z,2.0) ], lower_bound=None, upper_bound=1.0) # optional - SCIP
            sage: l.coefficients() # optional - SCIP
            [1.0, 2.0]
            
        """
        cdef int n = SCIPgetNVarsLinear(self._scip._scip, self._cons)
        cdef double *c = SCIPgetValsLinear(self._scip._scip, self._cons)

        coeffs = []
        for i in range(n):
            coeffs.append(c[i])
        return coeffs

    def variables(self):
        """
        Return the list of variables.

        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import LinearConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable(integer=True)) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: l = LinearConstraint(scip, [ (x,1.0), (z,2.0) ], lower_bound=None, upper_bound=1.0) # optional - SCIP
            sage: l.variables() # optional - SCIP
            [x0, x2]
        """
        cdef int n = SCIPgetNVarsLinear(self._scip._scip, self._cons)
        cdef SCIP_VAR **v = SCIPgetVarsLinear(self._scip._scip, self._cons)

        variables = []
        for i in range(n):
            variables.append(self._scip._variables[v[i].index])
        return variables

    def bounds(self):
        """
        Return the lower and upper bound of this constraint.

        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import LinearConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable(integer=True)) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: l = LinearConstraint(scip, [ (x,1.0), (z,2.0) ], lower_bound=None, upper_bound=1.0) # optional - SCIP
            sage: l.bounds() # optional - SCIP
            (None, 1.0)            
        """
        cdef double lhs = SCIPgetLhsLinear(self._scip._scip, self._cons)
        cdef double rhs = SCIPgetRhsLinear(self._scip._scip, self._cons)

        ret = lhs,rhs

        if lhs == -SCIPinfinity(self._scip._scip):
            ret = None,ret[1]

        if rhs == SCIPinfinity(self._scip._scip):
            ret = ret[0],None
            
        return ret

cdef class ANDConstraint(Constraint):
    """
    SCIP's AND constraint (res = a & b & c)
    """
    def __cinit__(self, SCIP scip, Variable res, variables, name=None, 
                  initial=True, 
                  separate=True, 
                  enforce=True, 
                  check=True,
                  propagate=True,
                  local=False,
                  modifiable=False,
                  dynamic=False,
                  removable=False,
                  stickingatnode=False):
        """
        Construct a new AND constraint

          `res = a & b & c` for `a,b,c \in variables.`

        INPUT:

        - ``scip`` - the SCIP instance this constraint will live
          in. The constraint is *not* added to the SCIP problem.

        - ``res`` - a SCIP variable

        - ``variables`` - list of SCIP variables

        - ``name`` - some name (default: ``None``)

        .. note::
        
            For the remaining parameters see :class:`LinearConstraint`.

        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import ANDConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: ANDConstraint(scip, x, [y,z])          # optional - SCIP
            x0 = x1 & x2
            
            sage: ANDConstraint(scip, x, [])             # optional - SCIP
            Traceback (most recent call last):
            ...
            ValueError: Parameter 'variables' must have at least one element.
        """
        self._scip = scip

        if name is None:
            name = "<some constraint>"
        self._name = name

        cdef int i
        cdef int n = len(variables)
        if n < 1:
            raise ValueError("Parameter 'variables' must have at least one element.")
        for i,v in enumerate(variables):
            if not (PY_TYPE_CHECK(v, Variable) and (<Variable>v)._scip is self._scip):
                raise TypeError("Variable must be of type SCIP Variable of this SCIP instance")
            
        cdef SCIP_VAR **_vars = <SCIP_VAR **>sage_malloc(sizeof(SCIP_VAR)*n)

        for i,v in enumerate(variables):
            _vars[i] = (<Variable>v)._var

        cdef SCIP_RETCODE _retstat_
        _retstat_ = SCIPcreateConsAnd(self._scip._scip, &self._cons, name, res._var, n, _vars,
                                      initial, separate, enforce, check, propagate, local, 
                                      modifiable, dynamic, removable, stickingatnode)

        sage_free(_vars)

        if _retstat_ != SCIP_OKAY:
            raise RuntimeError("Error %d initialising SCIP AND constraint."%_retstat_)

    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import ANDConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: ANDConstraint(scip, x, [y,z]) # optional - SCIP, indirect doctest
            x0 = x1 & x2
        """
        cdef int i
        cdef int n = SCIPgetNVarsAnd (self._scip._scip, self._cons)
        cdef SCIP_VAR **_vars = SCIPgetVarsAnd (self._scip._scip, self._cons)
        cdef SCIP_VAR *_var = SCIPgetResultantAnd (self._scip._scip, self._cons)

        l = []
        for i in range(n):
            l.append(_vars[i].name)
        return "%s = %s"%(_var.name," & ".join(l))


cdef class ORConstraint(Constraint):
    """
    SCIP's OR constraint (res = a | b | c)
    """
    def __cinit__(self, SCIP scip, Variable res, variables, name=None, 
                  initial=True, 
                  separate=True, 
                  enforce=True, 
                  check=True,
                  propagate=True,
                  local=False,
                  modifiable=False,
                  dynamic=False,
                  removable=False,
                  stickingatnode=False):
        """
        Construct a new OR constraint

          `res = a | b | c` for `a,b,c \in variables.`

        INPUT:

        - ``scip`` - the SCIP instance this constraint will live
          in. The constraint is *not* added to the SCIP problem.

        - ``res`` - a SCIP variable

        - ``variables`` - list of SCIP variables

        - ``name`` - some name (default: ``None``)

        .. note::
        
            For the remaining parameters see :class:`LinearConstraint`.

        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import ORConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: ORConstraint(scip, x, [y,z])          # optional - SCIP
            x0 = x1 | x2
            
            sage: ORConstraint(scip, x, [])             # optional - SCIP
            Traceback (most recent call last):
            ...
            ValueError: Parameter 'variables' must have at least one element.
        """
        self._scip = scip

        if name is None:
            name = "<some constraint>"
        self._name = name

        cdef int i
        cdef int n = len(variables)
        if n < 1:
            raise ValueError("Parameter 'variables' must have at least one element.")
        for i,v in enumerate(variables):
            if not (PY_TYPE_CHECK(v, Variable) and (<Variable>v)._scip is self._scip):
                raise TypeError("Variable must be of type SCIP Variable of this SCIP instance")
            
        cdef SCIP_VAR **_vars = <SCIP_VAR **>sage_malloc(sizeof(SCIP_VAR)*n)

        for i,v in enumerate(variables):
            _vars[i] = (<Variable>v)._var

        cdef SCIP_RETCODE _retstat_
        _retstat_ = SCIPcreateConsOr(self._scip._scip, &self._cons, name, res._var, n, _vars,
                                      initial, separate, enforce, check, propagate, local, 
                                      modifiable, dynamic, removable, stickingatnode)

        sage_free(_vars)

        if _retstat_ != SCIP_OKAY:
            raise RuntimeError("Error %d initialising SCIP OR constraint."%_retstat_)

    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import ORConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: ORConstraint(scip, x, [y,z]) # optional - SCIP, indirect doctest
            x0 = x1 | x2
        """
        cdef int i
        cdef int n = SCIPgetNVarsOr (self._scip._scip, self._cons)
        cdef SCIP_VAR **_vars = SCIPgetVarsOr (self._scip._scip, self._cons)
        cdef SCIP_VAR *_var = SCIPgetResultantAnd(self._scip._scip, self._cons)

        l = []
        for i in range(n):
            l.append(_vars[i].name)
        return "%s = %s"%(_var.name," | ".join(l))


cdef class XORConstraint(Constraint):
    """
    SCIP's XOR constraint (rhs = a ^ b ^ c)
    """
    def __cinit__(self, SCIP scip, rhs, variables, name=None, 
                  initial=True, 
                  separate=True, 
                  enforce=True, 
                  check=True,
                  propagate=True,
                  local=False,
                  modifiable=False,
                  dynamic=False,
                  removable=False,
                  stickingatnode=False):
        """
        Construct a new XOR constraint

          `rhs = a ^ b ^ c` for `a,b,c \in variables.`

        INPUT:

        - ``scip`` - the SCIP instance this constraint will live
          in. The constraint is *not* added to the SCIP problem.

        - ``rhs`` - a boolean value

        - ``variables`` - list of SCIP variables

        - ``name`` - some name (default: ``None``)

        .. note::
        
            For the remaining parameters see :class:`LinearConstraint`.

        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import XORConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: XORConstraint(scip, 0, [y,z])          # optional - SCIP
            0 = x1 ^ x2
            
            sage: XORConstraint(scip, x, [])             # optional - SCIP
            Traceback (most recent call last):
            ...
            ValueError: Parameter 'variables' must have at least one element.
        """
        self._scip = scip

        if name is None:
            name = "<some constraint>"
        self._name = name

        cdef int i
        cdef int n = len(variables)
        if n < 1:
            raise ValueError("Parameter 'variables' must have at least one element.")
        for i,v in enumerate(variables):
            if not (PY_TYPE_CHECK(v, Variable) and (<Variable>v)._scip is self._scip):
                raise TypeError("Variable must be of type SCIP Variable of this SCIP instance")

            
        cdef SCIP_VAR **_vars = <SCIP_VAR **>sage_malloc(sizeof(SCIP_VAR*)*n)

        for i,v in enumerate(variables):
            _vars[i] = (<Variable>v)._var

        cdef SCIP_RETCODE _retstat_
        _retstat_ = SCIPcreateConsXor(self._scip._scip, &self._cons, name, rhs, n, _vars,
                                      initial, separate, enforce, check, propagate, local, 
                                      modifiable, dynamic, removable, stickingatnode)

        sage_free(_vars)

        if _retstat_ != SCIP_OKAY:
            raise RuntimeError("Error %d initialising SCIP XOR constraint."%_retstat_)

    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import XORConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: XORConstraint(scip, 0, [y,z]) # optional - SCIP, indirect doctest
            0 = x1 ^ x2
        """
        cdef SCIP_RETCODE _status_
        cdef int n = SCIPgetNVarsXor(self._scip._scip, self._cons)
        cdef SCIP_VAR **_vars = SCIPgetVarsXor(self._scip._scip, self._cons)
        cdef bint rhs = SCIPgetRhsXor(self._scip._scip, self._cons)
        
        l = []
        for i in range(n):
            l.append(_vars[i].name)
        return "%d = %s"%(rhs," ^ ".join(l))

cdef class LogicORConstraint(Constraint):
    """
    SCIP's Logical OR operator.
    """
    def __cinit__(self, SCIP scip, variables, name=None, 
                  initial=True, 
                  separate=True, 
                  enforce=True, 
                  check=True,
                  propagate=True,
                  local=False,
                  modifiable=False,
                  dynamic=False,
                  removable=False,
                  stickingatnode=False):
        """
        Construct a new LogicOR constraint

          `a` or `b` or `c` for `a,b,c \in variables.`

        INPUT:

        - ``scip`` - the SCIP instance this constraint will live
          in. The constraint is *not* added to the SCIP problem.

        - ``variables`` - list of SCIP variables

        - ``name`` - some name (default: ``None``)

        .. note::
        
            For the remaining parameters see :class:`LinearConstraint`.

        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import LogicORConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: LogicORConstraint(scip, [y,z])          # optional - SCIP
            x1 or x2
            
            sage: LogicORConstraint(scip, [])             # optional - SCIP
            Traceback (most recent call last):
            ...
            ValueError: Parameter 'variables' must have at least one element.
        """

        self._scip = scip

        if name is None:
            name = "<some constraint>"
        self._name = name

        cdef int i
        cdef int n = len(variables)
        if n < 1:
            raise ValueError("Parameter 'variables' must have at least one element.")
        for i,v in enumerate(variables):
            if not (PY_TYPE_CHECK(v, Variable) and (<Variable>v)._scip is self._scip):
                raise TypeError("Variable must be of type SCIP Variable of this SCIP instance")
            
        cdef SCIP_VAR **_vars = <SCIP_VAR **>sage_malloc(sizeof(SCIP_VAR)*n)

        for i,v in enumerate(variables):
            _vars[i] = (<Variable>v)._var

        cdef SCIP_RETCODE _retstat_
        _retstat_ = SCIPcreateConsLogicor(self._scip._scip, &self._cons, name, n, _vars,
                                          initial, separate, enforce, check, propagate, local, 
                                          modifiable, dynamic, removable, stickingatnode)
        
        sage_free(_vars)

        if _retstat_ != SCIP_OKAY:
            raise RuntimeError("Error %d initialising SCIP XOR constraint."%_retstat_)

    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import LogicORConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: LogicORConstraint(scip, [y,z]) # optional - SCIP, indirect doctest
            x1 or x2
        """
        cdef SCIP_RETCODE _status_
        cdef int n = SCIPgetNVarsLogicor(self._scip._scip, self._cons)
        cdef SCIP_VAR **_vars = SCIPgetVarsLogicor(self._scip._scip, self._cons)
        
        l = []
        for i in range(n):
            l.append(_vars[i].name)
        return " or ".join(l)

cdef class SetPPCConstraint(Constraint):
    def __cinit__(self, SCIP scip, variables, type, name=None, 
                  initial=True, 
                  separate=True, 
                  enforce=True, 
                  check=True,
                  propagate=True,
                  local=False,
                  modifiable=False,
                  dynamic=False,
                  removable=False,
                  stickingatnode=False):
        """
        Construct a new partition, packing or covering constraint

          sum(a,b,c) <=/==/=> 1

        INPUT:

        - ``scip`` - the SCIP instance this constraint will live
          in. The constraint is *not* added to the SCIP problem.

        - ``variables`` - list of SCIP variables
        
        - ``type`` - either ``partition``, ``packing`` or ``covering``.

        - ``name`` - some name (default: ``None``)

        .. note::
        
            For the remaining parameters see :class:`LinearConstraint`.

        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import SetPPCConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: SetPPCConstraint(scip, [x,y,z], "partition") # optional - SCIP
            sum(x0, x1, x2) == 1

            sage: SetPPCConstraint(scip, [x,y,z], "packing")   # optional - SCIP
            sum(x0, x1, x2) <= 1

            sage: SetPPCConstraint(scip, [x,y,z], "covering")  # optional - SCIP
            sum(x0, x1, x2) >= 1

            sage: SetPPCConstraint(scip, [], "covering")  # optional - SCIP
            Traceback (most recent call last):
            ...
            ValueError: Parameter 'variables' must have at least one element.
        """

        self._scip = scip

        if name is None:
            name = "<some constraint>"
        self._name = name

        cdef int i
        cdef int n = len(variables)
        if n < 1:
            raise ValueError("Parameter 'variables' must have at least one element.")
        for i,v in enumerate(variables):
            if not (PY_TYPE_CHECK(v, Variable) and (<Variable>v)._scip is self._scip):
                raise TypeError("Variable must be of type SCIP Variable of this SCIP instance")

            
        cdef SCIP_VAR **_vars = <SCIP_VAR **>sage_malloc(sizeof(SCIP_VAR*)*n)

        for i,v in enumerate(variables):
            _vars[i] = (<Variable>v)._var

        cdef SCIP_RETCODE _retstat_
        if type == 'partition':
            _retstat_ = SCIPcreateConsSetpart(self._scip._scip, &self._cons, name, n, _vars,
                                              initial, separate, enforce, check, propagate, local, 
                                              modifiable, dynamic, removable, stickingatnode)
        elif type == 'packing':
            _retstat_ = SCIPcreateConsSetpack(self._scip._scip, &self._cons, name, n, _vars,
                                              initial, separate, enforce, check, propagate, local, 
                                              modifiable, dynamic, removable, stickingatnode)
        elif type == 'covering':
            _retstat_ = SCIPcreateConsSetcover(self._scip._scip, &self._cons, name, n, _vars,
                                               initial, separate, enforce, check, propagate, local, 
                                               modifiable, dynamic, removable, stickingatnode)
        else:
            raise ValueError("Type '%s' unknown."%type)

        sage_free(_vars)

        if _retstat_ != SCIP_OKAY:
            raise RuntimeError("Error %d initialising SCIP XOR constraint."%_retstat_)

    def _repr_(self):
        """
        EXAMPLE::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.constraint import SetPPCConstraint # optional - SCIP
            sage: scip = SCIP() # optional - SCIP
            sage: x = scip.variable(scip.add_variable()) # optional - SCIP
            sage: y = scip.variable(scip.add_variable()) # optional - SCIP
            sage: z = scip.variable(scip.add_variable()) # optional - SCIP
            sage: SetPPCConstraint(scip, [x,y,z], "partition") # optional - SCIP, indirect doctest
            sum(x0, x1, x2) == 1

            sage: SetPPCConstraint(scip, [x,y,z], "packing")   # optional - SCIP, indirect doctest
            sum(x0, x1, x2) <= 1

            sage: SetPPCConstraint(scip, [x,y,z], "covering")  # optional - SCIP, indirect doctest
            sum(x0, x1, x2) >= 1
        """
        cdef SCIP_RETCODE _status_
        cdef int n = SCIPgetNVarsSetppc(self._scip._scip, self._cons)
        cdef SCIP_VAR **_vars = SCIPgetVarsSetppc(self._scip._scip, self._cons)
        cdef SCIP_SetppcType _type = SCIPgetTypeSetppc(self._scip._scip, self._cons)
      
        l = []
        for i in range(n):
            l.append(_vars[i].name)
        variables = ", ".join(l)

        if _type == SCIP_SETPPCTYPE_PARTITIONING:
            return "sum(%s) == 1"%variables
        elif _type == SCIP_SETPPCTYPE_PACKING:
            return "sum(%s) <= 1"%variables
        elif _type == SCIP_SETPPCTYPE_COVERING:
            return "sum(%s) >= 1"%variables
 
cdef class QuadraticConstraint(Constraint):
    """
    SCIP's quadratic constraints"
    """
    def __cinit__(self, SCIP scip, linterms, quadterms, name, lhs, rhs,
                  initial=True, 
                  separate=True, 
                  enforce=True, 
                  check=True,
                  propagate=True,
                  local=False,
                  modifiable=False,
                  dynamic=False,
                  removable=False,
                  stickingatnode=False):
        """
        Create a new quadratic constraint.

        `'lhs <= \sum_{v,c \in linterms} c*v + \sum_{v1,v2,c12 \in quadterms} c12*v1*v2 <= rhs`

        INPUT:

        - ``scip`` - the SCIP instance this constraint will live
          in. The constraint is *not* added to the SCIP problem.

        - ``linterms`` - list of ``(v,c)`` tuples

        - ``qyadterms`` - list of ``(v1,v2,c12)`` tuples

        - ``name`` - some name (default: ``None``)

        - ``lhs`` - the lower bound for this constraint

        - ``rhs`` - the upper bound for this constraint

        .. note::
        
            For the remaining parameters see :class:`LinearConstraint`.
        """
        cdef double _lhs, _rhs
        cdef int i, nlinterms, nquadterms
        cdef SCIP_VAR **linvars, **quadvars1, **quadvars2
        cdef double *lincoefs, *quadcoefs

        self._scip = scip

        if name is None:
            name = "<quadratic constraint>"
        self._name = name


        _lhs = -SCIPinfinity(self._scip._scip) if lhs is None else float(lhs)
        _rhs = SCIPinfinity(self._scip._scip) if rhs is None else float(rhs)

        cdef SCIP_RETCODE _retstat_

        nlinterms = len(linterms)
        nquadterms = len(quadterms)

        for v,c in linterms:
            if not (PY_TYPE_CHECK(v, Variable) and (<Variable>v)._scip is self._scip):
                raise TypeError("Variable must be of type SCIP Variable of this SCIP instance")

        for v1,v2,c in quadterms:
            if not (PY_TYPE_CHECK(v1, Variable) and (<Variable>v1)._scip is self._scip):
                raise TypeError("Variable must be of type SCIP Variable of this SCIP instance")
            if not (PY_TYPE_CHECK(v2, Variable) and (<Variable>v2)._scip is self._scip):
                raise TypeError("Variable must be of type SCIP Variable of this SCIP instance")

        if nlinterms > 0:
            linvars = <SCIP_VAR**>sage_malloc(sizeof(SCIP_VAR*)*nlinterms)
            lincoefs = <double*>sage_malloc(sizeof(double)*nlinterms)

            for i,(v,c) in enumerate(linterms):
                lincoefs[i] = float(c)
                linvars[i] = (<Variable>v)._var
        else:
            linvars = NULL
            lincoefs = NULL


        if nquadterms > 0:
            quadvars1 = <SCIP_VAR**>sage_malloc(sizeof(SCIP_VAR*) * nquadterms)
            quadvars2 = <SCIP_VAR**>sage_malloc(sizeof(SCIP_VAR*) * nquadterms)
            quadcoefs = <double*>sage_malloc(sizeof(double) * nquadterms)

            for i,(v1,v2,c) in enumerate(quadterms):
                quadcoefs[i] = float(c)
                quadvars1[i] = (<Variable>v1)._var
                quadvars2[i] = (<Variable>v2)._var
        else:
            quadvars1 = NULL
            quadvars2 = NULL
            quadcoefs = NULL
            
        _retstat_ =  SCIPcreateConsQuadratic(self._scip._scip, &self._cons, name,
                                             nlinterms, linvars, lincoefs, 
                                             nquadterms, quadvars1, quadvars2, quadcoefs, 
                                             _lhs, _rhs, 
                                             initial, separate, enforce, check, propagate, 
                                             local, modifiable, dynamic, removable)

        if nlinterms:
            sage_free(linvars)
            sage_free(lincoefs)
    
        if nquadterms:
            sage_free(quadvars1)
            sage_free(quadvars2)
            sage_free(quadcoefs)

        if _retstat_ != SCIP_OKAY:
            raise RuntimeError("Error %d initialising SCIP quadratic constraint."%_retstat_)

    def _repr_(self):
        return "QuadraticConstraint"
    #     cdef double *lhs = SCIPgetLhsQuadratic(self._scip._scip, self._cons)
    #     cdef double *rhs = SCIPgetRhsQuadratic(self._scip._scip, self._cons)
    #     cdef int i
    #     cdef int nlinterms = SCIPgetNLinearVarsQuadratic(self._scip._scip, self._cons)
    #     cdef int nquadterms = SCIPgetNQuadVarsQuadratic(self._scip._scip, self._cons)

    #     cdef SCIP_VAR **linvars = SCIPgetLinearVarsQuadratic (self._scip._scip, self._cons)
    #     cdef SCIP_VAR **quadvars = SCIPgetQuadVarsQuadratic (self._scip._scip, self._cons)
    #     cdef double *lincoefs = SCIPgetCoefsLinearVarsQuadratic (self._scip._scip, self._cons)
    #     cdef double *quadcoefs = SCIPgetLinearCoefsQuadVarsQuadratic (self._scip._scip, self._cons)
    #     cdef double *sqrcoefs = SCIPgetSqrCoefsQuadVarsQuadratic (self._scip._scip, self._cons)
        

    #     if lhs[0] == -SCIPinfinity(self._scip._scip):
    #         lhs_s = "-inf <= "
    #     else:
    #         lhs_s = "%f <= "%(lhs[0],)

    #     if rhs[0] == SCIPinfinity(self._scip._scip):
    #         rhs_s = " <= inf"
    #     else:
    #         rhs_s = " <= %f"%(rhs[0],)

    #     c_s = []
    #     for i in range(n):
    #         c_s.append("%f*%s"%(c[i],v[i].name))
    #     c_s = " + ".join(c_s)

    #     return "%s%s%s"%(lhs_s, c_s, rhs_s)
            
