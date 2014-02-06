"""
Interface to 'Solving Constraint Integer Programs' (SCIP).

SCIP is described by its developers as: "SCIP is a framework for Constraint
Integer Programming oriented towards the needs of Mathematical Programming
experts who want to have total control of the solution process and access
detailed information down to the guts of the solver. SCIP can also be used as a
pure MIP solver or as a framework for branch-cut-and-price.

SCIP is implemented as C callable library and provides C++ wrapper classes for
user plugins. It can also be used as a standalone program to solve mixed integer
programs given in MPS Format, in ILOG LP Format or as ZIMPL model."

SCIP is available at http://scip.zib.de which also hosts the optional Sage
package required for this interface.

AUTHOR:

- Martin Albrecht (2010-11, initial version)
"""

################################################################################
#    Copyright (C) 2010-2013 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
################################################################################

include "../../ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from decl cimport *
from variable cimport Variable
from constraint cimport *
from sage.misc.updown import UpDown

from sage.misc.misc import get_verbose
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.numerical.mip import MIPSolverException
from sage.rings.all import RR

cdef class SCIP(GenericBackend):
    """
    An instance of a SCIP integer constraint solver.

    EXAMPLE::

        sage: from sage.libs.scip import SCIP                        # optional - SCIP
        sage: scip = SCIP(maximization=True)                         # optional - SCIP
        sage: x = scip.add_variable(obj=1.0, integer=True)           # optional - SCIP
        sage: y = scip.add_variable(obj=1.0, integer=True)           # optional - SCIP
        sage: scip.add_linear_constraint([(x,0.5),(y,0.5)],-3.0,3.0) # optional - SCIP
        sage: scip.solve()                                           # optional - SCIP
        sage: scip.get_variable_value(0)                             # optional - SCIP
        6.0
        sage: scip.get_variable_value(1)                             # optional - SCIP
        -0.0

    .. note::

        In the context of this interface we identify a SCIP solver with a
        particular problem. Thus, at creation also an empty new SCIP problem is
        created.
    """
    VARTYPE_INTEGER = SCIP_VARTYPE_INTEGER
    VARTYPE_CONTINUOUS = SCIP_VARTYPE_CONTINUOUS
    VARTYPE_BINARY = SCIP_VARTYPE_BINARY

    def __cinit__(self, maximization=True, parameters=None, name=""):
        """
        Create a new SCIP solver instance.

        INPUT:

        - ``maximization`` - if `True` the objective sense is to maximize
          (default: ``True``)

        - ``parameters`` - filename of a SCIP parameters file (default:
          ``None``)

        - ``name`` - some name for the problem (default: '')

        EXAMPLES:

        One can either import the SCIP solver directly::

            sage: from sage.libs.scip import SCIP # optional - SCIP
            sage: SCIP()                          # optional - SCIP
            SCIP Constraint Integer Program "" ( maximization, 0 variables, 0 constraints )

        Alternatively, one may use :func:`get_solver`::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP") # optional - SCIP
            sage: p                               # optional - SCIP
            SCIP Constraint Integer Program "" ( maximization, 0 variables, 0 constraints )
        """
        cdef SCIP_RETCODE _status_ = SCIPcreate(&self._scip)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d initialising SCIP."%_status_)

        _status_ = SCIPincludeDefaultPlugins(self._scip)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d loading default SCIP plugins."%_status_)

        _restat_ = SCIPcreateProb(self._scip, name, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d initialising SCIP Problem."%_status_)

        if maximization:
            self.set_sense(1)
        else:
            self.set_sense(-1)

        if parameters is not None:
            self.read_parameters(parameters)

        self._variables = []
        self._constraints = []

        self.set_verbosity(get_verbose())

    def __dealloc__(self):
        """
        TESTS::

            sage: from sage.libs.scip.scip import SCIP            # optional - SCIP
            sage: for i in range(100):                            # optional - SCIP
            ...       s = SCIP()                                  # optional - SCIP
            ...       x = s.add_variable(name='x')                # optional - SCIP
            ...       c = s.add_linear_constraint([(x,1.0)],0,0)  # optional - SCIP
            ...       del x                                       # optional - SCIP
            ...       del c                                       # optional - SCIP
            ...       del s                                       # optional - SCIP
        """
        if not self._scip:
            return

        cdef SCIP_RETCODE _status_ = SCIPfree(& self._scip)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d finalizing SCIP."%_status_)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.libs.scip.scip import SCIP # optional - SCIP
            sage: from sage.libs.scip.polynomials import boolean_polynomials as boolean_polynomials_scip # optional - SCIP
            sage: SCIP()                               # optional - SCIP
            SCIP Constraint Integer Program "" ( maximization, 0 variables, 0 constraints )

            sage: SCIP(maximization=False)             # optional - SCIP
            SCIP Constraint Integer Program "" ( minimization, 0 variables, 0 constraints )

            sage: sr = mq.SR(1,2,2,4,gf2=True,polybori=True,allow_zero_inversions=True)
            sage: F,s = sr.polynomial_system()
            sage: scip = SCIP(maximization=False)                # optional - SCIP
            sage: m1 = boolean_polynomials_scip(scip, F); scip # optional - SCIP
            SCIP Constraint Integer Program "" ( minimization, 168 variables, 216 constraints )
        """
        name = SCIPgetProbName(self._scip)
        maximization = SCIPgetObjsense(self._scip) == SCIP_OBJSENSE_MAXIMIZE
        nvars = SCIPgetNVars(self._scip)
        ncons = SCIPgetNConss(self._scip)

        return "SCIP Constraint Integer Program "+ "\"%s\" "%name+"( " + \
             ( "maximization" if maximization else "minimization" ) + \
             ", " + str(nvars) + " variables, " +  str(ncons) + " constraints )"

# #########################################
# MIP Backend Interface
# #########################################

    cpdef int add_variable(self, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, name=None):
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default, the
        variable is both positive and real.

        INPUT:

        - ``lower_bound`` - the lower bound of the variable (default: 0.0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default:
          ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of this variable in the objective
          function (default: 0.0)

        - ``name`` - an optional name for the newly added variable, this
          function picks a new name of the form "x%d" if no name is given
          (default: ``None``).

        OUTPUT: The index of the newly created variable

        EXAMPLE::

            sage: from sage.libs.scip.scip import SCIP # optional - SCIP
            sage: scip = SCIP()                        # optional - SCIP
            sage: scip.add_variable()                  # optional - SCIP
            0
            sage: scip.add_variable()                  # optional - SCIP
            1
            sage: x = scip.add_variable(upper_bound=2.0, integer=True)  # optional - SCIP
            sage: scip.is_variable_integer(x)          # optional - SCIP
            True
            sage: scip.col_bounds(x)                   # optional - SCIP
            (0.0, 2.0)

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")      # optional - SCIP
            sage: p.ncols()                            # optional - SCIP
            0
            sage: p.add_variable()                     # optional - SCIP
            0
            sage: p.ncols()                            # optional - SCIP
            1

        The returned variable index can be used to get access to the underlying
        SCIP Variable object::

            sage: from sage.libs.scip.scip import SCIP # optional - SCIP
            sage: scip = SCIP()                        # optional - SCIP
            sage: x = scip.add_variable()              # optional - SCIP
            sage: scip.variable(x)                     # optional - SCIP
            x0
        """
        cdef SCIP_RETCODE _status_
        if SCIPisTransformed((<SCIP>self)._scip):
            _status_ = SCIPfreeTransform((<SCIP>self)._scip)
            if _status_ != SCIP_OKAY:
                raise MIPSolverException("Problem is already transformed, but could not revert this.")

        if name is None:
            name = "x%d"%(SCIPgetNVars((<SCIP>self)._scip))

        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if  vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        if binary:
            vtype = SCIP.VARTYPE_BINARY
        elif continuous:
            vtype = SCIP.VARTYPE_CONTINUOUS
        elif integer:
            vtype = SCIP.VARTYPE_INTEGER

        v = Variable(self, name, lower_bound, upper_bound, obj, vtype)

        (<SCIP>self)._variables.append(v)
        return SCIPgetNVars((<SCIP>self)._scip) - 1

    cpdef int add_variables(self, int number, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, names=None) except -1:
        """
        Add ``number`` new variables.

        This amounts to adding new columns to the matrix. By default, the
        variables are both positive and real.

        INPUT:

        - ``n`` - the number of new variables (must be > 0)

        - ``lower_bound`` - the lower bound of the variable (default: 0.0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default:
          ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of all variables in the objective
          function (default: 0.0)

        - ``names`` - optional list of names (default: ``None``)

        OUTPUT: The index of the variable created last.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")    # optional - SCIP
            sage: p.ncols()                          # optional - SCIP
            0
            sage: p.add_variables(5)                 # optional - SCIP, indirect doctest
            4
            sage: p.ncols()                          # optional - SCIP
            5
        """
        ncols = SCIPgetNVars((<SCIP>self)._scip)

        if names is None:
            names = [None for _ in range(number)]
        for i in range(number):
            name = names[i]
            (<SCIP>self).add_variable(lower_bound=lower_bound,
                                      upper_bound=upper_bound,
                                      obj=obj,
                                      binary=binary,
                                      integer=integer,
                                      continuous=continuous,
                                      name=name)
        return SCIPgetNVars((<SCIP>self)._scip) - 1

    cpdef set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer) :

            *  1 Integer
            *  0 Binary
            * -1 Real

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")   # optional - SCIP
            sage: p.ncols()                         # optional - SCIP
            0
            sage: p.add_variable()                  # optional - SCIP
            0
            sage: p.set_variable_type(0,1)          # optional - SCIP
            sage: p.is_variable_integer(0)          # optional - SCIP
            True

            sage: p.set_variable_type(0,2)          # optional - SCIP
            Traceback (most recent call last):
            ...
            ValueError: Variable type '2' unknown to SCIP backend.
        """
        cdef Variable var = self._variables[variable]
        cdef SCIP_RETCODE _status_
        cdef bint infeasible

        if vtype == 1:
            _status_ = SCIPchgVarType(self._scip, var._var, SCIP_VARTYPE_INTEGER, &infeasible)
        elif vtype == 0:
            _status_ = SCIPchgVarType(self._scip, var._var, SCIP_VARTYPE_BINARY, &infeasible)
            # We must enforce the boundaries here
            lb, ub = var.bounds()
            if lb is None or lb != 1:
                lb = 0
            self.variable_lower_bound(variable, lb)
            if ub is None or ub != 0:
                ub = 1
            self.variable_upper_bound(variable, ub)
        elif vtype == -1:
            _status_ = SCIPchgVarType(self._scip, var._var, SCIP_VARTYPE_CONTINUOUS, &infeasible)
        else:
            raise ValueError("Variable type '%s' unknown to SCIP backend."%vtype)

        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d setting SCIP variable type."%_status_)

        if infeasible:
            if get_verbose() > 0:
                print "Problem became infeasible due to variable type change."

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")  # optional - SCIP
            sage: p.is_maximization()              # optional - SCIP
            True
            sage: p.set_sense(-1)                  # optional - SCIP
            sage: p.is_maximization()              # optional - SCIP
            False
        """
        cdef SCIP_RETCODE _status_
        if sense == 1:
            _status_ = SCIPsetObjsense(self._scip, SCIP_OBJSENSE_MAXIMIZE)
        elif sense == -1:
            _status_ = SCIPsetObjsense(self._scip, SCIP_OBJSENSE_MINIMIZE)
        else:
            raise ValueError("Sense '%s' unknown to SCIP backend."%sense)

        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d setting the objective sense."%_status_)

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient or ``None`` for
          reading (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")  # optional - SCIP
            sage: p.add_variable()                 # optional - SCIP
            0
            sage: p.objective_coefficient(0)   # optional - SCIP
            0.0
            sage: p.objective_coefficient(0,2) # optional - SCIP
            sage: p.objective_coefficient(0)   # optional - SCIP
            2.0
        """
        cdef SCIP_RETCODE _status_
        cdef Variable var = self._variables[variable]

        if SCIPisTransformed(self._scip):
            _status_ = SCIPfreeTransform(self._scip)
            if _status_ != SCIP_OKAY:
                raise MIPSolverException("Problem is already transformed, but could not revert this.")

        if coeff is not None:
            _status_ = SCIPchgVarObj(self._scip, var._var, coeff)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d setting the objective value of variable %d."%(_status_,variable))

        else:
            return var._var.obj

    cpdef set_objective(self, list coeff, d=0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                     # optional - SCIP
            sage: p.add_variables(5)                                  # optional - SCIP
            4
            sage: p.set_objective([1, 1, 2, 1, 3])                    # optional - SCIP
            sage: [p.objective_coefficient(x) for x in range(5)]  # optional - SCIP
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """
        if not coeff:
            for i in range(len(self._variables)):
                self.objective_coefficient(i, 0.0)

        for i,c in enumerate(coeff):
            self.objective_coefficient(i,c)

    cpdef set_verbosity(self, int level):
        """
        Sets the log (verbosity) level.

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        .. note::

            To get the standard SCIP verbosity, choose ``level=3``.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")   # optional - SCIP
            sage: p.set_verbosity(2)                # optional - SCIP
        """
        # the standard interface only supports levels up to three.
        if level == 3:
            level = 4
        self['display/verblevel'] = level

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None,
                                initial=True, separate=True, enforce=True, check=True, propagate=True,
                                local=False, modifiable=False, dynamic=False, removable=False, stickingatnode=False):
        r"""
        Add a linear constraint.

          `lower_bound <= \sum_{v,c \in coefficients} c*v <= upper_bound`

        INPUT:

        - ``coefficients`` an iterable with ``(c,v)`` pairs where ``c`` is a
          variable index (integer) and ``v`` is a value (real value).

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``name`` - an optional name for this row (default: ``None``)

        - ``initial`` - should the LP relaxation of constraint be in the initial
          LP? Usually set to ``True``. Set to `False`` for 'lazy
          constraints'. (default: ``True``, SCIP specific)

        - ``separate`` - should the constraint be separated during LP
          processing? (default: ``True``, SCIP specific)

        - ``enforce`` - should the constraint be enforced during node
          processing? ``True`` for model constraints, ``False`` for additional,
          redundant constraints.  (default: ``True``, SCIP specific)

        - ``check`` - should the constraint be checked for feasibility? ``True``
          for model constraints, ``False`` for additional, redundant
          constraints.  (default: ``True``, SCIP specific)

        - ``propagate`` - should the constraint be propagated during node
          processing? (default: `True``, SCIP specific)

        - ``local`` - is constraint only valid locally? Usually set to
          ``False``. Has to be set to ``True``, e.g., for branching constraints.
          (default: ``False``, SCIP specific)

        - ``modifiable`` - is constraint modifiable (subject to column
          generation)? Usually set to ``False``. In column generation
          applications, set to ``True`` if pricing adds coefficients to this
          constraint.  (default: ``False``, SCIP specific)

        - ``dynamic`` - Is constraint subject to aging? Usually set to
          ``False``. Set to ``True`` for own cuts which are seperated as
          constraints.  (default: ``False``, SCIP specific)

        - ``removable`` - should the relaxation be removed from the LP due to
          aging or cleanup? Usually set to ``False``. Set to ``True`` for 'lazy
          constraints' and 'user cuts'. (default: ``False``, SCIP specific)

        - ``stickingatnode`` - should the constraint always be kept at the node
          where it was added, even if it may be moved to a more global node?
          Usually set to ``False``. Set to ``True`` to for constraints that
          represent node data. (default: ``False``, SCIP specific)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                         # optional - SCIP
            sage: p.add_variables(5)                                      # optional - SCIP
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 0, 2) # optional - SCIP
            sage: p.row(0)                                                # optional - SCIP
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                                         # optional - SCIP
            (0.0, 2.0)
            sage: p.constraint(-1)                                        # optional - SCIP
            0.000000 <= 1.000000*x1 + 2.000000*x2 + 3.000000*x3 + 4.000000*x4 <= 2.000000
        """
        cdef LinearConstraint c
        V = (<SCIP>self)._variables

        coefficients = [(V[i],coeff) for i,coeff in coefficients]

        if name is None:
            name = ""

        c = LinearConstraint(self, coefficients, name, lower_bound, upper_bound,
                             initial, separate, enforce, check,
                             propagate, local, modifiable, dynamic, removable,
                             stickingatnode)
        self.add_scip_constraint(c)

    def solve(self, known_solutions=None):
        """
        Solve this instance.

        .. note::

            This method raises ``MIPSolverException`` exceptions when the
            solution can not be computed for any reason (none exists, or the LP
            solver was not able to find it, etc...)

        EXAMPLES::

            sage: from sage.libs.scip.scip import SCIP                                                   # optional - SCIP
            sage: from sage.libs.scip.polynomials import boolean_polynomials as boolean_polynomials_scip # optional - SCIP

            sage: B.<a,b,c> = BooleanPolynomialRing()
            sage: F = Sequence([a*b + c + 1, c+1]); F
            [a*b + c + 1, c + 1]

            sage: scip = SCIP(name='bool')                                                               # optional - SCIP
            sage: m1 = boolean_polynomials_scip(scip, F); scip                                           # optional - SCIP
            SCIP Constraint Integer Program "bool" ( maximization, 4 variables, 3 constraints )

            sage: scip = SCIP(name='bool')                                                               # optional - SCIP
            sage: m1 = boolean_polynomials_scip(scip, F, use_xor=False); scip                            # optional - SCIP
            SCIP Constraint Integer Program "bool" ( maximization, 6 variables, 3 constraints )

            sage: scip.solve()                                                                       # optional - SCIP
        """
        cdef SCIP_SOL  *sol
        cdef SCIP_SOL **sols
        cdef double obj
        cdef bint success

        verblevel = self['display/verblevel']
        if verblevel == 0:
            SCIPsetMessagehdlrQuiet((<SCIP>self)._scip, True)

        cdef SCIP_RETCODE _status_ = SCIPsolve((<SCIP>self)._scip)

        if verblevel == 0:
            SCIPsetMessagehdlrQuiet((<SCIP>self)._scip, False)

        if _status_ != SCIP_OKAY:
            raise MIPSolverException("Error %d solving SCIP instance."%_status_)

        cdef int n = SCIPgetNSols((<SCIP>self)._scip)

        if n == 0:
            raise MIPSolverException("SCIP: No solution was found.")

    cpdef get_objective_value(self):
        """
        Return the value of the objective function for the best solution found
        so far.

        .. note::

           Raises an `MIPSolverException` if ``solve`` has not been called yet.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                # optional - SCIP
            sage: p.add_variables(2)                             # optional - SCIP
            1
            sage: p.add_linear_constraint([(0,1), (1,2)], 0, 3)  # optional - SCIP
            sage: p.set_objective([2, 5])                        # optional - SCIP
            sage: p.solve()                                      # optional - SCIP

            sage: p.get_objective_value()                        # optional - SCIP
            7.5
            sage: p.get_variable_value(0)                        # optional - SCIP
            0.0
            sage: p.get_variable_value(1)                        # optional - SCIP
            1.5
        """
        cdef int n = SCIPgetNSols(self._scip)
        if n == 0:
            raise MIPSolverException("SCIP: No solution was found so far.")

        cdef SCIP_SOL  *sol = SCIPgetBestSol(self._scip)
        return SCIPgetSolOrigObj(self._scip, sol)

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. note::

           Raises a ``MIPSolverException`` if ``solve`` has not been called yet.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")         # optional - SCIP
            sage: p.add_variables(2)                      # optional - SCIP
            1
            sage: p.add_linear_constraint([(0, 1), (1, 2)], 0, 3) # optional - SCIP
            sage: p.set_objective([2, 5])                 # optional - SCIP
            sage: p.solve()                               # optional - SCIP

            sage: p.get_objective_value()                 # optional - SCIP
            7.5
            sage: p.get_variable_value(0)                 # optional - SCIP
            0.0
            sage: p.get_variable_value(1)                 # optional - SCIP
            1.5
        """
        cdef int n = SCIPgetNSols(self._scip)
        if n == 0:
            raise MIPSolverException("SCIP: No solution was found so far.")

        cdef SCIP_SOL *sol = SCIPgetBestSol(self._scip)
        cdef Variable var = self._variables[variable]
        return SCIPgetSolVal(self._scip, sol, var._var)

    cpdef int ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")  # optional - SCIP
            sage: p.ncols()                        # optional - SCIP
            0
            sage: p.add_variables(2)               # optional - SCIP
            1
            sage: p.ncols()                        # optional - SCIP
            2
        """
        return SCIPgetNVars(self._scip)

    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP") # optional - SCIP
            sage: p.nrows()                       # optional - SCIP
            0
            sage: p.add_variables(2)              # optional - SCIP
            1
            sage: p.add_linear_constraint([(0, 1), (1, 2)], 0, 3) # optional - SCIP
            sage: p.nrows()                         # optional - SCIP
            1
        """
        return SCIPgetNConss(self._scip)

    cpdef bint is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP") # optional - SCIP
            sage: p.is_maximization()             # optional - SCIP
            True
            sage: p.set_sense(-1)                 # optional - SCIP
            sage: p.is_maximization()             # optional - SCIP
            False
        """
        return SCIPgetObjsense(self._scip) == SCIP_OBJSENSE_MAXIMIZE

    cpdef problem_name(self, char *name = NULL):
        """
        Return or defines the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")               # optional - SCIP
            sage: p.problem_name("There once was a french fry") # optional - SCIP
            sage: print p.problem_name()                        # optional - SCIP
            There once was a french fry
        """
        cdef SCIP_RETCODE _status_

        if name is NULL:
            name = SCIPgetProbName(self._scip)
            return name
        else:
            _status_ = SCIPsetProbName(self._scip, name)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d setting SCIP problem name."%_status_)

    cpdef write_lp(self, char * name):
        """
        Write the problem to a .lp file

        INPUT:

        - ``name`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")         # optional - SCIP
            sage: p.add_variables(2)                      # optional - SCIP
            1
            sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3) # optional - SCIP
            sage: p.set_objective([2, 5])                 # optional - SCIP
            sage: p.write_lp(SAGE_TMP+"/lp_problem.lp")   # optional - SCIP
        """
        if not name.endswith(".lp"):
            _name = name + ".lp"
        else:
            _name = name
        SCIPwriteOrigProblem(self._scip, _name, NULL, False)

    cpdef write_mps(self, char * name, int modern):
        """
        Write the problem to a .mps file

        INPUT:

        - ``name`` -- the filename (string)
        - ``modern`` (ignored)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")         # optional - SCIP
            sage: p.add_variables(2)                      # optional - SCIP
            1
            sage: p.add_linear_constraint([(0, -1), (1, 0.5)], None, 3) # optional - SCIP
            sage: p.set_objective([2, 5])                 # optional - SCIP
            sage: p.write_mps(SAGE_TMP+"/lp_problem.mps", True) # optional - SCIP
            WARNING: At least one name of a constraint is empty, so file will be written with generic names.
            WARNING: write original problem with generic variable and constraint names
        """
        if not name.endswith(".mps"):
            _name = name + ".mps"
        else:
            _name = name
        SCIPwriteOrigProblem(self._scip, _name, NULL, False)

    cpdef row(self, int i):
        """
        Return a row.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_scip_constraint`` method.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")  # optional - SCIP
            sage: p.add_variables(5)                                # optional - SCIP
            4
            sage: p.add_linear_constraint( [(i,i) for i in range(5)], 2, 2) # optional - SCIP
            sage: p.row(0)                                          # optional - SCIP
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])

            sage: p.row_bounds(0)                                   # optional - SCIP
            (2.0, 2.0)
            sage: p.add_and_constraint(0, range(1,5)) # optional - SCIP
            sage: p.row(1)                            # optional - SCIP
            Traceback (most recent call last):
            ...
            TypeError: Expected linear constraint but got type `<type 'sage.libs.scip.constraint.ANDConstraint'>' instead.
        """
        cdef Constraint c = self._constraints[i]
        if not PY_TYPE_CHECK(c, LinearConstraint):
            raise TypeError("Expected linear constraint but got type `%s' instead."%type(c))
        return [v.index() for v in c.variables()], c.coefficients()

    cpdef row_bounds(self, int index):
        """
        Return the bounds of a specific constraint.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be
        ``None`` if the constraint is not bounded in the corresponding
        direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")            # optional - SCIP
            sage: p.add_variables(5)                         # optional - SCIP
            4
            sage: p.add_linear_constraint([(i,i) for i in range(5)], 2, 2) # optional - SCIP
            sage: p.row(0)                                   # optional - SCIP
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                            # optional - SCIP
            (2.0, 2.0)
        """
        cdef Constraint c = self._constraints[index]
        if not PY_TYPE_CHECK(c, LinearConstraint):
            raise TypeError("Expected linear constraint but got type `%s' instead."%type(c))
        return c.bounds()

    cpdef col_bounds(self, int index):
        """
        Return the bounds of a specific variable.

        INPUT:

        - ``index`` (integer) -- the variable's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be
        ``None`` if the variable is not bounded in the corresponding
        direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")  # optional - SCIP
            sage: p.add_variable()                 # optional - SCIP
            0
            sage: p.col_bounds(0)                  # optional - SCIP
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)             # optional - SCIP
            sage: p.col_bounds(0)                  # optional - SCIP
            (0.0, 5.0)
        """
        cdef Variable var = self._variables[index]
        return var.bounds()

    cpdef bint is_variable_binary(self, int index):
        """
        Return whether the given variable is binary.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")  # optional - SCIP
            sage: p.ncols()                                       # optional - SCIP
            0
            sage: p.add_variable()                                 # optional - SCIP
            0
            sage: p.set_variable_type(0,0)                         # optional - SCIP
            sage: p.is_variable_binary(0)                          # optional - SCIP
            True

        """
        cdef Variable var = self._variables[index]
        return var.is_binary()

    cpdef bint is_variable_integer(self, int index):
        """
        Return whether the given variable is integral.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")  # optional - SCIP
            sage: p.ncols()                        # optional - SCIP
            0
            sage: p.add_variable()                 # optional - SCIP
            0
            sage: p.set_variable_type(0,1)         # optional - SCIP
            sage: p.is_variable_integer(0)         # optional - SCIP
            True
        """
        cdef Variable var = self._variables[index]
        return var.is_integer()

    cpdef bint is_variable_continuous(self, int index):
        """
        Return whether the given variable is of continuous.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")  # optional - SCIP
            sage: p.ncols()                                        # optional - SCIP
            0
            sage: p.add_variable()                                 # optional - SCIP
            0
            sage: p.is_variable_continuous(0)                      # optional - SCIP
            True
            sage: p.set_variable_type(0,1)                         # optional - SCIP
            sage: p.is_variable_continuous(0)                      # optional - SCIP
            False

        """
        cdef Variable var = self._variables[index]
        return var.is_continuous()

    cpdef row_name(self, int index, char * name = NULL):
        """
        Return or define the ``index`` th row name.

        INPUT:

        - ``index`` (integer) -- the row's id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")        # optional - SCIP
            sage: x = p.add_variable()                   # optional - SCIP
            sage: p.add_linear_constraint([(0,1.0)],0,0) # optional - SCIP
            sage: p.row_name(0, "Some constraint")       # optional - SCIP
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if index > SCIPgetNConss((<SCIP>self)._scip):
            raise IndexError("Constraint index %d > %d."%(index, SCIPgetNConss((<SCIP>self)._scip)))
        cdef SCIP_CONS **conss = SCIPgetConss((<SCIP>self)._scip)
        if name is NULL:
            name = SCIPconsGetName(conss[index])
            return name
        else:
            raise NotImplementedError()

    cpdef col_name(self, int index):
        """
        Returns the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")          # optional - SCIP
            sage: p.add_variable(name='I am a variable')   # optional - SCIP
            0
            sage: print p.col_name(0)                      # optional - SCIP
            <BLANKLINE>
        """
        return ""

    cpdef variable_upper_bound(self, int index, value = False):
        """
        Return or define the upper bound on a variable.

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` if the variable has no
          upper bound. When set to ``False`` (default), the method
          returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")  # optional - SCIP
            sage: p.add_variable()                 # optional - SCIP
            0
            sage: p.col_bounds(0)                  # optional - SCIP
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)     # optional - SCIP
            sage: p.col_bounds(0)                  # optional - SCIP
            (0.0, 5.0)
        """
        cdef SCIP_RETCODE _status_
        cdef Variable var = self._variables[index]
        if value is False:
            return var.bounds()[1]
        else:
            if value is None:
                value = SCIPinfinity(self._scip)
            _status_ = SCIPchgVarUb(self._scip, var._var, value)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d when setting SCIP variable upper bound."%_status_)


    cpdef  variable_lower_bound(self, int index, value = False):
        """
        Return or define the lower bound on a variable.

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")  # optional - SCIP
            sage: p.add_variable()                 # optional - SCIP
            0
            sage: p.col_bounds(0)                  # optional - SCIP
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)     # optional - SCIP
            sage: p.col_bounds(0)                  # optional - SCIP
            (5.0, None)
        """
        cdef Variable var = self._variables[index]
        if value is False:
            return var.bounds()[0]
        else:
            if value is None:
                value = -SCIPinfinity(self._scip)
            _status_ = SCIPchgVarLb(self._scip, var._var, value)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d when setting SCIP variable upper bound."%_status_)

    ####################################################
    #
    #       SCIP specific functions
    #
    ###################################################

    def variable(self, int i):
        """
        Return the ``i``-th variable.

        OUTPUT:An object representing a variable.

        .. note::

            This function differs from :meth:`column`: it returns a
            SCIP Variable object.

        EXAMPLE::

            sage: from sage.libs.scip.scip import SCIP      # optional - SCIP
            sage: scip = SCIP(maximization=False)           # optional - SCIP
            sage: z = scip.add_variable(name='z', integer=True)  # optional - SCIP
            sage: scip.variable(-1)                         # optional - SCIP
            z
        """
        return self._variables[i]

    def constraint(self, int i):
        """
        Return the ``i``-th contraint.

        OUTPUT:

            An object representing a (possibly non-linear) constraint.

        .. note::

            This function differs from :meth:`row`: it returns a SCIP
            Constraint object and it supports non-linear constraints.

        EXAMPLE::

            sage: from sage.libs.scip.scip import SCIP                 # optional - SCIP
            sage: scip = SCIP(maximization=False)                      # optional - SCIP
            sage: x = scip.add_variable(name='x', integer=True)        # optional - SCIP
            sage: y = scip.add_variable(name='y', integer=True)        # optional - SCIP
            sage: z = scip.add_variable(name='z', integer=True)        # optional - SCIP
            sage: scip.add_quadratic_constraint([(x,-2.0)],[(x,y,1.0)],0.0, 1.0) # optional - SCIP
            sage: scip.constraint(-1)                                  # optional - SCIP
            0.000000 <= 1.000000*x*y + -2.000000*x <= 1.000000
        
        """
        return self._constraints[i]

    def get_all_solutions(self):
        """
        Returns a list of lists, containing all solutions found by the
        solver so far and their corresponding objective values.

        OUTPUT:

            ``[(s0,obj0),(s1,obj1),...]`` where ``si`` is a dictionary
            mapping variable indices to values and ``obji`` is the
            value of the objective function corresponding to ``si``.

        .. note::

            This function does not return all solutions for the
            problem, but all solutions found so far by the solver.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP") # optional - SCIP
            sage: p.add_variables(2)              # optional - SCIP
            1
            sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3) # optional - SCIP
            sage: p.set_objective([2, 5])         # optional - SCIP
            sage: _ = p.solve()                   # optional - SCIP
            sage: p.get_all_solutions()           # optional - SCIP
            [({0: 0.0, 1: 1.5}, 7.5), ({0: 0.0, 1: 0.0}, 0.0)]
        """
        cdef int n = SCIPgetNSols(self._scip)
        if n == 0:
            raise MIPSolverException("SCIP: No solution was found so far.")
        cdef SCIP_SOL **sols = SCIPgetSols(self._scip)
        S = []
        for i in range(n):
            s = {}
            for j,var in enumerate(self._variables):
                val = SCIPgetSolVal(self._scip, sols[i], (<Variable>var)._var)
                s[j] = val
            obj = SCIPgetSolOrigObj(self._scip, sols[i])
            S.append( (s, obj) )
        return S

    def add_quadratic_constraint(self, linterms, quadterms, lower_bound, upper_bound, name=None,
                                 initial=True, separate=True, enforce=True, check=True, propagate=True,
                                 local=False, modifiable=False, dynamic=False, removable=False, stickingatnode=False):
        r"""
        Add a quadratic constraint.

        That is, a constraint of the form

          `lower_bound <= \sum_{v,c \in linterms} c*v + \sum_{v1,v2,c \in quadterms} c*v1*v2 <= upper_bound`

        INPUT:

        - ``linterms`` - list of ``(variable,value)`` pairs

        - ``quadterms`` - list of ``(variable,variable,value)`` tuples

        - ``lower_bound`` - a lower bound, either a real value or
          ``None``

        - ``upper_bound`` - an upper bound, either a real value or
          ``None``

        - ``name`` - an optional name for this row (default: ``None``)

        For the remain parameters see documentation for
        :meth:`add_linear_constraint`.

        EXAMPLE::

            sage: from sage.libs.scip.scip import SCIP # optional - SCIP
            sage: scip = SCIP(maximization=False)      # optional - SCIP
            sage: x = scip.add_variable(name='x', integer=True) # optional - SCIP
            sage: y = scip.add_variable(name='y', integer=True) # optional - SCIP
            sage: z = scip.add_variable(name='z', integer=True) # optional - SCIP
            sage: scip.add_quadratic_constraint([],[(x,y,1.0)], 1.0, 1.0) # optional - SCIP
            sage: _ = scip.solve()                              # optional - SCIP
            sage: [scip.get_variable_value(i) for i in (x,y,z)] # optional - SCIP
            [1.0, 1.0, -0.0]
        """
        if name is None:
            name = ''
        V = self._variables
        linterms = [(V[v],coeff) for (v,coeff) in linterms]
        quadterms = [(V[v],V[w],coeff) for (v,w,coeff) in quadterms]

        cdef QuadraticConstraint c =  QuadraticConstraint(self, linterms, quadterms, name, lower_bound, upper_bound,
                                                          initial, separate, enforce, check,
                                                          propagate, local, modifiable, dynamic, removable,
                                                          stickingatnode)
        self.add_scip_constraint(c)

    def add_and_constraint(self, res, variables, name=None,
                              initial=True, separate=True, enforce=True, check=True, propagate=True,
                              local=False, modifiable=False, dynamic=False, removable=False, stickingatnode=False):
        """
        Add an AND constraint.

        If ``variables`` is ``[v_1,v_2,...,v_n]`` then this constraint
        restricts ``res`` to ``res = v_1 AND v_2 AND ... AND v_n``.

        INPUT:

        - ``res`` - a variable

        - ``variables`` - a list of variables

        For the remain parameters see documentation for
        :meth:`add_linear_constraint`.

        EXAMPLE::

            sage: from sage.libs.scip.scip import SCIP         # optional - SCIP
            sage: scip = SCIP(maximization=False)              # optional - SCIP
            sage: x = scip.add_variable(name='x', binary=True) # optional - SCIP
            sage: y = scip.add_variable(name='y', binary=True) # optional - SCIP
            sage: z = scip.add_variable(name='z', binary=True) # optional - SCIP
            sage: scip.add_and_constraint(x,[y,z])             # optional - SCIP
            sage: scip.constraint(-1)                          # optional - SCIP
            x = y & z
        """
        if name is None:
            name = ''
        res = self._variables[res]
        variables = [self._variables[i] for i in variables]
        cdef ANDConstraint c=  ANDConstraint(self, res, variables, name,
                                             initial, separate, enforce, check,
                                             propagate, local, modifiable, dynamic, removable,
                                             stickingatnode)
        self.add_scip_constraint(c)

    def add_xor_constraint(self, rhs, variables, name=None,
                              initial=True, separate=True, enforce=True, check=True, propagate=True,
                              local=False, modifiable=False, dynamic=False, removable=False, stickingatnode=False):
        """
        Add a XOR constraint.

        If ``variables`` is ``[v_1,v_2,...,v_n]`` then this expresses
        the constraint ``rhs = v_1 XOR v_2 XOR ... XOR v_n``.

        INPUT:

        - ``rhs`` - a bit

        - ``variables`` - a list of variables

        For the remain parameters see documentation for
        :meth:`add_linear_constraint`.

        EXAMPLE::

            sage: from sage.libs.scip.scip import SCIP         # optional - SCIP
            sage: scip = SCIP(maximization=False)              # optional - SCIP
            sage: x = scip.add_variable(name='x', binary=True) # optional - SCIP
            sage: y = scip.add_variable(name='y', binary=True) # optional - SCIP
            sage: z = scip.add_variable(name='z', binary=True) # optional - SCIP
            sage: scip.add_xor_constraint(1, [x,y,z])          # optional - SCIP
            sage: scip.constraint(-1)                          # optional - SCIP
            1 = x ^ y ^ z
        """
        if name is None:
            name = ''
        variables = [self._variables[i] for i in variables]
        cdef XORConstraint c = XORConstraint(self, rhs, variables, name,
                                             initial, separate, enforce, check,
                                             propagate, local, modifiable, dynamic, removable,
                                             stickingatnode)
        self.add_scip_constraint(c)

    def add_logicor_constraint(self, variables, name=None,
                               initial=True, separate=True, enforce=True, check=True, propagate=True,
                               local=False, modifiable=False, dynamic=False, removable=False, stickingatnode=False):
        """
        Add a logical OR constraint.

        If ``variables`` is ``[v_1,v_2,...,v_n]`` then this constraint
        modells that at at least one of them must be non-zero.

        INPUT:

        - ``variables`` - a list of variables

        For the remain parameters see documentation for
        :meth:`add_linear_constraint`.

        EXAMPLE::

            sage: from sage.libs.scip.scip import SCIP         # optional - SCIP
            sage: scip = SCIP(maximization=False)              # optional - SCIP
            sage: x = scip.add_variable(name='x', binary=True) # optional - SCIP
            sage: y = scip.add_variable(name='y', binary=True) # optional - SCIP
            sage: z = scip.add_variable(name='z', binary=True) # optional - SCIP
            sage: scip.add_logicor_constraint([x,y,z])         # optional - SCIP
            sage: scip.constraint(-1)                          # optional - SCIP
            x or y or z
            sage: _ = scip.solve()                             # optional - SCIP
        """
        if name is None:
            name = ''
        variables = [self._variables[i] for i in variables]
        cdef LogicORConstraint c = LogicORConstraint(self, variables, name,
                                                     initial, separate, enforce, check,
                                                     propagate, local, modifiable, dynamic, removable,
                                                     stickingatnode)
        self.add_scip_constraint(c)

    def add_or_constraint(self, res, variables, name=None,
                          initial=True, separate=True, enforce=True, check=True, propagate=True,
                          local=False, modifiable=False, dynamic=False, removable=False, stickingatnode=False):
        """
        Add an OR constraint to this SCIP instance.

        If ``variables`` is ``[v_1,v_2,...,v_n]`` then this constraint
        modells ``res = v_1 OR v_2 OR ... OR v_n``

        For the remain parameters see documentation for
        :meth:`add_linear_constraint`.

        INPUT:

        - ``res`` - a variable

        - ``variables`` - a list of variables

        EXAMPLE::

            sage: from sage.libs.scip.scip import SCIP         # optional - SCIP
            sage: scip = SCIP(maximization=False)              # optional - SCIP
            sage: x = scip.add_variable(name='x', binary=True) # optional - SCIP
            sage: y = scip.add_variable(name='y', binary=True) # optional - SCIP
            sage: z = scip.add_variable(name='z', binary=True) # optional - SCIP
            sage: scip.add_or_constraint(x,[y,z])              # optional - SCIP
            sage: scip.constraint(-1)                          # optional - SCIP
            x = y | z
        """
        if name is None:
            name = ''
        res = self._variables[res]
        variables = [self._variables[i] for i in variables]
        cdef ORConstraint c=  ORConstraint(self, res, variables, name,
                                           initial, separate, enforce, check,
                                           propagate, local, modifiable, dynamic, removable,
                                           stickingatnode)
        self.add_scip_constraint(c)

    def add_setppc_constraint(self, variables, type, name=None,
                              initial=True, separate=True, enforce=True, check=True, propagate=True,
                              local=False, modifiable=False, dynamic=False, removable=False, stickingatnode=False):
        """
        Add a set partition, set packing or set covering constraint.

        INPUT:

        - ``variables`` - a list of variables

        For the remain parameters see documentation for
        :meth:`add_linear_constraint`.

        EXAMPLE::

            sage: from sage.libs.scip.scip import SCIP                  # optional - SCIP
            sage: scip = SCIP(maximization=False)                       # optional - SCIP
            sage: a = scip.add_variable(name='a', binary=True, obj=1.0) # optional - SCIP
            sage: b = scip.add_variable(name='b', binary=True, obj=1.0) # optional - SCIP
            sage: c = scip.add_variable(name='c', binary=True, obj=1.0) # optional - SCIP
            sage: scip.add_setppc_constraint([a,b,c],'partition')       # optional - SCIP
            sage: scip.constraint(-1)                                   # optional - SCIP
            sum(a, b, c) == 1
            sage: scip.solve()                                          # optional - SCIP
            sage: [scip.get_variable_value(i) for i in (a,b,c)]         # optional - SCIP
            [1.0, 0.0, 0.0]

            sage: scip = SCIP(maximization=False)                       # optional - SCIP
            sage: a = scip.add_variable(name='a', binary=True, obj=1.0) # optional - SCIP
            sage: b = scip.add_variable(name='b', binary=True, obj=1.0) # optional - SCIP
            sage: c = scip.add_variable(name='c', binary=True, obj=1.0) # optional - SCIP
            sage: scip.add_setppc_constraint([a,b,c],'packing')         # optional - SCIP
            sage: scip.constraint(-1)                                   # optional - SCIP
            sum(a, b, c) <= 1
            sage: scip.solve()                                          # optional - SCIP
            sage: [scip.get_variable_value(i) for i in (a,b,c)]         # optional - SCIP
            [-0.0, -0.0, -0.0]

            sage: scip = SCIP(maximization=False)                       # optional - SCIP
            sage: a = scip.add_variable(name='a', binary=True, obj=1.0) # optional - SCIP
            sage: b = scip.add_variable(name='b', binary=True, obj=1.0) # optional - SCIP
            sage: c = scip.add_variable(name='c', binary=True, obj=1.0) # optional - SCIP
            sage: scip.add_setppc_constraint([a,b,c],'covering')        # optional - SCIP
            sage: scip.constraint(-1)                                   # optional - SCIP
            sum(a, b, c) >= 1
            sage: scip.solve()                                          # optional - SCIP
            sage: [scip.get_variable_value(i) for i in (a,b,c)]         # optional - SCIP
            [1.0, -0.0, -0.0]
        """
        if name is None:
            name = ''
        variables = [self._variables[i] for i in variables]
        cdef SetPPCConstraint c=  SetPPCConstraint(self, variables, type, name,
                                                   initial, separate, enforce, check,
                                                   propagate, local, modifiable, dynamic, removable,
                                                   stickingatnode)
        self.add_scip_constraint(c)

    def add_scip_constraint(self, Constraint c):
        """
        Adds an existing SCIP constraint object ``o`` to this problem.

        INPUT:

        - ``o`` - SCIP constraint object

        EXAMPLE::

            sage: from sage.libs.scip import SCIP                # optional - SCIP
            sage: from sage.libs.scip.variable import Variable   # optional - SCIP
            sage: from sage.libs.scip.constraint import QuadraticConstraint # optional - SCIP
            sage: scip = SCIP()                                  # optional - SCIP

            sage: x = Variable(scip, "x", 0.0, 1.0, 0.0, scip.VARTYPE_BINARY) # optional - SCIP
            sage: y = Variable(scip, "y", 0.0, 1.0, 0.0, scip.VARTYPE_BINARY) # optional - SCIP
            sage: q = [(x,y,1),(x,x,1)]                          # optional - SCIP
            sage: l = [(x,-3),(y,2)]                             # optional - SCIP
            sage: Q = QuadraticConstraint(scip, l, q, "quad", 0, 0) # optional - SCIP
            sage: scip.add_scip_constraint(Q)                       # optional - SCIP
        """
        cdef SCIP_RETCODE _status_
        if SCIPisTransformed(self._scip):
            _status_ = SCIPfreeTransform(self._scip)
            if _status_ != SCIP_OKAY:
                raise MIPSolverException("Problem is already transformed, but could not revert this.")

        _status_ = SCIPaddCons(self._scip, c._cons)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d adding constraint to SCIP instance."%_status_)
        self._constraints.append(c)

    def add_solution(self, s):
        """
        Add a known solution to this problem.

        This could, for example, be a solution which is sub-optimal
        such that the solver can improve upon it.

        EXAMPLE::

            sage: from sage.libs.scip import SCIP                # optional - SCIP
            sage: scip = SCIP(maximization=False)                # optional - SCIP
            sage: x = scip.add_variable(name='x',lower_bound=-1.0,upper_bound=1.0, obj=1.0) # optional - SCIP
            sage: y = scip.add_variable(name='y',lower_bound=-1.0,upper_bound=1.0, obj=1.0) # optional - SCIP
            sage: scip.add_quadratic_constraint([], [(x,y,1)], lower_bound=1.0, upper_bound=1.0) # optional - SCIP
            sage: scip.add_solution({x:1.0,y:1.0}) # optional - SCIP

            sage: scip.add_solution({x:1.0,y:0.5}) # optional - SCIP
            Traceback (most recent call last):
            ...
            ValueError: Proposed solution not stored since it violates at least one constraint.

            sage: scip.solve()                     # optional - SCIP
            sage: scip.get_all_solutions()         # optional - SCIP
            [({0: -1.0, 1: -1.0}, -2.0), ({0: 1.0, 1: 1.0}, 2.0)]
        """
        cdef SCIP_RETCODE _status_
        cdef SCIP_SOL *sol

        s2 = {}
        for k,v in s.iteritems():
            try:
                k = self._variables[k]
            except TypeError:
                assert(PY_TYPE_CHECK(k, Variable))
            assert((<Variable>k)._scip is self)
            s2[k] = v
        s = s2

        # First you have to build your problem, then you have to
        # transform your problem (SCIP only accepts solutions if it is
        # at least in the transformed stage, see here) via calling
        # SCIPtransformProb().

        if not SCIPisTransformed (self._scip):
            _status_  = SCIPtransformProb(self._scip)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d transforming SCIP Problem."%_status_)

        # Next, you create a new SCIP primal solution
        # by calling SCIPcreateSol()

        _status_ = SCIPcreateSol(self._scip, &sol, NULL)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d creating SCIP solution."%_status_)

        # and set all nonzero values by calling SCIPsetSolVal().

        for k,v in s.iteritems():
            _status_ = SCIPsetSolVal(self._scip, sol, (<Variable>k)._var, float(v))
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d setting SCIP solution value for variable %s."%(_status_,k))

        # After that, you add this solution by calling SCIPtrySol()
        # (the variable success should be true afterwards, if your
        # solution was correct)

        cdef bint stored = False
        cdef bint printreason = False

        if get_verbose() > 1:
            printreason = True

        _status_ = SCIPtrySol(self._scip, sol, printreason, True, True, True, &stored)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d trying SCIP solution."%(_status_))

        # and then release it by calling SCIPsolFree().
        _status = SCIPfreeSol(self._scip, &sol)
        if _status_ != SCIP_OKAY:
            raise RuntimeError("Error %d freeing SCIP solution."%(_status_))

        if not stored:
            raise ValueError("Proposed solution not stored since it violates at least one constraint.")

##########################
# Parameters
##########################

    def read_parameters(self, filename):
        r"""
        Read SCIP parameters from the file ``filename``

        INPUT:

        - ``filename`` - a file with SCIP parameters.

        EXAMPLE::

            sage: fn = tmp_filename()
            sage: fh = open(fn,'w')
            sage: fh.write("display/verblevel=4\n")
            sage: fh.close()

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")   # optional - SCIP
            sage: p['display/verblevel']            # optional - SCIP
            0
            sage: p.read_parameters(fn)             # optional - SCIP
            sage: p['display/verblevel']            # optional - SCIP
            4
        """
        cdef SCIP_RETCODE _status_ = SCIPreadParams(self._scip, filename)
        if _status_ != SCIP_OKAY:
            raise IOError("Error %d reading parameters from file '%s'"%(_status_,filename))

    def write_parameters(self, filename):
        """
        Write current SCIP parameters to file ``filename``

        INPUT:

        - ``filename`` - a writable path + filename

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")   # optional - SCIP
            sage: p['display/verblevel'] = 3        # optional - SCIP
            sage: fn = tmp_filename()
            sage: p.write_parameters(fn)            # optional - SCIP
            sage: p['display/verblevel'] = 2        # optional - SCIP
            sage: p.read_parameters(fn)             # optional - SCIP
            sage: p['display/verblevel']            # optional - SCIP
            3
        """
        cdef SCIP_RETCODE _status_ = SCIPwriteParams(self._scip, filename, 1, 0)
        if _status_ != SCIP_OKAY:
            raise IOError("Error %d writing parameters to file '%s'"%(_status_,filename))

    def reset_parameters(self):
        """
        Reset parameters to default values.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")   # optional - SCIP
            sage: p['display/verblevel'] = 3        # optional - SCIP
            sage: p.reset_parameters()              # optional - SCIP
            sage: p['display/verblevel']            # optional - SCIP
            4
        """
        cdef SCIP_RETCODE _status_ = SCIPresetParams(self._scip)
        if _status_ != SCIP_OKAY:
            raise IOError("Error %d resetting paramters'"%(_status_))

    def __getitem__(self, name):
        """
        Read the value for the parameter ``name``

        INPUT:

        - ``name`` - a parameter for a SCIP solver

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")   # optional - SCIP
            sage: p['branching/relpscost/maxreliable'] # optional - SCIP
            8.0
            sage: p['branching/relpscost/inititer']    # optional - SCIP
            0
            sage: p['separating/cmir/freq']            # optional - SCIP
            0

        .. note::

            Consult the SCIP documentation for supported parameters.
        """
        cdef SCIP_RETCODE _status_
        cdef int i
        cdef SCIP_Param **params = SCIPgetParams(self._scip)
        cdef int n = SCIPgetNParams (self._scip)
        cdef bint v0
        cdef int v1
        cdef long v2
        cdef double v3
        cdef char v4
        cdef char *v5

        for i in range(n):
            if name == params[i].name:
                break
        else:
            raise KeyError("Name '%s' unknown."%name)

        if params[i].paramtype == SCIP_PARAMTYPE_BOOL:
            _status_ = SCIPgetBoolParam (self._scip, name, &v0)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d retrieving parameter."%_status_)
            return v0

        if params[i].paramtype == SCIP_PARAMTYPE_INT:
            _status_ = SCIPgetIntParam(self._scip, name, &v1)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d retrieving parameter."%_status_)
            return v1

        if params[i].paramtype == SCIP_PARAMTYPE_LONGINT:
            _status_ = SCIPgetLongintParam(self._scip, name, &v2)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d retrieving parameter."%_status_)
            return v2

        if params[i].paramtype == SCIP_PARAMTYPE_REAL:
            _status_ = SCIPgetRealParam(self._scip, name, &v3)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d retrieving parameter."%_status_)
            return v3

        if params[i].paramtype == SCIP_PARAMTYPE_CHAR:
            _status_ = SCIPgetCharParam(self._scip, name, &v4)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d retrieving parameter."%_status_)
            return v4

        if params[i].paramtype == SCIP_PARAMTYPE_STRING:
            _status_ = SCIPgetStringParam(self._scip, name, &v5)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d retrieving parameter."%_status_)
            return v5

    def __setitem__(self, name, value):
        """
        Read the value for the parameter ``name``

        INPUT:

        - ``name`` - a parameter for a SCIP solver

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")   # optional - SCIP
            sage: p['branching/relpscost/maxreliable'] = 1 # optional - SCIP
            sage: p['branching/relpscost/inititer'] = 1    # optional - SCIP
            sage: p['separating/cmir/freq']  = -1          # optional - SCIP
            sage: p['branching/relpscost/maxreliable'] # optional - SCIP
            1.0
            sage: p['branching/relpscost/inititer']    # optional - SCIP
            1
            sage: p['separating/cmir/freq']            # optional - SCIP
            -1

        .. note::

            Consult the SCIP documentation for supported parameters.
        """
        cdef SCIP_RETCODE _status_
        cdef int i
        cdef SCIP_Param **params = SCIPgetParams(self._scip)
        cdef int n = SCIPgetNParams (self._scip)

        for i in range(n):
            if name == params[i].name:
                break
        if i == n:
            raise KeyError("Name '%s' unknown."%name)

        if params[i].paramtype == SCIP_PARAMTYPE_BOOL:
            _status_ = SCIPsetBoolParam (self._scip, name, value)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d setting parameter."%_status_)

        if params[i].paramtype == SCIP_PARAMTYPE_INT:
            _status_ = SCIPsetIntParam(self._scip, name, value)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d setting parameter."%_status_)

        if params[i].paramtype == SCIP_PARAMTYPE_LONGINT:
            _status_ = SCIPsetLongintParam(self._scip, name, value)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d setting parameter."%_status_)

        if params[i].paramtype == SCIP_PARAMTYPE_REAL:
            _status_ = SCIPsetRealParam(self._scip, name, value)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d setting parameter."%_status_)

        if params[i].paramtype == SCIP_PARAMTYPE_CHAR:
            _status_ = SCIPsetCharParam(self._scip, name, value)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d setting parameter."%_status_)

        if params[i].paramtype == SCIP_PARAMTYPE_STRING:
            _status_ = SCIPsetStringParam(self._scip, name, value)
            if _status_ != SCIP_OKAY:
                raise RuntimeError("Error %d setting parameter."%_status_)




