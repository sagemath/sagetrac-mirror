"""
SCIP Backend

AUTHORS:

- Nathann Cohen (2010-10): generic backend
- Matthias Koeppe (2017): stubs
- Moritz Firsching (2018-04): rest
"""

#*****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from os import sys
from os.path import splitext
from sage.numerical.mip import MIPSolverException
from sage.numerical.backends.generic_backend import GenericBackend
from pyscipopt import Model


class SCIPBackend(GenericBackend):

    """
    MIP Backend that uses the SCIP solver.

    TESTS:

    General backend testsuite::

        sage: p = MixedIntegerLinearProgram(solver="SCIP")                      # optional - pyscipopt
        sage: TestSuite(p.get_backend()).run(skip="_test_pickling")             # optional - pyscipopt
    """

    def __init__(self, maximization=True):
        """
        Constructor

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="SCIP")                  # optional - pyscipopt
        """
        self.model = Model('')
        if maximization:
            self.set_sense(1)
        else:
            self.set_sense(-1)
        self.model.hideOutput()
        self.obj_constant_term = 0.0

    def _get_model(self):
        """
        Get the model as a pyscipopt Model.

        EXAMPLE::
        sage: from sage.numerical.backends.generic_backend import get_solver    # optional - pyscipopt
        sage: p = get_solver(solver = "SCIP")                                   # optional - pyscipopt
        sage: p._get_model()                                                    # optional - pyscipopt
        <pyscipopt.scip.Model object at ...
        """
        return self.model

    def add_variable(self, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, name=None):
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive, real and the coefficient in the
        objective function is 0.0.

        INPUT:

        - ``lower_bound`` - the lower bound of the variable (default: 0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default: ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` - an optional name for the newly added variable (default: ``None``).

        OUTPUT: The index of the newly created variable

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.ncols()                                                     # optional - pyscipopt
            0
            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.ncols()                                                     # optional - pyscipopt
            1
            sage: p.add_variable(binary=True)                                   # optional - pyscipopt
            1
            sage: p.add_variable(lower_bound=-2.0, integer=True)                # optional - pyscipopt
            2
            sage: p.add_variable(continuous=True, integer=True)                 # optional - pyscipopt
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.add_variable(name='x', obj=1.0)                             # optional - pyscipopt
            3
            sage: p.col_name(3)                                                 # optional - pyscipopt
            u'x'
            sage: p.objective_coefficient(3)                                    # optional - pyscipopt
            1.0
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        if name == None:
            vname=''
        else:
            assert(type(name) in [str, unicode])
            vname = name

        if continuous:
            vtypestr = 'C'
        if binary:
            vtypestr = 'B'
        if integer:
            vtypestr = 'I'


        self.model.addVar(name=vname, vtype=vtypestr, ub=upper_bound, lb=lower_bound, obj=obj, pricedVar=False)
        return self.ncols() - 1

    def set_variable_type(self, variable, vtype):
        """
        Set the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer) :

            *  1  Integer
            *  0  Binary
            * -1  Real

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.ncols()                                                     # optional - pyscipopt
            0
            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.set_variable_type(0,1)                                      # optional - pyscipopt
            sage: p.is_variable_integer(0)                                      # optional - pyscipopt
            True
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        vtypenames = {1: 'I', 0: 'B', -1: 'C'}
        self.model.chgVarType(var=self.model.getVars()[variable], vtype=vtypenames[vtype])

    def set_sense(self, sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.is_maximization()                                           # optional - pyscipopt
            True
            sage: p.set_sense(-1)                                               # optional - pyscipopt
            sage: p.is_maximization()                                           # optional - pyscipopt
            False
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        if sense == 1:
            self.model.setMaximize()
        elif sense == -1:
            self.model.setMinimize()
        else:
            raise AssertionError("sense must be either 1 or -1")

    def objective_coefficient(self, variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient or ``None`` for
          reading (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.objective_coefficient(0)                                    # optional - pyscipopt
            0.0
            sage: p.objective_coefficient(0,2)                                  # optional - pyscipopt
            sage: p.objective_coefficient(0)                                    # optional - pyscipopt
            2.0
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        if coeff is None:
            return self.model.getVars()[variable].getObj()
        else:
            objexpr = self.model.getObjective()
            var = self.model.getVars()[variable]
            from pyscipopt.scip import quicksum
            linfun = quicksum([e*c for e, c in objexpr.terms.iteritems() if e != var]) + var*coeff
            self.model.setObjective(linfun, sense=self.model.getObjectiveSense())

    def problem_name(self, name=None):
        """
        Return or define the problem's name

        INPUT:

        - ``name`` -- the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.problem_name("Nomen est omen")                              # optional - pyscipopt
            sage: p.problem_name()                                              # optional - pyscipopt
            u'Nomen est omen'
        """
        if name is None:
            return self.model.getProbName()
        else:
            self.model.setProbName(name)

    def set_objective(self, coeff, d=0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` - a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        - ``d`` (double) -- the constant term in the linear function (set to `0` by default)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variables(5)                                            # optional - pyscipopt
            4
            sage: p.set_objective([1, 1, 2, 1, 3])                              # optional - pyscipopt
            sage: map(lambda x :p.objective_coefficient(x), range(5))           # optional - pyscipopt
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        from pyscipopt.scip import quicksum
        linfun = quicksum([c*x for c, x in zip(coeff, self.model.getVars())]) + d
        self.model.setObjective(linfun, sense=self.model.getObjectiveSense())


    def set_verbosity(self, level):
        """
        Set the verbosity level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 1.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.set_verbosity(1)                                            # optional - pyscipopt

        TODOs::

            - Currently, the output is written to stdout, even when running
            Jupyter: https://github.com/SCIP-Interfaces/PySCIPOpt/issues/116 .
            This should be fixed upstream
            - get access to more verbosity levels (e.g. via parameter settings)
        """
        if level == 0:
            self.model.hideOutput()
        elif level == 1:
            self.model.hideOutput(False)
        else:
            raise AssertionError('level must be "0" or "1"')


    def remove_constraint(self, i):
        r"""
        Remove a constraint from self.

        INPUT:

        - ``i`` -- index of the constraint to remove

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver='SCIP')                  # optional - pyscipopt
            sage: x, y = p['x'], p['y']                                         # optional - pyscipopt
            sage: p.add_constraint(2*x + 3*y <= 6)                              # optional - pyscipopt
            sage: p.add_constraint(3*x + 2*y <= 6)                              # optional - pyscipopt
            sage: p.add_constraint(x >= 0)                                      # optional - pyscipopt
            sage: p.set_objective(x + y + 7)                                    # optional - pyscipopt
            sage: p.set_integer(x); p.set_integer(y)                            # optional - pyscipopt
            sage: p.solve()                                                     # optional - pyscipopt
            9.0
            sage: p.remove_constraint(0)                                        # optional - pyscipopt
            sage: p.solve()                                                     # optional - pyscipopt
            10.0

        Removing fancy constraints does not make Sage crash::

            sage: MixedIntegerLinearProgram(solver = "SCIP").remove_constraint(-2) # optional - pyscipopt
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        if i < 0 or i >= self.nrows():
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        self.model.delCons(self.model.getConss()[i])

    def add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (real
          value).

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``name`` - an optional name for this row (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variables(5)                                            # optional - pyscipopt
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2.0, 2.0)   # optional - pyscipopt
            sage: p.row_bounds(0)                                               # optional - pyscipopt
            (2.0, 2.0)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo') # optional - pyscipopt
            sage: p.row_name(1)                                                 # optional - pyscipopt
            u'foo'
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        mvars = self.model.getVars()
        from pyscipopt.scip import quicksum
        linfun = quicksum([v*mvars[c] for c, v in coefficients])
        # we introduced patch 0001 for pyscipopt, in order to handle the case
        # when linfun is an empty expression.
        if name is None:
            name = ''

        if lower_bound is None:
            lower_bound = -self.model.infinity()
        if upper_bound is None:
            upper_bound = self.model.infinity()

        cons = lower_bound <= (linfun <= upper_bound)
        self.model.addCons(cons, name=name)

    def row(self, index):
        r"""
        Return a row

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_linear_constraint`` method.

        EXAMPLE::


            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variables(5)                                            # optional - pyscipopt
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)        # optional - pyscipopt
            sage: p.row(0)                                                      # optional - pyscipopt
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                                               # optional - pyscipopt
            (2.0, 2.0)
        """
        namedvars = [_.name for _ in self.model.getVars()]
        valslinear = self.model.getValsLinear(self.model.getConss()[index])
        pair = zip(*[[namedvars.index(_), v]
                       for _, v in valslinear.iteritems()])
        return (list(pair[0]), list(pair[1])) if pair != [] else ([], [])

    def row_bounds(self, index):
        """
        Return the bounds of a specific constraint.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variables(5)                                            # optional - pyscipopt
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)        # optional - pyscipopt
            sage: p.row_bounds(0)                                               # optional - pyscipopt
            (2.0, 2.0)
        """
        cons = self.model.getConss()[index]
        lhs = self.model.getLhs(cons)
        rhs = self.model.getRhs(cons)
        if lhs == -self.model.infinity():
            lhs = None
        if rhs == self.model.infinity():
            rhs = None
        return (lhs, rhs)

    def col_bounds(self, index):
        """
        Return the bounds of a specific variable.

        INPUT:

        - ``index`` (integer) -- the variable's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.col_bounds(0)                                               # optional - pyscipopt
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                                  # optional - pyscipopt
            sage: p.col_bounds(0)                                               # optional - pyscipopt
            (0.0, 5.0)
        """
        var = self.model.getVars()[index]
        lb = var.getLbOriginal()
        if lb == -self.model.infinity():
            lb = None
        ub = var.getUbOriginal()
        if ub == self.model.infinity():
            ub = None
        return (lb, ub)


    def add_col(self, indices, coeffs):
        """
        Add a column.

        INPUT:

        - ``indices`` (list of integers) -- this list constains the
          indices of the constraints in which the variable's
          coefficient is nonzero

        - ``coeffs`` (list of real values) -- associates a coefficient
          to the variable in each of the constraints in which it
          appears. Namely, the ith entry of ``coeffs`` corresponds to
          the coefficient of the variable in the constraint
          represented by the ith entry in ``indices``.

        .. NOTE::

            ``indices`` and ``coeffs`` are expected to be of the same
            length.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.ncols()                                                     # optional - pyscipopt
            0
            sage: p.nrows()                                                     # optional - pyscipopt
            0
            sage: p.add_linear_constraints(5, 0, None)                          # optional - pyscipopt
            sage: p.add_col(range(5), range(5))                                 # optional - pyscipopt
            sage: p.nrows()                                                     # optional - pyscipopt
            5
        """
        mcons = self.model.getConss()
        index = self.add_variable(lower_bound=-self.model.infinity())
        var = self.model.getVars()[index]

        for i, coeff in zip(indices, coeffs):
            self.model.addConsCoeff(var=var, cons=mcons[i], coeff=coeff)

    def solve(self):
        """
        Solve the problem.

        Sage uses SCIP's implementation of the branch-and-cut algorithm
        to solve the mixed-integer linear program.
        (If all variables are continuous, the algorithm reduces to solving the
        linear program by the simplex method.)

        EXAMPLE::

            sage: lp = MixedIntegerLinearProgram(solver = 'SCIP', maximization = False)
            sage: x, y = lp[0], lp[1]                                           # optional - pyscipopt
            sage: lp.add_constraint(-2*x + y <= 1)                              # optional - pyscipopt
            sage: lp.add_constraint(x - y <= 1)                                 # optional - pyscipopt
            sage: lp.add_constraint(x + y >= 2)                                 # optional - pyscipopt
            sage: lp.set_objective(x + y)                                       # optional - pyscipopt
            sage: lp.set_integer(x)                                             # optional - pyscipopt
            sage: lp.set_integer(y)                                             # optional - pyscipopt
            sage: lp.solve()                                                    # optional - pyscipopt
            2.0
            sage: lp.get_values([x, y])                                         # optional - pyscipopt
            [1.0, 1.0]

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_linear_constraints(5, 0, None)                          # optional - pyscipopt
            sage: p.add_col(range(5), range(5))                                 # optional - pyscipopt
            sage: p.solve()                                                     # optional - pyscipopt
            0
            sage: p.objective_coefficient(0,1)                                  # optional - pyscipopt
            sage: p.solve()                                                     # optional - pyscipopt
            Traceback (most recent call last):
            ...
            MIPSolverException: ...


            sage: lp = MixedIntegerLinearProgram(solver = "SCIP")               # optional - pyscipopt
            sage: v = lp.new_variable(nonnegative=True)                         # optional - pyscipopt
            sage: lp.add_constraint(v[1] +v[2] -2.0 *v[3], max=-1.0)            # optional - pyscipopt
            sage: lp.add_constraint(v[0] -4.0/3 *v[1] +1.0/3 *v[2], max=-1.0/3) # optional - pyscipopt
            sage: lp.add_constraint(v[0] +0.5 *v[1] -0.5 *v[2] +0.25 *v[3], max=-0.25) # optional - pyscipopt
            sage: lp.solve()                                                    # optional - pyscipopt
            0.0

        Solving a LP within the acceptable gap. No exception is raised, even if
        the result is not optimal. To do this, we try to compute the maximum
        number of disjoint balls (of diameter 1) in a hypercube::

            sage: g = graphs.CubeGraph(9)
            sage: p = MixedIntegerLinearProgram(solver = "SCIP")                # optional - pyscipopt

            sage: b = p.new_variable(binary=True)                               # optional - pyscipopt
            sage: p.set_objective(p.sum(b[v] for v in g))                       # optional - pyscipopt
            sage: for v in g:                                                   # optional - pyscipopt
            ....:     p.add_constraint(b[v]+p.sum(b[u] for u in g.neighbors(v)) <= 1) # optional - pyscipopt
            sage: p.add_constraint(b[v] == 1) # Force an easy non-0 solution    # optional - pyscipopt
            sage: p.solver_parameter("limits/absgap", 100)                      # optional - pyscipopt
            sage: p.solve() # rel tol 100                                       # optional - pyscipopt
            1

        """
        if (self.model.getStatus() != 'unknown') or (self.model.getStage() > 1):
            self.model.freeTransform()
        self.model.optimize()

        status = self.model.getStatus()

        if status == 'unbounded':
            raise MIPSolverException("SCIP : Solution is unbounded")
        elif status == 'infeasible':
            raise MIPSolverException("SCIP : There is no feasible solution")

        return 0

    def get_objective_value(self):
        """
        Returns the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variables(2)                                            # optional - pyscipopt
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)            # optional - pyscipopt
            sage: p.set_objective([2, 5])                                       # optional - pyscipopt
            sage: p.solve()                                                     # optional - pyscipopt
            0
            sage: p.get_objective_value()                                       # optional - pyscipopt
            7.5
            sage: p.get_variable_value(0) # abs tol 1e-15                       # optional - pyscipopt
            0.0
            sage: p.get_variable_value(1)                                       # optional - pyscipopt
            1.5
        """
        return self.model.getObjVal()

    def best_known_objective_bound(self):
        r"""
        Return the value of the currently best known bound.

        This method returns the current best upper (resp. lower) bound on the
        optimal value of the objective function in a maximization
        (resp. minimization) problem. It is equal to the output of
        :meth:`get_objective_value` if the MILP found an optimal solution, but
        it can differ if it was interrupted manually or after a time limit (cf
        :meth:`solver_parameter`).

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

        """
        return self.model.getPrimalbound()

    def get_relative_objective_gap(self):
        r"""
        Return the relative objective gap of the best known solution.

        For a minimization problem, this value is computed by
        `(\texttt{bestinteger} - \texttt{bestobjective}) / (1e-10 +
        |\texttt{bestobjective}|)`, where ``bestinteger`` is the value returned
        by :meth:`get_objective_value` and ``bestobjective`` is the value
        returned by :meth:`best_known_objective_bound`. For a maximization
        problem, the value is computed by `(\texttt{bestobjective} -
        \texttt{bestinteger}) / (1e-10 + |\texttt{bestobjective}|)`.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::


        TESTS:

        Just make sure that the variable *has* been defined, and is not just
        undefined::

        """
        return self.model.getGap()

    def get_variable_value(self, variable):
        """
        Returns the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variables(2)                                            # optional - pyscipopt
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)            # optional - pyscipopt
            sage: p.set_objective([2, 5])                                       # optional - pyscipopt
            sage: p.solve()                                                     # optional - pyscipopt
            0
            sage: p.get_objective_value()                                       # optional - pyscipopt
            7.5
            sage: p.get_variable_value(0) # abs tol 1e-15                       # optional - pyscipopt
            0.0
            sage: p.get_variable_value(1)                                       # optional - pyscipopt
            1.5
        """
        return self.model.getVal(self.model.getVars()[variable])

    def get_row_prim(self, i):
        r"""
        Returns the value of the auxiliary variable associated with i-th row.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "SCIP")                              # optional - pyscipopt
            sage: lp.add_variables(3)                                           # optional - pyscipopt
            2
            sage: lp.add_linear_constraint(zip([0, 1, 2], [8, 6, 1]), None, 48) # optional - pyscipopt
            sage: lp.add_linear_constraint(zip([0, 1, 2], [4, 2, 1.5]), None, 20) # optional - pyscipopt
            sage: lp.add_linear_constraint(zip([0, 1, 2], [2, 1.5, 0.5]), None, 8) # optional - pyscipopt
            sage: lp.set_objective([60, 30, 20])                                # optional - pyscipopt
            sage: lp.solver_parameter('presolving/maxrounds',0)                 # optional - pyscipopt
            sage: lp.solve()                                                    # optional - pyscipopt
            0
            sage: lp.get_objective_value()                                      # optional - pyscipopt
            280.0
            sage: lp.get_row_prim(0)                                            # optional - pyscipopt
            24.0
            sage: lp.get_row_prim(1)                                            # optional - pyscipopt
            20.0
            sage: lp.get_row_prim(2)                                            # optional - pyscipopt
            8.0
        """
        return self.model.getActivity(self.model.getConss()[i])


    def ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.ncols()                                                     # optional - pyscipopt
            0
            sage: p.add_variables(2)                                            # optional - pyscipopt
            1
            sage: p.ncols()                                                     # optional - pyscipopt
            2
        """
        return len(self.model.getVars())

    def nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.nrows()                                                     # optional - pyscipopt
            0
            sage: p.add_linear_constraints(2, 2, None)                          # optional - pyscipopt
            sage: p.nrows()                                                     # optional - pyscipopt
            2
        """
        return len(self.model.getConss())

    def col_name(self, index):
        """
        Return the ``index``th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variable(name='I am a variable')                        # optional - pyscipopt
            0
            sage: p.col_name(0)                                                 # optional - pyscipopt
            u'I am a variable'
        """
        return self.model.getVars()[index].name

    def row_name(self, index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1']) # optional - pyscipopt
            sage: p.row_name(0)                                                 # optional - pyscipopt
            u'Empty constraint 1'
        """
        return self.model.getConss()[index].name

    def is_variable_binary(self, index):
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.ncols()                                                     # optional - pyscipopt
            0
            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.set_variable_type(0,0)                                      # optional - pyscipopt
            sage: p.is_variable_binary(0)                                       # optional - pyscipopt
            True

        """
        return self.model.getVars()[index].vtype() == 'BINARY'

    def is_variable_integer(self, index):
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.ncols()                                                     # optional - pyscipopt
            0
            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.set_variable_type(0,1)                                      # optional - pyscipopt
            sage: p.is_variable_integer(0)                                      # optional - pyscipopt
            True
        """
        return self.model.getVars()[index].vtype() == 'INTEGER'


    def is_variable_continuous(self, index):
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.ncols()                                                     # optional - pyscipopt
            0
            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.is_variable_continuous(0)                                   # optional - pyscipopt
            True
            sage: p.set_variable_type(0,1)                                      # optional - pyscipopt
            sage: p.is_variable_continuous(0)                                   # optional - pyscipopt
            False

        """
        return self.model.getVars()[index].vtype() == 'CONTINUOUS'


    def is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.is_maximization()                                           # optional - pyscipopt
            True
            sage: p.set_sense(-1)                                               # optional - pyscipopt
            sage: p.is_maximization()                                           # optional - pyscipopt
            False
        """
        return self.model.getObjectiveSense() != 'minimize'

    def variable_upper_bound(self, index, value = False):
        """
        Return or define the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.col_bounds(0)                                               # optional - pyscipopt
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                                  # optional - pyscipopt
            sage: p.col_bounds(0)                                               # optional - pyscipopt
            (0.0, 5.0)

        TESTS:

        :trac:`14581`::

            sage: P = MixedIntegerLinearProgram(solver="SCIP")                  # optional - pyscipopt
            sage: x = P["x"]                                                    # optional - pyscipopt
            sage: P.set_max(x, 0)                                               # optional - pyscipopt
            sage: P.get_max(x)                                                  # optional - pyscipopt
            0.0

        Check that :trac:`10232` is fixed::

            sage: p = get_solver(solver="SCIP")                                 # optional - pyscipopt
            sage: p.variable_upper_bound(2)                                     # optional - pyscipopt
            Traceback (most recent call last):
            ...
            ValueError: The variable's id must satisfy 0 <= id < number_of_variables
            sage: p.variable_upper_bound(3, 5)                                  # optional - pyscipopt
            Traceback (most recent call last):
            ...
            ValueError: The variable's id must satisfy 0 <= id < number_of_variables

            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.variable_upper_bound(0, 'hey!')                             # optional - pyscipopt
            Traceback (most recent call last):
            ...
            TypeError: a float is required
        """
        if index < 0 or index >= self.ncols():
            raise ValueError("The variable's id must satisfy 0 <= id < number_of_variables")
        var = self.model.getVars()[index]
        if value is False:
            return var.getUbOriginal()
        else:
            self.model.chgVarUb(var=var, ub=value)

    def variable_lower_bound(self, index, value = False):
        """
        Return or define the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.col_bounds(0)                                               # optional - pyscipopt
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)                                  # optional - pyscipopt
            sage: p.col_bounds(0)                                               # optional - pyscipopt
            (5.0, None)

        TESTS:

        :trac:`14581`::

            sage: P = MixedIntegerLinearProgram(solver="SCIP")                  # optional - pyscipopt
            sage: x = P["x"]                                                    # optional - pyscipopt
            sage: P.set_min(x, 5)                                               # optional - pyscipopt
            sage: P.set_min(x, 0)                                               # optional - pyscipopt
            sage: P.get_min(x)                                                  # optional - pyscipopt
            0.0

        Check that :trac:`10232` is fixed::

            sage: p = get_solver(solver="SCIP")                                 # optional - pyscipopt
            sage: p.variable_lower_bound(2)                                     # optional - pyscipopt
            Traceback (most recent call last):
            ...
            ValueError: The variable's id must satisfy 0 <= id < number_of_variables
            sage: p.variable_lower_bound(3, 5)                                  # optional - pyscipopt
            Traceback (most recent call last):
            ...
            ValueError: The variable's id must satisfy 0 <= id < number_of_variables

            sage: p.add_variable()                                              # optional - pyscipopt
            0
            sage: p.variable_lower_bound(0, 'hey!')                             # optional - pyscipopt
            Traceback (most recent call last):
            ...
            TypeError: a float is required
        """
        if index < 0 or index >= self.ncols():
            raise ValueError("The variable's id must satisfy 0 <= id < number_of_variables")
        var = self.model.getVars()[index]
        if value is False:
            return var.getLbOriginal()
        else:
            self.model.chgVarLb(var=var, lb=value)

    def write_cip(self, filename):
        """
        Write the problem to a .cip file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variables(2)                                            # optional - pyscipopt
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)            # optional - pyscipopt
            sage: p.set_objective([2, 5])                                       # optional - pyscipopt
            sage: p.write_cip(os.path.join(SAGE_TMP, "lp_problem.cip"))         # optional - pyscipopt
            wrote problem to file ...
        """
        self.model.writeProblem(filename)

    def write_lp(self, filename):
        """
        Write the problem to a .lp file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variables(2)                                            # optional - pyscipopt
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)            # optional - pyscipopt
            sage: p.set_objective([2, 5])                                       # optional - pyscipopt
            sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))           # optional - pyscipopt
            wrote problem to file ...
        """
        filenamestr = filename
        fname, fext = splitext(filenamestr)

        if fext.lower() != '.lp':
            filenamestr = filename + '.lp'

        self.model.writeProblem(filenamestr)


    def write_mps(self, filename, modern):
        """
        Write the problem to a .mps file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")                               # optional - pyscipopt
            sage: p.add_variables(2)                                            # optional - pyscipopt
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)            # optional - pyscipopt
            sage: p.set_objective([2, 5])                                       # optional - pyscipopt
            sage: p.write_mps(os.path.join(SAGE_TMP, "lp_problem.mps"), 2)      # optional - pyscipopt
            wrote problem to file ...
        """
        filenamestr = filename
        fname, fext = splitext(filenamestr)

        if fext.lower() != '.mps':
            filenamestr = filename + '.mps'

        self.model.writeProblem(filenamestr)

    def __copy__(self):
        """
        Returns a copy of self.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "SCIP")                # optional - pyscipopt
            sage: b = p.new_variable()                                          # optional - pyscipopt
            sage: p.add_constraint(b[1] + b[2] <= 6)                            # optional - pyscipopt
            sage: p.set_objective(b[1] + b[2])                                  # optional - pyscipopt
            sage: copy(p).solve()                                               # optional - pyscipopt
            6.0
        """
        cp = type(self)(maximization=self.is_maximization())
        cp.problem_name(self.problem_name())
        for i, v in enumerate(self.model.getVars()):
            vtype = v.vtype
            cp.add_variable(self.variable_lower_bound(i),
                            self.variable_upper_bound(i),
                            binary=vtype == 'BINARY',
                            continuous=vtype == 'CONTINUOUS',
                            integer=vtype == 'INTEGER',
                            obj=self.objective_coefficient(i),
                            name=self.col_name(i))
        assert(self.ncols() == cp.ncols())

        for i in range(self.nrows()):
            cp.add_linear_constraint(zip(*self.row(i)),
                                     self.row_bounds(i)[0],
                                     self.row_bounds(i)[1],
                                     name=self.row_name(i))
        assert(self.nrows() == cp.nrows())
        return cp

    def solver_parameter(self, name, value=None):
        """
        Return or define a solver parameter

        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

        """
        if value is not None:
            if name.lower() == 'timelimit':
                self.model.setRealParam("limits/time", float(value))
            else:
                self.model.setParam(name, value)
        else:
            return self.model.getParam(name)



