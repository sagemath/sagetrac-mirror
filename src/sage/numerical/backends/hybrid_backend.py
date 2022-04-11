r"""
Hybrid Backend

AUTHORS:

- Nathann Cohen (2010-10)      : generic_backend template
- Matthias Koeppe (2016-03)    : this backend
- Yuan Zhou (2020-02)          : this backend

"""

#*****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#       Copyright (C) 2016 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.numerical.mip import MIPSolverException
from sage.modules.all import vector
from sage.numerical.backends.generic_backend import GenericBackend

class HybridBackend(GenericBackend):

    """
    A backend that delegates to a sequence of other backends.

    All setters delegate to all backends.
    All getters delegate to the last backend (the exact backend).

    """

    def __init__(self, maximization = True, backends = []):
        """
        Constructor

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))

        """

        if backends is []:
            raise ValueError("HybridBackend needs at least one backend to work with")
        self.backends = backends
        self.set_verbosity(0)

        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

    def base_ring(self):
        """
        Return the base ring.

        OUTPUT:

        A ring. The coefficients that the chosen solver supports.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.base_ring()
            Rational Field
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"), base_ring=AA)
            sage: p.base_ring()
            Algebraic Real Field
        """
        return self.backends[-1].base_ring()

    def add_variable(self, lower_bound=0, upper_bound=None,
                           binary=False, continuous=True, integer=False,
                           obj=None, name=None):
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

        INPUT:

        - ``lower_bound`` - the lower bound of the variable (default: 0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default: ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of this variable in the objective function (default: 0)

        - ``name`` - an optional name for the newly added variable (default: ``None``).

        OUTPUT: The index of the newly created variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable(continuous=True, integer=True)
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.add_variable(name='x',obj=1)
            1
            sage: p.col_name(1)
            'x'
            sage: p.objective_coefficient(1)
            1
        """
        for b in self.backends:
            result = b.add_variable(lower_bound, upper_bound, binary, continuous, integer, obj, name)
        return result

    def add_variables(self, n, lower_bound=0, upper_bound=None, binary=False, continuous=True, integer=False, obj=None, names=None):
        """
        Add ``n`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        INPUT:

        - ``n`` - the number of new variables (must be > 0)

        - ``lower_bound`` - the lower bound of the variable (default: 0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default: ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of all variables in the objective function (default: 0)

        - ``names`` - optional list of names (default: ``None``)

        OUTPUT: The index of the variable created last.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.ncols()
            0
            sage: p.add_variables(5)
            4
            sage: p.ncols()
            5
            sage: p.add_variables(2, names=['a','b'])
            6
        """
        for b in self.backends:
            result = b.add_variables(n, lower_bound, upper_bound, binary, continuous, integer, obj, names)
        return result

    def set_variable_type(self, variable, vtype):
        """
        Set the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer)

            *  1  Integer
            *  0  Binary
            *  -1  Continuous

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,-1)
            sage: p.is_variable_continuous(0)
            True
        """
        for b in self.backends:
            b.set_variable_type(variable, vtype)

    def set_sense(self, sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer)

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        for b in self.backends:
            b.set_sense(sense)

    def objective_coefficient(self, variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective
        function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2
        """
        if coeff is None:
            return self.backends[-1].objective_coefficient(variable)
        else:
            for b in self.backends:
                b.objective_coefficient(variable, coeff)

    def objective_constant_term(self, d=None):
        """
        Set or get the constant term in the objective function

        INPUT:

        - ``d`` (double) -- its coefficient.  If `None` (default), return the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.objective_constant_term()
            0
            sage: p.objective_constant_term(42)
            sage: p.objective_constant_term()
            42
        """
        if d is None:
            return self.backends[-1].objective_constant_term()
        else:
            for b in self.backends:
                b.objective_constant_term(d)

    def set_objective(self, coeff, d=0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose i-th element is the
          coefficient of the i-th variable in the objective function.

        - ``d`` (real) -- the constant term in the linear function (set to `0` by default)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1, 1, 2, 1, 3]

        Constants in the objective function are respected::

            sage: p = MixedIntegerLinearProgram(solver = ("GLPK", "InteractiveLP"))
            sage: x,y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.solve()
            47/5
        """
        for b in self.backends:
            b.set_objective(coeff, d)

    def set_verbosity(self, level):
        """
        Set the log (verbosity) level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.set_verbosity(2)
        """
        for b in self.backends:
            b.set_verbosity(level)

    def remove_constraint(self, i):
        r"""
        Remove a constraint.

        INPUT:

        - ``i`` -- index of the constraint to remove.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver=("GLPK", "InteractiveLP"))
            sage: v = p.new_variable(nonnegative=True)
            sage: x,y = v[0], v[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.solve()
            47/5
            sage: p.remove_constraint(0)
            sage: p.solve()
            10
            sage: p.get_values([x,y])
            [0, 3]
        """
        for b in self.backends:
            b.remove_constraint(i)

    ## def remove_constraints(self, constraints): # inherite from generic_backend

    def add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` -- an iterable of pairs ``(i, v)``. In each
          pair, ``i`` is a variable index (integer) and ``v`` is a
          value (element of :meth:`base_ring`).

        - ``lower_bound`` -- element of :meth:`base_ring` or
          ``None``. The lower bound.

        - ``upper_bound`` -- element of :meth:`base_ring` or
          ``None``. The upper bound.

        - ``name`` -- string or ``None``. Optional name for this row.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)
            (2, 2)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1, 1, name='foo')
            sage: p.row_name(1)
            'foo'
        """
        coefficients = list(coefficients)
        for b in self.backends:
            b.add_linear_constraint(coefficients, lower_bound, upper_bound, name)

    ## def add_linear_constraint_vector(self, degree, coefficients, lower_bound, upper_bound, name=None):  # inherite from generic_backend

    def add_col(self, indices, coeffs):
        """
        Add a column.

        INPUT:

        - ``indices`` (list of integers) -- this list contains the
          indices of the constraints in which the variable's
          coefficient is nonzero

        - ``coeffs`` (list of real values) -- associates a coefficient
          to the variable in each of the constraints in which it
          appears. Namely, the i-th entry of ``coeffs`` corresponds to
          the coefficient of the variable in the constraint
          represented by the i-th entry in ``indices``.

        .. NOTE::

            ``indices`` and ``coeffs`` are expected to be of the same
            length.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5
        """
        for b in self.backends:
            b.add_col(indices, coeffs)

    def add_linear_constraints(self, number, lower_bound, upper_bound, names=None):
        """
        Add constraints.

        INPUT:

        - ``number`` (integer) -- the number of constraints to add.

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``names`` - an optional list of names (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.row(4)
            ([], [])
            sage: p.row_bounds(4)
            (0, None)
        """
        for b in self.backends:
            b.add_linear_constraints(number, lower_bound, upper_bound, names)

    def solve(self):
        """
        Solve the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.solve()
            0
            sage: p.objective_coefficient(0,1)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...

        Test on LPs with non-negative variables.
        The first LP has a unique optimal solution; the second one is unbounded; the third one is infeasible::

            sage: p = MixedIntegerLinearProgram(maximization=True, solver=("GLPK","InteractiveLP"))
            sage: v = p.new_variable(nonnegative=True)
            sage: x, y = v[0], v[1]
            sage: p.add_constraint(-2*x + y <= 1)
            sage: p.add_constraint(x - y <= 1)
            sage: p.add_constraint(x + y >= 2)
            sage: p.set_objective(-y)
            sage: p.solve()
            -1/2
            sage: p.set_objective(y)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: InteractiveLP: Problem is unbounded
            sage: p.add_constraint(x + y <= 0)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...

        Test on LP that has free varialbe or equality constraint::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: h = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: h.add_variables(1, lower_bound=None, upper_bound=None)
            0
            sage: h.add_variables(1, lower_bound=0, upper_bound=None)
            1
            sage: h.add_linear_constraint([(0,2),(1,-1)],-1,None)
            sage: h.add_linear_constraint([(0,1),(1,-1)],None, 1)
            sage: h.add_linear_constraint([(0,1),(1,1)],2,2)
            sage: h.set_objective([0,-1])
            sage: h.solve()
            0
            sage: h.get_objective_value()
            -1/2
        """
        from sage.numerical.backends.interactivelp_backend import InteractiveLPBackend
        if self.backends[-1].parent() is InteractiveLPBackend:
            for b in self.backends[:-1]:
                try:
                    b.solve()
                except MIPSolverException:
                    # infeasilbe or unbounded LP.
                    return self.backends[-1].solve()
                try:
                    basic_variables = []
                    k = 1
                    for i in range(b.ncols()):
                        is_free_variable = bool(b.col_bounds(i) == (None, None))
                        if b.is_variable_basic(i):
                            if is_free_variable and (b.get_variable_value(i) < 0):
                                basic_variables.append(k+1)
                            else:
                                basic_variables.append(k)
                        if is_free_variable:
                            k += 2
                        else:
                            k += 1
                    for i in range(b.nrows()):
                        is_equality_constraint = not ((b.row_bounds(i)[0] == None) or (b.row_bounds(i)[1] == None))
                        if b.is_slack_variable_basic(i):
                            basic_variables.append(k)
                        if is_equality_constraint:
                            k += 1
                            basic_variables.append(k)
                        k += 1
                except NotImplementedError:
                    # no basis information is available in this backend. continue to the next backend.
                    continue
                else:
                    self.backends[-1].set_dictionary(basic_variables=basic_variables)
                    return self.backends[-1].solve()
        return self.backends[-1].solve()

    def get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behavior is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            15/2
            sage: p.get_variable_value(0)
            0
            sage: p.get_variable_value(1)
            3/2
        """
        return self.backends[-1].get_objective_value()

    def get_variable_value(self, variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behavior is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0,1), (1, 2)], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            15/2
            sage: p.get_variable_value(0)
            0
            sage: p.get_variable_value(1)
            3/2
        """
        return self.backends[-1].get_variable_value(variable)

    def ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.ncols()
            0
            sage: p.add_variables(2)
            1
            sage: p.ncols()
            2
        """
        return self.backends[-1].ncols()

    def nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(2, 0, None)
            sage: p.nrows()
            2
        """
        return self.backends[-1].nrows()

    def is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        return self.backends[-1].is_maximization()

    def problem_name(self, name = None):
        """
        Return or define the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.problem_name("There_once_was_a_french_fry")
            sage: print(p.problem_name())
            There_once_was_a_french_fry
        """
        if name == None:
            return self.backends[-1].problem_name()
        else:
            for b in self.backends:
                b.problem_name(name)

    ## def best_known_objective_bound(self)
    ## def get_relative_objective_gap(self)
    ## def write_lp(self, name):
    ##     """
    ##     Write the problem to a ``.lp`` file

    ##     INPUT:

    ##     - ``filename`` (string)

    ##     EXAMPLES::

    ##         sage: from sage.numerical.backends.generic_backend import get_solver
    ##         sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
    ##         sage: p.add_variables(2)
    ##         2
    ##         sage: p.add_linear_constraint([(0, 1], (1, 2)], None, 3)
    ##         sage: p.set_objective([2, 5])
    ##         sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))
    ##     """
    ##     raise NotImplementedError()

    ## def write_mps(self, name, modern):
    ##     """
    ##     Write the problem to a ``.mps`` file

    ##     INPUT:

    ##     - ``filename`` (string)

    ##     EXAMPLES::

    ##         sage: from sage.numerical.backends.generic_backend import get_solver
    ##         sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
    ##         sage: p.add_variables(2)
    ##         2
    ##         sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3)
    ##         sage: p.set_objective([2, 5])
    ##         sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))
    ##     """
    ##     raise NotImplementedError()

    def row(self, i):
        """
        Return a row

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_linear_constraint`` method.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 0, None)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
        """
        return self.backends[-1].row(i)

    def row_bounds(self, index):
        """
        Return the bounds of a specific constraint.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row_bounds(0)
            (2, 2)
        """
        return self.backends[-1].row_bounds(index)

    def col_bounds(self, index):
        """
        Return the bounds of a specific variable.

        INPUT:

        - ``index`` (integer) -- the variable's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variable(lower_bound=None)
            0
            sage: p.col_bounds(0)
            (None, None)
            sage: p.variable_lower_bound(0, 0)
            sage: p.col_bounds(0)
            (0, None)
        """
        return self.backends[-1].col_bounds(index)

    def is_variable_binary(self, index):
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_binary(0)
            False

        """
        return self.backends[-1].is_variable_binary(index)

    def is_variable_integer(self, index):
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_integer(0)
            False
        """
        return self.backends[-1].is_variable_integer(index)

    def is_variable_continuous(self, index):
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True

        """
        return self.backends[-1].is_variable_continuous(index)

    def row_name(self, index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1'])
            sage: p.row_name(0)
            'Empty constraint 1'

        """
        return self.backends[-1].row_name(index)

    def col_name(self, index):
        """
        Return the ``index``-th column name

        INPUT:

        - ``index`` (integer) -- the column id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variable(name="I_am_a_variable")
            0
            sage: p.col_name(0)
            'I_am_a_variable'
        """
        return self.backends[-1].col_name(index)

    def variable_upper_bound(self, index, value = False):
        """
        Return or define the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variable(lower_bound=None)
            0
            sage: p.col_bounds(0)
            (None, None)
            sage: p.variable_upper_bound(0) is None
            True
            sage: p.variable_upper_bound(0, 0)
            sage: p.col_bounds(0)
            (None, 0)
            sage: p.variable_upper_bound(0)
            0
            sage: p.variable_upper_bound(0, None)
            sage: p.variable_upper_bound(0) is None
            True
        """
        if value is False:
            return self.backends[-1].variable_upper_bound(index, value)
        else:
            for b in self.backends:
                b.variable_upper_bound(index, value)

    def variable_lower_bound(self, index, value = False):
        """
        Return or define the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has no lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
            sage: p.add_variable(lower_bound=None)
            0
            sage: p.col_bounds(0)
            (None, None)
            sage: p.variable_lower_bound(0) is None
            True
            sage: p.variable_lower_bound(0, 0)
            sage: p.col_bounds(0)
            (0, None)
            sage: p.variable_lower_bound(0)
            0
            sage: p.variable_lower_bound(0, None)
            sage: p.variable_lower_bound(0) is None
            True
        """
        if value is False:
            return self.backends[-1].variable_lower_bound(index, value)
        else:
            for b in self.backends:
                b.variable_lower_bound(index, value)

    ## def solver_parameter(self, name, value = None):
    ##     """
    ##     Return or define a solver parameter

    ##     INPUT:

    ##     - ``name`` (string) -- the parameter

    ##     - ``value`` -- the parameter's value if it is to be defined,
    ##       or ``None`` (default) to obtain its current value.

    ##     .. NOTE::

    ##        The list of available parameters is available at
    ##        :meth:`~sage.numerical.mip.MixedIntegerLinearProgram.solver_parameter`.

    ##     EXAMPLES::

    ##         sage: from sage.numerical.backends.generic_backend import get_solver
    ##         sage: p = get_solver(solver = ("GLPK", "InteractiveLP"))
    ##         sage: p.solver_parameter("timelimit")
    ##         sage: p.solver_parameter("timelimit", 60)
    ##         sage: p.solver_parameter("timelimit")
    ##     """
    ##     raise NotImplementedError()

    def is_variable_basic(self, index):
        """
        Test whether the given variable is basic.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver=("GLPK", "InteractiveLP"))
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_variable_basic(0)
            True
            sage: b.is_variable_basic(1)
            False
        """
        return self.backends[-1].is_variable_basic(index)

    def is_variable_nonbasic_at_lower_bound(self, index):
        """
        Test whether the given variable is nonbasic at lower bound.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver=("GLPK", "InteractiveLP"))
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_variable_nonbasic_at_lower_bound(0)
            False
            sage: b.is_variable_nonbasic_at_lower_bound(1)
            True
        """
        return self.backends[-1].is_variable_nonbasic_at_lower_bound(index)

    def is_slack_variable_basic(self, index):
        """
        Test whether the slack variable of the given row is basic.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver=("GLPK", "InteractiveLP"))
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_slack_variable_basic(0)
            True
            sage: b.is_slack_variable_basic(1)
            False
        """
        return self.backends[-1].is_slack_variable_basic(index)

    def is_slack_variable_nonbasic_at_lower_bound(self, index):
        """
        Test whether the given variable is nonbasic at lower bound.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver=("GLPK", "InteractiveLP"))
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_slack_variable_nonbasic_at_lower_bound(0)
            False
            sage: b.is_slack_variable_nonbasic_at_lower_bound(1)
            True
        """
        return self.backends[-1].is_slack_variable_nonbasic_at_lower_bound(index)
