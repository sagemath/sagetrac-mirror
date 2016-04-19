import sage.numerical.backends.glpk_backend as backend
from sage.numerical.backends.glpk_backend \
    import glp_bs, glp_nl, glp_nu
from sage.modules.all import vector
from sage.numerical.backends.abstract_backend_dictionary \
    import LPAbstractBackendDictionary


class LPGLPKBackendDictionary(LPAbstractBackendDictionary):
    r"""
    Construct a dictionary for an LP problem from an backend.

    INPUT:

        - ``backend`` -- the backend where the dictionary is
            constructed from

    OUTPUT:

       - a :class:`backend dictionary for an LP problem
       <LPGLPKBackendDictionary>`

    EXAMPLES:

    One needs an instance of :class:`GLPKBackend` to initialize
    this class::

        sage: from sage.numerical.backends.glpk_backend_dictionary \
              import LPGLPKBackendDictionary
        sage: p = MixedIntegerLinearProgram(maximization=True, solver="GLPK")
        sage: x = p.new_variable(nonnegative=True)
        sage: p.add_constraint(-x[0] + x[1] <= 2)
        sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
        sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
        sage: b = p.get_backend()
        sage: d = LPGLPKBackendDictionary(b)
        sage: d
        LP problem dictionary (use typeset mode to see details)
    """
    def __init__(self, backend):
        r"""
        See :class:`LPGLPKBackendDictionary` for documentation.

        TESTS::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: d = LPGLPKBackendDictionary(b)
            sage: TestSuite(d).run(skip=['_test_pickling'])

        An exception will be raised if the problem is not in standard form
        i.e. with <= constraints and >= 0 variable bounds::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(8 * x[0] + 2 * x[1], min=17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: d = LPGLPKBackendDictionary(b)
            Traceback (most recent call last):
            ...
            AttributeError: Problem constraints not in standard form.
        """
        super(LPGLPKBackendDictionary, self).__init__(backend)

    def basic_variables(self):
        r"""
        Return the basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES:

        Setting up the problem::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(\
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0

        Use function in :class:`LPGLPKBackendDictionary`:

            sage: d = LPGLPKBackendDictionary(b)

        Use function in :class:`InteractiveLPProblem`:

            sage: lp, basis = p.interactive_linear_program()
            sage: lpd = lp.dictionary(*basis)

        Compare results:

            sage: d.basic_variables()
            (x_0, x_1)
            sage: lpd.basic_variables()
            (x_0, x_1)
        """
        col_basics = tuple(
            self._x[i]
            for i in range(self._backend.ncols())
            if self._backend.get_col_stat(i) == glp_bs
        )
        row_basics = tuple(
            self._x[i + self._backend.ncols()]
            for i in range(self._backend.nrows())
            if self._backend.get_row_stat(i) == glp_bs
        )
        return vector(col_basics + row_basics)

    def constant_terms(self):
        r"""
        Return the constant terms of relations of ``self``.

        OUTPUT:

        - a vector.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(\
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: d = LPGLPKBackendDictionary(b)
            sage: d.constant_terms()
            (1.3, 3.3)
        """
        col_const = tuple(
            self._backend.get_variable_value(i)
            for i in range(self._backend.ncols())
            if self._backend.get_col_stat(i) == glp_bs
        )
        row_const = tuple(
            self._backend.row_bounds(i)[1] - self._backend.get_row_prim(i)
            for i in range(self._backend.nrows())
            if self._backend.get_row_stat(i) == glp_bs
        )
        return vector(col_const + row_const)

    def entering_coefficients(self):
        r"""
        Return coefficients of the entering variable.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(x[0] + x[1] - 7*x[2] + x[3] <= 22)
            sage: p.add_constraint(x[1] + 2*x[2] - x[3] <= 13)
            sage: p.add_constraint(5*x[0] + x[2] <= 11)
            sage: p.set_objective(2*x[0] + 3*x[1] + 4*x[2] + 13*x[3])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(\
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: d = LPGLPKBackendDictionary(b)
            sage: vars = d.nonbasic_variables()
            sage: vars
            (x_0, x_1, w_0, w_2)
            sage: d.enter(vars[0])
            sage: d.entering_coefficients()
            (5.0, 36.0, 26.0)
            sage: d.enter(vars[1])
            sage: d.entering_coefficients()
            (0.0, 1.0, 2.0)
        """
        if self._entering is None:
            raise ValueError("entering variable must be chosen to compute "
                             "its coefficients")

        index = tuple(self._x).index(self._entering)

        # Reverse signs for auxiliary variables
        if index < self._backend.ncols():
            tab_col = map(lambda (i, v):
                          (i, v) if i < self._backend.nrows() else (i, -v),
                          zip(*self._backend.eval_tab_col(
                              index + self._backend.nrows())))
        else:
            tab_col = map(lambda (i, v):
                          (i, v) if i < self._backend.nrows() else (i, -v),
                          zip(*self._backend.eval_tab_col(
                              index - self._backend.ncols())))

        # Sort the coefficients so coefficients of
        # problem variables comes first
        l = [0] * (self._backend.nrows())
        for (i, v) in tab_col:
            if i < self._backend.nrows():
                symbol = self._x[i + self._backend.ncols()]
            else:
                symbol = self._x[i - self._backend.nrows()]
            pos = tuple(self.basic_variables()).index(symbol)
            l[pos] = v

        return vector(l)

    def leaving_coefficients(self):
        r"""
        Return coefficients of the leaving variable.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(x[0] + x[1] - 7*x[2] + x[3] <= 22)
            sage: p.add_constraint(x[1] + 2*x[2] - x[3] <= 13)
            sage: p.add_constraint(5*x[0] + x[2] <= 11)
            sage: p.set_objective(2*x[0] + 3*x[1] + 4*x[2] + 13*x[3])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(\
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: d = LPGLPKBackendDictionary(b)
            sage: vars = d.basic_variables()
            sage: vars
            (x_2, x_3, w_1)
            sage: d.leave(vars[0])
            sage: d.leaving_coefficients()
            (5.0, 0.0, 0.0, 1.0)
            sage: d.leave(vars[1])
            sage: d.leaving_coefficients()
            (36.0, 1.0, 1.0, 7.0)
        """
        if self._leaving is None:
            raise ValueError("leaving variable must be chosen to compute "
                             "its coefficients")

        index = tuple(self._x).index(self._leaving)

        # Reverse signs for auxiliary variables
        if index < self._backend.ncols():
            raw_row = self._backend.eval_tab_row(
                index + self._backend.nrows())
            tab_row = map(lambda (i, v):
                          (i, v) if i < self._backend.nrows() else (i, -v),
                          zip(*raw_row))
        else:
            raw_row = self._backend.eval_tab_row(
                index - self._backend.ncols())
            tab_row = map(lambda (i, v):
                          (i, -v) if i < self._backend.nrows() else (i, v),
                          zip(*raw_row))

        l = [0] * (self._backend.ncols())
        for (i, v) in tab_row:
            if i < self._backend.nrows():
                symbol = self._x[i + self._backend.ncols()]
            else:
                symbol = self._x[i - self._backend.nrows()]
            pos = tuple(self.nonbasic_variables()).index(symbol)
            l[pos] = v

        return vector(l)

    def nonbasic_variables(self):
        r"""
        Return non-basic variables of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES:

        Setting up the problem::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(\
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0

        Use function in :class:`LPGLPKBackendDictionary`:

            sage: d = LPGLPKBackendDictionary(b)

        Use function in :class:`InteractiveLPProblem`:

            sage: lp, basis = p.interactive_linear_program()
            sage: lpd = lp.dictionary(*basis)

        Compare results:

            sage: d.nonbasic_variables()
            (w_0, w_1)
            sage: lpd.nonbasic_variables()
            (w_0, w_1)
        """
        col_nonbasics = tuple(
            self._x[i]
            for i in range(self._backend.ncols())
            if self._backend.get_col_stat(i) != glp_bs
        )
        row_nonbasics = tuple(
            self._x[i + self._backend.ncols()]
            for i in range(self._backend.nrows())
            if self._backend.get_row_stat(i) != glp_bs
        )
        return vector(col_nonbasics + row_nonbasics)

    def objective_coefficients(self):
        r"""
        Return coefficients of the objective of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES:

        Setting up the problem::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(\
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0

        Use function in :class:`LPGLPKBackendDictionary`:

            sage: d = LPGLPKBackendDictionary(b)

        Use function in :class:`InteractiveLPProblem`:

            sage: lp, basis = p.interactive_linear_program()
            sage: lpd = lp.dictionary(*basis)

        Compare results:

            sage: d.objective_coefficients()
            (-0.58, -0.76)
            sage: lpd.objective_coefficients() # rel tol 1e-9
            (-0.5800000000000001, -0.76)
        """
        col_coefs = tuple(
            self._backend.get_col_dual(i)
            for i in range(self._backend.ncols())
            if self._backend.get_col_stat(i) != glp_bs
        )
        row_coefs = tuple(
            -self._backend.get_row_dual(i)
            for i in range(self._backend.nrows())
            if self._backend.get_row_stat(i) != glp_bs
        )
        return vector(col_coefs + row_coefs)

    def update(self):
        r"""
        Update ``self`` using previously set entering and leaving variables.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(x[0] + x[1] - 7*x[2] + x[3] <= 22)
            sage: p.add_constraint(x[1] + 2*x[2] - x[3] <= 13)
            sage: p.add_constraint(5*x[0] + x[2] <= 11)
            sage: p.set_objective(2*x[0] + 3*x[1] + 4*x[2] + 13*x[3])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(\
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: d = LPGLPKBackendDictionary(b)
            sage: d.objective_value()
            1331.0
            sage: d.nonbasic_variables()
            (x_0, x_1, w_0, w_2)
            sage: d.enter(d.nonbasic_variables()[0])
            sage: d.basic_variables()
            (x_2, x_3, w_1)
            sage: d.leave(d.basic_variables()[0])
            sage: d.objective_value()
            1331.0
            sage: d.update()
            sage: d.basic_variables()
            (x_0, x_3, w_1)
            sage: d.nonbasic_variables()
            (x_1, x_2, w_0, w_2)
            sage: d.objective_value()
            261.8

        TESTS:

        An error will be raised if the pivot selected is zero::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(x[0] + x[1] - 7*x[2] + x[3] <= 22)
            sage: p.add_constraint(x[1] + 2*x[2] - x[3] <= 13)
            sage: p.add_constraint(5*x[0] + x[2] <= 11)
            sage: p.set_objective(2*x[0] + 3*x[1] + 4*x[2] + 13*x[3])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(\
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: d = LPGLPKBackendDictionary(b)
            sage: d.leave(d.basic_variables()[0])
            sage: d.leaving_coefficients()
            (5.0, 0.0, 0.0, 1.0)
            sage: d.enter(d.nonbasic_variables()[1])
            sage: d.leave(d.basic_variables()[0])
            sage: d.update()
            Traceback (most recent call last):
            ...
            ValueError: incompatible choice of entering and leaving variables
        """
        entering = self._entering
        if entering is None:
            raise ValueError("entering variable must be set before updating")

        leaving = self._leaving
        if leaving is None:
            raise ValueError("leaving variable must be set before updating")

        matching_index = tuple(self.basic_variables()).index(leaving)
        coef = self.entering_coefficients()[matching_index]
        if coef == 0:
            raise ValueError("incompatible choice of entering and leaving "
                             "variables")

        entering_index = tuple(self._x).index(entering)
        if entering_index < self._backend.ncols():
            self._backend.set_col_stat(entering_index, glp_bs)
        else:
            self._backend.set_row_stat(entering_index - self._backend.ncols(),
                                       glp_bs)

        leaving_index = tuple(self._x).index(leaving)
        if leaving_index < self._backend.ncols():
            self._backend.set_col_stat(leaving_index, glp_nl)
        else:
            self._backend.set_row_stat(leaving_index - self._backend.ncols(),
                                       glp_nu)

        if self._backend.warm_up() != 0:
            raise AttributeError("Warm up failed.")

    def add_row(self, nonbasic_coef, constant, slack_variable,
                integer_slack_variable=False):
        r"""
        Update a dictionary with an additional row based on a given dictionary.

        INPUT:

        - ``nonbasic_coef``-- a list of nonbasic coefficients for the new row

        - ``constant``-- a number of the constant term for the new row

        - ``slack_variable``-- a string of the name for the new slack variable

        - ``integer_slack_variable``-- (default: False) a boolean value
        indicating if the new slack variable is integer or not.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_backend_dictionary \
                  import LPGLPKBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True, \
                                                solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(x[0] + x[1] - 7*x[2] + x[3] <= 22)
            sage: p.add_constraint(x[1] + 2*x[2] - x[3] <= 13)
            sage: p.add_constraint(5*x[0] + x[2] <= 11)
            sage: p.set_objective(2*x[0] + 3*x[1] + 4*x[2] + 13*x[3])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(\
                backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: d = LPGLPKBackendDictionary(b)
            sage: d.basic_variables()
            (x_2, x_3, w_1)
            sage: d.nonbasic_variables()
            (x_0, x_1, w_0, w_2)
            sage: d.objective_coefficients()
            (-486.0, -10.0, -13.0, -95.0)
            sage: d.add_row(range(3,7), 2, 'z_0')
            sage: d.objective_coefficients()
            (-486.0, -10.0, -13.0, -95.0)
            sage: d.basic_variables()
            (x_2, x_3, w_1, z_0)
            sage: d.leave(d.basic_variables()[3])
            sage: d.leaving_coefficients()
            (3.0, 4.0, 5.0, 6.0)
            sage: b.solve()
            0
            sage: d.basic_variables()
            (x_2, x_3, w_1, z_0)
            sage: d.nonbasic_variables()
            (x_0, x_1, w_0, w_2)

        Variables have 0 as their coefficient will not show up in the
        tableau:

            sage: d.add_row(range(-2, 2), 5, 'z_1')
            sage: d.get_backend().row(4)
            ([2, 1, 0], [-1.0, -1.0, -7.0])
        """
        if len(nonbasic_coef) != self._backend.ncols():
            raise ValueError("Length of nonbasic coefficients incompatible")

        # Convert to problem variable coefficients
        coef_pairs, constant = (
            self._nonbasic_coef_to_problem_coef_(nonbasic_coef, constant)
        )
        self._backend.add_linear_constraint(
            coef_pairs, None, constant, slack_variable
        )

        # Update buffered variables
        row_index = self._backend.nrows() - 1
        self._load_new_variable_(
            index=row_index,
            name=self._backend.row_name(row_index),
            auxiliary=True
        )

        # Update basis status in the backend
        self._backend.set_row_stat(self._backend.nrows() - 1, glp_bs)
        if self._backend.warm_up() != 0:
            raise AttributeError("Warm up failed.")
