from sage.numerical.interactive_simplex_method import *


class LPAbstractBackendDictionary(LPAbstractDictionary):
    r"""
    Construct an abstract dictionary for an LP problem from a backend.

    INPUT:

        - ``backend`` -- the backend that the dictionary is
            constructed from

    OUTPUT:

       - a :class:`backend dictionary for an LP problem
       <LPAbstractBackendDictionary>`

    EXAMPLES:

    One needs an instance of :class:`MixedIntegerLinearProgram` to initialize
    this class::

        sage: from sage.numerical.backends.abstract_backend_dictionary \
              import LPAbstractBackendDictionary
        sage: p = MixedIntegerLinearProgram(maximization=True)
        sage: x = p.new_variable(nonnegative=True)
        sage: p.add_constraint(-x[0] + x[1] <= 2)
        sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
        sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
        sage: b = p.get_backend()
        sage: d = LPAbstractBackendDictionary(b)
        sage: d
        LP problem dictionary (use typeset mode to see details)
    """
    def __init__(self, backend):
        super(LPAbstractBackendDictionary, self).__init__()
        self._backend = backend
        for i in range(self._backend.nrows()):
            if self._backend.row_bounds(i)[0] != None \
               or self._backend.row_bounds(i)[1] == None:
                raise AttributeError("Problem constraints "
                                     "not in standard form.")

        for i in range(self._backend.ncols()):
            if self._backend.variable_lower_bound(i) == None:
                raise AttributeError("Problem variables "
                                     "not in standard form.")

        col_vars = tuple(
            LPAbstractBackendDictionary._format_(
                name=self._backend.col_name(i),
                symbol='x', index=i
            )
            for i in range(self._backend.ncols())
        )
        row_vars = tuple(
            LPAbstractBackendDictionary._format_(
                name=self._backend.row_name(i),
                symbol='w', index=i
            )
            for i in range(self._backend.nrows())
        )
        self._names = ", ".join(col_vars + row_vars)
        self._R = PolynomialRing(self._backend.base_ring(),
                                 self._names, order="neglex")
        self._x = self._R.gens()

    def __eq__(self, other):
        r"""
        Check if two LP problem dictionaries have the same
        reference.

        INPUT:

        - ``other`` -- anything

        OUTPUT:

        - ``True`` if ``other`` is an :class:`LPDictionary` with its
          reference the same as ``self``, ``False`` otherwise.

        TESTS:

        Setting up the problem::

            sage: from sage.numerical.backends.abstract_backend_dictionary \
                import LPAbstractBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: d = LPAbstractBackendDictionary(b)

        Test when two problems have the same reference::

            sage: d2 = d
            sage: d2 == d
            True

        Test when two problems have the same construction::

            sage: d3 = LPAbstractBackendDictionary(copy(p).get_backend())
            sage: d3 == d
            False
        """
        return (isinstance(other, LPAbstractBackendDictionary) and
                self._backend == other._backend)

    @staticmethod
    def _format_(name='', symbol='x', infix='_', index='0'):
        r"""
        Returns a proper name for a given parameter.

        INPUT::

        - ``name`` -- (defualt: '') the original name of the variable

        -``symbol`` -- (defualt: 'x') the symbol for the new name

        -``infix`` -- (default: '_') the character separate the symbol
        and its index for the new name

        -``index`` -- (default: '0') the index following the infix for
        the new name

        TESTS:

        If a name is given, then returns the name itself::

            sage: from sage.numerical.backends.abstract_backend_dictionary \
                  import LPAbstractBackendDictionary
            sage: LPAbstractBackendDictionary._format_('Name')
            'Name'

        However, if the name given is in the form 'symbol[index]', then it
        will be converted to 'symbol+infix+index' i.e. 'symbol_index'
        by default::

            sage: LPAbstractBackendDictionary._format_('x[3]')
            'x_3'
            sage: LPAbstractBackendDictionary._format_('x[2]', infix='~')
            'x~2'

        If no name is given, then a newly create name will be returned in the
        form of 'symbol+infix+index'::

            sage: LPAbstractBackendDictionary._format_(symbol='w', index='7')
            'w_7'
        """
        if name:
            return name.replace('[', infix).strip(']')
        else:
            return symbol + infix + str(index)

    def objective_value(self):
        r"""
        Return the value of the objective value.

        OUTPUT:

        - a number

        EXAMPLES::

            sage: from sage.numerical.backends.abstract_backend_dictionary \
                  import LPAbstractBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solve()
            0
            sage: d = LPAbstractBackendDictionary(b)
            sage: d.objective_value()
            14.08
        """
        return self._backend.get_objective_value()

    def get_backend(self):
        r"""
        Return the backend used to create the dictionary.

        OUTPUT:

        - The corresponding backend associated with self

        EXAMPLES::

            sage: from sage.numerical.backends.abstract_backend_dictionary \
                  import LPAbstractBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solve()
            0
            sage: d = LPAbstractBackendDictionary(b)
            sage: d.get_backend()
            <sage.numerical.backends.coin_backend.CoinBackend object at ...>
        """
        return self._backend

    def dictionary(self):
        r"""
        Return a regular LP dictionary matching ``self``.

        OUTPUT:

        - an :class:`LP dictionary <LPDictionary>`

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
            sage: view(d.dictionary()) # not tested

        Zeta is used as default problem name, and it can be changed::

            sage: b.problem_name("beta")
            sage: view(d.dictionary()) # not tested
        """
        rows = []
        for i in range(self._backend.nrows()):
            self.leave(self.basic_variables()[i])
            rows.append(self.leaving_coefficients())
        import sage.matrix.constructor as construc
        m = construc.matrix(rows)
        name = self._backend.problem_name()
        if not name:
            name = 'zeta'
        D = LPDictionary(m,
                         self.constant_terms(),
                         self.objective_coefficients(),
                         self.objective_value(),
                         self.basic_variables(),
                         self.nonbasic_variables(),
                         name)
        D._entering = self._entering
        D._leaving = self._leaving
        return D

    def _load_new_variable_(self, index, name, auxiliary=True):
        r"""
        Load a new variable to buffer.

        INPUT:

        - ``index`` -- the index of the new variable

        - ``name`` -- the name of the new variable

        - ``auxiliary`` -- (default: True) if true, the symbol of the
        new varible will be 'w', 'x' otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.abstract_backend_dictionary \
                  import LPAbstractBackendDictionary
            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Coin")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] + 2.1 * x[1])
            sage: b = p.get_backend()
            sage: b.solve()
            0
            sage: d = LPAbstractBackendDictionary(b)
            sage: d._x
            (x_0, x_1, w_0, w_1)
        """
        if auxiliary:
            symbol = 'w'
        else:
            symbol = 'x'

        self._names += ', '
        self._names += LPAbstractBackendDictionary._format_(
            name=name, symbol=symbol, index=index
        )
        self._R = PolynomialRing(self._backend.base_ring(),
                                 self._names, order="neglex")
        self._x = self._R.gens()

    def _nonbasic_coef_to_problem_coef_(self, nonbasic_coef, constant):
        r"""
        Returns coefficients of nonbasic variables rewritten in terms of
        problem variables.

        INPUT:

        - ``nonbasic_coef`` -- a list of nonbasic coefficients for the new row

        - ``constant`` -- a number of the constant term for the new row

        OUTPUT:

        - ``coef_pairs`` -- a list of tuples that indicate problem variable
        indices and coefficients

        - ``constant`` -- the re-calculated row bound

        EXAMPLE::

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
            sage: nonbasic_coef = [3, 4, 5, 6]
            sage: constant = 11
            sage: pairs, const = \
            d._nonbasic_coef_to_problem_coef_(nonbasic_coef, constant)
            sage: pairs
            [(0, -29.0), (1, -11.0)]
            sage: const
            -63.0
        """
        coefs = [0] * self._backend.ncols()
        for i, var in enumerate(self.nonbasic_variables()):
            index = self._x.index(var)
            if index < self._backend.ncols():
                coefs[index] += nonbasic_coef[i]
            else:
                row_pos = index - self._backend.ncols()
                row = self._backend.row(row_pos)
                for j, v in zip(*row):
                    coefs[j] -= nonbasic_coef[i] * v
                upper_bound = self._backend.row_bounds(row_pos)[1]
                constant -= nonbasic_coef[i] * upper_bound

        coef_pairs = [(i, coefs[i]) for i in range(self._backend.ncols())
                      if coefs[i] != 0]
        return coef_pairs, constant

