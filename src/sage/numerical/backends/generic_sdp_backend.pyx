r"""
Generic Backend for SDP solvers

This class only lists the methods that should be defined by any
interface with a SDP Solver. All these methods immediately raise
``NotImplementedError`` exceptions when called, and are obviously
meant to be replaced by the solver-specific method. This file can also
be used as a template to create a new interface : one would only need
to replace the occurences of ``"Nonexistent_SDP_solver"`` by the
solver's name, and replace ``GenericSDPBackend`` by
``SolverName(GenericSDPBackend)`` so that the new solver extends this
class.

AUTHORS:

- Ingolfur Edvardsson (2014-07): initial implementation

"""

#*****************************************************************************
#       Copyright (C) 2014 Ingolfur Edvardsson <ingolfured@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef class GenericSDPBackend:

    cpdef base_ring(self):
        """
        The base ring

        TESTS::

            sage: from sage.numerical.backends.generic_sdp_backend import GenericSDPBackend
            sage: GenericSDPBackend().base_ring()
            Real Double Field
        """
        from sage.rings.all import RDF
        return RDF

    cpdef zero(self):
        """
        Zero of the base ring

        TESTS::

            sage: from sage.numerical.backends.generic_sdp_backend import GenericSDPBackend
            sage: GenericSDPBackend().zero()
            0.0
        """
        return self.base_ring().zero()

    cpdef int add_variable(self, obj=0.0, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

        INPUT:

        - ``obj`` - (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` - an optional name for the newly added variable (default: ``None``).

        OUTPUT: The index of the newly created variable

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")    # optional - Nonexistent_LP_solver
            sage: p.ncols()                                           # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                    # optional - Nonexistent_LP_solver
            0
            sage: p.ncols()                                           # optional - Nonexistent_LP_solver
            1
            sage: p.add_variable(name='x',obj=1.0)                    # optional - Nonexistent_LP_solver
            3
            sage: p.col_name(3)                                       # optional - Nonexistent_LP_solver
            'x'
            sage: p.objective_coefficient(3)                          # optional - Nonexistent_LP_solver
            1.0
        """
        raise NotImplementedError()

    cpdef int add_variables(self, int n, names=None) except -1:
        """
        Add ``n`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        INPUT:

        - ``n`` - the number of new variables (must be > 0)

        - ``obj`` - (optional) coefficient of all variables in the objective function (default: 0.0)

        - ``names`` - optional list of names (default: ``None``)

        OUTPUT: The index of the variable created last.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")    # optional - Nonexistent_LP_solver
            sage: p.ncols()                                           # optional - Nonexistent_LP_solver
            0
            sage: p.add_variables(5)                                  # optional - Nonexistent_LP_solver
            4
            sage: p.ncols()                                           # optional - Nonexistent_LP_solver
            5
            sage: p.add_variables(2, lower_bound=-2.0, integer=True, names=['a','b']) # optional - Nonexistent_LP_solver
            6
        """
        raise NotImplementedError()

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.is_maximization()                              # optional - Nonexistent_LP_solver
            True
            sage: p.set_sense(-1)                              # optional - Nonexistent_LP_solver
            sage: p.is_maximization()                              # optional - Nonexistent_LP_solver
            False
        """
        raise NotImplementedError()

    cpdef  objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective
        function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            1
            sage: p.objective_coefficient(0)                         # optional - Nonexistent_LP_solver
            0.0
            sage: p.objective_coefficient(0,2)                       # optional - Nonexistent_LP_solver
            sage: p.objective_coefficient(0)                         # optional - Nonexistent_LP_solver
            2.0
        """
        raise NotImplementedError()

    cpdef  set_objective(self, list coeff, d=0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        - ``d`` (double) -- the constant term in the linear function (set to `0` by default)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")    # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                                 # optional - Nonexistent_LP_solver
            5
            sage: p.set_objective([1, 1, 2, 1, 3])                   # optional - Nonexistent_LP_solver
            sage: map(lambda x :p.objective_coefficient(x), range(5))  # optional - Nonexistent_LP_solver
            [1.0, 1.0, 2.0, 1.0, 3.0]

        Constants in the objective function are respected::

        """
        raise NotImplementedError()


    cpdef add_linear_constraint(self, coefficients, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (real
          value).

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``name`` - an optional name for this row (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                               # optional - Nonexistent_LP_solver
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2.0, 2.0) # optional - Nonexistent_LP_solver
            sage: p.row(0)                                         # optional - Nonexistent_LP_solver
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])                   # optional - Nonexistent_LP_solver
            sage: p.row_bounds(0)                                  # optional - Nonexistent_LP_solver
            (2.0, 2.0)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo') # optional - Nonexistent_LP_solver
            sage: p.row_name(-1)                                                          # optional - Nonexistent_LP_solver
            "foo"
        """
        raise NotImplementedError()


    cpdef add_linear_constraints(self, int number, names=None):
        """
        Add constraints.

        INPUT:

        - ``number`` (integer) -- the number of constraints to add.

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``names`` - an optional list of names (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")   # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                                # optional - Nonexistent_LP_solver
            5
            sage: p.add_linear_constraints(5, None, 2)          # optional - Nonexistent_LP_solver
            sage: p.row(4)                                      # optional - Nonexistent_LP_solver
            ([], [])
            sage: p.row_bounds(4)                               # optional - Nonexistent_LP_solver
            (None, 2.0)
        """
        raise NotImplementedError()

    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises ``SDPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.add_linear_constraints(5, 0, None)             # optional - Nonexistent_LP_solver
            sage: p.add_col(range(5), range(5))                    # optional - Nonexistent_LP_solver
            sage: p.solve()                                        # optional - Nonexistent_LP_solver
            0
            sage: p.objective_coefficient(0,1)                 # optional - Nonexistent_LP_solver
            sage: p.solve()                                        # optional - Nonexistent_LP_solver
            Traceback (most recent call last):
            ...
            SDPSolverException: ...
        """
        raise NotImplementedError()

    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.add_variables(2)                               # optional - Nonexistent_LP_solver
            2
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3) # optional - Nonexistent_LP_solver
            sage: p.set_objective([2, 5])                          # optional - Nonexistent_LP_solver
            sage: p.solve()                                        # optional - Nonexistent_LP_solver
            0
            sage: p.get_objective_value()                          # optional - Nonexistent_LP_solver
            7.5
            sage: p.get_variable_value(0)                          # optional - Nonexistent_LP_solver
            0.0
            sage: p.get_variable_value(1)                          # optional - Nonexistent_LP_solver
            1.5
        """

        raise NotImplementedError()

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.add_variables(2)                              # optional - Nonexistent_LP_solver
            2
            sage: p.add_linear_constraint([(0,1), (1, 2)], None, 3) # optional - Nonexistent_LP_solver
            sage: p.set_objective([2, 5])                         # optional - Nonexistent_LP_solver
            sage: p.solve()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.get_objective_value()                         # optional - Nonexistent_LP_solver
            7.5
            sage: p.get_variable_value(0)                         # optional - Nonexistent_LP_solver
            0.0
            sage: p.get_variable_value(1)                         # optional - Nonexistent_LP_solver
            1.5
        """

        raise NotImplementedError()

    cpdef int ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.add_variables(2)                               # optional - Nonexistent_LP_solver
            2
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            2
        """

        raise NotImplementedError()

    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.nrows()                                        # optional - Nonexistent_LP_solver
            0
            sage: p.add_linear_constraints(2, 2.0, None)         # optional - Nonexistent_LP_solver
            sage: p.nrows()                                      # optional - Nonexistent_LP_solver
            2
        """

        raise NotImplementedError()

    cpdef bint is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.is_maximization()                             # optional - Nonexistent_LP_solver
            True
            sage: p.set_sense(-1)                             # optional - Nonexistent_LP_solver
            sage: p.is_maximization()                             # optional - Nonexistent_LP_solver
            False
        """
        raise NotImplementedError()

    cpdef problem_name(self, char * name = NULL):
        """
        Return or define the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")   # optional - Nonexistent_LP_solver
            sage: p.problem_name("There once was a french fry") # optional - Nonexistent_LP_solver
            sage: print p.get_problem_name()                        # optional - Nonexistent_LP_solver
            There once was a french fry
        """

        raise NotImplementedError()

    cpdef row(self, int i):
        """
        Return a row

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_linear_constraint`` method.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                               # optional - Nonexistent_LP_solver
            5
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2) # optional - Nonexistent_LP_solver
            sage: p.row(0)                                     # optional - Nonexistent_LP_solver
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)                              # optional - Nonexistent_LP_solver
            (2.0, 2.0)
        """
        raise NotImplementedError()



    cpdef row_name(self, int index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_linear_constraints(1, 2, None, name="Empty constraint 1")  # optional - Nonexistent_LP_solver
            sage: p.row_name(0)                                     # optional - Nonexistent_LP_solver
            'Empty constraint 1'

        """
        raise NotImplementedError()

    cpdef col_name(self, int index):
        """
        Return the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable(name="I am a variable")            # optional - Nonexistent_LP_solver
            1
            sage: p.col_name(0)                                     # optional - Nonexistent_LP_solver
            'I am a variable'
        """
        raise NotImplementedError()


    cpdef solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

        .. NOTE::

           The list of available parameters is available at
           :meth:`~sage.numerical.sdp.SemidefiniteProgram.solver_parameter`.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.solver_parameter("timelimit")                   # optional - Nonexistent_LP_solver
            sage: p.solver_parameter("timelimit", 60)               # optional - Nonexistent_LP_solver
            sage: p.solver_parameter("timelimit")                   # optional - Nonexistent_LP_solver
        """
        raise NotImplementedError()



default_solver = None

def default_sdp_solver(solver = None):
    """
    Returns/Sets the default SDP Solver used by Sage

    INPUT:

    - ``solver`` -- defines the solver to use:


        - CVXOPT (``solver="CVXOPT"``). See the `CVXOPT
          <http://cvxopt.org/>`_ web site.

        ``solver`` should then be equal to one of ``"CVXOPT"``.

        - If ``solver=None`` (default), the current default solver's name is
          returned.

    OUTPUT:

    This function returns the current default solver's name if ``solver = None``
    (default). Otherwise, it sets the default solver to the one given. If this
    solver does not exist, or is not available, a ``ValueError`` exception is
    raised.

    EXAMPLE::

        sage: former_solver = default_sdp_solver()
        sage: default_sdp_solver("Cvxopt")
        sage: default_sdp_solver()
        'Cvxopt'
        sage: default_sdp_solver("Yeahhhhhhhhhhh")
        Traceback (most recent call last):
        ...
        ValueError: 'solver' should be set to 'CVXOPT' or None.
        sage: default_sdp_solver(former_solver)
    """
    global default_solver

    if solver is None:

        if default_solver is not None:
            return default_solver

        else:
            for s in ["Cvxopt"]:
                try:
                    default_sdp_solver(s)
                    return s
                except ValueError:
                    pass

    solver = solver.capitalize()


    if solver == "Cvxopt":
        try:
            from sage.numerical.backends.cvxopt_sdp_backend import CVXOPTSDPBackend
            default_solver = solver
        except ImportError:
            raise ValueError("CVXOPT is not available. Please refer to the documentation to install it.")

    else:
        raise ValueError("'solver' should be set to 'CVXOPT' or None.")

cpdef GenericSDPBackend get_solver(solver = None):
    """
    Return a solver according to the given preferences

    INPUT:

    - ``solver`` -- 1 solver should be available through this class:

        - CVXOPT (``solver="CVXOPT"``). See the `CVXOPT
          <http://cvxopt.org/>`_ web site.

        ``solver`` should then be equal to one of ``"CVXOPT"`` or ``None``.
          If ``solver=None`` (default), the default solver is used (see ``default_sdp_solver`` method.

    .. SEEALSO::

    - :func:`default_sdp_solver` -- Returns/Sets the default SDP solver.

    EXAMPLE::

        sage: from sage.numerical.backends.generic_sdp_backend import get_solver
        sage: p = get_solver()
    """
    if solver is None:
        solver = default_sdp_solver()

    else:
        solver = solver.capitalize()


    if solver == "Cvxopt":
        from sage.numerical.backends.cvxopt_sdp_backend import CVXOPTSDPBackend
        return CVXOPTSDPBackend()

    else:
        raise ValueError("'solver' should be set to 'CVXOPT' or None (in which case the default one is used).")
