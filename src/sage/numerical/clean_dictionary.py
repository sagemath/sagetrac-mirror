r"""
Interactive Simplex Method floating-point helper class

This module provides a "clean" dictionary view on a dictionary with
floating-point numbers.  Cleaning means to change almost-zeros to exact
zeros, allowing the Interactive Simplex Method to recognize primal and
dual feasibility and to avoid pivoting on zero pivot elements.  By
extension, it also changes almost-integers to exact integers, for the
benefit of mixed-integer programming.
"""

from sage.numerical.interactive_simplex_method \
    import LPDictionary, LPAbstractDictionary
from sage.modules.all import vector
from sage.rings.all import ZZ
from sage.matrix.all import matrix

class LPCleanDictionary(LPAbstractDictionary):
    r"""
    Construct a clean dictionary for an LP problem from a dictionary.

    INPUT:

        - ``dictionary`` -- the dictionary to be cleaned

        - ``epsilon`` -- (default: 1e-9) the tolerance of the cleaning process

    OUTPUT:

       - a :class:`dictionary for an LP problem
       <LPCleanDictionary>`

    EXAMPLES:

    One needs an instance inherited from :class:`LPAbstractDictionary`
    to initialize this class:

    Setting up a backend dictionary as the input::

        sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
        sage: x = p.new_variable(nonnegative=True)
        sage: p.add_constraint(-x[0] + x[1] <= 2)
        sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
        sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
        sage: lp, basis = p.interactive_lp_problem()
        sage: d = lp.dictionary(*basis)
        sage: view(d)  # not tested

    Create the clean dictionary using the backend dictionary::

        sage: from sage.numerical.clean_dictionary \
              import LPCleanDictionary
        sage: clean = LPCleanDictionary(d, epsilon=1e-5)
        sage: clean
        LP problem dictionary (use typeset mode to see details)
        sage: view(clean.dictionary())  # not tested
    """
    def __init__(self, dictionary, epsilon=1e-9):
        r"""
        See :class:`LPGLPKBackendDictionary` for documentation.

        TESTS::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
            sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: view(d)  # not tested
            sage: clean = LPCleanDictionary(d, epsilon=1e-5)
            sage: clean
            LP problem dictionary (use typeset mode to see details)
            sage: view(clean.dictionary())  # not tested
        """
        super(LPCleanDictionary, self).__init__()
        self._dict = dictionary
        self._epsilon = epsilon

    def _round_(self, var):
        r"""
        Round a number based on the epsilon.

        INPUT:

            - ``var`` -- the variable to be rounded

        EXAMPLE::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: clean._round_(0.1)
            0.100000000000000
            sage: clean._round_(0.00001)
            0
        """
        nearest_int = ZZ(round(var))
        if abs(var - nearest_int) <= self._epsilon:
            var = nearest_int
        return var

    def leave(self, v):
        r"""
        Set ``v`` as the leaving variable for the associated dictionary.

        INPUT:

        - ``v`` -- a basic variable of the associated dictionary,
        can be given as a string, an actual variable, or an integer
        interpreted as the index of a variable

        EXAMPLE::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
            sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: clean.basic_variables()
            (w_0, w_1)
            sage: clean.leave(clean.basic_variables()[0])
            sage: d.leaving()
            w_0
        """
        self._dict.leave(v)

    def enter(self, v):
        r"""
        Set ``v`` as the entering variable for the associated dictionary.

        INPUT:

        - ``v`` -- a basic variable of the associated dictionary,
        can be given as a string, an actual variable, or an integer
        interpreted as the index of a variable

        EXAMPLE::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
            sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: clean.nonbasic_variables()
            (m_0, m_1)
            sage: clean.enter(clean.nonbasic_variables()[0])
            sage: d._entering
            m_0
        """
        self._dict.enter(v)

    def _round_all_(self, iterable):
        r"""
        Round all element in an iterable.

        INPUT:

            - ``iterable`` -- the list to tranverse through

        EXAMPLE::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: l = [0.01, 0.000001]
            sage: clean._round_all_(l)
            (0.0100000000000000, 0.000000000000000)
        """
        return vector([self._round_(e) for e in iterable])

    def leaving_coefficients(self):
        r"""
        Return coefficients of the leaving variable.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
            sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: clean.leave(clean.basic_variables()[0])
            sage: clean.leaving_coefficients()
            (-1, 1)
        """
        coefs = self._dict.leaving_coefficients()
        return self._round_all_(coefs)

    def constant_terms(self):
        r"""
        Return the constant terms of relations of ``self``.

        OUTPUT:

        - a vector.

        EXAMPLES::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
            sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: d.constant_terms()
            (2.0, 3e-07)
            sage: clean.constant_terms()
            (2, 0)
        """
        terms = self._dict.constant_terms()
        return self._round_all_(terms)

    def objective_coefficients(self):
        r"""
        Return coefficients of the objective of ``self``.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
            sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: clean.objective_coefficients()
            (-5.5, -2.1)
        """
        coefs = self._dict.objective_coefficients()
        return self._round_all_(coefs)

    def objective_name(self):
        return self._dict.objective_name()

    def objective_value(self):
        r"""
        Return the value of the objective value.

        OUTPUT:

        - a number

        EXAMPLES::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
            sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: clean.objective_value()
            0
        """
        val = self._dict.objective_value()
        return self._round_(val)

    def basic_variables(self):
        r"""
        Return the basic variables of the associated dictionary.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
            sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: clean.basic_variables()
            (w_0, w_1)
        """
        return self._dict.basic_variables()

    def nonbasic_variables(self):
        r"""
        Return the non-basic variables of the associated dictionary.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
            sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: clean.nonbasic_variables()
            (m_0, m_1)
        """
        return self._dict.nonbasic_variables()

    def dictionary(self):
        r"""
        Return a cleaned regular LP dictionary of the associated dictionary.

        OUTPUT:

        - an :class:`LP dictionary <LPDictionary>`

        EXAMPLES::

            sage: from sage.numerical.clean_dictionary \
                  import LPCleanDictionary
            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(0.001 * x[0] + 0.0005 * x[1] <= 0.0000003)
            sage: p.set_objective(-5.5 * x[0] + -2.1 * x[1])
            sage: lp, basis = p.interactive_lp_problem()
            sage: d = lp.dictionary(*basis)
            sage: clean = LPCleanDictionary(d, epsilon=1e-4)
            sage: clean_dict = clean.dictionary()
            sage: view(clean_dict)  #not tested
        """
        old_leaving = self._dict._leaving
        rows = []
        for i, var in enumerate(self.basic_variables()):
            self.leave(self.basic_variables()[i])
            rows.append(self.leaving_coefficients())

        A = matrix(rows)
        D = LPDictionary(A,
                         self.constant_terms(),
                         self.objective_coefficients(),
                         self.objective_value(),
                         self.basic_variables(),
                         self.nonbasic_variables(),
                         self.objective_name())
        self.leave(old_leaving)
        D.enter(self._dict._entering)
        D.leave(old_leaving)
        return D