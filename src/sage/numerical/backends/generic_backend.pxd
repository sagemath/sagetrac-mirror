##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.structure.sage_object cimport SageObject

# We inherit from SageObject to make some testing infrastructure available.

cdef class GenericBackend (SageObject):
    cdef int add_variable(self, lower_bound=*, upper_bound=*, binary=*, continuous=*, integer=*, obj=*, name=*) except -1
    cdef int add_variables(self, int, lower_bound=*, upper_bound=*, binary=*, continuous=*, integer=*, obj=*, names=*) except -1
    cdef set_variable_type(self, int variable, int vtype)
    cdef set_sense(self, int sense)
    cdef objective_coefficient(self, int variable, coeff=*)
    cdef objective_constant_term(self, d=*)
    cdef set_objective(self, list coeff, d=*)
    cdef set_verbosity(self, int level)
    cdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=*)
    cdef add_linear_constraint_vector(self, degree, coefficients, lower_bound, upper_bound, name=*)
    cdef remove_constraint(self, int)
    cdef remove_constraints(self, constraints)
    cdef add_col(self, indices, coeffs)
    cdef add_linear_constraints(self, int number, lower_bound, upper_bound, names=*)
    cdef int solve(self) except -1
    cdef get_objective_value(self)
    cdef best_known_objective_bound(self)
    cdef get_relative_objective_gap(self)
    cdef get_variable_value(self, int variable)
    cdef bint is_maximization(self)
    cdef write_lp(self, name)
    cdef write_mps(self, name, int modern)
    cdef row(self, int i)
    cdef int ncols(self)
    cdef int nrows(self)
    cdef bint is_variable_binary(self, int)
    cdef bint is_variable_integer(self, int)
    cdef bint is_variable_continuous(self, int)
    cdef problem_name(self, name = *)
    cdef row_bounds(self, int index)
    cdef col_bounds(self, int index)
    cdef row_name(self, int index)
    cdef col_name(self, int index)
    cdef variable_upper_bound(self, int index, value = *)
    cdef variable_lower_bound(self, int index, value = *)
    cdef solver_parameter(self, name, value=*)
    cdef zero(self)
    cdef base_ring(self)
    cdef __copy__(self)
    cdef copy(self)
    cdef bint is_variable_basic(self, int index)
    cdef bint is_variable_nonbasic_at_lower_bound(self, int index)
    cdef bint is_slack_variable_basic(self, int index)
    cdef bint is_slack_variable_nonbasic_at_lower_bound(self, int index)

    cdef object obj_constant_term

cdef GenericBackend get_solver(constraint_generation = ?, solver = ?, base_ring = ?)
