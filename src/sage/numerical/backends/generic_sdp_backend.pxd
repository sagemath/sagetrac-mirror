##############################################################################
#       Copyright (C) 2014 Ingolfur Edvardsson <ingolfured@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

cdef class GenericSDPBackend:
    cpdef int add_variable(self, obj=*, name=*) except -1
    cpdef int add_variables(self, int, names=*) except -1
    cpdef set_sense(self, int sense)
    cpdef objective_coefficient(self, int variable, coeff=*)
    cpdef set_objective(self, list coeff, d=*)
    cpdef add_linear_constraint(self, constraints, name=*)
    cpdef add_linear_constraints(self, int number, names=*)
    cpdef int solve(self) except -1
    cpdef get_objective_value(self)
    cpdef get_variable_value(self, int variable)
    cpdef bint is_maximization(self)
    cpdef row(self, int i)
    cpdef int ncols(self)
    cpdef int nrows(self)
    cpdef problem_name(self, char * name = *)
    cpdef row_name(self, int index)
    cpdef col_name(self, int index)
    cpdef solver_parameter(self, name, value=*)
    cpdef zero(self)
    cpdef base_ring(self)

    cpdef obj_constant_term
    cdef dict matrices_dim

cpdef GenericSDPBackend get_solver(solver = ?)
