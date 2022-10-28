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
    cdef int add_variable(self, obj=*, name=*) except -1
    cdef int add_variables(self, int, names=*) except -1
    cdef set_sense(self, int sense)
    cdef objective_coefficient(self, int variable, coeff=*)
    cdef set_objective(self, list coeff, d=*)
    cdef add_linear_constraint(self, constraints, name=*)
    cdef add_linear_constraints(self, int number, names=*)
    cdef int solve(self) except -1
    cdef get_objective_value(self)
    cdef get_variable_value(self, int variable)
    cdef dual_variable(self, int variable, sparse=*)
    cdef slack(self, int variable, sparse=*)
    cdef bint is_maximization(self)
    cdef row(self, int i)
    cdef int ncols(self)
    cdef int nrows(self)
    cdef problem_name(self, name=*)
    cdef row_name(self, int index)
    cdef col_name(self, int index)
    cdef solver_parameter(self, name, value=*)
    cdef zero(self)
    cdef base_ring(self)

    cdef obj_constant_term
    cdef dict matrices_dim

cdef GenericSDPBackend get_solver(solver=?, base_ring=?)
