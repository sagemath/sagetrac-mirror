from cpython.array cimport array

cdef void reset_swap(int n, int *c, int *o)
cdef int next_swap(int n, int *c, int *o)
cpdef bint next_perm(array l)
cpdef map_to_list(array l, tuple values, int n)
cpdef list left_action_same_n(list l, list r)
cpdef list right_action_same_n(list l, list r)
cpdef list left_action_product(list l, list r)
cpdef list right_action_product(list l, list r)

