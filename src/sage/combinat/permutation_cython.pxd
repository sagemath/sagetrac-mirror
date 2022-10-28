from cpython.array cimport array

cdef void reset_swap(int n, int *c, int *o)
cdef int next_swap(int n, int *c, int *o)
cdef bint next_perm(array l)
cdef map_to_list(array l, tuple values, int n)
cdef list left_action_same_n(list l, list r)
cdef list right_action_same_n(list l, list r)
cdef list left_action_product(list l, list r)
cdef list right_action_product(list l, list r)

