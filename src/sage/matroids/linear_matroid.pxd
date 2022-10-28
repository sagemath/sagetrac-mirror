from sage.data_structures.bitset cimport bitset_t

from .matroid cimport Matroid
from .basis_exchange_matroid cimport BasisExchangeMatroid
from .lean_matrix cimport LeanMatrix, GenericMatrix, BinaryMatrix, TernaryMatrix, QuaternaryMatrix


cdef class LinearMatroid(BasisExchangeMatroid):
    cdef LeanMatrix _A, _representation
    cdef long *_prow
    cdef object _zero, _one

    cdef _forget(self)
    cdef base_ring(self)
    cdef characteristic(self)

    cdef list _setup_internal_representation(self, matrix, reduced_matrix, ring, keep_initial_representation)
    cdef __exchange_value(self, long x, long y)

    cdef representation(self, B=*, reduced=*, labels=*, order=*, lift_map=*)
    cdef _current_rows_cols(self, B=*)
    cdef representation_vectors(self)
    cdef LeanMatrix _basic_representation(self, B=*)
    cdef LeanMatrix _reduced_representation(self, B=*)
    
    
    cdef bint _is_field_isomorphism(self, LinearMatroid other, morphism)
    cdef is_field_equivalent(self, other)
    cdef is_field_isomorphism(self, other, morphism)
    # cdef is_field_isomorphic(self, other)  # TODO: currently only works as ``def``
    cdef _fast_isom_test(self, other)

    cdef _minor(self, contractions, deletions)
    cdef dual(self)
    cdef has_line_minor(self, k, hyperlines=*, certificate=*)
    cdef has_field_minor(self, N)

    cdef _exchange_value(self, e, f)
    cdef fundamental_cycle(self, B, e)
    cdef fundamental_cocycle(self, B, e)

    cdef _line_ratios(self, F)
    cdef _line_length(self, F)

    cdef _line_cross_ratios(self, F)
    cdef cross_ratios(self, hyperlines=*)
    cdef cross_ratio(self, F, a, b, c, d)
    cdef _line_cross_ratio_test(self, F, x, fundamentals)
    cdef _cross_ratio_test(self, x, fundamentals, hyperlines=*)

    cdef linear_extension(self, element, chain=*, col=*)
    cdef linear_coextension(self, element, cochain=*, row=*)
    cdef _linear_extensions(self, element, chains)
    cdef _linear_coextensions(self, element, cochains)
    cdef _extend_chains(self, C, f, fundamentals=*)
    cdef _linear_extension_chains(self, F, fundamentals=*)
    cdef linear_extension_chains(self, F=*, simple=*, fundamentals=*)
    cdef linear_coextension_cochains(self, F=*, cosimple=*, fundamentals=*)
    cdef linear_extensions(self, element=*, F=*, simple=*, fundamentals=*)
    cdef linear_coextensions(self, element=*, F=*, cosimple=*, fundamentals=*)

    cdef _is_3connected_shifting(self, certificate=*)
    cdef _is_4connected_shifting(self, certificate=*)

    cdef is_valid(self)

cdef class BinaryMatroid(LinearMatroid):
    cdef tuple _b_invariant, _b_partition
    cdef BinaryMatrix _b_projection, _eq_part

    cdef base_ring(self)
    cdef characteristic(self)

    cdef _current_rows_cols(self, B=*)
    cdef LeanMatrix _basic_representation(self, B=*)
    cdef LeanMatrix _reduced_representation(self, B=*)

    cdef  __fundamental_cocircuit(self, bitset_t, long x)

    cdef _is_isomorphic(self, other, certificate=*)

    cdef _minor(self, contractions, deletions)

    cdef _make_invariant(self)
    cdef _invariant(self)
    cdef bicycle_dimension(self)
    cdef brown_invariant(self)
    cdef _principal_tripartition(self)
    cdef BinaryMatrix _projection(self)
    cdef BinaryMatrix _projection_partition(self)
    cdef _fast_isom_test(self, other)

    cdef is_graphic(self)
    cdef is_valid(self)


cdef class TernaryMatroid(LinearMatroid):
    cdef object _two
    cdef tuple _t_invariant, _t_partition
    cdef TernaryMatrix _t_projection

    cdef base_ring(self)
    cdef characteristic(self)

    cdef _current_rows_cols(self, B=*)
    cdef LeanMatrix _basic_representation(self, B=*)
    cdef LeanMatrix _reduced_representation(self, B=*)

    cdef  __fundamental_cocircuit(self, bitset_t, long x)

    cdef _is_isomorphic(self, other, certificate=*)

    cdef _minor(self, contractions, deletions)

    cdef _make_invariant(self)
    cdef _invariant(self)
    cdef bicycle_dimension(self)
    cdef character(self)
    cdef _principal_quadripartition(self)
    cdef TernaryMatrix _projection(self)
    cdef _fast_isom_test(self, other)

    cdef is_valid(self)

cdef class QuaternaryMatroid(LinearMatroid):
    cdef object _x_zero, _x_one
    cdef tuple _q_invariant, _q_partition
    cdef QuaternaryMatrix _q_projection

    cdef base_ring(self)
    cdef characteristic(self)

    cdef _current_rows_cols(self, B=*)
    cdef LeanMatrix _basic_representation(self, B=*)
    cdef LeanMatrix _reduced_representation(self, B=*)

    cdef  __fundamental_cocircuit(self, bitset_t, long x)

    cdef _is_isomorphic(self, other, certificate=*)

    cdef _minor(self, contractions, deletions)

    cdef _make_invariant(self)
    cdef _invariant(self)
    cdef bicycle_dimension(self)
    cdef _principal_tripartition(self)
    cdef _fast_isom_test(self, other)

    cdef is_valid(self)

cdef class RegularMatroid(LinearMatroid):
    cdef _bases_count, _r_invariant
    cdef _r_projection, _r_hypergraph
    cdef _hypergraph_vertex_partition, _hypergraph_tuples

    cdef base_ring(self)
    cdef characteristic(self)

    cdef _is_isomorphic(self, other, certificate=*)

    cdef _invariant(self)
    cdef _fast_isom_test(self, other)

    cdef bases_count(self)
    cdef _projection(self)
    cdef _hypergraph(self)
    cdef _hypertest(self, other)
    cdef has_line_minor(self, k, hyperlines=*, certificate=*)
    cdef _linear_extension_chains(self, F, fundamentals=*)

    cdef is_graphic(self)
    cdef is_valid(self)
