from sage.structure.sage_object cimport SageObject

cdef class Matroid(SageObject):
    cdef public __custom_name
    cdef public _custom_name
    cdef public _cached_info
    cdef int _stored_full_rank
    cdef int _stored_size

    # virtual methods
    cdef groundset(self)
    cdef _rank(self, X)

    # internal methods, assuming verified input
    cdef _max_independent(self, X)
    cdef _circuit(self, X)
    cdef _fundamental_circuit(self, B, e)
    cdef _closure(self, X)
    cdef _corank(self, X)
    cdef _max_coindependent(self, X)
    cdef _cocircuit(self, X)
    cdef _fundamental_cocircuit(self, B, e)
    cdef _coclosure(self, X)
    cdef _augment(self, X, Y)

    cdef _is_independent(self, X)
    cdef _is_basis(self, X)
    cdef _is_circuit(self, X)
    cdef _is_closed(self, X)
    cdef _is_coindependent(self, X)
    cdef _is_cobasis(self, X)
    cdef _is_cocircuit(self, X)
    cdef _is_coclosed(self, X)

    cdef _minor(self, contractions, deletions)
    cdef _has_minor(self, N, bint certificate=*)
    cdef _line_length(self, F)
    cdef _extension(self, element, hyperplanes)

    cdef inline __subset(self, X):
        """
        Convert ``X`` to a ``frozenset`` and check that it is a subset
        of the groundset.

        See ``_subset`` for the corresponding Python method.
        """
        S = frozenset(X)
        if not self.groundset().issuperset(S):
            raise ValueError(f"{X!r} is not a subset of the groundset")
        return S

    cdef inline __subset_all(self, X):
        """
        If ``X`` is ``None``, return the groundset.

        Otherwise, do like ``_subset``:
        convert ``X`` to a ``frozenset`` and check that it is a subset
        of the groundset.

        See ``_subset_all`` for the corresponding Python method.
        """
        if X is None:
            return self.groundset()
        S = frozenset(X)
        if not self.groundset().issuperset(S):
            raise ValueError(f"{X!r} is not a subset of the groundset")
        return S

    # ** user-facing methods **
    cdef size(self)

    # matroid oracle
    cdef rank(self, X=*)
    cdef full_rank(self)
    cdef basis(self)
    cdef max_independent(self, X)
    cdef circuit(self, X=*)
    cdef fundamental_circuit(self, B, e)
    cdef closure(self, X)
    cdef k_closure(self, X, k)

    cdef augment(self, X, Y=*)

    cdef corank(self, X=*)
    cdef full_corank(self)
    cdef cobasis(self)
    cdef max_coindependent(self, X)
    cdef cocircuit(self, X=*)
    cdef fundamental_cocircuit(self, B, e)
    cdef coclosure(self, X)

    cdef loops(self)
    cdef is_independent(self, X)
    cdef is_dependent(self, X)
    cdef is_basis(self, X)
    cdef is_circuit(self, X)
    cdef is_closed(self, X)
    cdef is_subset_k_closed(self, X, int k)

    cdef coloops(self)
    cdef is_coindependent(self, X)
    cdef is_codependent(self, X)
    cdef is_cobasis(self, X)
    cdef is_cocircuit(self, X)
    cdef is_coclosed(self, X)

    # verification
    cdef is_valid(self)

    # enumeration
    cdef circuits(self)
    cdef nonspanning_circuits(self)
    cdef cocircuits(self)
    cdef noncospanning_cocircuits(self)
    cdef circuit_closures(self)
    cdef nonspanning_circuit_closures(self)
    cdef bases(self)
    cdef independent_sets(self)
    cdef independent_r_sets(self, long r)
    cdef nonbases(self)
    cdef dependent_r_sets(self, long r)
    cdef _extend_flags(self, flags)
    cdef _flags(self, r)
    cdef flats(self, r)
    cdef coflats(self, r)
    cdef hyperplanes(self)
    cdef f_vector(self)
    cdef broken_circuits(self, ordering=*)
    cdef no_broken_circuits_sets(self, ordering=*)

    # isomorphism
    cdef is_isomorphic(self, other, certificate=*)
    cdef _is_isomorphic(self, other, certificate=*)
    cdef isomorphism(self, other)
    cdef _isomorphism(self, other)
    cdef equals(self, other)
    cdef is_isomorphism(self, other, morphism)
    cdef _is_isomorphism(self, other, morphism)

    # minors, dual, truncation
    cdef minor(self, contractions=*, deletions=*)
    cdef contract(self, X)
    cdef delete(self, X)
    cdef _backslash_(self, X)
    cdef dual(self)
    cdef truncation(self)
    cdef has_minor(self, N, bint certificate=*)
    cdef has_line_minor(self, k, hyperlines=*, certificate=*)
    cdef _has_line_minor(self, k, hyperlines, certificate=*)

    # extension
    cdef extension(self, element=*, subsets=*)
    cdef coextension(self, element=*, subsets=*)
    cdef modular_cut(self, subsets)
    cdef linear_subclasses(self, line_length=*, subsets=*)
    cdef extensions(self, element=*, line_length=*, subsets=*)

    # connectivity
    cdef simplify(self)
    cdef cosimplify(self)
    cdef is_simple(self)
    cdef is_cosimple(self)
    cdef components(self)
    cdef is_connected(self, certificate=*)
    cdef connectivity(self, S, T=*)
    cdef _connectivity(self, S, T)
    cdef is_kconnected(self, k, certificate=*)
    cdef link(self, S, T)
    cdef _link(self, S, T)
    cdef _is_3connected_shifting(self, certificate=*)
    cdef _is_4connected_shifting(self, certificate=*)
    cdef _shifting_all(self, X, P_rows, P_cols, Q_rows, Q_cols, m)
    cdef _shifting(self, X, X_1, Y_2, X_2, Y_1, m)
    cdef is_3connected(self, certificate=*, algorithm=*)
    cdef is_4connected(self, certificate=*, algorithm=*)
    cdef _is_3connected_CE(self, certificate=*)
    cdef _is_3connected_BC(self, certificate=*)
    cdef _is_3connected_BC_recursion(self, basis, fund_cocircuits)

    # representability
    cdef _local_binary_matroid(self, basis=*)
    cdef binary_matroid(self, randomized_tests=*, verify=*)
    cdef is_binary(self, randomized_tests=*)
    cdef _local_ternary_matroid(self, basis=*)
    cdef ternary_matroid(self, randomized_tests=*, verify=*)
    cdef is_ternary(self, randomized_tests=*)

    # matroid k-closed
    cdef is_k_closed(self, int k)

    # matroid chordality
    cdef _is_circuit_chordal(self, frozenset C, bint certificate=*)
    cdef is_circuit_chordal(self, C, bint certificate=*)
    cdef is_chordal(self, k1=*, k2=*, bint certificate=*)
    cdef chordality(self)

    # optimization
    cdef max_weight_independent(self, X=*, weights=*)
    cdef max_weight_coindependent(self, X=*, weights=*)
    cdef is_max_weight_independent_generic(self, X=*, weights=*)
    cdef is_max_weight_coindependent_generic(self, X=*, weights=*)
    cdef intersection(self, other, weights=*)
    cdef _intersection(self, other, weights)
    cdef _intersection_augmentation(self, other, weights, Y)
    cdef intersection_unweighted(self, other)
    cdef _intersection_unweighted(self, other)
    cdef _intersection_augmentation_unweighted(self, other, Y)
    cdef partition(self)

    # invariants
    cdef _internal(self, B)
    cdef _external(self, B)
    cdef tutte_polynomial(self, x=*, y=*)
    cdef flat_cover(self, solver=*, verbose=*, integrality_tolerance=*)

    # misc
    cdef bergman_complex(self)
    cdef augmented_bergman_complex(self)

    # visualization
    cdef plot(self,B=*,lineorders=*,pos_method=*,pos_dict=*,save_pos=*)
    cdef show(self,B=*,lineorders=*,pos_method=*,pos_dict=*,save_pos=*,lims=*)
    cdef _fix_positions(self,pos_dict=*,lineorders=*)
