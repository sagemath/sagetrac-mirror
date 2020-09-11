#ifndef COMBINATORIAL_FACE_STRUCTURE
#define COMBINATORIAL_FACE_STRUCTURE

#include "bitsets.h"
/*
This structure is subject to change.
Use the access functions instead.
*/
struct face_struct{
    // All atoms that are contained in the face.
    uint64_t* atom_rep;

    // Coatoms such that their meet produces the face.
    // Note that in a simple polytope this is unique except for the empty face.
    uint64_t* coatom_gen;

    // Length of the above in terms of ``uint64_t``.
    size_t face_length;
    size_t face_length_coatom_gen;
};

enum algorithm_variant { standard, simple };

/*
#############################################################################
# Face Initalization
#############################################################################
*/

void face_clear(face_struct& face);

void face_copy(face_struct& dst, face_struct& src);

int initalize_face_with_allocate_instructions(\
        face_struct& face, size_t n_facets, size_t n_atoms, size_t step,
        void* location, size_t* alignment, size_t* size);

/*
#############################################################################
# Basic Access
#############################################################################
*/

int face_add_atom_safe(face_struct& face, size_t n);

/*
#############################################################################
# Underlying bitset functions access
#############################################################################
*/

inline int face_cmp(face_struct& a, face_struct& b){
    /*
    Return ``0`` if the faces are equal and consistently
    ``-1`` and ``1`` if not.
    */
    return bitset_cmp(a.atom_rep, b.atom_rep, a.face_length);
}

inline void face_add_atom(face_struct& face, size_t n){
    /*
    Add atom ``n`` to the face.

    This is unsafe and assumes enought memory.
    */
    bitset_add(face.atom_rep, n);
}

inline void face_discard_atom(face_struct& face, size_t n){
    bitset_discard(face.atom_rep, n);
}

inline void facet_set_coatom(face_struct& face, size_t n){
    /*
    Shall only be used for facets to specify their index.
    */
    bitset_clear(face.coatom_gen, face.face_length_coatom_gen);
    bitset_add(face.coatom_gen, n);
}

inline void face_set_first_n_atoms(face_struct& face, size_t n){
    bitset_set_first_n(face.atom_rep, face.face_length, n);
}

inline int atom_in_face(face_struct& face, int n){
    /*
    Return whether ``n`` is an atom of ``face``.
    */
    return bitset_in(face.atom_rep, n);
}

inline int is_zero(face_struct& face){
    /*
    Return whether ``face`` contains no coatoms.
    */
    return bitset_isempty(face.atom_rep, face.face_length);
}

inline size_t face_next_atom(face_struct& face, size_t n){
    /*
    Return the index of the next coatom in ``face`` with index >= ``n``.

    In case there are none, return ``-1``.
    */
    return bitset_next(face.atom_rep, face.face_length, n);
}

inline size_t count_atoms(face_struct& face){
    return count_atoms(face.atom_rep, face.face_length);
}

inline int test_alignment(face_struct& face){
    /*
    Return whether the data is correctly aligned.
    */
    return (test_alignment(face.atom_rep) && test_alignment(face.coatom_gen));
}

/*
#############################################################################
# Arithmetic
#############################################################################
*/

template <algorithm_variant N>
inline int is_subset(face_struct& A, face_struct& B){
    /*
    Return whether ``A`` is a subface of ``B``.
    */
    return is_subset(A.atom_rep, B.atom_rep, A.face_length);
}

template <>
inline int is_subset<simple>(face_struct& A, face_struct& B){
    // In the simple case we can use the coatom
    // generators.
    // Those are unique unless the face is empty.
    return is_subset(B.coatom_gen, A.coatom_gen, A.face_length_coatom_gen);
}

inline int is_subset(face_struct& A, face_struct& B){
    return is_subset<standard>(A,B);
}

template <algorithm_variant N>
inline void intersection(face_struct& dest, face_struct& A, face_struct& B){
    /*
    Set ``dest to the intersection of ``A`` and ``B``.
    */
    intersection(dest.atom_rep, A.atom_rep, B.atom_rep, dest.face_length);
}

template <>
inline void intersection<simple>(face_struct& dest, face_struct& A, face_struct& B){
    // In the simple case we must unite the coatom
    // generators as well.
    intersection(dest.atom_rep, A.atom_rep, B.atom_rep, dest.face_length);
    unite(dest.coatom_gen, A.coatom_gen, B.coatom_gen, dest.face_length_coatom_gen);
}

inline void intersection(face_struct& dest, face_struct& A, face_struct& B){
    intersection<standard>(dest, A, B);
}

#endif
