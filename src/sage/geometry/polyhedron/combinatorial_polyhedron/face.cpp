/*
This file defines the data structure for the face of a combinatorial polyhedron
and functions to access it.

#*****************************************************************************
#       Copyright (C) 2020 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
*/

#include <cstdint>
#include <cstdio>
#include "face.h"
using namespace std;

/*
#############################################################################
# Face Initalization
#############################################################################
*/

void face_clear(face_struct& face){
    bitset_clear(face.atom_rep, face.face_length);
    bitset_clear(face.coatom_gen, face.face_length_coatom_gen);
}

void face_copy(face_struct& dst, face_struct& src){
    /*
    Tries to copy as good as possible.
    Will set additional chunks to zero.
    Respect the maximal size of the destination.
    */
    bitset_copy(dst.atom_rep, src.atom_rep, dst.face_length, src.face_length);
    bitset_copy(dst.coatom_gen, src.coatom_gen, dst.face_length_coatom_gen, src.face_length_coatom_gen);
}

int initalize_face_with_allocate_instructions(\
        face_struct& face, size_t n_coatoms, size_t n_atoms, size_t step,
        void* location, size_t* alignment, size_t* size){
    /*
    Initialize a face step by step.

    ``n_coatoms`` and ``n_coatoms`` must much the corresponding numbers
    of Vrepresentatives and facets of the polyhedron.

    Return ``1`` as long as this function requests an aligned malloc.

    The caller must then allocate ``size`` memory to ``location[0]``,
    which needs to be ``alignment``-byte aligned.
    */
    switch (step) {
        case 0:
            face.face_length = get_face_length(n_atoms);
            face.face_length_coatom_gen = get_face_length(n_coatoms);
            size[0] = face.face_length*LIMB_BITS/8;
            alignment[0] = chunksize/8;
            return 1;
        case 1:
            // Assign the previous allocation instruction.
            face.atom_rep = (uint64_t*) location;

            size[0] = face.face_length_coatom_gen*LIMB_BITS/8;
            alignment[0] = chunksize/8;
            return 1;
        case 2:
            // Assign the previous allocation instruction.
            face.coatom_gen = (uint64_t*) location;
            return 0;
        default:
            return 0;
    }
}

/*
#############################################################################
# Basic Access
#############################################################################
*/

int face_add_atom_safe(face_struct& face, size_t n){
    /*
    Add atom `n` to the face.

    Return ``-1`` if the face is too small to do so.
    */
    if (n > face.face_length*LIMB_BITS)
        return -1;
    bitset_add(face.atom_rep, n);
    return 0;
}
