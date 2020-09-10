/*
#*****************************************************************************
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
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
#include <algorithm>
#include <cassert>
#include "face.h"
#include "face_list.h"
using namespace std;

void initialize_faces(face_list_struct& faces, size_t n_faces, size_t n_atoms){
    faces.n_faces = n_faces;
    faces.total_n_faces = n_faces;
    faces.n_atoms = n_atoms;
    faces.polyhedron_is_simple = 0;
}

void faces_copy(face_list_struct& dst, face_list_struct& src){
    /*
    Asserts that ``dst`` is allocated and fits enought faces.
    */
    face_struct* data = dst.faces;
    dst = src;
    dst.faces = data;
    for(size_t i=0; i<src.n_faces; i++){
        face_copy(dst.faces[i], src.faces[i]);
    }
}

void faces_shallow_copy(face_list_struct& dst, face_list_struct& src){
    /*
    Asserts that ``dst`` is allocated and fits enought faces.

    The faces itself are not copied, just the pointers set accordingly.

    Thus ``dest`` can be sorted and faces can be removed, but may not be altered.
    */
    face_struct* data = dst.faces;
    dst = src;
    dst.faces = data;
    for(size_t i=0; i<src.n_faces; i++){
        dst.faces[i] = src.faces[i];
    }
}

inline void intersection(face_list_struct& dest, face_list_struct& A, face_struct& B){
    /*
    Intersect any face in ``A`` with ``B`` to obtain ``dest``.
    */
    assert(dest.total_n_faces >= A.n_faces);
    assert(dest.n_atoms >= A.n_atoms);
    dest.n_faces = A.n_faces;

    int is_simple = A.polyhedron_is_simple;
    dest.polyhedron_is_simple = is_simple;

    for(size_t i=0; i < A.n_faces; i++){
        intersection(dest.faces[i], A.faces[i], B);
        dest.faces[i].coatom_gen_is_maximal = is_simple;
    }
}

inline int is_contained_in_one(face_struct& face, face_list_struct& faces){
    /*
    Return whether ``face`` is contained in one of ``faces``.
    */
    for(size_t i = 0; i < faces.n_faces; i++){
        if (is_subset(face, faces.faces[i]))
            return 1;
    }
    return 0;
}

inline int is_contained_in_one(face_struct& face, face_list_struct& faces, size_t skip){
    /*
    Return whether ``face`` is contained in one of ``faces``.

    Skips ``faces[skip]``.
    */
    face_list_struct faces_start = faces;
    face_list_struct faces_end = faces;

    faces_start.n_faces = skip;
    faces_end.n_faces -= skip+1;
    faces_end.faces += skip + 1;

    return is_contained_in_one(face, faces_start) || \
        is_contained_in_one(face, faces_end);
}

size_t get_next_level(\
        face_list_struct& faces, \
        face_list_struct& new_faces, \
        face_list_struct& visited_all){
    /*
    Set ``newfaces`` to be the facets of ``faces.faces[face.n_faces-1]``
    that are not contained in a face of ``visited_all``.

    Reduce the number of faces in ``faces`` by one.

    INPUT:

    - ``faces`` -- containing at least one face
    - ``new_faces`` -- needs to be of same size as ``faces``
    - ``visited_all`` -- the faces which have been visited before

    OUTPUT:

    - set ``newfaces`` to point to the new faces

    ALGORITHM:

    To get all facets of ``faces.faces[faces.n_faces-1]``, we would have to:
    - Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    - Add all the intersection of ``visited_all`` with the last face
    - Out of both the inclusion-maximal ones are of codimension one, i.e. facets.

    As we have visited all faces of ``visited_all``, we alter the algorithm
    to not revisit:
    Step 1: Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    Step 2: Out of those the inclusion-maximal ones are some of the facets.
            At least we obtain all of those, that we have not already visited.
            Maybe, we get some more.
    Step 3: Only keep those that we have not already visited.
            We obtain exactly the facets of ``faces[n_faces-1]`` that we have
            not visited yet.
    */

    // We keep track, which face in ``newfaces`` is a new face.
    size_t n_faces = faces.n_faces;
    int is_not_new_face[n_faces -1] = { 0 };

    // Step 1:
    faces.n_faces -= 1;
    intersection(new_faces, faces, faces.faces[n_faces-1]);

    for (size_t j = 0; j < n_faces-1; j++){
        if (new_faces.polyhedron_is_simple){
            // It suffices to check whether the face is non-empty and Step 3.
            if (is_zero(new_faces.faces[j]) ||
                    is_contained_in_one(new_faces.faces[j], visited_all))
                is_not_new_face[j] = 1;
        } else {
            // For each face we will do Step 2 and Step 3.
            if (is_contained_in_one(new_faces.faces[j], new_faces, j) || \
                    is_contained_in_one(new_faces.faces[j], visited_all))
                is_not_new_face[j] = 1;
        }
    }

    // Set ``new_faces`` to point to the correct ones.
    size_t n_new_faces = 0;
    for (size_t j = 0; j < n_faces -1; j++){
        if (is_not_new_face[j]) {
            // Not a new face of codimension 1.
            continue;
        }
        // It is a new face of codimension 1.
        // Either ``n_new_faces == j`` or ``new_faces.faces[n_new_faces]`` is not
        // a new face.
        swap(new_faces.faces[n_new_faces], new_faces.faces[j]);
        n_new_faces++;
    }
    new_faces.n_faces = n_new_faces;
    return n_new_faces;
}

size_t bit_rep_to_coatom_rep(face_struct& face, face_list_struct& coatoms, size_t *output){
    /*
    Write the coatom-representation of face in output. Return length.
    ``face_length`` is the length of ``face`` and ``coatoms[i]``
    in terms of uint64_t.
    ``n_coatoms`` length of ``coatoms``.
    */
    size_t count_length = 0;
    for (size_t i = 0; i < coatoms.n_faces; i++){
        if (is_subset(face, coatoms.faces[i])){
            // ``face`` is contain in ``coatoms[i]``,
            // then ``i`` is an element in the coatom-representation.
            output[count_length] = i;
            count_length++;
        }
    }
    return count_length;
}

int test_alignment(face_list_struct& faces){
   for(size_t i=0; i < faces.n_faces; i++){
       if (!test_alignment(faces.faces[i]))
           return 0;
   }
   return 1;
}
