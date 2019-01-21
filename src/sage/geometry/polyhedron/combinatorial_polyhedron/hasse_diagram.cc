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

#include <Python.h>
#include <math.h>
#include "hasse_diagram.h"
#include <cstdlib> //for aligned_alloc in C++17
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>


static uint64_t vertex_to_bit_dictionary[64];

//this dictionary helps storing a vector of 64 incidences as uint64_t, where each bit represents an incidence
void build_dictionary(){
    unsigned int i = 0;
    uint64_t count = 1;
    for (i=0; i< 64;i++){
        vertex_to_bit_dictionary[64-i-1] = count;
        count *= 2;
    }
}

//taken from https://codingforspeed.com/counting-the-number-of-leading-zeros-for-a-32-bit-integer-signed-or-unsigned/
//counts the number of leading zero bits of an uint64_t
inline unsigned int leading_zero_naive3(uint64_t x){
    unsigned n = 0;
    if (x == 0) return 64;
    while (1) {
        if (x > 0x7fffffffffffffff) break;
        n++;
        x <<= 1;
    }
    return n;
}

//counts the number of leading zero bits of a chunktype, where chunktype represents 1,2 or 4 uint64_t depending on the processor
inline unsigned int leading_zero_workaround(chunktype chunk){
    unsigned int i;
    unsigned int count = 0;
    uint64_t A[chunksize/64];
    store_register(A[0],chunk);
    for (i = 0;i < chunksize/64;i++){
        count += leading_zero_naive3(A[i]);
        if (count < 64*(i+1)){
            return count;
        }
    }
    return count;
}

//counts the number of trailing zero bits of an uint64_t
inline unsigned int trailing_zero_naive3(uint64_t x){
    unsigned n = 0;
    if (x == 0) return sizeof(x) * 8;
    while (1) {
        if (x % 2) break;
        n ++;
        x >>= 1;
    }
    return n;
}

//counts the number of trailing zero bits of a chunktype, where chunktype represents 1,2 or 4 uint64_t depending on the processor
inline unsigned int trailing_zero_workaround(chunktype chunk){
    unsigned int i;
    unsigned int count = 0;
    uint64_t A[chunksize/64];
    store_register(A[0],chunk);
    for (i = 0;i < chunksize/64;i++){
        count += trailing_zero_naive3(A[chunksize/64-i-1]);
        if (count < 64*(i+1)){
            return count;
        }
    }
    return count;
}


inline unsigned int naive_popcount(uint64_t A){
    unsigned int count = 0;
    while (A){
        count += A & 1;
        A >>= 1;
    }
    return count;
}


//initialization with a tuple of facets (each facet a tuple of vertices, vertices labeled 0,1,...)
CombinatorialPolyhedron::CombinatorialPolyhedron(unsigned int ** facets_pointer, unsigned int nr_facets_given, unsigned int *len_facets, unsigned int nr_vertices_given, int is_unbounded){
    unbounded = is_unbounded;
    if (is_unbounded)
        nr_lines = (unsigned int) is_unbounded - 1;
    else
        nr_lines = 0;
    build_dictionary();
    nr_vertices = nr_vertices_given;
    nr_facets = nr_facets_given;
    if ((!unbounded) && (nr_facets > nr_vertices)){//in this case the polar approach is actually better, so we will save the polar CombinatorialPolyhedron and compute accordingly
        //Polar_Init(py_tuple, nr_vertices_given);
        polar = 1;
        nr_vertices = nr_facets;
        nr_facets = nr_vertices_given;
    }
    get_facets_bitrep_from_facets_pointer(facets_pointer, len_facets);
    get_vertices_bitrep_from_facets_pointer(facets_pointer, len_facets);
}


//initialization with an incidence matrix given as tuple of tuples
CombinatorialPolyhedron::CombinatorialPolyhedron(unsigned int ** incidence_matrix, unsigned int nr_facets_given, unsigned int nr_vertices_given, int is_unbounded){
    unbounded = is_unbounded;
    if (is_unbounded)
        nr_lines = (unsigned int) is_unbounded - 1;
    else
        nr_lines = 0;
    build_dictionary();
    nr_vertices = nr_vertices_given;
    nr_facets = nr_facets_given;
    if ((!unbounded) && (nr_facets > nr_vertices)){//in this case the polar approach is actually much better, so we will save the polar CombinatorialPolyhedron and compute accordingly
        polar = 1;
        unsigned int nr_vertices_given = nr_vertices;
        nr_vertices = nr_facets;
        nr_facets = nr_vertices_given;
    }
    get_facets_from_incidence_matrix(incidence_matrix);
    get_vertices_from_incidence_matrix(incidence_matrix);
}

CombinatorialPolyhedron::~CombinatorialPolyhedron(){
    unsigned int i;
    deallocate_facets();
    deallocate_newfaces();
    deallocate_vertices();
    deallocate_allfaces();//must be called before deleting f_vector
    if (f_vector){
        delete[] f_vector;
    }
    if (edges){
        for (i=0; i< maxnumberedges; i++){
            if (edges[i]){
                delete[] edges[i];
            }
        }
        delete[] edges;
    }
    if (ridges){
        for (i=0; i< maxnumberedges; i++){
            if (ridges[i]){
                delete[] ridges[i];
            }
        }
        delete[] ridges;
    }
    if (incidences){
        for (i=0; i< maxnumberincidences; i++){
            if (incidences[i]){
                delete[] incidences[i];
            }
        }
        delete[] incidences;
    }
    face_iterator_destroy();
}

unsigned int CombinatorialPolyhedron::get_dimension(){
    if (dimension){
        return dimension;
    }
    dimension = calculate_dimension(facets, nr_facets);
    return dimension;
}


inline void CombinatorialPolyhedron::get_f_vector(unsigned long *vector){
    unsigned int i;
    get_dimension();//this assigns dimension the correct value
    if (!f_vector){
        get_f_vector_and_edges();
    }
    if (!polar){
        for (i = 0; i < dimension + 2; i++){
            vector[i] = f_vector[i];
        }
    }
    else {
        for (i = 0; i < dimension + 2; i++){
            vector[dimension + 1 - i] = f_vector[i];
        }
    }
}

inline unsigned int ** CombinatorialPolyhedron::get_edges(){
    if (polar){
        if (!nr_ridges){
            calculate_ridges();
        }
        return ridges;
    }
    else {
        if (!nr_edges){
            edgemode = 1;
            get_f_vector_and_edges();
        }
        return edges;
    }
}

inline unsigned int **  CombinatorialPolyhedron::get_ridges(){
    if (!polar){
        if (!nr_ridges){
            calculate_ridges();
        }
        return ridges;
    }
    else {
        if (!nr_edges){
            edgemode = 1;
            get_f_vector_and_edges();
        }
        return edges;
    }
}

//get_faces fills faces_to_return with all faces in dimension face_dimension and length_of_faces with length of the faces, if facet_repr then faces will be given with facet_incidences, otherwise with vertex incidences
void CombinatorialPolyhedron::get_faces(int face_dimension, unsigned int facet_repr, unsigned int **faces_to_return, unsigned int *length_of_faces){
    get_dimension();//this assigns dimension the correct value
    if (!f_vector){
            get_f_vector_and_edges();
    }
    if (f_vector[1] < nr_vertices){
        unbounded = 1;//in case unbounded was not passed properly, this is done to prevent crashing
    }
    unsigned int i;
    if (polar){ //if the polar is stored we should return the facet_repr for the vertex_repr and vice_versa
        facet_repr = polar - facet_repr;
        face_dimension = ((int) dimension) -1 - face_dimension;
    }
    if (face_dimension == -1){
        if (facet_repr) {
            for (i = 0; i < nr_facets; i++){
                faces_to_return[0][i] = i;
            }
            length_of_faces[0] = nr_facets;
            return;
        }
        else {
            length_of_faces[0] = 0;
            return;
        }
    }
    if (face_dimension < -1)
        return;
    unsigned int face_dimension_unsigned = (unsigned int) face_dimension;
    if (face_dimension_unsigned > dimension)
        return;
    if (face_dimension_unsigned == dimension -1){
        if (!facet_repr)
            return bitrep_to_list(facets, nr_facets, faces_to_return, length_of_faces, 0);
        else {
            for (i = 0; i < nr_facets; i++){
                faces_to_return[i][0] = i;
                length_of_faces[i] = 1;
            }
            return;
        }
    }
    if ((face_dimension_unsigned == 0) & (!unbounded)){
        if (nr_lines)
            return;
        if (facet_repr)
            return bitrep_to_list(vertices, nr_vertices, faces_to_return, length_of_faces, 1);
        else {
            for (i = 0; i < nr_vertices; i++){
                faces_to_return[i][0] = i;
                length_of_faces[i] = 1;
            }
            return;
        }
    }
    if (face_dimension_unsigned == dimension){
        if (!facet_repr) {
            for (i = 0; i < nr_vertices; i++){
                faces_to_return[0][i] = i;
            }
            length_of_faces[0] = nr_vertices;
            return;
        }
        else {
            length_of_faces[0] = 0;
            return;
        }
    }
    if (!allfaces_are_allocated || (allfaces_are_allocated[face_dimension_unsigned] != 2)){
        allocate_allfaces(face_dimension_unsigned);
        record_faces(face_dimension_unsigned);
    }
    if (!facet_repr)
        return bitrep_to_list(allfaces[face_dimension_unsigned], allfaces_counter[face_dimension_unsigned], faces_to_return, length_of_faces, 0);
    else
        return bitrep_to_list(allfaces_facet_repr[face_dimension_unsigned], allfaces_counter[face_dimension_unsigned], faces_to_return, length_of_faces, 1);
}

//initializes the face_iterator
void CombinatorialPolyhedron::face_iterator_init(int record_dimension, unsigned int vertex_repr, unsigned int facet_repr){
    if (face_iterator_is_initialized){
        face_iterator_destroy();//in case there is some face_iterator around, that wasn't used up
    }
    const unsigned int dim = get_dimension();
    allocate_newfaces();
    face_iterator_is_initialized = 1;
    face_iterator_current_dimension = dimension - 1;
    if (polar){
        if (record_dimension < 0){
            face_iterator_record_dimension = record_dimension;
        }
        else {
            face_iterator_record_dimension = ((int) dimension) - 1 - record_dimension;
        }
        face_iterator_vertex_repr = facet_repr;
        face_iterator_facet_repr = vertex_repr;
    }
    else {
        face_iterator_record_dimension = record_dimension;
        face_iterator_vertex_repr = vertex_repr;
        face_iterator_facet_repr = facet_repr;
    }

    //in case we have already recorded all faces in the given dimension, we might save time, just returning the stored lists
    if ((0 < face_iterator_record_dimension) && (face_iterator_record_dimension < (int) dimension - 1) && allfaces_are_allocated && allfaces_are_allocated[face_iterator_record_dimension]){
        face_iterator_counter = 0;
        face_iterator_is_initialized = 2;
        return;
    }

    face_iterator_yet_to_yield = nr_facets;
    face_iterator_nr_faces = new unsigned int[dim]();
    face_iterator_nr_faces[dimension -1] = nr_facets;
    face_iterator_nr_forbidden = new unsigned int[dim]();
    face_iterator_nr_forbidden[dimension -1] = 0;
    face_iterator_first_time = new unsigned int[dim]();
    face_iterator_first_time[dimension - 1] = 1;
}


//returns one face at a time of all faces of record_dimension dimension
//returns all faces if dimension == -2
//returns not the faces in dimesion ``-1`` and dimesion ``dimension``
inline unsigned int CombinatorialPolyhedron::face_iterator_call(unsigned int *Vface_to_return, unsigned int *Vlength, unsigned int *Hface_to_return, unsigned int *Hlength){
    unsigned int return_length = 0;
    if (!face_iterator_is_initialized)
        return 0;

    if (face_iterator_is_initialized == 2){
        if (face_iterator_counter < allfaces_counter[face_iterator_record_dimension]){
            if (!polar){
                if (face_iterator_vertex_repr)
                    bitrep_to_list(allfaces[face_iterator_record_dimension][face_iterator_counter], Vface_to_return, Vlength, 0);
                if (face_iterator_facet_repr)
                    bitrep_to_list(allfaces_facet_repr[face_iterator_record_dimension][face_iterator_counter], Hface_to_return, Hlength, 1);
                face_iterator_counter++;
                return 1;
            }
            else {
                if (face_iterator_vertex_repr)
                    bitrep_to_list(allfaces[face_iterator_record_dimension][face_iterator_counter], Hface_to_return, Hlength, 0);
                if (face_iterator_facet_repr)
                    bitrep_to_list(allfaces_facet_repr[face_iterator_record_dimension][face_iterator_counter], Vface_to_return, Vlength, 1);
                face_iterator_counter++;
                return 1;
            }
        }
        return 0;
    }

    while ((!return_length) && (face_iterator_current_dimension != dimension)){
        if (polar)
            return_length = face_iterator(Hface_to_return, Hlength, Vface_to_return, Vlength);
        else
            return_length = face_iterator(Vface_to_return, Vlength, Hface_to_return, Hlength);
    }


    //in this case the iterator is used up and we should clean up
    if (!return_length) {
        face_iterator_destroy();
    }
    return return_length;
}

void CombinatorialPolyhedron::face_iterator_destroy(){
    face_iterator_is_initialized = 0;
    if (face_iterator_nr_faces)
        delete[] face_iterator_nr_faces;
    if (face_iterator_nr_forbidden)
        delete[] face_iterator_nr_forbidden;
    if (face_iterator_first_time)
        delete[] face_iterator_first_time;
    face_iterator_nr_faces = NULL;
    face_iterator_nr_forbidden = NULL;
    face_iterator_first_time = NULL;
}

void CombinatorialPolyhedron::record_all_faces(){
    allocate_allfaces();
    if (unbounded){
        record_faces(0);
    }
    else {
        record_faces(1);
    }
}

unsigned long ** CombinatorialPolyhedron::get_incidences(int dimension_one, int dimension_two, unsigned long * nr_incidences_to_return, unsigned int * twisted){
    nr_incidences = 0;
    get_dimension();//this assigns dimension the correct value
    if (!f_vector)
        get_f_vector_and_edges();
    unsigned long i,j;
    if (polar){
        dimension_one = dimension -1 - dimension_one;
        dimension_two = dimension -1  - dimension_two;
    }
    twisted[0] = 0;
    if (nr_lines && ((dimension_one == 0) || (dimension_two == 0))){//if the unbounded polyhedron contains lines, then there are no vertices
        nr_incidences_to_return[0] = nr_incidences;
        return incidences;
    }
    if (dimension_one < dimension_two){
        int foo = dimension_one;
        dimension_one = dimension_two;
        dimension_two = foo;
        twisted[0] = 1;
    }
    if (dimension_two == dimension_one){
        unsigned long nr_of_those_faces = f_vector[dimension_two + 1];
        for (j = 0; j < nr_of_those_faces; j++){
            add_incidence(j,j);
        }
        nr_incidences_to_return[0] = nr_incidences;
        return incidences;
    }
    unsigned int one = (unsigned int) dimension_one;//if dimension_one < 0, then dimension_two < -1 and the next statement will be triggered anyway)
    if ((dimension_two < -1) || (one > dimension)){//in this case there will be no incidences
        nr_incidences_to_return[0] = nr_incidences;
        return incidences;
    }
    if (dimension_two == -1){//the -1-dimensional face is contained in every face
        unsigned long nr_of_those_faces = f_vector[dimension_one + 1];
        for (j = 0; j < nr_of_those_faces; j++){
            add_incidence(j,0);
        }
        nr_incidences_to_return[0] = nr_incidences;
        return incidences;
    }
    unsigned int two = (unsigned int) dimension_two;
    if (one == dimension){ //the polytope contains every face
        unsigned long nr_of_those_faces = f_vector[dimension_two + 1];
        for (j = 0; j < nr_of_those_faces; j++){
            add_incidence(0,j);
        }
        nr_incidences_to_return[0] = nr_incidences;
        return incidences;
    }
    if (f_vector[1] != nr_vertices){
        unbounded = 1;
    }
    if ((one == dimension -1) && (two == 0) && (!unbounded)){//getting the vertex-facet incidence, in the unbounded case we cannot assume the vertices to be actual vertices
        vertex_facet_incidences();
        nr_incidences_to_return[0] = nr_incidences;
        return incidences;
    }
    if (one == dimension -1){
        if (!allfaces_are_allocated || (allfaces_are_allocated[two] != 2)){
            allocate_allfaces(two);
            record_faces(two);
        }
        for (i = 0; i < nr_facets; i++)
            for (j = 0; j < allfaces_counter[two]; j++)
                if (is_subset(allfaces[two][j],facets[i]))
                    add_incidence(i,j);
        nr_incidences_to_return[0] = nr_incidences;
        return incidences;
    }
    if ((two == 0) && (!unbounded)){//in the unbounded case, we need to figure out, what the vertices actually are
        if (!allfaces_are_allocated || (allfaces_are_allocated[one] != 2)){
            allocate_allfaces(one);
            record_faces(one);
        }
        for (i = 0; i < allfaces_counter[one]; i++)
            for (j = 0; j < nr_vertices; j++)
                if (is_subset_facet_repr(allfaces_facet_repr[one][i],vertices[j]))
                    add_incidence(i,j);
        nr_incidences_to_return[0] = nr_incidences;
        return incidences;
    }
    //at this point we know 0 < two < one < dimension - 1 (in the unbounded case 0<= two)
    if (!allfaces_are_allocated || (allfaces_are_allocated[one] != 2) || (allfaces_are_allocated[two] != 2)){
        allocate_allfaces(one);
        allocate_allfaces(two);
        record_faces(two);
    }
    for (i = 0; i < allfaces_counter[one]; i++)
        for (j = 0; j < allfaces_counter[two]; j++)
            if (is_subset_facet_repr(allfaces_facet_repr[one][i],allfaces_facet_repr[two][j]))
                add_incidence(i,j);
    nr_incidences_to_return[0] = nr_incidences;
    return incidences;
}

//expects a sorted array of size lenght with pairwise distinct integers at least 0 and at most dim
inline unsigned long CombinatorialPolyhedron::get_flag_number_init(unsigned int *flagarray, unsigned int length){
    unsigned int i;
    unsigned int counter = get_dimension()-1;
    if (!polar){
        if (nr_lines > flagarray[0])
            return 0;
    }
    else {//turning the array around
        for (i=0;i< (length + 1)/2;i++){
            unsigned int saver = flagarray[i];
            flagarray[i] = dimension - 1 - flagarray[length - 1 - i]; //in case i = (length + 1)/2 - 1 both might be the same, which works fine
            flagarray[length - 1 - i] = dimension - 1 - saver;
        }
    }
    for(i=0;i<length;i++){
    if ((!allfaces_are_allocated) || (allfaces_are_allocated[flagarray[i]] != 2)){
        if (flagarray[i] < counter){
        counter = flagarray[i];
        }
            allocate_allfaces(flagarray[i]);
        }
    }
    record_faces(counter);
    return get_flag_number(flagarray,length);
}


//will set C to be the intersection of A and B
inline void CombinatorialPolyhedron::intersection(chunktype *A, chunktype *B, chunktype *C){
    unsigned int i;
    for (i = 0; i < length_of_face; i++){
        C[i] = bitwise_intersection(A[i],B[i]);
    }
}

//returns 1 if A is a proper subset of B, otherwise returns 0, this is done by checking if there is an element in A, which is not in B
inline int CombinatorialPolyhedron::is_subset(chunktype *A, chunktype *B){
    unsigned int i;
    for (i = 0; i < length_of_face; i++){
        if (bitwise_is_not_subset(A[i],B[i])){
            return 0;
        }
    }
    return 1;
}

//as above just in facet_repr
inline int CombinatorialPolyhedron::is_subset_facet_repr(chunktype *A, chunktype *B){
    unsigned int i;
    for (i = 0; i < length_of_face_in_facet_repr; i++){
        if (bitwise_is_not_subset(A[i],B[i])){
            return 0;
        }
    }
    return 1;
}

//counts the number of vertices in a face by counting bits set to one
inline unsigned int CombinatorialPolyhedron::CountFaceBits(chunktype* A1) {
    //this function is not implemented for speed (it basically gets called dimension times and once to convert a face to a tuple
    unsigned int i,count = 0;
    const unsigned int length_of_conversion_face = length_of_face*chunksize/64;
    unsigned long A[length_of_conversion_face];
    for (i=0;i<length_of_face;i++){
        store_register(A[i*chunksize/64],A1[i]);
    }
    for (i=0;i<length_of_conversion_face;i++){
        count += popcount(A[i]);
    }
    return count;
}

inline unsigned int CombinatorialPolyhedron::CountFaceBits_facet_repr(chunktype* A1) {
    //this function is not implemented for speed (it basically gets called dimension times and once to convert a face to a tuple
    unsigned int i,count = 0;
    const unsigned int length_of_conversion_face = length_of_face_in_facet_repr*chunksize/64;
    unsigned long A[length_of_conversion_face];
    for (i=0;i<length_of_face_in_facet_repr;i++){
    store_register(A[i*chunksize/64],A1[i]);
    }
    for (i=0;i<length_of_conversion_face;i++){
    count += popcount(A[i]);
    }
    return count;
}

//adds an edge to the edges list
inline void CombinatorialPolyhedron::add_edge(chunktype *face){
    unsigned int i,one = 0,two = 0;
    for (i = 0; i < length_of_face; i++){
        one += leading_zero_count(face[i]);
        if (one < (i+1)*chunksize){
            break;
        }
    }
    for (i = 0; i < length_of_face; i++){
        two += trailing_zero_count(face[length_of_face-i-1]);
        if (two < (i+1)*chunksize){
            break;
        }
    }
    add_edge(one,length_of_face*chunksize - two - 1);
}

//adds an edge to the edges list given as its two vertices
inline void CombinatorialPolyhedron::add_edge(unsigned int one, unsigned int two){
            if (nr_edges >= maxnumberedges*maxnumberedges){
                return;
            }
            unsigned int position_one = nr_edges / maxnumberedges;
            unsigned int position_two = 2*(nr_edges % maxnumberedges);
            if (!position_two){
                edges[position_one] = new unsigned int [maxnumberedges*2];
            }
            edges[position_one][position_two] = one;
            edges[position_one][position_two + 1] = two;
            nr_edges += 1;
        }

//adds a ridge to the ridge list given as its two facets
inline void CombinatorialPolyhedron::add_ridge(unsigned int one, unsigned int two){
    if (nr_ridges >= maxnumberedges*maxnumberedges){
        return;
    }
    unsigned int position_one = nr_ridges / maxnumberedges;
    unsigned int position_two = 2*(nr_ridges % maxnumberedges);
    if (!position_two){
        ridges[position_one] = new unsigned int [maxnumberedges*2];
    }
    ridges[position_one][position_two] = one;
    ridges[position_one][position_two + 1] = two;
    nr_ridges += 1;
}

//adds an incidence to the list of incidences, where one and two correspond to the number of the faces according to allfaces resp. vertices/facets
inline void CombinatorialPolyhedron::add_incidence(unsigned long one, unsigned long two){
    if (nr_incidences >= maxnumberincidences*maxnumberincidences){
        return;
    }
    unsigned int position_one = nr_incidences / maxnumberincidences;
    unsigned int position_two = 2*(nr_incidences % maxnumberincidences);
    if (!incidences[position_one]){
        incidences[position_one] = new unsigned long [maxnumberincidences*2];
    }
    incidences[position_one][position_two] = one;
    incidences[position_one][position_two + 1] = two;
    nr_incidences += 1;
}

//intersects the first 'lenfaces' faces of 'faces' with the 'face_to_intersect'-th face of faces and stores the result in 'nextfaces'
//then determines which ones are exactly of one dimension less by considering containment
//newfaces2 will point at those of exactly one dimension less which are not contained in any of the faces in 'forbidden'
//returns the number of those faces
inline unsigned int CombinatorialPolyhedron::get_next_level(chunktype **faces, unsigned int lenfaces, unsigned int face_to_intersect, chunktype **nextfaces, chunktype **nextfaces2, unsigned int nr_forbidden){
    unsigned int j,k, addthisface;
    unsigned int newfacescounter = 0;
    unsigned int counter = 0;
    for (j = 0; j < lenfaces; j++){
        if (j != face_to_intersect){
            intersection(faces[j],faces[face_to_intersect],nextfaces[counter]);
            counter ++;
        }
    }
    //we have create all possible intersection with the i_th-face, but some of them might not be of exactly one dimension less
    const unsigned int constlenfaces = lenfaces;
    int addfacearray[constlenfaces] = { };
    for (j = 0; j< lenfaces-1; j++){
    addfacearray[j] = 1;
            for(k = 0; k < j; k++){
        if(is_subset(nextfaces[j],nextfaces[k])){
                    addfacearray[j] = 0;
                    break;
                }
            }
            if (!addfacearray[j]) {
                continue;
            }
    for(k = j+1; k < lenfaces-1; k++){
        if(is_subset(nextfaces[j],nextfaces[k])){
        addfacearray[j] = 0;
        break;
        }
    }
    if (!addfacearray[j]) {
        continue;
    }
    for (k = 0; k < nr_forbidden; k++){//we do not want to double count any faces, the faces in forbidden we have considered already
                if(is_subset(nextfaces[j],forbidden[k])){
                    addfacearray[j] = 0;
                    break;
                }
            }
    }
    for (j = 0; j < lenfaces -1; j++){
            if (!addfacearray[j]) {
                continue;
            }
            nextfaces2[newfacescounter] = nextfaces[j];
            newfacescounter++;
        }
    return newfacescounter;
}

unsigned int CombinatorialPolyhedron::calculate_dimension(chunktype **faces, unsigned int nr_faces){
    unsigned int i,j,k, newfacescounter, dim;
    if (nr_faces == 0){
        return 0;//this isn't supposed to happen, but maybe the data is malformed
    }
    unsigned int bitcount = CountFaceBits(faces[0]);
    if (bitcount == 0){//if a polyhedron contains only the empty face as facet, then it is of dimension 0
        return 0;
    }
    if (nr_faces == 1){//if there is only one facet and this contains bitcount # of vertices/rays/lines then the polyhedron is of dimension bitcount -1
        return bitcount;
    }
    if (bitcount == 1){
        return 1;
    }
    void *nextfaces_creator[nr_faces-1];
    chunktype  *nextfaces2[nr_faces-1], *nextfaces[nr_faces-1];
    for (i=0; i < (nr_faces-1); i++){
        nextfaces_creator[i] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
        nextfaces[i] = (chunktype*) nextfaces_creator[i];
    }
    newfacescounter = get_next_level(faces,nr_faces,nr_faces-1,nextfaces,nextfaces2,0);//calculates the ridges contained in one facet
    dim =  calculate_dimension(nextfaces2,newfacescounter) + 1;//calculates the dimension of that facet
    if (dim == 1){
        dim = bitcount;//our face should be a somewhat a vertex, but if the polyhedron is unbounded, than our face will have dimension equal to the number of 'vertices' it contains, where some of the vertices might represent lines
    }
    for (i=0; i < (nr_faces - 1); i++){
        free(nextfaces_creator[i]);
    }
    return dim;
}

void CombinatorialPolyhedron::calculate_ridges(){//this is a much simpler version of belows get_f_vector_and_edges
    if (nr_facets <= 1){
        return;
    }
    if ((nr_lines > 0) && (dimension - 1 == nr_lines)){
        return;//in this case it might look like there is a ridge, but actually this is not a face but a line
    }
    unsigned int i,j,counter, addthisface, nr_forbidden = 0;
    unsigned long newfacescounter;
    void *nextfaces_creator[nr_facets-1];
    const unsigned int const_facets = nr_facets;
    forbidden = new chunktype *[const_facets];
    chunktype  *nextfaces2[nr_facets-1], *nextfaces[nr_facets-1];
    for (i=0; i < (nr_facets-1); i++){
        nextfaces_creator[i] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
        nextfaces[i] = (chunktype*) nextfaces_creator[i];
    }
    i = nr_facets;
    while (i--){
        newfacescounter = get_next_level(facets,i+1,i,nextfaces,nextfaces2,nr_forbidden);
        if (newfacescounter){
            counter = 0;
            for (j = 0; j < newfacescounter; j++){
                while (nextfaces[j] != nextfaces2[counter]){
                    counter++;
                }
                add_ridge(j,i);
            }
        }
        forbidden[nr_forbidden] = facets[i];
        nr_forbidden++;
    }
    for (i=0; i < (nr_facets - 1); i++){
        free(nextfaces_creator[i]);
    }
}


void CombinatorialPolyhedron::get_f_vector_and_edges(){
    unsigned int i;
    if (!dimension){
        get_dimension();
    }
    allocate_newfaces();
    if (!f_vector){
        const unsigned int const_dimension = dimension;
        f_vector = new unsigned long[const_dimension + 2]();
    }
    else
    {
        for (i = 1; i < dimension; i++)
            f_vector[i] = 0;
    }
    f_vector[0] = 1;
    f_vector[dimension+1] = 1;
    f_vector[dimension] = nr_facets;
    //f_vector[1] = nr_vertices; //this is calculated explicetly in order to make calculations also work for unbounded polyhedra
    get_f_vector_and_edges(facets,dimension-1,nr_facets,0);
}

//will be called from the initial get_f_vector_and_edges()
void CombinatorialPolyhedron::get_f_vector_and_edges(chunktype **faces, unsigned int dim, unsigned int nr_faces, unsigned int nr_forbidden){
    unsigned int i;
    unsigned long newfacescounter;
    if (dim == 1){
        if (edgemode)//in this case we want to record the edges
        for (i = 0; i < nr_faces; i++){
            add_edge(faces[i]);
        }
    }
    if (nr_faces <= 1){
        return;//there will be no newfaces, in the case of nr_facets == 1, newfaces might not have been initialized properly, however with correct data, one of the above things should happen
    }
    if (dim == nr_lines){
        return;//it might look like faces contain a common vertex, but by user input we know that this is not a vertex but a line, i.e. not a face
    }
    i = nr_faces;
    while (i--){
        newfacescounter = get_next_level(faces,i+1,i,newfaces[dim-1],newfaces2[dim-1],nr_forbidden);//get the facets contained in faces[i] but not in any of the forbidden
        f_vector[dim] += newfacescounter;
        if (newfacescounter){
            get_f_vector_and_edges(newfaces2[dim-1],dim-1,newfacescounter,nr_forbidden);//add all face in faces[i] to the f_vector, but not those which we have counted already
        }
        forbidden[nr_forbidden] = faces[i];//we have counted all faces in faces[i], so we do not want to count them ever again
        nr_forbidden++;
    }
}

void CombinatorialPolyhedron::record_faces(unsigned int lowest_dimension){
    if (nr_lines > lowest_dimension){
        lowest_dimension = nr_lines;
    }
    record_faces(facets,dimension-1,nr_facets,0,lowest_dimension);
    for (unsigned int i = lowest_dimension; i < dimension - 1; i++)
        if (allfaces_are_allocated[i])
            allfaces_are_allocated[i] = 2;//making sure that we record the faces of each dimension at most once
}

void CombinatorialPolyhedron::record_faces(chunktype **faces, unsigned int current_dimension, unsigned int nr_faces, unsigned int nr_forbidden, unsigned int lowest_dimension){
    unsigned int i;
    unsigned long newfacescounter;
    if ((current_dimension < dimension - 1) && (allfaces_are_allocated[current_dimension] == 1)){
        for (i = 0; i < nr_faces; i++){
            record_face(faces[i], current_dimension);
        }
    }
    if (current_dimension == lowest_dimension){
        return;
    }
    if (current_dimension == nr_lines){
        return;//it might look like faces contain a common vertex, but by user input we know that this is not a vertex but a line, i.e. not a face
    }
    i = nr_faces;
    while (i--){
        newfacescounter = get_next_level(faces,i+1,i,newfaces[current_dimension-1],newfaces2[current_dimension-1],nr_forbidden);//get the facets contained in faces[i] but not in any of the forbidden
        if (newfacescounter){
            record_faces(newfaces2[current_dimension-1],current_dimension-1,(unsigned int) newfacescounter,nr_forbidden, lowest_dimension);//visit all faces in  faces[i], but not those which we visited already
        }
        forbidden[nr_forbidden] = faces[i];//we have visited all faces in faces[i], so we do not want to visit them again
        nr_forbidden++;
    }
}

inline void CombinatorialPolyhedron::record_face(chunktype *face, unsigned int current_dimension){
    for (unsigned int j = 0; j < length_of_face; j++){
        load_register(allfaces[current_dimension][allfaces_counter[current_dimension]][j],face[j]);
    }
    record_face_facet_repr(face, current_dimension);
    allfaces_counter[current_dimension]++;
}

inline void CombinatorialPolyhedron::record_face_facet_repr(chunktype *face, unsigned int current_dimension){
    unsigned int i;
    unsigned int position, value;
    const unsigned int size_array = length_of_face_in_facet_repr*chunksize/64;
    uint64_t *array = new uint64_t [size_array]();
    for (i = 0; i < nr_facets; i++){
        if (is_subset(face, facets[i])){
            value = i % 64;
            position = i/64;
            array[position] += vertex_to_bit_dictionary[value];
        }
    }
    for (i=0;i<length_of_face_in_facet_repr;i++){
        load_register(allfaces_facet_repr[current_dimension][allfaces_counter[current_dimension]][i],array[i*chunksize/64]);
    }
    delete[] array;
}

//copied mostly from record_faces
inline unsigned int CombinatorialPolyhedron::face_iterator(unsigned int *Vface_to_return, unsigned int *Vlength, unsigned int *Hface_to_return, unsigned int *Hlength){
    unsigned int current_dimension = face_iterator_current_dimension;
    if (current_dimension == dimension){
        return 0;//the function is not supposed to be called in this case
    }
    unsigned int nr_faces = face_iterator_nr_faces[current_dimension];
    unsigned int nr_forbidden = face_iterator_nr_forbidden[current_dimension];
    chunktype **faces;
    if (current_dimension == dimension -1)
        faces = facets;
    else
        faces = newfaces2[current_dimension];
    unsigned int i;
    unsigned long newfacescounter;
    if ((face_iterator_record_dimension > -2) && (face_iterator_record_dimension != (int) current_dimension))
        face_iterator_yet_to_yield = 0;
    if (face_iterator_yet_to_yield > 0){
        if (face_iterator_vertex_repr){
            bitrep_to_list(faces[face_iterator_yet_to_yield - 1], Vface_to_return, Vlength, 0);
        }
        if (face_iterator_facet_repr){
            unsigned int counter = 0;
            for (i = 0; i < nr_facets; i++){
                if (is_subset(faces[face_iterator_yet_to_yield - 1], facets[i])){
                    Hface_to_return[counter] = i;
                    counter++;
                }
            }
            Hlength[0] = counter;
        }
        face_iterator_yet_to_yield--;
        return 1;
    }
    if ((int) current_dimension <= face_iterator_record_dimension){
        face_iterator_current_dimension++;
        return 0;
    }
    if (current_dimension == 0){
        face_iterator_current_dimension++;
        return 0;
    }
    if (nr_faces <= 1){
        face_iterator_current_dimension++;
        return 0;
    }
    if (current_dimension == nr_lines){
        face_iterator_current_dimension++;
        return 0;//it might look like faces contain a common vertex, but by user input we know that this is not a vertex but a line, i.e. not a face
    }
    i = nr_faces - 1;
    face_iterator_nr_faces[current_dimension]--;
    if (!face_iterator_first_time[current_dimension]){ //this is necessary, as we cannot at the new forbidden face at the end of the function
        forbidden[nr_forbidden] = faces[i+1];
        face_iterator_nr_forbidden[current_dimension]++;
        nr_forbidden = face_iterator_nr_forbidden[current_dimension];
    }
    else {
        face_iterator_first_time[current_dimension] = 0;
    }
    newfacescounter = get_next_level(faces,i+1,i,newfaces[current_dimension-1],newfaces2[current_dimension-1],nr_forbidden);//get the facets contained in faces[i] but not in any of the forbidden
    if (newfacescounter){
        face_iterator_first_time[current_dimension - 1] = 1;
        face_iterator_nr_faces[current_dimension - 1] = (unsigned int) newfacescounter;//newfacescounter is a small number, I had it be a long in order to fit addition to the f_vector
        face_iterator_nr_forbidden[current_dimension - 1] = nr_forbidden;
        face_iterator_yet_to_yield = (unsigned int) newfacescounter;
        face_iterator_current_dimension--;
        return 0;
    }
    return 0;
}

void CombinatorialPolyhedron::vertex_facet_incidences(){
    for (unsigned int i = 0; i < nr_facets; i++)
        vertex_facet_incidences(facets[i],i);
}

void CombinatorialPolyhedron::vertex_facet_incidences(chunktype *array1, unsigned int nr_facet){
    unsigned int i,j;
    unsigned int counter = 0;
    const unsigned int size_array = length_of_face*chunksize/64;
    uint64_t *array = new uint64_t [size_array]();
    for (i = 0; i < length_of_face;i++){
        store_register(array[i*chunksize/64],array1[i]);
    }
    for (i = 0; i < size_array;i++){
        for (j = 0; j < 64; j++){
            if (array[i] >= vertex_to_bit_dictionary[j]){
                add_incidence(nr_facet, i*64+j);
                counter++;
                array[i] -= vertex_to_bit_dictionary[j];
            }
        }
    }
}

unsigned long CombinatorialPolyhedron::get_flag_number(unsigned int *array, unsigned int len){
    if (len == 0){
        return 0;
    }
    if (!f_vector)
        get_f_vector_and_edges();
    if (len == 1){
        return f_vector[array[0]+1];
    }
    unsigned long i,j,k;
    const unsigned int const_len = len;
    unsigned long **saverarray = new unsigned long *[const_len];
    for (i= 0; i < len; i++){
        const unsigned int len_array = f_vector[array[i]+1];
        saverarray[i] = new unsigned long [len_array]();
    }
    unsigned int counter = 1;
    if (array[0] == 0){
        if (array[1] < dimension - 1){
            if (len > 2){
        for (i=0; i < f_vector[array[1]+1];i++){
            saverarray[1][i] = CountFaceBits(allfaces[array[1]][i]);
        }
        }
        else{
        unsigned long sum = 0;
        for (i= 0;i < f_vector[array[1]+1]; i++){
            sum += CountFaceBits(allfaces[array[1]][i]);
        }
        delete[] saverarray;
        return sum;
        }
        }
        else {
            unsigned long sum = 0;
            for (i = 0; i < nr_facets; i++){
                sum += CountFaceBits(facets[i]);
            }
            delete[] saverarray;
            return sum;
        }
        counter = 2;
    }
    else {
        for (i= 0; i < f_vector[array[0]+1]; i++){
            saverarray[0][i] = 1;
        }
    }
    while (counter < len -1){
        for (j= 0; j < f_vector[array[counter-1]+1]; j ++){
            for (k = 0; k < f_vector[array[counter]+1]; k++){
                if (is_subset_facet_repr(allfaces_facet_repr[array[counter]][k],allfaces_facet_repr[array[counter-1]][j])){
                    saverarray[counter][k] += saverarray[counter -1][j];
                }
            }
        }
        counter++;
    }
    unsigned long sum = 0;
    if (array[counter] == dimension -1){
        for (i = 0; i < f_vector[array[counter-1]+1]; i++){
            sum += saverarray[counter - 1][i]*CountFaceBits_facet_repr(allfaces_facet_repr[array[counter-1]][i]);
        }
    }
    else {
        for (j= 0; j < f_vector[array[counter-1]+1]; j ++)
            for (k = 0; k < f_vector[array[counter]+1]; k++){
                if (is_subset_facet_repr(allfaces_facet_repr[array[counter]][k],allfaces_facet_repr[array[counter-1]][j])){
                    sum += saverarray[counter -1][j];
                }
            }
    }
    delete[] saverarray;
    return sum;
}


//initialization

void CombinatorialPolyhedron::get_facets_bitrep_from_facets_pointer(unsigned int ** facets_pointer, unsigned int *len_facets){
    length_of_face = ((nr_vertices - 1)/chunksize + 1);//this determines the length of the face in terms of chunktype
    allocate_facets();
    if (polar){
        get_vertices_or_facets_bitrep_from_facets_pointer(facets_pointer, len_facets, facets, 1, nr_facets, nr_vertices, 0);
    }
    else {
        get_vertices_or_facets_bitrep_from_facets_pointer(facets_pointer, len_facets, facets, 0, nr_vertices, nr_facets, 0);
    }
}

void CombinatorialPolyhedron::get_vertices_bitrep_from_facets_pointer(unsigned int ** facets_pointer, unsigned int *len_facets){
    length_of_face_in_facet_repr = ((nr_facets - 1)/chunksize + 1);//this determines the length in facet representation in terms of chunktype
    allocate_vertices();
    if (!polar){
        get_vertices_or_facets_bitrep_from_facets_pointer(facets_pointer, len_facets, vertices, 1, nr_vertices, nr_facets, 1);
    }
    else {
        get_vertices_or_facets_bitrep_from_facets_pointer(facets_pointer, len_facets, vertices, 0, nr_facets, nr_vertices, 1);
    }
}

void CombinatorialPolyhedron::get_vertices_or_facets_bitrep_from_facets_pointer(unsigned int ** facets_pointer, unsigned int *len_facets, chunktype** facets_or_vertices, unsigned int flip, unsigned int nr_vertices_given, unsigned int nr_facets_given, unsigned int facet_repr){
    unsigned int i,j;
    if (flip){//getting the vertices from the original polytope, those will be facets or vertices depending on wether we consider polar or not
        const int size_three = nr_facets_given;
        unsigned int *old_facets_walker = new unsigned int [size_three]();
        unsigned int new_facets_array[size_three];
        unsigned int length_that_face;
        for(i = 0;i<nr_vertices_given;i++){
            length_that_face = 0;
            for (j=0;j < nr_facets_given;j++){
                if (i == facets_pointer[j][old_facets_walker[j]]){//testing if vertex i is contained in the j-th facet
                    new_facets_array[length_that_face] = j;
                    old_facets_walker[j]++;
                    length_that_face++;
                    if (old_facets_walker[j] >= len_facets[j])
                        old_facets_walker[j]--;
                }
            }
            char_from_array(new_facets_array, length_that_face, facets_or_vertices[i], facet_repr);
        }
        delete[] old_facets_walker;
    }
    else {//getting the facets from the original polytope, those will be facets or vertices depending on wether we consider polar or not
        for(i = 0;i<nr_facets_given;i++){
            char_from_array(facets_pointer[i], len_facets[i], facets_or_vertices[i], facet_repr);
        }
    }
}

void CombinatorialPolyhedron::get_facets_from_incidence_matrix(unsigned int **incidence_matrix){
    length_of_face = ((nr_vertices - 1)/chunksize + 1);//this determines the length of the face in terms of chunktype
    allocate_facets();
    if (polar){
        get_facets_or_vertices_from_incidence_matrix(incidence_matrix, facets, 1, nr_facets, nr_vertices, 0);
    }
    else {
        get_facets_or_vertices_from_incidence_matrix(incidence_matrix, facets, 0, nr_vertices, nr_facets, 0);
    }
}

void CombinatorialPolyhedron::get_vertices_from_incidence_matrix(unsigned int **incidence_matrix){
    length_of_face_in_facet_repr = ((nr_facets - 1)/chunksize + 1);//this determines the length in facet representation in terms of chunktype
    allocate_vertices();
    if (!polar){
        get_facets_or_vertices_from_incidence_matrix(incidence_matrix, vertices, 1, nr_vertices, nr_facets, 1);
    }
    else {
        get_facets_or_vertices_from_incidence_matrix(incidence_matrix, vertices, 0, nr_facets, nr_vertices, 1);
    }
}

void CombinatorialPolyhedron::get_facets_or_vertices_from_incidence_matrix(unsigned int **incidence_matrix, chunktype** facets_or_vertices, unsigned int flip, unsigned int nr_vertices_given, unsigned int nr_facets_given, unsigned int facet_repr){
    unsigned int i,j;
    if (flip){//getting the vertices from the original polytope, those will be facets or vertices depending on wether we consider polar or not
        const int size_three = nr_facets_given;
        unsigned int new_facets_array[size_three];
        unsigned int length_that_face;
        for(i = 0;i<nr_vertices_given;i++){
            length_that_face = 0;
            for (j=0;j < nr_facets_given;j++){
                if (incidence_matrix[j][i]){
                    new_facets_array[length_that_face] = j;
                    length_that_face++;
                }
            }
            char_from_array(new_facets_array, length_that_face, facets_or_vertices[i], facet_repr);
        }
    }
    else {//getting the facets from the original polytope, those will be facets or vertices depending on wether we consider polar or not
        for(i = 0;i<nr_facets_given;i++){
            char_from_incidence_list(incidence_matrix[i], nr_vertices_given, facets_or_vertices[i], facet_repr);
        }
    }
}

//conversions ----------------------------------------------------

inline void CombinatorialPolyhedron::bitrep_to_list(chunktype *array1, unsigned int *face_to_return, unsigned int *length_of_faces, unsigned int facet_repr){
    unsigned int i,j;
    unsigned int face_length = length_of_face;
    if (facet_repr){
        face_length = length_of_face_in_facet_repr;
    }
    unsigned int counter = 0;
    const unsigned int size_array = face_length*chunksize/64;
    uint64_t *array = new uint64_t [size_array]();
    for (i = 0; i < face_length;i++){
        store_register(array[i*chunksize/64],array1[i]);
    }
    for (i = 0; i < size_array;i++){
        for (j = 0; j < 64; j++){
            if (array[i] >= vertex_to_bit_dictionary[j]){
                face_to_return[counter] = i*64 + j;
                counter++;
                array[i] -= vertex_to_bit_dictionary[j];
            }
        }
    }
    length_of_faces[0] = counter;
    delete[] array;
    return;
}

inline void CombinatorialPolyhedron::bitrep_to_list(chunktype **array1, unsigned int len, unsigned int **faces_to_return, unsigned int *length_of_faces, unsigned int facet_repr){
    for(unsigned int i = 0;i < len; i++){
        bitrep_to_list(array1[i], faces_to_return[i], &length_of_faces[i], facet_repr);
    }
    return;
}

void CombinatorialPolyhedron::char_from_incidence_list(unsigned int *incidence_list, unsigned int nr_vertices_given, chunktype *array1, unsigned int facet_repr){
    unsigned int face_length;
    if (facet_repr){
        face_length = length_of_face_in_facet_repr;
    }
    else {
        face_length = length_of_face;
    }
    unsigned int entry, position, value,i ;
    const unsigned int size_array = face_length*chunksize/64;
    uint64_t *array = new uint64_t [size_array]();
    while (nr_vertices_given--) {
        entry = incidence_list[nr_vertices_given];
        if (entry){
            value = nr_vertices_given % 64;
            position = nr_vertices_given/64;
            array[position] += vertex_to_bit_dictionary[value];
        }
    }
    for (i=0;i<face_length;i++){
        load_register(array1[i],array[i*chunksize/64]);
    }
    delete[] array;
}

void CombinatorialPolyhedron::char_from_array(unsigned int* input, unsigned int len, chunktype *array1, unsigned int facet_repr){
    unsigned int face_length;
    if (facet_repr){
        face_length = length_of_face_in_facet_repr;
    }
    else {
        face_length = length_of_face;
    }
    unsigned int entry, position, value,i ;
    const unsigned int size_array = face_length*chunksize/64;
    uint64_t *array = new uint64_t [size_array]();
    while (len--) {
        entry = input[len];
        value = entry % 64;
        position = entry/64;
        array[position] += vertex_to_bit_dictionary[value];
    }
    for (i=0;i<face_length;i++){
        load_register(array1[i],array[i*chunksize/64]);
    }
    delete[] array;
}




//allocation and deallocation -----------------------------------

void CombinatorialPolyhedron::allocate_facets(){
    if (facets_are_allocated){
        return;
    }
    unsigned int i;
    const int size_one = nr_facets;
    const int size_two = length_of_face;
    facets = new chunktype *[size_one];
    facets_allocator = new void *[size_one];
    for(i = 0;i<nr_facets;i++){
        facets_allocator[i] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
        facets[i] = (chunktype*) facets_allocator[i];
    }
    facets_are_allocated = 1;
}
void CombinatorialPolyhedron::deallocate_facets(){
    unsigned int i;
    if (!facets_are_allocated)
        return;
    for (i=0;i <nr_facets;i++){
        free(facets_allocator[i]);
    }
    delete[] facets;
    delete[] facets_allocator;
    facets_are_allocated = 0;
}
void CombinatorialPolyhedron::allocate_vertices(){
    if (vertices_are_allocated)
        return;
    unsigned int i;
    const int size_one = nr_vertices;
    const int size_two = length_of_face_in_facet_repr;
    vertices = new chunktype *[size_one];
    vertices_allocator = new void *[size_one];
    for(i = 0;i<nr_vertices;i++){
        vertices_allocator[i] = aligned_alloc(chunksize/8,length_of_face_in_facet_repr*chunksize/8);
        vertices[i] = (chunktype*) vertices_allocator[i];
    }
    vertices_are_allocated = 1;
}
void CombinatorialPolyhedron::deallocate_vertices(){
    unsigned int i;
    if (!vertices_are_allocated)
        return;
    for (i=0;i <nr_vertices;i++){
        free(vertices_allocator[i]);
    }
    delete[] vertices;
    delete[] vertices_allocator;
    vertices_are_allocated = 0;
}
void CombinatorialPolyhedron::allocate_newfaces(){
    if (newfaces_are_allocated)
        return;
    unsigned int i,j;
    const unsigned int const_dimension = dimension;
    const unsigned int const_facets = nr_facets;

    forbidden = new chunktype *[const_facets]();
    newfaces = new chunktype **[const_dimension -1 ]();
    newfaces2 = new chunktype **[const_dimension -1 ]();
    newfaces_allocator = new void **[const_dimension -1 ];
    for (i=0; i < dimension -1; i++){
        newfaces[i] = new chunktype*[const_facets-1]();
        newfaces2[i] = new chunktype*[const_facets-1]();
        newfaces_allocator[i] = new void*[const_facets-1];
        for (j=0; j < nr_facets-1; j++){
            newfaces_allocator[i][j] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
            newfaces[i][j] = (chunktype*) newfaces_allocator[i][j];
        }
    }
    newfaces_are_allocated = 1;
}
void CombinatorialPolyhedron::deallocate_newfaces(){
    if (!newfaces_are_allocated)
        return;
    unsigned int i,j;
    for (i=0;i <dimension -1;i++){
        for (j=0; j < nr_facets-1; j++){
            free(newfaces_allocator[i][j]);
        }
    }
    delete[] forbidden;
    delete[] newfaces2;
    delete[] newfaces;
    delete[] newfaces_allocator;
    newfaces_are_allocated = 0;
}
//allocates allfaces in a certain dimension, must be smaller than dimension and at least 1, if dimension is 0 will be allocated in all dimensions
void CombinatorialPolyhedron::allocate_allfaces(){
    unsigned int i;
    if (!dimension){
        get_dimension();
    }
    for (i = 1; i < dimension; i++){
        allocate_allfaces(i);
    }
    if (unbounded)
        allocate_allfaces(0);
}

void CombinatorialPolyhedron::allocate_allfaces(unsigned int dimension_to_allocate){
    unsigned int i;
    if (!f_vector){
        get_f_vector_and_edges();
    }
    if (!allfaces_are_allocated){
        const unsigned int const_dimension = dimension;
        allfaces_counter = new unsigned long[const_dimension]();
        allfaces_are_allocated = new unsigned int[const_dimension]();
        allfaces_allocator = new void **[const_dimension]();
        allfaces = new chunktype **[const_dimension]();
        allfaces_facet_repr_allocator = new void **[const_dimension]();
        allfaces_facet_repr = new chunktype **[const_dimension]();
    }
    if ((0 <= dimension_to_allocate) && (dimension_to_allocate < dimension))
        {
        if (!allfaces_are_allocated[dimension_to_allocate]){
            const unsigned long size_one = f_vector[dimension_to_allocate+1];
            allfaces[dimension_to_allocate] = new chunktype *[size_one]();
            allfaces_allocator[dimension_to_allocate] = new void *[size_one]();
            allfaces_facet_repr[dimension_to_allocate] = new chunktype *[size_one]();
            allfaces_facet_repr_allocator[dimension_to_allocate] = new void *[size_one]();
            for(i = 0;i<f_vector[dimension_to_allocate+1];i++){
                allfaces_allocator[dimension_to_allocate][i] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
                allfaces[dimension_to_allocate][i] = (chunktype*) allfaces_allocator[dimension_to_allocate][i];
                allfaces_facet_repr_allocator[dimension_to_allocate][i] = aligned_alloc(chunksize/8,length_of_face_in_facet_repr*chunksize/8);
                allfaces_facet_repr[dimension_to_allocate][i] = (chunktype*) allfaces_facet_repr_allocator[dimension_to_allocate][i];
            }
            allfaces_are_allocated[dimension_to_allocate] = 1;
        }
    }
}

void CombinatorialPolyhedron::deallocate_allfaces(){
    if (!allfaces_are_allocated)
        return;
    unsigned int i,j;
    for (j = 0; j < dimension; j++){
        if (allfaces_are_allocated[j]){
            for (i = 0; i < f_vector[j+1]; i++){
                free(allfaces_allocator[j][i]);
                free(allfaces_facet_repr_allocator[j][i]);
            }
        }
    }
    delete[] allfaces_counter;
    delete[] allfaces;
    delete[] allfaces_allocator;
    delete[] allfaces_facet_repr;
    delete[] allfaces_facet_repr_allocator;
    delete[] allfaces_are_allocated;
}




CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(unsigned int ** facets_pointer, unsigned int nr_facets, unsigned int *len_facets, unsigned int nr_vertices, int is_unbounded) {
    CombinatorialPolyhedron_ptr C = new CombinatorialPolyhedron(facets_pointer, nr_facets, len_facets, nr_vertices, is_unbounded);
    return C;
}

CombinatorialPolyhedron_ptr init_CombinatorialPolyhedron(unsigned int ** incidence_matrix, unsigned int nr_facets, unsigned int nr_vertices, int is_unbounded) {
    CombinatorialPolyhedron_ptr C = new CombinatorialPolyhedron(incidence_matrix, nr_facets, nr_vertices, is_unbounded);
    return C;
}

unsigned int dimension(CombinatorialPolyhedron_ptr C){
  return (*C).get_dimension();
}

unsigned int ** edges(CombinatorialPolyhedron_ptr C){
  return (*C).get_edges();
}


void f_vector(CombinatorialPolyhedron_ptr C, unsigned long *vector){
  return (*C).get_f_vector(vector);
}

unsigned int ** ridges(CombinatorialPolyhedron_ptr C){
  return (*C).get_ridges();
}

unsigned long ** incidences(CombinatorialPolyhedron_ptr C, int dimension_one, int dimension_two, unsigned long * nr_incidences, unsigned int * twisted){
    return (*C).get_incidences(dimension_one, dimension_two, nr_incidences, twisted);
}

void record_all_faces(CombinatorialPolyhedron_ptr C){
    (*C).record_all_faces();
}

void get_faces(CombinatorialPolyhedron_ptr C, int dimension, unsigned int facet_repr, unsigned int **faces_to_return, unsigned int *length_of_faces){
    return (*C).get_faces(dimension, facet_repr, faces_to_return, length_of_faces);
}

unsigned int face_iterator(CombinatorialPolyhedron_ptr C, unsigned int *Vface_to_return, unsigned int *Vlength, unsigned int *Hface_to_return, unsigned int *Hlength){
    return (*C).face_iterator_call(Vface_to_return, Vlength, Hface_to_return, Hlength);
}

void face_iterator_init(CombinatorialPolyhedron_ptr C, int dimension, unsigned int vertex_repr, unsigned int facet_repr){
    return (*C).face_iterator_init(dimension, vertex_repr, facet_repr);
}

unsigned long get_flag(CombinatorialPolyhedron_ptr C, unsigned int * flagarray, unsigned int length){
    return (*C).get_flag_number_init(flagarray, length);
}


void delete_CombinatorialPolyhedron(CombinatorialPolyhedron_ptr C){
  delete(C);
}

unsigned long get_maxnumberedges(){
    return (unsigned long) maxnumberedges;
}

unsigned long get_maxnumberincidences(){
    return (unsigned long) maxnumberincidences;
}
