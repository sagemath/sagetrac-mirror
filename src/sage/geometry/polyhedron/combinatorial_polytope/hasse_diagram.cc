//Copyright: see base.pyx

#include <Python.h>
#include <math.h>
#include "hasse_diagram.h"



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
CombinatorialPolytope::CombinatorialPolytope(PyObject* py_tuple, unsigned int nr_vertices_given){
    build_dictionary();
    nr_vertices = nr_vertices_given;
    nr_facets = PyTuple_Size(py_tuple);
    if (nr_facets > nr_vertices){//in this case the polar approach is actually better, so we will save the polar CombinatorialPolytope and compute accordingly
        //Polar_Init(py_tuple, nr_vertices_given);
        polar = 1;
        nr_vertices = nr_facets;
        nr_facets = nr_vertices_given;
    }
    get_facets_from_tuple(py_tuple);
    get_vertices_from_tuple(py_tuple);
}

//initialization with an incidence matrix given as tuple of tuples
CombinatorialPolytope::CombinatorialPolytope(PyObject* py_tuple){
    build_dictionary();
    nr_vertices = PyTuple_Size(PyTuple_GetItem(py_tuple,0));
    nr_facets = PyTuple_Size(py_tuple);
    if (nr_facets > nr_vertices){//in this case the polar approach is actually much better, so we will save the polar CombinatorialPolytope and compute accordingly
        polar = 1;
        unsigned int nr_vertices_given = nr_vertices;
        nr_vertices = nr_facets;
        nr_facets = nr_vertices_given;
    }
    get_facets_from_incidence_matrix(py_tuple);
    get_vertices_from_incidence_matrix(py_tuple);
}

CombinatorialPolytope::~CombinatorialPolytope(){
    unsigned int i,j;
    deallocate_facets();
    deallocate_newfaces();
    deallocate_vertices();
    deallocate_allfaces();//must be called before deleting f_vector
    delete(f_vector);
    delete(edges);
    delete(ridges);
}

unsigned int CombinatorialPolytope::get_dimension(){
    if (dimension){
        return dimension;
    }
    dimension =  calculate_dimension(facets, nr_facets);
    return dimension;
}

inline PyObject* CombinatorialPolytope::get_f_vector(){
    if (!f_vector){
        calculate_f_vector();
    }
    return tuple_from_f_vector();
}

inline PyObject* CombinatorialPolytope::get_edges(){
    if (polar){
        if (!nr_ridges){
            calculate_ridges();
        }
        return tuple_from_ridges();
    }
    else {
        if (!nr_edges){
            edgemode = 1;
            calculate_f_vector();
        }
        return tuple_from_edges();
    }
}

inline PyObject* CombinatorialPolytope::get_ridges(){
    if (!polar){
        if (!nr_ridges){
            calculate_ridges();
        }
        return tuple_from_ridges();
    }
    else {
        if (!nr_edges){
            edgemode = 1;
            calculate_f_vector();
        }
        return tuple_from_edges();
    }
}

inline PyObject* CombinatorialPolytope::get_faces(unsigned int face_dimension){
    //TODO, special cases, dimension to low, face_dimension == 1, etc.
    if (face_dimension == dimension -1){
        return tuple_from_faces(facets, nr_facets,0);
    }
    if (!allfaces_are_allocated || !allfaces_are_allocated[face_dimension]){
        allocate_allfaces(face_dimension);
        record_faces(facets,dimension-1,nr_facets,0,face_dimension);
    }
    return tuple_from_faces(allfaces[face_dimension],allfaces_counter[face_dimension],0);
}

void CombinatorialPolytope::record_all_faces(){
    allocate_allfaces(0);
    record_faces(facets,dimension-1,nr_facets,0,1);
}

//will set C to be the intersection of A and B
inline void CombinatorialPolytope::intersection(chunktype *A, chunktype *B, chunktype *C){
    unsigned int i;
    for (i = 0; i < length_of_face; i++){
        C[i] = bitwise_intersection(A[i],B[i]);
    }
}

//returns 1 if A is a proper subset of B, otherwise returns 0, this is done by checking if there is an element in A, which is not in B
inline int CombinatorialPolytope::is_subset(chunktype *A, chunktype *B){
    unsigned int i;
    for (i = 0; i < length_of_face; i++){
        if (bitwise_is_not_subset(A[i],B[i])){
            return 0;
        }
    }
    return 1;
}

//counts the number of vertices in a face by counting bits set to one
inline unsigned int CombinatorialPolytope::CountFaceBits(chunktype* A1) {
    //this function is not implemented for speed (it basically gets called dimension times and once to convert a face to a tuple
    unsigned int i,count = 0;
    const unsigned int length_of_conversion_face = length_of_face*chunksize/64;
    unsigned long A[length_of_conversion_face];
    unsigned long n;
    for (i=0;i<length_of_face;i++){
        store_register(A[i*chunksize/64],A1[i]);
    }
    for (i=0;i<length_of_conversion_face;i++){
        count += popcount(A[i]);
    }
    return count;
}

//adds an edge to the edges list
inline void CombinatorialPolytope::add_edge(chunktype *face){
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
inline void CombinatorialPolytope::add_edge(unsigned int one, unsigned int two){
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
inline void CombinatorialPolytope::add_ridge(unsigned int one, unsigned int two){
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

//intersects the first 'lenfaces' faces of 'faces' with the 'face_to_intersect'-th face of faces and stores the result in 'nextfaces'
//then determines which ones are exactly of one dimension less by considering containment
//newfaces2 will point at those of exactly one dimension less which are not contained in any of the faces in 'forbidden'
//returns the number of those faces
inline unsigned int CombinatorialPolytope::get_next_level(chunktype **faces, unsigned int lenfaces, unsigned int face_to_intersect, chunktype **nextfaces, chunktype **nextfaces2, unsigned int nr_forbidden){
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
    int addfacearray[lenfaces] = { };
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

unsigned int CombinatorialPolytope::calculate_dimension(chunktype **faces, unsigned int nr_faces){
    unsigned int i,j,k, newfacescounter, dimension;
    unsigned int bitcount = CountFaceBits(faces[0]);
    if (bitcount == 1){//we have encountered vertex, if the facets have dimension 0, the polytope has dimension 1
        return 1;
    }
    if (bitcount == 2){//we have encountered edge
        return 2;
    }
    if (bitcount == 0){//the polytope has only the empty face as facet
        return 0;
    }
    void *nextfaces_creator[nr_faces-1];
    chunktype  *nextfaces2[nr_faces-1], *nextfaces[nr_faces-1];
    for (i=0; i < (nr_faces-1); i++){
        nextfaces_creator[i] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
        nextfaces[i] = (chunktype*) nextfaces_creator[i];
    }
    newfacescounter = get_next_level(faces,nr_faces,nr_faces-1,nextfaces,nextfaces2,0);//calculates the ridges contained in one facet
    dimension =  calculate_dimension(nextfaces2,newfacescounter) + 1;//calculates the dimension of that facet
    for (i=0; i < (nr_faces - 1); i++){
        free(nextfaces_creator[i]);
    }
    return dimension;
}

void CombinatorialPolytope::calculate_ridges(){//this is a much simpler version of belows get_f_vector_and_edges
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

void CombinatorialPolytope::get_f_vector_and_edges(chunktype **faces, unsigned int dimension, unsigned int nr_faces, unsigned int nr_forbidden){
    unsigned int i,j;
    unsigned long newfacescounter;
    if (dimension == 1){
        if (edgemode)//in this case we want to record the edges
        for (i = 0; i < nr_faces; i++){
            add_edge(faces[i]);
        }
    }
    if (dimension == 0){
        return;
    }
    i = nr_faces;
    while (i--){
        newfacescounter = get_next_level(faces,i+1,i,newfaces[dimension-1],newfaces2[dimension-1],nr_forbidden);//get the facets contained in faces[i] but not in any of the forbidden
        f_vector[dimension] += newfacescounter;
        if (newfacescounter){
            get_f_vector_and_edges(newfaces2[dimension-1],dimension-1,newfacescounter,nr_forbidden);//add all face in faces[i] to the f_vector, but not those which we have counted already
        }
        forbidden[nr_forbidden] = faces[i];//we have counted all faces in faces[i], so we do not want to count them ever again
        nr_forbidden++;
    }
}


void CombinatorialPolytope::record_faces(chunktype **faces, unsigned int current_dimension, unsigned int nr_faces, unsigned int nr_forbidden, unsigned int lowest_dimension){
    unsigned int i,j;
    unsigned long newfacescounter;
    if ((current_dimension < dimension - 1) && allfaces_are_allocated[current_dimension]){
        for (i = 0; i < nr_faces; i++){
            record_face(faces[i], current_dimension);
        }
    }
    if (current_dimension == lowest_dimension){
        return;
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

inline void CombinatorialPolytope::record_face(chunktype *face, unsigned int current_dimension){
    for (unsigned int j = 0; j < length_of_face; j++){
        load_register(allfaces[current_dimension][allfaces_counter[current_dimension]][j],face[j]);
    }
    allfaces_counter[current_dimension]++;
}

 inline void CombinatorialPolytope::calculate_f_vector(){
    unsigned int i,j;
    if (!dimension){
        get_dimension();
    }
    allocate_newfaces();
    const unsigned int const_dimension = dimension;
    f_vector = new unsigned long[const_dimension + 2]();
    f_vector[0] = 1;
    f_vector[dimension+1] = 1;
    f_vector[dimension] = nr_facets;
    //f_vector[1] = nr_vertices; //this is commented out in order to make calculations also work for unbounded polyhedra
    get_f_vector_and_edges(facets,dimension-1,nr_facets,0);
}


void CombinatorialPolytope::get_facets_from_tuple(PyObject* py_tuple){
    length_of_face = ((nr_vertices - 1)/chunksize + 1);//this determines the length of the face in terms of chunktype
    allocate_facets();
    if (polar){
        get_facets_or_vertices_from_tuple(py_tuple, facets, 1, nr_facets, nr_vertices, 0);
    }
    else {
        get_facets_or_vertices_from_tuple(py_tuple, facets, 0, nr_vertices, nr_facets, 0);
    }
}

void CombinatorialPolytope::get_vertices_from_tuple(PyObject* py_tuple){
    length_of_face_in_facet_repr = ((nr_facets - 1)/chunksize + 1);//this determines the length in facet representation in terms of chunktype
    allocate_vertices();
    if (!polar){
        get_facets_or_vertices_from_tuple(py_tuple, vertices, 1, nr_vertices, nr_facets, 1);
    }
    else {
        get_facets_or_vertices_from_tuple(py_tuple, vertices, 0, nr_facets, nr_vertices, 1);
    }
}

void CombinatorialPolytope::get_facets_or_vertices_from_tuple(PyObject* py_tuple, chunktype** facets_or_vertices, unsigned int flip, unsigned int nr_vertices_given, unsigned int nr_facets_given, unsigned int facet_repr){
    unsigned int i,j;
    if (flip){//getting the vertices from the original polytope, those will be facets or vertices depending on wether we consider polar or not
        const int size_three = nr_facets_given;
        unsigned int *old_facets_walker = new unsigned int [size_three]();
        unsigned int new_facets_array[size_three];
        unsigned int length_that_face;
        for(i = 0;i<nr_vertices_given;i++){
            length_that_face = 0;
            for (j=0;j < nr_facets_given;j++){
                if (i == (unsigned int) PyInt_AsLong(PyTuple_GetItem(PyTuple_GetItem(py_tuple, j),old_facets_walker[j]))){
                    new_facets_array[length_that_face] = j;
                    old_facets_walker[j]++;
                    length_that_face++;
                    if (old_facets_walker[j] >= PyTuple_Size(PyTuple_GetItem(py_tuple, j)))
                        old_facets_walker[j]--;
                }
            }
            char_from_array(new_facets_array, length_that_face, facets_or_vertices[i], facet_repr);
        }
        delete(old_facets_walker);
    }
    else {//getting the facets from the original polytope, those will be facets or vertices depending on wether we consider polar or not
        for(i = 0;i<nr_facets_given;i++){
            char_from_tuple(PyTuple_GetItem(py_tuple,i), facets_or_vertices[i], facet_repr);
        }
    }
}

void CombinatorialPolytope::get_facets_from_incidence_matrix(PyObject* py_tuple){
    length_of_face = ((nr_vertices - 1)/chunksize + 1);//this determines the length of the face in terms of chunktype
    allocate_facets();
    if (polar){
        get_facets_or_vertices_from_incidence_matrix(py_tuple, facets, 1, nr_facets, nr_vertices, 0);
    }
    else {
        get_facets_or_vertices_from_incidence_matrix(py_tuple, facets, 0, nr_vertices, nr_facets, 0);
    }
}

void CombinatorialPolytope::get_vertices_from_incidence_matrix(PyObject* py_tuple){
    length_of_face_in_facet_repr = ((nr_facets - 1)/chunksize + 1);//this determines the length in facet representation in terms of chunktype
    allocate_vertices();
    if (!polar){
        get_facets_or_vertices_from_incidence_matrix(py_tuple, vertices, 1, nr_vertices, nr_facets, 1);
    }
    else {
        get_facets_or_vertices_from_incidence_matrix(py_tuple, vertices, 0, nr_facets, nr_vertices, 1);
    }
}

void CombinatorialPolytope::get_facets_or_vertices_from_incidence_matrix(PyObject* py_tuple, chunktype** facets_or_vertices, unsigned int flip, unsigned int nr_vertices_given, unsigned int nr_facets_given, unsigned int facet_repr){
    unsigned int i,j;
    if (flip){//getting the vertices from the original polytope, those will be facets or vertices depending on wether we consider polar or not
        const int size_three = nr_facets_given;
        unsigned int new_facets_array[size_three];
        unsigned int length_that_face;
        for(i = 0;i<nr_vertices_given;i++){
            length_that_face = 0;
            for (j=0;j < nr_facets_given;j++){
                if (PyInt_AsLong(PyTuple_GetItem(PyTuple_GetItem(py_tuple, j),i))){
                    new_facets_array[length_that_face] = j;
                    length_that_face++;
                }
            }
            char_from_array(new_facets_array, length_that_face, facets_or_vertices[i], facet_repr);
        }
    }
    else {//getting the facets from the original polytope, those will be facets or vertices depending on wether we consider polar or not
        for(i = 0;i<nr_facets_given;i++){
            char_from_incidence_tuple(PyTuple_GetItem(py_tuple,i), facets_or_vertices[i], facet_repr);
        }
    }
}

//conversions ----------------------------------------------------

void CombinatorialPolytope::char_from_tuple(PyObject* py_tuple, chunktype *array1, unsigned int facet_repr){
    unsigned int face_length;
    if (facet_repr){
        face_length = length_of_face_in_facet_repr;
    }
    else {
        face_length = length_of_face;
    }
    unsigned int len, entry, position, value,i ;
    const unsigned int size_array = face_length*chunksize/64;
    uint64_t *array = new uint64_t [size_array]();
    len = PyTuple_Size(py_tuple);
    while (len--) {
        entry = (unsigned int) PyInt_AsLong(PyTuple_GetItem(py_tuple, len));
        value = entry % 64;
        position = entry/64;
        array[position] += vertex_to_bit_dictionary[value];
    }
    for (i=0;i<face_length;i++){
        load_register(array1[i],array[i*chunksize/64]);
    }
    delete(array);
}

inline PyObject* CombinatorialPolytope::tuple_from_char(chunktype *array1, unsigned int facet_repr){
    unsigned int i,j;
    unsigned int face_length = length_of_face;
    if (facet_repr){
        face_length = length_of_face_in_facet_repr;
    }
    unsigned int len = CountFaceBits(array1);
    unsigned int counter = 0;
    PyObject* py_tuple = PyTuple_New(len);
    const unsigned int size_array = face_length*chunksize/64;
    uint64_t *array = new uint64_t [size_array]();
    for (i = 0; i < face_length;i++){
        store_register(array[i*chunksize/64],array1[i]);
    }
    for (i = 0; i < size_array;i++){
        for (j = 0; j < 64; j++){
            if (array[i] >= vertex_to_bit_dictionary[j]){
                PyTuple_SET_ITEM(py_tuple,counter,PyLong_FromUnsignedLong(i*64+j));
                counter++;
                array[i] -= vertex_to_bit_dictionary[j];
            }
        }
    }
    return py_tuple;
}

inline PyObject* CombinatorialPolytope::tuple_from_faces(chunktype **array1, unsigned int len, unsigned int facet_repr){
    PyObject* py_tuple = PyTuple_New(len);
    for(unsigned int i = 0;i < len; i++){
        PyTuple_SET_ITEM(py_tuple,i,tuple_from_char(array1[i],0));
    }
    return py_tuple;
}

void CombinatorialPolytope::char_from_incidence_tuple(PyObject* py_tuple, chunktype *array1, unsigned int facet_repr){
    unsigned int face_length;
    if (facet_repr){
        face_length = length_of_face_in_facet_repr;
    }
    else {
        face_length = length_of_face;
    }
    unsigned int len, entry, position, value,i ;
    const unsigned int size_array = face_length*chunksize/64;
    uint64_t *array = new uint64_t [size_array]();
    len = PyTuple_Size(py_tuple);
    while (len--) {
        entry = (unsigned int) PyInt_AsLong(PyTuple_GetItem(py_tuple, len));
        if (entry){
            value = len % 64;
            position = len/64;
            array[position] += vertex_to_bit_dictionary[value];
        }
    }
    for (i=0;i<face_length;i++){
        load_register(array1[i],array[i*chunksize/64]);
    }
    delete(array);
}

void CombinatorialPolytope::char_from_array(unsigned int* input, unsigned int len, chunktype *array1, unsigned int facet_repr){
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
    delete(array);
}

inline PyObject* CombinatorialPolytope::tuple_from_f_vector(){
    PyObject *py_tuple;
    unsigned int i;
    py_tuple = PyTuple_New(dimension + 2);
    if (polar){
        for(i=0;i< dimension + 2;i++){
            PyTuple_SET_ITEM(py_tuple,dimension + 1 - i,PyLong_FromUnsignedLong(f_vector[i]));
        }
        return py_tuple;
    }
    for(i=0;i< dimension + 2;i++){
        PyTuple_SET_ITEM(py_tuple,i,PyLong_FromUnsignedLong(f_vector[i]));
    }
    return py_tuple;
}

inline PyObject* CombinatorialPolytope::tuple_from_edges(){
    PyObject *py_tuple;
    py_tuple = PyTuple_New(nr_edges);
    unsigned int i;
    unsigned int position_one, position_two;
    for (i=0;i < nr_edges;i++){
        position_one = i / maxnumberedges;
        position_two = 2*(i % maxnumberedges);
        PyObject *py_edge = PyTuple_New(2);
        PyTuple_SET_ITEM(py_edge,0,Py_BuildValue("i",edges[position_one][position_two]));
        PyTuple_SET_ITEM(py_edge,1,Py_BuildValue("i",edges[position_one][position_two+1]));
        PyTuple_SET_ITEM(py_tuple,i,py_edge);
    }
    return py_tuple;
}

inline PyObject* CombinatorialPolytope::tuple_from_ridges(){
    PyObject *py_tuple;
    py_tuple = PyTuple_New(nr_ridges);
    unsigned int i;
    unsigned int position_one, position_two;
    for (i=0;i < nr_ridges;i++){
        position_one = i / maxnumberedges;
        position_two = 2*(i % maxnumberedges);
        PyObject *py_ridge = PyTuple_New(2);
        PyTuple_SET_ITEM(py_ridge,0,Py_BuildValue("i",ridges[position_one][position_two]));
        PyTuple_SET_ITEM(py_ridge,1,Py_BuildValue("i",ridges[position_one][position_two+1]));
        PyTuple_SET_ITEM(py_tuple,i,py_ridge);
    }
    return py_tuple;
}


//allocation and deallocation -----------------------------------

void CombinatorialPolytope::allocate_facets(){
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
void CombinatorialPolytope::deallocate_facets(){
    unsigned int i;
    if (!facets_are_allocated)
        return;
    for (i=0;i <nr_facets;i++){
        free(facets_allocator[i]);
    }
    delete(facets);
    delete(facets_allocator);
    facets_are_allocated = 0;
}
void CombinatorialPolytope::allocate_vertices(){
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
void CombinatorialPolytope::deallocate_vertices(){
    unsigned int i;
    if (!vertices_are_allocated)
        return;
    for (i=0;i <nr_vertices;i++){
        free(vertices_allocator[i]);
    }
    delete(vertices);
    delete(vertices_allocator);
    vertices_are_allocated = 0;
}
void CombinatorialPolytope::allocate_newfaces(){
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
void CombinatorialPolytope::deallocate_newfaces(){
    if (!newfaces_are_allocated)
        return;
    unsigned int i,j;
    for (i=0;i <dimension -1;i++){
        for (j=0; j < nr_facets-1; j++){
            free(newfaces_allocator[i][j]);
        }
    }
    delete(forbidden);
    delete(newfaces2);
    delete(newfaces);
    delete(newfaces_allocator);
    newfaces_are_allocated = 0;
}
//allocates allfaces in a certain dimension, must be smaller than dimension and at least 1, if dimension is 0 will be allocated in all dimensions
void CombinatorialPolytope::allocate_allfaces(unsigned int dimension_to_allocate){
    unsigned int i;
    if (!f_vector){
        calculate_f_vector();
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
    if (dimension_to_allocate == 0){
        for (i = 1; i < dimension; i++){
            allocate_allfaces(i);
        }
    }
    if ((1 <= dimension_to_allocate) && (dimension_to_allocate < dimension))
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

void CombinatorialPolytope::deallocate_allfaces(){
    if (!allfaces_are_allocated)
        return;
    unsigned int i,j;
    for (j = 1; j < dimension; j++){
        if (allfaces_are_allocated[j]){
            for (i = 0; i < f_vector[j+1]; i++){
                free(allfaces_allocator[j][i]);
                free(allfaces_facet_repr_allocator[j][i]);
            }
        }
    }
    delete(allfaces_counter);
    delete(allfaces);
    delete(allfaces_allocator);
    delete(allfaces_facet_repr);
    delete(allfaces_facet_repr_allocator);
    delete(allfaces_are_allocated);
}




CombinatorialPolytope_ptr init_CombinatorialPolytope(PyObject* py_tuple, unsigned int nr_vertices) {
    CombinatorialPolytope_ptr C = new CombinatorialPolytope(py_tuple,nr_vertices);
    return C;
}

CombinatorialPolytope_ptr init_CombinatorialPolytope(PyObject* py_tuple) {
    CombinatorialPolytope_ptr C = new CombinatorialPolytope(py_tuple);
    return C;
}

unsigned int dimension(CombinatorialPolytope_ptr C){
  return (*C).get_dimension();
}

PyObject* edges(CombinatorialPolytope_ptr C){
  return (*C).get_edges();
}

PyObject* f_vector(CombinatorialPolytope_ptr C){
  return (*C).get_f_vector();
}

PyObject* ridges(CombinatorialPolytope_ptr C){
  return (*C).get_ridges();
}

void record_all_faces(CombinatorialPolytope_ptr C){
    (*C).record_all_faces();
}

PyObject* get_faces(CombinatorialPolytope_ptr C, unsigned int dimension){
    return (*C).get_faces(dimension);
}

void delete_CombinatorialPolytope(CombinatorialPolytope_ptr C){
  delete(C);
}
