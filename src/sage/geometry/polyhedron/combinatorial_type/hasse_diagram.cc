//Copyright: see base.pyx

#include <Python.h>
#include <math.h>
#include "hasse_diagram.h"



static uint64_t vertex_to_bit_dictionary[64];

void build_dictionary(){
    unsigned int i = 0;
    uint64_t count = 1;
    for (i=0; i< 64;i++){
        vertex_to_bit_dictionary[64-i-1] = count;
        count *= 2;
    }
}

//taken from https://codingforspeed.com/counting-the-number-of-leading-zeros-for-a-32-bit-integer-signed-or-unsigned/
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






CombinatorialType::CombinatorialType(PyObject* py_tuple, unsigned int nr_vertices_given){//this is initalizing CombinatorialType with a tuple of facets, each facet given as a tuple of vertices, the vertices must be exactly 0,...,nr_vertices-1
            build_dictionary();
            unsigned int i;
            nr_vertices = nr_vertices_given;
            nr_facets = PyTuple_Size(py_tuple);
            if (nr_facets > nr_vertices){//in this case the polar approach is actually much better, so we will save the polar CombinatorialType and compute accordingly
                Polar_Init(py_tuple, nr_vertices_given);
                return;
            }
            length_of_face = ((nr_vertices - 1)/chunksize + 1);//this determines the length of the face in terms of chunktype
            const int size_one = nr_facets;
            const int size_two = length_of_face;
            facets = new chunktype *[size_one];
            facets_allocator = new void *[size_one];
            for(i = 0;i<nr_facets;i++){
                facets_allocator[i] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
                facets[i] = (chunktype*) facets_allocator[i];
                //chunktype saver[size_two];
                //facets[i] = saver;
                //facets[i] = new chunktype [size_two];
                char_from_tuple(PyTuple_GetItem(py_tuple,i), facets[i]);
            }
            facets_are_allocated = 1;
        }
CombinatorialType::CombinatorialType(PyObject* py_tuple){//this is initalizing CombinatorialType incidence matrix as tuple of tuples
            build_dictionary();
            unsigned int i;
            nr_vertices = PyTuple_Size(PyTuple_GetItem(py_tuple,0));
            nr_facets = PyTuple_Size(py_tuple);
            if (nr_facets > nr_vertices){//in this case the polar approach is actually much better, so we will save the polar CombinatorialType and compute accordingly
                Polar_Init(py_tuple);
                return;
            }
            length_of_face = ((nr_vertices - 1)/chunksize + 1);//this determines the length of the face in terms of chunktype
            const int size_one = nr_facets;
            const int size_two = length_of_face;
            facets = new chunktype *[size_one];
            facets_allocator = new void *[size_one];
            for(i = 0;i<nr_facets;i++){
                facets_allocator[i] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
                facets[i] = (chunktype*) facets_allocator[i];
                //chunktype saver[size_two];
                //facets[i] = saver;
                //facets[i] = new chunktype [size_two];
                char_from_incidence_tuple(PyTuple_GetItem(py_tuple,i), facets[i]);
            }
            facets_are_allocated = 1;
        }
CombinatorialType::~CombinatorialType(){
            unsigned int i,j;
            if (facets_are_allocated){
                for (i=0;i <nr_facets;i++){
                    free(facets_allocator[i]);
                }
            }
            if (newfaces_are_allocated){
                for (i=0;i <dimension -1;i++){
                    for (j=0; j < nr_facets-1; j++){
                        free(newfaces_allocator[i][j]);
                    }
                }
            }
            delete(edges);
            delete(ridges);
        }
unsigned int CombinatorialType::get_dimension(){
            if (dimension){
                return dimension;
            }
            dimension =  calculate_dimension(facets, nr_facets);
            return dimension;
        }
inline PyObject* CombinatorialType::get_f_vector(){
            //need to take care of the polar case yet
            if (!f_vector){
                calculate_f_vector();
            }
            return tuple_from_f_vector();
        }
inline PyObject* CombinatorialType::get_edges(){
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
inline PyObject* CombinatorialType::get_ridges(){
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
inline void CombinatorialType::intersection(chunktype *A, chunktype *B, chunktype *C){
            unsigned int i;
            for (i = 0; i < length_of_face; i++){
                C[i] = bitwise_intersection(A[i],B[i]);
            }
        }
inline int CombinatorialType::is_subset(chunktype *A, chunktype *B){
            unsigned int i;
            for (i = 0; i < length_of_face; i++){
                if (bitwise_is_not_subset(A[i],B[i])){
                    return 0;
                }
            }
            return 1;
        }
inline unsigned int CombinatorialType::CountFaceBits(chunktype* A1) {
            //this function is not implemented for speed (it basically gets called dimension times and once to convert a face to a tuple
            unsigned int i,count = 0;
            const unsigned int length_of_conversion_face = length_of_face*chunksize/64;
            unsigned long A[length_of_conversion_face];
            unsigned long n;
            for (i=0;i<length_of_face;i++){
                store_register(A[i*chunksize/64],A1[i]);
            }
            for (i=0;i<length_of_conversion_face;i++){
                n = A[i];
                while (n) { 
                    count += n & 1; 
                    n >>= 1; 
                }
            }
            return count;
        }
inline void CombinatorialType::add_edge(chunktype *face){
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
            //first_two_vertices_of_face(face,&one,&two);
            //add_edge_helper(one,two);
        }
inline void CombinatorialType::add_edge(unsigned int one, unsigned int two){
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
inline void CombinatorialType::add_ridge(unsigned int one, unsigned int two){
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
inline unsigned int CombinatorialType::get_next_level(chunktype **faces, unsigned int lenfaces, unsigned int face_to_intersect, chunktype **nextfaces, chunktype **nextfaces2, unsigned int nr_forbidden){
            unsigned int j, k, addthisface;
            unsigned int newfacescounter = 0;
            unsigned int counter = 0;
            for (j = 0; j < lenfaces; j++){
                if (j != face_to_intersect){
                    intersection(faces[j],faces[face_to_intersect],nextfaces[counter]);
                    counter ++;
                }
            }
            //we have create all possible intersection with the i_th-face, but some might not be proper i_th-faces
            for (j = 0; j< lenfaces-1; j++){
                addthisface = 1;
                for(k = 0; k < newfacescounter; k++){
                    if(is_subset(nextfaces[j],nextfaces2[k])){
                        addthisface = 0;
                        break;
                    }
                }
                if (!addthisface) {
                    continue;
                }
                for(k = j+1; k < lenfaces-1; k++){
                    if(is_subset(nextfaces[j],nextfaces[k])){
                        addthisface = 0;
                        break;
                    }
                }
                if (!addthisface) {
                    continue;
                }
                //and some of them we might have counted already
                for (k = 0; k < nr_forbidden; k++){
                    if(is_subset(nextfaces[j],forbidden[k])){
                        addthisface = 0;
                        break;
                    }
                }
                if (!addthisface) {
                    continue;
                }
                nextfaces2[newfacescounter] = nextfaces[j];
                newfacescounter++;
            }
            return newfacescounter;
        }
unsigned int CombinatorialType::calculate_dimension(chunktype **faces, unsigned int nr_faces){
            unsigned int i,j,k, newfacescounter, dimension;
            unsigned int bitcount = CountFaceBits(faces[0]);
            if (bitcount == 1){
                return 1;
            }
            if (bitcount == 2){
                return 2;
            }
            void *nextfaces_creator[nr_faces-1];
            chunktype  *nextfaces2[nr_faces-1], *nextfaces[nr_faces-1];
            for (i=0; i < (nr_faces-1); i++){
                nextfaces_creator[i] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
                nextfaces[i] = (chunktype*) nextfaces_creator[i];
            }
            newfacescounter = get_next_level(faces,nr_faces,nr_faces-1,nextfaces,nextfaces2,0);
            dimension =  calculate_dimension(nextfaces2,newfacescounter) + 1;
            for (i=0; i < (nr_faces - 1); i++){
                free(nextfaces_creator[i]);
            }
            return dimension;
        }
void CombinatorialType::calculate_ridges(){//this is a much simpler version of belows get_f_vector_and_edges
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
void CombinatorialType::get_f_vector_and_edges(chunktype **faces, unsigned int dimension, unsigned int nr_faces, unsigned int nr_forbidden){
            unsigned int i,j,counter, addthisface;
            unsigned long newfacescounter;
            if (dimension == 1){
                if (edgemode && !polar)
                for (i = 0; i < nr_faces; i++){
                    add_edge(faces[i]);
                }
            }
            if (dimension == 0){
                return;
            }
            i = nr_faces;
            while (i--){
                newfacescounter = get_next_level(faces,i+1,i,newfaces[dimension-1],newfaces2[dimension-1],nr_forbidden);
                f_vector[dimension] += newfacescounter;
                if (newfacescounter){
                    get_f_vector_and_edges(newfaces2[dimension-1],dimension-1,newfacescounter,nr_forbidden);
                }
                forbidden[nr_forbidden] = faces[i];
                nr_forbidden++;
            }
        }
inline void CombinatorialType::calculate_f_vector(){
            unsigned int i,j;
            if (!dimension){
                get_dimension();
            }
            const unsigned int const_dimension = dimension;
            const unsigned int const_facets = nr_facets;
            f_vector = new unsigned long[const_dimension + 2]();
            forbidden = new chunktype *[const_facets];
            newfaces = new chunktype **[const_dimension -1 ];
            newfaces2 = new chunktype **[const_dimension -1 ];
            newfaces_allocator = new void **[const_dimension -1 ];
            for (i=0; i < dimension -1; i++){
                newfaces[i] = new chunktype*[const_facets-1];
                newfaces2[i] = new chunktype*[const_facets-1];
                newfaces_allocator[i] = new void*[const_facets-1];
                for (j=0; j < nr_facets-1; j++){
                    newfaces_allocator[i][j] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
                    newfaces[i][j] = (chunktype*) newfaces_allocator[i][j];
                }
            }
            newfaces_are_allocated = 1;
            f_vector[0] = 1;
            f_vector[dimension+1] = 1;
            f_vector[dimension] = nr_facets;
            //f_vector[1] = nr_vertices; //this is commented out in order to make calculations also work for unbounded polyhedra
            get_f_vector_and_edges(facets,dimension-1,nr_facets,0);
        }
void CombinatorialType::Polar_Init(PyObject* py_tuple, unsigned int nr_vertices_given){
            unsigned int i,j;
            nr_facets = nr_vertices_given;
            nr_vertices = PyTuple_Size(py_tuple);
            polar = 1;
            length_of_face = ((nr_vertices - 1)/chunksize + 1);//this determines the length of the face in terms of chunktype
            const int size_one = nr_facets;
            const int size_two = length_of_face;
            const int size_three = nr_vertices;
            facets = new chunktype *[size_one];
            facets_allocator = new void *[size_one];
            unsigned int *old_facets_walker = new unsigned int [size_three]();
            unsigned int new_facets_array[size_three];
            unsigned int length_that_face;
            for(i = 0;i<nr_facets;i++){
                facets_allocator[i] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
                facets[i] = (chunktype*) facets_allocator[i];
                length_that_face = 0;
                for (j=0;j < nr_vertices;j++){
                    if (i == (unsigned int) PyInt_AsLong(PyTuple_GetItem(PyTuple_GetItem(py_tuple, j),old_facets_walker[j]))){
                        new_facets_array[length_that_face] = j;
                        old_facets_walker[j]++;
                        length_that_face++;
                        if (old_facets_walker[j] >= PyTuple_Size(PyTuple_GetItem(py_tuple, j)))
                            old_facets_walker[j]--;
                    }
                }
                char_from_array(new_facets_array, length_that_face, facets[i]);
            }
            facets_are_allocated = 1;
        }
void CombinatorialType::Polar_Init(PyObject* py_tuple){
            unsigned int i,j;
            nr_facets = PyTuple_Size(PyTuple_GetItem(py_tuple,0));
            nr_vertices = PyTuple_Size(py_tuple);
            polar = 1;
            length_of_face = ((nr_vertices - 1)/chunksize + 1);//this determines the length of the face in terms of chunktype
            const int size_one = nr_facets;
            const int size_two = length_of_face;
            const int size_three = nr_vertices;
            facets = new chunktype *[size_one];
            facets_allocator = new void *[size_one];
            unsigned int new_facets_array[size_three];
            unsigned int length_that_face;
            for(i = 0;i<nr_facets;i++){
                facets_allocator[i] = aligned_alloc(chunksize/8,length_of_face*chunksize/8);
                facets[i] = (chunktype*) facets_allocator[i];
                length_that_face = 0;
                for (j=0;j < nr_vertices;j++){
                    if (PyInt_AsLong(PyTuple_GetItem(PyTuple_GetItem(py_tuple, j),i))){
                        new_facets_array[length_that_face] = j;
                        length_that_face++;
                    }
                }
                char_from_array(new_facets_array, length_that_face, facets[i]);
            }
            facets_are_allocated = 1;
        }
        //conversions
void CombinatorialType::char_from_tuple(PyObject* py_tuple, chunktype *array1){//transforms a tuple to an array, the first entry gives the length of the tuple (and therefore the expected length of the array)
            unsigned int len, entry, position, value,i ;
            const unsigned int size_array = length_of_face*chunksize/64;
            uint64_t *array = new uint64_t [size_array]();
            len = PyTuple_Size(py_tuple);
            while (len--) {
                entry = (unsigned int) PyInt_AsLong(PyTuple_GetItem(py_tuple, len));
                value = entry % 64;
                position = entry/64;
                array[position] += vertex_to_bit_dictionary[value];
            }
            for (i=0;i<length_of_face;i++){
                load_register(array1[i],array[i*chunksize/64]);
            }
            delete(array);
        }
void CombinatorialType::char_from_incidence_tuple(PyObject* py_tuple, chunktype *array1){
            unsigned int len, entry, position, value,i ;
            const unsigned int size_array = length_of_face*chunksize/64;
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
            for (i=0;i<length_of_face;i++){
                load_register(array1[i],array[i*chunksize/64]);
            }
            delete(array);
        }
void CombinatorialType::char_from_array(unsigned int* input, unsigned int len, chunktype *array1){
            unsigned int entry, position, value,i ;
            const unsigned int size_array = length_of_face*chunksize/64;
            uint64_t *array = new uint64_t [size_array]();
            while (len--) {
                entry = input[len];
                value = entry % 64;
                position = entry/64;
                array[position] += vertex_to_bit_dictionary[value];
            }
            for (i=0;i<length_of_face;i++){
                load_register(array1[i],array[i*chunksize/64]);
            }
            delete(array);
        }
inline PyObject* CombinatorialType::tuple_from_f_vector(){
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
inline PyObject* CombinatorialType::tuple_from_edges(){
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
inline PyObject* CombinatorialType::tuple_from_ridges(){
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

CombinatorialType_ptr init_CombinatorialType(PyObject* py_tuple, unsigned int nr_vertices) {
    CombinatorialType_ptr C = new CombinatorialType(py_tuple,nr_vertices);
    return C;
}

CombinatorialType_ptr init_CombinatorialType(PyObject* py_tuple) {
    CombinatorialType_ptr C = new CombinatorialType(py_tuple);
    return C;
}

unsigned int dimension(CombinatorialType_ptr C){
  return (*C).get_dimension();
}

PyObject* edges(CombinatorialType_ptr C){
  return (*C).get_edges();
}

PyObject* f_vector(CombinatorialType_ptr C){
  return (*C).get_f_vector();
}

PyObject* ridges(CombinatorialType_ptr C){
  return (*C).get_ridges();
}

void delete_CombinatorialType(CombinatorialType_ptr C)
{
  delete(C);
}
