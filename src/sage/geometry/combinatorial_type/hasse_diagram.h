#ifndef HASSEDIAGRAM__H
#define HASSEDIAGRAM__H

#include <immintrin.h>
#define chunktype __m256i
const unsigned int maxnumberedges = 16348;

class CombinatorialType {
    public:
        CombinatorialType(PyObject* py_tuple, unsigned int nr_vertices_given);
        CombinatorialType(PyObject* py_tuple);
        ~CombinatorialType();
        unsigned int get_dimension();
        inline PyObject* get_f_vector();
        inline PyObject *get_edges();
        inline PyObject *get_ridges();
    private:
        int polar = 0;
        void **facets_allocator;
        int facets_are_allocated = 0;
        chunktype **facets = NULL, ***newfaces = NULL, ***newfaces2 = NULL, **forbidden = NULL;
        void ***newfaces_allocator;
        int newfaces_are_allocated = 0;
        unsigned long *f_vector = NULL;
        unsigned int nr_vertices, nr_facets, length_of_face;
        unsigned int dimension = 0;
        unsigned long nr_edges = 0, nr_ridges = 0;
        unsigned int **edges = new unsigned int *[maxnumberedges];
        unsigned int **ridges = new unsigned int *[maxnumberedges];
        unsigned int edgemode = 0;
        inline void intersection(chunktype *A, chunktype *B, chunktype *C);
        inline int is_subset(chunktype *A, chunktype *B);
        inline unsigned int CountFaceBits(chunktype* A1);
        inline void add_edge(chunktype *face);
        inline void add_edge(unsigned int one, unsigned int two);
        inline void add_ridge(unsigned int one, unsigned int two);
        inline unsigned int get_next_level(chunktype **faces, unsigned int lenfaces, unsigned int face_to_intersect, chunktype **nextfaces, chunktype **nextfaces2, unsigned int nr_forbidden);
        unsigned int calculate_dimension(chunktype **faces, unsigned int nr_faces);
        void calculate_ridges();
        void get_f_vector_and_edges(chunktype **faces, unsigned int dimension, unsigned int nr_faces, unsigned int nr_forbidden);
        void calculate_f_vector();
        void Polar_Init(PyObject* py_tuple, unsigned int nr_vertices_given);
        void Polar_Init(PyObject* py_tuple);
        //conversions
        void char_from_tuple(PyObject* py_tuple, chunktype *array1);
        void char_from_incidence_tuple(PyObject* py_tuple, chunktype *array1);
        void char_from_array(unsigned int* input, unsigned int len, chunktype *array1);
        inline PyObject* tuple_from_f_vector();
        inline PyObject* tuple_from_edges();
        inline PyObject* tuple_from_ridges();
};

typedef CombinatorialType* CombinatorialType_ptr;

CombinatorialType_ptr init_CombinatorialType(PyObject* py_tuple, unsigned int nr_vertices);//initialize by facets as tuples of vertices
CombinatorialType_ptr init_CombinatorialType(PyObject* py_tuple);//initialize by incidence_matrix

unsigned int dimension(CombinatorialType_ptr C);
PyObject* f_vector(CombinatorialType_ptr C);
PyObject* edges(CombinatorialType_ptr C);
PyObject* ridges(CombinatorialType_ptr C);

void delete_CombinatorialType(CombinatorialType_ptr);

#endif
