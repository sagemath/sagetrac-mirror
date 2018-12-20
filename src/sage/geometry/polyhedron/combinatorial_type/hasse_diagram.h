#ifndef HASSEDIAGRAM__H
#define HASSEDIAGRAM__H


//as of now, 512bit does not have something like _mm256_testc_si256, which is the bottle neck of this function, so it does not make sene to implement it
//512 bit commands
//#if __AVX512F__
//#define chunktype __m512i
//#define bitwise_intersection(one,two) _mm512_and_si512((one),(two))

#if __AVX2__
//256 bit commands
#include <immintrin.h>
#define chunktype __m256i
const unsigned int chunksize = 256;
#define bitwise_intersection(one,two) _mm256_and_si256((one),(two)) // this is supposed to something as (one) & (two)
#define bitwise_is_not_subset(one,two) !_mm256_testc_si256((two),(one)) // this is supposed to something as (one) & ~(two)
#define store_register(one,two) _mm256_storeu_si256((__m256i*)&(one),(two)) //this is supposed to be something as one = two, where two is a register
#define load_register(one,two) (one) = _mm256_loadu_si256((const __m256i*)&(two)) //this is supposed to be somethign as one = two, where one is a register
#define leading_zero_count(one) leading_zero_workaround(one) //the workaround is not extremely fast, but it is not needed often (only for edges)
#define trailing_zero_count(one) trailing_zero_workaround(one) //the workaround is not extremely fast, but it is not needed often

#else
#if __SSE4_1__
//128 bit commands
#include <emmintrin.h>
#include <smmintrin.h>
#define chunktype __m128i
const unsigned int chunksize = 128;
#define bitwise_intersection(one,two) _mm_and_si128((one),(two))
#define bitwise_is_not_subset(one,two) !_mm_testc_si128((two),(one))
#define store_register(one,two) _mm_storeu_si128((__m128i*)&(one),(two))
#define load_register(one,two) (one) = _mm_loadu_si128((const __m128i*)&(two))
#define leading_zero_count(one) leading_zero_workaround(one) 
#define trailing_zero_count(one) trailing_zero_workaround(one) 

#else
//64 bit commands
#define chunktype uint64_t
const unsigned int chunksize = 64;
#define bitwise_intersection(one,two) (one) & (two)
#define bitwise_is_not_subset(one,two) (one) & ~(two)
#define store_register(one,two) one = two
#define load_register(one,two) one = two
#define leading_zero_count(one) leading_zero_naive3(one)
#define trailing_zero_count(one) trailing_zero_naive3(one)

#endif
#endif

const unsigned int maxnumberedges = 16348;//^2 (the edges will be build as an array of arrays, such that I can save up to maxnumberedges*maxnumberedges edges, the number should contain a high power of two


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
