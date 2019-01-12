//Copyright: see base.pyx

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

#if __POPCNT__
#define popcount(A) _mm_popcnt_u64(A)
#else
#define popcount(A) naive_popcount(A)
#endif

const unsigned int maxnumberedges = 16348;//^2 (the edges will be build as an array of arrays, such that we can save up to maxnumberedges*maxnumberedges edges, the number should contain a high power of two
const unsigned int maxnumberincidences = 16348;//^2 the maximal number of incidences between l-faces and k-faces

class CombinatorialPolytope {
    public:
        CombinatorialPolytope(PyObject* py_tuple, unsigned int nr_vertices_given);//initialization with a tuple of facets (each facet a tuple of vertices, vertices labeled 0,1,...)
        CombinatorialPolytope(PyObject* py_tuple);//initialization with an incidence matrix given as tuple of tuples
        ~CombinatorialPolytope();//cleanup to avoid memory leak
        unsigned int get_dimension();
        inline PyObject* get_f_vector();
        inline PyObject* get_edges();//returns a tuple of edges, each edges as tuple of vertices
        inline PyObject* get_ridges();//returns a tuple of ridges, each ridge as tuple of its to facets
        inline PyObject* get_faces(int face_dimension, unsigned int facet_repr);
        void record_all_faces();
        inline PyObject* get_incidences(int dimension_one, int dimension_two);
        inline unsigned long get_flag_number(PyObject* py_tuple); 
    private:
        int polar = 0;//in order to speed things up, we will consider the dual/polar whenever the number of vertices is smaller than the number of facets
        void **facets_allocator;
        int facets_are_allocated = 0;
        void **vertices_allocator = NULL;
        int vertices_are_allocated = 0;
        chunktype **facets = NULL;//facets as incidences of vertices
        chunktype **vertices = NULL;//vertices as incidenes of facets
        void ***newfaces_allocator = NULL;
        int newfaces_are_allocated = 0;
        chunktype ***newfaces = NULL, ***newfaces2 = NULL, **forbidden = NULL;
        unsigned long *f_vector = NULL;
        unsigned int nr_vertices, nr_facets, length_of_face, length_of_face_in_facet_repr;
        unsigned int dimension = 0;
        unsigned long nr_edges = 0, nr_ridges = 0, nr_incidences = 0;
        unsigned int **edges = new unsigned int *[maxnumberedges];
        unsigned int **ridges = new unsigned int *[maxnumberedges];
        unsigned long **incidences = new unsigned long *[maxnumberincidences]();
        unsigned int edgemode = 0;
        
        unsigned int *allfaces_are_allocated = NULL;//this should be set to 1 if it is just allocated and to 2 if the corresponding faces have been recorded already
        void ***allfaces_allocator = NULL;
        chunktype ***allfaces = NULL;
        void ***allfaces_facet_repr_allocator = NULL;
        chunktype ***allfaces_facet_repr = NULL;
        unsigned long *allfaces_counter = NULL;

        inline void intersection(chunktype *A, chunktype *B, chunktype *C);//will set C to be the intersection of A and B
        inline int is_subset(chunktype *A, chunktype *B);//returns 1 if A is a proper subset of B, otherwise returns 0
        inline int is_subset_facet_repr(chunktype *A, chunktype *B);//as above just in facet_repr
        inline unsigned int CountFaceBits(chunktype* A1);//counts the number of vertices in a face by counting bits set to one
        inline unsigned int CountFaceBits_facet_repr(chunktype* A1);
	inline void add_edge(chunktype *face);//adds an edge to the edges list
        inline void add_edge(unsigned int one, unsigned int two);//adds an edge to the edges list given as its two vertices
        inline void add_ridge(unsigned int one, unsigned int two);//adds a ridge to the ridge list given as its two facets
        inline void add_incidence(unsigned long one, unsigned long two);//adds an incidence to the list of incidences, where one and two correspond to the number of the faces according to allfaces resp. vertices/facets
        
        //get_next_level intersects the first 'lenfaces' faces of 'faces' with the 'face_to_intersect'-th face of faces and stores the result in 'nextfaces'
        //then determines which ones are exactly of one dimension less by considering containment
        //newfaces2 will point at those of exactly one dimension less which are not contained in any of the faces in 'forbidden'
        //returns the number of those faces
        inline unsigned int get_next_level(chunktype **faces, unsigned int lenfaces, unsigned int face_to_intersect, chunktype **nextfaces, chunktype **nextfaces2, unsigned int nr_forbidden);
        unsigned int calculate_dimension(chunktype **faces, unsigned int nr_faces);
        void calculate_ridges();
        void get_f_vector_and_edges();
        void get_f_vector_and_edges(chunktype **faces, unsigned int dimension, unsigned int nr_faces, unsigned int nr_forbidden);
        void record_faces(unsigned int lowest_dimension);
        void record_faces(chunktype **faces, unsigned int current_dimension, unsigned int nr_faces, unsigned int nr_forbidden, unsigned int lowest_dimension);
        inline void record_face(chunktype *face, unsigned int current_dimension);
        inline void record_face_facet_repr(chunktype *face, unsigned int current_dimension);
        void vertex_facet_incidences();
        void vertex_facet_incidences(chunktype *array1, unsigned int nr_facet);
        unsigned long get_flag_number(unsigned int *array, unsigned int len);

        //initialization
        void get_facets_from_tuple(PyObject* py_tuple);
        void get_vertices_from_tuple(PyObject* py_tuple);
        void get_facets_or_vertices_from_tuple(PyObject* py_tuple, chunktype** facets_or_vertices, unsigned int flip, unsigned int nr_vertices_given, unsigned int nr_facets_given, unsigned int facet_repr);
        void get_facets_from_incidence_matrix(PyObject* py_tuple);
        void get_vertices_from_incidence_matrix(PyObject* py_tuple);
        void get_facets_or_vertices_from_incidence_matrix(PyObject* py_tuple, chunktype** facets_or_vertices, unsigned int flip, unsigned int nr_vertices_given, unsigned int nr_facets_given, unsigned int facet_repr);
        
        
        //conversions
        void char_from_tuple(PyObject* py_tuple, chunktype *array1, unsigned int facet_repr);
        inline PyObject* tuple_from_char(chunktype *array1, unsigned int facet_repr);
        inline PyObject* tuple_from_faces(chunktype **array1, unsigned int len, unsigned int facet_repr);
        void char_from_incidence_tuple(PyObject* py_tuple, chunktype *array1, unsigned int facet_repr);
        void char_from_array(unsigned int* input, unsigned int len, chunktype *array1, unsigned int facet_repr);
        inline PyObject* tuple_from_f_vector();
        inline PyObject* tuple_from_edges();
        inline PyObject* tuple_from_ridges();
        inline PyObject* tuple_from_incidences(unsigned int twisted);
        inline PyObject* n_trivial_tuples(unsigned int nr_elements); //Builds a tuple of the form ((1,), (2,), (3,), ... , (nr_elements,))
        inline PyObject* one_trivial_tuple(unsigned int nr_elements); //Builds a tuple of the form ((1, 2, 3, ... , nr_elements),)
        
        //allocation and deallocation
        void allocate_facets();
        void deallocate_facets();
        void allocate_vertices();
        void deallocate_vertices();
        void allocate_newfaces();
        void deallocate_newfaces();
        void allocate_allfaces(unsigned int dimension_to_allocate);//allocates allfaces in a certain dimension, must be smaller than dimension and at least 1, if dimension is 0 will allocate all dimensions
        void deallocate_allfaces();
        
};

typedef CombinatorialPolytope* CombinatorialPolytope_ptr;

CombinatorialPolytope_ptr init_CombinatorialPolytope(PyObject* py_tuple, unsigned int nr_vertices);//initialize by facets as tuples of vertices
CombinatorialPolytope_ptr init_CombinatorialPolytope(PyObject* py_tuple);//initialize by incidence_matrix

unsigned int dimension(CombinatorialPolytope_ptr C);
PyObject* f_vector(CombinatorialPolytope_ptr C);
PyObject* edges(CombinatorialPolytope_ptr C);
PyObject* ridges(CombinatorialPolytope_ptr C);
PyObject* incidences(CombinatorialPolytope_ptr C, int dimension_one, int dimension_two);
void record_all_faces(CombinatorialPolytope_ptr C);
PyObject* get_faces(CombinatorialPolytope_ptr C, int dimension, unsigned int facet_repr);
unsigned long get_flag(CombinatorialPolytope_ptr C, PyObject* py_tuple);

void delete_CombinatorialPolytope(CombinatorialPolytope_ptr);

#endif
