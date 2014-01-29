"""
    sage: import sage.libs.decoding
    sage: sage.libs.decoding.f()
    1
"""

#include <stdlib.h>
# #include <decoding/decoding.h>
cdef extern from "decoding/decoding.h":
    pass

# This file is decoding.c from the decoding sources. We include it
# instead of linking against decoding.a for now, because there is some
# linker incompatibility preventing linking this static library into a
# dynamically loaded library
cdef extern from "decoding_base.c":
    pass

# #define RING_NAME ring
# #define PRIME 101
# #include <decoding/rings/gfp_word.c>
# #include <decoding/algos.c>
#
# This .h file contains the above four lines. Just because I don't
# know how to insert a #define in a Cython file
cdef extern from "decoding_gfp_word.h":
    ctypedef struct ring_vec:
        pass
    ctypedef struct ring_rs_code:
        pass

    void ring_ring_init()
    void ring_vec_inits(ring_vec, ...)
    void ring_rs_code_init_with_support(ring_rs_code, long, long)
    void ring_rs_code_random_codeword(ring_vec, ring_rs_code)
    void ring_vec_random_hamming_weight(ring_vec, long, long)
    void ring_vec_add(ring_vec, ring_vec, ring_vec)
    long ring_rs_code_sudan_koetter(ring_vec **, ring_rs_code, ring_vec, long)


def f():
    # /* This exemple shows how to use the
    #    Guruswami-sudan algorithm. */

    # /* Generate code for GF(PRIME) where PRIME
    #    holds within a machine word. */




    # int main(void) {


    #   /* The ``uma" type (Unaligned Memory Access) shall be
    #      used everywhere you need indices. In particular,
    #      all function in decoding returning polynomials
    #      degrees, valuations, number of elements of a certain
    #      set are all of type ``uma". It usually corresponds
    #      to long, but do not assume that. You can safely
    #      assume that all int can be casted to uma without
    #      problem. */

    #   uma nc,i;
    cdef long nc, i

    #   /* We need three vectors: the message -------> m 
    #                             the error ---------> e
    #                             the received word -> y */
    # ring_vec m,e,y;
    cdef ring_vec m,e,y

    #   /* We need an array of vectors. It will receive
    #      the list of vectors returned by the Guruswami-Sudan
    #      algorithm. */
    # ring_vec *c;
    cdef ring_vec *c


    #   /* We need a Reed-Solomon code over ring. */
    #   ring_rs_code rs;
    cdef ring_rs_code rs

    #   /* We initialize the ring. */
    # ring_ring_init();
    ring_ring_init()

    #   /* We initialize the vectors. */
    #   ring_vec_inits(m,e,y,NULL);
    ring_vec_inits(m,e,y,NULL)

    #   /* We make a Reed-Solomon code of parameters [50,10,41] over 
    #      our ring. Note that all the functions are prefixed by
    #      ``ring_". */
    #   ring_rs_code_init_with_support(rs,50,10);
    ring_rs_code_init_with_support(rs,50,10)

    #   /* We generate a random codeword. This is our message m. */
    #   ring_rs_code_random_codeword(m,rs);
    ring_rs_code_random_codeword(m,rs)

    #   /* We generate an error vector e of length 50 and weight 22. */
    #   ring_vec_random_hamming_weight(e,50,22);
    ring_vec_random_hamming_weight(e,50,22)

    #   /* We build the received word y = m + e. */
    #   ring_vec_add(y,m,e);
    ring_vec_add(y,m,e)

    #   /* We call the Sudan algorithm to correct y. We tell the algorithm
    #      that we want to correct the maximum possible number of errors.
    #      The list of codewords returned by our algorithm will be put
    #      in c. */
    #   nc = ring_rs_code_sudan_koetter(&c,rs,y,SUDAN_MAX_ERRORS);
    nc = ring_rs_code_sudan_koetter(&c,rs,y,-1)
    return nc

    #   /* nc contains the number of codewords in the array c. Then
    #      we print the original message m and the list returned by
    #      the algorithm. */
    #   ring_decoding_printf(stdout,"Original message: %v\n",m);
    #   for( i = 0 ; i < nc ; i++ ) {
    #     ring_decoding_printf(stdout,"Codeword %d: %v\n",i,c[i]);
    #   }

    #   return 0;

    # }


