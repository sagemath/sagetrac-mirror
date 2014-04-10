from libc.stdint cimport uint8_t, uint_fast64_t

from libcpp.vector cimport vector

cdef extern from "borie.hpp":
    const unsigned N

    cdef cppclass SGroup nogil:
        cppclass type
    unsigned SGroup_N "SGroup::N"
    unsigned SGroup_vect_len "SGroup::vect_len"
    cppclass SGroup_type "SGroup::type":
        uint8_t &operator[](int)
        bint operator==(const SGroup_type &vp)
    SGroup_type SGroup_one "SGroup::one" ()
    void SGroup_mult "SGroup::mult" (SGroup_type &a, SGroup_type &a, SGroup_type &b)

    bint is_canonical(const vector[vector[SGroup_type] ] sgs,
                      const SGroup_type v)
