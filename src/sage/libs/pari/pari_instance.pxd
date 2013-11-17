include 'decl.pxi'

from sage.libs.pari.gen cimport gen
from sage.structure.parent cimport Parent

cdef class PariInstance(Parent):
    cdef unsigned long num_primes
    cdef unsigned long default_prec
    cdef int init_stack(self, size_t size) except -1
    cdef void clear_stack(self)
    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address)
    cdef long get_var(self, v)
    cdef gen _new_gen(self, GEN x)
    cdef gen new_gen(self, GEN x)
    cdef gen new_ref(self, GEN g, gen x)
    cdef gen new_gen_from_int(self, int value)
    cdef gen new_gen_from_double(self, double value)
    cdef gen new_gen_from_mpz_t(self, mpz_t value)
    cdef gen new_gen_from_mpq_t(self, mpq_t value)
    cdef gen new_gen_from_padic(self, long ordp, long relprec, mpz_t prime, mpz_t p_pow, mpz_t unit)
    cdef gen new_gen_from_mpz_t_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf)
    cdef gen new_gen_from_mpq_t_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)
    cdef gen new_t_POL_from_int_star(self, int *vals, int length, long varnum)
    cdef GEN toGEN(self, x) except NULL
    cdef inline GEN _new_GEN_from_double(self, double)
    cdef inline GEN _new_GEN_from_mpz_t(self, mpz_t value)
    cdef inline GEN _new_GEN_from_mpq_t(self, mpq_t value)
    cdef GEN _new_GEN_from_mpz_t_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc)
    cdef GEN _new_GEN_from_mpz_t_matrix_rotate90(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc)
    cdef GEN _new_GEN_from_mpq_t_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)

cpdef int prec_bits_to_dec(int prec_in_bits)
cpdef int prec_dec_to_bits(int prec_in_dec)
cpdef int prec_bits_to_words(int prec_in_bits)
cpdef int prec_words_to_bits(int prec_in_words)
cpdef int prec_dec_to_words(int prec_in_dec)
cpdef int prec_words_to_dec(int prec_in_words)
