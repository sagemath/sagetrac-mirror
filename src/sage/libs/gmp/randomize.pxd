"""
Generate random rationals in Sage
"""

from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *
from sage.misc.randstate cimport randstate, current_randstate, SAGE_RAND_MAX

###########################

cdef inline void mpq_randomize_entry(mpq_t x, mpz_t num_bound, mpz_t den_bound):
    cdef randstate rstate = current_randstate()
    mpz_urandomm(mpq_numref(x), rstate.gmp_state, num_bound)
    mpz_urandomm(mpq_denref(x), rstate.gmp_state, den_bound)
    if mpz_sgn(mpq_denref(x)) == 0:
        mpz_set_si(mpq_denref(x),1)
    if rstate.c_random() % 2:
        mpz_mul_si(mpq_numref(x), mpq_numref(x), -1)
    mpq_canonicalize(x)

cdef inline void mpq_randomize_entry_as_int(mpq_t x, mpz_t bound):
    cdef randstate rstate = current_randstate()
    mpz_urandomm(mpq_numref(x), rstate.gmp_state, bound)
    mpz_set_si(mpq_denref(x), 1)
    if rstate.c_random() % 2:
        mpz_mul_si(mpq_numref(x), mpq_numref(x), -1)

cdef inline void mpq_randomize_entry_recip_uniform(mpq_t x):
    cdef randstate rstate = current_randstate()
    # Numerator is selected the same way as ZZ.random_element();
    # denominator is selected in a similar way, but
    # modified to give only positive integers.  (The corresponding
    # probability distribution is $X = \mbox{trunc}(1/R)$, where R
    # varies uniformly between 0 and 1.)
    cdef int den = rstate.c_random() - SAGE_RAND_MAX/2
    if den == 0: den = 1
    mpz_set_si(mpq_numref(x), (SAGE_RAND_MAX/5*2) / den)
    den = rstate.c_random()
    if den == 0: den = 1
    mpz_set_si(mpq_denref(x), SAGE_RAND_MAX / den)
    mpq_canonicalize(x)
