include 'sage/ext/interrupt.pxi'
include 'sage/ext/stdsage.pxi'

cdef str fmpz_get_python_str(const fmpz_t f, int base=10):
    cdef char *s
    cdef str res

    s = <char *>sage_malloc(fmpz_sizeinbase(f, base) + 2)
    if s == NULL:
        raise MemoryError("unable to allocate enough memory for the string "
                "representation of an integer")

    if not sig_on_no_except():
        sage_free(s)
        cython_check_exception()
    fmpz_get_str(s, base, f)
    sig_off()

    res = str(s)

    sage_free(s)

    return res
