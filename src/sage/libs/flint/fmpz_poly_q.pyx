from fmpz_poly cimport fmpz_poly_get_python_str, fmpz_poly_is_one

cdef str fmpz_poly_q_get_python_str(
        const fmpz_poly_q_t func, str name='x', bint latex=False, int base=10):
    cdef str num, den

    num = fmpz_poly_get_python_str(
            fmpz_poly_q_numref(func), name, latex, base)

    if fmpz_poly_is_one(fmpz_poly_q_denref(func)):
        return num

    den = fmpz_poly_get_python_str(
            fmpz_poly_q_denref(func), name, latex, base)

    if latex:
        return r'\frac{{{num}}}{{{den}}}'.format(num=num, den=den)

    if ' ' in num:
        num = '(' + num + ')'
    if ' ' in den or '*' in den:
        den = '(' + den + ')'

    return num + '/' + den
