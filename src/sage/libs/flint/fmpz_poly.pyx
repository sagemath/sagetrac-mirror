from fmpz cimport fmpz_sgn, fmpz_get_python_str

import itertools

cdef str fmpz_poly_get_python_str(
        const fmpz_poly_t poly, str name='x', bint latex=False, int base=10):

    if fmpz_poly_is_zero(poly):
        return '0'

    cdef slong i, sign
    cdef fmpz *coeff
    cdef list str_list = []
    cdef str coeff_str, sign_str

    for i in range(fmpz_poly_length(poly)):
        coeff = fmpz_poly_get_coeff_ptr(poly, i)
        sign = fmpz_sgn(coeff)
        if sign == 0:
            continue
        elif sign > 0:
            sign_str = '+'
            coeff_str = fmpz_get_python_str(coeff, base)
        else:
            sign_str = '-'
            coeff_str = fmpz_get_python_str(coeff, base)[1:]
        if i > 0:
            if coeff_str == '1':
                coeff_str = ''
            elif not latex:
                coeff_str += '*'
            coeff_str += name
        if i > 1:
            if latex and i >= 10:
                coeff_str += '^{{{i}}}'.format(i=i)
            else:
                coeff_str += '^{i}'.format(i=i)
        str_list.append((sign_str, coeff_str))

    itr = itertools.chain.from_iterable(reversed(str_list))

    sign_str = itr.next()
    if sign_str == '+':
        return ' '.join(itr)
    else:
        return sign_str + ' '.join(itr)
