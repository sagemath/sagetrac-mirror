"""
Some utility routines for the hypergeometric motives package that benefit
significantly from Cythonization.
"""

from cpython cimport array

cpdef hgm_coeffs(long p, int f, gamma, array.array m, int D,
                 array.array gtab0, gtab1):
    r"""
    Compute coefficients for the hypergeometric trace formula.

    This function is not intended for direct user access.

    TESTS::

        sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
        sage: import array
        sage: from sage.modular.hypergeometric_misc import hgm_coeffs
        sage: H = Hyp(cyclotomic=([3],[4]))
        sage: H.euler_factor(2, 7, cache_p=True)
        7*T^2 - 3*T + 1
        sage: gamma = H.gamma_array()
        sage: gtable = H.gauss_table(7, 1, 2)
        sage: m = array.array('i', [0]*6)
        sage: D = 1
        sage: hgm_coeffs(7, 1, gamma, m, D, gtable)
        [7, 2*7, 6*7, 7, 6, 4*7]
    """
    cdef int gl, j, k, l, r, r1, v, gv, prec
    cdef long i, q1

    q1 = p ** f - 1
    gl = len(gamma)
    cdef array.array gamma_array1 = array.array('i', gamma.keys())
    cdef array.array gamma_array2 = array.array('i', gamma.values())
    cdef array.array r_array = array.array('i', [0]) * gl

    R = gtab1[0].parent()
    prec = R.precision_cap()
    ans = []
    for r in range(q1):
        # First determine whether this term is forced to be zero
        # for divisibility reasons. If so, skip the p-adic arithmetic.
        i = 0
        for k in range(gl):
            v = gamma_array1[k]
            gv = gamma_array2[k]
            r1 = v * r % q1
            r_array[k] = r1
            i += gtab0[r1] * gv
        i //= (p - 1)
        l = i + f * (D + m[0] - m[r])
        if l >= prec:
            ans.append(R.zero())
            continue
        u = R.one()
        u1 = R.one()
        if i % 2: u = -u
        for k in range(gl):
            gv = gamma_array2[k]
            r1 = r_array[k]
            if gv > 0:
                for j in range(gv): u *= gtab1[r1]
            else:
                for j in range(-gv): u1 *= gtab1[r1]
        ans.append((u / u1) << l)
    return ans
