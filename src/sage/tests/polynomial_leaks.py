r"""
This file gather leaks that appeared in :trac:`27261`.
"""
import resource
import gc

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

def get_resources():
    gc.collect()
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

def test_leak(f, args, repeat, n_success):
    r"""
    Test that there is no memory when calling ``f(*args)``.

    Note that (lib)Singular has pre-allocated buckets, so we have to run a lot
    of iterations to fill those up first.
    """
    zeros = 0
    for i in range(repeat):
        n = f(*args)
        if n == 0:
            zeros += 1
            if zeros >= n_success:
                return
        else:
            zeros = 0
    raise RuntimeError('leak f={} args={} repeat={} n_success={}'.format(f, args, repeat, n_success))

def leak_call(base, nvars, repeat):
    r"""
    TESTS::

        sage: from sage.tests.polynomial_leaks import test_leak, leak_call
        sage: for base in [ZZ, QQ, GF(7), GF(49), AA, ZZ['y']]:
        ....:     test_leak(leak_call, (base, 2, 10), 40, 10)
        ....:     test_leak(leak_call, (base, 50, 10), 40, 10)
    """
    assert nvars >= 2
    R = PolynomialRing(base, 'x', nvars)
    gens = R.gens()
    p1 = R.zero()
    p2 = R.one()
    p3 = sum(gens)
    p4 = (gens[0] + gens[1])**10
    polys = [p1, p2, p3, p4]
    v1 = list(gens)
    v2 = [v1[0] + v1[1]] + v1[1:]
    v3 = [base.zero()] * nvars
    v4 = [(-1)**i * base.one() for i in range(nvars)]
    v5 = [base.gen() + gens[0]] + gens[1:]
    values = [v1, v2, v3, v4, v5]
    before = get_resources()
    for i in range(repeat):
        for p in polys:
            for v in values:
                _ = p(*v)
                _ = p(v)
    after = get_resources()
    return after - before

def leak_subs(base, nvars, repeat):
    r"""
    TESTS::

        sage: from sage.tests.polynomial_leaks import test_leak, leak_subs
        sage: for base in [ZZ, QQ, GF(7), GF(49), AA, ZZ['y']]:
        ....:     test_leak(leak_subs, (base, 2, 10), 50, 10)
        ....:     test_leak(leak_subs, (base, 50, 10), 50, 10)
    """
    R = PolynomialRing(base, 'x', nvars)
    gens = R.gens()
    p1 = R.zero()
    p2 = R.one()
    p3 = sum(gens)
    p4 = (gens[0] + gens[1])**10
    polys = [p1, p2, p3, p4]
    d1 = {str(g): g for g in R.gens()}
    d2 = {str(g): R.zero() for g in R.gens()}
    d3 = {str(g): base.one() for g in R.gens()}
    values = [d1, d2, d3]
    before = get_resources()
    for i in range(repeat):
        for p in polys:
            for d in values:
                _ = p.subs(**d)
                _ = p.subs(d)
    after = get_resources()
    return after - before
