from sage.misc.misc import cputime
from sage.tests.benchmark import Benchmark

from sage.all import *

def benchmark(n=-1):
    """
    Run a well-chosen range of Sage commands and record the time it
    takes for each to run.

    INPUT:
        n -- int (default: -1) the benchmark number; the default
             of -1 runs all the benchmarks.
    OUTPUT:
        list -- summary of timings for each benchmark.
        int -- if n == -1, also return the total time

    EXAMPLE:
        sage: from sage.tests.benchmarks.misc import *
        sage: _ = benchmark()
        Running benchmark 0
        Benchmark 0: Factor the following polynomial over
            the rational numbers: (x^97+19*x+1)*(x^103-19*x^97+14)*(x^100-1)
        Time: ... seconds
        Running benchmark 1
        Find the Mordell-Weil group of the elliptic curve 5077A using mwrank
        Time: ... seconds
        Running benchmark 2
        Some basic arithmetic with very large Integer numbers: '3^1000001 * 19^100001
        Time: ... seconds
        Running benchmark 3
        Some basic arithmetic with very large Rational numbers: '(2/3)^100001 * (17/19)^100001
        Time: ... seconds
        Running benchmark 4
        Rational polynomial arithmetic using Sage. Compute (x^29+17*x-5)^200.
        Time: ... seconds
        Running benchmark 5
        Rational polynomial arithmetic using Sage. Compute (x^19 - 18*x + 1)^50 one hundred times.
        Time: ... seconds
        Running benchmark 6
        Compute the p-division polynomials of y^2 = x^3 + 37*x - 997 for primes p < 40.
        Time: ... seconds
        Running benchmark 7
        Compute the Mordell-Weil group of y^2 = x^3 + 37*x - 997.
        Time: ... seconds
        Running benchmark 8

    """

    if isinstance(n, list):
        t = cputime()
        v = [benchmark(m) for m in n]
        return v, cputime(t)

    if n != -1:
        print "Running benchmark %s"%n
        try:
            b = eval("Bench%d()"%n)
        except NameError:
            raise RuntimeError, "no benchmark %s"%n
        t = b.sage()
        print b
        print "Time: %s seconds"%t
        return (n, t, str(b))

    t = cputime()
    m = 0
    v = []
    while True:
        try:
            v.append(benchmark(m))
            m += 1
        except RuntimeError:
            break
    return v, cputime(t)

class Bench0(Benchmark):
    """
    EXAMPLE:
        sage: from sage.tests.benchmarks.misc import *
        sage: print Bench0()
        Benchmark 0: Factor the following polynomial over
            the rational numbers: (x^97+19*x+1)*(x^103-19*x^97+14)*(x^100-1)

    """

    def __init__(self):
        self.repr_str = """Benchmark 0: Factor the following polynomial over
        the rational numbers: (x^97+19*x+1)*(x^103-19*x^97+14)*(x^100-1)"""

    def sage(self):
        x = polygen(QQ,"x")
        f = (x**97+19*x+1)*(x**103-19*x**97+14)*(x**100-1)
        t = cputime()
        F = f.factor()
        return cputime(t)

class Bench1(Benchmark):
    """
    EXAMPLE:
        sage: from sage.tests.benchmarks.misc import *
        sage: print Bench1()
        Find the Mordell-Weil group of the elliptic curve 5077A using mwrank

    """

    def __init__(self):
        self.repr_str = """Find the Mordell-Weil group of the elliptic curve 5077A using mwrank"""

    def sage(self):
        E = mwrank_EllipticCurve([0, 0, 1, -7, 6])
        t = cputime()
        g = E.gens()
        return cputime(t)

class Bench2(Benchmark):
    """
    EXAMPLE:
        sage: from sage.tests.benchmarks.misc import *
        sage: print Bench2()
        Some basic arithmetic with very large Integer numbers: '3^1000001 * 19^100001

    """

    def __init__(self):
        self.repr_str = """Some basic arithmetic with very large Integer numbers: '3^1000001 * 19^100001"""

    def sage(self):
        t = cputime()
        a = ZZ(3)**1000001 * ZZ(19)**100001
        return cputime(t)

class Bench3(Benchmark):
    """
    EXAMPLE:
        sage: from sage.tests.benchmarks.misc import *
        sage: print Bench3()
        Some basic arithmetic with very large Rational numbers: '(2/3)^100001 * (17/19)^100001

    """

    def __init__(self):
        self.repr_str = """Some basic arithmetic with very large Rational numbers: '(2/3)^100001 * (17/19)^100001"""

    def sage(self):
        t = cputime()
        a = QQ('2/3')**100001 * QQ('17/19')**100001
        return cputime(t)

class Bench4(Benchmark):
    """
    EXAMPLE:
        sage: from sage.tests.benchmarks.misc import *
        sage: print Bench4()
        Rational polynomial arithmetic using Sage. Compute (x^29+17*x-5)^200.

    """

    def __init__(self):
        self.repr_str = """Rational polynomial arithmetic using Sage. Compute (x^29+17*x-5)^200."""

    def sage(self):
        x = PolynomialRing(QQ, 'x').gen()
        t = cputime()
        f = x**29 + 17*x-5
        a = f**200
        return cputime(t)

class Bench5(Benchmark):
    """
    EXAMPLE:
        sage: from sage.tests.benchmarks.misc import *
        sage: print Bench5()
        Rational polynomial arithmetic using Sage. Compute (x^19 - 18*x + 1)^50 one hundred times.

    """

    def __init__(self):
        self.repr_str = """Rational polynomial arithmetic using Sage. Compute (x^19 - 18*x + 1)^50 one hundred times."""

    def sage(self):
        x = PolynomialRing(QQ, 'x').gen()
        t = cputime()
        f = x**19 - 18*x + 1
        w = [f**50 for _ in range(100)]
        return cputime(t)

class Bench6(Benchmark):
    """
    EXAMPLE:
        sage: from sage.tests.benchmarks.misc import *
        sage: print Bench6()
        Compute the p-division polynomials of y^2 = x^3 + 37*x - 997 for primes p < 40.

    """

    def __init__(self):
        self.repr_str = """Compute the p-division polynomials of y^2 = x^3 + 37*x - 997 for primes p < 40."""

    def sage(self):
        E = EllipticCurve([0,0,0,37,-997])
        t = cputime()
        for p in [2,3,5,7,11,13,17,19,23,29,31,37]:
            f = E.division_polynomial(p)
        return cputime(t)

class Bench7(Benchmark):
    """
    EXAMPLE:
        sage: from sage.tests.benchmarks.misc import *
        sage: print Bench7()
        Compute the Mordell-Weil group of y^2 = x^3 + 37*x - 997.

    """

    def __init__(self):
        self.repr_str = """Compute the Mordell-Weil group of y^2 = x^3 + 37*x - 997."""

    def sage(self):
        E = EllipticCurve([0,0,0,37,-997])
        t = cputime()
        G = E.gens()
        return cputime(t)
