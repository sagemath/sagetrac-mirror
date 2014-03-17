"""
Benchmarks for matrices

This file has many functions for computing timing benchmarks
of various methods for random matrices with given bounds for
the entries.  The systems supported are Sage and Magma.

The basic command syntax is as follows::

    sage: import sage.matrix.benchmark as b
    sage: print "starting"; import sys; sys.stdout.flush(); b.report([b.det_ZZ], 'Test', systems=['sage'])
    starting...
    ======================================================================
              Test
    ======================================================================
    ...
    ======================================================================
"""
from sage.tests.benchmarks.benchmark import Benchmark
from sage.interfaces.magma import magma

from constructor import random_matrix, Matrix
from sage.rings.all import ZZ, QQ, GF
from sage.misc.misc import alarm, cancel_alarm, cputime
from sage.ext.c_lib import AlarmInterrupt

verbose = False

timeout = 60

def report(F, title, systems = ['sage', 'magma'], **kwds):
    """
    Run benchmarks with default arguments for each function in the list F.

    INPUT:

    - ``F`` - a list of callables used for benchmarking
    - ``title`` - a string describing this report
    - ``systems`` - a list of systems (supported entries are 'sage' and 'magma')
    - ``**kwds`` - keyword arguments passed to all functions in ``F``

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: print "starting"; import sys; sys.stdout.flush(); b.report([b.det_ZZ], 'Test', systems=['sage'])
        starting...
        ======================================================================
                  Test
        ======================================================================
        ...
        ======================================================================
    """
    import os
    if len(systems) > 2:
        raise NotImplementedError, "at most two systems ('sage' or 'magma')"
    print '='*70
    print ' '*10 + title
    print '='*70
    os.system('uname -a')
    print '\n'
    for f in F:
        print "-"*70
        print f.__doc__.strip()
        print ('%15s'*len(systems))%tuple(systems)
        w = []
        for s in systems:
            alarm(timeout)
            try:
                # t = f(**kwds).run(systems=[s])
                t = eval('f(**kwds).%s()'%s)
            except AlarmInterrupt:
                t = -timeout
            cancel_alarm()
            w.append(float(t))
        if len(w) > 1:
            if w[1] == 0:
                w.append(0.0)
            else:
                w.append(w[0]/w[1])

        w = tuple(w)
        print ('%15.3f'*len(w))%w
    print '='*70


#######################################################################
# Dense Benchmarks over ZZ
#######################################################################

def report_ZZ(**kwds):
    """
    Reports all the benchmarks for integer matrices and few
    rational matrices.

    INPUT:

    - ``**kwds`` - passed through to :func:`report`

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: print "starting"; import sys; sys.stdout.flush(); b.report_ZZ(systems=['sage'])  # long time (15s on sage.math, 2012)
        starting...
        ======================================================================
        Dense benchmarks over ZZ
        ======================================================================
        ...
        ======================================================================
    """
    F = [Vecmat_ZZ, Rank_ZZ, Rank2_ZZ, Charpoly_ZZ, Smithform_ZZ, Det_ZZ,
        Det_QQ, MatrixMultiply_ZZ, MatrixAdd_ZZ, MatrixAdd_ZZ_2, Nullspace_ZZ]

    title = 'Dense benchmarks over ZZ'
    report(F, title, **kwds)

# Integer Nullspace

class Nullspace_ZZ(Benchmark):
    """
    Nullspace over ZZ:
    Given a n+1 x n matrix over ZZ with random entries
    between min and max, compute the nullspace.

    INPUT:

    - ``n`` - matrix dimension (default: ``200``)
    - ``min`` - minimal value for entries of matrix (default: ``0``)
    - ``max`` - maximal value for entries of matrix (default: ``2**32``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.nullspace_ZZ(200)
        sage: tm = b.nullspace_ZZ(200, system='magma')  # optional - magma
    """

    def __init__(self, n=200, min=0, max=2**32):
        self.n, self.min, self.max = n, min, max

    def sage(self):
        A = random_matrix(ZZ, self.n+1, self.n, x=self.min, y=self.max+1).change_ring(QQ)
        t = cputime()
        v = A.kernel()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := RMatrixSpace(RationalField(), n+1,n)![Random(%s,%s) : i in [1..n*(n+1)]];
t := Cputime();
K := Kernel(A);
s := Cputime(t);
"""%(self.n,self.min,self.max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

class Charpoly_ZZ(Benchmark):
    """
    Characteristic polynomial over ZZ:
    Given a n x n matrix over ZZ with random entries between min and
    max, compute the charpoly.

    INPUT:

    - ``n`` - matrix dimension (default: ``100``)
    - ``min`` - minimal value for entries of matrix (default: ``0``)
    - ``max`` - maximal value for entries of matrix (default: ``9``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.charpoly_ZZ(100)
        sage: tm = b.charpoly_ZZ(100, system='magma')  # optional - magma
    """

    def __init__(self, n=100, min=0, max=9):
        self.n, self.min, self.max = n, min, max

    def sage(self):
        A = random_matrix(ZZ, self.n, self.n, x=self.min, y=self.max+1)
        t = cputime()
        v = A.charpoly()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
t := Cputime();
K := CharacteristicPolynomial(A);
s := Cputime(t);
"""%(self.n,self.min,self.max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

class Rank_ZZ(Benchmark):
    """
    Rank over ZZ:
    Given a n x (n+10) matrix over ZZ with random entries
    between min and max, compute the rank.

    INPUT:

    - ``n`` - matrix dimension (default: ``700``)
    - ``min`` - minimal value for entries of matrix (default: ``0``)
    - ``max`` - maximal value for entries of matrix (default: ``9``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.rank_ZZ(300)
        sage: tm = b.rank_ZZ(300, system='magma')  # optional - magma
    """

    def __init__(self, n=700, min=0, max=9):
        self.n, self.min, self.max = n, min, max

    def sage(self):
        A = random_matrix(ZZ, self.n, self.n+10, x=self.min, y=self.max+1)
        t = cputime()
        v = A.rank()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := RMatrixSpace(IntegerRing(), n, n+10)![Random(%s,%s) : i in [1..n*(n+10)]];
t := Cputime();
K := Rank(A);
s := Cputime(t);
"""%(self.n,self.min,self.max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

class Rank2_ZZ(Benchmark):
    """
    Rank 2 over ZZ:
    Given a (n + 10) x n matrix over ZZ with random entries
    between min and max, compute the rank.

    INPUT:

    - ``n`` - matrix dimension (default: ``400``)
    - ``min`` - minimal value for entries of matrix (default: ``0``)
    - ``max`` - maximal value for entries of matrix (default: ``2**64``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.rank2_ZZ(300)
        sage: tm = b.rank2_ZZ(300, system='magma')  # optional - magma
    """
    def __init__(self, n=400, min=0, max=2**64):
        self.n, self.min, self.max = n, min, max
    def sage(self):
        A = random_matrix(ZZ, self.n+10, self.n, x=self.min, y=self.max+1)
        t = cputime()
        v = A.rank()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := RMatrixSpace(IntegerRing(), n+10, n)![Random(%s,%s) : i in [1..n*(n+10)]];
t := Cputime();
K := Rank(A);
s := Cputime(t);
"""%(self.n,self.min,self.max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

# Smith Form

class Smithform_ZZ(Benchmark):
    """
    Smith Form over ZZ:
    Given a n x n matrix over ZZ with random entries
    between min and max, compute the Smith normal form.

    INPUT:

    - ``n`` - matrix dimension (default: ``128``)
    - ``min`` - minimal value for entries of matrix (default: ``0``)
    - ``max`` - maximal value for entries of matrix (default: ``9``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.smithform_ZZ(100)
        sage: tm = b.smithform_ZZ(100, system='magma')  # optional - magma
    """

    def __init__(self, n=128, min=0, max=9):
        self.n, self.min, self.max = n, min, max

    def sage(self):
        A = random_matrix(ZZ, self.n, self.n, x=self.min, y=self.max+1)
        t = cputime()
        v = A.elementary_divisors()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
t := Cputime();
K := ElementaryDivisors(A);
s := Cputime(t);
"""%(self.n,self.min,self.max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

class MatrixMultiply_ZZ(Benchmark):
    """
    Matrix multiplication over ZZ
    Given an n x n matrix A over ZZ with random entries
    between min and max, inclusive, compute A * (A+1).

    INPUT:

    - ``n`` - matrix dimension (default: ``300``)
    - ``min`` - minimal value for entries of matrix (default: ``-9``)
    - ``max`` - maximal value for entries of matrix (default: ``9``)
    - ``times`` - number of experiments (default: ``1``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_multiply_ZZ(200)
        sage: tm = b.matrix_multiply_ZZ(200, system='magma')  # optional - magma
    """

    def __init__(self, n=300, min=-9, max=9, times=1):
        self.n, self.min, self.max = n, min, max
        self.times = times

    def sage(self):
        A = random_matrix(ZZ, self.n, self.n, x=self.min, y=self.max+1)
        B = A + 1
        t = cputime()
        for z in range(self.times):
            v = A * B
        return cputime(t)/self.times

    def magma(self):
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
B := A + 1;
t := Cputime();
for z in [1..%s] do
    K := A * B;
end for;
s := Cputime(t);
"""%(self.n,self.min,self.max,self.times)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))/times

class MatrixAdd_ZZ(Benchmark):
    """
    Matrix addition over ZZ
    Given an n x n matrix A and B over ZZ with random entries between
    ``min`` and ``max``, inclusive, compute A + B ``times`` times.

    INPUT:

    - ``n`` - matrix dimension (default: ``200``)
    - ``min`` - minimal value for entries of matrix (default: ``-9``)
    - ``max`` - maximal value for entries of matrix (default: ``9``)
    - ``system`` - either 'sage' or 'magma' (default: 'sage')
    - ``times`` - number of experiments (default: ``50``)

    EXAMPLES::
    """
    def __init__(self, n=200, min=-9, max=9, times=50):
        self.n, self.min, self.max = n, min, max
        self.times = times

    def sage(self):
        A = random_matrix(ZZ, self.n, self.n, x=self.min, y=self.max+1)
        B = random_matrix(ZZ, self.n, self.n, x=self.min, y=self.max+1)
        t = cputime()
        for z in range(self.times):
            v = A + B
        return cputime(t)/self.times

    def magma(self):
        code = """
n := %s;
min := %s;
max := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(min,max) : i in [1..n^2]];
B := MatrixAlgebra(IntegerRing(), n)![Random(min,max) : i in [1..n^2]];
t := Cputime();
for z in [1..%s] do
    K := A + B;
end for;
s := Cputime(t);
"""%(self.n,self.min,self.max,self.times)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))/times

class MatrixAdd_ZZ_2(MatrixAdd_ZZ):
    """
    Matrix addition over ZZ.
    Given an n x n matrix A and B over ZZ with random ``bits``-bit
    entries, compute A + B.

    INPUT:

    - ``n`` - matrix dimension (default: ``200``)
    - ``bits`` - bitsize of entries
    - ``times`` - number of experiments (default: ``50``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_add_ZZ_2(200)
        sage: tm = b.matrix_add_ZZ_2(200, system='magma')  # optional - magma
    """

    def __init__(self, n=200, bits=16, times=50):
        b = 2**bits
        MatrixAdd_ZZ.__init__(self, n=n, min=-b, max=b, times=times)


class Det_ZZ(Benchmark):
    """
    Dense integer determinant over ZZ.
    Given an n x n matrix A over ZZ with random entries
    between min and max, inclusive, compute det(A).

    INPUT:

    - ``n`` - matrix dimension (default: ``200``)
    - ``min`` - minimal value for entries of matrix (default: ``1``)
    - ``max`` - maximal value for entries of matrix (default: ``100``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.det_ZZ(200)
        sage: tm = b.det_ZZ(200, system='magma')  # optional - magma
    """

    def __init__(self, n=200, min=1, max=100):
        self.n, self.min, self.max = n, min, max

    def sage(self):
        A = random_matrix(ZZ, self.n, self.n, x=self.min, y=self.max+1)
        t = cputime()
        d = A.determinant()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
t := Cputime();
d := Determinant(A);
s := Cputime(t);
"""%(self.n,self.min,self.max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

class Det_QQ(Benchmark):
    """
    Dense rational determinant over QQ.
    Given an n x n matrix A over QQ with random entries
    with numerator bound and denominator bound, compute det(A).

    INPUT:

    - ``n`` - matrix dimension (default: ``200``)
    - ``num_bound`` - numerator bound, inclusive (default: ``10``)
    - ``den_bound`` - denominator bound, inclusive (default: ``10``)
    - ``system`` - either 'sage' or 'magma' (default: 'sage')

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.det_QQ(200)
        sage: ts = b.det_QQ(10, num_bound=100000, den_bound=10000)
        sage: tm = b.det_QQ(200, system='magma')  # optional - magma
    """
    def __init__(self, n=300, num_bound=10, den_bound=10):
        self.n = n
        self.num_bound, self.den_bound = num_bound, den_bound

    def sage(self):
        A = random_matrix(QQ, self.n, self.n, num_bound=self.num_bound, den_bound=self.den_bound)
        t = cputime()
        d = A.determinant()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := MatrixAlgebra(RationalField(), n)![Random(%s,%s)/Random(1,%s) : i in [1..n^2]];
t := Cputime();
d := Determinant(A);
s := Cputime(t);
"""%(self.n,-self.num_bound,self.num_bound,self.den_bound)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

class Vecmat_ZZ(Benchmark):
    """
    Vector matrix multiplication over ZZ.

    Given an n x n  matrix A over ZZ with random entries
    between min and max, inclusive, and v the first row of A,
    compute the product v * A.

    INPUT:

    - ``n`` - matrix dimension (default: ``300``)
    - ``min`` - minimal value for entries of matrix (default: ``-9``)
    - ``max`` - maximal value for entries of matrix (default: ``9``)
    - ``times`` - number of runs (default: ``200``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.vecmat_ZZ(300)  # long time
        sage: tm = b.vecmat_ZZ(300, system='magma')  # optional - magma
    """

    def __init__(self, n=300, min=-9, max=9, times=200):
        self.n, self.min, self.max = n, min, max
        self.times = times

    def sage(self):
        A = random_matrix(ZZ, self.n, self.n, x=self.min, y=self.max+1)
        v = A.row(0)
        t = cputime()
        for z in range(self.times):
            w = v * A
        return cputime(t)/self.times

    def magma(self):
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
v := A[1];
t := Cputime();
for z in [1..%s] do
    K := v * A;
end for;
s := Cputime(t);
"""%(self.n,self.min,self.max,self.times)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))/times

#######################################################################
# Dense Benchmarks over GF(p), for small p.
#######################################################################

def report_GF(p=16411, **kwds):
    """
    Runs all the reports for finite field matrix operations, for
    prime p=16411.

    INPUT:

    - ``p`` - ignored
    - ``**kwds`` - passed through to :func:`report`

    .. note::

        right now, even though p is an input, it is being ignored!  If
        you need to check the performance for other primes, you can
        call individual benchmark functions.

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: print "starting"; import sys; sys.stdout.flush(); b.report_GF(systems=['sage'])
        starting...
        ======================================================================
        Dense benchmarks over GF with prime 16411
        ======================================================================
        ...
        ======================================================================
    """
    F = [Rank_GF, Rank2_GF, Nullspace_GF, Charpoly_GF,
         MatrixMultiply_GF, Det_GF]
    title = 'Dense benchmarks over GF with prime %i' % p
    report(F, title, **kwds)

# Nullspace over GF

class Nullspace_GF(Benchmark):
    """
    Given a n+1 x n  matrix over GF(p) with random
    entries, compute the nullspace.

    INPUT:

    - ``n`` - matrix dimension (default: 300)
    - ``p`` - prime number (default: ``16411``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.nullspace_GF(300)
        sage: tm = b.nullspace_GF(300, system='magma')  # optional - magma
    """

    def __init__(self, n=300, p=16411):
        self.n, self.p = n, p

    def sage(self):
        A = random_matrix(GF(self.p), self.n, self.n+1)
        t = cputime()
        v = A.kernel()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := Random(RMatrixSpace(GF(%s), n, n+1));
t := Cputime();
K := Kernel(A);
s := Cputime(t);
"""%(self.n,self.p)
        if verbose: print code
        magma.eval(code)
        return magma.eval('s')

# Characteristic Polynomial over GF

class Charpoly_GF(Benchmark):
    """
    Given a n x n matrix over GF with random entries, compute the
    charpoly.

    INPUT:

    - ``n`` - matrix dimension (default: 100)
    - ``p`` - prime number (default: ``16411``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.charpoly_GF(100)
        sage: tm = b.charpoly_GF(100, system='magma')  # optional - magma
    """

    def __init__(self, n=100, p=16411):
        self.n, self.p = n, p

    def sage(self):
        A = random_matrix(GF(self.p), self.n, self.n)
        t = cputime()
        v = A.charpoly()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
K := CharacteristicPolynomial(A);
s := Cputime(t);
"""%(self.n,self.p)
        if verbose: print code
        magma.eval(code)
        return magma.eval('s')

class MatrixAdd_GF(Benchmark):
    """
    Given two n x n matrix over GF(p) with random entries, add them.

    INPUT:

    - ``n`` - matrix dimension (default: 300)
    - ``p`` - prime number (default: ``16411``)
    - ``times`` - number of experiments (default: ``100``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_add_GF(500, p=19)
        sage: tm = b.matrix_add_GF(500, p=19, system='magma')  # optional - magma
    """

    def __init__(self, n=1000, p=16411,times=100):
        self.n, self.p = n, p
        self.times = times

    def sage(self):
        A = random_matrix(GF(self.p), self.n, self.n)
        B = random_matrix(GF(self.p), self.n, self.n)
        t = cputime()
        for self.n in range(self.times):
            v = A + B
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
B := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
for z in [1..%s] do
    K := A + B;
end for;
s := Cputime(t);
"""%(self.n,self.p,self.p,self.times)
        if verbose: print code
        magma.eval(code)
        return magma.eval('s')

# Matrix multiplication over GF(p)

class MatrixMultiply_GF(Benchmark):
    """
    Given an n x n matrix A over GF(p) with random entries, compute
    A * (A+1).

    INPUT:

    - ``n`` - matrix dimension (default: 100)
    - ``p`` - prime number (default: ``16411``)
    - ``times`` - number of experiments (default: ``3``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_multiply_GF(100, p=19)
        sage: tm = b.matrix_multiply_GF(100, p=19, system='magma')  # optional - magma
    """

    def __init__(self, n=100, p=16411, times=3):
        self.n, self.p = n, p
        self.times = times

    def sage(self):
        A = random_matrix(GF(self.p), self.n)
        B = A + 1
        t = cputime()
        for self.n in range(self.times):
            v = A * B
        return cputime(t) / self.times

    def magma(self):
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
B := A + 1;
t := Cputime();
for z in [1..%s] do
    K := A * B;
end for;
s := Cputime(t);
"""%(self.n,self.p,self.times)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))/times

class Rank_GF(Benchmark):
    """
    Rank over GF(p):
    Given a n x (n+10) matrix over GF(p) with random entries, compute the rank.

    INPUT:

    - ``n`` - matrix dimension (default: 300)
    - ``p`` - prime number (default: ``16411``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.rank_GF(1000)
        sage: tm = b.rank_GF(1000, system='magma')  # optional - magma
    """

    def __init__(self, n=500, p=16411):
        self.n, self.p = n, p

    def sage(self):
        A = random_matrix(GF(self.p), self.n, self.n+10)
        t = cputime()
        v = A.rank()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
K := Rank(A);
s := Cputime(t);
"""%(self.n,self.p)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

class Rank2_GF(Benchmark):
    """
    Rank over GF(p): Given a (n + 10) x n matrix over GF(p) with
    random entries, compute the rank.

    INPUT:

    - ``n`` - matrix dimension (default: 300)
    - ``p`` - prime number (default: ``16411``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.rank2_GF(500)
        sage: tm = b.rank2_GF(500, system='magma')  # optional - magma
    """

    def __init__(self, n=500, p=16411):
        self.n, self.p = n, p

    def sage(self):
        A = random_matrix(GF(self.p), self.n+10, self.n)
        t = cputime()
        v = A.rank()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
K := Rank(A);
s := Cputime(t);
"""%(self.n,self.p)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

class Det_GF(Benchmark):
    """
    Dense determinant over GF(p).
    Given an n x n matrix A over GF with random entries compute
    det(A).

    INPUT:

    - ``n`` - matrix dimension (default: 300)
    - ``p`` - prime number (default: ``16411``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.det_GF(1000)
        sage: tm = b.det_GF(1000, system='magma')  # optional - magma
    """

    def __init__(self, n=400, p=16411):
        self.n, self.p = n, p

    def sage(self):
        A = random_matrix(GF(self.p), self.n, self.n)
        t = cputime()
        d = A.determinant()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
d := Determinant(A);
s := Cputime(t);
"""%(self.n,self.p)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

#######################################################################
# Dense Benchmarks over QQ
#######################################################################

def hilbert_matrix(n):
    """
    Returns the Hilbert matrix of size n over rationals.

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: b.hilbert_matrix(3)
        [  1 1/2 1/3]
        [1/2 1/3 1/4]
        [1/3 1/4 1/5]
    """
    A = Matrix(QQ,n,n)
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            A[i,j] =  QQ(1)/((i+1)+(j+1)-1)
    return A

# Reduced row echelon form over QQ

class Echelon_QQ(Benchmark):
    """
    Given a n x (2*n) matrix over QQ with random integer entries
    between min and max, compute the reduced row echelon form.

    INPUT:

    - ``n`` - matrix dimension (default: ``300``)
    - ``min`` - minimal value for entries of matrix (default: ``-9``)
    - ``max`` - maximal value for entries of matrix (default: ``9``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.echelon_QQ(100)
        sage: tm = b.echelon_QQ(100, system='magma')  # optional - magma
    """

    def __init__(self, n=100, min=0, max=9):
        self.n, self.min, self.max = n, min, max

    def sage(self):
        A = random_matrix(ZZ, self.n, 2*self.n, x=self.min, y=self.max+1).change_ring(QQ)
        t = cputime()
        v = A.echelon_form()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := RMatrixSpace(RationalField(), n, 2*n)![Random(%s,%s) : i in [1..n*2*n]];
t := Cputime();
K := EchelonForm(A);
s := Cputime(t);
"""%(self.n,self.min,self.max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

# Invert a matrix over QQ.

class Inverse_QQ(Benchmark):
    """
    Given a n x n matrix over QQ with random integer entries
    between min and max, compute the reduced row echelon form.

    INPUT:

    - ``n`` - matrix dimension (default: ``300``)
    - ``min`` - minimal value for entries of matrix (default: ``-9``)
    - ``max`` - maximal value for entries of matrix (default: ``9``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.inverse_QQ(100)
        sage: tm = b.inverse_QQ(100, system='magma')  # optional - magma
    """

    def __init__(self, n=100, min=0, max=9):
        self.n, self.min, self.max = n, min, max

    def sage(self):
        A = random_matrix(ZZ, self.n, self.n, x=self.min, y=self.max+1).change_ring(QQ)
        t = cputime()
        v = ~A
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := MatrixAlgebra(RationalField(), n)![Random(%s,%s) : i in [1..n*n]];
t := Cputime();
K := A^(-1);
s := Cputime(t);
"""%(self.n,self.min,self.max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

# Matrix multiplication over QQ
class MatrixMultiply_QQ(Benchmark):
    """
    Given an n x n matrix A over QQ with random entries
    whose numerators and denominators are bounded by bnd,
    compute A * (A+1).

    INPUT:

    - ``n`` - matrix dimension (default: ``300``)
    - ``bnd`` - numerator and denominator bound (default: ``bnd``)
    - ``times`` - number of experiments (default: ``1``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_multiply_QQ(100)
        sage: tm = b.matrix_multiply_QQ(100, system='magma')  # optional - magma
    """

    def __init__(self, n=100, bnd=2, times=1):
        self.n, self.bnd, self.times = n, bnd, times

    def sage(self):
        A = random_matrix(QQ, self.n, self.n, num_bound=self.bnd, den_bound=self.bnd)
        B = A + 1
        t = cputime()
        for z in range(self.times):
            v = A * B
        return cputime(t)/self.times

    def magma(self):
        A = magma(random_matrix(QQ, n, n, num_bound=bnd, den_bound=bnd))
        code = """
n := %s;
A := %s;
B := A + 1;
t := Cputime();
for z in [1..%s] do
    K := A * B;
end for;
s := Cputime(t);
"""%(self.n, A.name(), self.times)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))/times

# Determinant of Hilbert matrix
class DetHilbert_QQ(Benchmark):
    """
    Runs the benchmark for calculating the determinant of the hilbert
    matrix over rationals of dimension n.

    INPUT:

    - ``n`` - matrix dimension (default: ``300``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.det_hilbert_QQ(50)
        sage: tm = b.det_hilbert_QQ(50, system='magma')  # optional - magma
    """

    def __init__(self, n=80):
        self.n = n

    def sage(self):
        A = hilbert_matrix(self.n)
        t = cputime()
        d = A.determinant()
        return cputime(t)

    def magma(self):
        code = """
h := HilbertMatrix(%s);
tinit := Cputime();
d := Determinant(h);
s := Cputime(tinit);
delete h;
"""%self.n
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

# inverse of Hilbert matrix
class InvertHilbert_QQ(Benchmark):
    """
    Runs the benchmark for calculating the inverse of the hilbert
    matrix over rationals of dimension n.

    INPUT:

    - ``n`` - matrix dimension (default: ``300``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.invert_hilbert_QQ(30)
        sage: tm = b.invert_hilbert_QQ(30, system='magma')  # optional - magma
    """

    def __init__(self, n=40):
        self.n = n

    def sage(self):
        A = hilbert_matrix(self.n)
        t = cputime()
        d = A**(-1)
        return cputime(t)

    def magma(self):
        code = """
h := HilbertMatrix(%s);
tinit := Cputime();
d := h^(-1);
s := Cputime(tinit);
delete h;
"""%self.n
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

class MatrixVector_QQ(Benchmark):
    """
    Compute product of square ``n`` matrix by random vector with num and
    denom bounded by ``h`` the given number of ``times``.

    INPUT:

    - ``n`` - matrix dimension (default: ``300``)
    - ``h`` - numerator and denominator bound (default: ``bnd``)
    - ``times`` - number of experiments (default: ``1``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.MatrixVector_QQ(500)
        sage: tm = b.MatrixVector_QQ(500, system='magma')  # optional - magma
    """

    def __init__(self, n=1000,h=100,times=1):
        self.n, self.h, self.times = n, h, times

    def sage(self):
        V=QQ**self.n
        v=V.random_element(self.h)
        M=random_matrix(QQ,self.n)
        t=cputime()
        for i in range(self.times):
            w=M*v
        return cputime(t)

    def magma(self):
        code = """
            n:=%s;
            h:=%s;
            times:=%s;
            v:=VectorSpace(RationalField(),n)![Random(h)/(Random(h)+1) : i in [1..n]];
            M:=MatrixAlgebra(RationalField(),n)![Random(h)/(Random(h)+1) : i in [1..n^2]];
            t := Cputime();
            for z in [1..times] do
                W:=v*M;
            end for;
            s := Cputime(t);
        """%(self.n,self.h,times)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

#######################################################################
# Dense Benchmarks over machine reals
# Note that the precision in reals for MAGMA is base 10, while in
# sage it is in base 2
#######################################################################

# Real Nullspace

class Nullspace_RR(Benchmark):
    """
    Nullspace over RR:
    Given a n+1 x n matrix over RR with random entries
    between min and max, compute the nullspace.

    INPUT:

    - ``n`` - matrix dimension (default: ``300``)
    - ``min`` - minimal value for entries of matrix (default: ``0``)
    - ``max`` - maximal value for entries of matrix (default: ``10``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.nullspace_RR(100)
        sage: tm = b.nullspace_RR(100, system='magma')  # optional - magma
    """

    def __init__(self, n=300, min=0, max=10):
        self.n, self.min, self.max = n, min, max

    def sage(self):
        from sage.rings.real_mpfr import RR
        A = random_matrix(ZZ, self.n+1, self.n, x=self.min, y=self.max+1).change_ring(RR)
        t = cputime()
        v = A.kernel()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := RMatrixSpace(RealField(16), n+1,n)![Random(%s,%s) : i in [1..n*(n+1)]];
t := Cputime();
K := Kernel(A);
s := Cputime(t);
"""%(self.n,self.min,self.max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

class Nullspace_RDF(Benchmark):
    """
    Nullspace over RDF:
    Given a n+1 x n  matrix over RDF with random entries
    between min and max, compute the nullspace.

    INPUT:

    - ``n`` - matrix dimension (default: ``300``)
    - ``min`` - minimal value for entries of matrix (default: ``0``)
    - ``max`` - maximal value for entries of matrix (default: `10``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.nullspace_RDF(100)  # long time
        sage: tm = b.nullspace_RDF(100, system='magma')  # optional - magma
    """

    def __init__(self, n=300, min=0, max=10):
        self.n, self.min, self.max = n, min, max

    def sage(self):
        from sage.rings.real_double import RDF
        A = random_matrix(ZZ, self.n+1, self.n, x=self.min, y=self.max+1).change_ring(RDF)
        t = cputime()
        v = A.kernel()
        return cputime(t)

    def magma(self):
        code = """
n := %s;
A := RMatrixSpace(RealField(16), n+1,n)![Random(%s,%s) : i in [1..n*(n+1)]];
t := Cputime();
K := Kernel(A);
s := Cputime(t);
"""%(self.n,self.min,self.max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
