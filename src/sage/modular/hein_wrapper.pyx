# Sage wrapper for Jeffery Hain's code to compute spaces of newforms by Birch's method
# Copyright (C) 2017-8 by Kiran S. Kedlaya <kskedl@gmail.com>
#
#encoding=utf8
#distutils: language = c++
#distutils: libraries = ['gmpxx', 'gmp', 'm', 'pthread']
#distutils: sources = ['AutomorphismZZ.cpp', 'CharacterZZ.cpp', 'Eigenvector.cpp', 'GenusZZ.cpp', 'IsometryZZ.cpp', 'MathZZ.cpp', 'NeighborIteratorZZ.cpp', 'QuadFormZZ.cpp', 'SparseMatrix.cpp']
#distutils: include_dirs = ['/usr/local/include', '.']
#distutils: extra_compile_args = ['-g', '-Wall', '-O3', '-Wextra', '-Werror', '-pedantic', '-std=c++11']

import sys
reload(sys)
sys.setdefaultencoding('utf8')

from cython.operator cimport dereference as deref

from libcpp.map cimport map
from libcpp.unordered_map cimport unordered_map
from libcpp.memory cimport shared_ptr, make_shared

from sage.rings.integer cimport Integer
from sage.libs.gmp.types cimport mpz_t

from sage.quadratic_forms.ternary_qf import find_a_ternary_qf_by_level_disc
from sage.arith.misc import squarefree_divisors
from sage.matrix.constructor import matrix

cdef extern from "gmpxx.h":
    cdef cppclass mpz_class:
        mpz_class(mpz_t a)

    cdef cppclass mpq_class:
        pass

cdef extern from "QuadForm.h":
    cdef cppclass QuadForm[R,F]:
        QuadForm(const R& a, const R& b, const R& c,
            const R& f, const R& g, const R& h, long reduced=false)

cdef extern from "SparseMatrix.h":
    cdef cppclass SparseMatrix:
        ctypedef map[long, map[long, long]] dataMap
        const dataMap& data() const;
        long num_rows() const;

cdef extern from "Character.h":
    cdef cppclass Character[R,F]:
        Character(const R& cond)

cdef extern from "Genus.h":
    cdef cppclass Genus[R,F]:
        ctypedef shared_ptr[QuadForm[R,F]] QuadFormPtr
        ctypedef shared_ptr[HeckeOperator[R,F]] HeckePtr

        Genus(const QuadForm[R,F]& q, long numThreads=0)
        void add_character(const Character[R,F]& chi)
        void compute_hecke_operators(const R& p, long numThreads)
        HeckePtr hecke_operator(const R& p, const R& cond)

cdef extern from "HeckeOperator.h":
    cdef cppclass HeckeOperator[R,F]:
        SparseMatrix& matrix() const;

ctypedef Genus[mpz_class, mpq_class] GenusZZ
ctypedef QuadForm[mpz_class, mpq_class] QuadFormZZ
ctypedef Character[mpz_class, mpq_class] CharacterZZ

ctypedef shared_ptr[QuadFormZZ] QuadFormZZPtr
ctypedef shared_ptr[GenusZZ] GenusZZPtr

cdef GenusZZPtr create_genus_from_coeffs(Integer a, Integer b,
                                          Integer c, Integer f,
                                          Integer g, Integer h):
    cdef QuadFormZZPtr q
    cdef GenusZZPtr genus
    q = make_shared[QuadFormZZ](mpz_class(a.value), mpz_class(b.value),
                                mpz_class(c.value), mpz_class(f.value),
                                mpz_class(g.value), mpz_class(h.value), False)
    genus = make_shared[GenusZZ](deref(q))
    return(genus)

cpdef hecke_birch(N, l):
    r"""
    Compute Hecke operators at level N for each prime in l via Birch's method.

    INPUT::

    - ``N`` -- a squarefree positive integer

    - ``l`` -- a list of prime positive integers not dividing ``N``

    OUTPUT::

    A dict with one entry for each pair (d, p) where d is a squarefree
    divisor of N and p is an entry of l. This entry is a sparse integer matrix
    which, when N has an odd number of prime factors, represents the action of the 
    Hecke operator T_p on the subspace of the space M_2(Gamma_0(N), QQ)^{new} on which 
    the Atkin-Lehner involutions act via the Kronecker character of level d. (Note that
    the Eisenstein series is included for d=1.)

    When N has an even number of prime factors, one gets a similar answer except that
    the space includes forms that are old at the smallest prime factor of N.

    EXAMPLES::

    In this example, there are no forms that are old at 7 because X_0(7) has genus 0::

        sage: h = hecke_birch(91, [2, 3])
        sage: h[(13, 3)]
        [-2]

    Note the presence of the Eisenstein series for the trivial character::

        sage: h[(1,2)].charpoly().factor()
        (x - 3) * (x^2 - 2)
        sage: h[(1,3)].charpoly().factor()
        (x - 4) * (x^2 - 2)

    The following example demonstrates the presence of oldforms::

        sage: hecke_birch(11, [19])[1,19].charpoly().factor()
        (x - 20) * x
        sage: hecke_birch(11*13, [19])[1,19].charpoly().factor()
        (x - 20) * x * (x^4 - 8*x^3 - 25*x^2 + 154*x + 387)
        sage: V = ModularSymbols(11*13, sign=-1).new_subspace()
        sage: V.hecke_polynomial(19).factor()   # No factor of x
        (x - 2) * (x^4 - 8*x^3 - 25*x^2 + 154*x + 387) * (x^6 + 10*x^5 + 3*x^4 - 196*x^3 - 561*x^2 - 454*x - 104)

    TESTS::

    Compare with Brandt symbols::

        sage: B = BrandtModule(151)
        sage: h = hecke_birch(151, [2])
        sage: P = prod(h[i].charpoly() for i in h)
        sage: P == B.hecke_polynomial(2)
        True
    """
    cdef Integer d, p
    cdef SparseMatrix mat
    cdef SparseMatrix.dataMap dat

    if not N.is_squarefree():
        raise ValueError("level " + str(N) +
                         " is not squarefree")
    for p in l:
        if not p.is_prime:
            raise ValueError("index " + str(p) + " is not prime")
        if N % p == 0:
            raise ValueError("prime " + str(p) +
                             " divides level " + str(N))
    ans = {}
    Q = find_a_ternary_qf_by_level_disc(4 * N, N)
    (a, b, c, f, g, h) = Q.coefficients()
    cdef GenusZZPtr genus = create_genus_from_coeffs(a, b, c, f, g, h)
    for d in squarefree_divisors(N):
        deref(genus).add_character(CharacterZZ(mpz_class(d.value)))

    for p in l:
        deref(genus).compute_hecke_operators(mpz_class(p.value), 0)
        for d in squarefree_divisors(N):
            hecke = deref(genus).hecke_operator(mpz_class(p.value),
                                                      mpz_class(d.value))
            mat = deref(hecke).matrix()
            dat = mat.data()
            datdict = {}
            for (i, dict1) in dat:
                for j in dict1:
                    datdict[(i, j)] = Integer(dat[i][j])
            n = mat.num_rows()
            ans[(d, p)] = matrix(n, n, datdict)
    return ans
