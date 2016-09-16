"""
Linbox interface
"""

## NOTE: The sig_on()/sig_off() stuff can't go in here -- it has to be in the
## code that calls these functions.  Otherwise strangely objects get left
## in an incorrect state.

from sage.misc.misc import verbose, get_verbose, cputime, UNAME

from cython.operator import dereference
from libcpp.vector cimport vector

from sage.libs.linbox.modular cimport *
from sage.libs.gmp.mpz cimport *
from sage.rings.integer cimport Integer

from sage.modules.vector_modn_sparse cimport set_entry

##################### Data Structures

#### integer ring
cdef extern from "givaro/zring.h":
    cdef cppclass GivaroIntegerRing "Givaro::ZRing<Givaro::Integer>":
        GivaroIntegerRing()

# for c++ debugging purposes
cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream() except +
        ostream& write (char* s, int n) except +
    cdef ostream cout


cdef extern from "givaro/givinteger.h" namespace "Givaro":
    #cdef cppclass mpz_class:
    #    mpz_class()
    #    mpz_class(mpz_t a)

    # LinBox integers are actualy Givaro Integers, 
    # which are themselves thin wrappers around GMPXX's mpz_class, 
    # which are themsemves thin wrappers around GMP's mpz_t...

    # the Givaro version currently bundled with sage has no constructor of Givaro::Integer from mpz_t....
    # Later Givaro version do have it. In the meanwhile, we have to use the horrible SpyInteger hack.
    cdef cppclass GivaroInteger "Givaro::Integer":
        GivaroInteger()
        # GivaroInteger(const mpz_class &a)    --- in more recent Givaro
        void output "print" (ostream &)

    cdef cppclass LinBoxSpyInteger "LinBox::SpyInteger":   # necessary evil
        mpz_t get_mpz(GivaroInteger &t)                    #    for now

#### ZZ[X] elements
ctypedef vector[GivaroInteger] LinBox_ZZX_element

#### sparse matrices mod small prime
cdef extern from "linbox/matrix/sparse-matrix.h"  namespace "LinBox":
    cdef cppclass SparseMatrix[Field]:
        const Field & field()
        SparseMatrix (const ModIntField &F, size_t m, size_t n)            # FixME Field
        void setEntry (size_t i, size_t j, const ModIntFieldElement & x)   # FixME Field

#### sparse matrices mod small prime
cdef extern from "linbox/vector/vector.h" namespace "LinBox":
    cdef cppclass DenseVector[Field]:
        DenseVector()
        DenseVector(const Field &F, size_t n)
        void setEntry (size_t i, const ModIntFieldElement & x)   # FixME Field
        ModIntFieldElement & getEntry (size_t i)                 # FixME Field

#### dense integer matrices
cdef extern from "linbox/matrix/dense-matrix.h" namespace "LinBox":
    cdef cppclass BlasMatrix[Ring]:
        BlasMatrix(const GivaroIntegerRing &F, int &nrows, int &ncols)
        void setEntry (size_t i, size_t j, const GivaroInteger & x)
        const GivaroInteger getEntry (size_t i, size_t j)
        #void write(ostream)

ctypedef BlasMatrix[GivaroIntegerRing] DenseMatrixZZ
ctypedef DenseVector[ModIntField] DenseVectorGFp
ctypedef SparseMatrix[ModIntField] SparseMatrixGFp


######### algorithms
cdef extern from "linbox/solutions/methods.h" namespace "LinBox":
    cdef cppclass Method:
        Method()
    cdef cppclass SparseEliminationTraits(Method):
        SparseEliminationTraits()
    cdef cppclass WiedemannTraits(Method):
        WiedemannTraits()


# this typedef works around the unability of cython to parse references correctly
# ---> https://mail.python.org/pipermail/cython-devel/2016-February/004673.html
ctypedef unsigned long & long_ref; 

cdef extern from "linbox/solutions/rank.h" namespace "LinBox":
    cdef long_ref rank(unsigned long &, const SparseMatrixGFp&, const Method&) 
    cdef long_ref rank(unsigned long &, const DenseMatrixZZ&) 

ctypedef vector[ModIntFieldElement] vec_mod_int

cdef extern from "linbox/solutions/solve.h" namespace "LinBox":
    cdef long_ref solve(DenseVectorGFp &, SparseMatrixGFp&, const DenseVectorGFp&, const Method&) 

cdef extern from "linbox/solutions/det.h" namespace "LinBox":
    cdef GivaroInteger & det(GivaroInteger &, DenseMatrixZZ &) 

cdef extern from "linbox/solutions/minpoly.h" namespace "LinBox":
    cdef LinBox_ZZX_element & minpoly(LinBox_ZZX_element &, DenseMatrixZZ &) 

cdef extern from "linbox/solutions/charpoly.h" namespace "LinBox":
    cdef LinBox_ZZX_element & charpoly(LinBox_ZZX_element &, DenseMatrixZZ &) 


##########################################################################
## Sparse matrices mod n
##########################################################################

cdef class Linbox_matrix_modn_sparse:
    cdef set(self, int modulus, size_t nrows, size_t ncols, c_vector_modint *rows):
        cdef ModIntField * F_ptr = new ModIntField(modulus)
        cdef SparseMatrixGFp *linbox_M = new SparseMatrixGFp(dereference(F_ptr), nrows, ncols)
        
        cdef int i, j
        for i in range(nrows):
            for j in range(rows[i].num_nonzero):
                linbox_M.setEntry(i, rows[i].positions[j], <ModIntFieldElement> rows[i].entries[j])
        
        self.nrows = nrows
        self.ncols = ncols
        self._M = linbox_M
    
    def __deallocate__(self):
        cdef SparseMatrixGFp *M = <SparseMatrixGFp *> self._M
        del M

    cdef size_t rank(self, SparseAlgorithm algorithm):
        cdef unsigned long r
        cdef SparseMatrixGFp *M = <SparseMatrixGFp *> self._M

        rank(r, dereference(M), SparseEliminationTraits())
        return r

    cdef void solve(self, c_vector_modint **x, c_vector_modint *b, SparseAlgorithm algorithm):
        cdef SparseMatrixGFp *M = <SparseMatrixGFp *> self._M
        cdef SparseEliminationTraits method = SparseEliminationTraits()
        
        cdef ModIntField F = M.field()
        cdef DenseVectorGFp B = DenseVectorGFp(F, self.ncols)
        cdef DenseVectorGFp X = DenseVectorGFp(F, self.nrows)

        #B = DenseVectorGFp(F, self.ncols)
        #X = DenseVectorGFp(F, self.nrows)

        # scatter b into a dense vector
        cdef size_t i
        for i in range(b.num_nonzero):
            B.setEntry(b.positions[i], <ModIntFieldElement> b.entries[i])
        solve(X, dereference(M), B, method)
        
        # linbox returns a dense solution, and we have to make it sparse somehow...
        for i in range(self.ncols):
            set_entry(dereference(x), i, <int> X.getEntry(i))

##########################################################################
## Sparse matrices over ZZ
##########################################################################


##########################################################################
## Dense matrices over ZZ
##########################################################################

######### interface ##########

cdef class Linbox_integer_dense:
    def __deallocate__(self):
        cdef SparseMatrixGFp *M = <SparseMatrixGFp *> self._M
        if M:
            del M

    cdef set(self, size_t nrows, size_t ncols, mpz_t ** matrix):
        cdef size_t i, j
        cdef GivaroIntegerRing *ZZ = new GivaroIntegerRing()
        cdef LinBoxSpyInteger spy   # this crap just have static methods
        cdef GivaroInteger t, u
        cdef DenseMatrixZZ *M = new DenseMatrixZZ(dereference(ZZ), nrows, ncols)    

        self.nrows = nrows
        self.ncols = ncols
        # copy the input matrix into a (dense) LinBox matrix. 
        # IMPROVE-ME: the copying is not necessary, because Linbox
        # often copies it again. However, we would need to deal with
        # the various wrappers around mpz_t.
        for i in range(nrows):
            for j in range(ncols):
                mpz_set(spy.get_mpz(t), matrix[i][j])
                M.setEntry(i, j, t)

        self._M = M

    cdef size_t rank(self):
        cdef size_t r
        cdef DenseMatrixZZ *M = <DenseMatrixZZ *> self._M
        rank(r, dereference(M))
        return r

    cdef det(self):
        cdef DenseMatrixZZ *M = <DenseMatrixZZ *> self._M
        cdef GivaroInteger d
        cdef Integer z  # Sage Integer
        cdef LinBoxSpyInteger spy   # this crap just have static methods

        det(d, dereference(M)) # we get a LinBox integer, and we need to convert it to a Sage integer..
        z = Integer()
        mpz_set(z.value, spy.get_mpz(d));
        return z

    cdef minpoly(self):
        cdef DenseMatrixZZ *M = <DenseMatrixZZ *> self._M
        cdef LinBox_ZZX_element pM
        cdef LinBoxSpyInteger spy   # this crap just have static methods

        minpoly(pM, dereference(M))
        # make a list of sage integers from the vector of LinBox integers
        l = []
        for i in range(pM.size()):
            k = Integer()
            mpz_set(k.value, spy.get_mpz(pM[i]))
            l.append(k)
        return l
    
    cdef charpoly(self):
        cdef DenseMatrixZZ *M = <DenseMatrixZZ *> self._M
        cdef LinBox_ZZX_element pM = LinBox_ZZX_element()
        cdef LinBoxSpyInteger spy   # this crap just have static methods

        # print("going in")
        charpoly(pM, dereference(M))
        # print("gone out")
        # make a list of sage integers from the vector of LinBox integers
        l = []
        for i in range(pM.size()):
            k = Integer()
            mpz_set(k.value, spy.get_mpz(pM[i]))
            l.append(k)
        return l
