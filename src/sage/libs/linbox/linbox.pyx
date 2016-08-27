"""
Linbox interface
"""

## NOTE: The sig_on()/sig_off() stuff can't go in here -- it has to be in the
## code that calls these functions.  Otherwise strangely objects get left
## in an incorrect state.

from sage.libs.gmp.mpz cimport *
from sage.rings.integer cimport Integer
from sage.misc.misc import verbose, get_verbose, cputime, UNAME

from cython.operator import dereference
from libcpp.vector cimport vector

##################### Data Structures


#### sparse matrices mod small prime

# cdef extern from "linbox/matrix/sparse-matrix.h"  namespace "LinBox": in v1.4.2
cdef extern from "linbox/blackbox/sparse.h" namespace "LinBox":  # v1.3.2
    cdef cppclass SparseMatrix[Field,Row]:
        SparseMatrix (const ModIntField &F, size_t m, size_t n)            # FixME Field
        void setEntry (size_t i, size_t j, const ModIntFieldElement & x)   # FixME Field

ctypedef pair[size_t, ModIntFieldElement] sparse_vector_pair
ctypedef vector[sparse_vector_pair] SparseSeqVectorGFp
ctypedef SparseMatrix[ModIntField, SparseSeqVectorGFp] SparseMatrixGFp

#### integer ring
cdef extern from "linbox/field/PID-integer.h" namespace "LinBox":
    cdef cppclass LinBoxIntegerRing "LinBox::PID_integer":
        LinBoxIntegerRing()

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream() except +
        ostream& write (char* s, int n) except +
    cdef ostream cout

cdef extern from "linbox/integer.h":
    #cdef cppclass mpz_class:
    #    mpz_class()
    #    mpz_class(mpz_t a)

    # LinBox integers are actualy Givaro Integers, 
    # which are themselves thin wrappers around GMPXX's mpz_class, 
    # which are themsemves thin wrappers around GMP's mpz_t...
    #
    # the Givaro version currently bundled with sage has no constructor of Givaro::Integer from mpz_t....
    # Later Givaro version do have it. In the meanwhile, we have to use the horrible SpyInteger hack.
    cdef cppclass LinBoxInteger "LinBox::Integer":
        LinBoxInteger()
        # LinBoxInteger(const mpz_class &a)    --- in more recent Givaro
        void output "print" (ostream &)

    cdef cppclass LinBoxSpyInteger "LinBox::SpyInteger":   # necessary evil
        mpz_t get_mpz(LinBoxInteger &t)                    #    for now

#### ZZ[X] elements
ctypedef vector[LinBoxInteger] LinBox_ZZX_element

#### dense integer matrices
cdef extern from "linbox/matrix/blas-matrix.h" namespace "LinBox":
    cdef cppclass BlasMatrix[Ring]:
        BlasMatrix(const LinBoxIntegerRing &F, int &nrows, int &ncols)
        void setEntry (size_t i, size_t j, const LinBoxInteger & x)
        const LinBoxInteger getEntry (size_t i, size_t j)
        #void write(ostream)

ctypedef BlasMatrix[LinBoxIntegerRing] DenseMatrixZZ


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
    cdef long_ref solve(vec_mod_int &, SparseMatrixGFp&, const vec_mod_int&, const Method&) 

cdef extern from "linbox/solutions/det.h" namespace "LinBox":
    cdef LinBoxInteger & det(LinBoxInteger &, DenseMatrixZZ &) 

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
        cdef vec_mod_int B = vec_mod_int(self.ncols)
        cdef vec_mod_int X = vec_mod_int(self.nrows)

        # scatter b into a dense vector
        cdef size_t i
        for i in range(b.num_nonzero):
            B[b.positions[i]] = <ModIntFieldElement> b.entries[i]
        solve(X, dereference(M), B, method)
        
        # linbox returns a dense solution, and we have to make it sparse somehow...
        for i in range(self.ncols):
            set_entry(dereference(x), i, <int> X[i])

##########################################################################
## Sparse matrices over ZZ
##########################################################################


##########################################################################
## Dense matrices over ZZ
##########################################################################

cdef extern from "linbox/linbox-sage.h":
    void linbox_integer_dense_minpoly(mpz_t* &minpoly, size_t &degree, size_t n, mpz_t** matrix)

    void linbox_integer_dense_charpoly(mpz_t* &charpoly, size_t &degree, size_t n, mpz_t** matrix)

    void linbox_integer_dense_delete_array(mpz_t* f)

    int linbox_integer_dense_matrix_matrix_multiply(mpz_t** ans, mpz_t **A, mpz_t **B, size_t A_nr, size_t A_nc, size_t B_nc)

    unsigned long linbox_integer_dense_rank(mpz_t** matrix, size_t nrows,
                                            size_t ncols)

    void linbox_integer_dense_det(mpz_t ans, mpz_t** matrix,
                             size_t nrows, size_t ncols)

    void linbox_integer_dense_smithform(mpz_t **v,
                                        mpz_t **matrix,
                                        size_t nrows, size_t ncols)

######### interface ##########

cdef class Linbox_integer_dense:
    def __deallocate__(self):
        cdef SparseMatrixGFp *M = <SparseMatrixGFp *> self._M
        if M:
            del M

    cdef set(self, size_t nrows, size_t ncols, mpz_t ** matrix):
        cdef size_t i, j
        cdef LinBoxIntegerRing *ZZ = new LinBoxIntegerRing()
        cdef LinBoxSpyInteger spy   # this crap just have static methods
        cdef LinBoxInteger t, u
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
        cdef LinBoxInteger d
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
