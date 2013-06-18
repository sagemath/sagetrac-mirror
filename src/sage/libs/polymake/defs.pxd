# distutils: language = c++

from libcpp.string cimport string
from sage.libs.gmp.types cimport mpz_t, mpq_t

cdef extern from "wrap.h" namespace "polymake":
    pass

cdef extern from "polymake/Main.h" namespace "polymake":
    cdef cppclass Main:
        Main(char*)
        void set_application(char*)
        void set_preference(char*)
        void set_custom(char*, int) except +ValueError

cdef extern from "polymake/Main.h":
    cdef cppclass PropertyValue:
        pass


cdef extern from "polymake/Rational.h" namespace 'polymake':
    cdef cppclass Integer:
        mpz_t get_rep()
        Py_ssize_t strsize(int)
        int compare(int)

    cdef cppclass Rational:
        Rational()
        Rational(mpq_t)
        mpq_t get_rep()
        Rational set(mpq_t)

cdef extern from "polymake/client.h":
    cdef extern PerlObject load "perl::Object::load" (char*) except +
    cdef cppclass PerlObject "perl::Object":
        PerlObject()
        PerlObject(char*) except +ValueError
        void VoidCallPolymakeMethod(char*) except +ValueError
        void save(char*)
        PropertyValue take(char*)
        PropertyValue give(char*) # do not add except here, see pm_get for why

    PerlObject CallPolymakeFunction (char*) except +ValueError
    PerlObject CallPolymakeFunction1 "CallPolymakeFunction" \
            (char*, int) except +ValueError
    PerlObject CallPolymakeFunction2 "CallPolymakeFunction" \
            (char*, int, int) except +ValueError
    PerlObject CallPolymakeFunction3 "CallPolymakeFunction" \
            (char*, int, int, int) except +ValueError
    PerlObject* new_PerlObject_from_PerlObject "new perl::Object" (PerlObject)
    PerlObject CallPolymakeFunction_PerlObject2 "CallPolymakeFunction" (char*, PerlObject, PerlObject) except +ValueError
    bint BoolCallPolymakeFunction_PerlObject2 "CallPolymakeFunction" (char*, PerlObject, PerlObject) except +ValueError

cdef extern from "polymake/SparseMatrix.h" namespace "polymake":
    pass

cdef extern from "polymake/Matrix.h" namespace 'polymake':
    cdef cppclass MatrixRational "Matrix<Rational>":
        MatrixRational()
        MatrixRational(int nr, int nc)
        void assign(int r, int c, Rational val)
        MatrixRational operator|(MatrixRational)
        Py_ssize_t rows()
        Py_ssize_t cols()

    Rational get_element "WRAP_CALL"(MatrixRational, int i, int j)

#cdef extern from "polymake/GenericVector.h" namespace 'polymake':
    MatrixRational ones_vector_Rational "ones_vector<Rational>" ()

#cdef extern from "polymake/GenericMatrix.h" namespace 'polymake':
    MatrixRational unit_matrix_Rational "unit_matrix<Rational>" (int dim)


    void pm_assign "WRAP_OUT" (PropertyValue, MatrixRational)

    # the except clause below is fake
    # it is used to catch errors in PerlObject.give(), however adding
    # the except statement to the declaration of give() makes cython
    # split lines like
    #        pm_get(self.pm_obj.give(prop), pm_res)
    # and store the result of give() first. This causes problems since
    # PropertyValue doesn't have a default constructor.
    void pm_get "WRAP_IN" (PropertyValue, Integer) except +ValueError
    void pm_get_String "WRAP_IN" (PropertyValue, string) except +ValueError
    void pm_get_Rational "WRAP_IN" (PropertyValue, Rational) except +ValueError
    void pm_get_MatrixRational "WRAP_IN" (PropertyValue, MatrixRational) \
            except +ValueError
    void pm_get_VectorInteger "WRAP_IN" (PropertyValue, VectorInteger) \
            except +ValueError
    void pm_get_PerlObject "WRAP_IN" (PropertyValue, PerlObject) \
            except +ValueError

cdef extern from "polymake/Vector.h" namespace 'polymake':
    cdef cppclass VectorInteger "Vector<Integer>":
        VectorInteger()
        VectorInteger(int nr)
        #void assign(int r, int val)
        Integer get "operator[]" (int i)
        int size()
