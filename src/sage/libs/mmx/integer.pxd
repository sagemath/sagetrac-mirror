from common cimport *

cdef extern from "numerix/integer.hpp" namespace "mmx":
     cdef cppclass integer:
          integer* integer()
          integer* integer(int)
          integer* integer(integer&)

     integer operator + (const integer& x1, const integer& x2)
     integer operator + (const integer& x1, const unsigned long int& x2)
     integer operator - (const integer& x1, const integer& x2)
     integer operator - (const integer& x1, const unsigned long int& x2)
     integer operator * (const integer& x1, const integer& x2)
     integer operator * (const integer& x1, const unsigned long int& x2)
     integer operator / (const integer& x1, const integer& x2)
     integer operator / (const integer& x1, const unsigned long int& x2)

     void add (integer& r, const integer& x1,const integer& x2)
     void add (integer& r, const integer& x1,unsigned long int x2)

     port operator << (const port& out, const integer& x)

cdef extern from "numerix/integer.hpp":
     string as_string (const integer& x)

cdef class MMXint:
     cdef integer* x
