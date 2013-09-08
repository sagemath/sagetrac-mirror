cdef extern from "<cstdio>":
     int printf (char * format, ... )

cdef extern from "<iostream>" namespace "std":
     cdef cppclass ostream:
          pass
     ostream cout
     ostream& operator<< (ostream& os, const char* s)

cdef extern from "basix/port.hpp" namespace "mmx":
     cdef cppclass port:
          pass
     cdef port mmout
     port operator << (const port& out, char* x)

cdef extern from "basix/string.hpp" namespace "mmx":
     cdef cppclass string:
          pass

     char* as_charp (const string& s)
     void free_charp (char* s)

cdef extern from "basix/generic.hpp" namespace "mmx":
     cdef cppclass generic:
          pass

cdef extern from "basix/mmx_syntax.hpp" namespace "mmx":
     string as_mmx (const generic& g)
     string flatten_as_mmx (const generic& g)

cdef extern from "basix/syntactic.hpp" namespace "mmx":
     cdef cppclass syntactic:
          pass

     port operator << (const port& out, const syntactic& x)
     string as_string (const syntactic& g)
     generic as_generic (const syntactic& g)
