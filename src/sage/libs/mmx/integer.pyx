cdef class MMXint:
     def __cinit__(self):
         self.x = new integer()

     def __init__(self, other):
         if isinstance(other, MMXint):
             self.x[0] = <integer> self.x[0]
         else:
             self.x[0] = integer(<int> other)

     def __dealloc__(self):
         del self.x

     def my_print(self):
         cdef char* s
         s = as_charp(as_string(self.x[0]))
         printf("%s\n", s)
         print(s)
         free_charp(s)
         mmout << self.x[0]

     def __repr__(self):
         return <bytes> as_charp(as_string(self.x[0]))

     def __add__(MMXint self, other):
         cdef MMXint res = MMXint.__new__(MMXint)
         if isinstance(other, MMXint):
             res.x[0] = <integer> self.x[0] + <integer> (<MMXint>other).x[0]
         else:
             res.x[0] = <integer> self.x[0] + <int> other
         return res

     def __sub__(MMXint self, other):
         cdef MMXint res = MMXint.__new__(MMXint)
         if isinstance(other, MMXint):
             res.x[0] = <integer> self.x[0] - <integer> (<MMXint>other).x[0]
         else:
             res.x[0] = <integer> self.x[0] - <int> other
         return res

     def __mul__(MMXint self, other):
         cdef MMXint res = MMXint.__new__(MMXint)
         if isinstance(other, MMXint):
             res.x[0] = <integer> self.x[0] * <integer> (<MMXint>other).x[0]
         else:
             res.x[0] = <integer> self.x[0] * <int> other
         return res

     def __div__(MMXint self, other):
         cdef MMXint res = MMXint.__new__(MMXint)
         if isinstance(other, MMXint):
             res.x[0] = <integer> self.x[0] / <integer> (<MMXint>other).x[0]
         else:
             res.x[0] = <integer> self.x[0] / <int> other
         return res
