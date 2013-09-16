include "sage/ext/stdsage.pxi"

cdef class MMXpadic_naive:
     def __cinit__(self):
         self.x = new mmx_padic_naive()

     def __init__(self, int val, int p):
         cdef integer c_val = integer(val)
         self.p = p
         mmx_set_modulus(p)
         self.x[0] = itsc_naive(c_val)
#         self.x = new mmx_padic_naive(itsc(c_val))

     def __dealloc__(self):
         del self.x

     def __add__(MMXpadic_naive self, other):
         cdef MMXpadic_naive res = MMXpadic_naive.__new__(MMXpadic_naive)
         if isinstance(other, MMXpadic_naive):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic_naive> self.x[0] + <mmx_padic_naive> (<MMXpadic_naive>other).x[0]
         return res

     def __sub__(MMXpadic_naive self, other):
         cdef MMXpadic_naive res = MMXpadic_naive.__new__(MMXpadic_naive)
         if isinstance(other, MMXpadic_naive):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic_naive> self.x[0] - <mmx_padic_naive> (<MMXpadic_naive>other).x[0]
         return res

     def __mul__(MMXpadic_naive self, other):
         cdef MMXpadic_naive res = MMXpadic_naive.__new__(MMXpadic_naive)
         if isinstance(other, MMXpadic_naive):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic_naive> self.x[0] * <mmx_padic_naive> (<MMXpadic_naive>other).x[0]
         return res

     def __div__(MMXpadic_naive self, other):
         cdef MMXpadic_naive res = MMXpadic_naive.__new__(MMXpadic_naive)
         if isinstance(other, MMXpadic_naive):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic_naive> self.x[0] / <mmx_padic_naive> (<MMXpadic_naive>other).x[0]
         return res

     def __repr__(self):
         return <bytes> as_charp(as_mmx(as_generic(flatten_naive(self.x[0]))))

     def my_print(MMXpadic_naive self):
         mmout << flatten_naive(self.x[0])
         print

cdef class MMXpadic_relaxed:
     def __cinit__(self):
         self.x = new mmx_padic_relaxed()

     def __init__(self, int val, int p):
         cdef integer c_val = integer(val)
         self.p = p
         mmx_set_modulus(p)
         self.x[0] = itsc_relaxed(c_val)
#         self.x = new mmx_padic_relaxed(itsc(c_val))

     def __dealloc__(self):
         del self.x

     def __add__(MMXpadic_relaxed self, other):
         cdef MMXpadic_relaxed res = MMXpadic_relaxed.__new__(MMXpadic_relaxed)
         if isinstance(other, MMXpadic_relaxed):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic_relaxed> self.x[0] + <mmx_padic_relaxed> (<MMXpadic_relaxed>other).x[0]
         return res

     def __sub__(MMXpadic_relaxed self, other):
         cdef MMXpadic_relaxed res = MMXpadic_relaxed.__new__(MMXpadic_relaxed)
         if isinstance(other, MMXpadic_relaxed):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic_relaxed> self.x[0] - <mmx_padic_relaxed> (<MMXpadic_relaxed>other).x[0]
         return res

     def __mul__(MMXpadic_relaxed self, other):
         cdef MMXpadic_relaxed res = MMXpadic_relaxed.__new__(MMXpadic_relaxed)
         if isinstance(other, MMXpadic_relaxed):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic_relaxed> self.x[0] * <mmx_padic_relaxed> (<MMXpadic_relaxed>other).x[0]
         return res

     def __div__(MMXpadic_relaxed self, other):
         cdef MMXpadic_relaxed res = MMXpadic_relaxed.__new__(MMXpadic_relaxed)
         if isinstance(other, MMXpadic_relaxed):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic_relaxed> self.x[0] / <mmx_padic_relaxed> (<MMXpadic_relaxed>other).x[0]
         return res

     def __repr__(self):
         return <bytes> as_charp(as_mmx(as_generic(flatten_relaxed(self.x[0]))))

     def my_print(MMXpadic_relaxed self):
         mmout << flatten_relaxed(self.x[0])
         print
