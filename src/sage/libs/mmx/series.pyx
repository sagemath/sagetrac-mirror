include "sage/ext/stdsage.pxi"

cdef class MMXpadic:
     def __cinit__(self):
         self.x = new mmx_padic()

     def __init__(self, int val, int p):
         cdef integer c_val = integer(val)
         self.p = p
         mmx_set_modulus(p)
         self.x[0] = itsc(c_val)
#         self.x = create_padic(itsc(c_val))

     def __dealloc__(self):
         del self.x

     def __add__(MMXpadic self, other):
         cdef MMXpadic res = MMXpadic.__new__(MMXpadic)
         if isinstance(other, MMXpadic):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic> self.x[0] + <mmx_padic> (<MMXpadic>other).x[0]
         return res

     def __sub__(MMXpadic self, other):
         cdef MMXpadic res = MMXpadic.__new__(MMXpadic)
         if isinstance(other, MMXpadic):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic> self.x[0] - <mmx_padic> (<MMXpadic>other).x[0]
         return res

     def __mul__(MMXpadic self, other):
         cdef MMXpadic res = MMXpadic.__new__(MMXpadic)
         if isinstance(other, MMXpadic):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic> self.x[0] * <mmx_padic> (<MMXpadic>other).x[0]
         return res

     def __div__(MMXpadic self, other):
         cdef MMXpadic res = MMXpadic.__new__(MMXpadic)
         if isinstance(other, MMXpadic):
             mmx_set_modulus(self.p)
             res.p = self.p
             res.x[0] = <mmx_padic> self.x[0] / <mmx_padic> (<MMXpadic>other).x[0]
         return res

     def __repr__(self):
         return <bytes> as_charp(as_mmx(as_generic(flatten(self.x[0]))))

     def my_print(MMXpadic self):
         mmout << flatten(self.x[0])
         print
