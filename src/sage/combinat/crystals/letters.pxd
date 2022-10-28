from sage.structure.element cimport Element

cdef class Letter(Element):
    cdef readonly int value

cdef class EmptyLetter(Element):
    cdef readonly str value
    cdef e(self, int i)
    cdef f(self, int i)
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)

cdef class Crystal_of_letters_type_A_element(Letter):
    cdef Letter e(self, int i)
    cdef Letter f(self, int i)
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)

cdef class Crystal_of_letters_type_B_element(Letter):
    cdef Letter e(self, int i)
    cdef Letter f(self, int i)
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)

cdef class Crystal_of_letters_type_C_element(Letter):
    cdef Letter e(self, int i)
    cdef Letter f(self, int i)
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)

cdef class Crystal_of_letters_type_D_element(Letter):
    cdef Letter e(self, int i)
    cdef Letter f(self, int i)
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)

cdef class Crystal_of_letters_type_G_element(Letter):
    cdef Letter e(self, int i)
    cdef Letter f(self, int i)
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)

cdef class LetterTuple(Element):
    cdef readonly tuple value
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)

cdef class Crystal_of_letters_type_E6_element(LetterTuple):
    cdef LetterTuple e(self, int i)
    cdef LetterTuple f(self, int i)

cdef class Crystal_of_letters_type_E6_element_dual(LetterTuple):
    cdef LetterTuple lift(self)
    cdef LetterTuple retract(self, LetterTuple p)
    cdef LetterTuple e(self, int i)
    cdef LetterTuple f(self, int i)

cdef class Crystal_of_letters_type_E7_element(LetterTuple):
    cdef LetterTuple e(self, int i)
    cdef LetterTuple f(self, int i)

cdef class BKKLetter(Letter):
    cdef Letter e(self, int i)
    cdef Letter f(self, int i)

cdef class QueerLetter_element(Letter):
    cdef Letter e(self, int i)
    cdef Letter f(self, int i)
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)

cdef class LetterWrapped(Element):
    cdef readonly Element value
    cdef tuple _to_tuple(self)
    cdef LetterWrapped e(self, int i)
    cdef LetterWrapped f(self, int i)
    cdef int epsilon(self, int i)
    cdef int phi(self, int i)
