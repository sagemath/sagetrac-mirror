
from sage.structure.element cimport CommutativeAlgebraElement
from sage.structure.element cimport Element


cdef class AlgebraFMElement(CommutativeAlgebraElement):
    r"""
    Generic class for elements lying in ring extensions

    AUTHOR:

    - Xavier Caruso (2016)
    """
    def __init__(self, parent, x, *args, **kwds):
        from sage.rings.algebra_from_morphism import AlgebraFromMorphism
        if not isinstance(parent, AlgebraFromMorphism):
            raise TypeError("%s is not an instance of AlgebraFromMorphism" % parent)
        if isinstance(x, AlgebraFMElement):
            x = x._backend()
        try:
            parentx = x.parent()
            if parent.base().has_coerce_map_from(parentx):
                x = parent.base().coerce_map_from(parentx)(x)
                x = parent.defining_morphism()(x)
        except AttributeError:
            pass
        Element.__init__(self, parent)
        ring = parent._backend()
        self._element = ring(x, *args, **kwds)

    def _repr_(self):
        r"""
        Return a string representation of this element

        By default, it is the same as the strict representation
        of this element viewed as an element of the underlying
        ring.

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = L.gen()
            sage: E(x)._repr_()
            'z4'
            sage: x._repr_()
            'z4'
        """
        return str(self._element)

    def _latex_(self):
        r"""
        Return a latex representation of this element

        By default, it is the same as the latex representation
        of this element viewed as an element of the underlying
        ring.

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = L.gen()
            sage: E(x)._latex_()
            'z_{4}'
            sage: x._latex_()
            'z_{4}'
        """
        from sage.misc.latex import latex
        return str(latex(self._element))

    def _backend(self):
        return self._element

    cpdef _richcmp_(left, right, int op):
        if isinstance(right, AlgebraFMElement):
            right = right._backend()
        return left._element._richcmp_(right, op)

    cpdef _add_(self,other):
        r"""
        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = E.random_element()
            sage: y = E.random_element()
            sage: (x+y).parent() is E
            True
        """
        return self.__class__(self._parent, self._element + other._backend())

    cpdef _neg_(self):
        return self.__class__(self._parent, -self._element)

    cpdef _sub_(self,other):
        r"""
        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = E.random_element()
            sage: y = E.random_element()
            sage: (x-y).parent() is E
            True
        """
        return self.__class__(self._parent, self._element - other._backend())

    cpdef _mul_(self,other):
        r"""
        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = E.random_element()
            sage: y = E.random_element()
            sage: (x*y).parent() is E
            True
        """
        return self.__class__(self._parent, self._element * other._backend())

    #cpdef _div_(self,other)
    #cpdef _floordiv_(self,other)

    def additive_order(self):
        r"""
        Return the additive order of this element

        EXAMPLES::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = E.gen()
            sage: x.additive_order()
            5
        """
        return self._element.additive_order()

    def multiplicative_order(self):
        r"""
        Return the multiplicite order of this element

        EXAMPLES::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)

            sage: x = E.gen()
            sage: x.multiplicative_order()
            624
        """
        return self._element.multiplicative_order()

    def is_unit(self):
        r"""
        Return ``True`` if this element is a unit in the ring
        of the parent extension, ``False`` otherwise

        EXAMPLES::

            sage: A.<x> = PolynomialRing(QQ)
            sage: E = RingExtension(A,QQ)
            sage: E(4).is_unit()
            True
            sage: E(x).is_unit()
            False
        """
        return self._element.is_unit()

    def is_nilpotent(self):
        r"""
        Return ``True`` if this element is nilpotent in the ring
        of the parent extension, ``False`` otherwise

        EXAMPLES::

            sage: A.<x> = PolynomialRing(QQ)
            sage: E = RingExtension(A,QQ)
            sage: E(0).is_nilpotent()
            True
            sage: E(x).is_nilpotent()
            False
        """
        return self._element.is_nilpotent()

    def is_prime(self):
        r"""
        Return ``True`` if this element is a prime element 
        in the ring of the parent extension, ``False`` otherwise

        EXAMPLES::

            sage: A.<x> = PolynomialRing(QQ)
            sage: E = RingExtension(A,QQ)
            sage: E(x^2+1).is_prime()
            True
            sage: E(x^2-1).is_prime()
            False
        """
        return self._element.is_prime()


cdef class RingExtensionWithBasisElement(AlgebraFMElement):
    def _repr_(self):
        parent = self._parent
        names = parent._names
        _, _, j = parent.vector_space()
        coeffs = j(self)
        s = ""
        for i in range(len(names)):
            if coeffs[i].is_zero(): continue
            c = coeffs[i]
            sign = 1
            if (-c)._is_atomic():
                c = -c
                sign = -sign
            if s == "":
                if sign == -1: s = "-"
            else:
                s += " + " if sign == 1 else " - "
            ss = ""
            if c != 1:
                if c._is_atomic():
                    ss += "%s" % c
                else:
                    ss += "(%s)" % c
                if names[i] != "": ss += "*"
            ss += names[i]
            if ss == "": ss += "1"
            s += ss
        if s == "": return "0"
        return s

    def matrix(self, base=None):
        from sage.matrix.matrix_space import MatrixSpace
        parent = self._parent
        if base is None:
            base = parent._base
        _, _, j = parent.vector_space(base)
        x = self._backend()
        M = [ j(x * b._backend()) for b in parent.basis(base) ]
        return MatrixSpace(base, len(M))(M)

    def trace(self, base=None):
        if base is None or base is self._parent._base:
            return self.matrix().trace()
        t = self
        b = self._parent
        while b is not base:
            t = t.trace()
            if b is b.base_ring():
                raise ValueError("(%s) is not defined over (%s)" % (self, base))
            b = b.base_ring()
        return t

    def norm(self, base=None):
        if base is None or base is self._parent._base:
            return self.matrix().determinant()
        n = self
        b = self._parent
        while b is not base:
            n = n.norm()
            if b is b.base_ring():
                raise ValueError("(%s) is not defined over (%s)" % (self, base))
            b = b.base_ring()
        return n

    def charpoly(self, base=None, var='x'):
        return self.matrix(base).charpoly(var)
