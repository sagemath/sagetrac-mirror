"""
Lie Algebra Elements

AUTHORS:

- Travis Scrimshaw (2005-05-04): Initial implementation
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#from sage.misc.abstract_method import abstract_method
#from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.misc import repr_lincomb
from copy import copy
#from functools import total_ordering
#from sage.structure.element import ModuleElement, RingElement, coerce_binop
#from sage.structure.sage_object import SageObject
from sage.combinat.free_module import CombinatorialFreeModuleElement
from sage.structure.element cimport have_same_parent, coercion_model
from sage.structure.element_wrapper cimport ElementWrapper

# TODO: Have the other classes inherit from this?
# TODO: Should this be a mixin class (or moved to the category)?
#class LieAlgebraElement_generic(ModuleElement):
#    """
#    Generic methods for all Lie algebra elements.
#    """
#    def __mul__(self, other):
#        """
#        If we are multiplying two non-zero elements, automatically
#        lift up to the universal enveloping algebra.
#
#        EXAMPLES::
#        """
#        if self == 0 or other == 0:
#            return self.parent().zero()
#        # Otherwise we lift to the UEA
#        return self.lift() * other

# TODO: Factor out parts of CombinatorialFreeModuleElement into a SparseFreeModuleElement?
# TODO: Do we want a dense version?
class LieAlgebraElement(CombinatorialFreeModuleElement):
    """
    A Lie algebra element.
    """
    # Need to bypass the coercion model
    def __mul__(self, y):
        """
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: y*x
            x*y - z
        """
        if self == 0 or y == 0:
            return self.parent().zero()
        if y in self.base_ring():
            return y * self
        # Otherwise we lift to the UEA
        return self.lift() * y

    #def _im_gens_(self, codomain, im_gens):
    #    """
    #    Return the image of ``self`` in ``codomain`` under the map that sends
    #    the images of the generators of the parent of ``self`` to the
    #    tuple of elements of ``im_gens``.
    #
    #    EXAMPLES::
    #    """
    #    s = codomain.zero()
    #    if not self: # If we are 0
    #        return s
    #    names = self.parent().variable_names()
    #    return codomain.sum(c * t._im_gens_(codomain, im_gens, names)
    #                        for t, c in self._monomial_coefficients.iteritems())

    # TODO: Move to the category/lift morphism?
    def lift(self):
        """
        Lift ``self`` to the universal enveloping algebra.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: x.lift().parent() == L.universal_enveloping_algebra()
            True
        """
        UEA = self.parent().universal_enveloping_algebra()
        gen_dict = UEA.gens_dict()
        s = UEA.zero()
        if not self:
            return s
        for t, c in self._monomial_coefficients.iteritems():
            s += c * gen_dict[t]
        return s

    def is_constant(self):
        """
        Check if ``self`` is a constant (i.e. zero).

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: a = x + y
            sage: a.is_constant()
            False
            sage: L.zero().is_constant()
            True
        """
        return not self._monomial_coefficients

    def dict(self):
        """
        Return ``self`` as a dictionary mapping monomials to coefficients.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: a = 3*x - 1/2*z
            sage: a.dict()
            {'x': 3, 'z': -1/2}
        """
        return copy(self._monomial_coefficients)

    def list(self):
        """
        Return ``self`` as a list of pairs ``(m, c)`` where ``m`` is a
        monomial and ``c`` is the coefficient.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: a = 3*x - 1/2*z
            sage: a.list()
            [('x', 3), ('z', -1/2)]
        """
        return sorted(self._monomial_coefficients.items())

cdef class LieAlgebraElementWrapper(ElementWrapper):
    """
    Wrap an element as a Lie algebra element.
    """
    def __richcmp__(self, right, int op):
        """
        Perform a rich comparison.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
            sage: L.bracket(L.gen(0), L.gen(1)) == -L.bracket(L.gen(1), L.gen(0))
            True

        The next doctests show similar behavior, although on elements of
        other classes::

            sage: L = lie_algebras.three_dimensional_by_rank(QQ, 3)
            sage: L.bracket(L.gen(0), L.gen(1)) == -L.bracket(L.gen(1), L.gen(0))
            True

            sage: L = lie_algebras.three_dimensional_by_rank(QQ, 1)
            sage: L.bracket(L.gen(0), L.gen(1)) == -L.bracket(L.gen(1), L.gen(0))
            True

        Check inequality::

            sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
            sage: L.bracket(L.gen(0), L.gen(1)) != -L.bracket(L.gen(1), L.gen(0))
            False
            sage: L.zero() == 0
            True
            sage: L.zero() != 0
            False

        The next doctests show similar behavior, although on elements of
        other classes::

            sage: L = lie_algebras.three_dimensional_by_rank(QQ, 3)
            sage: L.bracket(L.gen(0), L.gen(1)) != -L.bracket(L.gen(1), L.gen(0))
            False
            sage: L.an_element()
            sage: L.an_element() == 0
            False
            sage: L.an_element() != 0
            True

            sage: L = lie_algebras.three_dimensional_by_rank(QQ, 1)
            sage: L.bracket(L.gen(0), L.gen(1)) != -L.bracket(L.gen(1), L.gen(0))
            False
            sage: L.zero() == 0
            True
            sage: L.zero() != 0
            False
            sage: L.zero() >= 0
            True
            sage: L.zero() < 0
            False
        """
        if right == self.parent().base_ring().zero():
            if op == 3: # !=
                return self.__nonzero__()
            if op in [1,2,5]: # <=, ==, >=
                return not self.__nonzero__()
            return False # <, >
        if not have_same_parent(self, right):
            try:
                self, right = coercion_model.canonical_coercion(self, right)
            except (TypeError, ValueError):
                return op == 3
        if op == 3: # !=
            return self.value != right.value
        if op in [1,2,5]: # <=, ==, >=
            return self.value == right.value
        return False # <, >

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x + y
            x + y
        """
        return repr(self.value)

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x')
            sage: L.<x0,x1,x2> = LieAlgebra(associative=R.gens())
            sage: latex(x0 + x1)
            x_{0} + x_{1}
        """
        from sage.misc.latex import latex
        return latex(self.value)

    def _ascii_art_(self):
        """
        Return an ascii art representation of ``self``.
        """
        from sage.typeset.ascii_art import ascii_art
        return ascii_art(self.value)

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.
        """
        from sage.typeset.unicode_art import unicode_art
        return unicode_art(self.value)

    def __nonzero__(self):
        """
        Return if ``self`` is non-zero.
        """
        return bool(self.value)

    cpdef _add_(self, right):
        """
        Add ``self`` and ``rhs``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x + y
            x + y
        """
        return self.__class__(self.parent(), self.value + right.value)

    cpdef _sub_(self, right):
        """
        Subtract ``self`` and ``rhs``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: x - y
            x - y
        """
        return self.__class__(self.parent(), self.value - right.value)

    # We need to bypass the coercion framework
    # We let the universal enveloping algebra handle the rest if both
    #   arguments are non-zero
    def __mul__(self, x):
        """
        If we are multiplying two non-zero elements, automatically
        lift up to the universal enveloping algebra.

        .. TODO::

            Write tests for this method once :trac:`16822` is
            implemented.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: S = GroupAlgebra(G, QQ)
            sage: L.<x,y> = LieAlgebra(associative=S.gens())
            sage: u = x*3; u
            3*(1,2,3)
            sage: parent(u) == L
            True
            sage: u = x*(3/2); u
            3/2*(1,2,3)
            sage: parent(u) == L
            True
            sage: elt = x*y - y*x; elt  # not tested: needs #16822
            sage: S(elt)  # not tested: needs #16822
            (2,3) - (1,3)
        """
        if self.value == 0 or x == 0:
            return self.parent().zero()
        if x in self.base_ring():
            return x * self
        # Otherwise we lift to the UEA
        return self.lift() * x

    def __div__(self, x):
        """
        Division by coefficients.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 3)
            sage: x = L.an_element(); x
            p1 + p2 + p3 + q1 + q2 + q3 + z
            sage: x / 2
            1/2*p1 + 1/2*p2 + 1/2*p3 + 1/2*q1 + 1/2*q2 + 1/2*q3 + 1/2*z
        """
        return self * (~x)

    def _acted_upon_(self, scalar, self_on_left=False):
        """
        Return the action of a scalar on ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: 3*x
            3*x
            sage: parent(3*x) == parent(x)
            True
            sage: x / 2
            1/2*x
            sage: y * (1/2)
            1/2*y
            sage: y * 1/2
            1/2*y
            sage: 1/2 * y
            1/2*y
            sage: QQ(1/2) * y
            1/2*y
        """
        # This was copied and IDK if it still applies (TCS):
        # With the current design, the coercion model does not have
        # enough information to detect apriori that this method only
        # accepts scalars; so it tries on some elements(), and we need
        # to make sure to report an error.
        if hasattr( scalar, 'parent' ) and scalar.parent() != self.base_ring():
            # Temporary needed by coercion (see Polynomial/FractionField tests).
            if self.base_ring().has_coerce_map_from(scalar.parent()):
                scalar = self.base_ring()( scalar )
            else:
                return None
        if self_on_left:
            return self.__class__(self.parent(), self.value * scalar)
        return self.__class__(self.parent(), scalar * self.value)

    def __neg__(self):
        """
        Return the negation of ``self``.

        EXAMPLES::

            sage: R = FreeAlgebra(QQ, 3, 'x,y,z')
            sage: L.<x,y,z> = LieAlgebra(associative=R.gens())
            sage: -x
            -x
        """
        return self.__class__(self.parent(), -self.value)

    def __getitem__(self, i):
        """
        Redirect the ``__getitem__()`` to the wrapped element.

        EXAMPLES::

            sage: L = lie_algebras.sl(QQ, 2, representation='matrix')
            sage: m = L.gen(0)
            sage: m[0,0]
            0
            sage: m[0][1]
            1
        """
        return self.value.__getitem__(i)

# TODO: Also used for vectors, find a better name
cdef class LieAlgebraMatrixWrapper(LieAlgebraElementWrapper):
    def __init__(self, parent, value):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = lie_algebras.Heisenberg(QQ, 1, representation="matrix")
            sage: z = L.z()
            sage: z.value.is_immutable()
            True
        """
        value.set_immutable() # Make the matrix immutable for hashing
        LieAlgebraElementWrapper.__init__(self, parent, value)

cdef class StructureCoefficientsElement(LieAlgebraMatrixWrapper):
    """
    An element of a Lie algebra given by structure coefficients.
    """
    def _repr_(self):
        """
        EXAMPLES::
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult=self.parent()._print_options['scalar_mult'],
                            repr_monomial=self.parent()._repr_generator,
                            strip_one=True)

    def _latex_(self):
        """
        EXAMPLES::
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult=self.parent()._print_options['scalar_mult'],
                            latex_scalar_mult=self.parent()._print_options['latex_scalar_mult'],
                            repr_monomial=self.parent()._latex_term,
                            is_latex=True, strip_one=True)

    cpdef bracket(self, right):
        """
        Return the Lie bracket ``[self, right]``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}, ('y','z'): {'x':1}, ('z','x'): {'y':1}})
            sage: x.bracket(y)
            z
            sage: y.bracket(x)
            -z
            sage: (x + y - z).bracket(x - y + z)
            -2*y - 2*z
        """
        if not have_same_parent(self, right):
            self, right = coercion_model.canonical_coercion(self, right)
        return self._bracket_(right)

    # We need this method because the LieAlgebra.bracket method (from the
    #   category) calls this, where we are guaranteed to have the same parent.
    cpdef _bracket_(self, right):
        """
        Return the Lie bracket ``[self, right]``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}, ('y','z'): {'x':1}, ('z','x'): {'y':1}})
            sage: x._bracket_(y)
            z
            sage: y._bracket_(x)
            -z
        """
        P = self.parent()
        cdef dict s_coeff = P._s_coeff
        d = P.dimension()
        cdef list ret = [P.base_ring().zero()]*d
        cdef int i1, i2, i3
        for i1 in range(d):
            c1 = self.value[i1]
            if not c1:
                pass
            for i2 in range(d):
                c2 = right.value[i2]
                if not c2:
                    pass
                if (i1, i2) in s_coeff:
                    for i3 in range(d):
                        ret[i3] += c1 * c2 * s_coeff[i1, i2][i3]
                elif (i2, i1) in s_coeff:
                    for i3 in range(d):
                        ret[i3] -= c1 * c2 * s_coeff[i2, i1][i3]
        return self.__class__(P, P._M(ret))

    def __iter__(self):
        """
        Iterate over ``self``.
        """
        zero = self.base_ring().zero()
        I = self.parent()._indices
        cdef int i
        for i,v in enumerate(self.value):
            if v != zero:
                yield (I[i], v)

    cpdef to_vector(self):
        """
        Return ``self`` as a vector.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: a = x + 3*y - z/2
            sage: a.to_vector()
            (1, 3, -1/2)
        """
        return self.value

    def lift(self):
        """
        Return the lift of ``self`` to the universal enveloping algebra.

        EXAMPLES::
        """
        UEA = self.parent().universal_enveloping_algebra()
        gens = UEA.gens()
        return UEA.sum(c * gens[i] for i, c in self.value.iteritems())

    cpdef dict monomial_coefficients(self, bint copy=True):
        """
        Return the monomial coefficients of ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: a = 2*x - 3/2*y + z
            sage: list(a)
            [('x', 2), ('y', -3/2), ('z', 1)]
            sage: a = 2*x - 3/2*z
            sage: list(a)
            [('x', 2), ('z', -3/2)]
        """
        I = self.parent()._indices
        return {I[i]: v for i,v in self.value.monomial_coefficients()}

    def __getitem__(self, i):
        """
        Return the coefficient of the basis element indexed by ``i``.

        EXAMPLES::
        """
        return self.value[self.parent()._indices.index(i)]

