r"""
Elements of finitely presented graded modules

This class implements construction and basic manipulation of elements of the
Sage parent :class:`sage.modules.fp_over_steenrod_algebra.fp_module.FP_Module`, which models
finitely presented modules over connected graded algebras.

.. NOTE:: This class is used by the derived class
    :class:`sage.modules.fp_over_steenrod_algebra.fpa_element.FPA_Element`.

AUTHORS:

    - Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
    - Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
      original software to Sage version 8.9.
    - Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added 
      new documentation and tests.

"""

#*****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu> and
#                          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.element import ModuleElement as SageModuleElement

from .free_element import FreeModuleElement

from .timing import g_timings

cdef class FP_Element():

    cdef object _free_element
    cdef object _parent

    def __init__(self, module, coefficients):
        r"""
        Create a module element of a finitely presented graded module over
        a connected graded algebra.
        """
        # Store the free representation of the element.
        self._free_element = FreeModuleElement(module.j.codomain(), coefficients)

        self._parent = module

#        SageModuleElement.__init__(self, parent=module)

    def parent(self):
        return self._parent


    def free_element(self):
        return self._free_element


    def coefficients(self):
        r"""
        The coefficients of this module element.

        OUTPUT:: A tuple of elements of the algebra over which this module is
        defined.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: M = FP_Module([0,1], SteenrodAlgebra(2), [[Sq(4), Sq(3)]])
            sage: x = M.element_from_coordinates((0,0,1), 5)

            sage: x
            <0, Sq(4)>
            sage: x.coefficients()
            (0, Sq(4))

            sage: y = M.element_from_coordinates((0,0,0), 5)
            sage: y
            <0, 0>
            sage: y.coefficients()
            (0, 0)

        """
        return self._free_element.coefficients()


    @cached_method
    def degree(self):
        r"""
        The degree of this element.

        OUTPUT:: The integer degree of this element, or ``None`` if this is the
        zero element.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: M = FP_Module([0,1], SteenrodAlgebra(2), [[Sq(4), Sq(3)]])
            sage: x = M.an_element(7)

            sage: x
            <Sq(0,0,1), Sq(3,1)>
            sage: x.degree()
            7
    
            sage: # The zero element has no degree::
            sage: (x-x).degree() is None
            True

        TESTS:

            sage: N = FP_Module([0], SteenrodAlgebra(2), [[Sq(2)]])
            sage: y = Sq(2)*N.generator(0)
            sage: y == 0
            True
            sage: y.degree() is None
            True

        """
        return self._free_element.degree() if self.nonzero() else None


    def __repr__(self):
        r"""
        Return a string representation of this element.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: M = FP_Module([0,1], SteenrodAlgebra(2), [[Sq(4), Sq(3)]])
            sage: [M.an_element(n) for n in range(1,10)]
            [<Sq(1), 1>,
             <Sq(2), Sq(1)>,
             <Sq(0,1), Sq(2)>,
             <Sq(1,1), Sq(3)>,
             <Sq(2,1), Sq(4)>,
             <Sq(0,2), Sq(5)>,
             <Sq(0,0,1), Sq(3,1)>,
             <Sq(1,0,1), Sq(1,2)>,
             <Sq(2,0,1), Sq(2,2)>]

        """
        return self._free_element.__repr__()


    def _lmul_(self, a):
        r"""
        Act by left multiplication on this element by ``a``.

        INPUT::

        - ``a`` -- an element of the algebra this module is defined over.

        OUTPUT:: the module element `a\cdot x` where `x` is this module element.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FP_Module([0,3], A2, [[Sq(2)*Sq(4), Sq(3)]])
            sage: A2.Sq(2)*M.generator(1)
            <0, Sq(2)>
            sage: A2.Sq(2)*(A2.Sq(1)*A2.Sq(2)*M.generator(0) + M.generator(1))
            <Sq(2,1), Sq(2)>

        TESTS:

            sage: elements = [M.an_element(n) for n in range(1,10)]
            sage: a = A2.Sq(3)
            sage: [a*x for x in elements]
            [<Sq(1,1), 0>,
             <0, 0>,
             <Sq(3,1), Sq(3)>,
             <0, Sq(1,1)>,
             <0, 0>,
             <Sq(3,2), Sq(3,1)>,
             <Sq(3,0,1), Sq(7)>,
             <Sq(1,1,1), Sq(5,1)>,
             <0, Sq(3,2)>]

        """
        global g_timings

        xxx = self._free_element._lmul_(a)

        g_timings.Start('fp_element_arithmetic_')
        res = self.parent()(xxx)
        g_timings.End()

        return res


    def __neg__(self):
        r"""
        Return the additive inverse of this element.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FP_Module([0], A2)

            sage: x = M.an_element(6);x
            <Sq(0,2)>

            sage: -x
            <Sq(0,2)>

            sage: x + (-x) == 0
            True

        """
        global g_timings

        xxx = -self._free_element

        g_timings.Start('fp_element_arithmetic_')
        res = self.parent()(xxx)
        g_timings.End()

        return res


    def __add__(self, other):
        r"""
        Return the module sum of this and the given module element.

        Implementation of this function allows Sage to make sense of the +
        operator for instances of this class.

        INPUT::

        - ``other`` -- another element of this element's module.  Only elements
          of the same degree are allowed to be added together.

        OUTPUT:: the module sum of this element and the given element ``other``.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FP_Module([0], A2)

            sage: x = M.an_element(6);x
            <Sq(0,2)>

            sage: -x
            <Sq(0,2)>

            sage: x + (-x) == 0
            True

        TESTS:

            sage: x = M.an_element(4)
            sage: y = M.an_element(5)
            sage: x+y
            Traceback (most recent call last):
            ...
            ValueError: can not add element of degree 4 and 5
            sage: z = M.zero()
            sage: x+z == x
            True
            sage: z+x
            <Sq(1,1)>
            sage: y+z
            <Sq(2,1)>

        """
        global g_timings
        xxx = self._free_element + other.free_element()

        g_timings.Start('fp_element_arithmetic_')
        res = self.parent()(xxx)
        g_timings.End()

        return res


    def __eq__(self, other):
        r"""
        Compare this element with ``other``.

        Implementation of this function allows Sage to make sense of the ==
        operator for instances of this class.

        INPUT::

        - ``other`` -- An instance of this class.

        - ``op`` -- An integer specifying the comparison operation to be
          carried out: If ``op`` == 2, then return ``True`` if and only if the
          elements are equal.  If ``op`` == 3, then return ``True `` if and
          only if the elements are not equal.  Otherwise, return ``False``.

        OUTPUT:: A boolean.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FP_Module((0,1), A2)
            sage: x = M([Sq(1), 1]); x
            <Sq(1), 1>
            sage: y = M([0, Sq(1)]); y
            <0, Sq(1)>

        Multiplying by Sq(1) takes x to the element y::

            sage: A2.Sq(1)*x == y
            True

        TESTS:

            sage: N = FP_Module([0], A2)
            sage: x._richcmp_(M.an_element(4), 2)  # Elements of different degrees aren't equal
            False
            sage: w = N.an_element(1)
            sage: x._richcmp_(w, 2) # Elements of different modules aren't equal.
            False
            sage: x._richcmp_(w, 3) # Elements of different modules aren't equal.
            True
            sage: z = M.zero()
            sage: x._richcmp_(z, 2) # Compare the non-trivial x to the zero element.
            False
            sage: x._richcmp_(z, 3) # Compare the non-trivial x to the zero element.
            True
            sage: z._richcmp_(z, 2) # Compare the zero element to itself.
            True
            sage: z._richcmp_(z, 3) # Compare the zero element to itself.
            False

        """


#        if type(other) is int:
#            return self._coefficients == (other,)
#        if self._coefficients == other.coefficients():
#            return True
#        else:
        return not (self.__add__(other.__neg__())).nonzero()

#
#        same = True
#        if self.parent() != other.parent() or\
#            self.degree() != other.degree() or\
#            (self._add_(other._neg_())).nonzero():
#            same = False
#
#        # Equality
#        if op == 2:
#            return same
#
#        # Non-equality
#        if op == 3:
#            return not same
#
#        return False


    def vector_presentation(self):
        r"""
        A coordinate vector representing this module element when it is non-zero.

        These are coordinates with respect to the basis chosen by
        :meth:`sage.modules.fp_over_steenrod_algebra.fp_module.FP_Module.basis_elements`.
        When the element is zero, it has no well defined degree, and this
        function returns ``None``.

        OUTPUT:: A vector of elements in the ground field of the algebra for
        this module when this element is non-zero.  Otherwise, the value
        ``None``.

        .. SEEALSO::

            :meth:`sage.modules.fp_over_steenrod_algebra.fp_module.FP_Module.vector_presentation`
            :meth:`sage.modules.fp_over_steenrod_algebra.fp_module.FP_Module.basis_elements`
            :meth:`sage.modules.fp_over_steenrod_algebra.fp_module.FP_Module.element_from_coordinates`

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FP_Module((0,1), A2)

            sage: x = M.an_element(7)
            sage: v = x.vector_presentation(); v
            (1, 0, 0, 0, 0, 1, 0)
            sage: type(v)
            <type 'sage.modules.vector_mod2_dense.Vector_mod2_dense'>

            sage: V = M.vector_presentation(7)
            sage: v in V
            True

            sage: M.element_from_coordinates(v, 7) == x
            True

        We can use the basis for the module elements in the degree of `x`,
        together with the coefficients `v` to recreate the element `x`::

            sage: basis = M.basis_elements(7)
            sage: x_ = sum( [c*b for (c,b) in zip(v, basis)] ); x_
            <Sq(0,0,1), Sq(3,1)>
            sage: x == x_
            True

        TESTS:

            sage: M.zero().vector_presentation() is None
            True

        """
        global g_timings

        # We cannot represent the zero element since it does not have a degree,
        # and we therefore do not know which vectorspace it belongs to.
        #
        # In this case, we could return the integer value 0 since coercion would
        # place it inside any vectorspace.  However, this will not work for
        # homomorphisms, so we we return None to be consistent.
        if self._free_element.degree() is None:
            return None

        F_n = self.parent().vector_presentation(self._free_element.degree())

        v = self._free_element.vector_presentation()

        g_timings.Start('lin_alg')
        qv = F_n.quotient_map()(v)
        g_timings.End()

        return qv


    def nonzero(self):
        r"""
        Determine if this element is non-zero.

        OUTPUT:: The boolean value ``True`` if this element is non-zero, and ``False``
        otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: M = FP_Module([0,2,4], SteenrodAlgebra(2), [[Sq(4),Sq(2),0]])
            sage: M(0)._nonzero_()
            False
            sage: M((Sq(6), 0, Sq(2)))._nonzero_()
            True
            sage: a = M((Sq(1)*Sq(2)*Sq(1)*Sq(4), 0, 0))
            sage: b = M((0, Sq(2)*Sq(2)*Sq(2), 0))
            sage: a._nonzero_()
            True
            sage: b._nonzero_()
            True
            sage: (a + b)._nonzero_()
            False

        """
        pres = self.vector_presentation()
        return False if pres is None else (pres != 0)


    def is_zero(self):
        return not self.nonzero()


    def normalize(self):
        r"""
        A normalized form of ``self``.

        OUTPUT:: An instance of this element class representing the same
        module element as this element.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: M = FP_Module([0,2,4], SteenrodAlgebra(2), [[Sq(4),Sq(2),0]])

            sage: m = M((Sq(6), 0, Sq(2))); m
            <Sq(6), 0, Sq(2)>
            sage: m.normalize()
            <Sq(6), 0, Sq(2)>
            sage: m == m.normalize()
            True

            sage: n = M((Sq(4), Sq(2), 0)); n
            <Sq(4), Sq(2), 0>
            sage: n.normalize()
            <0, 0, 0>
            sage: n == n.normalize()
            True

        """
        if not self.nonzero():
            return self.parent().zero()

        v = self.vector_presentation()
        return self.parent().element_from_coordinates(v, self.degree())


    def __hash__(self):
        r"""
        A hash value representing this element.

        TESTS:

            sage: from sage.modules.fp_over_steenrod_algebra.fp_module import FP_Module
            sage: M = FP_Module([0,1], SteenrodAlgebra(2), [[Sq(4),Sq(3)]])
            sage: M([Sq(3), Sq(2)]).__hash__() == M([Sq(1)*Sq(2), Sq(2)]).__hash__()
            True

        """
        return hash(self.coefficients())

