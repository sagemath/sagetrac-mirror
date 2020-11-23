r"""
Elements of finitely generated free graded modules

This class implements construction and basic manipulation of 
elements of the Sage parent
:class:`sage.modules.fp_over_steenrod_algebra.free_module.FreeModule`, which models
free graded modules over connected algebras.

.. NOTE:: This class is intended for private use by
    :class:`sage.modules.fp_over_steenrod_algebra.fp_module.FP_Module`.

For an overview of the free module API, see :doc:`free_module`.

AUTHORS:

    - Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
    - Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
      original software to Sage version 8.9.
    - Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added 
      new documentation and tests.

"""

#*****************************************************************************
#       Copyright (C) 2019 Robert R. Bruner <rrb@math.wayne.edu>
#                     and  Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.element import ModuleElement as SageModuleElement


class FreeModuleElement(SageModuleElement):

    def __init__(self, module, coefficients):
        r"""
        Create a module element of a finitely generated free graded module over
        a connected graded algebra.

        INPUT::

        - ``module`` -- the parent instance of this module element.

        - ``coefficients`` -- a tuple of homogeneous algebra coefficients.

        OUTPUT:: The module element given by the coefficients.

        .. NOTE:: This constructor should not be used explicitly, instead use
              the parent's call method.  The reason for this is that the
              dynamic type of the element class changes as a consequence of the
              category system.

        EXAMPLES::  

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import FreeModule
            sage: M = FreeModule((0, 1), SteenrodAlgebra(2))

            sage: M([0, 0])
            <0, 0>

            sage: M([1, 0])
            <1, 0>

            sage: M([0, 1])
            <0, 1>

            sage: M([Sq(1), 1])
            <Sq(1), 1>

        """
        if isinstance(coefficients, FreeModuleElement):
            self._coefficients = coefficients._coefficients
        else:
            self._coefficients = tuple([module.base_ring()(x) for x in coefficients])

        if len(self._coefficients) != len(module.generator_degrees()):
            raise ValueError('the number of coefficients must match the '
                'number of module generators: %d' % len(module.generator_degrees()))

        # Check homogenity and store the degree of the element.
        self._degree = None
        for g, c in zip(module.generator_degrees(), self._coefficients):
            if not c.is_zero():
                d = g + c.degree()

                # XXX todo: Measure how much time is spent in this loop.  Since
                #           construction of free module elements is done
                #           frequently e.g. when computing resolutions, we could
                #           potentially speed things up by dropping this
                #           test.  Since this class constructor is for internal
                #           use only, we could justify commenting in the
                #           following break statement:
                # self._degree = d
                # break

                if self._degree == None:
                    self._degree = d
                else:
                    if self._degree != d:
                        raise ValueError('non-homogeneous element defined')

        SageModuleElement.__init__(self, parent=module)


    def coefficients(self):
        r"""
        The coefficients of this module element.

        OUTPUT:: A tuple of elements of the algebra over which this module is
        defined.

        EXAMPLES::  

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import FreeModule
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeModule((0,1), A)
            sage: x = M.element_from_coordinates((0,0,0,1), 5); x
            <0, Sq(4)>
            sage: x.coefficients()
            (0, Sq(4))
            sage: y = M.element_from_coordinates((0,0,0,0), 5); y
            <0, 0>
            sage: y.coefficients()
            (0, 0)

        """
        return self._coefficients


    def degree(self):
        r"""
        The degree of this element.

        OUTPUT:: the integer degree of this element, or None if this is the zero
        element.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeModule((0,1), A)
            sage: x = M.an_element(7); x
            <Sq(0,0,1), Sq(3,1)>
            sage: x.degree()
            7

        The zero element has no degree::

            sage: (x-x).degree() is None
            True

        """
        return self._degree


    def _repr_(self):
        r"""
        Return a string representation of this element.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A = SteenrodAlgebra(2)
            sage: M = FreeModule((0,1), A)
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
        return '<%s>' % ', '.join(['%s' % c for c in self._coefficients])


    def _lmul_(self, a):
        r"""
        Act by left multiplication on this element by ``a``.

        INPUT::

        - ``a`` -- an element of the algebra this module is defined over.

        OUTPUT:: the module element `a\cdot x` where `x` is this module element.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FreeModule((0,0,3), A2)
            sage: A2.Sq(2)*M.generator(1)
            <0, Sq(2), 0>
            sage: A2.Sq(2)*(A2.Sq(1)*A2.Sq(2)*M.generator(1) + M.generator(2))
            <0, Sq(2,1), Sq(2)>

        TESTS:

            sage: elements = [M.an_element(n) for n in range(1,10)]
            sage: a = A2.Sq(3)
            sage: [a*x for x in elements]
            [<Sq(1,1), Sq(1,1), 0>,
             <0, 0, 0>,
             <Sq(3,1), Sq(3,1), Sq(3)>,
             <0, 0, Sq(1,1)>,
             <0, 0, 0>,
             <Sq(3,2), Sq(3,2), Sq(3,1)>,
             <Sq(3,0,1), Sq(3,0,1), Sq(7)>,
             <Sq(1,1,1), Sq(1,1,1), Sq(5,1)>,
             <0, 0, Sq(3,2)>]

        """

        return self.parent()((a*c for c in self._coefficients))


    def _neg_(self):
        r"""
        Return the additive inverse of this element.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FreeModule((0,), A2)
            sage: x = M.an_element(6);x
            <Sq(0,2)>
            sage: -x
            <Sq(0,2)>
            sage: x + (-x) == 0
            True

        """
        return self.parent()([-c for c in self._coefficients])


    def _add_(self, other):
        r"""
        Return the module sum of this and the given module element.

        Implementation of this function allows Sage to make sense of the +
        operator for instances of this class.

        INPUT::

        - ``other`` -- another element of this element's module.  Only elements
          of the same degree are allowed to be added together.

        OUTPUT:: the module sum of this element and the given element ``other``.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FreeModule((0,), A2)
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

        if self.parent() != other.parent():
            raise TypeError('can not add element in different modules')
        elif self._degree == None: # if self = 0, degree is None
            return self.parent()(other.coefficients())
        elif other._degree == None:   # if other = 0, degree is None
            return self.parent()(self._coefficients)
        elif self._degree != other._degree:
            raise ValueError('can not add element of degree %s and %s'\
                  %(self._degree, other._degree))
        else:
            return self.parent()(
                [x + y for x,y in zip(self._coefficients, other.coefficients())])


    def _richcmp_(self, other, op):
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

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FreeModule((0,1), A2)
            sage: x = M([Sq(1), 1]); x
            <Sq(1), 1>
            sage: y = M([0, Sq(1)]);y
            <0, Sq(1)>

        Multiplying by Sq(1) takes x to the element y::

            sage: A2.Sq(1)*x == y
            True

        TESTS:

            sage: N = FreeModule((0,), A2)
            sage: x._richcmp_(M.an_element(4), op=2)  # Elements of different degrees aren't equal
            False
            sage: w = N.an_element(1)
            sage: x._richcmp_(w, op=2) # Elements of different modules aren't equal.
            False
            sage: x._richcmp_(w, op=3) # Elements of different modules aren't equal.
            True
            sage: z = M.zero()
            sage: x._richcmp_(z, op=2) # Compare the non-trivial x to the zero element.
            False
            sage: x._richcmp_(z, op=3) # Compare the non-trivial x to the zero element.
            True
            sage: z._richcmp_(z, op=2) # Compare the zero element to itself.
            True
            sage: z._richcmp_(z, op=3) # Compare the zero element to itself.
            False

        """

        same = True
        if self.parent() != other.parent():
            same = False
        elif self._degree != other._degree and self._degree != None and other._degree != None:
            same = False
        elif (self._add_(other._neg_()))._nonzero_():
            same = False

        # Equality
        if op == 2:
            return same

        # Non-equality
        if op == 3:
            return not same

        return False

    @cached_method
    def vector_presentation(self):
        r"""
        A coordinate vector representing this module element when it is non-zero.

        These are coordinates with respect to the basis chosen by
        :meth:`sage.modules.fp_over_steenrod_algebra.free_module.FreeModule.basis_elements`.
        When the element is zero, it has no well defined degree, and this
        function returns ``None``.

        OUTPUT:: A vector of elements in the ground field of the algebra for
        this module when this element is non-zero.  Otherwise, the value
        ``None``.

        .. SEEALSO::

            :meth:`sage.modules.fp_over_steenrod_algebra.free_module.FreeModule.vector_presentation`
            :meth:`sage.modules.fp_over_steenrod_algebra.free_module.FreeModule.basis_elements`
            :meth:`sage.modules.fp_over_steenrod_algebra.free_module.FreeModule.element_from_coordinates`

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FreeModule((0,1), A2)
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

        # We cannot represent the zero element since it does not have a degree,
        # and we therefore do not know which vectorspace it belongs to.
        # 
        # In this case, we could return the integer value 0 since coercion would
        # place it inside any vectorspace.  However, this will not work for
        # homomorphisms, so we we return None to be consistent.
        if self._degree is None:
             return None

        bas_gen = self.parent().basis_elements(self._degree)
        base_vec = self.parent().vector_presentation(self._degree)

        base_dict = dict(zip(bas_gen, base_vec.basis()))

        # Create a sparse representation of the element.
        sparse_coeffs = [x for x in enumerate(self._coefficients) if not x[1].is_zero()]

        vector = base_vec.zero()
        for summand_index, algebra_element in sparse_coeffs:
            for scalar_coefficient, monomial in zip(algebra_element.coefficients(), algebra_element.monomials()):
                vector += scalar_coefficient*base_dict[monomial*self.parent().generator(summand_index)]

        return vector


    def _nonzero_(self):
        r"""
        Determine if this element is non-zero.

        OUTPUT:: The boolean value True if this element is non-zero, and False
        otherwise.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FreeModule((0,1), A2)
            sage: y = M.an_element(12); y
            <Sq(2,1,1), Sq(4,0,1)>
            sage: y._nonzero_()
            True
            sage: x = M([0, Sq(1)])
            sage: x._nonzero_()
            True
            sage: (A2.Sq(1)*x)._nonzero_()
            False

        """

        if self._degree == None:
            return False

        return not all(c == 0 for c in self._coefficients)


    def __hash__(self):
        r"""
        A hash value representing this element.

        TESTS:

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import *
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: M = FreeModule((0,1), A2)
            sage: M([Sq(3), Sq(2)]).__hash__() == M([Sq(1)*Sq(2), Sq(2)]).__hash__()
            True

        """
        return hash(self._coefficients)

