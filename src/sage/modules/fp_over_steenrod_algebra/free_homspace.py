r"""
The set of homomorphisms of finitely generated free graded modules

This class implements methods for construction and basic
manipulation of homsets of finitely generated free graded modules over a connected
graded `k`-algebra, where `k` is a field.

.. NOTE:: This class is intended for private use by
    :class:`sage.modules.fp_over_steenrod_algebra.fp_homspace.FP_ModuleHomspace`.

For an overview of the free module API, see :doc:`free_module`.

TESTS:

    sage: from sage.modules.fp_over_steenrod_algebra.free_module import FreeModule
    sage: from sage.misc.sage_unittest import TestSuite
    sage: A = SteenrodAlgebra(2)
    sage: F1 = FreeModule((1,3), A);
    sage: F2 = FreeModule((2,3), A);
    sage: homset = Hom(F1, F2); homset
    Set of Morphisms from Finitely presented free module on 2 generators ...
    sage: homset([F2((Sq(1), 1)), F2((0, Sq(2)))])
    Module homomorphism of degree 2 defined by sending the generators
      [<1, 0>, <0, 1>]
    to
      [<Sq(1), 1>, <0, Sq(2)>]
    sage: TestSuite(homset).run(verbose=True)
    running ._test_additive_associativity() . . . pass
    running ._test_an_element() . . . pass
    running ._test_cardinality() . . . pass
    running ._test_category() . . . pass
    running ._test_elements() . . .
      Running the test suite of self.an_element()
      running ._test_category() . . . pass
      running ._test_eq() . . . pass
      running ._test_new() . . . pass
      running ._test_nonzero_equal() . . . pass
      running ._test_not_implemented_methods() . . . pass
      running ._test_pickling() . . . pass
      pass
    running ._test_elements_eq_reflexive() . . . pass
    running ._test_elements_eq_symmetric() . . . pass
    running ._test_elements_eq_transitive() . . . pass
    running ._test_elements_neq() . . . pass
    running ._test_eq() . . . pass
    running ._test_new() . . . pass
    running ._test_not_implemented_methods() . . . pass
    running ._test_pickling() . . . pass
    running ._test_some_elements() . . . pass
    running ._test_zero() . . . pass

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

from __future__ import absolute_import

from sage.categories.homset import Homset
from sage.misc.cachefunc import cached_method


def is_FreeModuleHomspace(x):
    r"""
    Check if the given object is of type FreeModuleHomspace.

    OUTPUT:: The boolean ``True`` if and only if ``x`` is of type
    FreeModuleHomspace, and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.modules.fp_over_steenrod_algebra.free_module import FreeModule
        sage: from sage.modules.fp_over_steenrod_algebra.free_homspace import is_FreeModuleHomspace
        sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
        sage: F = FreeModule((1,3), A2)
        sage: L = FreeModule((2,3), A2)
        sage: is_FreeModuleHomspace(Hom(F, L))
        True

    TESTS:

        sage: is_FreeModuleHomspace(0)
        False

    """
    return isinstance(x, FreeModuleHomspace)

from .free_morphism import FreeModuleMorphism


class FreeModuleHomspace(Homset):
    # In the category framework, Elements of the class FP_Module are of the
    # class FP_Element, see
    # http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html#implementing-the-category-framework-for-the-elements
    Element = FreeModuleMorphism

    def _element_constructor_(self, values):
        r"""
        Construct any element of this homset.

        This function is used internally by the ()-method when creating
        homomorphisms.

        INPUT::

        - ``values`` -- A tuple of values (i.e. elements of the
        codomain for this homset) corresponding bijectively to the generators
        of the domain of this homset, or the zero integer constant.

        OUTPUT:: An instance of the morphism class.  The returned morphism is
        defined by mapping the module generators in the domain to the given
        values.

        OUTPUT:: A module homomorphism.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import FreeModule
            sage: from sage.modules.fp_over_steenrod_algebra.free_homspace import is_FreeModuleHomspace
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = FreeModule((1,3), A2)
            sage: L = FreeModule((2,5), A2)
            sage: H = Hom(F, L)

            sage: values = (A2.Sq(4)*L.generator(0), A2.Sq(3)*L.generator(1))
            sage: f = H(values); f
            Module homomorphism of degree 5 defined by sending the generators
              [<1, 0>, <0, 1>]
            to
              [<Sq(4), 0>, <0, Sq(3)>]

            sage: H(0)
            The trivial homomorphism.

        """
        if isinstance(values, FreeModuleMorphism):
            return values
        elif values == 0:
            return self.zero()
        else:
            return self.element_class(self, values)


    def _an_element_(self):
        r"""
        Return a morphism belonging to this homspace.

        OUTPUT:: A morphism in this homspace.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import FreeModule
            sage: from sage.modules.fp_over_steenrod_algebra.free_homspace import is_FreeModuleHomspace
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = FreeModule((1,3), A2)
            sage: L = FreeModule((2,3), A2)
            sage: H = Hom(F, L)
            sage: H._an_element_()
            The trivial homomorphism.

        """
        return self.zero()


    @cached_method
    def zero(self):
        r"""
        Return the trivial morphism of this homspace.

        OUTPUT:: The morphism evaluating to the zero element for any element in
        the domain.

        EXAMPLES::

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import FreeModule
            sage: from sage.modules.fp_over_steenrod_algebra.free_homspace import is_FreeModuleHomspace
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: F = FreeModule((1,3), A2)
            sage: L = FreeModule((2,3), A2)
            sage: H = Hom(F, L)
            sage: H.zero()
            The trivial homomorphism.

        """
        return self.element_class(self, self.codomain().zero())


    def identity(self):
        r"""
        Return the identity morphism, if this is an endomorphism set.

        OUTPUT:: The identity endomorphism.

        TESTS:

            sage: from sage.modules.fp_over_steenrod_algebra.free_module import FreeModule
            sage: from sage.modules.fp_over_steenrod_algebra.free_homspace import is_FreeModuleHomspace
            sage: A2 = SteenrodAlgebra(2, profile=(3,2,1))
            sage: L = FreeModule((2,3), A2)
            sage: H = Hom(L, L)
            sage: H.identity()
            The identity homomorphism.

        TESTS:

            sage: F = FreeModule((1,3), A2)
            sage: H = Hom(F, L)
            sage: H.identity()
            Traceback (most recent call last):
            ...
            TypeError: This homspace does not consist of endomorphisms.

        """
        if self.is_endomorphism_set():
            return self.element_class(self, self.codomain().generators())
        else:
            raise TypeError("This homspace does not consist of endomorphisms.")

