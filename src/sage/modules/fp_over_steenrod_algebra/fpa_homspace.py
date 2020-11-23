r"""
The set of homomorphisms of finitely presented graded modules over the Steenrod algebra

This class implements methods for construction and basic
manipulation of homsets of finitely presented graded modules over the
`\operatorname{mod} p` Steenrod algebra.  

For an overview of the API, see :doc:`fpa_module`.

TESTS:

    sage: from sage.modules.fp_over_steenrod_algebra.fpa_module import FPA_Module
    sage: from sage.misc.sage_unittest import TestSuite
    sage: A = SteenrodAlgebra(2, profile=(3,2,1))
    sage: F = FPA_Module([1,3], A)
    sage: L = FPA_Module([2,3], A, [[Sq(2),Sq(1)], [0,Sq(2)]])
    sage: homset = Hom(F, L); homset
    Set of Morphisms from Finitely presented module on 2 generators ...
    sage: homset.an_element()
    Module homomorphism of degree 0 defined by sending the generators
      [<1, 0>, <0, 1>]
    to
      [<0, 0>, <Sq(1), 0>]
    sage: homset([L((Sq(1), 1)), L((0, Sq(2)))])
    Module homomorphism of degree 2 defined by sending the generators
      [<1, 0>, <0, 1>]
    to
      [<Sq(1), 1>, <0, Sq(2)>]
    sage: Hom(F, L) ([L((Sq(1), 1)), L((0, Sq(2)))]).kernel()
    Module homomorphism of degree 0 defined by sending the generators
      [<1, 0>, <0, 1>]
    to
      [<0, 1>, <Sq(0,1), 0>]
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

from .fpa_morphism import FPA_ModuleMorphism
from .fp_homspace import FP_ModuleHomspace


class FPA_ModuleHomspace(FP_ModuleHomspace):
    # In the category framework, Elements of the class FPA_ModuleHomspace are of the
    # class FPA_ModuleMorphism, see
    # http://doc.sagemath.org/html/en/thematic_tutorials/coercion_and_categories.html#implementing-the-category-framework-for-the-elements
    Element = FPA_ModuleMorphism

