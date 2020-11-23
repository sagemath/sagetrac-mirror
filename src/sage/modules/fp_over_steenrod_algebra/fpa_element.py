r"""
Elements of finitely presented modules over the Steenrod algebra

This class implements construction and basic manipulation of elements of the
Sage parent :class:`sage.modules.fp_over_steenrod_algebra.fpa_module.FPA_Module`, which models
finitely presented modules over the `\operatorname{mod} p` Steenrod algebra.

For an overview of the API, see :doc:`fpa_module`.

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

from .fp_element import FP_Element


class FPA_Element(FP_Element):

    def __init__(self, module, coefficients):
        r"""
        Create a module element of a finitely presented graded module over
        the Steenrod algebra.

        INPUT:

        - ``module`` -- The parent instance of this module element.

        - ``coefficients`` -- A tuple of homogeneous elements of the algebra
          over which the module is defined.

        OUTPUT: The module element given by the coefficients.

        .. NOTE:: Never use this constructor explicitly, but rather the parent's
            call method, or this class' __call__ method.  The reason for this
            is that the dynamic type of the element class changes as a
            consequence of the category system.

        TESTS:

            sage: from sage.modules.fp_over_steenrod_algebra.fpa_module import FPA_Module
            sage: from sage.modules.fp_over_steenrod_algebra.fpa_element import FPA_Element
            sage: FPA_Element(FPA_Module([0], SteenrodAlgebra(2)), [Sq(2)])
            <Sq(2)>

        """
        FP_Element.__init__(self, module, coefficients)


