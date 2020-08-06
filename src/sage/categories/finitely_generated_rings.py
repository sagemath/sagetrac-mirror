r"""
Finitely generated rings


A ring is `R` *finitely generated* if there exists a finite subset
`A \subset R` such that every element of `R` can be written as a sum of
elements of the form:

.. MATH::

    a_1a_2 \cdots a_n

with `a_i \in A` for `i=1, \ldots, n`.

AUTHORS:

- Reimundo Heluani (08-06-2020): Initial implementation
"""

#*****************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.rings import Rings
from sage.misc.abstract_method import abstract_method

class FinitelyGeneratedRings(CategoryWithAxiom):
    """
    The category of finitely generated rings.

    EXAMPLES::

        sage: C = Rings().FinitelyGenerated(); C
        Category of finitely generated rings

    TESTS::

        sage: TestSuite(C).run()
    """

    _base_category_class_and_axiom = (Rings, "FinitelyGeneratedAsDistributiveMagma")

    class ParentMethods:

        @abstract_method
        def gens(self):
            """
            The generators of this ring.

            EXAMPLES::

                sage: R.<x,y> = ZZ[]
                sage: R.gens()
                (x, y)
            """

        def ngens(self):
            """
            The number of generators of this ring.

            EXAMPLES::

                sage: R.<x,y> = ZZ[]
                sage: R.ngens()
                2
            """
            return len(self.gens())


