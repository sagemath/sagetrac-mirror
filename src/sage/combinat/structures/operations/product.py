# -*- coding: utf-8 -*-
"""
Cauchy product of classes of combinatorial structures.

Let `F` and `G` be both classes of combinatorial structures.

We denote `H = F \times G` the (Cauchy) product of `F` and `G`.

AUTHOR:

- Jean-Baptiste Priez (2014)
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import product, imap, izip
from sage.misc.ascii_art import ascii_art_list
from sage.structure.list_clone import ClonableArray
from sage.combinat.integer_vector import IntegerVectors
from sage.misc.misc_c import prod
from sage.combinat.structures import Structures, Structure
from sage.combinat.structures.operations import _Operations


class CauchyProduct(_Operations):
    """
    Cauchy product of classes of combinatorial structures.

    Let `F` and `G` be both classes of combinatorial structures.

    We denote `H = F \times G` the (Cauchy) product of `F` and `G`.
    """

    def _repr_(self):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: C = Compositions()
            sage: CB = C*B; CB
            Product of structures : `Compositions of non-negative integers`, `Binary trees`
        """

        return "Product of structures : `" + "`, `".join(map(repr, self._structures)) + "`"

    def generating_series(self):
        """
        The generating series `h` of `H`: the product of `F` and `G` is
        defined by the product of its generating series:

        MATH::

            h(t) = f(t) \times g(t)

        """
        return prod(map(lambda F: F.generating_series(), self._structures))


    class GradedComponent(Structures.GradedComponent):

        def __iter__(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: C = Compositions()
                sage: CB = C*B
                sage: for obj in CB.graded_component(2).list(): obj
                [[1, 1], .]
                [[2], .]
                [[1], [., .]]
                [[], [., [., .]]]
                [[], [[., .], .]]

            """
            structs = self.ambient()._structures
            for I in IntegerVectors(self.grading(), len(structs)):
                for tup in product(*imap(
                        lambda (F, ci): F.graded_component(ci),
                        izip(structs, I)
                )):
                    yield self._element_constructor_(tup)
        
    class Element(Structure, ClonableArray):
        """
        TESTS::

            sage: B = BinaryTrees()
            sage: C = Compositions()
            sage: CB = C*B
            sage: obj = CB.graded_component(2).first(); obj
            [[1, 1], .]
            sage: obj.parent() is CB
            True

        """

        def check(self):
            pass

        def _ascii_art_(self):
            """
            TESTS::

                sage: B = BinaryTrees()
                sage: C = Compositions()
                sage: CB = C*B
                sage: ascii_art(CB.graded_component(2).list())
                [                              [       ]  [       ] ]
                [                              [ , o   ]  [ ,   o ] ]
                [ [ *   ]  [      ]  [      ]  [    \  ]  [    /  ] ]
                [ [ *,  ], [ **,  ], [ *, o ], [     o ], [   o   ] ]
            """
            return ascii_art_list(self)