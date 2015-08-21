# -*- coding: utf-8 -*-
r"""
Shuffle product of trees

AUTHOR:

- Jean-Baptiste Priez
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.combinat import CombinatorialClass


class TreeShuffleProduct(CombinatorialClass):
    """
    TESTS::

        sage: from sage.combinat.hopf_algebras.pbt.tree_shuffle_product import TreeShuffleProduct
        sage: b1 = BinaryTree([[],[]]); ascii_art(b1)
          o
         / \
        o   o
        sage: b2 = BinaryTree([[], None]); ascii_art(b2)
          o
         /
        o
        sage: TreeShuffleProduct(b1, b2)
        Tree shuffle product of [[., .], [., .]] and [[., .], .]
        sage: ascii_art(list(TreeShuffleProduct(b1, b2)))
        [       o       o       o      __o__      _o_      o     ]
        [      /       /       /      /     \    /   \    / \    ]
        [     o      _o_      o      o       o  o     o  o   o   ]
        [    /      /   \    / \            /        /        \  ]
        [   o      o     o  o   o          o        o          o ]
        [  / \          /        \        /          \        /  ]
        [ o   o  ,     o  ,       o,     o    ,       o,     o   ]
        sage: b3 = BinaryTree([None, []]); ascii_art(b3)
        o
         \
          o
        sage: ascii_art(list(TreeShuffleProduct(b1, b3)))
        [     _o_      _o_        o       ]
        [    /   \    /   \      / \      ]
        [   o     o  o     o    o   o     ]
        [  / \            / \        \    ]
        [ o   o    ,     o   o,       o   ]
        [                              \  ]
        [                               o ]
    """

    def __init__(self, tree1, tree2):
        self._t1 = tree1
        self._t2 = tree2

    def __repr__(self):
        return "Tree shuffle product of " + str(self._t1) + \
               " and " + str(self._t2)

    def __contains__(self, x):
        # TODO:: test the instance
        # FIXME:: implement
        return False

    def cardinality(self):
        from sage.rings.arith import binomial

        def des_lr(bt, i, j=0):
            if bt.is_empty():
                return j
            return des_lr(bt[i], i, j + 1)

        i = des_lr(self._t1, 1)
        j = des_lr(self._t2, 0)
        return binomial(i + j, j)

    def __iter__(self):
        """
        TESTS::

            sage: from sage.combinat.hopf_algebras.pbt.tree_shuffle_product import TreeShuffleProduct
            sage: tsp = TreeShuffleProduct(BinaryTree([None, [[],[]]]),BinaryTree([[],[]])); tsp
            Tree shuffle product of [., [[., .], [., .]]] and [[., .], [., .]]
            sage: for bt in list(tsp): print bt
            [[[., [[., .], [., .]]], .], [., .]]
            [[., [[[., .], [., .]], .]], [., .]]
            [[., [[., .], [[., .], .]]], [., .]]
            [[., [[., .], [., [., .]]]], [., .]]
            [., [[[[., .], [., .]], .], [., .]]]
            [., [[[., .], [[., .], .]], [., .]]]
            [., [[[., .], [., [., .]]], [., .]]]
            [., [[., .], [[[., .], .], [., .]]]]
            [., [[., .], [[., [., .]], [., .]]]]
            [., [[., .], [., [[., .], [., .]]]]]
        """
        if self._t1.is_empty():
            yield self._t2
        elif self._t2.is_empty():
            yield self._t1
        else:
            bt1 = self._t1
            bt2 = self._t2
            B = bt1.parent()._element_constructor_
            if hasattr(bt1, "label"):
                for bt in TreeShuffleProduct(bt1, bt2[0]):
                    yield B([bt, bt2[1]], label=bt2.label())
                for bt in TreeShuffleProduct(bt1[1], bt2):
                    yield B([bt1[0], bt], label=bt1.label())
            else:
                for bt in TreeShuffleProduct(bt1, bt2[0]):
                    yield B([bt, bt2[1]])
                for bt in TreeShuffleProduct(bt1[1], bt2):
                    yield B([bt1[0], bt])
