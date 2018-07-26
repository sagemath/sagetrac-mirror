# -*- coding: utf-8 -*-
"""
Patches d'algos classiques sur les arbres.

TRAC : 14498

# TODO:: À ajouter dans le preamble ::

+- Jean-Baptiste Priez (2013):
+      -> several classical algorithms
+      -> some research algorithms from _[LodayRonco]
+
+References
+----------
+
+.. [LodayRonco] Hopf algebra of the planar binary trees,
+    Jean-Louis Loday and
+    Maria O. Ronco
 """
from sage.combinat.cha.patches.monkey_patching import MonkeyPatch
from sage.combinat.binary_tree import BinaryTree

class _(MonkeyPatch, BinaryTree):

    def sylvestrohedron_greater(self):
        '''
        The list of all trees greater than ''self''.
        This is the transitive ideal of its successors.

        The tree
            __o__
           /     \
          o       o
         / \     /
        o   o   o

        has these trees greater than it:

        To          , o        , o        , o        ,  o       ,   o      ,
        | \            \          \          \           \           \
        |  o            o          o           o         _o_        __o__
        |   \            \          \           \       /   \      /     \
        |    o            o          o          _o_    o     o    o       o
        |     \            \        / \        /   \    \     \    \     /
        |      o            o      o   o      o     o    o     o    o   o
        |       \            \          \          /
        |        o            o          o        o
        |         \          /
        L          o        o
           o        ,   o      ,   _o_      ,   _o__     ,   __o__    ,   ___o___  ,
          / \          / \        /   \        /    \       /     \      /       \  
         o   o        o   o      o     o      o     _o_    o       o    o         o 
              \            \          / \          /   \    \       \    \       /  
               o            o        o   o        o     o    o       o    o     o   
                \            \            \            /      \            \
                 o            o            o          o        o            o
                  \          /
                   o        o


             _o_    ,     __o__  T
            /   \        /     \ |
           o     o      o       o|
          / \     \    / \     / |
         o   o     o  o   o   o  |
                                 |
                                 |
                                 |
                                 |
                                 |
                                 J


        TESTS::

            sage: B = BinaryTree
            sage: b = B([None, B([None, B([None, B([])])])]);b
            [., [., [., [., .]]]]
            sage: b.sylvestrohedron_greater()
            [[., [., [., [., .]]]]]
            sage: b = B([B([B([B([]), None]), None]), None]);b
            [[[[., .], .], .], .]
            sage: b.sylvestrohedron_greater()
            [[., [., [., [., .]]]], [., [., [[., .], .]]], [., [[., .], [., .]]], [., [[., [., .]], .]], [., [[[., .], .], .]], [[., .], [., [., .]]], [[., .], [[., .], .]], [[., [., .]], [., .]], [[., [., [., .]]], .], [[., [[., .], .]], .], [[[., .], .], [., .]], [[[., .], [., .]], .], [[[., [., .]], .], .], [[[[., .], .], .], .]]
        '''
        from sage.combinat.tools import transitive_ideal
        return transitive_ideal(lambda x: x.sylvestrohedron_succ(), self)

    def sylvestrohedron_pred(self):
        '''
        Compute the list of predecessor of ''self'' in the
        sylvestrohedron.
        This list is computed by all left rotate possible on
        its nodes.

        For this tree

            __o__
           /     \
          o       o
         / \     /
        o   o   o

        the list is

        T       o ,       _o_  T
        |      /         /   \ |
        |    _o_        o     o|
        |   /   \      /     / |
        |  o     o    o     o  |
        | / \        /         |
        Lo   o      o          J

        TESTS::

            sage: B = BinaryTree
            sage: b = B([B([B([B([]), None]), None]), None]);b
            [[[[., .], .], .], .]
            sage: b.sylvestrohedron_pred()
            []
            sage: b = B([None, B([None, B([None, B([])])])]);b
            [., [., [., [., .]]]]
            sage: b.sylvestrohedron_pred()
            [[[., .], [., [., .]]], [., [[., .], [., .]]], [., [., [[., .], .]]]]
        '''
        res = []
        if self.is_empty():
            return []
        if not self[1].is_empty():
            res.append(self.left_rotate())
        B = self.parent()._element_constructor_
        return (res +
                [B([g, self[1]]) for g in self[0].sylvestrohedron_pred()] +
                [B([self[0], d]) for d in self[1].sylvestrohedron_pred()])

    def sylvestrohedron_smaller(self):
        '''
        The list of all trees smaller than ''self''.
        This is the transitive ideal of its predecessors.

        The tree
            __o__
           /     \
          o       o
         / \     /
        o   o   o

        has these trees smaller than it:

        T    __o__  ,       _o_  ,        o ,         o,         o,           oT
        |   /     \        /   \         /           /          /            / |
        |  o       o      o     o      _o_          o          o            o  |
        | / \     /      /     /      /   \        / \        /            /   |
        |o   o   o      o     o      o     o      o   o      o            o    |
        |              /            / \          /          /            /     |
        |             o            o   o        o          o            o      |
        |                                      /          / \          /       |
        |                                     o          o   o        o        |
        |                                                            /         |
        L                                                           o          J

        TESTS::

            sage: B = BinaryTree
            sage: b = B([None, B([None, B([None, B([])])])]);b
            [., [., [., [., .]]]]
            sage: b.sylvestrohedron_smaller()
            [[., [., [., [., .]]]], [., [., [[., .], .]]], [., [[., .], [., .]]], [., [[., [., .]], .]], [., [[[., .], .], .]], [[., .], [., [., .]]], [[., .], [[., .], .]], [[., [., .]], [., .]], [[., [., [., .]]], .], [[., [[., .], .]], .], [[[., .], .], [., .]], [[[., .], [., .]], .], [[[., [., .]], .], .], [[[[., .], .], .], .]]
            sage: b = B([B([B([B([]), None]), None]), None]);b
            [[[[., .], .], .], .]
            sage: b.sylvestrohedron_smaller()
            [[[[[., .], .], .], .]]
        '''
        from sage.combinat.tools import transitive_ideal
        return transitive_ideal(lambda x: x.sylvestrohedron_pred(), self)

    def sylvestrohedron_succ(self):
        '''
        Compute the list of successors of ''self'' in the sylvestrohedron.
        There is the list of all trees obtains by a right rotate of
        one of its nodes.

        The list of successor of

            __o__
           /     \
          o       o
         / \     /
        o   o   o

        is

        T  _o__     ,   ___o___  ,     _o_    T
        | /    \       /       \      /   \   |
        |o     _o_    o         o    o     o  |
        |     /   \    \       /    / \     \ |
        |    o     o    o     o    o   o     o|
        |         /      \                    |
        L        o        o                   J

        TESTS::

            sage: B = BinaryTree
            sage: b = B([B([B([B([]), None]), None]), None]);b
            [[[[., .], .], .], .]
            sage: b.sylvestrohedron_succ()
            [[[[., .], .], [., .]], [[[., .], [., .]], .], [[[., [., .]], .], .]]
        '''
        res = []
        if self.is_empty():
            return []
        B = self.parent()._element_constructor_
        if not self[0].is_empty():
            res.append(self.right_rotate())
        return (res +
             [B([g, self[1]]) for g in self[0].sylvestrohedron_succ()] +
             [B([self[0], d]) for d in self[1].sylvestrohedron_succ()])

    def q_hook_length_formula(self):
        '''
        Compute the number of permutations which give
        by binary search insertion algorithm the same
        shape tree (`self`).

        .. MATH::

            f_{eq} (T) = \frac{\mid T\mid !}{\prod_{t\in T} \mid t\mid}

        where `\mid T\mid` is the node number of `T` and `t\in T` the set
        of all subtree of `T`.

        There is 20 permutations which give this shape binary tree:

            __o__
           /     \
          o       o
         / \     /
        o   o   o

        by the binary search insertion algorithm.

        TESTS::

            sage: b = BinaryTree([[[],[]],[[],None]]); b
            [[[., .], [., .]], [[., .], .]]
            sage: b.q_hook_length_formula()(q=1)
            20
            sage: BinaryTree([[],[]]).q_hook_length_formula()
            (((q + 2)*q + 2)*q + 1)*q/((q + 1)*q + 1)
        '''
        from sage.combinat.q_analogues import q_factorial, q_int
        from sage.symbolic.ring import SymbolicRing

        q = SymbolicRing().var('q')

        def product_of_subtrees(b):
            if b.is_empty():
                return q ** 0
            return (q ** (-b[0].node_number()) * q_int(b.node_number())) * \
                product_of_subtrees(b[0]) * product_of_subtrees(b[1])

        return q_factorial(self.node_number()) / product_of_subtrees(self)

    def over(self, bt):
        '''
        The ``over`` (`/`) operation defined by Loday-Ronco [LodayRonco]_::

              o       __o__       _o_
             / \  /  /     \  =  /   \
            o   o   o       o   o     o
                     \     /           \
                      o   o           __o__
                                     /     \
                                    o       o
                                     \     /
                                      o   o

        TESTS::

            sage: b1 = BinaryTree([[],[]])
            sage: b2 = BinaryTree([[None,[]],[[],None]])
            sage: b1.over(b2)
            [[., .], [., [[., [., .]], [[., .], .]]]]
        '''
        B = self.parent()._element_constructor_
        if self.is_empty():
            return bt
        lab = None
        if hasattr(self, "label"):
            lab = self.label()
        else:
            return B([self[0], self[1].over(bt)], lab)

    def __div__(self, bt):
        '''
        The ``over`` operation on trees.
        .. see ::method::*over*.

        TESTS::

            sage: b1 = BinaryTree([[],[]])
            sage: b2 = BinaryTree([[None,[]],[[],None]])
            sage: b1/b2
            [[., .], [., [[., [., .]], [[., .], .]]]]
        '''
        return self.over(bt)

    def under(self, bt):
        '''
        The ``under`` (`\`) operation defined by Loday-Ronco [LodayRonco]_::

        TESTS::

            sage: b1 = BinaryTree([[],[[None,[]],None]])
            sage: b2 = BinaryTree([[],[None,[]]])
            sage: b1.under(b2)
            [[[[., .], [[., [., .]], .]], .], [., [., .]]]
        '''
        B = self.parent()._element_constructor_
        if bt.is_empty():
            return self
        lab = None
        if hasattr(bt, "label"):
            lab = bt.label()
        else:
            return B([self.under(bt[0]), bt[1]], lab)

    def _backslash_(self, bt):
        '''
        The ``under`` operation on trees.
        .. see ::method::*under*.

        TESTS::

            sage: b1 = BinaryTree([[],[[None,[]],None]])
            sage: b2 = BinaryTree([[],[None,[]]])
            sage: b1\b2
            [[[[., .], [., .]], [., .]], [[., .], .]]
        '''
        return self.under(bt)
