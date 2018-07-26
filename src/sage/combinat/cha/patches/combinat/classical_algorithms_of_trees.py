"""
This patch implement a several classical algorithms on trees

TRAC : 14498

"""
from sage.combinat.cha.patches.monkey_patching import MonkeyPatch

from sage.combinat.abstract_tree import AbstractTree
from sage.combinat.binary_tree import BinaryTree, LabelledBinaryTree


class _(MonkeyPatch, AbstractTree):

    def pre_order_traversal(self, action=lambda node: None):
        '''
        The depth first pre-order traversal algorithm.

        For example on the following binary tree `b`::

              ___3____
             /        \
            1         _7_
             \       /   \
              2     5     8
                   / \
                  4   6

        the ``depth first pre-order traversal algorithm`` explores `b` in the
        following order of nodes `3,1,2,7,5,4,6,8`.

        The algorithm is::

            manipulate the root
            then explore each subtrees (by the algorithm)


        An other example::

                __1____
               /  /   /
              2  6   8_
              |  |  / /
              3_ 7 9 10
             / /
            4 5

        The algorithm explores this tree in the following order:
        `1,2,3,4,5,6,7,8,9,10`.

        TESTS::

            sage: l = []
            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).\
            ....:     canonical_labelling(); b
            3[1[., 2[., .]], 7[5[4[., .], 6[., .]], 8[., .]]]
            sage: b.pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [3, 1, 2, 7, 5, 4, 6, 8]
            sage: t = OrderedTree([[[[],[]]],[[]],[[],[]]]).\
            ....:     canonical_labelling(); t
            1[2[3[4[], 5[]]], 6[7[]], 8[9[], 10[]]]
            sage: l = []
            sage: t.pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            sage: l = []
            sage: BinaryTree().canonical_labelling().\
            ....:    pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            []
            sage: OrderedTree([]).canonical_labelling().\
            ....:    pre_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [1]
        '''
        if self.is_empty():
            return
        stack = []
        stack.append(self)
        while len(stack) > 0:
            node = stack.pop()
            action(node)
            for i in range(len(node)):
                subtree = node[-i - 1]
                if not subtree.is_empty():
                    stack.append(subtree)

    def post_order_traversal(self, action=lambda node: None):
        '''
        The depth first post-order traversal algorithm.

        For example on the following binary tree `b`::

              ___3____
             /        \
            1         _7_
             \       /   \
              2     5     8
                   / \
                  4   6

        the ``depth first post-order traversal algorithm`` explores `b` in
        the following order of nodes `2,1,4,6,5,8,7,3`.

        The algorithm is::

            explore each subtrees (by the algorithm)
            then manipulate the root

        An other example::

                __1____
               /  /   /
              2  6   8_
              |  |  / /
              3_ 7 9 10
             / /
            4 5

        The algorithm explores this tree in the following order:
        `4,5,3,2,7,6,9,10,8,1`.

        TESTS::

            sage: l = []
            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).\
            ....:    canonical_labelling(); b
            3[1[., 2[., .]], 7[5[4[., .], 6[., .]], 8[., .]]]
            sage: b.post_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [2, 1, 4, 6, 5, 8, 7, 3]
            sage: t = OrderedTree([[[[],[]]],[[]],[[],[]]]).\
            ....:    canonical_labelling(); t
            1[2[3[4[], 5[]]], 6[7[]], 8[9[], 10[]]]
            sage: l = []
            sage: t.post_order_traversal(lambda node: l.append(node.label()))
            sage: l
            [4, 5, 3, 2, 7, 6, 9, 10, 8, 1]
            sage: l = []
            sage: BinaryTree().canonical_labelling().\
            ....:    post_order_traversal(lambda node: l.append(
            ....:            node.label()
            ....: ))
            sage: l
            []
            sage: OrderedTree([]).canonical_labelling().post_order_traversal(
            ....:    lambda node: l.append(node.label())
            ....: )
            sage: l
            [1]
        '''
        if self.is_empty():
            return
        for subtree in self:
            subtree.post_order_traversal(action)
        action(self)

    def breadth_first_order_traversal(self, action=lambda node: None):
        '''
        The breadth first order traversal algorithm.

        For example on the following binary tree `b`::

              ___3____
             /        \
            1         _7_
             \       /   \
              2     5     8
                   / \
                  4   6

        the ``breadth first order traversal algorithm`` explores `b` in the
        following order of nodes `3,1,7,2,5,8,4,6`.

        The algorithm is::

            queue <- ( root )
            while the queue is not empty:
                node <- pop( queue )
                manipulate the node
                append in the queue all subtrees of the node

        TESTS::

            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).\
            ....:    canonical_labelling()
            sage: l = []
            sage: b.breadth_first_order_traversal(
            ....:    lambda node: l.append(node.label())
            ....: ); l
            [3, 1, 7, 2, 5, 8, 4, 6]
            sage: t = OrderedTree([[[[],[]]],[[]],[[],[]]]).\
            ....:    canonical_labelling(); t
            1[2[3[4[], 5[]]], 6[7[]], 8[9[], 10[]]]
            sage: l = []
            sage: t.breadth_first_order_traversal(
            ....:    lambda node: l.append(node.label())
            ....: ); l
            [1, 2, 6, 8, 3, 7, 9, 10, 4, 5]
            sage: l = []
            sage: BinaryTree().canonical_labelling().\
            ....:     breadth_first_order_traversal(
            ....:        lambda node: l.append(node.label())
            ....:     ); l
            []
            sage: OrderedTree([]).canonical_labelling().\
            ....:     breadth_first_order_traversal(
            ....:         lambda node: l.append(node.label())
            ....:     ); l
            [1]
        '''
        if self.is_empty():
            return
        queue = []
        queue.append(self)
        while len(queue) > 0:
            node = queue.pop()
            action(node)
            for subtree in node:
                if not subtree.is_empty():
                    queue.insert(0, subtree)


class _(MonkeyPatch, BinaryTree):

    def in_order_traversal(self,
        node_action=lambda node: None,
        leaf_action=lambda leaf: None
    ):
        '''
        The depth first infix-order traversal algorithm.

        For example on the following binary tree `b` where we denote leafs by
        `a,b,c...` and nodes by `1,2,3...`::

              ____3____
             /         \
            1          __7__
           / \        /     \
          a   2      _5_     8
             / \    /   \   / \
            b   c  4     6 h   i
                  / \   / \
                 d   e f   g

        the ``depth first infixe-order traversal algorithm`` explores `b` in
        the following order of nodes `a,1,b,2,c,3,d,4,e,5,f,6,g,7,h,8,i`.

        The algorithm is::

            explore the left subtree
            manipulate the root
            explore the right subtree

        TESTS::

            sage: nb_leaf = 0
            sage: def l_action(_):
            ....:    global nb_leaf
            ....:    nb_leaf += 1
            sage: nb_node = 0
            sage: def n_action(_):
            ....:    global nb_node
            ....:    nb_node += 1
            sage: BinaryTree().in_order_traversal(n_action, l_action)
            sage: nb_leaf, nb_node
            (1, 0)
            sage: nb_leaf, nb_node = 0, 0
            sage: b = BinaryTree([[],[[],[]]]); b
            [[., .], [[., .], [., .]]]
            sage: b.in_order_traversal(n_action, l_action)
            sage: nb_leaf, nb_node
            (6, 5)
            sage: nb_leaf, nb_node = 0, 0
            sage: b = b.canonical_labelling()
            sage: b.in_order_traversal(n_action, l_action)
            sage: nb_leaf, nb_node
            (6, 5)
            sage: l = []
            sage: b.in_order_traversal(lambda node: l.append( node.label() ))
            [1, 2, 3, 4, 5]
            sage: leaf = 'a'
            sage: l = []
            sage: def l_action(_):
            ....:    global leaf, l
            ....:    l.append(leaf)
            ....:    leaf = chr( ord(leaf)+1 )
            sage: n_action = lambda node: l.append( node.label() )
            sage: b = BinaryTree([[None,[]],[[[],[]],[]]]).\
            ....:    canonical_labelling()
            sage: b.in_order_traversal(n_action, l_action)
            sage: l
            ['a', 1, 'b', 2, 'c', 3, 'd', 4, 'e', 5, 'f', 6, 'g', 7, 'h', 8, 'i']
        '''
        if self.is_empty():
            leaf_action(self)
            return
        self[0].in_order_traversal(node_action, leaf_action)
        node_action(self)
        self[1].in_order_traversal(node_action, leaf_action)

    def canonical_permutation(self, left_to_right=True):
        '''
        Compute the canonical permutation associated to the binary search tree
        insertion from the `right` or the `left` of the permutation.

        TESTS::

            sage: b = BinaryTree([[[],[]],[[],None]])
            sage: b.canonical_permutation(False)
            [1, 3, 2, 5, 6, 4]
            sage: b.canonical_permutation(True)
            [4, 2, 1, 3, 6, 5]
            sage: b.canonical_permutation()
            [4, 2, 1, 3, 6, 5]
            sage: b.canonical_permutation().binary_search_tree().shape() == b
            True
            sage: b.canonical_permutation(False).\
            ....:    binary_search_tree(False).shape() == b
            True
            sage: b.canonical_permutation(False).\
            ....:    binary_search_tree().shape() == b
            False
            sage: b.canonical_permutation().\
            ....:    binary_search_tree(False).shape() == b
            False
        '''
        from sage.combinat.permutation import Permutation
        assert(isinstance(left_to_right, bool))
        l = []
        if left_to_right:
            self.canonical_labelling().pre_order_traversal(
                lambda node: l.append(node.label())
            )
        else:
            self.canonical_labelling().post_order_traversal(
                lambda node: l.append(node.label())
            )
        return Permutation(l)

    def right_rotate(self):
        '''
        Right rotation operation of tree:

            o                     _o_
           /                     /   \
          o    -right-rotate->  o     o
         / \                         /
        o   o                       o

              __o__                         _o__
             /     \                       /    \
            o       o  -right-rotate->    o     _o_
           / \                           /     /   \
          o   o                         o     o     o
         /     \                               \
        o       o                               o

        TESTS::

            sage: b = BinaryTree([[[],[]], None]); b
            [[[., .], [., .]], .]
            sage: b.right_rotate()
            [[., .], [[., .], .]]
            sage: b = BinaryTree([[[[],None],[None,[]]], []]);b
            [[[., .], .], [., [., .]]], [., .]]
            sage: b.right_rotate()
            [[[., .], .], [[., [., .]], [., .]]]
        '''
        B = self.parent()._element_constructor_
        return B([self[0][0], B([self[0][1], self[1]])])

    def left_rotate(self):
        '''
        Right rotation operation of tree:

          _o_                        o
         /   \                      /
        o     o  -left-rotate->    o
             /                    / \
            o                    o   o

              __o__                            o
             /     \                          /
            o       o  -left-rotate->        o
           / \                              /
          o   o                            o
         /     \                          / \
        o       o                        o   o
                                        /     \
                                       o       o


        TESTS::

            sage: b = BinaryTree([[],[[],None]]); b
            [[., .], [[., .], .]]
            sage: b.left_rotate()
            [[[., .], [., .]], .]
            sage: b.left_rotate().right_rotate() == b
            True
        '''
        B = self.parent()._element_constructor_
        return B([B([self[0], self[1][0]]), self[1][1]])


class _(MonkeyPatch, LabelledBinaryTree):

    def right_rotate(self):
        '''
        Right rotation operation of tree:

                y                      x
               / \                    / \
              x   C -right-rotate->  A   y
             / \                        / \
            A   B                      B   C

        TESTS::

            sage: LB = LabelledBinaryTree
            sage: b = LB([LB([LB([],"A"), LB([],"B")],"x"),LB([],"C")], "y"); b
            y[x[A[., .], B[., .]], C[., .]]
            sage: b.right_rotate()
            x[A[., .], y[B[., .], C[., .]]]
        '''
        B = self.parent()._element_constructor_
        return B([
            self[0][0],
            B([self[0][1], self[1]], self.label())
        ], self[0].label())

    def left_rotate(self):
        '''
        Left rotation operation of tree:

                y                    x
               / \                  / \
              x   C <-left-rotate- A   y
             / \                      / \
            A   B                    B   C


        TESTS::

            sage: LB = LabelledBinaryTree
            sage: b = LB([LB([LB([],"A"), LB([],"B")],"x"),LB([],"C")], "y"); b
            y[x[A[., .], B[., .]], C[., .]]
            sage: b == b.right_rotate().left_rotate()
            True
        '''
        B = self.parent()._element_constructor_
        return B([
            B([self[0], self[1][0]], self.label()),
            self[1][1]
        ], self[1].label())
