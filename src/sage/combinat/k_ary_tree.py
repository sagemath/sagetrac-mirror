# -*- coding: utf-8 -*-
"""
K-ary Trees

This module deals with k-ary trees as mathematical (in particular immutable)
objects.

AUTHORS:

- Adrien Boussicault, Bérénice Delcroix-Oger (2015): initial implementation.

"""
#*****************************************************************************
#       Copyright (C) 2015
#           Adrien Boussicault <boussica@labri.fr>,
#           Bérénice Delcroix-Oger <berenice.delcroix@math.univ-toulouse.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.list_clone import ClonableArray
from sage.combinat.abstract_tree import (AbstractClonableTree,
                                         AbstractLabelledClonableTree)
from sage.combinat.ordered_tree import LabelledOrderedTrees
from sage.rings.integer import Integer
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.combinat.combinatorial_map import combinatorial_map
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.misc.cachefunc import cached_method
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.cartesian_product import CartesianProduct
from sage.functions.other import binomial


class KAryTree(AbstractClonableTree, ClonableArray):
    """
    K-ary trees.

    K-ary trees here mean ordered (a.k.a. plane) finite k-ary
    trees, where "ordered" means that the children of each node are
    ordered.

    K-ary trees contain nodes and leaves, where each node has k
    children while each leaf has no children. The number of leaves
    of a k-ary tree always equals k-1 times the number of nodes, plus `1`.

    INPUT:

    - ``children`` -- ``None`` (default) or a list, tuple or iterable of
      length k of k-ary trees or convertible objects. This corresponds
      to the standard recursive definition of a k-ary tree as either a
      leaf or a k-tuple of k-ary trees. Syntactic sugar allows leaving out
      all but the outermost calls of the ``KAryTree()`` constructor, so
      that, e. g., ``KAryTree([KAryTree(None),KAryTree(None)])`` can
      be shortened to ``KAryTree([None,None])``. The arity is given by the
      size of the list given in parameter.
      It is also allowed to abbreviate ``[None, None, ...]`` by ``[]``,
      if the arity is given in parameter or if the arity can be deduced
      by the constructor.

    - ``arity`` -- ``None`` (default) or a positive integer. This corresponds
      to the arity of the tree. If ``None`` is given then the constructor will
      try to deduce the arity from the size of ``children``.

    - ``check`` -- (default: ``True``) whether check for k-arity should be
      performed or not.

    EXAMPLES::

        sage: KAryTree(None)
        .
        sage: KAryTree(None,3)
        .
        sage: KAryTree([None, None, None])
        [., ., .]
        sage: KAryTree([], 0)
        .
        sage: KAryTree([], 1)
        [.]
        sage: KAryTree([], 4)
        [., ., ., .]
        sage: KAryTree([None, None, None], 3)
        [., ., .]
        sage: KAryTree([None, [None, None, None], None])
        [., [., ., .], .]
        sage: KAryTree([None, [], None])
        [., [., ., .], .]
        sage: KAryTree([[None, None], None])
        [[., .], .]
        sage: KAryTree("[[., ., .], ., .]")
        [[., ., .], ., .]
        sage: KAryTree([None, KAryTree([None, None])])
        [., [., .]]
        sage: KAryTree([KAryTree([None])])
        [[.]]

        sage: KAryTree([[None, None], None, None])
        Traceback (most recent call last):
        ...
        ValueError: this is not a 3-ary tree

        sage: KAryTree([KAryTree([None, None, None]), None])
        Traceback (most recent call last):
        ...
        ValueError: this is not a 2-ary tree

    TESTS::

        sage: t1 = KAryTree([[None, [[],[[], None]]],[[],[]]])
        sage: t2 = KAryTree([[[],[]],[]])
        sage: with t1.clone() as t1c:
        ....:     t1c[1,1,1] = t2
        sage: t1 == t1c
        False
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        TODO DOC
        Ensure that k-ary trees created by the enumerated sets and directly
        are the same and that they are instances of :class:`KAryTree`.

        TESTS::

            sage: from sage.combinat.k_ary_tree import KAryTrees_all
            sage: issubclass(KAryTrees_all().element_class, KAryTree)
            True
            sage: t0 = KAryTree([[],[[], None]])
            sage: t0.parent()
            k-ary trees
            sage: type(t0)
            <class 'sage.combinat.k_ary_tree.KAryTrees_all_with_category.element_class'>

            sage: t1 = KAryTrees()([[],[[], None]])
            sage: t1.parent() is t0.parent()
            True
            sage: type(t1) is type(t0)
            True

            sage: t1 = KAryTrees(2,4)([[],[[], None]])
            sage: t1.parent() is t0.parent()
            True
            sage: type(t1) is type(t0)
            True
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        The automatic parent of the elements of this class.

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: KAryTree._auto_parent
            k-ary trees
            sage: KAryTree([None, None]).parent()
            k-ary trees
        """
        return KAryTrees_all()

    def __init__(self, parent, children = None, arity = None, check = True):
        """
        TESTS::

            sage: KAryTree([None, None]).parent()
            k-ary trees
            sage: KAryTree("[., [., [., .]]]")
            [., [., [., .]]]
            sage: KAryTree("[., [., ., .], .]")
            [., [., ., .], .]
            sage: KAryTree("[.,.,.]", 4)
            Traceback (most recent call last):
            ...
            ValueError: this is not a 4-ary tree
            sage: all(
            ....:     KAryTree(repr(tree)) == tree
            ....:     for arity in range(1, 4)
            ....:     for size in range(5)
            ....:     for tree in KAryTrees(arity, size)
            ....: )
            True
        """
        if (isinstance(children, str)):  # if the input is the repr of a binary tree
            children = children.replace(".","None")
            from ast import literal_eval
            children = literal_eval(children)
        if children is None:
            children = []
        elif (children == [] or children == ()) and not arity is None:
            children = [None for i in range(arity)]
        if (children.__class__ is self.__class__ and
            children.parent() == parent):
            children = list(children)
        else:
            children = [self.__class__(parent, x, arity=len( children ) ) for x in children]
        if arity is None:
            self._arity = len( children )
        else:
            self._arity = arity
        ClonableArray.__init__(self, parent, children, check=check)

    def arity(self):
        r"""
        Return the arity of a k-ary tree.

        If the tree is empty, the arity is 0.

        Examples::

            sage: KAryTree("[.]").arity()
            1
            sage: KAryTree("[.,.]").arity()
            2
            sage: KAryTree("[., [., ., .], .]").arity()
            3
            sage: KAryTree(".").arity()
            0
        """
        return len(self)

    def check(self):
        """
        Check that ``self`` is a k-ary tree.

        EXAMPLES::

            sage: KAryTree([[None, None], [None, None]])     # indirect doctest
            [[., .], [., .]]
            sage: KAryTree(None) # indirect doctest
            .
        """
        if self and len(self) != self._arity :
            raise ValueError("this is not a %d-ary tree"%(self._arity))
            for tree in self:
                if tree and tree.arity() != self.arity():
                    raise ValueError("this is not a %d-ary tree"%(self._arity))

    def _repr_(self):
        """
        TESTS::

            sage: t1 = KAryTree([[], None]); t1  # indirect doctest
            [[., .], .]
            sage: KAryTree([[None, t1], None])   # indirect doctest
            [[., [[., .], .]], .]
        """
        if not self:
            return "."
        else:
            return super(KAryTree, self)._repr_()

    def is_empty(self):
        """
        Return whether ``self`` is empty.

        The notion of emptiness employed here is the one which defines
        a binary tree to be empty if its root is a leaf. There is
        precisely one empty k-ary tree.

        EXAMPLES::

            sage: KAryTree().is_empty()
            True
            sage: KAryTree(None).is_empty()
            True
            sage: KAryTree(".").is_empty()
            True
            sage: KAryTree([]).is_empty()
            True
            sage: KAryTree([], 1).is_empty()
            False
            sage: KAryTree([[], None]).is_empty()
            False
        """
        return not self

    def make_leaf(self):
        """
        Modify ``self`` so that it becomes a leaf (i. e., an empty tree).

        .. NOTE:: ``self`` must be in a mutable state.

        .. SEEALSO::

            :meth:`make_node <sage.combinat.binary_tree.BinaryTree.make_node>`

        EXAMPLES::

            sage: t = KAryTree([None, None])
            sage: t.make_leaf()
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with t.clone() as t1:
            ....:     t1.make_leaf()
            sage: t, t1
            ([., .], .)
        """
        self._require_mutable()
        self.__init__(self.parent(), None)

    def comb(self, side=0):
        r"""
        Return the comb of a k-ary trees : there are k combs in a k-ary
        tree (one for every direction). A comb is defined as the list of
        children of the nodes on branches whose direction does not change.

        INPUT:
        - ``side`` -- set to 0 to obtain the leftmost comb, to i+1 to
        obtain the leftmost comb in the set of combs which are on the
        right side of the ith comb and to k-1 to obtain the rightmost comb
        of a k-ary tree.

        OUTPUT:

        A list of `k` tuples of k-ary trees.

        EXAMPLES::

            sage: T = KAryTree([[None,None, None],[[None,None, [None, None, None]], None, None],None])
            sage: T.comb(0)
            [[None, ., .]]
            sage: T.comb(1)
            [[[., ., [., ., .]], None, .]]
            sage: T.comb(2)
            []
            sage: T = KAryTree([[[[None]]]])
            sage: T.comb(0)
            [[None], [None], [None]]


        """
        if self.is_empty():
            return []
        d=self.arity()
        if not side < d :
            raise ValueError("Value %d is a wrong side value : it must be strictly smaller than the arity %d of the tree"%(side,d))
        tree=self[side]
        res=[]
        fc=[]
        while not tree.is_empty():
            for i in range(d):
                if i==side:
                    fc.append(None)
                else:
                    fc.append(tree[i])
            res.append(fc)
            fc=[]
            tree=tree[side]
        return res

    def hook_number(self):
        r"""
        Return the number of hooks in a k-ary trees.

        The hook of a vertex v is the union of {v}, and all the
        branches from {v} in which the direction does not change.

        There is a unique way to partition the vertices in hooks.
        The number of hooks in such the partition is the hook number
        of the tree.

        We can obtain this partition recursively by extracting the root's
        hook and iterating the processus on each tree of the remaining
        forest.

        EXAMPLES::

            sage: T = KAryTree(None)
            sage: T.hook_number()
            0
            sage: T = KAryTree( [None,None,None] )
            sage: T.hook_number()
            1
            sage: T = KAryTree([[None, [None, None]], [[None, None], None]])
            sage: T.hook_number()
            3
            sage: T = KAryTree([None,[[None,None,None],None,[None,None,None]],[None,None,[None,None,None]]] )
            sage: T.hook_number()
            3
        """
        if self.is_empty() or self==None:
            return 0
        s=1
        for i in range(self.arity()):
            for h in self.comb(i):
                if len(h)>0:
                    for el in h:
                        if not(el==None) and not(el.is_empty()):
                            s+=el.hook_number()
        return s

    def twisting_number(self):
        r"""
        Return a k-tuple where the ith element of the tuple is the number
        of straight branches in the k-ary tree tree

        OUTPUT :

        A list of size `k` of non negative integers.

        EXAMPLES::

            sage: T = KAryTree(None)
            sage: T.twisting_number()
            []
            sage: T = KAryTree( [None,None,None] )
            sage: T.twisting_number()
            [0, 0, 0]
            sage: T = KAryTree([[None, [None, None]], [[None, None], None]])
            sage: T.twisting_number()
            [2, 2]
            sage: T = KAryTree([None,[[None,None,None],None,[None,None,None]],[None,None,[None,None,None]]] )
            sage: T.twisting_number()
            [1, 1, 2]
        """
        wn=[]
        d=self.arity()
        for i in range(d):
            wn.append(0)
        if self.is_empty() or self==None:
            return wn
        for i in range(self.arity()):
            if len(self.comb(i))>0:
                wn[i]=wn[i]+1
            for h in self.comb(i):
                for j in range(len(h)):
                    el=h[j]
                    if not(el==None) and not(el.is_empty()):
                        partres=el.twisting_number()
                        for k in range(d):
                            wn[k]=wn[k]+partres[k]
                        if el[j].is_empty:
                            wn[j]=wn[j]+1
        return wn


# Abstract class to serve as a Factory no instance are created.
class KAryTrees(UniqueRepresentation, Parent):
    r"""
    Factory for k-ary trees.

    INPUT:

    - ``arity`` -- (optional) an integer
    - ``size`` -- (optional) an integer

    OUPUT:

    - the set of all k-ary trees (of the given ``arity`` and ``size`` if
      specified )

    EXAMPLES::

        sage: KAryTrees()
        k-ary trees

        sage: KAryTrees(3)
        3-ary trees

        sage: KAryTrees(3, 2)
        3-ary trees of size 2

    .. NOTE:: this is a factory class whose constructor returns instances of
              subclasses.

    .. NOTE:: the fact that KAryTrees is a class instead of a simple callable
              is an implementation detail. It could be changed in the future
              and one should not rely on it.
    """
    @staticmethod
    def __classcall_private__(cls, arity=None, size=None):
        """
        TESTS::

            sage: from sage.combinat.k_ary_tree import KAryTrees_all, KAryTrees_size, KAryTrees_arity
            sage: isinstance(KAryTrees(3, 2), KAryTrees)
            True
            sage: isinstance(KAryTrees(3), KAryTrees)
            True
            sage: isinstance(KAryTrees(), KAryTrees)
            True
            sage: KAryTrees(3, 2) is KAryTrees_size(3,2)
            True
            sage: KAryTrees(3) is KAryTrees_arity(3)
            True
            sage: KAryTrees(2, 5).cardinality()
            42
            sage: KAryTrees() is KAryTrees_all()
            True
        """
        if size is None and arity is None:
            return KAryTrees_all()
        else:
            if not (isinstance(arity, (Integer, int)) and arity >= 0):
                raise ValueError("arity must be a positive integer")
            if size is None:
                return KAryTrees_arity(Integer(arity))
            else:
                if not (isinstance(size, (Integer, int)) and size >= 0):
                    raise ValueError("size must be a nonnegative integer")
                return KAryTrees_size(Integer(arity), Integer(size))

    @cached_method
    def leaf(self):
        """
        Return a leaf tree with ``self`` as parent.

        EXAMPLES::

            sage: KAryTrees().leaf()
            .

        TEST::

            sage: (KAryTrees().leaf() is
            ....:  sage.combinat.k_ary_tree.KAryTrees_all().leaf())
            True
        """
        return self(None)


from sage.structure.parent import Parent
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer import Integer
from sage.structure.list_clone import ClonableArray


class IntegerPair(ClonableArray):
    def __init__(self, parent, v, check = True):
        ClonableArray.__init__(self, parent, v, check=check)
    def check(self):
        if not (not self or len(self) == 2):
            raise ValueError("this is not a pair o integer")
    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)
    @lazy_class_attribute
    def _auto_parent(cls):
        return IntegerPairs()


class IntegerPairs(UniqueRepresentation, Parent):
    def __init__(self):
        Parent.__init__(self, category = InfiniteEnumeratedSets())

    def _repr_(self):
        return "Pair of non negative integers"

    def __contains__(self, elt):
        return (
            ( Integer(elt[0]) == Integer(0) ) and
            ( Integer(elt[1]) == Integer(0) )
        ) or (
            Integer(elt[0]) >= Integer(1) and Integer(elt[0]) >= Integer(1)
        )

    def __iter__(self):
        yield self._element_constructor_([0, 0])
        n = 2
        while True:
            for i in range(1, n):
                yield self._element_constructor_([i, n-i])
            n += 1

    def __call__(self, elt):
        if elt in self:
            return self._element_constructor_( elt )
        else:
            raise ValueError("Value %s is not a pair of non negative integer."%(elt))

    def an_element(self):
        return self._element_constructor_([Integer(42),Integer(3)])

    def next(self, v):
        if v[1] == 1:
            return self._element_constructor_( [1, v[0]+1] )
        return self._element_constructor_( [v[0]+1, v[0]-1] )

    def _element_constructor_(self, *args, **keywords):
        return self.element_class(self,*args, **keywords)

    Element = IntegerPair




#################################################################
# Enumerated set of all kary trees
#################################################################
class KAryTrees_all(DisjointUnionEnumeratedSets, KAryTrees):

    def __init__(self):
        """
        TODO

        TESTS::

            sage: from sage.combinat.k_ary_tree import KAryTrees_all
            sage: K = KAryTrees_all()
            sage: K.cardinality()
            +Infinity

            sage: it = iter(K)
            sage: (next(it), next(it), next(it), next(it), next(it))
            (., [.], [[.]], [., .], [[[.]]])
            sage: next(it).parent()
            k-ary trees
            sage: K([None, None, None])
            [., ., .]

            TOTO : est-ce que l'on permet : K([]) et K(None) qui seraient égal
            respectivement à : K([None, None, None]) et K(None, 3).

            sage: K is KAryTrees_all()
            True

            #sage: TestSuite(K).run() # long time
            """

        DisjointUnionEnumeratedSets.__init__(
            self, Family(
                IntegerPairs(),
                lambda x : KAryTrees_size(arity=x[0], size=x[1])
            ),
            facade=True, keepkey = False
        )

    def _repr_(self):
        """
        TEST::

            sage: KAryTrees()   # indirect doctest
            k-ary trees
        """
        return "k-ary trees"

    def __contains__(self, x):
        """
        TESTS::

            sage: K = KAryTrees()
            sage: 1 in K
            False
            sage: K([None, None, None]) in K
            True
        """
        return isinstance(x, self.element_class)

    def __call__(self, x=None, *args, **keywords):
        """
        Ensure that ``None`` instead of ``0`` is passed by default.

        TESTS::

            sage: K = KAryTrees()
            sage: K([None, None])
            [., .]
        """
        return super(KAryTrees, self).__call__(x, *args, **keywords)

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        EXAMPLES::

            sage: KAryTrees().unlabelled_trees()
            k-ary trees
        """
        return self

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        EXAMPLES::

            sage: KAryTrees().labelled_trees()
            Labelled k-ary trees
        """
        return LabelledKAryTrees()

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: K = KAryTrees()
            sage: K._element_constructor_([None, None, None])
            [., ., .]
            sage: K([[None, None, None],None, None])
            [[., ., .], ., .]
            sage: K(None)
            .
        """
        return self.element_class(self, *args, **keywords)

    Element = KAryTree


#################################################################
# Enumerated set of all kary trees with a fixed arity
#################################################################
class KAryTrees_arity(DisjointUnionEnumeratedSets, KAryTrees):

    def __init__(self, arity):
        """
        TODO

        TESTS::

            sage: from sage.combinat.k_ary_tree import KAryTrees_all, KAryTrees_arity
            sage: K = KAryTrees_arity(3)
            sage: K.cardinality()
            +Infinity

            sage: it = iter(K)
            sage: (next(it), next(it), next(it), next(it), next(it))
            (., [., ., .], [[., ., .], ., .], [., [., ., .], .], [., ., [., ., .]])
            sage: next(it).parent()
            k-ary trees
            sage: K([None, None, None])
            [., ., .]

            TOTO : est-ce que l'on permet : K([]) et K(None) qui seraient égal
            respectivement à : K([None, None, None]) et K(None, 3).

            sage: K is KAryTrees_arity(3)
            True

            #sage: TestSuite(K).run() # long time
            """
        self._arity = arity
        DisjointUnionEnumeratedSets.__init__(
            self, Family(
                NonNegativeIntegers(), lambda x: KAryTrees_size(self._arity, x)
            ), facade=True, keepkey = False
        )

    def _repr_(self):
        """
        TEST::

            sage: KAryTrees(3)   # indirect doctest
            3-ary trees
        """
        return "%d-ary trees"%(self._arity)

    def __contains__(self, x):
        """
        TESTS::

            sage: from sage.combinat.k_ary_tree import KAryTrees_arity
            sage: K = KAryTrees_arity(3)
            sage: 1 in K
            False
            sage: K([None, None, None]) in K
            True
            sage: K1 = KAryTrees()
            sage: K1([None, None]) in K
            False
        """
        return isinstance(x, self.element_class) and (
            x.is_empty() or x.arity() == self._arity
        )

    def __call__(self, x=None, *args, **keywords):
        """
        Ensure that ``None`` instead of ``0`` is passed by default.

        TODO

        TESTS::

            sage: from sage.combinat.k_ary_tree import KAryTrees_arity
            sage: K = KAryTrees_arity(4)
            sage: K([None, None, None, None])
            [., ., ., .]
            sage: K([])
            [., ., ., .]
        """
        return super(KAryTrees, self).__call__(
            x, arity=self._arity, *args, **keywords
        )

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        EXAMPLES::

            sage: from sage.combinat.k_ary_tree import KAryTrees_arity
            sage: KAryTrees_arity(3).unlabelled_trees()
            3-ary trees
        """
        return self

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        EXAMPLES::

            sage: from sage.combinat.k_ary_tree import KAryTrees_arity

            TODO
            #sage: KAryTrees_arity(3).labelled_trees()
            #Labelled k-ary trees
        """
        return LabelledKAryTrees(self._arity)

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: from sage.combinat.k_ary_tree import KAryTrees_arity
            sage: K = KAryTrees_arity(3)
            sage: K._element_constructor_([None, None, None])
            [., ., .]
            sage: K([[None, None, None],None, None]) # indirect doctest
            [[., ., .], ., .]
            sage: K(None)    # indirect doctest
            .
        """
        return self.element_class(self, *args, **keywords)

    Element = KAryTree



from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
#################################################################
# Enumerated set of k-ary trees of a given size
#################################################################
class KAryTrees_size(KAryTrees):
    """
    The enumerated sets of k-ary trees of given size

    TESTS::

    TODO

        #sage: from sage.combinat.k_ary_tree import KAryTrees_size
        #sage: for i in range(6): TestSuite(KAryTrees_size(3,i)).run()
    """
    def __init__(self, arity, size):
        """
        TESTS::

            sage: S = KAryTrees(2, 3)
            sage: S == loads(dumps(S))
            True

            sage: S is KAryTrees(2, 3)
            True
        """
        super(KAryTrees_size, self).__init__(
            category = FiniteEnumeratedSets()
        )
        self._arity = arity
        self._size = size

    def _repr_(self):
        """
        TESTS::

            sage: KAryTrees(2, 3)   # indirect doctest
            2-ary trees of size 3
        """
        return "%d-ary trees of size %s" % (self._arity, self._size)

    def __contains__(self, x):
        """
        TESTS::

            sage: K = KAryTrees(2, 3)
            sage: 1 in K
            False
            sage: KS = KAryTrees()
            sage: KS([[],[]]) in K
            True
            sage: KS([None, []]) in K
            False
            sage: KS([None, [], []]) in K
            False
        """
        return isinstance(x, self.element_class) and (
            x.is_empty() or (
                x.node_number() == self._size and x.arity() == self._arity
            )
        )

    def _an_element_(self):
        """
        TESTS::

            sage: KAryTrees(2, 0).an_element()  # indirect doctest
            .
        """
        return self.first()

    def cardinality(self):
        """
        The cardinality of ``self``

        TODO
        This is a ?? number.

        TESTS::

            TODO
            sage: KAryTrees(2, 0).cardinality()
            1
            sage: KAryTrees(2, 5).cardinality()
            42
        """
        return binomial(self._arity*self._size, self._size)/(
            (self._arity-1)*self._size + 1
        )

    def __iter__(self):
        """
        A basic generator.

        .. TODO:: could be optimized.

        TESTS::

            sage: KAryTrees(2, 0).list()
            [.]
            sage: KAryTrees(2, 1).list()
            [[., .]]
            sage: KAryTrees(2, 3).list()
            [[[[., .], .], .],
             [[., [., .]], .],
             [[., .], [., .]],
             [., [[., .], .]],
             [., [., [., .]]]]
        """
        if self._size == 0:
            yield self._element_constructor_(None)
        else:
            for v in IntegerVectors( self._size-1, length=self._arity ):
                cp = CartesianProduct(
                    * map( lambda x:self.__class__(self._arity, x), v )
                )
                for children in cp:
                    yield self._element_constructor_( children )

    @lazy_attribute
    def _parent_for(self):
        """
        The parent of the elements generated by ``self``.

        TESTS::

            sage: S = KAryTrees(2, 3)
            sage: S._parent_for
            k-ary trees
        """
        return KAryTrees_all()

    @lazy_attribute
    def element_class(self):
        """
        TESTS::

            sage: K = KAryTrees(2, 3)
            sage: K.element_class
            <class 'sage.combinat.k_ary_tree.KAryTrees_all_with_category.element_class'>
            sage: K.first().__class__ == KAryTrees().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: S = KAryTrees(2, 0)
            sage: S([None])   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: this is not a 2-ary tree
            sage: S(None)   # indirect doctest
            .

            sage: S = KAryTrees(2, 1)   # indirect doctest
            sage: S([None, None])
            [., .]
            sage: S([None, []])
            Traceback (most recent call last):
            ...
            ValueError: wrong number of nodes
        """
        res = self.element_class(
            self._parent_for, arity=self._arity, *args, **keywords
        )
        if res.node_number() != self._size:
            raise ValueError("wrong number of nodes")
        return res

class LabelledKAryTree(AbstractLabelledClonableTree, KAryTree):
    """
    Labelled k-ary trees.

    A labelled k-ary tree is a k-ary tree (see :class:`KAryTree` for
    the meaning of this) with a label assigned to each node.
    The labels need not be integers, nor are they required to be distinct.
    ``None`` can be used as a label.

    .. WARNING::

        While it is possible to assign values to leaves (not just nodes)
        using this class, these labels are disregarded by various
        methods such as
        :meth:`~sage.combinat.abstract_tree.AbstractLabelledTree.labels`,
        :meth:`~sage.combinat.abstract_tree.AbstractLabelledTree.map_labels`,
        and (ironically)
        :meth:`~sage.combinat.abstract_tree.AbstractLabelledTree.leaf_labels`.

    INPUT:

    - ``children`` -- ``None`` (default) or a list, tuple or iterable of
      length `k` of labelled `k`-ary trees or convertible objects. This
      corresponds to the standard recursive definition of a labelled
      k-ary tree as being either a leaf, or a pair of:

      - a k-tuple of labelled binary trees,
      - and a label.

      (The label is specified in the keyword variable ``label``; see
      below.)

      Syntactic sugar allows leaving out all but the outermost calls
      of the ``LabelledKAryTree()`` constructor, so that, e. g.,
      ``LabelledKAryTree([LabelledKAryTree(None),LabelledKAryTree(None)])``
      can be shortened to ``LabelledKAryTree([None,None])``. However,
      using this shorthand, it is impossible to label any vertex of
      the tree other than the root (because there is no way to pass a
      ``label`` variable without calling ``LabelledKAryTree``
      explicitly).

      It is also allowed to abbreviate ``[None, ...]`` by ``[]`` by using
      the arity parameter if one does not want to label the leaves
      (which one should not do anyway!).

    - `̀̀̀`arity`` -- ``None`` (default) or a positive integer. This corresponds
      to the arity of the tree. If ``None`` is given then the constructor will
      try to deduce the arity from the size of ``children``.

    - ``label`` -- (default: ``None``) the label to be put on the root
      of this tree.

    - ``check`` -- (default: ``True``) whether checks should be
      performed or not.

    .. TODO::

        It is currently not possible to use ``LabelledKAryTree()``
        as a shorthand for ``LabelledKAryTree(None)`` (in analogy to
        similar syntax in the ``KAryTree`` class).

    EXAMPLES::

        sage: LabelledKAryTree(None)
        .
        sage: LabelledKAryTree(None, label="ae")    # not well supported
        'ae'
        sage: LabelledKAryTree([])
        .
        sage: LabelledKAryTree([], arity=2, label=3)    # not well supported
        3[., .]
        sage: LabelledKAryTree([None, None])
        None[., .]
        sage: LabelledKAryTree([None, None], label=5)
        5[., .]
        sage: LabelledKAryTree([None, []])
        None[., None[., .]]
        sage: LabelledKAryTree([None, [], None], label=4)
        4[., None[., ., .], .]
        sage: LabelledKAryTree([[], None])
        None[None[., .], .]
        sage: LabelledKAryTree("[[], .]", label=False)
        False[None[., .], .]
        sage: LabelledKAryTree([None, LabelledKAryTree([None, None], label=4)], label=3)
        3[., 4[., .]]
        sage: LabelledKAryTree([None, KAryTree([None, None])], label=3)
        3[., None[., .]]

        sage: LabelledKAryTree([[None, None], None, []])
        Traceback (most recent call last):
        ...
        ValueError: this is not a 3-ary tree

        sage: LBT = LabelledKAryTree
        sage: t1 = LBT([[LBT([], arity=2, label=2), None], None], label=4); t1
        4[None[2[., .], .], .]

    TESTS::

        sage: t1 = LabelledKAryTree([[None, [[],[[], None]]],[[],[]]])
        sage: t2 = LabelledKAryTree([[[],[]],[]])
        sage: with t1.clone() as t1c:
        ....:     t1c[1,1,1] = t2
        sage: t1 == t1c
        False

    We check for :trac:`16314`::

        sage: t1 = LBT([ LBT([LBT([], arity=2, label=2),
        ....:                 LBT([], arity=2, label=5)], label=6),
        ....:            None], label=4); t1
        4[6[2[., .], 5[., .]], .]
        sage: class Foo(LabelledKAryTree):
        ....:     pass
        sage: t2 = Foo(t1.parent(), t1); t2
        4[6[2[., .], 5[., .]], .]
        sage: t2.label()
        4
        sage: t2[0].label()
        6
        sage: t2.__class__, t2[0].__class__
        (<class '__main__.Foo'>, <class '__main__.Foo'>)
    """
    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        Ensure that trees created by the sets and directly are the same and
        that they are instances of :class:`LabelledTree`.

        TESTS::

            sage: issubclass(LabelledKAryTrees().element_class, LabelledKAryTree)
            True
            sage: t0 = LabelledKAryTree([[],[[], None]], label = 3)
            sage: t0.parent()
            Labelled k-ary trees
            sage: type(t0)
            <class 'sage.combinat.k_ary_tree.LabelledKAryTrees_with_category.element_class'>
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)


    @lazy_class_attribute
    def _auto_parent(cls):
        """
        The automatic parent of the elements of this class.

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: LabelledKAryTree._auto_parent
            Labelled k-ary trees
            sage: LabelledKAryTree([], arity=3, label = 3).parent()
            Labelled k-ary trees
        """
        return LabelledKAryTrees()

    def _repr_(self):
        """
        TESTS::

            sage: LBT = LabelledKAryTree
            sage: t1 = LBT([[LBT([], arity=2, label=2), None], None], label=4); t1
            4[None[2[., .], .], .]
            sage: LBT([[],[[], None]], label = 3)   # indirect doctest
            3[None[., .], None[None[., .], .]]
        """
        if not self:
            if self._label is not None:
                return repr(self._label)
            else:
                return "."
        else:
            return "%s%s" % (self._label, self[:])

    _UnLabelled = KAryTree


class LabelledKAryTrees(LabelledOrderedTrees):
    """
    This is a parent stub to serve as a factory class for trees with various
    labels constraints.
    """
    def _repr_(self):
        """
        TESTS::

            sage: LabelledKAryTrees()   # indirect doctest
            Labelled k-ary trees
        """
        return "Labelled k-ary trees"

    def _an_element_(self):
        """
        Return a labelled k-ary tree.

        EXAMPLE::

            sage: LabelledKAryTrees().an_element()   # indirect doctest
            toto[42[3[., .], 3[., .]], 5[None[., .], None[., .]]]
        """
        LT = self._element_constructor_
        t  = LT([None, None], label = 3)
        t1 = LT([t,t], label = 42)
        t2  = LT([[None, None], [None, None]], label = 5)
        return LT([t1,t2], label = "toto")

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        EXAMPLES::

            sage: LabelledKAryTrees().unlabelled_trees()
            k-ary trees

        This is used to compute the shape::

            sage: t = LabelledKAryTrees().an_element().shape(); t
            [[[., .], [., .]], [[., .], [., .]]]
            sage: t.parent()
            k-ary trees

        TESTS::

            TODO

            #sage: t = LabelledKAryTrees().an_element()
            #sage: t.canonical_labelling()
            #4[2[1[., .], 3[., .]], 6[5[., .], 7[., .]]]
        """
        return KAryTrees_all()

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        EXAMPLES::

            sage: LabelledKAryTrees().labelled_trees()
            Labelled k-ary trees
        """
        return self

    Element = LabelledKAryTree
