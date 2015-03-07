from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.ordered_tree import LabelledOrderedTrees


class FreeShuffleOperad(CombinatorialFreeModule):
    r"""
    The free shuffle operad over any given set of generators

    Generators are given as a list of labelled ordered trees
    """

    def __init__(self, R, gens, order='pathlex'):
        """
        EXAMPLES::

            sage: from sage.operads.free_shuffle_operad import FreeShuffleOperad
            sage: L = LabelledOrderedTrees()
            sage: lf = lambda A: L([],label=A) # leaf
            sage: g = L([lf(1),lf(2)],label="a")
            sage: A = FreeShuffleOperad(QQ,(g,)); A
            The Free shuffle operad over Rational Field with generators (a[1[], 2[]],)
            sage: TestSuite(A).run()
        """
        CombinatorialFreeModule.__init__(self, R, LabelledOrderedTrees(),
                                         monomial_cmp=self.cmp_fun(order))
        self._gens = gens
        self._term_order = order

    def cmp_fun(self, order):
        """
        ordering according to a term order

        possible values so far:

        - 'pathlex'
        """
        if order == 'pathlex':
            return lambda x, y: self.pathlex_cmp(x, y)
        return lambda x, y: cmp(x, y)

    def pathlex_cmp(self, x, y):
        """
        term order according to Dotsenko-Koroshkin

        first one compares the number of leaves

        then one compares the words one by one in deglex order

        then one compares the permutations in revlex order

        EXAMPLES::

            sage: from sage.operads.free_shuffle_operad import FreeShuffleOperad
            sage: L = LabelledOrderedTrees()
            sage: lf = lambda A: L([],label=A) # leaf
            sage: g = L([lf(1),lf(2)],label="a")
            sage: A = FreeShuffleOperad(QQ,(g,))

            sage: u = L([L([lf(1),lf(3)],label="a"),lf(2)],label="b")
            sage: u.pathlex_key()
            (['ba', 'b', 'ba'], [1, 3, 2])

            sage: v = L([L([lf(1),lf(2)],label="a"),lf(3)],label="b")
            sage: v.pathlex_key()
            (['ba', 'ba', 'b'], [1, 2, 3])

            sage: w = L([lf(1),L([lf(2),lf(3)],label="a")],label="b")
            sage: w.pathlex_key()
            (['b', 'ba', 'ba'], [1, 2, 3])

            sage: A.pathlex_cmp(u,v)
            -1
        """
        lx, px = x.pathlex_key()
        ly, py = y.pathlex_key()
        if len(px) != len(py):
            return len(px) < len(py)  # ok
        elif lx != ly:
            dlx = []
            for u in lx:
                dlx += [len(u), u]
            dly = []
            for v in ly:
                dly += [len(v), v]
            return cmp(dlx, dly)  # ok ? degree-lex ?
        else:
            return cmp(px.complement(), py.complement())  # ok reverse-lex

    def gens(self):
        """
        Return the generators of the given free shuffle operad

        EXAMPLES::

            sage: from sage.operads.free_shuffle_operad import FreeShuffleOperad
            sage: L = LabelledOrderedTrees()
            sage: lf = lambda A: L([],label=A) # leaf
            sage: g = L([lf(1),lf(2)],label="a")
            sage: A = FreeShuffleOperad(QQ,(g,)); A.gens()
            (a[1[], 2[]],)
        """
        return self._gens

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.operads.free_shuffle_operad import FreeShuffleOperad
            sage: L = LabelledOrderedTrees()
            sage: lf = lambda A: L([],label=A) # leaf
            sage: g = L([lf(1),lf(2)],label="a")
            sage: FreeShuffleOperad(QQ,(g,))            # indirect doctest
            The Free shuffle operad over Rational Field with generators (a[1[], 2[]],)
        """
        return "The Free shuffle operad over %s with generators " % (self.base_ring()) + str(self.gens())

    def degree_on_basis(self, x):
        return len(x.leaf_labels())

    def shuffle_composition_on_basis_as_tree(self, x, y, i, sigma):
        """
        Returns the shuffle composition x o_{i,sigma} y of two planar
        trees in the free shuffle operad as a planar tree.

        INPUT:

        - x and y are labeled ordered trees

        - i is an integer between 1 and the number of leaves of x

        - sigma is a permutation that fixes 1,..,i and is a shuffle on the rest

        EXAMPLES::

            sage: from sage.operads.free_shuffle_operad import FreeShuffleOperad
            sage: L = LabelledOrderedTrees()
            sage: lf = lambda A: L([],label=A) # leaf
            sage: g = L([lf(1),lf(2)],label="a")
            sage: A = FreeShuffleOperad(QQ,(g,))
            sage: Trees = A.basis().keys()
            sage: s = Permutation([1,2,3])
            sage: A.shuffle_composition_on_basis_as_tree(g,g,1,s)
            a[a[1[], 2[]], 3[]]
            sage: s = Permutation([1,3,2])
            sage: A.shuffle_composition_on_basis_as_tree(g,g,1,s)
            a[a[1[], 3[]], 2[]]
            sage: s = Permutation([1,2,3])
            sage: A.shuffle_composition_on_basis_as_tree(g,g,2,s)
            a[1[], a[2[], 3[]]]
        """
        # FIRST SHIFT INDICES (ok)
        shifted_y = y.map_leaves(lambda z: z + i - 1)

        def shift(l):
            if l > i:
                return l + self.degree_on_basis(y) - 1
            else:
                return l
        shifted_x = x.map_leaves(shift)
        # then graft (ok)
        gr = self.graft(shifted_x, shifted_y, i)
        # then do the shuffle
        return gr.map_leaves(lambda z: sigma(z))

    def graft(self, x, y, i):
        """
        simple grafting of y on leaf i of x, in a recursive way

        assumes that all leaves have distinct labels

        EXAMPLES::
        """
        if x.node_number() == 1:
            return y
        else:
            with x.clone() as t:
                for j in range(len(t)):
                    if i in t[j].leaf_labels():
                        t[j] = self.graft(t[j], y, i)
            return t

    def shuffle_composition_on_basis(self, x, y, i, sigma):
        """
        This computes the shuffle composition x o_{i,sigma} y as a
        planar tree, for planar trees x and y.

        INPUT:

        - x and y are labeled ordered trees

        - i is an integer between 1 and the number of leaves of x

        - sigma is a permutation

        EXAMPLES::

            sage: from sage.operads.free_shuffle_operad import FreeShuffleOperad
            sage: L = LabelledOrderedTrees()
            sage: lf = lambda A: L([],label=A) # leaf
            sage: g = L([lf(1),lf(2)],label="a")
            sage: A = FreeShuffleOperad(QQ,(g,))
            sage: Trees = A.basis().keys()
            sage: s = Permutation([1,2,3])
            sage: A.shuffle_composition_on_basis(g,g,1,s)
            B[a[a[1[], 2[]], 3[]]]
            sage: s = Permutation([1,3,2])
            sage: A.shuffle_composition_on_basis(g,g,1,s)
            B[a[a[1[], 3[]], 2[]]]
            sage: s = Permutation([1,2,3])
            sage: A.shuffle_composition_on_basis(g,g,2,s)
            B[a[1[], a[2[], 3[]]]]
        """
        return self.basis()[self.shuffle_composition_on_basis_as_tree(x, y, i,
                                                                      sigma)]

    def shuffle_composition(self):
        """
        The composition of the free shuffle operad

        By default, this is implemented from
        :meth:`.shuffle_composition_on_basis`, if available.

        EXAMPLES::

            sage: from sage.operads.free_shuffle_operad import FreeShuffleOperad
            sage: L = LabelledOrderedTrees()
            sage: lf = lambda A: L([],label=A) # leaf
            sage: g = L([lf(1),lf(2)],label="a")
            sage: A = FreeShuffleOperad(QQ,(g,))
            sage: z = A(g)
            sage: s = Permutation([1,2,3])
            sage: A.shuffle_composition()(z,z,1,s)
            B[a[a[1[], 2[]], 3[]]]
            sage: s = Permutation([1,3,2])
            sage: A.shuffle_composition()(z,z,1,s)
            B[a[a[1[], 3[]], 2[]]]
            sage: s = Permutation([1,2,3])
            sage: A.shuffle_composition()(z,z,2,s)
            B[a[1[], a[2[], 3[]]]]
        """
        if self.shuffle_composition_on_basis is not NotImplemented:
            return self._module_morphism(self._module_morphism(self.shuffle_composition_on_basis, position=0, codomain=self), position=1)
        else:
            return NotImplemented
