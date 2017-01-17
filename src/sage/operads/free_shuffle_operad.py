"""
Shuffle operads

Basic implementation of shuffle operads.

REFERENCES:

.. [DoKo] V. Dotsenko and A. Khoroshkin, Anton, *Gröbner bases for operads*.
   Duke Math. J. 153 (2010), no. 2, 363-396.
"""
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.ordered_tree import LabelledOrderedTrees


class FreeShuffleOperad(CombinatorialFreeModule):
    r"""
    The free shuffle operad over any given set of generators.

    Generators are given as a list of labelled ordered trees.
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

    def cmp_fun(self, order=None):
        """
        ordering according to a term order

        IL FAUT UTILISER plutot une "key" !

        possible values so far for order:

        - 'pathlex'
        """
        if order is None:
            order = 'pathlex'
        if order == 'pathlex':
            return lambda x, y: self.pathlex_cmp(x, y)

    def pathlex_cmp(self, x, y):
        """
        term order according to Dotsenko-Koroshkin

        First one compares the number of leaves,

        then one compares the words one by one in deglex order,

        then one compares the permutations in revlex order.

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
        Return the generators of the given free shuffle operad.

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
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.operads.free_shuffle_operad import FreeShuffleOperad
            sage: L = LabelledOrderedTrees()
            sage: lf = lambda A: L([],label=A) # leaf
            sage: g = L([lf(1),lf(2)],label="a")
            sage: FreeShuffleOperad(QQ,(g,))            # indirect doctest
            The Free shuffle operad over Rational Field with generators (a[1[], 2[]],)
        """
        msg = "The Free shuffle operad over {} with generators {}"
        return msg.format(self.base_ring(), self.gens())

    def degree_on_basis(self, x):
        """
        Return the degree of a basis element ``x`` (a planar tree).

        This is the number of leaves.

        EXAMPLES::

            sage: from sage.operads.free_shuffle_operad import FreeShuffleOperad
            sage: L = LabelledOrderedTrees()
            sage: lf = lambda A: L([],label=A) # leaf
            sage: g = L([lf(1),lf(2)],label="a")
            sage: A = FreeShuffleOperad(QQ,(g,))
            sage: A.degree_on_basis(g)
            2
        """
        return len(x.leaf_labels())

    def shuffle_composition_on_basis_as_tree(self, x, y, i, sigma):
        r"""
        Return the shuffle composition of two planar
        trees in the free shuffle operad as a planar tree.

        The shuffle composition x o_{i,sigma} y depends on a position i and
        a permutation sigma.

        INPUT:

        - ``x`` and ``y`` are labeled ordered trees

        - ``i`` is an integer between 1 and the number of leaves of ``x``

        - ``sigma`` is a permutation that fixes 1,..,i and is a shuffle on
          the rest

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
        Return the simple grafting of ``y`` on leaf ``i`` of ``x``.

        This is done in a recursive way, and assumes that all leaves
        have distinct labels.

        INPUT:

        - ``x`` and ``y`` are labeled ordered trees

        - ``i`` is an integer between 1 and the number of leaves of ``x``

        OUTPUT:

        a labelled ordered tree

        EXAMPLES::

            sage: from sage.operads.free_shuffle_operad import FreeShuffleOperad
            sage: L = LabelledOrderedTrees()
            sage: lf = lambda A: L([],label=A) # leaf
            sage: g = L([lf(1),lf(2)],label="a")
            sage: h = L([lf(2),lf(3)],label="a")
            sage: A = FreeShuffleOperad(QQ,(g,))
            sage: A.graft(g, h, 2)
            ?
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
        Return the shuffle composition of two planar
        trees in the free shuffle operad.

        The shuffle composition x o_{i,sigma} y depends on a position i and
        a permutation sigma.

        INPUT:

        - ``x`` and ``y`` are labeled ordered trees

        - ``i`` is an integer between 1 and the number of leaves of ``x``

        - ``sigma`` is a permutation

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
        Return the composition operation of the free shuffle operad.

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
        return NotImplemented


def divisors_of_shape(alpha, beta):
    """
    Find all occurences of the tree beta as a divisor of the tree alpha.

    pas au point, il faudrait renvoyer des arbres avec trous ?
    """
    for t in alpha.subtrees():
        b = extract_divisor_at_root(t, beta)
        if b is False:
            pass
        else:
            yield b  # il faudrait renvoyer l'arbre total


def extract_divisor_at_root(alpha, beta):
    """
    Check if beta is a divisor of alpha at the root

    If yes, return the tree with beta collapsed.

    If no, return ``False``

    pas au point a cause des etiquettes ! penible..

    EXAMPLES::

        sage: L = LabelledOrderedTrees()
        sage: lf = lambda A: L([],label=A)  # leaf
        sage: ga = L([lf(1),lf(2)],label="a")
        sage: gb = L([lf(1),lf(2)],label="b")
        sage: h = L([ga,lf(3)],label="a")
        sage: i = L([L([lf(1),lf(2)],'b'),L([lf(3),lf(4)],'c')],label="a"
        sage: extract_divisor_at_root(ga,h)
        False

        sage: x0 = L([lf(1),lf(3)],label="a")
        sage: x1 = L([x0,lf(2)],label="a")
        sage: x2 = L([lf(1),lf(2)],label="a")
        sage: x3 = L([x2,lf(3)],label="a")
        sage: extract_divisor_at_root(x1, x3)
        False

    Large example::

        sage: x13 = L([lf(1),lf(3)],label="a")
        sage: x45 = L([lf(4),lf(5)],label="a")
        sage: x1345 = L([x13,x45], label="a")
        sage: large = L([x1345,lf(2)], label="a")
        sage: small = L([x13,lf(2)], label="a")
        sage: extract_divisor_at_root(large, small)
        False
    """
    if len(alpha) != len(beta) or alpha.label() != beta.label():
        return False

    list_match = []
    alpha_labels = []
    beta_labels = []
    for i, beta_i in enumerate(beta):
        alpha_i = alpha[i]
        if len(beta_i):
            match = extract_divisor_at_root(alpha_i, beta_i)
            if match is False:
                return False
            list_match.append(match)
            alpha_labels.append(leftmost_leaf_label(match))
        else:
            list_match.append(LabelledOrderedTree([alpha_i], label="H"))
            alpha_labels.append(leftmost_leaf_label(alpha_i))
        beta_labels.append(leftmost_leaf_label(beta_i))

    wa = Word(alpha_labels).standard_permutation()
    wb = Word(beta_labels).standard_permutation()
    print wa, wb
    if wa != wb:
        return False

    return LabelledOrderedTree([subtree for t in list_match for subtree in t],
                               label="H")



def leftmost_leaf_label(tree):
    """
    c'est aussi l'etiquette minimale
    """
    if not(len(tree)):
        return tree.label()
    return leftmost_leaf_label(tree[0])
