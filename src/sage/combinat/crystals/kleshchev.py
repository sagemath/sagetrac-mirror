"""
Kleshchev Partition (Tuple) Crystals
"""
#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.all import ZZ
from sage.combinat.partition_tuple import PartitionTuple
from sage.combinat.root_system.cartan_type import CartanType

class KleshchevCrystal(Parent, UniqueRepresentation):
    r"""
    The Kleshchev multipartition (or partition tuple) crystal.

    We consider type `A_n^{(1)}` crystals, and let `r = (r_i \mid r_i \in
    \ZZ / n \ZZ)` be a finite sequence of length `k` (often called the level)
    and `\lambda = \sum_i \Lambda_{r_i}`. We will model the highest weight
    `U_q(\mathfrak{g})`-crystal `B(\lambda)` by a particular subset of
    partition tuples of level `k`. We define a regular crystal structure
    as follows:

    We create a word in the letters `A_i` and `R_i` by taking the partitions
    (components) from right to left and reading from top to bottom in each
    partition (in English notation) where an addable `i`-cell is an `A_i` and
    a removable `i`-cell is an `R_i`. We then cancel pairs `R_i A_i`, so the
    reduced word is of the form `A_i \cdots A_i R_i \cdots R_i`. Now `e_i`
    (resp. `f_i`) acts by removing (resp. adding) the cell corresponding to
    the left-most `R_i` (resp. right-most `A_i`) in the reduced word.

    The Kleshchev crystal `B(\lambda)` is the crystal generated by the empty
    partition tuple. We can compute the weight of an element `\mu` by taking
    `\lambda - \sum_{i=0}^n c_i \alpha_i` where `c_i` is the number of cells
    of `n`-residue `i` in `\mu`. Partition tuples in the crystal are known
    as `r`-Kleshchev partition tuples, and if `r = (r_1)`, then the partitions
    are known as *Kleshchev* partitions.

    .. NOTE::

        We can describe `r`-Kleshchev partition tuples in `B(\lambda)` as
        partition tuples `\mu` such that `\mu^{(t)}_{r_t - r_{t+1} + x} <
        \mu^{(t+1)}_x` for all `x \geq 1` and `1 \leq t \leq k - 1`.

    INPUT:

    - ``n`` -- for type `A_n^{(1)}`
    - ``r`` -- the sequence `r`

    EXAMPLES:

    We first do an example of a level 1 crystal::

        sage: C = crystals.KleshchevPartitions(2, [0])
        sage: C
        The crystal of multipartitions of type ['A', 2, 1] and residues (0,)
        sage: mg = C.module_generators[0]
        sage: mg
        []
        sage: mg.f(0)
        [1]
        sage: mg.f(1)
        sage: mg.f(2)
        sage: mg.f_string([0,2,1,0])
        [1, 1, 1, 1]
        sage: mg.f_string([0,1,2,0])
        [2, 2]
        sage: S = C.subcrystal(max_depth=5)
        sage: G = C.digraph(subset=S)
        sage: B = crystals.LSPaths(['A',2,1], [1,0,0])
        sage: SB = B.subcrystal(max_depth=5)
        sage: GB = B.digraph(subset=SB)
        sage: G.is_isomorphic(GB, edge_labels=True)
        True

    Now a higher level crystal::

        sage: C = crystals.KleshchevPartitions(2, [0,2])
        sage: mg = C.module_generators[0]
        sage: mg
        ([], [])
        sage: mg.f(0)
        ([1], [])
        sage: mg.f(2)
        ([], [1])
        sage: mg.f_string([0,1,2,0])
        ([2, 2], [])
        sage: mg.f_string([0,2,1,0])
        ([1, 1, 1, 1], [])
        sage: mg.f_string([2,0,1,0])
        ([2], [2])
        sage: S = C.subcrystal(max_depth=3)
        sage: G = C.digraph(subset=S)
        sage: B = crystals.LSPaths(['A',2,1], [1,0,1])
        sage: SB = B.subcrystal(max_depth=3)
        sage: GB = B.digraph(subset=SB)
        sage: G.is_isomorphic(GB, edge_labels=True)
        True

    The ordering of the residues gives a different representation of the
    higher level crystals (but it is still isomorphic)::

        sage: C2 = crystals.KleshchevPartitions(2, [2,0])
        sage: mg2 = C2.highest_weight_vector()
        sage: mg2.f_string([0,1,2,0])
        ([2], [2])
        sage: mg2.f_string([0,2,1,0])
        ([1, 1, 1], [1])
        sage: mg2.f_string([2,0,1,0])
        ([2, 1], [1])
        sage: S2 = C2.subcrystal(max_depth=5)
        sage: G2 = C2.digraph(subset=S)
        sage: G.is_isomorphic(G2, edge_labels=True)
        True

    REFERENCES:

    .. [Ariki2001] Susumu Ariki. On the classification of simple modules for
       cyclotomic Hecke algebras of type `G(m,1,n)` and Kleshchev
       multipartitions. Osaka J. Math. **38** (2001). :arxiv:`9908004v2`.

    .. [Vazirani2002] Monica Vazirani. *Parameterizing Hecek algebra modules:
       Bernstein-Zelevinsky multisegments, Kleshchev multipartitions, and
       crystal graphs*. Transform. Groups **7** (2002). pp. 267-303.
       :arxiv:`0107052v1`, :doi:`10.1007/s00031-002-0014-1`.

    .. [TingleyLN] Peter Tingley. Explicit `\widehat{\mathfrak{sl}}_n` crystal
       maps between cylindric plane partitions, multi-partitions, and
       multi-segments. Lecture notes.
       http://webpages.math.luc.edu/~ptingley/lecturenotes/explicit_bijections.pdf

    .. [Tingley2007] Peter Tingley. Three combinatorial models for
       `\widehat{\mathfrak{sl}}_n` crystals, with applications to cylindric
       plane partitions. International Mathematics Research Notices. (2007).
       :arxiv:`0702062v3`.
    """
    @staticmethod
    def __classcall_private__(cls, n, r):
        """
        Normalize input to ensure a uniqure representation.

        EXAMPLES::

            sage: C1 = crystals.KleshchevPartitions(2, [0,2])
            sage: C2 = crystals.KleshchevPartitions(2, (0,2))
            sage: C1 is C2
            True
        """
        if r in ZZ:
            r = [r]
        M = IntegerModRing(n+1)
        r = tuple(map(M, r))
        return super(KleshchevCrystal, cls).__classcall__(cls, n, r)

    def __init__(self, n, r):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(2, [0,2])
            sage: TestSuite(C).run() # long time
        """
        self._cartan_type = CartanType(['A', n, 1])
        self._r = r
        Parent.__init__(self, category=(HighestWeightCrystals(), RegularCrystals()))
        self.module_generators = (self.element_class(self, PartitionTuple([[]]*len(r))),)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: crystals.KleshchevPartitions(2, [0,2])
            The crystal of multipartitions of type ['A', 2, 1] and residues (0, 2)
        """
        return "The crystal of multipartitions of type {} and residues {}".format(self._cartan_type, self._r)

    class Element(ElementWrapper):
        """
        An element in the multipartition crystal.
        """
        def e(self, i):
            r"""
            Return the action of `e_i` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: C = crystals.KleshchevPartitions(2, [0,2])
                sage: x = C(PartitionTuple([[5,4,1],[3,2,1,1]]))
                sage: x.e(2)
                ([5, 4, 1], [3, 1, 1, 1])
            """
            rem_cells = []
            rs = self.parent()._r
            rev_enum = lambda l: reversed(list(enumerate(l)))

            for j,p in rev_enum(self.value.components()):
                r = rs[j]
                for a,v in enumerate(p):
                    res = v - a + r
                    if res - 1 == i and (a >= len(p)-1 or p[a] > p[a+1]): # Removable
                        rem_cells.append((j,a,v-1))
                    elif res == i and len(rem_cells) > 0 and (a == 0 or p[a-1] > p[a]) > 0: # Addable
                        rem_cells.pop()
                if r - len(p) == i and len(rem_cells) > 0: # Last addable cell
                    rem_cells.pop()

            if len(rem_cells) == 0:
                return None
            c = rem_cells[0]
            if len(rs) == 1: # Special case when it is a single partition
                c = c[1:]
            return self.__class__(self.parent(), self.value.remove_cell(*c))

        def f(self, i):
            r"""
            Return the action of `f_i` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: C = crystals.KleshchevPartitions(2, [0,2])
                sage: x = C(PartitionTuple([[5,4,1],[3,2,1,1]]))
                sage: x.e(2)
                ([5, 4, 1], [3, 1, 1, 1])
                sage: x.f(2)
                ([5, 4, 1], [4, 2, 1, 1])
            """
            add_cells = []
            rs = self.parent()._r
            rev_enum = lambda l: reversed(list(enumerate(l)))

            for j,p in enumerate(self.value.components()):
                r = rs[j]
                if r - len(p) == i: # Last addable cell
                    add_cells.append((j,len(p),0))
                for a,v in rev_enum(p):
                    res = v - a + r
                    if res == i and (a == 0 or p[a-1] > p[a]): # Addable
                        add_cells.append((j,a,v))
                    elif res - 1 == i and len(add_cells) > 0 and (a >= len(p)-1 or p[a] > p[a+1]): # Removable
                        add_cells.pop()

            if len(add_cells) == 0:
                return None
            c = add_cells[0]
            if len(rs) == 1: # Special case when it is a single partition
                c = c[1:]
            return self.__class__(self.parent(), self.value.add_cell(*c))

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: C = crystals.KleshchevPartitions(2, [0,2])
                sage: x = C(PartitionTuple([[5,4,1],[3,2,1,1]]))
                sage: x.epsilon(2)
                1
            """
            ep = 0
            rs = self.parent()._r
            rev_enum = lambda l: reversed(list(enumerate(l)))
            for j,p in rev_enum(self.value.components()):
                r = rs[j]
                for a,v in enumerate(p):
                    res = v - a + r
                    if res - 1 == i and (a >= len(p)-1 or p[a] > p[a+1]): # Addable
                        ep += 1
                    elif res == i and ep > 0 and (a == 0 or p[a-1] > p[a]): # Removable
                        ep -= 1
                if r - len(p) == i and ep > 0: # Last addable cell
                    ep -= 1
            return ep

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: C = crystals.KleshchevPartitions(2, [0,2])
                sage: x = C(PartitionTuple([[5,4,1],[3,2,1,1]]))
                sage: x.phi(2)
                1
            """
            phi = 0
            rs = self.parent()._r
            rev_enum = lambda l: reversed(list(enumerate(l)))
            for j,p in enumerate(self.value.components()):
                r = rs[j]
                if r - len(p) == i: # Last addable cell
                    phi += 1
                for a,v in rev_enum(p):
                    res = v - a + r
                    if res == i and (a == 0 or p[a-1] > p[a]): # Removable
                        phi += 1
                    elif res - 1 == i and phi > 0 and (a >= len(p)-1 or p[a] > p[a+1]): # Addable
                        phi -= 1
            return phi

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: C = crystals.KleshchevPartitions(2, [0,2])
                sage: x = C(PartitionTuple([[5,4,1],[3,2,1,1]]))
                sage: x.weight()
                3*Lambda[0] - Lambda[1]
                sage: x.Phi() - x.Epsilon()
                3*Lambda[0] - Lambda[1]
            """
            WLR = self.parent().weight_lattice_realization()
            alpha = WLR.simple_roots()
            La = WLR.fundamental_weights()
            r = self.parent()._r
            pt = self.value
            wt = WLR.sum(La[ZZ(x)] for x in r)
            return wt - sum(alpha[pt.content(*c, multicharge=r)] for c in pt.cells())

