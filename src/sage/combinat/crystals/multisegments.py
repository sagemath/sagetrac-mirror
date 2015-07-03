r"""
Crystal of Bernstien-Zelevinsky Multisegments
"""

#*****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
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
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.all import ZZ

class InfinityCrystalOfMultisegments(Parent, UniqueRepresentation):
    r"""
    The type `A_n^{(1)}` crystal `B(\infty)` realized using
    Bernstein-Zelevinsky (BZ) multisegments.

    Using (a modified version of the) notation from [JL2009]_, for `\ell \in
    \ZZ_{>0}` and `i \in \ZZ /(n+1)\ZZ`, a segment of length `\ell` and head `i`
    is the sequence of consecutive residues `[i,i+1,\dots,i+\ell-1]`.  The
    notation  for a segment of length `\ell` and head `i` is simplified to
    `[i;\ell)`.  Similarly, a segment of length `\ell` and tail `i` is the
    sequence of consecutive residues `[i-\ell+1,\dots,i-1,i]`.  The latter is
    denoted simply by `(\ell;i]`.  Finally, a multisegment is a formal linear
    combination of segments, usually written in the form

    .. MATH::

        \psi =
        \sum_{\substack{i \in \ZZ/(n+1)\ZZ \\ \ell \in \ZZ_{>0}}}
        m_{(\ell;i]}(\ell;i].

    Such a multisegment is called aperiodic if, for every `\ell>0`, there exists
    some `i \in \ZZ/(n+1)\ZZ` such that `(\ell;i]` does not appear in `\psi`.
    Denote the set of all periodic multisegments, together with the empty
    multisegment `\varnothing`, by `\Psi`.  We define a crystal
    structure on multisegments as follows.  Set `S_{\ell,i} = \sum_{k\ge \ell}
    (m_{(k;i-1]} - m_{(k;i]})` and let `\ell_f` be the minimal `\ell` that
    attains the value `\min_{\ell>0} S_{\ell,i}`. Then `f_i\psi =
    \psi_{\ell_f,i}`, where

    .. MATH::

        \psi_{\ell_f,i} =
        \begin{cases}
         \psi + (1;i] & \text{ if } \ell_f = 1,\\
         \psi + (\ell_f;i] - (\ell_f-1;i-1] & \text{ if } \ell_f > 1.
        \end{cases}

    Similarly, let `\ell_e` be the maximal `\ell` that attains the value
    `\min_{\ell>0} S_{\ell,i}`.  Then `e_i\psi = \psi_{\ell_e,i}`, where

    .. MATH::

        \psi_{\ell_e,i} =
        \begin{cases}
         0 & \text{ if } \min_{\ell>0} S_{\ell,i} = 0, \\
         \psi + (1;i] & \text{ if } \ell_e = 1,\\
         \psi - (\ell_e;i] + (\ell_e-1;i-1] & \text{ if } \ell_e > 1.
        \end{cases}

    INPUT:

    - ``n`` -- for type `A_n^{(1)}`

    EXAMPLES::

        sage: B = crystals.infinity.Multisegments(2)
        sage: x = B([(8,1),(6,0),(5,1),(5,0),(4,0),(4,1),(4,1),(3,0),(3,0),(3,1),(3,1),(1,0),(1,2),(1,2)]); x
        {(8; 1], (6; 0], (5; 0], (5; 1], (4; 0], (4; 1], (4; 1],
         (3; 0], (3; 0], (3; 1], (3; 1], (1; 0], (1; 2], (1; 2]}
        sage: x.f(1)
        {(8; 1], (6; 0], (5; 0], (5; 1], (4; 0], (4; 1], (4; 1],
         (3; 0], (3; 0], (3; 1], (3; 1], (2; 1], (1; 2], (1; 2]}
        sage: x.f(1).f(1)
        {(8; 1], (6; 0], (5; 0], (5; 1], (4; 0], (4; 1], (4; 1], (3; 0],
         (3; 0], (3; 1], (3; 1], (2; 1], (1; 1], (1; 2], (1; 2]}
        sage: x.e(1)
        {(7; 0], (6; 0], (5; 0], (5; 1], (4; 0], (4; 1], (4; 1], (3; 0], (3; 0], (3; 1], (3; 1], (1; 0], (1; 2], (1; 2]}
        sage: x.e(1).e(1)
        sage: x.f(0)
        {(8; 1], (6; 0], (5; 0], (5; 1], (4; 0], (4; 1], (4; 1],
         (3; 0], (3; 0], (3; 1], (3; 1], (2; 0], (1; 0], (1; 2]}

    We check an `\widehat{\mathfrak{sl}}_2` example against the generalized
    Young walls::

        sage: B = crystals.infinity.Multisegments(1)
        sage: S = B.subcrystal(max_depth=4)
        sage: G = B.digraph(subset=S)
        sage: C = crystals.infinity.GeneralizedYoungWalls(1)
        sage: SC = C.subcrystal(max_depth=4)
        sage: GC = C.digraph(subset=SC)
        sage: G.is_isomorphic(GC, edge_labels=True)
        True

    REFERENCES:

    .. [JL2009] Nicolas Jacon and Cedric Lecouvey.
       *Kashiwara and Zelevinsky involutions in affine type A*.
       Pac. J. Math. 243(2):287-311 (2009).

    .. [LTV1999] Bernard Leclerc, Jean-Yves Thibon, and Eric Vasserot.
       *Zelevinsky's involution at roots of unity*.
       J. Reine Angew. Math. 513:33-51 (1999).
    """
    def __init__(self, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Multisegments(2)
            sage: TestSuite(B).run()
        """
        self._cartan_type = CartanType(['A', n, 1])
        Parent.__init__(self, category=(HighestWeightCrystals(), InfiniteEnumeratedSets()))
        self.module_generators = (self.module_generator(),)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: crystals.infinity.Multisegments(2)
            The infinity crystal of multisegments of type ['A', 2, 1]
        """
        return "The infinity crystal of multisegments of type {}".format(self._cartan_type)

    def module_generator(self):
        """
        Return the module generator of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Multisegments(2)
            sage: B.module_generator()
            {}
        """
        return self.element_class(self, ())

    def weight_lattice_realization(self):
        """
        Return a realization of the weight lattice of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.Multisegments(2)
            sage: B.weight_lattice_realization()
            Extended weight lattice of the Root system of type ['A', 2, 1]
        """
        return self._cartan_type.root_system().weight_lattice(extended=True)

    class Element(ElementWrapper):
        """
        An element in a BZ multisegments crystal.
        """
        def __init__(self, parent, value):
            """
            Initialize ``self``.

            EXAMPLES::

                sage: B = crystals.infinity.Multisegments(2)
                sage: mg = B.module_generator()
                sage: TestSuite(mg).run()
            """
            ZM = IntegerModRing(parent.cartan_type().rank())
            value = [(k, ZM(i)) for k,i in value]
            def sort_cmp(x,y):
                c = cmp(y[0],x[0])
                if not c:
                    c = cmp(ZZ(x[1]),ZZ(y[1]))
                return c
            ElementWrapper.__init__(self, parent, tuple(sorted(value, cmp=sort_cmp)))

        def _repr_(self):
            r"""
            Return a string representation of ``self``.

            EXAMPLES::

                sage: B = crystals.infinity.Multisegments(2)
                sage: B([(4,2), (3,0), (3,1), (1,1), (1,0)])
                {(4; 2], (3; 0], (3; 1], (1; 0], (1; 1]}
            """
            seg = lambda x: "({}; {}]".format(x[0], x[1])
            return "{" + ", ".join(seg(x) for x in self.value) + "}"

        def _latex_(self):
            r"""
            Return a LaTeX representation of ``self``.

            EXAMPLES::

                sage: B = crystals.infinity.Multisegments(2)
                sage: latex(B([(4,2), (3,0), (3,1), (1,1), (1,0)]))
                \bigl\{(4; 2], (3; 0], (3; 1], (1; 0], (1; 1]\bigr\}
            """
            seg = lambda x: "({}; {}]".format(x[0], x[1])
            return "\\bigl\\{" + ", ".join(seg(x) for x in self.value) + "\\bigr\\}"

        def _sig(self, i):
            r"""
            Return an `i`-signature of ``self``.

            INPUT:

            - ``i`` -- an element of the indexing set

            OUTPUT:

            A pair ``(m, p, ep)`` where ``m`` and ``p`` correspond to the
            block length of the unmatched `-` and `+` respectively or ``None``
            if there is no such block and ``ep`` is the number of unmatched
            `-`.

            EXAMPLES::

                sage: B = crystals.infinity.Multisegments(2)
                sage: b = B([(4,2), (3,0), (3,1), (1,1), (1,0)])
                sage: b._sig(0)
                (1, None, 1)
                sage: b._sig(1)
                (None, None, 0)
            """
            if not self.value:
                return (None, None)
            pos = []
            block = self.value[0][0]
            cur = 0
            for k,j in self.value:
                if k != block:
                    if cur != 0:
                        pos.append((block, cur))
                    cur = 0
                    block = k
                if j + 1 == i: # + or (
                    cur += 1
                elif j == i: # - or )
                    cur -= 1
            if cur != 0:
                pos.append((block, cur))
            # Now cancel all +- pairs
            cur = 0
            m = None
            p = None
            ep = 0
            for k,c in pos:
                old = cur
                cur += c
                if cur < 0:
                    m = k
                    p = None
                    ep -= cur
                    cur = 0
                elif not cur:
                    p = None
                elif cur > 0 and old <= 0:
                    p = k
            return (m, p, ep)

        def e(self, i):
            r"""
            Return the action of `e_i` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: B = crystals.infinity.Multisegments(2)
                sage: b = B([(4,2), (3,0), (3,1), (1,1), (1,0)])
                sage: b.e(0)
                {(4; 2], (3; 0], (3; 1], (1; 1]}
                sage: b.e(1)
                sage: b.e(2)
                {(3; 0], (3; 1], (3; 1], (1; 0], (1; 1]}
            """
            i = IntegerModRing(self.parent()._cartan_type.rank())(i)
            m = self._sig(i)[0]
            if m is None:
                return None

            M = self.value
            a = M.index((m, i))
            k = M[a][0]
            if k == 1:
                return self.__class__(self.parent(), M[:a] + M[a+1:])
            return self.__class__(self.parent(), M[:a] + ((k-1,i-1),) + M[a+1:])

        def f(self, i):
            r"""
            Return the action of `f_i` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: B = crystals.infinity.Multisegments(2)
                sage: b = B([(4,2), (3,0), (3,1), (1,1), (1,0)])
                sage: b.f(0)
                {(4; 2], (3; 0], (3; 1], (1; 0], (1; 0], (1; 1]}
                sage: b.f(1)
                {(4; 2], (3; 0], (3; 1], (1; 0], (1; 1], (1; 1]}
                sage: b.f(2)
                {(4; 2], (4; 2], (3; 0], (1; 0], (1; 1]}
            """
            i = IntegerModRing(self.parent()._cartan_type.rank())(i)
            p = self._sig(i)[1]
            M = self.value
            if p is None:
                return self.__class__(self.parent(), ((1, i),) + M)

            a = M.index((p, i-1))
            return self.__class__(self.parent(), M[:a] + ((M[a][0]+1,i),) + M[a+1:])

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: B = crystals.infinity.Multisegments(2)
                sage: b = B([(4,2), (3,0), (3,1), (1,1), (1,0)])
                sage: b.epsilon(0)
                1
                sage: b.epsilon(1)
                0
                sage: b.epsilon(2)
                1
            """
            i = IntegerModRing(self.parent()._cartan_type.rank())(i)
            return self._sig(i)[2]

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            Let `\psi \in \Psi`. Define `\varphi_i(\psi) :=
            \varepsilon_i(\psi) + \langle h_i, \mathrm{wt}(\psi) \rangle`,
            where `h_i` is the `i`-th simple coroot and `\mathrm{wt}(\psi)` is the
            :meth:`weight` of `\psi`.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::

                sage: B = crystals.infinity.Multisegments(2)
                sage: b = B([(4,2), (3,0), (3,1), (1,1), (1,0)])
                sage: b.phi(0)
                1
                sage: b.phi(1)
                0
                sage: mg = B.module_generator()
                sage: mg.f(1).phi(0)
                -1
            """
            h = self.parent().weight_lattice_realization().simple_coroots()
            return self.epsilon(i) + self.weight().scalar(h[i])

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: B = crystals.infinity.Multisegments(2)
                sage: b = B([(4,2), (3,0), (3,1), (1,1), (1,0)])
                sage: b.weight()
                4*delta
            """
            WLR = self.parent().weight_lattice_realization()
            alpha = WLR.simple_roots()
            n = self.parent()._cartan_type.rank()
            return WLR.sum(alpha[j % n] for k,i in self.value for j in range(ZZ(i),ZZ(i)+k))

