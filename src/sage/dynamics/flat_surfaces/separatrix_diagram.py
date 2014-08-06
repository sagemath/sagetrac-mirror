r"""
Saddle configuration, separatrix diagram and cylinder diagram.

A separatrix diagram is a couple of permutation ``(bot,top)`` such that in
there is the same number of cycles in the cycle decompositions of both
``bot`` and ``top``. A cylinder diagram is a separatrix diagram together
with a bijection between the cycles of ``bot`` and the one of ``top``.

A cylinder diagram encodes the combinatorics of cylinder decomposition of a
completely periodic direction in a translation surface. If we adjoin coordinates
to this combinatorial datum, we have a complete description of the underlying
surface. In the case of arithmetic curves, the coordinates can be taken to be
rational numbers.

This representation of a surface is used in various constructions:

 - square tiled surfaces
 - Thurston-Veech construction of pseudo-Anosov diffeomorphism
 - description of the cusp of Teichmueller curves

TODO::

* We need a more general structure to encode configurations of structure of
  saddle connections (which need not be completely periodic directions (see
  [EMZ], [MZ]_)

* Gray code for conjugacy classes of permutation in order to optimize the
  generation of separatrix and cylinder diagrams.


REFERENCES:

.. [EMZ] A. Eskin, H. Masur, A. Zorich "Principal boundary ... and Siegel-Veech
         constant"

.. [MZ]  H. Masur, A. Zorich "Multiple saddle connections on flat surfaces and
         the principal boundary of the moduli spaces of quadratic
         differentials"

.. [N]   Y. Naveh "Tight upper bounds on the number of invariant components on
         translation surfaces", Isr. J. Math. 165, 211-231 (2008)

EXAMPLES:

Separatrix diagrams::

    sage: s = SeparatrixDiagram('(0,1,2)(3,4)(5,6,7)','(0,4,1,2)(3,7)(5,6)')
    sage: s
    Separatrix diagram
     bot (0,1,2)(3,4)(5,6,7)
     top (0,4,1,2)(3,7)(5,6)
    sage: s.bot_cycle_tuples()
    [(0, 1, 2), (3, 4), (5, 6, 7)]
    sage: s.top_cycle_tuples()
    [(0, 4, 1, 2), (3, 7), (5, 6)]

Cylinder diagrams::

    sage: c = CylinderDiagram([((0,),(4,)),((1,2),(0,1,3)),((3,4),(2,))])
    sage: print c
    (0)-(4) (1,2)-(0,1,3) (3,4)-(2)
    sage: print c.separatrix_diagram()
    Separatrix diagram
     bot (0)(1,2)(3,4)
     top (0,1,3)(2)(4)

They can also be built from separatrix diagram::

    sage: s = SeparatrixDiagram('(0,1,2)(3,4)(5,6,7)','(0,4,1,2)(3,7)(5,6)')
    sage: s
    Separatrix diagram
     bot (0,1,2)(3,4)(5,6,7)
     top (0,4,1,2)(3,7)(5,6)
    sage: s.to_cylinder_diagram([(0,1),(1,0),(2,2)])
    (0,1,2)-(3,7) (3,4)-(0,4,1,2) (5,6,7)-(5,6)
"""
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from gray_codes import GrayCodeSwitch

import itertools
import sage.rings.arith as arith
from sage.rings.integer import Integer

from sage.misc.permutation import (perm_check, equalize_perms, init_perm,
        perm_cycle_tuples, perm_cycle_string, perm_compose, perm_compose_i,
        perm_orbit, perm_invert)

#
# Abelian and quadratic Separatrix Diagram
#

# main class

class SeparatrixDiagram(SageObject):
    r"""
    Separatrix diagram of oriented foliation.

    A separatrix diagram is a 2-tuple of permutations ``(bot,top)`` such that
    ``bot`` and ``top`` share the same number of cycles.

    bot (resp. top) has to be thought a bottom (resp. top) of a potential face
    as in the following

            -- bot -->
        -------------------
           <-- top --

    The order for bot and top is choosen in such a way that it cooresponds to
    the orientation of a face.

    EXAMPLES::

        sage: s = SeparatrixDiagram('(0,2)(1,3,4)','(0,4)(2,1,3)')
        sage: print s
        Separatrix diagram
         bot (0,2)(1,3,4)
         top (0,4)(1,3,2)
        sage: print s.stratum()
        H_3(4)
    """
    def __init__(self,data,top=None,check=True):
        r"""
        TESTS::

            sage: s = SeparatrixDiagram('(0,1)(2,3,4)','(0,2,4)(1,3)')
            sage: s == loads(dumps(s))
            True
            sage: s == SeparatrixDiagram(s)
            True
            sage: s == SeparatrixDiagram(str(s))
            True
        """
        if top is None:
            if isinstance(data,SeparatrixDiagram):
                bot = data.bot()
                top = data.top()
            elif isinstance(data,(list,tuple)) and len(data) == 2:
                bot,top = data
            elif isinstance(data,str):
                bot,top = data.split('-')
            else:
                raise ValueError, "the argument data is not valid"
        else:
            bot = data

        self._bot = init_perm(bot)
        self._top = init_perm(top)
        n = equalize_perms((self._bot, self._top))

        bot_seen = [True] * n
        top_seen = [True] * n
        bot_to_cycle = [None]*n
        top_to_cycle = [None]*n
        bot_cycles = []
        top_cycles = []
        for i in xrange(n):
            if bot_seen[i]:
                c=[]
                k = len(bot_cycles)
                while bot_seen[i]:
                    bot_to_cycle[i] = k
                    bot_seen[i] = False
                    c.append(i)
                    i = self._bot[i]
                bot_cycles.append(tuple(c))
                k += 1
            if top_seen[i]:
                c=[]
                k = len(top_cycles)
                while top_seen[i]:
                    top_to_cycle[i] = k
                    top_seen[i] = False
                    c.append(i)
                    i = self._top[i]
                top_cycles.append(tuple(c))
                k += 1

        self._bot_cycles = bot_cycles
        self._top_cycles = top_cycles
        self._bot_to_cycle = bot_to_cycle
        self._top_to_cycle = top_to_cycle

        if check:
            self._check()

    def _check(self):
        r"""
        Check that the data of self is valid, i.e.

          * self._bot, self._top are valid permutations
          * the number of cylces of bot and top are the same

        EXAMPLES::

            sage: c = SeparatrixDiagram('(0,1)(2,3)','(0,2,3,1)') #indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: bot has 2 cylinders whereas top has 1
        """
        perm_check(self._bot)
        perm_check(self._top)

        p_bot = self.bot_cycle_tuples()
        p_top = self.top_cycle_tuples()
        if len(p_top) != len(p_bot):
            raise ValueError, "bot has %d cylinders whereas top has %d" %(len(p_bot),len(p_top))

    def to_directed_graph(self):
        r"""
        Return a graph that encodes this separatrix diagram.

        The vertices correspond to separatrix and the edges are of two types

        - 'b' neighboor corresponds to the right neighbors on the bottom
          permutation
        - 't' edges correspond to the neighbor of the top permutation

        EXAMPLES::

            sage: S = SeparatrixDiagram('(0,1)(2,3,4)','(0,3,2)(1,4)')
            sage: G = S.to_directed_graph(); G
            Looped multi-digraph on 5 vertices
            sage: G.vertices()
            [0, 1, 2, 3, 4]
            sage: G.edges()
            [(0, 1, 'b'), (0, 3, 't'), (1, 0, 'b'), (1, 4, 't'), (2, 0, 't'), (2, 3, 'b'), (3, 2, 't'), (3, 4, 'b'), (4, 1, 't'), (4, 2, 'b')]
        """
        from sage.graphs.digraph import DiGraph

        G = DiGraph(multiedges=True,loops=True)
        for i in xrange(self.nseps()):
            G.add_edge(i,self._top[i],'t')
            G.add_edge(i,self._bot[i],'b')
        return G

    def _repr_(self):
        r"""
        String representation of self

        EXAMPLES::

            sage: d = SeparatrixDiagram('(0,1)(2)','(0)(1,2)')
            sage: repr(d) #indirect doctest
            '(0,1)(2)-(0)(1,2)'
        """
        return self.bot_cycle_string() + "-" + self.top_cycle_string()

    #TODO
    def _latex_(self):
        r"""
        """
        n = self._n
        if len(self.vertices_out().cycle_type()) == 1:
            v = self.vertices()[0]
            m = 360. / (2*n)
            d = dict([i,(v.index(i),v.index(-i))] for i in range(1,self._n+1))
            s = "\\begin{tikzpicture}\n"
            for i,(vout,vin) in d.iteritems():
                s += "    \draw [-triangle 45] (0,0) -- (%f:0.8cm);\n" %(vout*m)
                s += "    \draw (%f:0.8cm) -- (%f:1cm);\n" %(vout*m,vout*m)
                s += "    \draw (%f:1cm) \n" %(vout*m)
                v1 = Integer(vout+vin)/2
                if vout+vin > 2*n:
                    v2 = v1 - n
                else:
                    v2 = v1 + n
                d1 = min([abs(vout-v1), abs(vout-v1-2*n), abs(vout-v1+2*n)])
                d2 = min([abs(vout-v2), abs(vout-v2-2*n), abs(vout-v2+2*n)])
                if d1 < d2:
                    vint = v1
                    d = d1
                else:
                    vint = v2
                    d = d2

                dint = '%fcm' %(1.5+d/2.)
                ct1 = '%fcm' %(d/2.)

                if cyclic_direction(vout,vint,vin) == 1:
                    ct2 = '%fcm' %(-d/2.)
                    ct3 = '%fcm' %(d/2.)
                else:
                    ct2 = '%fcm' %(d/2.)
                    ct3 = '%fcm' %(-d/2.)

                s += "    ..controls +(%f:%s) and +(%f:%s) ..\n" %(vout*m,ct1,(vint+n/2.)*m,ct2)
                s += "    (%f:%s)\n" %(vint*m,dint)
                s += "    ..controls +(%f:%s) and +(%f:%s) ..\n" %((vint+n/2.)*m,ct3,vin*m,ct1)
                s += "    (%f:1cm);\n" %(vin*m)
                s += "    \draw [-open triangle 45] (%f:1cm) -- (%f:0.6cm);\n" %(vin*m,vin*m)
                s += "    \draw (%f:0.6cm) -- (0,0);\n" %(vin*m)
            s += "\\end{tikzpicture}"
            return s
        else:
            return ""

    #
    # Comparisons and canonic labels
    #

    def __eq__(self, other):
        r"""
        Equality test

        TESTS::

            sage: d1 = SeparatrixDiagram('(0)','(0)')
            sage: d2 = SeparatrixDiagram('(0,1)(2)','(0,1)(2)')
            sage: d3 = SeparatrixDiagram('(0,1)(2)','(0,2)(1)')
            sage: d1 == d1 and d2 == d2 and d3 == d3
            True
            sage: d1 == d2 or d1 == d3 or d2 == d3 or d3 == d2
            False
        """
        if not isinstance(other, SeparatrixDiagram):
            raise NotImplemented
        return self._bot == other._bot and self._top == other._top

    def __ne__(self,other):
        r"""
        Difference test

        TESTS::

            sage: d1 = SeparatrixDiagram('(0)','(0)')
            sage: d2 = SeparatrixDiagram('(0,1)(2)','(0,1)(2)')
            sage: d3 = SeparatrixDiagram('(0,1)(2)','(0,2)(1)')
            sage: d1 != d1 or d2 != d2 or d3 != d3
            False
            sage: d1 != d2 and d1 != d3 and d2 != d3 and d3 != d2
            True
        """
        if not isinstance(other,SeparatrixDiagram):
            raise NotImplemented
        return self._bot != other._bot or self._top != other._top

    def __cmp__(self,other):
        r"""
        Comparison

        TESTS::

            sage: s3_0 = SeparatrixDiagram([0,1,2],[0,1,2])
            sage: S3 = [s3_0]

            sage: s2_0 = SeparatrixDiagram([1,0,2],[1,0,2])
            sage: s2_1 = SeparatrixDiagram([1,0,2],[2,1,0])
            sage: s2_2 = SeparatrixDiagram([1,0,2],[0,2,1])
            sage: s2_3 = SeparatrixDiagram([2,1,0],[1,0,2])
            sage: s2_4 = SeparatrixDiagram([2,1,0],[2,1,0])
            sage: s2_5 = SeparatrixDiagram([2,1,0],[0,2,1])
            sage: s2_6 = SeparatrixDiagram([0,2,1],[1,0,2])
            sage: s2_7 = SeparatrixDiagram([0,2,1],[2,1,0])
            sage: s2_8 = SeparatrixDiagram([0,2,1],[0,2,1])
            sage: S2 = [s2_0,s2_1,s2_2,s2_3,s2_4,s2_5,s2_6,s2_7,s2_8]

            sage: s1_0 = SeparatrixDiagram([1,2,0],[1,2,0])
            sage: s1_1 = SeparatrixDiagram([1,2,0],[2,0,1])
            sage: s1_2 = SeparatrixDiagram([2,0,1],[1,2,0])
            sage: s1_3 = SeparatrixDiagram([2,0,1],[2,0,1])
            sage: S1 = [s1_0,s1_1,s1_2,s1_3]

            sage: all(cmp(s,s) == 0 for s in S1+S2+S3)
            True
            sage: all(cmp(s1,s2) == -1 for s1 in S1 for s2 in S2)
            True
            sage: all(cmp(s2,s1) == 1 for s1 in S1 for s2 in S2)
            True
            sage: all(cmp(s2,s3) == -1 for s2 in S2 for s3 in S3)
            True
            sage: all(cmp(s3,s2) == 1 for s2 in S2 for s3 in S3)
            True
            sage: all(cmp(s1,s3) == -1 for s1 in S1 for s3 in S3)
            True
            sage: all(cmp(s3,s1) == 1 for s1 in S1 for s3 in S3)
            True
        """
        if not isinstance(other, SeparatrixDiagram):
            raise NotImplemented

        test = cmp(self.nseps(), other.nseps())
        if test: return test

        test = cmp(self.ncyls(),other.ncyls())
        if test: return test

        test = cmp(map(len,self.bot_cycle_tuples()),map(len,other.bot_cycle_tuples()))
        if test: return test

        test = cmp(map(len,self.top_cycle_tuples()),map(len,other.top_cycle_tuples()))
        if test: return test

        test = cmp(self._bot,other._bot)
        if test: return test

        test = cmp(self._top,other._top)
        if test: return test

        return 0

    def is_isomorphic(self, other, return_map=False):
        r"""
        Test whether self is isomorphic to other.

        EXAMPLES::

            sage: bot = [1,2,0,3]
            sage: top = [1,0,3,2]
            sage: s = SeparatrixDiagram(bot,top); s
            Separatrix diagram
             bot (0,1,2)(3)
             top (0,1)(2,3)
            sage: m = [3,0,1,2]
            sage: bot2 = [0]*4
            sage: top2 = [0]*4
            sage: for i in xrange(4):
            ...     bot2[m[i]] = m[bot[i]]
            ...     top2[m[i]] = m[top[i]]
            sage: ss = SeparatrixDiagram(bot2,top2)
            sage: s.is_isomorphic(ss)
            True
            sage: m = [1,2,0,3]
            sage: for i in xrange(4):
            ...     bot2[m[i]] = m[bot[i]]
            ...     top2[m[i]] = m[top[i]]
            sage: ss = SeparatrixDiagram(bot2,top2)
            sage: s.is_isomorphic(ss)
            True
        """
        if not isinstance(other, SeparatrixDiagram):
            raise NotImplemented
        return self.relabel(inplace=False) == other.relabel(inplace=False)

    def relabel(self, perm=None, return_map=False, inplace=True):
        r"""
        New labels for self.


        INPUT:

        - ``perm`` - a permutation as a list of elements of [0,...,n] or None
          (default: None) - the permutation used to relabel the separatrix
          diagram. If None, then canonic labels are used.

        - ``return_map`` - boolean (default: False) - whether or not return the
          permutation used to relabel.

        - ``inplace`` - boolean (default: False) - if True modify self if not
          return a relabeled copy.

        EXAMPLES::

        Perform some renumbering of a separatrix diagram using a permutation::

            sage: S=SeparatrixDiagram('(0)(2,3,4)','(0,3,2)(1)'); S
            Separatrix diagram
             bot (0)(1)(2,3,4)
             top (0,3,2)(1)(4)
            sage: S.relabel(perm=[1,0,2,3,4],inplace=False)
            Separatrix diagram
             bot (0)(1)(2,3,4)
             top (0)(1,3,2)(4)
            sage: S.relabel(perm=[1,2,0,3,4],inplace=False)
            Separatrix diagram
             bot (0,3,4)(1)(2)
             top (0,1,3)(2)(4)

        Canonical labels::

            sage: bot = '(0,1,3,6,7,5)(2,4)(8)(9)'
            sage: top = '(0)(1,2)(3,4,5)(6,7,8,9)'
            sage: s = SeparatrixDiagram(bot,top)
            sage: s.relabel(inplace=False)
            Separatrix diagram
             bot (0)(1)(2,3,4,5,6,7)(8,9)
             top (0,1,2,3)(4,7,9)(5)(6,8)

        TESTS::

            sage: top = [1,0,3,4,2]
            sage: bot = [3,2,4,0,1]
            sage: t = [None]*5; b = [None]*5
            sage: for p in permutations([0,1,2,3,4]):
            ...       for i in xrange(5):
            ...           t[p[i]] = p[top[i]]
            ...           b[p[i]] = p[bot[i]]
            ...       s = SeparatrixDiagram(b,t)
            ...       s.relabel()
            ...       print s
            Separatrix diagram
             bot (0,1)(2,3,4)
             top (0,2,4)(1,3)
            Separatrix diagram
             bot (0,1)(2,3,4)
             top (0,2,4)(1,3)
            Separatrix diagram
             bot (0,1)(2,3,4)
             top (0,2,4)(1,3)
            Separatrix diagram
             bot (0,1)(2,3,4)
             top (0,2,4)(1,3)
            Separatrix diagram
             bot (0,1)(2,3,4)
             top (0,2,4)(1,3)
            Separatrix diagram
             bot (0,1)(2,3,4)
             top (0,2,4)(1,3)
            ...
            Separatrix diagram
             bot (0,1)(2,3,4)
             top (0,2,4)(1,3)
        """
        if perm is None:
            perm = self._compute_normal_form()
            if inplace:
                self.__dict__ = self._normal_form.__dict__
                if return_map:
                    raise NotImplementedError
                return None

            if return_map:
                raise NotImplementedError
            return self._normal_form

        n = self.degree()
        perm.extend(xrange(len(perm),n))
        bot = [None] * self.degree()
        top = [None] * self.degree()

        for i in xrange(n):
            bot[perm[i]] = perm[self._bot[i]]
            top[perm[i]] = perm[self._top[i]]

        S = SeparatrixDiagram(bot,top)
        if inplace:
            self._top = S._top
            self._bot = S._bot
            self._top_cycles = S._top_cycles
            self._bot_cycles = S._bot_cycles
            self._top_to_cycle = S._top_to_cycle
            self._bot_to_cycle = S._bot_to_cycle
            if hasattr('self','_normal_form') and self is self._normal_form:
                del self._normal_form

            if return_map:
                return perm
            return None
        if return_map:
            return perm, S
        return S

    def _compute_normal_form(self):
        r"""
        Returns a normal form for self

        1) first compute the orbit of G = <bot,top>

        2) for each of the connected component
          top is the lexicographic minimum (1)(2)(3,4)(5,6,7)
          bot is used to discover new guys

        3) sort (normal_top, normal_bot) and concatenate them

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0,1)(2,3,4)','(0,3,4)(1,2)')
            sage: s.relabel(inplace=False) #indirect doctest
            Separatrix diagram
             bot (0,1)(2,3,4)
             top (0,2)(1,3,4)

            sage: s = SeparatrixDiagram('(0,5,2)(1,3,4)(6,7,8)','(0,3,7,8)(1,5)(2,4,6)')
            sage: s
            Separatrix diagram
             bot (0,5,2)(1,3,4)(6,7,8)
             top (0,3,7,8)(1,5)(2,4,6)
            sage: s.relabel(inplace=False) #indirect doctest
            Separatrix diagram
             bot (0,1,2)(3,4,5)(6,7,8)
             top (0,1,3,6)(2,5,7)(4,8)
        """
        from sage.misc.permutation import perms_canonical_labels, perms_relabel

        if not hasattr(self, '_normal_form'):
            bot = self._bot
            top = self._top
            n = len(bot)

            G = self.to_directed_graph()
            cs = map(tuple,G.connected_components())

            cs_type_nb = {} # (bot,top) -> nb of them

            for c in cs:
                c_inv = dict((c[i],i) for i in xrange(len(c)))
                cbot = [None]*len(c)
                ctop = [None]*len(c)
                for i in c:
                    cbot[c_inv[i]] = c_inv[bot[i]]
                    ctop[c_inv[i]] = c_inv[top[i]]

                (cbot,ctop), _ = perms_canonical_labels([cbot,ctop])

                bt=(tuple(cbot),tuple(ctop))

                if bt in cs_type_nb:
                    cs_type_nb[bt] += 1
                else:
                    cs_type_nb[bt] = 1

            shift = 0
            bot = []
            top = []
            for key in sorted(cs_type_nb.keys()):
                for j in xrange(cs_type_nb[key]):
                    bot.extend(shift+i for i in key[0])
                    top.extend(shift+i for i in key[1])
                    shift += len(key[0])

            self._normal_form = SeparatrixDiagram(bot,top)
            self._normal_form._normal_form = self._normal_form

        return self._normal_form

    def is_in_normal_form(self):
        r"""
        Test normal form

        Return True if self is in normal form and False otherwise.

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0,1,2)(3,4,5)(6,7,8)','(0,3,7,8)(1,5)(2,4,6)')
            sage: s.is_in_normal_form()
            False
            sage: s.relabel(inplace=False).is_in_normal_form()
            True
        """
        self._compute_normal_form()
        return self == self._normal_form

    #
    # Attributes access
    #

    def degree(self):
        r"""
        Return the degree (number of separatrices) of this separatrix diagram.

        EXAMPLES::

            sage: S = SeparatrixDiagram('(0,1)(2,3)','(1,3,2)(0)')
            sage: S.degree()
            4
        """
        return len(self._top)

    nseps = degree

    def ncyls(self):
        r"""
        Return the number of cylinders of this separatrix diagram.

        EXAMPLES::

            sage: S = SeparatrixDiagram('(0,1)(2,3)','(1,3,2)(0)')
            sage: S.ncyls()
            2
        """
        return len(self._top_cycles)

    def profile(self):
        r"""
        Return the angles around each vertex

        EXAMPLES::

            sage: a = AbelianStratum(1,1,0)
            sage: s = a.separatrix_diagrams()[0]
            sage: s.profile()
            [2, 2, 1]
        """
        from sage.combinat.partition import Partition

        p = map(len,perm_cycle_tuples(self.outgoing_edges_perm(),singletons=True))
        return Partition(sorted(p, reverse=True))

    def euler_characteristic(self):
        r"""
        Return the Euler characteristic

        EXAMPLES::

            sage: SeparatrixDiagram('(0)','(0)').euler_characteristic()
            0

            sage: CylinderDiagram([((0,),(0,))]).euler_characteristic()
            0
            sage: CylinderDiagram([((0,1),(0,2)), ((2,),(1,))]).euler_characteristic()
            -2
        """
        p = self.profile()
        return Integer(len(p)-sum(p))

    def genus(self):
        r"""
        Return the genus

        EXAMPLES::

            sage: CylinderDiagram([((0,),(0,))]).genus()
            1
            sage: CylinderDiagram([((0,1),(0,1))]).genus()
            1
            sage: CylinderDiagram([((0,1,2),(0,1,2))]).genus()
            2
            sage: CylinderDiagram([((0,1,2,3),(0,1,2,3))]).genus()
            2
            sage: CylinderDiagram([((0,1,2,3,4),(0,1,2,3,4))]).genus()
            3
        """
        return Integer(1 - self.euler_characteristic()//2)

    def stratum(self):
        r"""
        Return the abelian stratum

        EXAMPLES::

            sage: SeparatrixDiagram('(0)(1)(2)','(0)(1)(2)').stratum()
            H_1(0^3)
            sage: SeparatrixDiagram('(0,1)(2)','(0,2)(1)').stratum()
            H_2(2)
        """
        from abelian_strata import AbelianStratum

        return AbelianStratum([i-1 for i in self.profile()])

    def bot(self):
        r"""
        The bot permutation as a list from 0 to nseps-1

        Warning: the output list should not be modified

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0)(1,2)','(0,1)(2)')
            sage: s.bot()
            [0, 2, 1]
        """
        return list(self._bot)

    def bot_perm(self):
        r"""
        Return the bot as a permutation (element of a group)

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0)(1,2)','(0,1)(2)')
            sage: s.bot_perm()
            (2,3)
        """
        from sage.groups.perm_gps.permgroup_element import PermutationGroupElement

        return PermutationGroupElement([i+1 for i in self._bot])

    def bot_orbit(self, i):
        r"""
        Return the orbit of i under the bot permutation

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0,1)(2,5)(3,4,6)','(0,1,5)(2,3,6)(4)')
            sage: s.bot_orbit(0)
            (0, 1)
            sage: s.bot_orbit(4)
            (3, 4, 6)
        """
        return self._bot_cycles[self._bot_to_cycle[i]]

    def bot_cycle_tuples(self):
        r"""
        Return the cycles of the bottom permutation as a list of tuples.

        EXAMPLES::

            sage: S = SeparatrixDiagram('(0,2)(3,4)','(0)(1,2,3)')
            sage: S.bot_cycle_tuples()
            [(0, 2), (1,), (3, 4)]
        """
        return self._bot_cycles

    def bot_cycle_string(self):
        r"""
        Return the cycles of the top permutation as a string.

        EXAMPLES::

            sage: S = SeparatrixDiagram('(0,2)(3,4)','(0)(1,2,3)')
            sage: S.bot_cycle_string()
            '(0,2)(1)(3,4)'
        """
        return ''.join('(' + ','.join(map(str,c)) +')' for c in self.bot_cycle_tuples())

    def top(self):
        r"""
        Return the top permutation of self as a list.

        Warning: the output should not be modified

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0,1,3)(2,4)','(0,4)(1,2,3)')
            sage: s.top()
            [4, 2, 3, 1, 0]
        """
        return self._top

    def top_perm(self):
        r"""
        Return the top as a permutation

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0)(1,2)','(1)(0,2)')
            sage: s.top_perm()
            (1,3)
        """
        from sage.groups.perm_gps.permgroup_element import PermutationGroupElement

        return PermutationGroupElement([i+1 for i in self._top])

    def top_orbit(self,i):
        r"""
        Return the orbit of ``i`` under the top permutation.

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0,1)(2,5)(3,4,6)','(0,1,5)(2,3,6)(4)')
            sage: s.top_orbit(0)
            (0, 1, 5)
            sage: s.top_orbit(6)
            (2, 3, 6)
        """
        return self._top_cycles[self._top_to_cycle[i]]

    def top_cycle_tuples(self):
        r"""
        Return the cycle of the top permutation as a list of tuples.

        EXAMPLES::

            sage: S = SeparatrixDiagram('(0,2)(3,4)','(0)(1,2,3)')
            sage: S.top_cycle_tuples()
            [(0,), (1, 2, 3), (4,)]
        """
        return self._top_cycles

    def top_cycle_string(self):
        r"""
        Return the cycle of the top permutation as a string.

        EXAMPLES::

            sage: S = SeparatrixDiagram('(0,2)(3,4)','(0)(1,2,3)')
            sage: S.top_cycle_string()
            '(0)(1,2,3)(4)'
        """
        return ''.join('(' + ','.join(map(str,c)) + ')' for c in self.top_cycle_tuples())

    def automorphism_group(self, implementation='graph'):
        r"""
        Return the automorphism group of self.

        That is the centralizer of the permutations top and bottom.

        INPUT:

        - ``implementation`` - either graph or gap

        EXAMPLES::

            sage: S = SeparatrixDiagram('(0,3,1,4,2)','(0,1,2,3,4)')
            sage: S.automorphism_group(implementation='graph')
            Permutation Group with generators [(1,2,3,4,5)]
            sage: S.automorphism_group(implementation='gap')
            Subgroup of (Symmetric group of order 5! as a permutation group) generated by [(1,2,3,4,5), (1,4,2,5,3)]
        """
        if implementation == 'graph':
            return self.to_directed_graph().automorphism_group(edge_labels=True)

        elif implementation == 'gap':
            from sage.groups.perm_gps.permgroup import PermutationGroup
            from sage.groups.perm_gps.permgroup_named import SymmetricGroup

            return SymmetricGroup(self.nseps()).centralizer(PermutationGroup([self.top_perm(),self.bot_perm()]))

        else:
            raise ValueError, "implementation should be either 'graph' or 'gap'"

    def homological_dimension_of_cylinders(self):
        r"""
        Returns the dimension in the first homology group of the span of waist
        curves of horizontal cylinders.

        EXAMPLES:

        Homological dimension in the stratum H(2)::

            sage: c = CylinderDiagram('(0,1,2)-(0,1,2)')
            sage: c.stratum()
            H_2(2)
            sage: c.homological_dimension_of_cylinders()
            1
            sage: c = CylinderDiagram('(0,1)-(1,2) (2)-(0)')
            sage: c.stratum()
            H_2(2)
            sage: c.homological_dimension_of_cylinders()
            2

        Homological dimensions for cylinder diagrams in H(1,1)::

            sage: c = CylinderDiagram('(0,1,2,3)-(0,1,2,3)')
            sage: c.stratum()
            H_2(1^2)
            sage: c.homological_dimension_of_cylinders()
            1
            sage: c = CylinderDiagram('(0,1)-(0,2) (2,3)-(1,3)')
            sage: c.stratum()
            H_2(1^2)
            sage: c.homological_dimension_of_cylinders()
            2
            sage: c = CylinderDiagram('(0,1,2)-(1,2,3) (3)-(0)')
            sage: c.stratum()
            H_2(1^2)
            sage: c.homological_dimension_of_cylinders()
            2
            sage: c = CylinderDiagram('(0,1)-(2,3) (2)-(0) (3)-(1)')
            sage: c.stratum()
            H_2(1^2)
            sage: c.homological_dimension_of_cylinders()
            2
        """
        return Integer(self.ncyls() - SeparatrixDiagram.to_directed_graph(self).connected_components_number() + 1)

    #
    # Vertices of the separatrix diagram
    #

    def outgoing_edges_perm(self):
        r"""
        Permutation associated to turning around vertices in trigonometric
        order.

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0,1)','(2,3)')
            sage: s.outgoing_edges_perm()
            [1, 0, 3, 2]

            sage: s = SeparatrixDiagram('(0,5,2)(1,3,4)(6,7,8)','(0,3,7,8)(1,5)(2,4,6)')
            sage: s.outgoing_edges_perm()
            [7, 0, 8, 2, 5, 4, 3, 1, 6]

        """
        return perm_compose_i(self._bot,self._top)

    def incoming_edges_perm(self):
        r"""
        Permutation associated to turning around vertices in trigonometric
        order.

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0,1)','(2,3)')
            sage: s.incoming_edges_perm()
            [1, 0, 3, 2]

            sage: s = SeparatrixDiagram('(0,5,2)(1,3,4)(6,7,8)','(0,3,7,8)(1,5)(2,4,6)')
            sage: s.incoming_edges_perm()
            [4, 2, 1, 8, 7, 3, 0, 6, 5]
        """
        return perm_compose(self._top,self._bot)

    #
    # to cylinder diagram
    #

    def to_cylinder_diagram(self, pairing):
        r"""
        Return a cylinder diagram with the given pairing

        The pairing should be a list of 2-tuples of integer.

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0,1,3)(2,4)','(0,2)(1,4,3)'); s
            Separatrix diagram
             bot (0,1,3)(2,4)
             top (0,2)(1,4,3)

            sage: s.to_cylinder_diagram([(0,0),(1,1)])
            (0,1,3)-(0,2) (2,4)-(1,4,3)
            sage: s.to_cylinder_diagram([(1,1),(0,0)])
            (0,1,3)-(0,2) (2,4)-(1,4,3)

            sage: s.to_cylinder_diagram([(0,1),(1,0)])
            (0,1,3)-(1,4,3) (2,4)-(0,2)
            sage: s.to_cylinder_diagram([(1,0),(0,1)])
            (0,1,3)-(1,4,3) (2,4)-(0,2)
        """
        from copy import copy

        other = copy(self)
        other.__class__ = CylinderDiagram

        bots = self.bot_cycle_tuples()
        tops = self.top_cycle_tuples()

        other._bot_to_cyl = [None]*self.nseps()
        other._top_to_cyl = [None]*self.nseps()

        for i in xrange(len(pairing)):
            b = bots[pairing[i][0]]
            t = tops[pairing[i][1]]
            cyl = (b[0],t[0])

            for j in b:
                other._bot_to_cyl[j] = cyl
            for j in t:
                other._top_to_cyl[j] = cyl

        return other

    def cylinder_diagram_iterator(self,connected=True,up_to_isomorphism=True):
        r"""
        Construct all cylinder diagrams from given separatrix diagram (i.e. a pair
        of permutations).

        INPUT:

        - ``connected`` - boolean (default: True) - if true, returns only
          connected cylinder diagrams.

        - ``up_to_isomorphism`` - boolean (default: True) - take care of
          isomorphism problem. It is memory efficient and probably faster to set
          this option to False.

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0,1)(2,3)(4,5)','(1,2)(3,4)(5,0)')
            sage: for c in s.cylinder_diagram_iterator(): print c
            (0,5)-(0,4) (1,4)-(1,3) (2,3)-(2,5)
            (0,3)-(0,5) (1,2)-(1,4) (4,5)-(2,3)
            (0,5)-(3,4) (1,4)-(0,2) (2,3)-(1,5)
            sage: G = s.automorphism_group(); G
            Permutation Group with generators [(1,3,5)(2,4,6), (1,6)(2,5)(3,4)]
            sage: G.order()
            6
            sage: sum(1 for _ in s.cylinder_diagram_iterator(up_to_isomorphism=False))
            6
        """
        cbot = self.bot_cycle_tuples()
        ctop0 = self.top_cycle_tuples()
        n = self.nseps()

        connected = not connected

        if up_to_isomorphism:
            s = set([])
            for ctop in itertools.permutations(ctop0):
                c = CylinderDiagram(zip(cbot,ctop),check=False)
                cc = c.relabel(inplace=False)
                if cc not in s:
                    s.add(cc)
                    if (connected or cc.is_connected()) and cc.smallest_integer_lengths():
                        yield cc

        else:
            for ctop in itertools.permutations(ctop0):
                c = CylinderDiagram(zip(cbot,ctop),check=False)
                if (connected or c.is_connected()) and c.smallest_integer_lengths():
                    yield c

    def cylinder_diagrams(self, connected=True,up_to_isomorphism=True):
        r"""
        Return the list of cylinder diagrams associated to this separatrix
        diagram.

        We warn that the cylinder diagram may be renumeroted in the output list
        (in order to prevent repetitions). If you care about numerotation the
        option ``up_to_isomorphism`` should be set to False.

        INPUT:

        - ``connected`` - boolean (default: True)

        - ``up_to_isomorphism`` - boolean (default: True)

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0)(1)(2)','(0)(1)(2)')
            sage: for c in s.cylinder_diagrams(connected=True): print c
            (0)-(2) (1)-(0) (2)-(1)
            sage: for c in s.cylinder_diagrams(connected=False): print c
            (0)-(0) (1)-(1) (2)-(2)
            (0)-(0) (1)-(2) (2)-(1)
            (0)-(2) (1)-(0) (2)-(1)

            sage: s = SeparatrixDiagram('(0,1)(2)','(0)(1,2)')
            sage: for c in s.cylinder_diagrams(): print c
            (0,1)-(0,2) (2)-(1)

        In the example below, there is no isomorphism problem for the cylinder
        diagram generation as the separatrix diagram admit no automorphism::

            sage: s = SeparatrixDiagram('(0,3)(1,4,5)(2)','(0)(1,2)(3,4,5)')
            sage: for c in s.cylinder_diagrams(): print c
            (0,1,2)-(0,1,5) (3,5)-(2,4) (4)-(3)
            (0,2,4)-(2,5) (1,3)-(0,1,4) (5)-(3)
            (0,3,1)-(5) (2,5)-(3,4) (4)-(0,2,1)
            sage: for c in s.cylinder_diagrams(up_to_isomorphism=False): print c
            (0,3)-(3,4,5) (1,4,5)-(1,2) (2)-(0)
            (0,3)-(1,2) (1,4,5)-(3,4,5) (2)-(0)
            (0,3)-(1,2) (1,4,5)-(0) (2)-(3,4,5)
            sage: s.automorphism_group()
            Permutation Group with generators [()]
        """
        return sorted(self.cylinder_diagram_iterator(
            connected=connected,
            up_to_isomorphism=up_to_isomorphism))

class QuadraticSeparatrixDiagram(SageObject):
    r"""
    A quadratic separatrix diagram

    sigma  : permutation around vertices
    twin   : an involution without fixed point
    cycles : a pair union of cycles

    EXAMPLES::

        sage: s = QuadraticSeparatrixDiagram([(1,1),(2,2)]); s
        (1,1)-(2,2)
        sage: s.stratum()
        Q_0(-1^4)

        sage: s = QuadraticSeparatrixDiagram([(1,1),(2,3,2,3)]); s
        (1,1)-(2,3,2,3)
        sage: s.stratum()
        Q_1(2, -1^2)
    """
    def __init__(self,bot,top):
        pass

    def stratum(self):
        pass



def cyclic_direction(x,y,z):
    r"""
    Returns 1 or -1 depending on the cyclic ordering of (x,y,z)
    """
    if x < y and y < z: return 1
    if y < z and z < x: return 1
    if z < x and x < y: return 1
    return -1

# iterators

def canonical_perm(part,i=0):
    r"""
    Return the canonical permutation with the given part

    EXAMPLES::

        sage: from sage.dynamics.flat_surfaces.separatrix_diagram import canonical_perm
        sage: canonical_perm([3,2])
        [1, 2, 0, 4, 3]
        sage: canonical_perm([2,2,2])
        [1, 0, 3, 2, 5, 4]
        sage: canonical_perm([1,1,3])
        [0, 1, 3, 4, 2]
    """
    res = []
    for p in part:
        res.extend(xrange(i+1,i+p))
        res.append(i)
        i += p
    return res

def canonical_perm_i(part,i=0):
    r"""
    Return the canonical permutation reversed

    EXAMPLES::

        sage: from sage.dynamics.flat_surfaces.separatrix_diagram import canonical_perm_i
        sage: canonical_perm_i([3,2])
        [2, 0, 1, 4, 3]
        sage: canonical_perm_i([2,2,2])
        [1, 0, 3, 2, 5, 4]
        sage: canonical_perm_i([1,1,3])
        [0, 1, 4, 2, 3]
    """
    res = []
    for p in part:
        res.append(i+p-1)
        res.extend(xrange(i,i+p-1))
        i += p
    return res

def separatrix_diagram_fast_iterator(profile,ncyls=None):
    r"""
    Iterator over separatrix diagram with given ``profile``

    Return a list of 3-tuples ``[bot, top, s]`` where ``bot`` and ``top`` are
    list on 0, ..., nseps-1 that corresponds to a separatrix diagram with
    profile ``profile`` while ``s`` is the element conjugacy class corresponding
    to the profile which equals ``bot * top``.

    If ncyls is not None, it should be a list of integers from which the number
    of cylinders is considered.

    Warning: each isomorphism class of separatrix diagram is output more than
    once in general. If you want a unique representative in each isomorphism
    class you may consider the method separatrix_diagram_iterator instead.

    EXAMPLES::

        sage: from sage.dynamics.flat_surfaces.separatrix_diagram import separatrix_diagram_fast_iterator
        sage: for s in separatrix_diagram_fast_iterator([3]): print s
        ([0, 2, 1], [1, 0, 2], [(0, 1, 2)])
        ([1, 2, 0], [1, 2, 0], [(0, 2, 1)])
        ([2, 1, 0], [1, 0, 2], [(0, 2, 1)])
        sage: for s in separatrix_diagram_fast_iterator([2,2]): print s
        ([0, 2, 3, 1], [1, 2, 0, 3], [(0, 1), (2, 3)])
        ([0, 1, 3, 2], [1, 0, 2, 3], [(0, 1), (2, 3)])
        ([1, 2, 3, 0], [1, 2, 3, 0], [(0, 2), (1, 3)])
        ([1, 3, 2, 0], [1, 2, 0, 3], [(0, 2), (1, 3)])
        ([3, 2, 1, 0], [1, 0, 3, 2], [(0, 2), (1, 3)])
        ([3, 1, 0, 2], [1, 2, 0, 3], [(0, 3), (1, 2)])
        ([2, 3, 0, 1], [1, 0, 3, 2], [(0, 3), (1, 2)])
    """
    from sage.combinat.partition import Partition,Partitions
    from sage.combinat.permutation_conjugacy_class_iterator import conjugacy_class_iterator

    part = Partition(profile)
    n = sum(part)
    d = (n+len(part))//2  # the maximum number of cylinders is known
                          # to be g+s-1 from a theorem of Y. Naveh
    res = set([])

    tops = [[]]
    if ncyls is None:
        ncyls = range(1,d+1)
    else:
        if isinstance(ncyls,(int,Integer)):
            ncyls = set([ncyls])
        else:
            ncyls = set(map(Integer,ncyls))
        for i in ncyls:
            if i < 1 or i > d:
                raise ValueError, "%d is not possible as number of cylinders"%i

    # build the list of admissible tops up to conjugacy class
    for k in xrange(1,d+1):
        tops.append([])
        if k in ncyls:
            for p in Partitions(n,length=k):
                tops[-1].append((canonical_perm(p),canonical_perm_i(p)))

    for s in conjugacy_class_iterator(part):
        for k in xrange(len(tops)):
            for top,top_i in tops[k]:
                bot = range(len(top_i))
                for cycle in s:
                    for i in xrange(len(cycle)-1):
                        bot[cycle[i]] = top_i[cycle[i+1]]
                    bot[cycle[-1]] = top_i[cycle[0]]

                seen = [True]*len(bot)
                nb_cycles = 0
                for i in xrange(len(bot)):
                    if seen[i]:
                        seen[i] = False
                        nb_cycles += 1
                        j = bot[i]
                        while seen[j]:
                            seen[j] = False
                            j = bot[j]

                if nb_cycles == k:
                    yield (bot,top,s)

def separatrix_diagram_iterator(profile,ncyls=None):
    r"""
    Iterator over separatrix diagram with given ``profile``

    Warning: to prevent isomorphism class to be output twice the function
    implement a cache mechanism. If you intend to iterate through a huge
    class of separatrix_diagram and do not care about isomorphism problem use
    separatrix_diagram_fast_iterator instead.

    EXAMPLES::

        sage: from sage.dynamics.flat_surfaces.separatrix_diagram import separatrix_diagram_iterator

        sage: for s in separatrix_diagram_iterator([1,1]): print s
        Separatrix diagram
         bot (0,1)
         top (0,1)
        Separatrix diagram
         bot (0)(1)
         top (0)(1)

        sage: for s in separatrix_diagram_iterator([3]): print s
        Separatrix diagram
         bot (0)(1,2)
         top (0,1)(2)
        Separatrix diagram
         bot (0,1,2)
         top (0,1,2)

        sage: for s in separatrix_diagram_iterator([2,2]): print s
        Separatrix diagram
         bot (0)(1,2,3)
         top (0,1,2)(3)
        Separatrix diagram
         bot (0)(1)(2,3)
         top (0,1)(2)(3)
        Separatrix diagram
         bot (0,1,2,3)
         top (0,1,2,3)
        Separatrix diagram
         bot (0,1)(2,3)
         top (0,2)(1,3)
    """
    res = set([])
    for bot,top,_ in separatrix_diagram_fast_iterator(profile,ncyls):
        s = SeparatrixDiagram(bot,top)
        s.relabel(inplace=True)
        if s not in res:
            res.add(s)
            yield s

def separatrix_diagrams(profile,ncyls=None):
    r"""
    Return the list of separatrix diagram with given ``profile``.

    EXAMPLES::

        sage: from sage.dynamics.flat_surfaces.separatrix_diagram import separatrix_diagrams
        sage: for s in separatrix_diagrams([3]): print s
        Separatrix diagram
         bot (0,1,2)
         top (0,1,2)
        Separatrix diagram
         bot (0)(1,2)
         top (0,1)(2)
    """
    res = set([])
    for bot,top,_ in separatrix_diagram_fast_iterator(profile,ncyls):
        s = SeparatrixDiagram(bot,top)
        s.relabel(inplace=True)
        res.add(s)
    return sorted(res)

#
# Cylinder diagram
#  (or completely periodic decomposition)
#

def string_to_cycle(s):
    if len(s) < 2:
        raise ValueError, "Wrong syntax"
    if s[0] != '(':
        raise ValueError, "A cycle string should start with an opening paranthesis"
    if s[-1] != ')':
        raise ValueError, "A cycle string should end with a closing paranthesis"
    return tuple(Integer(i) for i in s[1:-1].split(','))

def orientation_cover(alpha,phi,a,verbose=0):
    r"""
    Build the cylinder diagram of Abelian differentials that double covers it.

    A quadratic differrential separatrix diagram is given by three permutations

    - sigma: the permutation of 1/2-separatrices around vertices
    - alpha: the permutation of 1/2-separatrices that describe the separatrices
       (it is a fixed point free involution)
    - phi: the permutation of 1/2-separatrices that describe the cycles.
    """
    if verbose: print " orientation cover"
    cyls = []
    todo = [True]*a

    for i in xrange(a):
        if todo[i]:
            todo[i] = False
            b = [i]
            if alpha[i] >= a:
                t = [i]
            else:
                t = [alpha[i]]
            if verbose: print "  top from %d,  bot from %d"%(i,t[0])
            j = phi[i]
            if j >= a:
                j = phi[j]
            while j != i:
                todo[j] = False
                b.append(j)

                if alpha[j] >= a:
                    t.append(j)
                else:
                    t.append(alpha[j])
                if verbose: print "  add %d to bot,  add %d to top"%(j,b[-1])

                j = phi[j]
                if j >= a:
                    j = phi[j]

            cyls.append((b,t))

    return CylinderDiagram(cyls)

def hyperelliptic_cylinder_diagram_iterator(a,verbose=False):
    r"""
    Return an iterator over cylinder diagrams of Abelian differentials that
    double covers Q((a-2), -1^(a+2)).

    TODO:

    - An optimization could be obtained by considering the generation of
      k-subsets of {1,...,n} up to the cyclic symmetry of the tree.

    INPUT:

    - ``a`` - integer - angle of the conical singularity of the quadratic
      differential.

    - ``verbose`` - integer (default: 0) - output various information during the
      iteration (mainly for debug).
    """
    from sage.combinat.plane_tree import admissible_plane_tree_iterator
    from gray_codes import revolving_door_switch_iterator

    cyl_diags = set([])
    if verbose is True: verbose=1
    aa = a//2
    B = [False]*(2*a+2)  # open loops indicator
                         # if B[k] is not False, it is where loop k starts
    sigma = range(1,a) + [0] + range(a,2*a+2)
    for t,n,l in admissible_plane_tree_iterator(a):
        # Build the initial tree
        L = []
        p = 2*n-a               # the number of poles
        ll = 0                  # leaf counter
        s = 0                   # 1/2-separatrix counter
        sp = a                  # pole counter
        alpha = [None]*(2*a+2)  # edge permutation
        phi = [None]*(2*a+2)    # face permutation
        if verbose:
            print "n = %d,  l = %d,  p = %d"%(n,l,p)
            print "t =", t
        for k in xrange(1,n+2):
            if verbose: print " k = %d"%k
            for kk in xrange(t[k-1],t[k]-1,-1): # close current loops
                if not (B[kk] is False):
                    if verbose:
                        print " close loop from %d to %d"%(B[kk],s)
                    alpha[B[kk]] = s
                    alpha[s] = B[kk]
                    phi[s] = (B[kk]-1)%a
                    phi[B[kk]] = (s-1)%a
                    s += 1
                    if verbose > 2:
                        print " alpha =",alpha
                        print " phi   =",phi
            if ll < p and t[k] >= t[k+1]:
                L.append(s) # store the leaf
                # t[k] is a pole
                if verbose: print " pole at %d"%s
                alpha[s] = sp
                alpha[sp] = s
                phi[s] = sp
                phi[sp] = (s-1)%a
                s += 1
                sp += 1
                ll += 1
                B[t[k]] = False
                if verbose > 2:
                    print " alpha =",alpha
                    print " phi   =",phi

            elif k != n+1: # not at the end -> open a loop
                if t[k] >= t[k+1]: # store the leaf
                    L.append(s)
                if verbose: print " open loop at %d"%s
                B[t[k]] = s
                s += 1
                if verbose > 2:
                    print " alpha =",alpha
                    print " phi   =",phi

        if verbose:
            print " tree is over"
            print " alpha =", alpha
            print " phi =", phi

        for pp in xrange(a+p,2*a+2,2):
            if verbose: print " paired poles (%d,%d)"%(pp,pp+1)
            alpha[pp] = phi[pp] = pp+1
            alpha[pp+1] = phi[pp+1] = pp
            if verbose > 1:
                print " alpha =",alpha
                print " phi   =",phi

        assert len(L) == l, "This may not happen"

        # yield the canonical sepx. diag
        if verbose:
            print " ="*(3*a+7)
            print " sigma =",sigma
            print " alpha =",alpha
            print " phi   =",phi
            print " ="*(3*a+7)

        c = orientation_cover(alpha,phi,a,verbose=verbose)
        c.relabel(inplace=True)
        if c not in cyl_diags:
            cyl_diags.add(c)
            yield c

        # Make the poles vary among the leaves
        #TODO: optimization when tree has nontrivial cyclic symmetry
        if p != 0 and p != l:
            if verbose:
                print " start revolving door(%d,%d)"%(l,p)
                print " leaves are at separatrices", L
            for i,j in revolving_door_switch_iterator(l,p):
                i = L[i]
                j = L[j]
                if verbose > 1:
                    print " revolve i=%d j=%d"%(i,j)
                a_i = alpha[i]
                a_j = alpha[j]
                s_i = sigma[i]
                s_a_j = sigma[a_j]
                ss_i = phi[alpha[i]]  # sigma^-1(i)
                ss_j = phi[alpha[j]]  # sigma^-1(j)
                a_s_i = alpha[s_i]
                a_s_a_j = alpha[s_a_j]

                assert sigma[j] == a_j, "sigma[%d] = %d != alpha[%d]"%(j,sigma[j],a_j)
                assert phi[i] == a_i, "phi[%d] != alpha[i]"%(i,i)
                assert phi[a_s_i] == i, "phi[%d] != %d"%(a_s_i,i)
                assert phi[j] == j, "phi[%d] + %d"%(j,j)
                assert phi[a_s_a_j] == a_j, "phi[%d] != alpha[%d]"%(a_s_a_j,j)

                alpha[i]     = a_j
                alpha[a_j]   = i
                alpha[j]     = a_i
                alpha[a_i]   = j

                sigma[i]     = a_j    # old_sigma[j]
                sigma[a_j]   = s_i    # old_sigma[i]
                sigma[j]     = s_a_j  # old_sigma[a_j]

                phi[i]       = i      # old_phi[a_s_i]
                phi[j]       = a_i    # old_phi[i]

                if s_i != j: # and a_s_i == a_j
                    phi[a_s_i]   = a_j    # old_phi[a_s_a_j]
                    phi[a_i]     = ss_j   # old_phi[a_j]
                else:
                    phi[a_i] = a_j

                if s_a_j != i:
                    phi[a_j]     = ss_i   # old_phi[a_i]
                    phi[a_s_a_j] = j      # old_phi[j]
                else:
                    phi[a_j] = j

                if verbose:
                    print " ="*(3*a+7)
                    print " sigma =",sigma
                    print " alpha =",alpha
                    print " phi   =",phi
                    print " ="*(3*a+7)

                for i in xrange(2*a+2):
                    ii = phi[alpha[sigma[i]]]
                    assert ii == i, "f_a_s(%d) == %d != %d"%(i,ii,i)
                for i in xrange(a,2*a+2):
                    assert sigma[i] == i, "sigma[%d] = %d != %d"%(i,sigma[i],i)
                c = orientation_cover(alpha,phi,a,verbose=verbose)
                c.relabel(inplace=True)
                if c not in cyl_diags:
                    cyl_diags.add(c)
                    yield c

            # reinitialize sigma
            sigma = range(1,a) + [0] + range(a,2*a+2)

class CylinderDiagram(SeparatrixDiagram):
    r"""
    Separatrix diagram with pairing.

    Each cylinder is stored as a couple (bot,top) for which the orientation is
    as follows

     ---------------------
    |     <-- top --      |
    |                     |
    |     -- bot -->      |
     ---------------------

    INPUT:

    - ``data`` - list of 2-tuples - matching of bottom-top pairs

    EXAMPLES:

    We first build the simplest cylinder diagram which corresponds to a torus::

        sage: CylinderDiagram([((0,),(0,))])
        (0)-(0)

    The same initialized from a string::

        sage: CylinderDiagram('(0)-(0)')
        (0)-(0)

    The following initialize a cylinder diagram with two cylinder which gives a
    surface of genus 2 with one singularity of degree 2::

        sage: CylinderDiagram([((0,1),(0,2)),((2,),(1,))])
        (0,1)-(0,2) (2)-(1)

    ALGORITHM:

    A cylinder is represented by a couple (i,j) where i is the min in bot and j
    is the min in top. The data _top_to_cyl and _bot_to_cyl corresponds to the
    association of a separatrix to the corresponding 2-tuple. The twist
    coordinate correspond to the shift betwenn those two indices.
    """
    def __init__(self,data,check=True):
        r"""
        TESTS::

            sage: c = CylinderDiagram([((0,),(0,))])
            sage: CylinderDiagram(str(c)) == c
            True
            sage: loads(dumps(c)) == c
            True
        """
        bot = []
        top = []

        if isinstance(data,str):
            data = [(string_to_cycle(b),string_to_cycle(t)) for b,t in (w.split('-') for w in data.split(' '))]

        for b,t in data:
            bot.append(tuple(b))
            top.append(tuple(t))

        SeparatrixDiagram.__init__(self,tuple(bot),tuple(top))

        b2c = [None] * self.nseps() # bot separatrix -> cylinder (bot_min_index, top_min_index)
        t2c = [None] * self.nseps() # top separatrix -> cylinder (bot_min_index, top_min_index)
        for b,t in data:
            cyl = (min(b),min(t))
            for j in b: b2c[j] = cyl
            for j in t: t2c[j] = cyl

        self._bot_to_cyl = b2c
        self._top_to_cyl = t2c

        #from sage.misc.latex import latex
        #if latex.has_file("tikz.sty"):
            #latex.add_to_preamble('\\usepackage{tikz}')
            #latex.add_to_preamble('\\usetikzlibrary{arrows}')
            #latex.add_to_jsmath_avoid_list('\\begin{tikzpicture}')

    def _repr_(self):
        r"""
        String representation

        TESTS::

            sage: c = CylinderDiagram([((0,1),(1,2)),((2,),(0,))])
            sage: repr(c) #indirect doctest
            '(0,1)-(1,2) (2)-(0)'
        """
        l = []
        for b,t in self.cylinders():
            l.append('(' + ','.join(map(str,b)) + ')-(' + ','.join(map(str,t)) + ')')
        return ' '.join(l)

    def __cmp__(self,other):
        r"""
        Comparison

        TESTS::

            sage: s3_0 = SeparatrixDiagram([0,1,2],[0,1,2])
            sage: S3 = [s3_0]

            sage: s2_0 = SeparatrixDiagram([1,0,2],[1,0,2])
            sage: s2_1 = SeparatrixDiagram([1,0,2],[2,1,0])
            sage: s2_2 = SeparatrixDiagram([1,0,2],[0,2,1])
            sage: s2_3 = SeparatrixDiagram([2,1,0],[1,0,2])
            sage: s2_4 = SeparatrixDiagram([2,1,0],[2,1,0])
            sage: s2_5 = SeparatrixDiagram([2,1,0],[0,2,1])
            sage: s2_6 = SeparatrixDiagram([0,2,1],[1,0,2])
            sage: s2_7 = SeparatrixDiagram([0,2,1],[2,1,0])
            sage: s2_8 = SeparatrixDiagram([0,2,1],[0,2,1])
            sage: S2 = [s2_0,s2_1,s2_2,s2_3,s2_4,s2_5,s2_6,s2_7,s2_8]

            sage: s1_0 = SeparatrixDiagram([1,2,0],[1,2,0])
            sage: s1_1 = SeparatrixDiagram([1,2,0],[2,0,1])
            sage: s1_2 = SeparatrixDiagram([2,0,1],[1,2,0])
            sage: s1_3 = SeparatrixDiagram([2,0,1],[2,0,1])
            sage: S1 = [s1_0,s1_1,s1_2,s1_3]

            sage: all(cmp(s,s) == 0 for s in S1+S2+S3)
            True
            sage: all(cmp(s1,s2) == -1 for s1 in S1 for s2 in S2)
            True
            sage: all(cmp(s2,s1) == 1 for s1 in S1 for s2 in S2)
            True
            sage: all(cmp(s2,s3) == -1 for s2 in S2 for s3 in S3)
            True
            sage: all(cmp(s3,s2) == 1 for s2 in S2 for s3 in S3)
            True
            sage: all(cmp(s1,s3) == -1 for s1 in S1 for s3 in S3)
            True
            sage: all(cmp(s3,s1) == 1 for s1 in S1 for s3 in S3)
            True
        """
        if not isinstance(other, CylinderDiagram):
            raise NotImplemented

        test = SeparatrixDiagram.__cmp__(self,other)
        if test: return test

        test = cmp(self._bot_to_cyl,other._top_to_cyl)
        if test: return test

        return 0

    #
    # access to attribute
    #


    def to_directed_graph(self):
        r"""
        Return a labeled directed graph that encodes the cylinder diagram.

        EXAMPLES::

            sage: c = CylinderDiagram('(0,1,5)-(2,5) (2)-(0,1,3) (3,4)-(4)'); c
            (0,1,5)-(2,5) (2)-(0,1,3) (3,4)-(4)
            sage: G = c.to_directed_graph(); G
            Looped multi-digraph on 6 vertices
            sage: G.edges()
            [(0, 1, 'b'), (0, 1, 't'), (0, 2, 'c'), (0, 5, 'c'), (1, 2, 'c'), (1, 3, 't'), (1, 5, 'b'), (1, 5, 'c'), (2, 0, 'c'), (2, 1, 'c'), (2, 2, 'b'), (2, 3, 'c'), (2, 5, 't'), (3, 0, 't'), (3, 4, 'b'), (3, 4, 'c'), (4, 3, 'b'), (4, 4, 'c'), (4, 4, 't'), (5, 0, 'b'), (5, 2, 'c'), (5, 2, 't'), (5, 5, 'c')]
        """
        G = SeparatrixDiagram.to_directed_graph(self)

        for cb,ct in self.cylinders():
            for i in cb:
                for j in ct:
                    G.add_edge(i,j,'c')

        return G

    def relabel(self, inplace=True, return_map=False):
        r"""
        Return a cylinder diagram with canonical labels.

        EXAMPLES::

            sage: import itertools
            sage: for p in itertools.permutations([0,1,2,3]):
            ...      c = CylinderDiagram([((p[0],),(p[1],)),((p[1],p[2]),(p[0],p[3])),((p[3],),(p[2],))])
            ...      cc,m = c.relabel(inplace=False, return_map=True)
            ...      b  = c.bot() ; t  = c.top()
            ...      bb = cc.bot(); tt = cc.top()
            ...      print cc
            ...      print all(bb[m[i]] == m[b[i]] for i in xrange(c.nseps())),
            ...      print all(tt[m[i]] == m[t[i]] for i in xrange(c.nseps()))
            (0,1)-(2,3) (2)-(1) (3)-(0)
            True True
            (0,1)-(2,3) (2)-(1) (3)-(0)
            True True
            (0,1)-(2,3) (2)-(1) (3)-(0)
            True True
            (0,1)-(2,3) (2)-(1) (3)-(0)
            True True
            (0,1)-(2,3) (2)-(1) (3)-(0)
            True True
            (0,1)-(2,3) (2)-(1) (3)-(0)
            True True
            (0,1)-(2,3) (2)-(1) (3)-(0)
            True True
            (0,1)-(2,3) (2)-(1) (3)-(0)
            True True
            ...
            (0,1)-(2,3) (2)-(1) (3)-(0)
            True True

            sage: import itertools
            sage: for p in itertools.permutations([0,1,2,3,4,5]):
            ...      c1 = ((p[0],p[4]),(p[0],p[3]))
            ...      c2 = ((p[1],p[3]),(p[1],p[5]))
            ...      c3 = ((p[2],p[5]),(p[2],p[4]))
            ...      c = CylinderDiagram([c1,c2,c3])
            ...      cc,m = c.relabel(inplace=False, return_map=True)
            ...      b  = c.bot() ; t  = c.top()
            ...      bb = cc.bot(); tt = cc.top()
            ...      print cc
            ...      print all(bb[m[i]] == m[b[i]] for i in xrange(c.nseps())),
            ...      print all(tt[m[i]] == m[t[i]] for i in xrange(c.nseps()))
            (0,5)-(0,4) (1,4)-(1,3) (2,3)-(2,5)
            True True
            (0,5)-(0,4) (1,4)-(1,3) (2,3)-(2,5)
            True True
            (0,5)-(0,4) (1,4)-(1,3) (2,3)-(2,5)
            True True
            (0,5)-(0,4) (1,4)-(1,3) (2,3)-(2,5)
            True True
            (0,5)-(0,4) (1,4)-(1,3) (2,3)-(2,5)
            True True
            ...
            (0,5)-(0,4) (1,4)-(1,3) (2,3)-(2,5)
            True True
            (0,5)-(0,4) (1,4)-(1,3) (2,3)-(2,5)
            True True

        TESTS::

            sage: c = CylinderDiagram('(0,1)-(0,2) (3,5,4)-(1,4,6) (2,6)-(3,5)')
            sage: c is c.relabel(inplace=False)
            False
            sage: c.relabel(inplace=False) is c.relabel(inplace=False)
            True
            sage: c.relabel(inplace=False).relabel(inplace=False) is c.relabel(inplace=False)
            True
        """
        if not hasattr(self,'_normal_form'):
            G = self.to_directed_graph()
            _,m = G.canonical_label(certify=True,edge_labels=True)
            # m = [m[i] for i in xrange(self.nseps())]
            # GG the new digraph
            # m from the digraph to its canonic labels
            cyls = []
            for b,t in self.cylinders():
                cyls.append((tuple(m[i] for i in b),tuple(m[i] for i in t)))

            self._normal_form = CylinderDiagram(cyls,check=False)
            self._normal_labels = m

            self._normal_form._normal_form = self._normal_form
            self._normal_form._normal_labels = range(self.nseps())

        if inplace:
            self.__dict__ = self._normal_form.__dict__

        if return_map:
            return self._normal_form, self._normal_labels
        elif not inplace:
            return self._normal_form

    def separatrix_diagram(self):
        r"""
        Return the underlying separatrix diagram

        EXAMPLES::

            sage: s = SeparatrixDiagram('(0,1)(2,3,4)','(0,3)(1,4,2)'); s
            Separatrix diagram
             bot (0,1)(2,3,4)
             top (0,3)(1,4,2)
            sage: c = s.to_cylinder_diagram([(0,1),(1,0)]); c
            (0,1)-(1,4,2) (2,3,4)-(0,3)
            sage: c.separatrix_diagram() == s
            True
        """
        return SeparatrixDiagram(self._bot,self._top,check=False)

    def cylinders(self):
        r"""
        Cylinders of self

        EXAMPLES::

            sage: c = CylinderDiagram('(0,2,4)-(1,3,5) (1,5)-(0) (3)-(2,4)')
            sage: c
            (0,2,4)-(1,3,5) (1,5)-(0) (3)-(2,4)
            sage: c.cylinders()
            [((0, 2, 4), (1, 3, 5)), ((1, 5), (0,)), ((3,), (2, 4))]
        """
        return [(b,self.top_orbit(self._bot_to_cyl[b[0]][1])) for b in self.bot_cycle_tuples()]

    def bot_to_cyl(self, j):
        r"""
        Return the cylinder above the separatrix j

        EXAMPLES::

            sage: c = CylinderDiagram('(0,2,4)-(1,3,5) (1,5)-(0) (3)-(2,4)')
            sage: c
            (0,2,4)-(1,3,5) (1,5)-(0) (3)-(2,4)
            sage: c.bot_to_cyl(0)
            ((0, 2, 4), (1, 3, 5))
            sage: c.bot_to_cyl(1)
            ((1, 5), (0,))
            sage: c.bot_to_cyl(3)
            ((3,), (2, 4))
        """
        jb,jt = self._bot_to_cyl[j]
        return self.bot_orbit(jb), self.top_orbit(jt)

    def top_to_cyl(self, j):
        r"""
        Return the cylinder below the separatrix j
        """
        jb,jt = self._top_to_cyl[j]
        return self.bot_orbit(jb), self.top_orbit(jt)

    #
    # properties
    #

    def is_connected(self):
        r"""
        Check the connectedness of this cylinder diagram.

        TESTS::

            sage: CylinderDiagram('(0)-(1) (1)-(0)').is_connected()
            True
            sage: CylinderDiagram('(0,1)-(0) (2)-(1,2)').is_connected()
            True

            sage: CylinderDiagram('(0)-(0) (1)-(1)').is_connected()
            False
            sage: CylinderDiagram('(0,1)-(3) (2)-(2) (3)-(0,1)').is_connected()
            False
        """
        from sage.graphs.graph import Graph
        G = Graph()
        for b,t in self.cylinders():
            G.add_edges((b[0],b[j]) for j in xrange(1,len(b)))
            G.add_edges((t[0],t[j]) for j in xrange(1,len(t)))
            G.add_edge(b[0],t[0])
        return G.num_verts() == self.nseps() and G.is_connected()

    #
    # symmetries
    #

    def inverse(self):
        r"""
        Return the inverse cylinder diagram

        The inverse of a cylinder diagram is the cylinder diagram in which all
        cylinders have been reversed. It corresponds to the multiplication by
        `-1` on the underlying Abelian differential.

        Combinatorially the operation is b0-t0 ... bk-tk becomes t0-b0 ... tk-bk

        EXAMPLES::

            sage: c = CylinderDiagram('(0,1)-(0,2) (3,5,4)-(1,4,6) (2,6)-(3,5)')
            sage: c
            (0,1)-(0,2) (2,6)-(3,5) (3,5,4)-(1,4,6)
            sage: c.inverse()
            (0,2)-(0,1) (1,4,6)-(3,5,4) (3,5)-(2,6)

        The inversion is an involution on cylinder diagrams::

            sage: all(cc.inverse().inverse() == cc for cc in AbelianStratum(4).cylinder_diagrams()) # long time
            True
        """
        return CylinderDiagram([(t,b) for (b,t) in self.cylinders()])

    def vertical_symmetry(self):
        r"""
        Return the cylinder diagram obtained by reflecting the cylinder
        configuration along the vertical axis.

        EXAMPLES::

            sage: c = CylinderDiagram('(0,3,4)-(0,3,5) (1,2,5)-(1,2,4)')
            sage: c.vertical_symmetry()
            (0,3,5)-(0,3,4) (1,2,4)-(1,2,5)
            sage: A = AbelianStratum(2,2)
            sage: all(c.vertical_symmetry().stratum() == A for c in A.cylinder_diagrams())
            True
        """
        return CylinderDiagram(tuple((b[::-1],t[::-1]) for b,t in self.cylinders()))

    def horizontal_symmetry(self):
        r"""
        Return the cylinder diagram obtained by reflecting the cylinder
        configuration along the horizontal axis.

        EXAMPLES::

            sage: c = CylinderDiagram('(0,3,4)-(0,3,5) (1,2,5)-(1,2,4)')
            sage: c.horizontal_symmetry()
            sage: A = AbelianStratum(2,2)
            sage: all(c.horizontal_symmetry().stratum() == A for c in A.cylinder_diagrams())
            True
        """
        return CylinderDiagram(tuple((t[::-1],b[::-1]) for b,t in self.cylinders()))

    def automorphism_group(self, order=False):
        r"""
        Return the automorphism group

        INPUT:

        - ``order`` - boolean (default: False) - whether or not return the order
          of the group
        """
        return self.to_directed_graph().automorphism_group(edge_labels=True, order=order)

    # TODO
    def quotient(self, H=None):
        r"""
        Quotient by a subgroup.
        """
        if H is None:
            H = self.automorphism_group()

        elif not H.is_subgroup(self.automorphism_group()):
            raise ValueError, "H must be a subgroup of the automorphism group"

        e = H.orbits()
        e_inv = [None]*self.nseps()  # sep -> new nb
        for i in xrange(len(e)):
            for j in c:
                e_inv[j-1] = i

        raise NotImplementedError

    def is_hyperelliptic(self,verbose=False):
        r"""
        Test of hyperellipticity

        Each stratum of Abelian differentials as up to three connected
        components. For the strata H(2g-2) and H(g-1,g-1) there is a special
        component called *hyperelliptic* in which all translation surfaces
        `(X,\omega)` in that component are such that `X` is hyperelliptic.

        This function returns True if and only if the cylinder diagrams
        correspond to a decomposition of a surface associated to the
        hyperelliptic components in H(2g-2) or H(g-1,g-1).

        EXAMPLES:

        In genus 2, strata H(2) and H(1,1), all surfaces are hyperelliptic::

            sage: for c in AbelianStratum(2).cylinder_diagrams():
            ...      print c
            ...      print c.is_hyperelliptic()
            (0,1,2)-(0,1,2)
            True
            (0,1)-(0,2) (2)-(1)
            True

            sage: for c in AbelianStratum(1,1).cylinder_diagrams():
            ...      print c
            ...      print c.is_hyperelliptic()
            (0,1,3,2)-(0,1,3,2)
            True
            (0,1,2)-(0,1,3) (3)-(2)
            True
            (0,3)-(0,2) (1,2)-(1,3)
            True
            (0,1)-(2,3) (2)-(1) (3)-(0)
            True

        In higher genera, some of them are, some of them are not::

            sage: C = AbelianStratum(4).cylinder_diagrams()
            sage: len(C)
            22
            sage: len(filter(lambda c: c.is_hyperelliptic(), C))
            5

            sage: C = AbelianStratum(2,2).cylinder_diagrams()
            sage: len(C)
            50
            sage: len(filter(lambda c: c.is_hyperelliptic(), C))
            12
        """
        z = self.stratum().zeros()
        if z == [0] or z == [2] or z == [1,1]: return True
        if 0 in z:
            raise NotImplementedError, "is_hyperelliptic method not implemented for cylinder diagrams with fake zeros"
        ns = self.nseps()
        if len(z) == 1: # minimal stratum H(2g-2)
            for cy in self.cylinders():
                if len(cy[0]) != len(cy[1]): return False
            b = self.bot()
            t = self.top()
            # build list of seps in cyclic order around zero, starting by outgoing sep 0
            lout = [0]
            lin = []
            for _ in xrange(ns):
                lin.append(t[lout[-1]])
                lout.append(b[lin[-1]])
            if verbose: print 'lin  ', lin; print 'lout', lout
            # build involution on separatrices
            p = [None]*ns
            for a in xrange(ns):
                p[lout[a]] = lin[(a+ns//2)%ns]
            if verbose: print "involution on seps", p
            # wsep = counter of sepatrices with a wpt
            wsep = 0
            for cy in self.cylinders():
                for k in cy[0]:
                    # check that p(k) is on the top of the cyl that has k on its bottom
                    if p[k] not in cy[1]: return False
                    # check that if k is on bot and top of cyl, then p(k) = k
                    if k in cy[1]:
                        if k != p[k]: return False
                        wsep += 1
            if verbose: print "wsep", wsep
            # check number of w pts
            if wsep + 2*self.ncyls() != z[0] + 3: return False
            # check that cylinders are stable under involution
            if self != CylinderDiagram(
                [(cy[0],tuple(map(lambda x: p[x],cy[0]))) for cy in self.cylinders()]):
                return False
            return True
        elif len(z) == 2: # should be stratum H(g-1,g-1)
            if z[0] != z[1]: return False
            for cy in self.cylinders():
                if len(cy[0]) != len(cy[1]): return False
            b = self.bot()
            t = self.top()
            # build list of seps in cyclic order around first zero, starting by outgoing sep 0
            lout = [0]
            lin = []
            for _ in xrange(ns//2):
                lin.append(t[lout[-1]])
                lout.append(b[lin[-1]])
            if verbose: print 'lin  ', lin; print 'lout', lout
            # build list of seps in cyclic order around the other zero
            a = 0
            while a in lout: a += 1
            llout = [a]
            llin = []
            for _ in xrange(ns//2):
                llin.append(t[llout[-1]])
                llout.append(b[llin[-1]])
            if verbose: print 'llin  ', llin; print 'llout', llout
            # now, try each way the involution could send lout to llout
            for j in xrange(ns//2):
                test = True
                # build involution on separatrices
                p = [None]*ns
                for a in xrange(ns//2):
                    p[lout[a]] = llin[(j+a)%(ns//2)]
                    p[llout[a]] = lin[(a-j-1)%(ns//2)]
                if verbose: print "involution on seps", p
                wsep = 0
                for cy in self.cylinders():
                    for k in cy[0]:
                        # check that p(k) is on the top of the cyl that has k on its bottom
                        if p[k] not in cy[1]:
                            test = False
                            break
                        # check that if k is on bot and top of cyl, then p(k) = k
                        if k in cy[1]:
                            if k != p[k]:
                                test = False
                                break
                            wsep += 1
                    if test is False: break
                if test is False: continue # try next j
                if verbose: print "wsep", wsep
                # check number of w pts
                if wsep + 2*self.ncyls() != 2*z[0] + 4:
                    continue # try next j
                # check that cylinders are stable under involution
                if self != CylinderDiagram(
                    [(cy[0],tuple(map(lambda x: p[x],cy[0]))) for cy in self.cylinders()]):
                    continue # try next j
                return True
            return False

        else:
            return False

    #
    # construction
    #

    def dual_graph(self):
        r"""
        The dual graph of the stable curve at infinity in the horizontal
        direction.

        This graph is defines as follows. Cut each horizontal cylinder along a
        circumference, then the vertices are the equivalence class of half
        cylinder modulo the relation "linked by a saddle connection" and the
        edges are the circumferences.

        EXAMPLES:

        We consider the three diagrams of the stratum H(1,1)::

            sage: c1 = CylinderDiagram('(0,1,2,3)-(0,1,2,3)')
            sage: c1.stratum()
            H_2(1^2)
            sage: c1.dual_graph()
            Looped multi-graph on 1 vertex
            sage: c2 = CylinderDiagram('(0,1)-(1,2) (2,3)-(0,3)')
            sage: c2.stratum()
            H_2(1^2)
            sage: c2.dual_graph()
            Looped multi-graph on 1 vertex
            sage: c3 = CylinderDiagram('(0,1)-(2,3) (2)-(0) (3)-(1)')
            sage: c3.stratum()
            H_2(1^2)
            sage: c3.dual_graph()
            Looped multi-graph on 2 vertices
        """
        from sage.graphs.graph import Graph
        cb = self.bot_cycle_tuples()
        ct = self.top_cycle_tuples()

        # first compute the equivalence class of half cylinders (i.e. gives vertices)
        V = Graph()
        V.add_vertices('%db' %c[0] for c in cb)
        V.add_vertices('%dt' %c[0] for c in ct)

        for i in xrange(self.nseps()):
            V.add_edge(
                ('%db' %self._bot_to_cyl[i][0]),
                ('%dt' %self._top_to_cyl[i][1]))

        # the dual graph
        G = Graph(loops=True,multiedges=True)
        cc = map(tuple,V.connected_components())
        hc2cc = {} # half-cyl to conn comp
        for c in cc:
            for e in c:
                hc2cc[e] = c
        for c in self.cylinders():
            G.add_edge(hc2cc['%db' %c[0][0]],hc2cc['%dt' %c[1][0]],(c[0][0],c[1][0]))
        return G

    def matrix_relation(self):
        r"""
        Return the matrix of relation on lengths.

        The output matrix has size ncyls x nseps and
         m[cyl,sep] = (1 if sep in top of cyl 0 else) - (1 if sep in bot of cyl 0 else)
        """
        from sage.matrix.constructor import matrix

        m = matrix(self.ncyls(),self.nseps(),sparse=True)
        for i,(top,bot) in enumerate(self.cylinders()):
            for t in top:
                m[i,t] = 1
            for b in bot:
                m[i,b] += -1
        return m


    #
    # Abelian differentials / coordinates
    #

    def stratum_component(self):
        r"""
        Return the connected component of stratum in which the cylinder diagram
        belongs.

        EXAMPLES::

            sage: CylinderDiagram('(0,1)-(0,2) (2)-(1)').stratum_component()
            H_2(2)^hyp

            sage: c = CylinderDiagram('(0,3,2,1)-(1,4,3,2) (4,7,6,5)-(0,7,6,5)')
            sage: c.stratum_component()
            H_4(3^2)^hyp
            sage: c = CylinderDiagram('(0,1,4)-(1,6,7) (2,5,3)-(0,2,4) (6)-(5) (7)-(3)')
            sage: c.stratum_component()
            H_4(3^2)^nonhyp

            sage: c = CylinderDiagram('(0,6)-(1,7) (1,5,4,3,2)-(2,6,5,4,3) (7,9,8)-(0,9,8)')
            sage: c.stratum_component()
            H_5(4^2)^hyp
            sage: c = CylinderDiagram('(0,2,6,1)-(0,8,1,9,2,5,7,4) (3,7,4,8,9,5)-(3,6)')
            sage: c.stratum_component()
            H_5(4^2)^even
            sage: c = CylinderDiagram('(3,7,4,8,9,5)-(0,8,1,9,2,5,7,4) (0,2,6,1)-(3,6)')
            sage: c.stratum_component()
            H_5(4^2)^odd
        """
        stratum = self.stratum()
        cc = stratum._cc
        if len(cc) == 1:
            return cc[0](stratum)

        from abelian_strata import HypASC, OddASC, EvenASC, NonHypASC

        if cc[0] is HypASC:
            if self.is_hyperelliptic():
                return HypASC(stratum)
            elif len(cc) == 2:
                return cc[1](stratum)

        if self.spin_parity() == 0:
            return EvenASC(stratum)
        else:
            return OddASC(stratum)

    @cached_method
    def smallest_integer_lengths(self):
        r"""
        Check if there is a integer solution that satisfy the cylinder
        condition.

        If there is a solution, the function returns a list of lengths for
        separatrices that minimize the sum of lengths. Otherwise, returns False.

        EXAMPLES::

            sage: c = CylinderDiagram('(0,1)-(0,2) (2,3)-(1,3)')
            sage: c.smallest_integer_lengths()
            (4, [1, 1, 1, 1])
            sage: c = CylinderDiagram('(0,1,2)-(3) (3)-(0) (4)-(1,2,4)')
            sage: c.smallest_integer_lengths()
            False

            sage: c = CylinderDiagram('(0,1)-(0,5) (2)-(3) (3,6)-(8) (4,8)-(6,7) (5)-(2,4) (7)-(1)')
            sage: c.smallest_integer_lengths()
            (13, [1, 2, 1, 1, 1, 2, 1, 2, 2])
        """
        try:
            from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
        except ImportError:
            raise ImportError, "Must have MixedIntegerLinearProgram... install GLPK"

        n = self.nseps()
        bot = self.bot_cycle_tuples()
        top = [self.top_orbit(self._bot_to_cyl[b[0]][1]) for b in bot]

        p = MixedIntegerLinearProgram(maximization=False)
        scl = p.new_variable()
        p.set_objective(sum(scl[i] for i in xrange(n)))
        for i in xrange(n):
            p.add_constraint(scl[i],min=1)
        for b,t in itertools.izip(bot,top):
            p.add_constraint(
                    sum(scl[i] for i in set(b).difference(t)) ==
                    sum(scl[i] for i in set(t).difference(b))
                    )

        try:
            total = Integer(p.solve())
            lengths = [Integer(p.get_values(scl[i])) for i in xrange(n)]
            return total, lengths
        except MIPSolverException:
            return False

    #
    # homology
    #

    def to_ribbon_graph(self):
        r"""
        Return a ribbon graph

        A *ribbon graph* is a graph embedded in an oriented surface such that
        its complement is a union of topological discs. To a cylinder diagram we
        associate the graph which consists of separatrices together with a
        choice of one vertical edge in each cylinder.

        The edges of the ribbon graph are labeled by ``(i,nseps+i)`` for
        separatrices and by ``(2(nseps+j),2(nseps+j)+1)`` for vertical in
        cylinders.

        EXAMPLES::

            sage: C = CylinderDiagram([((0,1),(0,2)),((2,),(1,))])
            sage: C.stratum()
            H_2(2)
            sage: R = C.to_ribbon_graph(); R
            Ribbon graph with 1 vertex, 5 edges and 2 faces
            sage: l,m = R.cycle_basis(intersection=True)
            sage: m.rank() == 2 * C.genus()
            True

        TESTS::

            sage: f = lambda c: c.to_ribbon_graph().cycle_basis(intersection=True)[1]

            sage: a = AbelianStratum(2)
            sage: all(f(c).rank() == 4 for c in a.cylinder_diagrams())
            True
            sage: a = AbelianStratum(1,1)
            sage: all(f(c).rank() == 4 for c in a.cylinder_diagrams())
            True
        """
        from homology import RibbonGraphWithAngles

        n = self.nseps()
        m = self.ncyls()

        edges = [(i,n+i) for i in xrange(n)] + [(2*(n+i),2*(n+i)+1) for i in xrange(m)]
        faces = []
        angles = [1] * (2*(n+m))
        half = Integer(1)/Integer(2)
        for j,(b,t) in enumerate(self.cylinders()):
            face = [i for i in b] + [2*(n+j)] + [n+i for i in t] + [2*(n+j)+1]
            faces.append(tuple(face))
            angles[b[0]] = angles[n+t[0]] = angles[2*(n+j)] = angles[2*(n+j)+1] = half

        return RibbonGraphWithAngles(edges=edges,faces=faces,angles=angles)

    def spin_parity(self):
        r"""
        Return the spin parity of any surface that is built from this cylinder
        diagram.

        EXAMPLES::

            sage: c = CylinderDiagram('(0,1,2,3,4)-(0,1,2,3,4)')
            sage: c.spin_parity()
            0
            sage: c = CylinderDiagram('(0,1,2,3,4)-(0,1,4,2,3)')
            sage: c.spin_parity()
            1

            sage: c = CylinderDiagram('(0,2,6,1)-(0,8,1,9,2,5,7,4) (3,7,4,8,9,5)-(3,6)')
            sage: c.spin_parity()
            0
            sage: c = CylinderDiagram('(3,7,4,8,9,5)-(0,8,1,9,2,5,7,4) (0,2,6,1)-(3,6)')
            sage: c.spin_parity()
            1
        """
        if any(z%2 for z in self.stratum().zeros()):
            return None
        return self.to_ribbon_graph().spin_parity()

    def circumferences_of_cylinders(self,ring=None):
        r"""
        Return the set of circumferences of cylinders as cycles in the chain
        space.
        """
        from sage.all import vector
        if ring is None:
            from sage.rings.integer_ring import ZZ
            ring = ZZ

        g = self.to_ribbon_graph()
        C = g.chain_complex(ring)
        C1 = C.chain_space(1)
        Z1 = C.cycle_space(1)
        n = g.num_edges()

        l = []
        for (b,t) in self.cylinders():
            l.append(Z1(vector(ring,n,dict(g.dart_to_edge(i,orientation=True) for i in b))))
        return l

    #
    # build one or many origamis
    #

    def an_origami(self, verbose=False):
        r"""
        Return one origami with this diagram cylinder if any.
        """
        res = self.smallest_integer_lengths()
        if res is False:
            return Fasle
        m,lengths = res

        if verbose:
            print "objective value", objective
            for i, v in values.iteritems():
                print 'x_%s = %s' % (i, int(round(v)))

        widths = [sum(lengths[i] for i in bot) for bot in self.bot_cycle_tuples()]
        areas = [widths[i] for i in xrange(self.ncyls())]

        v = [0]
        for a in areas:
            v.append(v[-1] + a)

        # initialization of bottom squares: sep_i -> bottom position
        sep_bottom_pos = [None] * self.nseps()
        for i,(bot,_) in enumerate(self.cylinders()):
            w = 0
            for j in bot:
                sep_bottom_pos[j] = v[i] + w
                w += lengths[j]

        # initialization of sigma_h which remains constant
        lx = range(1, v[-1]+1)
        for i in xrange(self.ncyls()):
            for j in xrange(v[i], v[i+1], widths[i]):
                lx[j+widths[i]-1] = j

        # initialization of y except the top 
        ly = []
        for i in xrange(self.ncyls()):
            ly.extend([None]*widths[i])

        # build the top interval without twist
        for i,(_,top_seps) in enumerate(self.cylinders()):
            top = []
            for k in reversed(top_seps):
                top.extend(range(sep_bottom_pos[k],sep_bottom_pos[k]+lengths[k]))
            ly[v[i+1]-widths[i]:v[i+1]] = top

        # yield the origami without twist
        from origamis.origami import Origami_dense
        return Origami_dense(tuple(lx), tuple(ly))

    def origami_iterator(self,n):
        r"""
        Iteration over all origamis with n squares.
 
        INPUT:

        - ``n`` - positive integer - the number of squares
        """
        for w,h in self.widths_and_heights_iterator(n):
            for o in self.cylcoord_to_origami_iterator(w, h):
                yield o

    def origamis(self,n=None):
        r"""
        Return the set of origamis.

        If ``n`` is None then return the origamis with less number of squares.
        """
        if n is None:
            res = self.smallest_integer_lengths()
            if res is False:
                return False
            n = res[0]

        return list(self.origami_iterator(n))

    def widths_and_heights_iterator(self, n):
        """
        OUTPUT:

        -  ``l`` - the lengths (as many as separatrices)

        - ``h`` - the heights (as many as cylinders)
        """
        from sage.combinat.partition import OrderedPartitions
        from sage.rings.arith import divisors
        from sage.all import vector

        m = self.matrix_relation()

        min_lengths = [1] * self.nseps()
        for i in xrange(self.ncyls()):
            pos = m.nonzero_positions_in_row(i)
            pos_m = filter(lambda j: m[i,j] == -1, pos)
            pos_p = filter(lambda j: m[i,j] == 1, pos)
            if len(pos_m) == 1:
                min_lengths[pos_m[0]] = max(min_lengths[pos_m[0]], len(pos_p))
            if len(pos_p) == 1:
                min_lengths[pos_p[0]] = max(min_lengths[pos_m[0]], len(pos_m))

        min_widths = []
        for bot,top in self.cylinders():
            min_widths.append(max(
                sum(min_lengths[j] for j in top),
                sum(min_lengths[j] for j in bot)))


        for a in itertools.ifilter(
              lambda x: all(x[i] >= min_widths[i] for i in xrange(self.ncyls())),
              OrderedPartitions(n, self.ncyls())):
            area_div = tuple(filter(lambda d: d >= min_widths[i],divisors(a[i])) for i in xrange(self.ncyls()))
            for w in itertools.product(*area_div):
                h = [Integer(a[i]/w[i]) for i in xrange(self.ncyls())]

                # from here the resolution becomes linear and convex ...
                #TODO: program a linear and convex solution
                seps_b = [c[0] for c in self.cylinders()]
                nseps_b = map(len, seps_b)
                lengths = tuple(OrderedPartitions(w[i], nseps_b[i]) for i in xrange(self.ncyls()))
                for l_by_cyl in itertools.product(*lengths):
                    l = vector([0]*self.nseps())
                    for i in xrange(self.ncyls()):
                        for j in xrange(nseps_b[i]):
                            l[seps_b[i][j]] = l_by_cyl[i][j]
                    if not m*l:
                        yield l,h

    def cylcoord_to_origami_iterator(self, lengths, heights):
        r"""
        Convert coordinates of the cylinders into an origami.

        INPUT:

        - ``lengths`` - lengths of the separatrices

        - ``heights`` - heights of the cylinders

        OUTPUT:

        - iterator over all possible origamis with those lengths and heights...
        """
        from origamis.origami import Origami_dense

        _VERBOSE = False # for debug

        widths = [sum(lengths[i] for i in bot) for bot in self.bot_cycle_tuples()]
        areas = [heights[i]*widths[i] for i in xrange(self.ncyls())]

        if _VERBOSE:
            print "areas of cylinders", areas

        # intialization of partial volumes: the set of squares in cylinder i is range(v[i],v[i+1])
        v = [0]
        for a in areas:
            v.append(v[-1] + a)
        if _VERBOSE:
            print "partial volumes", v

        # initialization of bottom squares: sep_i -> bottom position
        sep_bottom_pos = [None] * self.nseps()
        for i,(bot,_) in enumerate(self.cylinders()):
            w = 0
            for j in bot:
                sep_bottom_pos[j] = v[i] + w
                w += lengths[j]

        if _VERBOSE:
            print "sep_bottom_pos", sep_bottom_pos

        # initialization of sigma_h which remains constant
        lx = range(1, v[-1]+1)
        for i in xrange(self.ncyls()):
            for j in xrange(v[i], v[i+1], widths[i]):
                lx[j+widths[i]-1] = j

        if _VERBOSE:
            print "permutation x", lx

        # initialization of y except the top 
        ly = []
        for i in xrange(self.ncyls()):
            ly.extend(range(v[i]+widths[i],v[i+1]))
            ly.extend([None]*widths[i])

        if _VERBOSE:
            print "permutation y without gluings", ly

        # build the top interval without twist
        for i,(_,top_seps) in enumerate(self.cylinders()):
            top = []
            for k in reversed(top_seps):
                top.extend(range(sep_bottom_pos[k],sep_bottom_pos[k]+lengths[k]))
            ly[v[i+1]-widths[i]:v[i+1]] = top

        # yield the one without twist
        yield Origami_dense(tuple(lx), tuple(ly))

        # yield the others using GrayCodeSwitch
        for i,o in GrayCodeSwitch(widths):
            if _VERBOSE: print i,o
            if o == 1:
                ly.insert(v[i+1]-widths[i],ly.pop(v[i+1]-1))
            else:
                ly.insert(v[i+1]-1,ly.pop(v[i+1]-widths[i]))
            if _VERBOSE: print ly
            yield Origami_dense(tuple(lx),tuple(ly))

    def cylcoord_to_origami(self, lengths, heights, twists=None):
        r"""
        Convert coordinates of the cylinders into an origami.

        INPUT:

        - ``lengths`` - lengths of the separatrices

        - ``heights`` - heights of the cylinders

        - ``twists`` - twists for cylinders


        EXAMPLES::

            sage: c = CylinderDiagram([((0,1),(1,2)),((2,),(0,))])
            sage: c.stratum()
            H(2)
            sage: c.cylcoord_to_origami([1,1,1],[1,1]).stratum()
            H(2)
            sage: o1 = c.cylcoord_to_origami([2,1,2],[1,1],[1,0])
            sage: o1 = o1.standard_form()
            sage: o2 = c.cylcoord_to_origami([2,1,2],[1,1],[0,1])
            sage: o2 = o2.standard_form()
            sage: o3 = c.cylcoord_to_origami([2,1,2],[1,1],[1,1])
            sage: o3 = o3.standard_form()
            sage: all(o.stratum() == AbelianStratum(2) for o in [o1,o2,o3])
            True
            sage: o1 == o2 or o1 == o3 or o3 == o1
            False

        If the lengths are not compatible with the cylinder diagram a ValueError
        is raised::

            sage: c.cylcoord_to_origami([1,2,3],[1,1])
            Traceback (most recent call last):
            ...
            ValueError: lengths are not compatible with cylinder equations

        TESTS::

            sage: c = CylinderDiagram([((0,),(1,)), ((1,2,3),(0,2,3))])
            sage: c
            Cylinder diagram (0)-(1) (1,2,3)-(0,2,3)
            sage: lengths = [1,1,1,1]
            sage: heights = [1,1]
            sage: c.cylcoord_to_origami(lengths,heights,[0,0])
            (1)(2,3,4)
            (1,2)(3,4)
            sage: c.cylcoord_to_origami(lengths,heights,[0,1])
            (1)(2,3,4)
            (1,2,3)(4)
            sage: c.cylcoord_to_origami(lengths,heights,[0,2])
            (1)(2,3,4)
            (1,2,4)(3)
        """
        from origamis.origami import Origami_dense

        widths = [sum(lengths[i] for i in bot) for bot,_ in self.cylinders()]

        if widths != [sum(lengths[i] for i in top) for _,top in self.cylinders()]:
            raise ValueError, "lengths are not compatible with cylinder equations"

        if twists is None:
            twists = [0] * len(widths)
        elif len(twists) != len(widths):
            raise ValueError, "not enough twists"
        else:
            twists = [(-twists[i])%widths[i] for i in xrange(len(widths))] 
        areas = [heights[i]*widths[i] for i in xrange(self.ncyls())]

        # intialization of partial volumes: the set of squares in cylinder i is range(v[i],v[i+1])
        v = [0]
        for a in areas:
            v.append(v[-1] + a)

        # initialization of bottom squares: sep_i -> bottom position
        sep_bottom_pos = [None] * self.nseps()
        for i,(bot,_) in enumerate(self.cylinders()):
            w = 0
            for j in bot:
                sep_bottom_pos[j] = v[i] + w
                w += lengths[j]

        # build the permutation r
        lx = range(1, v[-1]+1)
        for i in xrange(self.ncyls()):
            for j in xrange(v[i], v[i+1], widths[i]):
                lx[j+widths[i]-1] = j

        # build permutation u with the given twists
        ly = []
        for i,(_,top_seps) in enumerate(self.cylinders()):
            # everything excepted the top
            ly.extend(range(v[i]+widths[i],v[i+1]))

            # the top
            k = top_seps[0]
            top = range(sep_bottom_pos[k],sep_bottom_pos[k]+lengths[k])
            for k in reversed(top_seps[1:]):
                top.extend(range(sep_bottom_pos[k],sep_bottom_pos[k]+lengths[k]))
            ly.extend(top[twists[i]:] + top[:twists[i]])

        # yield the one without twist
        return Origami_dense(tuple(lx), tuple(ly))


    #TODO
#    def chain_complex_dual(self, ring=None):
#        r"""
#        Return a chain complex for the cylinder diagram
#
#        The vertices are in bijection with the cylinder of self
#        The edges are in bijection with separatrices and cylinders
#
#        """
#        from homology import TranslationSurfaceChainComplex
#        from sage.rings.integer import Integer
#
#        if ring is None:
#            from sage.rings.integer_ring import IntegerRing
#            ring = IntegerRing()
#
#        vertices = []   # list of list of vertex = integers from 0 to ncyls
#                        # (in or out, edge)
#        edges = {}      # label -> (start vertex,end)
#        angles = {}     # (in or out,label) -> angle to next edge
#
#        cyls = self.cylinders()
#        t2c = [None]*self.nseps()
#        b2c = [None]*self.nseps()
#
#        for k,(cb,ct) in enumerate(cyls):
#            for i in cb: b2c[i] = k
#            for i in ct: t2c[i] = k
#
#        for k,(cb,ct) in enumerate(cyls):
#            vertex = []
#            e = 'c%d' %k
#            edges[e] = (k,k)
#            angles[(1,e)] = Integer(1)/Integer(2)
#            angles[(-1,e)] = Integer(1)/Integer(2)
#            # the incoming edges from bottom
#            for i in cb:
#                e = 's%d' %i
#                edges[e] = (t2c[i],k)
#                vertex.append((-1,e))
#                angles[(-1,e)] = Integer(0)
#            angles[(-1,e)] = Integer(1)/Integer(2)
#
#            # the central edge (outgoing)
#            vertex.append((1, 'c%d' %k))
#
#            # the outgoing edges from top
#            for i in ct:
#                e = 's%d' %i
#                vertex.append((1,e))
#                angles[(1,e)] = Integer(0)
#            angles[(1,e)] = Integer(1)/Integer(2)
#
#            # the central edge (incoming)
#            vertex.append((-1, 'c%d' %k))
#
#            vertices.append(vertex)
#
#        return TranslationSurfaceChainComplex(ring,vertices,edges,angles)



