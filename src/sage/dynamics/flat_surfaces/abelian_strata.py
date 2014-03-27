r"""
Strata of differentials on Riemann surfaces

The space of Abelian differentials on Riemann surfaces of a given genus is
stratified by degrees of zeros. Each stratum has one, two or three connected
components each of which is associated to an extended Rauzy class. The
:meth:`~sage.combinat.iet.strata.AbelianStratum.components` method
lists connected components of a stratum.

The work for Abelian differentials was done by Maxim Kontsevich and Anton
Zorich in [KonZor03]_ and for quadratic differentials by Erwan Lanneau in
[Lan08]. Zorich gave an algorithm to pass from a connected component of a
stratum to the associated Rauzy class (for both interval exchange
transformations and linear involutions) in [Zor08]_ and is implemented for
Abelian stratum at different level (approximately one for each component):

- for connected stratum :meth:`~ConnectedComponentOfAbelianStratum.permutation_representative`

- for hyperellitic component
  :meth:`~HypConnectedComponentOfAbelianStratum.permutation_representative`

- for non hyperelliptic component, the algorithm is the same as for connected
  component

- for odd component :meth:`~OddConnectedComponentOfAbelianStratum.permutation_representative`

- for even component :meth:`~EvenConnectedComponentOfAbelianStratum.permutation_representative`

The inverse operation (pass from an interval exchange transformation to
the connected component) is partially written in [KonZor03]_ and
simply named here
:meth:`~sage.dynamics.interval_exchanges.template.PermutationIET.components`.

All the code here was first available on Mathematica [ZS]_.

REFERENCES:

.. [KonZor03] M. Kontsevich, A. Zorich "Connected components of the moduli space
   of Abelian differentials with prescripebd singularities" Invent. math. 153,
   631-678 (2003)

.. [Zor08] A. Zorich "Explicit Jenkins-Strebel representatives of all strata of
   Abelian and quadratic differentials", Journal of Modern Dynamics, vol. 2,
   no 1, 139-185 (2008) (http://www.math.psu.edu/jmd)

.. [ZS] Anton Zorich, "Generalized Permutation software"
   (http://perso.univ-rennes1.fr/anton.zorich/Software/software_en.html)

AUTHORS:

- Vincent Delecroix (2009-09-29): initial version


EXAMPLES:

Construction of a stratum from a list of singularity degrees::

    sage: a = AbelianStratum(1,1)
    sage: print a
    H_2(1^2)
    sage: print a.genus()
    2
    sage: print a.dimension()
    5

::

    sage: a = AbelianStratum(4,3,2,1)
    sage: print a
    H_6(4, 3, 2, 1)
    sage: print a.genus()
    6
    sage: print a.dimension()
    15

By convention, the degrees are always written in decreasing order::

    sage: a1 = AbelianStratum(4,3,2,1)
    sage: a1
    H_6(4, 3, 2, 1)
    sage: a2 = AbelianStratum(2,3,1,4)
    sage: a2
    H_6(4, 3, 2, 1)
    sage: a1 == a2
    True

It is possible to lis strata and their connected components::

    sage: AbelianStratum(10).components()
    [H_6(10)^hyp, H_6(10)^odd, H_6(10)^even]

Get a list of strata with constraints on genus or on the number of intervals
of a representative::

    sage: AbelianStrata(genus=3).list()
    [H_3(4), H_3(3, 1), H_3(2^2), H_3(2, 1^2), H_3(1^4)]

Obtains the connected components of a stratum::

    sage: a = AbelianStratum(0)
    sage: print a.components()
    [H_1(0)^hyp]

::

    sage: @cached_function
    ....: def nb_irred_perm(n):
    ....:     if n == 0 or n == 1: return 1
    ....:     return factorial(n) - sum(nb_irred_perm(k) * factorial(n - k) for k in xrange(1,n))
    sage: [nb_irred_perm(i) for i in xrange(10)]
    [1, 1, 1, 3, 13, 71, 461, 3447, 29093, 273343]

::

    sage: A = AbelianStrata(dimension=5, fake_zeros=True)
    sage: N = 0
    sage: for a in A:
    ....:    for cc in a.components():
    ....:       for z in set(a.zeros()):
    ....:           p = cc.permutation_representative(left_degree=z)
    ....:           n = p.rauzy_diagram().cardinality()
    ....:           print "%13s, %d  :  %d"%(cc, z, n)
    ....:           print p
    ....:           N += n
    H_2(2, 0)^hyp, 0  :  11
    0 1 2 3 4
    4 2 1 3 0
    H_2(2, 0)^hyp, 2  :  35
    0 1 2 3 4
    4 1 3 2 0
     H_2(1^2)^hyp, 1  :  15
    0 1 2 3 4
    4 3 2 1 0
     H_1(0^4)^hyp, 0  :  10
    0 1 2 3 4
    4 0 1 2 3
    sage: N
    71
    sage: nb_irred_perm(5)
    71

"""
#*****************************************************************************
#       Copyright (C) 2009 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets

from sage.combinat.partition import Partitions
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.infinity import Infinity

from strata import Stratum, StratumComponent, Strata


class AbelianStratum(Stratum):
    """
    Stratum of Abelian differentials.

    A stratum with a marked outgoing separatrix corresponds to Rauzy diagram
    with left induction, a stratum with marked incoming separatrix correspond
    to Rauzy diagram with right induction.
    If there is no marked separatrix, the associated Rauzy diagram is the
    extended Rauzy diagram (consideration of the
    :meth:`sage.dynamics.interval_exchanges.template.Permutation.symmetric`
    operation of Boissy-Lanneau).

    When you want to specify a marked separatrix, the degree on which it is is
    the first term of your degrees list.

    INPUT:

    - ``marked_separatrix`` - ``None`` (default) or 'in' (for incoming
      separatrix) or 'out' (for outgoing separatrix).

    EXAMPLES:

    Creation of an Abelian stratum and get its connected components::

        sage: a = AbelianStratum(2, 2)
        sage: print a
        H_3(2^2)
        sage: a.components()
        [H_3(2^2)^hyp, H_3(2^2)^odd]

    Get a permutation representative of a connected component::

        sage: a = AbelianStratum(2,2)
        sage: a_hyp, a_odd = a.components()
        sage: print a_hyp.permutation_representative()
        0 1 2 3 4 5 6
        6 5 4 3 2 1 0
        sage: print a_odd.permutation_representative()
        0 1 2 3 4 5 6
        3 2 4 6 5 1 0

    You can specify the alphabet::

        sage: print a_odd.permutation_representative(alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        A B C D E F G
        D C E G F B A
    """
    _name = 'H'
    _latex_name = '\\mathcal{H}'

    def __init__(self, *l, **d):
        """
        TESTS::

            sage: s = AbelianStratum(0)
            sage: s == loads(dumps(s))
            True
            sage: s = AbelianStratum(1,1,1,1)
            sage: s == loads(dumps(s))
            True
        """
        if len(l) == 1:
            try:
                Integer(l[0])
            except StandardError:
                l = l[0]

        if isinstance(l, dict):
            l = sum(([v]*e for v,e in l.iteritems()), [])

        zeros = self._zeros = map(Integer, filter(lambda x: x, l))
        if any(z < 0 for z in self._zeros):
            raise ValueError("the degrees must be non negative")
        self._nb_fake_zeros = sum(1 for x in l if not x)
        self._zeros.sort(reverse=True)

        s = sum(self._zeros)
        if s%2:
            raise ValueError("the sum of the degrees must be even")
        genus = s//2 + 1

        if genus == 1:
            self._cc = (HypASC,)

        elif genus == 2:
            self._cc = (HypASC,)

        elif genus == 3:
            if zeros == [2, 2] or zeros == [4]:
                self._cc = (HypASC, OddASC)
            else:
                self._cc = (ASC,)

        elif len(zeros) == 1:
            # just one zeros [2g-2]
            self._cc = (HypASC, OddASC, EvenASC)

        elif zeros == [genus-1, genus-1]:
            # two similar zeros [g-1, g-1]
            if genus % 2 == 0:
                self._cc = (HypASC, NonHypASC)

            else:
                self._cc = (HypASC, OddASC, EvenASC)

        elif all(x%2 == 0 for x in zeros):
            # even zeroes [2 l_1, 2 l_2, ..., 2 l_n]
            self._cc = (OddASC, EvenASC)

        else:
            # connected
            self._cc = (ASC, )

    def zeros(self, fake_zeros=True):
        r"""
        Returns the list of zeros of self.

        The return list *should not* be modified

        EXAMPLES::

            sage: AbelianStratum([1,2,3]).zeros()
            [3, 2, 1]
            sage: AbelianStratum({2:4}).zeros()
            [2, 2, 2, 2]
        """
        if fake_zeros:
            return self._zeros + [0]*self._nb_fake_zeros
        return self._zeros

    def nb_zeros(self, fake_zeros=True):
        r"""
        Returns the number of zeros of self.

        EXAMPLES::

            sage: AbelianStratum(0).nb_zeros()
            1
            sage: AbelianStratum({2:4,3:2}).nb_zeros()
            6
        """
        if fake_zeros:
            return len(self._zeros) + self._nb_fake_zeros
        return len(self._zeros)

    def nb_fake_zeros(self):
        r"""
        Return the number of fake zeros.

        EXAMPLES::

            sage: AbelianStratum(0).nb_fake_zeros()
            1
            sage: AbelianStratum(1,1,0,0).nb_fake_zeros()
            2

            sage: QuadraticStratum(0,4,2,2).nb_fake_zeros()
            1
         """
        return self._nb_fake_zeros

    def genus(self):
        r"""
        Return the genus of the stratum.

        OUTPUT:

        integer -- the genus

        EXAMPLES::

            sage: AbelianStratum(0).genus()
            1
            sage: AbelianStratum(1,1).genus()
            2
            sage: AbelianStratum(3,2,1).genus()
            4
        """
        return Integer(sum(self._zeros)//2+1)

    def dimension(self):
        r"""
        Return the complex dimension of this stratum.

        The dimension is `2g-2+s+1` where `g` is the genus of surfaces in the
        stratum, `s` the number of singularities. The complex dimension of a
        stratum is also the number of intervals of any interval exchange
        transformations associated to the strata.

        EXAMPLES::

            sage: AbelianStratum(0).dimension()
            2
            sage: AbelianStratum(0,0).dimension()
            3
            sage: AbelianStratum(2).dimension()
            4
            sage: AbelianStratum(1,1).dimension()
            5

        ::

            sage: a = AbelianStratum(4,3,2,1,0)
            sage: p = a.permutation_representative()
            sage: len(p) == a.dimension()
            True
        """
        return 2 * self.genus() + self.nb_zeros() - 1

    #
    # Connected component
    #

    def has_odd_component(self):
        return all(z%2 == 1 for z in self.zeros() and self.genus() != 2)

    def has_even_component(self):
        return all(z%2 == 0 for z in self.zeros() and self.genus() >= 4)

    def has_hyperelliptic_component(self):
        z = self.zeros()
        return (len(z) == 1) or (len(z) == 2 and z[0] == z[1])

    def has_non_hyperelliptic_component(self):
        z = self.zeros()
        return ((len(z) == 2) and (z[0] == z[1]) and (z[0]%2 == 1))

    def unique_component(self):
        r"""
        If the stratum is connected, return the unique component.

        EXAMPLES::

            sage: a = AbelianStratum(1,1); a
            H_2(1^2)
            sage: a.unique_component()
            H_2(1^2)^hyp

            sage: a = AbelianStratum(3,2,1); a
            H_4(3, 2, 1)
            sage: a.unique_component()
            H_4(3, 2, 1)^c
        """
        if self.is_connected(): return self._cc[0](self)
        raise ValueError, "The stratum is not connected"

    def odd_component(self):
        r"""
        Return the odd component of self (if any).

        EXAMPLES::

            sage: a = AbelianStratum([2,2]); a
            H_3(2^2)
            sage: a.odd_component()
            H_3(2^2)^odd
        """
        if OddASC in self._cc: return OddASC(self)
        raise ValueError, "No odd spin component in this stratum"

    def even_component(self):
        r"""
        Return the even component of self (if any)

        EXAMPLES::

            sage: a = AbelianStratum({2:4}); a
            H_5(2^4)
            sage: a.even_component()
            H_5(2^4)^even
        """
        if EvenASC in self._cc: return EvenASC(self)
        raise ValueError, "No even spin component in this stratum"

    def hyperelliptic_component(self):
        r"""
        Return the hyperelliptic component of self (if any)

        EXAMPLES::

            sage: a = AbelianStratum(10); a
            H_6(10)
            sage: a.hyperelliptic_component()
            H_6(10)^hyp
        """
        if HypASC in self._cc: return HypASC(self)
        raise ValueError, "No hyperelliptic component in this stratum"

    def non_hyperelliptic_component(self):
        r"""
        Return the non hyperelliptic component of self (if any)

        EXAMPLES::

            sage: a = AbelianStratum(3,3); a
            H_4(3^2)
            sage: a.non_hyperelliptic_component()
            H_4(3^2)^nonhyp
        """
        if NonHypASC in self._cc: return NonHypASC(self)
        raise ValueError, "No non hyperelliptic component in this stratum"


    #
    # Quadratic cover
    #

    def orientation_quotients(self,fake_zeros=False):
        r"""
        Return the list of quadratic strata such that their orientation cover
        are contained in this stratum.

        If ``fake_zeros`` (default: False) is True we do care about poles which
        becomes a marked zero.

        EXAMPLES:

        The stratum H(2g-2) has one conic singularities of angle `2(2g-1)pi`. The
        only way a surface in H(2g-2) covers a quadratic differential is that
        the quadratic differential has as unique zeros a conical singularity of
        angle `(2g-1) \pi`. The number of poles may vary and give a collection
        of possibilities::

            sage: AbelianStratum(2).orientation_quotients()
            [Q_0(1, -1^5)]
            sage: AbelianStratum(4).orientation_quotients()
            [Q_1(3, -1^3), Q_0(3, -1^7)]
            sage: AbelianStratum(6).orientation_quotients()
            [Q_2(5, -1), Q_1(5, -1^5), Q_0(5, -1^9)]

        A stratum with two zeros may or may not have orientation quotients::

            sage: AbelianStratum(1,1).orientation_quotients()
            [Q_1(2, -1^2), Q_0(2, -1^6)]
            sage: AbelianStratum(2,2).orientation_quotients()
            [Q_1(1^2, -1^2), Q_0(1^2, -1^6), Q_1(4, -1^4), Q_0(4, -1^8)]
            sage: AbelianStratum(3,1).orientation_quotients()
            []

        To impose that covering of poles are fake zeros, switch option
        ``fake_zeros`` to ``True``::

            sage: AbelianStratum(2,2,0,0).orientation_quotients(fake_zeros=True)
            [Q_1(1^2, -1^2)]
        """
        e = {}
        for i in self.zeros(fake_zeros=False):
            if i not in e: e[i] = 0
            e[i] += 1

        # the odd degrees (corresponding to angles 2((2m+1)+1) times pi should
        # be non ramified and hence come by pair.
        if any(e[i]%2 for i in e if i%2):
            return []

        pairings = []
        for d,m in e.iteritems():
            if d%2: # if the degree is odd it is necessarily non ramified
                pairings.append([(d,m//2)])
            else: # if the degree is even ramified and non ramified are possible
                pairings.append([(d,k) for k in range(m//2+1)])

        import itertools
        from quadratic_strata import QuadraticStratum
        res = []

        for p in itertools.product(*pairings):
            ee = dict((d-1,0) for d in e)
            ee.update((2*d,0) for d in e)
            for d,m in p:
                ee[d-1] += e[d]-2*m
                ee[2*d] += m

            degrees = []
            for d in ee: degrees.extend([d]*ee[d])

            s = sum(degrees)
            for nb_poles in xrange(s%4,s+5,4):
                q = QuadraticStratum(degrees + [-1]*nb_poles)
                if not q.is_empty() and (not fake_zeros or q.nb_poles() <= self.nb_fake_zeros()):
                    res.append(q)

        return res

class AbelianStratumComponent(StratumComponent):
    r"""
    Connected component of Abelian stratum.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'c'

    def spin(self):
        r"""
        Return ``None`` since surfaces in this component have no spin.

        EXAMPLES::

            sage: c = AbelianStratum([1,1,1,1]).unique_component(); c
            H_3(1^4)^c
            sage: c.spin() is None
            True
        """
        return None
    def permutation_representative(self, left_degree=None, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitely interval exchange
        transformations for each stratum in [Zor08]_.

        INPUT:

        - ``reduced`` - boolean (default: ``True``): whether you
          obtain a reduced or labelled permutation

        - ``alphabet`` - an alphabet or ``None``: whether you want to
          specify an alphabet for your permutation

        - ``left_degree`` - the degree of the singularity on the left of the
          interval.

        OUTPUT:

        permutation -- a permutation which lives in this component

        EXAMPLES::

            sage: c = AbelianStratum(1,1,1,1).unique_component()
            sage: p = c.permutation_representative(alphabet="abcdefghi")
            sage: p
            a b c d e f g h i
            i d c b e h g f a
            sage: print p.stratum_component()
            H_3(1^4)^c

            sage: cc = AbelianStratum(3,2,1,0).unique_component()
            sage: p = cc.permutation_representative(left_degree=3); p
            0 1 2 3 4 5 6 7 8 9 10
            10 1 3 2 4 7 6 5 9 8 0
            sage: p.stratum_component()
            H_4(3, 2, 1, 0)^c
            sage: p.marking().left()
            4
            sage: p.rauzy_diagram()  # long time
            Rauzy diagram with 1060774 permutations

            sage: p = cc.permutation_representative(left_degree=2); p
            0 1 2 3 4 5 6 7 8 9 10
            10 1 4 3 2 6 5 7 9 8 0
            sage: p.stratum_component()
            H_4(3, 2, 1, 0)^c
            sage: p.marking().left()
            3
            sage: p.rauzy_diagram()  # long time
            Rauzy diagram with 792066 permutations

            sage: p = cc.permutation_representative(left_degree=1); p
            0 1 2 3 4 5 6 7 8 9 10
            10 1 3 2 6 5 4 9 8 7 0
            sage: p.stratum_component()
            H_4(3, 2, 1, 0)^c
            sage: p.marking().left()
            2
            sage: p.rauzy_diagram()  # long time
            Rauzy diagram with 538494 permutations

            sage: p = cc.permutation_representative(left_degree=0); p
            0 1 2 3 4 5 6 7 8 9 10
            10 2 1 5 4 3 8 7 6 9 0
            sage: p.stratum_component()
            H_4(3, 2, 1, 0)^c
            sage: p.marking().left()
            1
            sage: p.rauzy_diagram()  # long time
            Rauzy diagram with 246914 permutations
        """
        z = self.stratum().zeros(fake_zeros=False)
        n = self.stratum().nb_fake_zeros()

        if left_degree is not None:
            if not isinstance(left_degree, (int,Integer)):
                raise ValueError, "left_degree (=%d) should be one of the degree"%left_degree
            if left_degree == 0:
                if n == 0:
                    raise ValueError, "left_degree (=%d) should be one of the degree"%left_degree
            elif left_degree not in z:
                raise ValueError, "left_degree (=%d) should be one of the degree"%left_degree
            else:
                z.remove(left_degree)
                z.insert(len(z),left_degree)


        k = 1
        l0 = [0]
        l1 = [1]
        for d in z:
            for i in xrange(d):
                l0.append(2*k)
                l1.append(2*(k+1) if k%2 else 2*(k-1))
                k += 1
            l0.append(2*k-1)
            l1.append(2*k-1)

        l0.pop()
        l1.pop()
        l0.append(1)
        l1.append(0)

        if n != 0:
            g = self.stratum().genus()
            interval = range(5*g, 5*g+n)

            if left_degree == 0:
                l0[-1:-1] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1[1:1] = interval


        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            p = ReducedPermutationIET([l0, l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

    def rauzy_diagram(self, reduced=True):
        r"""
        Returns the Rauzy diagram associated to this connected component.

        OUTPUT:

        rauzy diagram -- the Rauzy diagram associated to this stratum

        EXAMPLES::

            sage: c = AbelianStratum(0).components()[0]
            sage: r = c.rauzy_diagram()
        """
        return self.permutation_representative(reduced=reduced).rauzy_diagram()

ASC = AbelianStratumComponent

class HypAbelianStratumComponent(ASC):
    """
    Hyperelliptic component of Abelian stratum.
    """
    _name = 'hyp'

    def spin(self):
        r"""
        Return the spin parity of hyperelliptic stratum.

        EXAMPLES:

        For the strata `H(2g-2)`::

            sage: c = AbelianStratum(0).hyperelliptic_component()
            sage: c.spin()
            1
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            1

            sage: c = AbelianStratum(2).hyperelliptic_component()
            sage: c.spin()
            1
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            1

            sage: c = AbelianStratum(4).hyperelliptic_component()
            sage: c.spin()
            0
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            0

        For the strata `H(g-1,g-1)`::

            sage: c = AbelianStratum(2,2).hyperelliptic_component()
            sage: c.spin()
            0
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            0

            sage: c = AbelianStratum(4,4).hyperelliptic_component()
            sage: c.spin()
            1
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            1
        """
        z = self.stratum().zeros(fake_zeros=False)
        if len(z) == 0:
            return Integer(1)
        elif len(z) == 1:
            return Integer(((self.stratum().genus()+1)//2) % 2)
        elif len(z) == 2:
            if z[0] % 2:
                return None
            return Integer(((self.stratum().genus()+1)//2) %2)

    def permutation_representative(self, left_degree=None, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitely interval exchange
        transformations for each stratum in [Zor08]_.

        INPUT:

        - ``reduced`` - boolean (defaut: ``True``): whether you obtain
          a reduced or labelled permutation

        - ``alphabet`` - alphabet or ``None`` (defaut: ``None``):
          whether you want to specify an alphabet for your
          representative

        EXAMPLES::

            sage: c = AbelianStratum(0).hyperelliptic_component()
            sage: p = c.permutation_representative()
            sage: p
            0 1
            1 0
            sage: p.stratum_component()
            H_1(0)^hyp

            sage: c = AbelianStratum(0,0).hyperelliptic_component()
            sage: p = c.permutation_representative(alphabet="abc")
            sage: p
            a b c
            c b a
            sage: p.stratum_component()
            H_1(0^2)^hyp

            sage: c = AbelianStratum(2,2).hyperelliptic_component()
            sage: p = c.permutation_representative(alphabet="ABCDEFGHIJKL")
            sage: p
            A B C D E F G
            G F E D C B A
            sage: c = AbelianStratum(1,1,0).hyperelliptic_component()
            sage: p = c.permutation_representative(left_degree=1); p
            0 1 2 3 4 5
            5 1 4 3 2 0
            sage: p.marking().left()
            2
            sage: p.rauzy_diagram()
            Rauzy diagram with 90 permutations

            sage: p = c.permutation_representative(left_degree=0); p
            0 1 2 3 4 5
            5 3 2 1 4 0
            sage: p.marking().left()
            1
            sage: p.rauzy_diagram()
            Rauzy diagram with 20 permutations
        """
        g = self._stratum.genus()
        n = self._stratum.nb_fake_zeros()
        m = len(self._stratum.zeros(fake_zeros=False))

        if left_degree is not None:
            if not isinstance(left_degree, (int,Integer)) or left_degree not in self.stratum().zeros():
                raise ValueError("left_degree (=%d) should be one of the degree"%left_degree)

        if m == 0:  # on the torus
            if n == 1:
                l0 = [0, 1]
                l1 = [1, 0]
            elif n == 2:
                l0 = [0, 1, 2]
                l1 = [2, 1, 0]
            else:
                l0 = range(1, n+2)
                l1 = [n+1] + range(1, n+1)

        elif m == 1:  # H(2g-2,0^n) or H(0,2g-2,0^(n-1))
            l0 = range(1, 2*g+1)
            l1 = range(2*g, 0, -1)
            interval = range(2*g+1, 2*g+n+1)

            if left_degree == 0:
                l0[-1:-1] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1[1:1] = interval

        else:  # H(g-1,g-1,0^n) or H(0,g-1,g-1,0^(n-1))
            l0 = range(1, 2*g+2)
            l1 = range(2*g+1, 0, -1)
            interval = range(2*g+2, 2*g+n+2)

            if left_degree == 0:
                l0[-1:-1] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1[1:1] = interval

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            p = ReducedPermutationIET([l0, l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1])
        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

HypASC = HypAbelianStratumComponent


class NonHypAbelianStratumComponent(ASC):
    """
    Non hyperelliptic component of Abelian stratum.
    """
    _name = 'nonhyp'

NonHypASC = NonHypAbelianStratumComponent

class EvenAbelianStratumComponent(ASC):
    """
    Connected component of Abelian stratum with even spin structure.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'even'

    def spin(self):
        r"""
        Return ``0``.

        EXAMPLES::

            sage: c = AbelianStratum(4,2).even_component(); c
            H_4(4, 2)^even
            sage: c.spin()
            0
        """
        return Integer(0)

    def permutation_representative(self, left_degree=None, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitely interval exchange
        transformations for each stratum in [Zor08]_.

        INPUT:

        - ``reduced`` - boolean (defaut: True): whether you obtain a reduced or
          labelled permutation

        - ``left_degree`` - integer (optional) - a specified degree of zero at
          the left of the interval.

        - ``alphabet`` - alphabet or None (defaut: None): whether you want to
          specify an alphabet for your representative

        - ``relabel`` - boolean (default: True) - if False uses Zorich's natural
          numbering otherwise uses 0,1,...

        EXAMPLES::

            sage: c = AbelianStratum(6).even_component()
            sage: c
            H_4(6)^even
            sage: p = c.permutation_representative(alphabet=range(8))
            sage: p
            0 1 2 3 4 5 6 7
            5 4 3 2 7 6 1 0
            sage: p.stratum_component()
            H_4(6)^even

        ::

            sage: c = AbelianStratum(4,4).even_component()
            sage: c
            H_5(4^2)^even
            sage: p = c.permutation_representative(alphabet=range(11))
            sage: p
            0 1 2 3 4 5 6 7 8 9 10
            5 4 3 2 6 8 7 10 9 1 0
            sage: p.stratum_component()
            H_5(4^2)^even

        Different markings lead to different Rauzy diagrams::

            sage: c = AbelianStratum(4,2,0).even_component()
            sage: p = c.permutation_representative(left_degree=4); p
            0 1 2 3 4 5 6 7 8 9
            6 5 4 3 7 9 8 2 0 1
            sage: p.stratum_component()
            H_4(4, 2, 0)^even
            sage: p.marking().left()
            5
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 66140 permutations

            sage: p = c.permutation_representative(left_degree=2); p
            0 1 2 3 4 5 6 7 8 9
            7 6 5 4 3 9 8 2 0 1
            sage: p.stratum_component()
            H_4(4, 2, 0)^even
            sage: p.marking().left()
            3
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 39540 permutations

            sage: p = c.permutation_representative(left_degree=0); p
            0 1 2 3 4 5 6 7 8 9
            7 5 4 3 2 9 8 1 6 0
            sage: p.stratum_component()
            H_4(4, 2, 0)^even
            sage: p.marking().left()
            1
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 11792 permutations
        """
        z = self._stratum.zeros(fake_zeros=False)
        n = self._stratum.nb_fake_zeros()
        g = self._stratum.genus()

        if left_degree is not None:
            if not isinstance(left_degree, (int,Integer)):
                raise ValueError, "left_degree (=%d) should be one of the degree"%left_degree
            if left_degree == 0:
                if n == 0:
                    raise ValueError, "left_degree (=%d) should be one of the degree"%left_degree
            elif left_degree not in z:
                raise ValueError, "left_degree (=%d) should be one of the degree"%left_degree
            else:
                z.remove(left_degree)
                z.insert(0,left_degree)

        l0 = range(3*g-2)
        l1 = [6, 5, 4, 3, 2, 7, 9, 8]
        for k in range(10, 3*g-4, 3):
            l1 += [k, k+2, k+1]
        l1 += [1, 0]

        k = 4
        for d in z:
            for i in range(d/2-1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 3
            k += 3

        # if there are marked points we transform 0 in [3g-2, 3g-3, ...]
        if n != 0:
            interval = range(3*g-2, 3*g - 2 + n)

            if left_degree == 0:
                k = l0.index(6)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            p = ReducedPermutationIET([l0, l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

EvenASC = EvenAbelianStratumComponent

class OddAbelianStratumComponent(ASC):
    r"""
    Connected component of an Abelian stratum with odd spin parity.
    """
    _name = 'odd'
    def spin(self):
        r"""
        Returns 1 which is, by definition, the spin parity of this stratum component.

        EXAMPLES::

            sage: c = AbelianStratum(4).odd_component(); c
            H_3(4)^odd
            sage: c.spin()
            1
         """
        return 1

    def permutation_representative(self, left_degree=None, reduced=True, alphabet=None, relabel=True):
        """
        Returns the Zorich representative of this connected component.

        A. Zorich constructs explicitely interval exchange
        transformations for each stratum in [Zor08]_.

        EXAMPLES::

            sage: a = AbelianStratum(6).odd_component()
            sage: p = a.permutation_representative()
            sage: p
            0 1 2 3 4 5 6 7
            3 2 5 4 7 6 1 0
            sage: p.stratum_component()
            H_4(6)^odd

        ::

            sage: a = AbelianStratum(4,4).odd_component()
            sage: p = a.permutation_representative()
            sage: p
            0 1 2 3 4 5 6 7 8 9 10
            3 2 5 4 6 8 7 10 9 1 0
            sage: p.stratum_component()
            H_5(4^2)^odd

        Different markings lead to different Rauzy diagrams::

            sage: c = AbelianStratum(4,2,0).odd_component()
            sage: p = c.permutation_representative(left_degree=4); p
            0 1 2 3 4 5 6 7 8 9
            4 3 6 5 7 9 8 2 0 1
            sage: p.stratum_component()
            H_4(4, 2, 0)^odd
            sage: p.marking().left()
            5
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 147090 permutations

            sage: p = c.permutation_representative(left_degree=2); p
            0 1 2 3 4 5 6 7 8 9
            4 3 5 7 6 9 8 2 0 1
            sage: p.stratum_component()
            H_4(4, 2, 0)^odd
            sage: p.marking().left()
            3
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 87970 permutations

            sage: p = c.permutation_representative(left_degree=0); p
            0 1 2 3 4 5 6 7 8 9
            4 2 5 7 6 9 8 1 3 0
            sage: p.stratum_component()
            H_4(4, 2, 0)^odd
            sage: p.marking().left()
            1
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 27754 permutations
        """
        zeros = self.stratum().zeros(fake_zeros=False)
        n = self._stratum.nb_fake_zeros()
        g = self._stratum.genus()

        if left_degree is not None:
            if not isinstance(left_degree, (int,Integer)):
                raise ValueError, "left_degree (=%d) should be one of the degree"%left_degree
            if left_degree == 0:
                if n == 0:
                    raise ValueError, "left_degree (=%d) should be one of the degree"%left_degree
            elif left_degree not in zeros:
                raise ValueError, "left_degree (=%d) should be one of the degree"%left_degree
            else:
                zeros.remove(left_degree)
                zeros.insert(0,left_degree)

        z = [x//2 for x in zeros]

        l0 = range(3*g-2)
        l1 = [3, 2]
        for k in range(4, 3*g-4, 3):
            l1 += [k, k+2, k+1]
        l1 += [1, 0]

        k = 4
        for d in z:
            for i in range(d-1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 3
            k += 3

        # marked points
        if n != 0:
            interval = range(3*g-2, 3*g-2+n)

            if left_degree == 0:
                k = l0.index(3)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            p = ReducedPermutationIET([l0, l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

OddASC = OddAbelianStratumComponent





#
# iterators for Abelian strata with constraints on genus and dimension
#

def AbelianStrata(genus=None, dimension=None, fake_zeros=None):
    r"""
    Abelian strata.

    INPUT:

    - ``genus`` - a non negative integer or None

    - ``dimension`` - a non negative integer or None

    - ``fake_zeros`` - boolean

    EXAMPLES:

    Abelian strata with a given genus::

        sage: for s in AbelianStrata(genus=1): print s
        H_1(0)

    ::

        sage: for s in AbelianStrata(genus=2): print s
        H_2(2)
        H_2(1^2)

    ::

        sage: for s in AbelianStrata(genus=3): print s
        H_3(4)
        H_3(3, 1)
        H_3(2^2)
        H_3(2, 1^2)
        H_3(1^4)

    ::

        sage: for s in AbelianStrata(genus=4): print s
        H_4(6)
        H_4(5, 1)
        H_4(4, 2)
        H_4(4, 1^2)
        H_4(3^2)
        H_4(3, 2, 1)
        H_4(3, 1^3)
        H_4(2^3)
        H_4(2^2, 1^2)
        H_4(2, 1^4)
        H_4(1^6)

    Get outside of the tests.
    Abelian strata with a given number of intervals

    sage for s in AbelianStrata(dimension=2): print s
    H^out([0])

    sage for s in AbelianStrata(dimension=3): print s
    H^out([0], 0)

    sage for s in AbelianStrata(dimension=4): print s
    H^out([2])
    H^out([0], 0, 0)


    Get outisde of tests
    sage  for s in AbelianStrata(dimension=5): print s
    H^out(2, [0])
    H^out([2], 0)
    H^out([1], 1)
    H^out([0], 0, 0, 0)
    """
    fake_zeros = bool(fake_zeros)
    if genus is None:
        if dimension is None:
            return AbelianStrata_all()
        dimension = Integer(dimension)
        return AbelianStrata_d(dimension,fake_zeros)
    genus = Integer(genus)

    if dimension is None:
        return AbelianStrata_g(genus)
    dimension = Integer(dimension)
    return AbelianStrata_gd(genus,dimension,fake_zeros)

class AbelianStrata_class(Strata):
    r"""
    Generic class for abelian strata.
    """
    pass

class AbelianStrata_g(AbelianStrata_class):
    r"""
    Stratas of genus g surfaces without fake zeros.

    INPUT:

    - ``genus`` - a non negative integer

    EXAMPLES::

        sage: AbelianStrata(genus=2).list()
        [H_2(2), H_2(1^2)]
        sage: AbelianStrata(genus=3).list()
        [H_3(4), H_3(3, 1), H_3(2^2), H_3(2, 1^2), H_3(1^4)]
        sage: AbelianStrata(genus=4).random_element() #random
        H_4(4, 2)
    """
    def __init__(self,genus):
        r"""
        TESTS::

            sage: s = AbelianStrata(genus=3)
            sage: TestSuite(s).run()
            sage: loads(dumps(s)) == s
            True
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._genus = genus

    def __eq__(self, other):
        return isinstance(other, AbelianStrata_g) and self._genus == other._genus

    def __ne__(self, other):
        return not self.__eq__(other)

    def __reduce__(self):
        return (AbelianStrata_g,(self._genus,))

    def __contains__(self, c):
        r"""
        Containance test

        TESTS::

            sage: a = AbelianStrata(genus=3)
            sage: all(s in a for s in a)
            True

            sage: a = AbelianStrata(genus=3,fake_zeros=False)
            sage: all(s in a for s in a)
            True
        """
        if not isinstance(c, AbelianStratum):
            return False

        return c.genus() == self._genus

    def _repr_(self):
        r"""
        TESTS::

            sage: repr(AbelianStrata(genus=3))   #indirect doctest
            'Abelian strata of genus 3 surfaces'
        """
        return "Abelian strata of genus %d surfaces"%self._genus

    def cardinality(self):
        r"""
        Return the number of abelian strata with a given genus.

        EXAMPLES::

            sage: AbelianStrata(genus=1).cardinality()
            1
            sage: AbelianStrata(genus=2).cardinality()
            2
            sage: AbelianStrata(genus=3).cardinality()
            5
            sage: AbelianStrata(genus=4).cardinality()
            11
        """
        if self._genus == 0:
            return Integer(0)
        if self._genus == 1:
            return Integer(1)
        return Partitions(2*self._genus-2).cardinality()

    def __iter__(self):
        r"""
        TESTS::

            sage: list(AbelianStrata(genus=1))
            [H_1(0)]
        """
        if self._genus == 0:
            pass
        elif self._genus == 1:
            yield AbelianStratum([0])
        else:
            for p in Partitions(2*self._genus-2):
                yield AbelianStratum(p._list)

    def random_element(self):
        r"""
        Return a random stratum.
        """
        if self._genus == 0:
            raise ValueError, "No stratum with that genus"
        if self._genus == 1:
            return AbelianStratum([0])
        return AbelianStratum(Partitions(2*self._genus - 2).random_element())

    def first(self):
        r"""
        Return the first element of this list of strata.

        EXAMPLES::

            sage: AbelianStrata(genus=3).first()
            H_3(4)
            sage: AbelianStrata(genus=4).first()
            H_4(6)
        """
        return AbelianStratum([2*self._genus-2])

    an_element_ = first

    def last(self):
        r"""
        Return the last element of this list of strata.

        EXAMPLES::

            sage: AbelianStrata(genus=4).last()
            H_4(1^6)
            sage: AbelianStrata(genus=5).last()
            H_5(1^8)
        """
        return AbelianStratum({1:2*self._genus-2})

class AbelianStrata_d(AbelianStrata_class):
    r"""
    Strata with prescribed dimension.

    INPUT:

    - ``dimension`` - an integer greater than 1

    - ``fake_zeros`` - boolean (dafault: False) - allows or not fake zeros

    EXAMPLES::

        sage: for a in AbelianStrata(dimension=5,fake_zeros=True):
        ...      print a
        ...      print a.permutation_representative()
        H_2(2, 0)
        0 1 2 3 4
        4 1 3 2 0
        H_2(1^2)
        0 1 2 3 4
        4 3 2 1 0
        H_1(0^4)
        0 1 2 3 4
        4 0 1 2 3
    """
    def __init__(self,dimension,fake_zeros):
        r"""
        TESTS::

            sage: s = AbelianStrata(dimension=10,fake_zeros=True)
            sage: TestSuite(s).run()
            sage: loads(dumps(s)) == s
            True

            sage: s = AbelianStrata(dimension=10,fake_zeros=False)
            sage: TestSuite(s).run()
            sage: loads(dumps(s)) == s
            True
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._dimension = dimension
        self._fake_zeros = fake_zeros

    def __eq__(self, other):
        return (isinstance(other, AbelianStrata_d) and
                (self._dimension == other._dimension) and
                (self._fake_zeros == other._fake_zeros))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, c):
        r"""
        Containance test

        TESTS::

            sage: a = AbelianStrata(dimension=7,fake_zeros=True)
            sage: all(s in a for s in a)
            True

            sage: a = AbelianStrata(dimension=7,fake_zeros=False)
            sage: all(s in a for s in a)
            True
        """
        if not isinstance(c, AbelianStratum):
            return False

        if c.dimension() != self._dimension:
            return False

        if (self._fake_zeros is False) and c.nb_fake_zeros():
            return False

        return True

    def __reduce__(self):
        return (AbelianStrata_d, (self._dimension,self._fake_zeros))

    def _repr_(self):
        r"""
        TESTS::

            sage: repr(AbelianStrata(dimension=2))  #indirect doctest
            'Abelian strata of dimension 2'
        """
        return "Abelian strata of dimension %d" %(self._dimension)

    def first(self):
        r"""
        Returns the first stratum.

        EXAMPLES::

            sage: AbelianStrata(dimension=2).first()
            H_1(0)
            sage: AbelianStrata(dimension=3).first()
            H_1(0^2)
            sage: AbelianStrata(dimension=4).first()
            H_2(2)
        """
        n = self._dimension
        if n%2:
            return AbelianStratum([(n-3)//2,(n-3)//2])
        return AbelianStratum([n-2])

    an_element = first

    def last(self):
        r"""
        Return the last stratum.

        EXAMPLES::

            sage: AbelianStrata(dimension=9,fake_zeros=True).last()
            H_1(0^8)
            sage: AbelianStrata(dimension=9,fake_zeros=False).last()
            H_3(1^4)

            sage: AbelianStrata(dimension=10,fake_zeros=True).last()
            H_1(0^9)
            sage: AbelianStrata(dimension=10,fake_zeros=False).last()
            H_4(2^3)
        """
        n = self._dimension
        if self._fake_zeros:
            return AbelianStratum({0:n-1})
        else:
            if n == 4:
                return AbelianStratum([2])
            if n == 5:
                return AbelianStratum([1,1])
            elif n == 6:
                return AbelianStratum([4])
            else:
                nn = (n-2)%4
                return AbelianStratum({2:3-nn,1:2*((n-10)//4)+2*nn})

    def __iter__(self):
        r"""
        Iterator.

        TESTS::

            sage: for a in AbelianStrata(dimension=4,fake_zeros=True): print a
            H_2(2)
            H_1(0^3)
        """
        n = self._dimension
        if n < 2:
            pass
        elif self._fake_zeros:
            for s in xrange(1+n%2, n, 2):
                for p in Partitions(n-1, length=s):
                    yield AbelianStratum([k-1 for k in p])
        else:
            if n == 2:
                yield AbelianStratum([0])
            else:
                for s in xrange(1+n%2, n, 2):
                    for p in Partitions(n-1,length=s,min_part=2):
                        yield AbelianStratum([k-1 for k in p])

    def cardinality(self):
        r"""
        Return the number of Abelian strata with given dimension.

        EXAMPLES::

            sage: AbelianStrata(dimension=5,fake_zeros=True).cardinality()
            3
            sage: AbelianStrata(dimension=5,fake_zeros=False).cardinality()
            1

            sage: AbelianStrata(dimension=6,fake_zeros=True).cardinality()
            4
            sage: AbelianStrata(dimension=6,fake_zeros=False).cardinality()
            1

            sage: AbelianStrata(dimension=7,fake_zeros=True).cardinality()
            6
            sage: AbelianStrata(dimension=7,fake_zeros=False).cardinality()
            2

            sage: AbelianStrata(dimension=12,fake_zeros=True).cardinality()
            29
            sage: AbelianStrata(dimension=12,fake_zeros=False).cardinality()
            7

        TESTS::

            sage: for d in xrange(1,15):
            ...     A = AbelianStrata(dimension=d,fake_zeros=True)
            ...     assert len(A.list()) == A.cardinality()
            ...     A = AbelianStrata(dimension=d,fake_zeros=False)
            ...     assert len(A.list()) == A.cardinality()
        """
        n = self._dimension
        if n < 2:
            return Integer(0)

        if self._fake_zeros:
            return sum(Partitions(n-1,length=s).cardinality() for s in xrange(1+n%2,n,2))
        if n == 2:
            return Integer(1)
        return sum(Partitions(n-1,length=s,min_part=2).cardinality() for s in xrange(1+n%2,n,2))

#TODO: first, last, cardinality
class AbelianStrata_gd(AbelianStrata_class):
    r"""
    Abelian strata with presrcribed genus and dimension.

    INPUT:

    - ``genus`` - an integer - the genus of the surfaces

    - ``dimension`` - an integer - the dimension of strata

    - ``fake_zeros`` - boolean - allows or not fake zeros

    """
    def __init__(self,genus,dimension,fake_zeros):
        r"""
        TESTS::

            sage: s = AbelianStrata(genus=4,dimension=10)
            sage: loads(dumps(s)) == s
            True
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._genus = genus
        self._dimension = dimension
        self._fake_zeros = fake_zeros

    def __eq__(self, other):
        return (isinstance(other, AbelianStrata_class) and
                self._genus == other._genus and
                self._dimension == other._dimension and
                self._fake_zeros == other._fake_zeros)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, c):
        r"""
        Containance test

        TESTS::

            sage: a = AbelianStrata(dimension=7,fake_zeros=True)
            sage: all(s in a for s in a)
            True

            sage: a = AbelianStrata(dimension=7,fake_zeros=False)
            sage: all(s in a for s in a)
            True
        """
        if not isinstance(c, AbelianStratum):
            return False

        if c.dimension() != self._dimension:
            return False
        if c.genus() != self._genus:
            return False

        if not fake_zeros and c.nb_fake_zeros():
            return False

        return True

    def __reduce__(self):
        return (AbelianStrata_gd, (self._genus,self._dimension,self._fake_zeros))

    def _repr_(self):
        r"""
        TESTS::

            sage: a = AbelianStrata(genus=2,dimension=4)
            sage: repr(a)   #indirect doctest
            'Abelian strata of genus 2 surfaces and dimension 4'
        """
        return "Abelian strata of genus %d surfaces and dimension %d"  %(self._genus, self._dimension)

    def __iter__(self):
        r"""
        TESTS::

            sage: list(AbelianStrata(genus=2,dimension=4))
            [H_2(2)]
        """
        min_part = 1
        if self._fake_zeros is False:
            min_part = 2

        if self._genus == 0 or self._dimension <= 1:
            pass
        elif self._genus == 1:
            if self._fake_zeros is False:
                if self._dimension == 2:
                    yield AbelianStratum([0])
            elif self._dimension >= 2:
                    yield AbelianStratum([0]*(self._dimension-1))
        else:
            s = self._dimension - 2*self._genus + 1
            for p in Partitions(2*self._genus - 2 + s, length=s, min_part=min_part):
                yield AbelianStratum([k-1 for k in p])

class AbelianStrata_all(AbelianStrata_class):
    r"""
    Abelian strata without fake zeros.

    EXAMPLES::

        sage: A = AbelianStrata(); A
        Abelian strata
        sage: a = iter(A)
        sage: for _ in xrange(10): print a.next()
        H_1(0)
        H_2(2)
        H_2(1^2)
        H_3(4)
        H_3(3, 1)
        H_3(2^2)
        H_3(2, 1^2)
        H_3(1^4)
        H_4(6)
        H_4(5, 1)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: s = AbelianStrata()
            sage: TestSuite(s).run()
            sage: loads(dumps(s)) == s
            True
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())

    def __eq__(self, other):
        return isinstance(other,AbelianStrata_all)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self,c):
        r"""
        Contain test.
        """
        if not isinstance(c, AbelianStratum):
            return False

        if c.genus() == 1:
            return c.nb_fake_zeros() == 1

        return c.nb_fake_zeros() == 0

    def _repr_(self):
        r"""
        TESTS::

            sage: repr(AbelianStrata())   #indirect doctest
            'Abelian strata'
        """
        return "Abelian strata"

    def cardinality(self):
        r"""
        Return infinity.

        EXAMPLES::

            sage: AbelianStrata().cardinality()
            +Infinity
        """
        return Infinity

    def first(self):
        return AbelianStratum([0])

    an_element = first

    def __iter__(self):
        """
        Iterator.

        TESTS::

            sage: iter(AbelianStrata()).next()
            H_1(0)
        """
        from itertools import count
        for g in count(1):
            for a in AbelianStrata_g(g):
                yield a
