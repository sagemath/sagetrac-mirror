r"""
Strata of quadratic differentials on Riemann surfaces
More precisely, we are interested in meromorphic quadratic differentials
with at most simple poles on closed compact connected Riemann surfaces,
which are not globally the square of an abelian differential.

The moduli space of such quadratic differentials on Riemann surfaces
of a given genus is a complex orbifold, stratified by the degrees of
zeros (zeros of degree -1, or simple poles, being allowed). The strata
themselves are complex orbifolds. Most strata are connected but some
(infinitely many) are not.

A stratum corresponds to the Sage object
:class:`~sage.dynamics.flat_surfaces.strata.QuadraticStratum`.

The classification of connected components of strata of quadratic
differentials was established by Erwan Lanneau in [Lan08],
after a similar classification was established by Kontsevich
and Zorich in [KonZor03]_.

Each stratum has one or two connected components and each
component is associated to an extended Rauzy class. The
:meth:`~sage.dynamics.flat_surfaces.strata.QuadraticStratum.components`
method gives the decomposition of a stratum into its connected components.

A representative for each connected component of stratum is given by Zorich in [Zor08]_.

This is implemented here following [Zor08]_:

- genus zero stratum :meth:`~GZQSC.permutation_representative`

- genus one stratum :meth:`~GOQSC.permutation_representative`

- genus two hyperellitic component :meth:`~GTHQSC.permutation_representative`

- genus two non-hyperellitic component :meth:`~GTNQSC.permutation_representative`

- hyperelliptic component :meth:`~HQSC.permutation_representative`

- non-hyperelliptic component :meth:`~NQSC.permutation_representative`

- regular component of exceptional stratum :meth:`~REQSC.permutation_representative`

- irregular component of exceptional stratum :meth:`~IEQSC.permutation_representative`

The inverse operation, i.e., starting from a permutation, determine
the connected component it lives in, is partially written in [KonZor03]_.

See:
:meth:`~sage.dynamics.interval_exchanges.template.PermutationIET.connected_component`.

The code here implements the descriptions in [Zor08]_. Zorich already
implemented all this for Mathematica in [ZS]_.

See also :mod:`abelian_strata` for Abelian strata.

AUTHORS:

- Vincent Delecroix (2009-09-29): initial version

- Samuel Lelievre (2010-10-08): quadratic strata

EXAMPLES:

Construction of a stratum from a list of singularity degrees::

    sage: a = QuadraticStratum(2,2)
    sage: print a
    Q_2(2^2)
    sage: print a.genus()
    2

::

    sage: a = QuadraticStratum(4,3,2,2,1)
    sage: print a
    Q_4(4, 3, 2^2, 1)
    sage: print a.genus()
    4

By convention, the degrees are always written in decreasing order::

    sage: a1 = QuadraticStratum(7,5,3,1)
    sage: a1
    Q_5(7, 5, 3, 1)
    sage: a2 = QuadraticStratum(3,1,7,5)
    sage: a2
    Q_5(7, 5, 3, 1)
    sage: a1 == a2
    True

List the connected components of a stratum::

    sage: a = QuadraticStratum(6,2)
    sage: print a.components()
    [Q_3(6, 2)^hyp, Q_3(6, 2)^nonhyp]

::

    sage: a = QuadraticStratum(12)
    sage: cc = a.components()
    sage: print cc
    [Q_4(12)^reg, Q_4(12)^irr]
    sage: for c in cc:
    ...       print c; print c.permutation_representative()
    Q_4(12)^reg
    0 1 2 1 2 3 4 3 4 5
    5 6 7 6 7 0
    Q_4(12)^irr
    0 1 2 3 4 5 6 5
    7 6 4 7 3 2 1 0

::

    sage: a = QuadraticStratum(1, 1, 1, 1)
    sage: print a.components()
    [Q_2(1^4)^hyp]
    sage: c = a.components()[0]
    sage: p = c.permutation_representative(); p
    0 1 2 3 1 4 5
    2 6 5 4 6 3 0

.. TODO::

    implement connected components for generalized permutations
    iterators for quadratic strata with constraints on genus and/or dimension
    implement an "alphabet" option for the "permutation_representative" method
    include examples with fake zeros in documentation
    include an attribute "hyp" / "nonhyp" for components of strata?
    include __init__ for all classes of connected components? no

"""
#*****************************************************************************
#       Copyright (C) 2010 Vincent Delecroix <20100.delecroix@gmail.com>
#       Copyright (C) 2010 Samuel Lelievre <samuel.lelievre@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.infinity import Infinity
from sage.structure.parent import Parent

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.rings.rational import Rational


from strata import Stratum, StratumComponent,Strata

class QuadraticStratum(Stratum):
    r"""
    Stratum of quadratic differentials.

    EXAMPLES::

        sage: Q = QuadraticStratum(15,-1,-1,-1); Q
        Q_4(15, -1^3)
        sage: Q.components()
        [Q_4(15, -1^3)^c]

        sage: Q = QuadraticStratum(6,6); Q
        Q_4(6^2)
        sage: Q.components()
        [Q_4(6^2)^hyp, Q_4(6^2)^reg, Q_4(6^2)^irr]
    """
    _name = "Q"
    _latex_name = "\\mathcal{Q}"
    def __init__(self, *l, **kwds):
        """
        TESTS::

            sage: a = QuadraticStratum(-1,-1,-1,-1)
            sage: loads(dumps(a)) == a
            True
            sage: QuadraticStratum([])
            Traceback (most recent call last):
            ...
            ValueError: the list must be non empty
        """
        genus = kwds.get('genus', None)

        if len(l) == 1 and isinstance(l[0],(tuple,list,dict)):
            l = l[0]
        if not l:
            raise ValueError("the list must be nonempty")
        if isinstance(l, dict):
            l = sum(([v]*e for v,e in l.iteritems()),[])
        allzeros = map(Integer, l)
        allzeros.sort(reverse=True)
        self._nb_poles = allzeros.count(-1)
        self._nb_fake_zeros = allzeros.count(0)
        self._zeros = allzeros[:len(allzeros)-self._nb_fake_zeros-self._nb_poles]

        # nfz (for "non fake zeros")
        nfz = self._zeros + [-1] * self._nb_poles

        if genus is not None:
            g = sum(nfz) + 4
            if 4*genus > g:
                raise ValueError
            elif 4*genus < g:
                self._nb_poles += g - 4*genus
                nfz.extend([-1] * (g-4*genus))

        genus = sum(nfz) + 4
        if genus%4:
            raise ValueError("the sum of the zeros should sum up to an integer congruent to 0 mod 4")
        genus //= 4


        # Lanneau classification of connected components

        #TODO: check genus 2 components
        #TODO: in genus 2, decide between GTHQSC/GTNQSC and HQSC/NQSC
        if genus == 0:
            self._cc = (GZQSC,)

        #TODO: all genus 1 strata are connected, but two are hyperelliptic; give the component a different name then?
        elif genus == 1:
            if self._zeros == []:
                # print "The stratum Q(0) is empty!"
                self._cc = ()
            elif self._zeros == [1]:
                # print "The stratum Q(1,-1) is empty!"
                self._cc = ()
            else:
                self._cc = (GOQSC,)

        elif genus == 2:
            if self._zeros == [4]:
                # print "The stratum Q(4) is empty!"
                self._cc = ()
            elif self._zeros == [3,1]:
                # print "The stratum Q(3,1) is empty!"
                self._cc = ()
            elif self._zeros == [6]: self._cc = (HQSC,GTNQSC)
            elif self._zeros == [3,3]: self._cc = (HQSC,GTNQSC)
            elif self._zeros == [2,2]: self._cc = (GTHQSC,)
            elif self._zeros == [2,1,1]: self._cc = (GTHQSC,)
            elif self._zeros == [1,1,1,1]: self._cc = (GTHQSC,)
            else: self._cc = (GTNQSC,)

        elif genus == 3 and self._nb_poles == 1 and all((z == -1 or z%3 == 0) for z in self._zeros):
                self._cc = (REQSC, IEQSC)

        elif genus == 4 and self._nb_poles == 0 and all(z%3 == 0 for z in self._zeros):
            if self._zeros in [[12],[9,3]]:
                self._cc = (REQSC, IEQSC)
            elif self._zeros in [[6,6],[6,3,3],[3,3,3,3]]:
                self._cc = (HQSC, REQSC, IEQSC)
        else:
            if len(nfz) == 2 and nfz[0]%4 == 2 and nfz[1]%4 == 2:
                self._cc = (HQSC,NQSC)
            elif len(nfz) == 4 and nfz[0] == nfz[1] and nfz[2] == nfz[3] and nfz[0]%2 and nfz[2]%2:
                self._cc = (HQSC,NQSC)
            elif len(nfz) == 3 and nfz[0] == nfz[1] and nfz[0]%2 and nfz[2]%4 == 2:
                self._cc = (HQSC,NQSC)
            elif len(nfz) == 3 and nfz[1] == nfz[2] and nfz[1]%2 and nfz[0]%4 == 2:
                self._cc = (HQSC,NQSC)
            else:
                self._cc = (CQSC,)

    def zeros(self, poles=True, fake_zeros=True):
        r"""
        Returns the list of zeros of self.

        EXAMPLES::

            sage: AbelianStratum([1,2,3]).zeros()
            [3, 2, 1]
            sage: AbelianStratum({2:4}).zeros()
            [2, 2, 2, 2]

            sage: QuadraticStratum({-1:4}).zeros()
            [-1, -1, -1, -1]
            sage: QuadraticStratum({1:8}).zeros()
            [1, 1, 1, 1, 1, 1, 1, 1]
        """
        z = self._zeros[:]
        if fake_zeros:
            z.extend([0]*self._nb_fake_zeros)
        if poles:
            z.extend([-1]*self._nb_poles)
        return z


    def nb_zeros(self, poles=True, fake_zeros=True):
        r"""
        EXAMPLES::

            sage: QuadraticStratum({-1:4}).nb_zeros()
            4
            sage: QuadraticStratum({-1:4,1:4}).nb_zeros()
            8
        """
        n = Integer(len(self._zeros))
        if poles:
            n += self._nb_poles
        if fake_zeros:
            n += self._nb_fake_zeros
        return n

    def nb_poles(self):
        return self._nb_poles

    def nb_fake_zeros(self):
        return self._nb_fake_zeros

    def genus(self):
        r"""
        Returns the genus.

        EXAMPLES::

            sage: QuadraticStratum(-1,-1,-1,-1).genus()
            0
        """
        return Integer(sum(self.zeros(poles=True,fake_zeros=False))//4 + 1)

    def dimension(self):
        r"""
        Returns the complex dimension of this stratum.

        The complex dimension is the number of intervals minus one of any linear
        involution associated to this stratum.

        EXAMPLES::

            sage: QuadraticStratum({-1:4}).dimension()
            2
        """
        return 2 * self.genus() + self.nb_zeros() - 2

    def orientation_cover(self, fake_zeros=False):
        r"""
        Return the stratum of Abelian differentials which contains the set of
        orientation cover of quadratic differentials in this stratum.

        OPTIONS:

        - ``fake_zeros`` - boolean - if True, add fake zeros which corresponds
          to the double cover of poles.

        EXAMPLES::

            sage: q = QuadraticStratum({4:1,-1:4}); q
            Q_1(4, -1^4)
            sage: a1 = q.orientation_cover(); a1
            H_3(2^2)
            sage: q.orientation_cover(fake_zeros=True)
            H_3(2^2, 0^4)

        For hyperelliptic strata (orientation cover of quadratic strata of the
        form `Q(n,-1^{n+4})`) the dimension coincide. From [Lan08] we know that
        it only happens for those ones::

            sage: q = QuadraticStratum(4,genus=0); q
            Q_0(4, -1^8)
            sage: q.dimension()
            7
            sage: q.orientation_cover().dimension()
            7

            sage: q = QuadraticStratum(3,1,genus=0); q
            Q_0(3, 1, -1^8)
            sage: q.dimension()
            8
            sage: a = q.orientation_cover(); a
            H_4(4, 2)
            sage: a.dimension()
            9
        """
        from abelian_strata import AbelianStratum

        l = []
        for z in self.zeros():
            if z == -1:
                if fake_zeros:
                    l.append(0)
            elif z % 2:
                l.append(z+1)
            else:
                l.append(z//2)
                l.append(z//2)

        return AbelianStratum(l)

    def spin(self):
        r"""
        Return the spin structure (None, 0 or 1) of that component.

        Any quadratic differential has a canonic double cover (called the
        orientation cover) which is a surface, generally of higher genera, with
        an Abelian differential. The *spin structure* of the quadratic
        differential is the spin structure of that double cover.

        The spin is None if any degree of zero of the quadratic differential is
        congruent to 2 mod 4. Otherwise, denoting respectively k1 and k3 the
        number of degree of zero congruent to 1 and 3 modulo 4 we have

        ..MATH::

            spin = ((k1-k3)/4) mod 2

        The proof of that formula is the object of [La04]_.

        REFERENCES:

        .. [La04] E. Lanneau, "Spin"

        EXAMPLES::

            sage: QuadraticStratum(1,3).spin()
            0
            sage: QuadraticStratum(1,3,genus=1).spin()
            1
            sage: QuadraticStratum(2,2).spin() is None
            True
        """
        k = [0,0,0,0]
        for z in self.zeros():
            zz = z%4
            if zz == 2:
                return None
            k[zz] += 1

        return (abs(k[1]-k[3])//4) % 2

    def has_hyperelliptic_component(self):
        r"""
        Returns True if and only if self has a connected component which
        contains only hyperelliptic surfaces.

        EXAMPLES::

            sage: QuadraticStratum(2,2).has_hyperelliptic_component()
            True
            sage: QuadraticStratum(3,1).has_hyperelliptic_component()
            False
        """
        return HQSC in self._cc or GTHQSC in self._cc

    def hyperelliptic_component(self):
        r"""
        Returns the hyperelliptic component of self (or raise a ValueError).

        EXAMPLES::

            sage: QuadraticStratum(2,2).hyperelliptic_component()
            Q_2(2^2)^hyp

            sage: QuadraticStratum(3,1).hyperelliptic_component()
            Traceback (most recent call last):
            ...
            ValueError: the stratum has no hyperelliptic component
        """
        if HQSC in self._cc:
            return HQSC(self)
        if GTHQSC in self._cc:
            return GTHQSC(self)
        raise ValueError, "the stratum has no hyperelliptic component"

    def has_non_hyperelliptic_component(self):
        r"""
        Test whether this stratum has a non hyperelliptic component.

        EXAMPLES::

            sage: QuadraticStratum(10,10).has_non_hyperelliptic_component()
            True
            sage: QuadraticStratum(6,6).has_non_hyperelliptic_component()
            False
        """
        return NQSC in self._cc or GTNQSC in self._cc

    def non_hyperelliptic_component(self):
        r"""
        Returns the non hyperelliptic component of this stratum (or raise a
        ValueError).

        EXAMPLES::

            sage: QuadraticStratum(10,10).non_hyperelliptic_component()
            Q_6(10^2)^nonhyp
            sage: QuadraticStratum(2,2).non_hyperelliptic_component()
            Traceback (most recent call last):
            ...
            ValueError: no non hyperelliptic component
        """
        if NQSC in self._cc:
            return NQSC(self)
        elif GTNQSC in self._cc:
            return GTNQSC(self)
        raise ValueError, "no non hyperelliptic component"

    def has_regular_and_irregular_components(self):
        r"""
        Test whether this component has a pair of regular, irregular components.

        The list of strata of quadratic differentials that admit such a pair
        are:

        in genus 3: Q(9,-1), Q(6,3,-1), Q(3,3,3,-1)
        in genus 4: Q(12), Q(9,3), Q(6,6), Q(6,3,3), Q(3,3,3,3)

        EXAMPLES::

            sage: QuadraticStratum(9,-1).has_regular_and_irregular_components()
            True
            sage: QuadraticStratum(11,1).has_regular_and_irregular_components()
            False
        """
        return REQSC in self._cc

    def regular_component(self):
        r"""
        Returns the regular component of that stratum or raise a ValueError.

        EXAMPLES::

            sage: QuadraticStratum(12).regular_component()
            Q_4(12)^reg
            sage: QuadraticStratum(2,2,2,2).regular_component()
            Traceback (most recent call last):
            ...
            ValueError: no regular component for this stratum
        """
        if REQSC in self._cc:
            return REQSC(self)
        raise ValueError, "no regular component for this stratum"

    def irregular_component(self):
        r"""
        Returns the irregular component of that stratum or raise a ValueError.

        EXAMPLES::

            sage: QuadraticStratum(3,3,3,-1).irregular_component()
            Q_3(3^3, -1)^irr
            sage: QuadraticStratum(2,2).irregular_component()
            Traceback (most recent call last):
            ...
            ValueError: no irregular component for this stratum
        """
        if IEQSC in self._cc:
            return IEQSC(self)
        raise ValueError, "no irregular component for this stratum"

    def unique_component(self):
        r"""
        Returns the unique component of self or raise a ValueError.

        EXAMPLES::

            sage: QuadraticStratum({1:1, -1:5}).unique_component()
            Q_0(1, -1^5)^c
            sage: QuadraticStratum(3,2,-1).unique_component()
            Q_2(3, 2, -1)^nonhyp

            sage: QuadraticStratum(12).unique_component()
            Traceback (most recent call last):
            ...
            ValueError: several components for this stratum
        """
        if len(self._cc) != 1:
            raise ValueError, "several components for this stratum"
        return self._cc[0](self)

    def random_cylindric_permutation(self):
        r"""
        Return a random cylindric permutation that belongs to this stratum.

        EXAMPLES::

            sage: Q = QuadraticStratum(4,4)
            sage: Q.random_cylindric_permutation()  # random
            0 1 2 3 4 5 1 5
            2 6 3 6 4 0
        """
        return self.random_component().random_cylindric_permutation()


    def random_cylindric_permutation(self):
        r"""
        Return a random cylindric permutation that belongs to this stratum.
        """
        return self.random_component().random_cylindric_permutation()

class QuadraticStratumComponent(StratumComponent):
    r"""
    Generic class for component of quadratic stratum.
    """
    _name = 'c'

    def random_cylindric_permutation(self, nsteps=64):
        r"""
        Return a cylindric permutation of the form ``p = ((0,...),(..., 0))``
        where 0 can be any label.

        EXAMPLES::

            sage: Q = QuadraticStratum({4:1,-1:4}); Q
            Q_1(4, -1^4)
            sage: c = Q.components()[0]
            sage: p = c.random_cylindric_permutation(); p   # random
            0 1 2 3 3 1
            2 4 4 5 5 0
            sage: p.stratum()
            Q_1(4, -1^4)
            sage: p = c.random_cylindric_permutation(); p
            0 1 2 2
            3 4 4 1 5 5 3 0
            sage: p.stratum()
            Q_1(4, -1^4)

            sage: Q = QuadraticStratum(6,6); Q
            Q_4(6^2)
            sage: c_hyp, c_reg, c_irr = Q.components()

            sage: for _ in xrange(3):
            ...      p = c_reg.random_cylindric_permutation()
            ...      print p
            ...      print p.stratum()
            0 1 2 3 4 5 6 7 3
            1 4 6 2 8 7 8 5 0
            Q_4(6^2)
            0 1 2 3 4 2 5 6
            6 4 7 8 3 7 5 1 8 0
            Q_4(6^2)
            0 1 2 3 4 5 6 3 7 6
            4 1 8 2 5 7 8 0
            Q_4(6^2)
        """
        import sage.misc.prandom as prandom

        p = self.permutation_representative()

        for _ in xrange(nsteps):
            while not p.has_rauzy_move(0):
                p = p.rauzy_move(1)
            while not p.has_rauzy_move(1):
                p = p.rauzy_move(0)

            rd = prandom.random()
            if rd < 0.1:  # inplace (symmetric)
                p._inversed_twin()
                p._reversed_twin()
            elif rd < 0.55:
                p = p.rauzy_move(0)
            else:
                p = p.rauzy_move(1)

        while p[0][0] != p[1][-1]:
            rd = prandom.random()
            if rd < 0.1:  # inplace (symmetric)
                p._inversed_twin()
                p._reversed_twin()
            elif not p.has_rauzy_move(0):
                p = p.rauzy_move(1)
            elif not p.has_rauzy_move(1) or rd < .55:
                p = p.rauzy_move(0)
            else:
                p = p.rauzy_move(1)

        return p

    def orientation_cover_component(self, fake_zeros=False):
        r"""
        Return the connected component of Abelian stratum component which
        contains the set of orientation cover of quadratic differentials in this
        connected component.

        ALGORITHM:

        The spin only depends on the component and the only components for which
        the double cover belongs to a hyperelliptic component are
        Q(k,-1^k+4) and Q(2g-1,2g-1,-1,-1)^hyp

        OPTIONS:

        - ``fake_zeros`` - boolean - if True, add fake zeros which corresponds
          to the double cover of poles.

        EXAMPLES::

            sage: cc = QuadraticStratum(5,genus=0).unique_component()
            sage: cc.orientation_cover_component()
            H_4(6)^hyp

            sage: cc = QuadraticStratum(5,genus=1).unique_component()
            sage: cc.orientation_cover_component()
            H_4(6)^odd

            sage: cc = QuadraticStratum(5,genus=2).unique_component()
            sage: cc.orientation_cover_component()
            H_4(6)^even

            sage: cc = QuadraticStratum(1,1,genus=0).unique_component()
            sage: cc.orientation_cover_component()
            H_3(2^2)^odd

            sage: cc = QuadraticStratum(4,genus=0).unique_component()
            sage: cc.orientation_cover_component()
            H_3(2^2)^hyp

            sage: cc = QuadraticStratum(1,1,genus=1).unique_component()
            sage: cc.orientation_cover_component()
            H_3(2^2)^hyp
        """
        stratum = self.stratum()
        zeros = stratum.zeros(poles=False, fake_zeros=False)
        astratum = stratum.orientation_cover(fake_zeros=fake_zeros)

        # 0) easy case: the stratum is connected
        if astratum.is_connected():
            return astratum.unique_component()

        # 1) canonical construction of Abelian hyperelliptic components
        if len(zeros) == 1 and stratum.genus() == 0:
                return astratum.hyperelliptic_component()

        # 2) other possibility for the hyperelliptic components (low genus
        # special case for Q(1,1,-1,-1))
        if zeros == [1,1] and stratum.nb_poles() == 2:
            return astratum.hyperelliptic_component()

        # 3) composant with spin
        spin = stratum.spin()
        if spin == 0:
            return astratum.even_component()
        elif spin == 1:
            return astratum.odd_component()

        # 4) as connected components are treated in 0 there remains only non
        # hyperelliptic components
        return astratum.non_hyperelliptic_component()

    def random_cylindric_permutation(self, nsteps=64):
        r"""
        Return a cylindric permutation of the form ``p = ((0,...),(..., 0))``
        where 0 can be any label.

        EXAMPLES::

            sage: Q = QuadraticStratum({4:1,-1:4}); Q
            Q_1(4, -1^4)
            sage: c = Q.components()[0]
            sage: p = c.random_cylindric_permutation()
            XX
            YY
            sage: p.stratum()
            H(2^2)
            sage: p = c.random_cylindric_permutation()
            XX
            YY
            sage: p.stratum()
            H(2^2)

            sage: Q = QuadraticStratum(12); Q
            Q_4(12)
            sage: c_hyp, c_reg, c_irr = Q.components()

            sage: for _ in xrange(3):
            ...      p = c_reg.random_cylindric_permutation()
            ...      print p
            ...      print p.stratum()
            XX
            YY
            Q_4(12)
            XX
            YY
            Q_4(12)
            XX
            YY
            Q_4(12)
        """
        import sage.misc.prandom as prandom

        p = self.permutation_representative()

        for _ in xrange(nsteps):
            while not p.has_rauzy_move(0):
                p = p.rauzy_move(1)
            while not p.has_rauzy_move(1):
                p = p.rauzy_move(0)

            rd = prandom.random()
            if rd < 0.1:  # inplace (symmetric)
                p._inversed_twin()
                p._reversed_twin()
            elif rd < 0.55:
                p = p.rauzy_move(0)
            else:
                p = p.rauzy_move(1)

        while p[0][0] != p[1][-1]:
            rd = prandom.random()
            if rd < 0.1:  # inplace (symmetric)
                p._inversed_twin()
                p._reversed_twin()
            elif not p.has_rauzy_move(0):
                p = p.rauzy_move(1)
            elif not p.has_rauzy_move(1) or rd < .55:
                p = p.rauzy_move(0)
            else:
                p = p.rauzy_move(1)

        return p

QSC = QuadraticStratumComponent

class GenusZeroQuadraticStratumComponent(QSC):
    r"""
    This class is intended to be called internally rather than directly.
    Call only with appropriate parameters, in particular correct genus: no
    consistency check inside, no prediction as to what may happen otherwise.
    """
    def permutation_representative(self, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns a generalized permutation representative.

        EXAMPLES::

            sage: QuadraticStratum(-1,-1,-1,-1).permutation_representative()
            0 1 1
            2 2 0
            sage: QuadraticStratum(1,genus=0).permutation_representative()
            0 1 1
            2 2 3 3 4 4 0
            sage: QuadraticStratum(2,genus=0).permutation_representative()
            0 1 1
            2 2 3 3 4 4 5 5 0
            sage: QuadraticStratum([1,1],genus=0).permutation_representative()
            0 1 1
            2 2 3 4 4 5 5 3 6 6 0
            sage: QuadraticStratum([2,1],genus=0).permutation_representative()
            0 1 1
            2 2 3 3 4 5 5 6 6 4 7 7 0
        """
        p = self._stratum.nb_poles()
        f = self._stratum.nb_fake_zeros()
        z = self._stratum.zeros(fake_zeros=False,poles=False)

        ll = [0]
        if f: ll = map(lambda x:'0'+str(x),range(f+1))

        l0 = ll + [2*p-6,2*p-6]
        l1 = []
        for k in xrange(1,p-3):
            l1.extend([2*k,2*k,2*k+1])
        for k in xrange(len(z)):
            for i in xrange(sum(z[:k])+1,sum(z[:k+1])):
                del l1[2*i+k]
        l2 = [2 * sum(z[:k+1]) + 1 for k in xrange(len(z))]
        l2.reverse()
        l1.extend(l2)
        l1.extend([1,1])
        l1.extend(ll)

        if reduced == True:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationLI
            p = ReducedPermutationLI([l0,l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationLI
            p = LabelledPermutationLI([l0,l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

GZQSC = GenusZeroQuadraticStratumComponent

class GenusOneQuadraticStratumComponent(QSC):
    r"""
    This class is intended to be called internally rather than directly.
    Call only with appropriate parameters, in particular correct genus: no
    consistency check inside, no prediction as to what may happen otherwise.
    """
    def permutation_representative(self, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns a generalized permutation representative.

        EXAMPLES::

            sage: QuadraticStratum([2],genus=1).permutation_representative()
            0 1 2 2
            1 3 3 0
            sage: QuadraticStratum([3],genus=1).permutation_representative()
            0 1 2 2 3 3
            1 4 4 0
            sage: QuadraticStratum([1,1],genus=1).permutation_representative()
            0 1 2 3 3
            2 1 4 4 0
            sage: Q = QuadraticStratum([2,1],genus=1)
            sage: Q.permutation_representative(alphabet='abcdef')
            a b c c d e e
            d b f f a

        TESTS::

            sage: QuadraticStratum([1],genus=1).permutation_representative()
            Traceback (most recent call last):
            ...
            EmptySetError: The stratum is empty
        """
        p = self._stratum.nb_poles()
        f = self._stratum.nb_fake_zeros()
        z = self._stratum.zeros(fake_zeros=False,poles=False)

        if f: ll = map(lambda x:'0'+str(x),range(f+1))
        else: ll = [0]

        l0 = ll + [1]
        for k in xrange(1,p):
            l0.extend([2*k,2*k+1,2*k+1])
        l1 = [2 * (sum(z[:k])) for k in xrange(1,len(z))]
        l1.extend([1, 2*p, 2*p])
        l1.extend(ll)

        for k in xrange(len(z)):
            for i in xrange(sum(z[:k])+1,sum(z[:k+1])):
                del l0[2*i+k]

        if reduced == True:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationLI
            p = ReducedPermutationLI([l0,l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationLI
            p = LabelledPermutationLI([l0,l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))

        return p


GOQSC = GenusOneQuadraticStratumComponent

class GenusTwoHyperellipticQuadraticStratumComponent(QSC):
    r"""
    This class is intended to be called internally rather than directly.
    Call only with appropriate parameters, in particular correct genus: no
    consistency check inside, no prediction as to what may happen otherwise.
    """
    _name = 'hyp'

    def permutation_representative(self, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns a generalized permutation representative.

        TESTS:

            sage: Q = QuadraticStratum([1,1,1,1],genus=2)
            sage: H = Q.hyperelliptic_component()
            sage: H.permutation_representative()
            0 1 2 3 1 4 5
            2 6 5 4 6 3 0

            sage: Q = QuadraticStratum([2,1,1],genus=2)
            sage: H = Q.hyperelliptic_component()
            sage: H.permutation_representative(alphabet='abcdef')
            a b c b d e
            f e d f c a


            sage: Q = QuadraticStratum([2,2],genus=2)
            sage: H = Q.hyperelliptic_component()
            sage: H.permutation_representative()
            0 1 2 1 3
            4 3 4 2 0
        """
        p = self._stratum.nb_poles()
        f = self._stratum.nb_fake_zeros()
        z = self._stratum.zeros(fake_zeros=False,poles=False)

        if f: ll = map(lambda x:'0'+str(x),range(f+1))
        else: ll = [0]

        if z == [2,2]:
            l0 = ll + [6, 5, 6, 4]
            l1 = [2, 4, 2, 5]
        elif z == [2,1,1]:
            l0 = ll + [6, 5, 6, 4, 3]
            l1 = [2, 3, 4, 2, 5]
        elif z == [1,1,1,1]:
            l0 = ll + [6, 1, 5, 6, 4, 3]
            l1 = [1, 2, 3, 4, 2, 5]
        else:
            raise ValueError, "Wrong stratum for GenusTwoHyperellipticQuadraticStratumComponent"

        l1.extend(ll)

        if reduced == True:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationLI
            p = ReducedPermutationLI([l0,l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationLI
            p = LabelledPermutationLI([l0,l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))

        return p



GTHQSC = GenusTwoHyperellipticQuadraticStratumComponent

#TODO: optimize construction of l0 and l1
class GenusTwoNonhyperellipticQuadraticStratumComponent(QSC):
    r"""
    This class is intended to be called internally rather than directly.
    Call only with appropriate parameters, in particular correct genus: no
    consistency check inside, no prediction as to what may happen otherwise.
    """
    _name = 'nonhyp'

    def permutation_representative(self, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns a generalized permutation representative.

        TESTS::

            sage: Q = QuadraticStratum([6],genus=2)
            sage: N = Q.non_hyperelliptic_component()
            sage: N.permutation_representative()
            0 1 2 1 3
            3 4 4 5 5 2 0
            sage: Q = QuadraticStratum([3,3],genus=2)
            sage: N = Q.non_hyperelliptic_component()
            sage: N.permutation_representative()
            0 1 2 1 3
            4 3 5 5 4 6 6 2 0
        """
        p = self._stratum.nb_poles()
        f = self._stratum.nb_fake_zeros()
        z = self._stratum.zeros(fake_zeros=False,poles=False)

        if z == [4]:
            raise ValueError("The stratum Q(4) is empty!")

        if z == [3,1]:
            raise ValueError("The stratum Q(3,1) is empty!")

        if z in [[2,2],[2,1,1],[1,1,1,1]]:
            raise ValueError("This stratum has no non-hyperelliptic component!")

        if f: ll = map(lambda x:'0'+str(x),range(f+1))
        else: ll = [0]

        l0 = ll + [2*p+6,1,2*p+5,2*p+6,p+4,p+3]
        l1 = range(1,p+5)
        for k in xrange(p+5,2*p+5):
            l1.extend([2*p+7-k,k,k])
        l1.extend([2,2*p+5])
        l1.extend(ll)

        for k in xrange(len(z)):
            for i in xrange(sum(z[:k])+1,sum(z[:k+1])):
                if i in l0: del l0[l0.index(i)]
                if i in l0: del l0[l0.index(i)]
                if i in l1: del l1[l1.index(i)]
                if i in l1: del l1[l1.index(i)]

        if reduced == True:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationLI
            p = ReducedPermutationLI([l0,l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationLI
            p = LabelledPermutationLI([l0,l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

GTNQSC = GenusTwoNonhyperellipticQuadraticStratumComponent

class HyperellipticQuadraticStratumComponent(QSC):
    r"""
    This class is intended to be called internally rather than directly.
    Call only with appropriate parameters, in particular correct genus: no
    consistency check inside, no prediction as to what may happen otherwise.
    """
    _name = 'hyp'

    def permutation_representative(self, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns a generalized permutation representative.

        EXAMPLES::

            sage: cc = QuadraticStratum([6],genus=2).hyperelliptic_component()
            sage: cc.permutation_representative()
            0 1 2 3 4 1
            5 4 3 2 5 0
            sage: cc = QuadraticStratum([3,3],genus=2).hyperelliptic_component()
            sage: cc.permutation_representative()
            0 1 2 3 4 5 1
            6 5 4 3 2 6 0
            sage: cc = QuadraticStratum(10,10).hyperelliptic_component()
            sage: cc.permutation_representative()
            0 1 2 3 4 5 6 1 7 8 9 10 11
            11 10 9 8 7 12 6 5 4 3 2 12 0
        """
        p = self._stratum.nb_poles()
        f = self._stratum.nb_fake_zeros()
        z = self._stratum.zeros(fake_zeros=False,poles=True)

        if f: ll = map(lambda x:'0'+str(x),range(f+1))
        else: ll = [0]

        if len(z) == 2: r = z[0]//2; s = z[1]//2
        elif len(z) == 4: r = z[0]+1; s = z[2]+1
        elif len(z) == 3 and z[0]%2: r = z[0]+1; s = z[2]//2
        elif len(z) == 3: r = z[0]//2; s = z[2]+1
        else: raise ValueError, "This stratum has no hyperelliptic component!"

        l0 = ll + ['A']
        l0.extend(range(1,r+1))
        l0.append('A')
        l0.extend(range(r+1,r+s+1))
        l1 = range(r+s,r,-1)
        l1.append('B')
        l1.extend(range(r,0,-1))
        l1.append('B')
        l1.extend(ll)

        if reduced == True:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationLI
            p = ReducedPermutationLI([l0,l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationLI
            p = LabelledPermutationLI([l0,l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

HQSC = HyperellipticQuadraticStratumComponent

class ConnectedQuadraticStratumComponent(QSC):
    r"""
    Connected component of stratum of quadratic differentials.


    This class is intended to be called internally rather than directly.
    Call only with appropriate parameters, in particular correct genus: no
    consistency check inside, no prediction as to what may happen otherwise.
    """
    _name = 'c'

    def permutation_representative(self, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns a generalized permutation representative.

        NOTES:

        The representative is made by constructing two lists l0 and l1
        which correspond to a generalized permutation representative
        for the stratum with simple zeros, and then erasing some
        elements from l0 and l1 (this corresponds to collapsing saddle
        connections to merge zeros). It may be possible to find a
        faster way to obtain the desired l0 and l1 than by constructing
        the long versions and erasing symbols; the current implementation
        has a loop with "del l0[l0.index(i)]" and "del l1[l1.index(i)]".

        EXAMPLES::

            sage: cc = QuadraticStratum(6,-1,-1).non_hyperelliptic_component()
            sage: p = cc.permutation_representative(); p
            0 1 2 1 3
            3 4 4 5 5 2 0
            sage: p.stratum_component()
            Q_2(6, -1^2)^nonhyp

            sage: cc = QuadraticStratum([3,3],genus=2).non_hyperelliptic_component()
            sage: p = cc.permutation_representative(); p
            0 1 2 1 3
            4 3 5 5 4 6 6 2 0
            sage: p.stratum_component()
            Q_2(3^2, -1^2)^nonhyp

            sage: cc = QuadraticStratum(8).unique_component()
            sage: p = cc.permutation_representative(); p
            0 1 2 1 2 3
            4 5 4 5 3 0
            sage: p.stratum_component()
            Q_3(8)^c

            sage: Q = QuadraticStratum(4,4).unique_component()
            sage: p = Q.permutation_representative(); p
            0 1 2 3 2 3 4
            5 6 5 6 1 4 0
            sage: p.stratum_component()
            Q_3(4^2)^c

            sage: Q = QuadraticStratum({12:1,-1:4}).unique_component()
            sage: p = Q.permutation_representative()
            sage: p
            0 1 2 1 2 3 3 4 4 5 5 6 6 7
            8 9 8 9 7 0
            sage: p.stratum()
            Q_3(12, -1^4)
        """
        p = self._stratum.nb_poles()
        f = self._stratum.nb_fake_zeros()
        z = self._stratum.zeros(fake_zeros=False,poles=False)
        g = (sum(z)-p)//4+1

        if f: ll = map(lambda x:'0'+str(x),range(f+1))
        else: ll = [0]

        l0 = ll + [1,2,3]
        l1 = [4*g-3+p,3,4*g-2+p,2,4*g-3+p,1,4*g-2+p]

        for k in xrange(g-2):  # V(k) in Zorich 2008
            l0.extend([4+4*k,4*g-1+p+2*k,7+4*k,4*g+p+2*k,6+4*k,4*g-1+p+2*k,5+4*k,4*g+p+2*k])
        for l in xrange(p):    # W(l) in Zorich 2008
            l0.extend([4*g-4+l,6*g+p-5+l,6*g+p-5+l])
        l0.append(4*g-4+p)
        l1.extend(xrange(4,4*g-3+p))
        l1.extend(ll)

        for k in xrange(len(z)):
            for i in xrange(sum(z[:k])+1,sum(z[:k+1])):
                del l0[l0.index(i)]
                del l1[l1.index(i)]

        if reduced == True:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationLI
            p = ReducedPermutationLI([l0,l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationLI
            p = LabelledPermutationLI([l0,l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

CQSC = ConnectedQuadraticStratumComponent

class NonhyperellipticQuadraticStratumComponent(CQSC):
    r"""
    Non hyperelliptic component of stratum of quadratic differentials.
    """
    _name = 'nonhyp'

NQSC = NonhyperellipticQuadraticStratumComponent

class RegularExceptionalQuadraticStratumComponent(QSC):
    r"""
    This class is intended to be called internally rather than directly.
    Call only with appropriate parameters, in particular correct genus: no
    consistency check inside, no prediction as to what may happen otherwise.
    """
    _name = 'reg'

    def permutation_representative(self, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns a generalized permutation representative.

        EXAMPLES::

            sage: cc = QuadraticStratum(9,-1).regular_component()
            sage: p = cc.permutation_representative(); p
            0 1 2 1 2 3 3 4
            5 6 5 6 4 0
            sage: p.stratum_component()
            Q_3(9, -1)^reg

            sage: cc = QuadraticStratum(6,3,-1).regular_component()
            sage: p = cc.permutation_representative(); p
            0 1 2 3 1 2 4 4 5
            6 7 6 7 3 5 0
            sage: p.stratum_component()
            Q_3(6, 3, -1)^reg

            sage: cc = QuadraticStratum(12).regular_component()
            sage: p = cc.permutation_representative(); p
            0 1 2 1 2 3 4 3 4 5
            5 6 7 6 7 0
            sage: p.stratum_component()
            Q_4(12)^reg
        """

        p = self._stratum.nb_poles()
        f = self._stratum.nb_fake_zeros()
        z = self._stratum.zeros(fake_zeros=False,poles=True)

        if f: ll = map(lambda x:'0'+str(x),range(f+1))
        else: ll = [0]

        if z == [12]:
            l0 = ll + [1, 2, 1, 2, 3, 4, 3, 4, 5]
            l1 = [5, 6, 7, 6, 7]
        elif z == [9, 3]:
            l0 = ll + [1, 2, 3, 4, 2, 'A', 5, 6]
            l1 = [1, 4, 5, 7, 6, 7, 'A', 3]
        elif z == [6, 6]:
            l0 = ll +[1, 2, 3, 4, 'A', 2, 5, 6, 'A']
            l1 = [1, 4, 5, 7, 6, 7, 3]
        elif z == [6, 3, 3]:
            l0 = ll + [1, 2, 3, 'B', 4, 2, 'A', 5, 6]
            l1 = [1, 4, 5, 7, 6, 7, 'B', 'A', 3]
        elif z == [3, 3, 3, 3]:
            l0 = ll + [1, 2, 3, 'B', 4, 2, 'A', 5, 6]
            l1 = [1, 4, 5, 7, 6, 'C', 7, 'C', 'B', 'A', 3]
        elif z == [9, -1]:
            l0 = ll + [1, 2, 1, 2, 3, 3, 4]
            l1 = [5, 6, 5, 6, 4]
        elif z == [6, 3, -1]:
            l0 = ll + [1, 2, 3, 1, 2, 4, 4, 5]
            l1 = [6, 7, 6, 7, 3, 5]
        elif z == [3, 3, 3, -1]:
            l0 = ll + [1, 2, 3, 4, 2, 3, 5, 5, 6]
            l1 = [7, 1, 8, 7, 8, 4, 6]
        else:
            raise ValueError, "RegularExceptionalQuadraticStratumComponent applies only to the 4 exceptional strata"

        l1.extend(ll)

        if reduced == True:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationLI
            p = ReducedPermutationLI([l0,l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationLI
            p = LabelledPermutationLI([l0,l1])

        if alphabet:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))

        return p

REQSC = RegularExceptionalQuadraticStratumComponent

class IrregularExceptionalQuadraticStratumComponent(QSC):
    r"""
    This class is intended to be called internally rather than directly.
    Call only with appropriate parameters, in particular correct genus: no
    consistency check inside, no prediction as to what may happen otherwise.
    """
    _name = 'irr'

    def permutation_representative(self, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns a generalized permutation representative.

        EXAMPLES::

            sage: cc = QuadraticStratum(9,-1).irregular_component()
            sage: p = cc.permutation_representative(); p
            0 1 2 3 4 1 2 3 4 5
            5 6 6 0
            sage: p.stratum_component()  # optional
            Q_3(9, -1)^irr

            sage: cc = QuadraticStratum(6,3,-1).irregular_component()
            sage: p = cc.permutation_representative(); p
            0 1 2 3 4 5 1 2 3 4 5 6
            6 7 7 0
            sage: p.stratum_component()
            Q_3(6, 3, -1)^irr

            sage: cc = QuadraticStratum(12).irregular_component()
            sage: p = cc.permutation_representative(); p
            0 1 2 3 4 5 6 5
            7 6 4 7 3 2 1 0
            sage: p.stratum_component()
            Q_4(12)^irr
        """

        p = self._stratum.nb_poles()
        f = self._stratum.nb_fake_zeros()
        z = self._stratum.zeros(fake_zeros=False,poles=True)

        if f: ll = map(lambda x:'0'+str(x),range(f+1))
        else: ll = [0]

        if z == [12]:
            l0 = ll + [1, 2, 3, 4, 5, 6, 5]
            l1 = [7, 6, 4, 7, 3, 2, 1]
        elif z == [9, 3]:
            l0 = ll + [1, 2, 3, 4, 'A', 3, 'A', 5, 6]
            l1 = [1, 5, 7, 4, 2, 6, 7]
        elif z == [6, 6]:
            l0 = ll + [1, 2, 3, 4, 3, 'A', 5, 6]
            l1 = [1, 5, 7, 4, 2, 6, 'A', 7]
        elif z == [6, 3, 3]:
            l0 = ll + [1, 2, 'B', 3, 4, 'A', 3, 'A', 5, 6]
            l1 = [1, 5, 7, 'B', 4, 2, 6, 7]
        elif z == [3, 3, 3, 3]:
            l0 = ll + [1, 2, 'B', 3, 4, 'A', 3, 'A', 5, 'C', 6]
            l1 = [1, 5, 7, 'B', 4, 2, 6, 'C', 7]
        elif z == [9, -1]:
            l0 = ll + [1, 2, 3, 4, 1, 2, 3, 4, 5]
            l1 = [5, 6, 6]
        elif z == [6, 3, -1]:
            l0 = ll + [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6]
            l1 = [6, 7, 7]
        elif z == [3, 3, 3, -1]:
            l0 = ll + [1, 2, 3, 4, 5, 1, 6, 2, 3, 4, 5, 6, 7]
            l1 = [7, 8, 8]
        else:
            raise ValueError, "IrregularExceptionalQuadraticStratumComponent applies only to the 4 exceptional strata"

        l1.extend(ll)

        if reduced == True:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationLI
            p = ReducedPermutationLI([l0,l1])

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationLI
            p = LabelledPermutationLI([l0,l1])

        if alphabet:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

IEQSC = IrregularExceptionalQuadraticStratumComponent


def QuadraticStrata(genus=None, dimension=None, min_nb_poles=None, max_nb_poles=None,nb_poles=None):
    r"""
    Quadratic strata.

    INPUT:

    - ``genus`` - a non negative integer or None

    - ``dimension`` - a non negative integer or None

    - ``min_nb_poles``, ``max_nb_poles`` - the minimum and maximum number of
      poles allowed

    - ``nb_poles`` - integer - the number of poles (if the option is set then
      the options ``min_nb_poles`` and ``max_nb_poles`` are ignored)

    EXAMPLES::

        sage: Q = QuadraticStrata(genus=2); Q
        Quadratic strata of genus 2 surfaces
        sage: Q.cardinality()
        +Infinity
        sage: i = iter(Q)
        sage: i.next()
        Q_2(2^2)
        sage: i.next()
        Q_2(2, 1^2)
        sage: i.next()
        Q_2(1^4)
        sage: i.next()
        Q_2(5, -1)

        sage: Q = QuadraticStrata(dimension=5); Q
        Quadratic strata of dimension 5

        sage: Q = QuadraticStrata(genus=3,max_nb_poles=6); Q
        Quadratic strata of genus 3 surfaces with at most 6 poles
        sage: Q.cardinality()
        463
        sage: for q in Q: print q
        Q_3(8)
        Q_3(7, 1)
        Q_3(6, 2)
        ...
        Q_3(2^2, 1^10, -1^6)
        Q_3(2, 1^12, -1^6)
        Q_3(1^14, -1^6)

        sage: Q = QuadraticStrata(genus=2,nb_poles=0); Q
        Quadratic strata of genus 2 surfaces with no pole
        sage: for q in Q: print q
        Q_2(2^2)
        Q_2(2, 1^2)
        Q_2(1^4)


        sage: Q = QuadraticStrata(dimension=7,min_nb_poles=1,max_nb_poles=3); Q
        Quadratic strata of dimension 7 with at least 1 and at most 3 poles
        sage: for q in Q: print q
        Q_3(8, 1, -1)
        Q_3(7, 2, -1)
        Q_3(6, 3, -1)
        Q_3(5, 4, -1)
        Q_2(2, 1^3, -1)
        Q_3(10, -1^2)
        Q_2(4, 1^2, -1^2)
        Q_2(3, 2, 1, -1^2)
        Q_2(2^3, -1^2)
        Q_2(6, 1, -1^3)
        Q_2(5, 2, -1^3)
        Q_2(4, 3, -1^3)

        sage: QuadraticStrata(dimension=6,genus=0)
        Quadratic strata of genus 0 surfaces and dimension 6
    """
    if nb_poles is not None:
        min_nb_poles = max_nb_poles = Integer(nb_poles)
    else:
        if min_nb_poles is None:
            min_nb_poles = Integer(0)
        else:
            min_nb_poles = Integer(min_nb_poles)

        if max_nb_poles is None:
            max_nb_poles = Infinity
        else:
            if max_nb_poles == Infinity:
                max_nb_poles = Infinity
            else:
                max_nb_poles = Integer(max_nb_poles)

        if max_nb_poles < min_nb_poles:
            raise ValueError, "min_nb_poles should be less or equal than max_nb_poles"

    if genus is None:
        if dimension is None:
            return QuadraticStrata_all()
        dimension = Integer(dimension)
        return QuadraticStrata_d(dimension,min_nb_poles,max_nb_poles)

    genus = Integer(genus)

    if dimension is None:
        return QuadraticStrata_g(genus,min_nb_poles,max_nb_poles)

    dimension = Integer(dimension)

    return QuadraticStrata_gd(genus,dimension,min_nb_poles,max_nb_poles)

class QuadraticStrata_class(Strata):
    r"""
    Generic class for strata of quadratic differentials.
    """
    def _repr_(self):
        r"""
        TESTS::

            sage: repr(AbelianStrata(genus=3))   #indirect doctest
            'Abelian strata of genus 3 surfaces'
        """
        s = self._repr_base_()

        if self._min_nb_poles != 0 and self._max_nb_poles != Infinity:
            if self._min_nb_poles == self._max_nb_poles:
                s += " with %d pole"%(self._min_nb_poles)
                if self._min_nb_poles != 1:
                    s += 's'
            else:
                s += " with at least %d and at most %d poles"%(self._min_nb_poles,self._max_nb_poles)

        elif self._min_nb_poles != 0:
            s += " with at least %d pole"%(self._min_nb_poles)
            if self._min_nb_poles != 1:
                s += 's'

        elif self._max_nb_poles == 0:
            s += " with no pole"

        elif self._max_nb_poles != Infinity:
            s += " with at most %d pole"%(self._max_nb_poles)
            if self._max_nb_poles != 1:
                s += 's'

        return s

    def __ne__(self, other):
        return not self.__eq__(other)

class QuadraticStrata_g(QuadraticStrata_class):
    r"""
    Stratas of genus g surfaces.

    INPUT:

    - ``genus`` - a non negative integer

    - ``min_nb_poles``, ``max_nb_poles`` - the number of poles

    EXAMPLES::

        sage: Q = QuadraticStrata(genus=3); Q
        Quadratic strata of genus 3 surfaces
        sage: Q.cardinality()
        +Infinity
        sage: i = iter(Q)
        sage: i.next()
        Q_3(8)

        sage: Q = QuadraticStrata(genus=2,max_nb_poles=1); Q
        Quadratic strata of genus 2 surfaces with at most 1 pole
        sage: Q.list()
        [Q_2(2^2), Q_2(2, 1^2), Q_2(1^4), Q_2(5, -1), Q_2(4, 1, -1), Q_2(3, 2, -1), Q_2(3, 1^2, -1), Q_2(2^2, 1, -1), Q_2(2, 1^3, -1), Q_2(1^5, -1)]

        sage: Q = QuadraticStrata(genus=4,nb_poles=3); Q
        Quadratic strata of genus 4 surfaces with 3 poles
        sage: Q.cardinality()
        176
    """
    def __init__(self,genus,min_nb_poles=None,max_nb_poles=None):
        r"""
        TESTS::

            sage: s = QuadraticStrata(genus=3)
            sage: TestSuite(s).run()
            sage: loads(dumps(s)) == s
            True

            sage: s = QuadraticStrata(genus=3,nb_poles=5)
            sage: TestSuite(s).run()
            sage: loads(dumps(s)) == s
            True

            sage: s = QuadraticStrata(genus=3,min_nb_poles=2,max_nb_poles=6)
            sage: TestSuite(s).run()
            sage: loads(dumps(s)) == s
            True
        """
        self._genus = genus
        self._min_nb_poles = min_nb_poles
        self._max_nb_poles = max_nb_poles

        if max_nb_poles != Infinity:
            Parent.__init__(self, category=FiniteEnumeratedSets())
        else:
            Parent.__init__(self,category=InfiniteEnumeratedSets())

    def __eq__(self, other):
        r"""
        Equality test.

        TESTS::

            sage: QuadraticStrata(genus=1) == QuadraticStrata(genus=1)
            True
            sage: QuadraticStrata(genus=1) == QuadraticStrata(genus=2)
            False
            sage: QuadraticStrata(genus=1) == QuadraticStrata(genus=1,min_nb_poles=5)
            False
            sage: QuadraticStrata(genus=1) == QuadraticStrata(genus=1,max_nb_poles= 3)
            False
        """
        return (isinstance(other, QuadraticStrata_g) and
                self._genus == other._genus and
                self._min_nb_poles == other._min_nb_poles and
                self._max_nb_poles == other._max_nb_poles)

    def __reduce__(self):
        r"""
        Support for pickling.

        EXAMPLES::

            sage: QuadraticStrata(genus=1).__reduce__()
            (<class 'sage.dynamics.flat_surfaces.quadratic_strata.QuadraticStrata_g'>, (1, 0, +Infinity))
        """
        return (QuadraticStrata_g,(self._genus,self._min_nb_poles, self._max_nb_poles))

    def __contains__(self, c):
        r"""
        Containance test

        TESTS::

            sage: a = QuadraticStrata(genus=3,max_nb_poles=4)
            sage: all(s in a for s in a)
            True
        """
        if not isinstance(c, QuadraticStratum):
            return False

        return (c.genus() == self._genus and
                c.nb_poles() >= self._min_nb_poles and
                c.nb_poles() <= self._max_nb_poles)

    def _repr_base_(self):
        r"""
        Base string for the representation.

        EXAMPLES::

            sage: QuadraticStrata(genus=4)._repr_base_()
            'Quadratic strata of genus 4 surfaces'
        """
        return "Quadratic strata of genus %d surfaces"%self._genus

    def __iter__(self):
        r"""
        TESTS::

            sage: list(AbelianStrata(genus=1))
            [H_1(0)]
        """
        from itertools import count
        from sage.combinat.partition import Partitions

        g = self._genus
        nb_poles = self._min_nb_poles

        while nb_poles <= self._max_nb_poles:
            s = 4*g-4+nb_poles
            if s == 0:
                nb_poles += 1
                continue
            for p in Partitions(s):
                Q = QuadraticStratum(list(p)+[-1]*nb_poles)
                if not Q.is_empty():
                    yield Q
            nb_poles += 1

    def first(self):
        r"""
        Return the first element of this list of strata.

        EXAMPLES::

            sage: Q = QuadraticStrata(genus=4); Q
            Quadratic strata of genus 4 surfaces
            sage: Q.first()
            Q_4(12)

            sage: Q = QuadraticStrata(genus=3,nb_poles=1); Q
            Quadratic strata of genus 3 surfaces with 1 pole
            sage: Q.first()
            Q_3(9, -1)

            sage: Q = QuadraticStrata(genus=3,min_nb_poles=2); Q
            Quadratic strata of genus 3 surfaces with at least 2 poles
            sage: Q.first()
            Q_3(10, -1^2)
        """
        p = self._min_nb_poles
        g = self._genus
        if g == 1: # Q(0) and Q(1,-1) are empty
            if p == 0 or p == 1:
                raise NotImplementedError, "empty list"
        elif g == 2 and p == 0:  # Q(4) and Q(3,1) are empty
            return QuadraticStratum(2,2)
        return QuadraticStratum([4*g-4+p]+[-1]*p)

    an_element_ = first

    def last(self):
        r"""
        Return the last element of this list of strata.

        EXAMPLES::

            sage: Q = QuadraticStrata(genus=2, nb_poles=0); Q
            Quadratic strata of genus 2 surfaces with no pole
            sage: Q.last()
            Q_2(1^4)

            sage: Q = QuadraticStrata(genus=0); Q
            Quadratic strata of genus 0 surfaces
            sage: Q.last()
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite list

        TESTS::

            sage: Q = QuadraticStrata(genus=1,nb_poles=2)
            sage: Q.list()[-1] == Q.last()
            True

            sage: Q = QuadraticStrata(genus=1,max_nb_poles=2)
            sage: Q.list()[-1] == Q.last()
            True

            sage: Q = QuadraticStrata(genus=0,nb_poles=10)
            sage: Q.list()[-1] == Q.last()
            True
        """
        if self._max_nb_poles == Infinity:
            raise NotImplementedError("infinite list")
        p = self._max_nb_poles
        g = self._genus
        if g == 0:
            if p < 4:
                raise NotImplementedError("empty list")
        if g == 1: # Q(0) and Q(1,-1) are empty
            if p == 0 or p == 1:
                raise NotImplementedError("empty list")
        return QuadraticStratum([1]*(4*g-4+p) + [-1]*p)

class QuadraticStrata_d(QuadraticStrata_class):
    r"""
    Strata with prescribed dimension.

    INPUT:

    - ``dimension`` - an integer greater than 1

    - ``min_nb_poles``, ``max_nb_poles`` - integers (optional) - the min or max
      number of poles

    EXAMPLES::

        sage: for q in QuadraticStrata(dimension=5): print q
        Q_3(8)
        Q_2(2, 1^2)
        Q_2(4, 1, -1)
        Q_2(3, 2, -1)
        Q_2(6, -1^2)
        Q_1(2, 1, -1^3)
        Q_1(4, -1^4)
        Q_0(2, -1^6)

        sage: Q = QuadraticStrata(dimension=6,nb_poles=1); Q
        Quadratic strata of dimension 6 with 1 pole
        sage: for q in Q: print q
        Q_3(9, -1)
        Q_2(3, 1^2, -1)
        Q_2(2^2, 1, -1)
    """
    def __init__(self,dimension,min_nb_poles,max_nb_poles):
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
        self._min_nb_poles = min_nb_poles
        self._max_nb_poles = max_nb_poles

    def __eq__(self, other):
        return (isinstance(other, QuadraticStrata_g) and
                self._dimension == other._dimension and
                self._min_nb_poles == other._min_nb_poles and
                self._max_nb_poles == other._max_nb_poles)

    def __reduce__(self):
        r"""
        Support for pickling.

        EXAMPLES::

            sage: QuadraticStrata(dimension=12).__reduce__()
            (<class 'sage.dynamics.flat_surfaces.quadratic_strata.QuadraticStrata_d'>, (12, 0, +Infinity))
        """
        return (QuadraticStrata_d,(self._dimension,self._min_nb_poles,self._max_nb_poles))

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
        if not isinstance(c, QuadraticStratum):
            return False

        return (c.dimension() == self._dimension and
                c.nb_poles() >= self._min_nb_zeros and
                c.nb_poles() <= self._max_nb_zeros)

    def _repr_base_(self):
        r"""
        Base string for the representation.

        EXAMPLES::

            sage: QuadraticStrata(dimension=12)._repr_base_()
            'Quadratic strata of dimension 12'
        """
        return "Quadratic strata of dimension %d" %(self._dimension)

    def __iter__(self):
        r"""
        Iterator.

        TESTS::

            sage: for q in QuadraticStrata(dimension=6): print q
            Q_3(7, 1)
            Q_3(6, 2)
            Q_3(5, 3)
            Q_3(4^2)
            Q_2(1^4)
            Q_3(9, -1)
            Q_2(3, 1^2, -1)
            Q_2(2^2, 1, -1)
            Q_2(5, 1, -1^2)
            Q_2(4, 2, -1^2)
            Q_2(3^2, -1^2)
            Q_2(7, -1^3)
            Q_1(1^3, -1^3)
            Q_1(3, 1, -1^4)
            Q_1(2^2, -1^4)
            Q_1(5, -1^5)
            Q_0(1^2, -1^6)
            Q_0(3, -1^7)
        """
        from sage.combinat.partition import Partitions

        d = self._dimension
        if d == 2:
            yield QuadraticStratum(-1,-1,-1,-1)
        else:
            m = max(0,self._min_nb_poles)
            M = min(2*d-2,self._max_nb_poles+1)
            for p in xrange(m,M):
                for z in xrange((d+p)%2,min(d+3-p,(2*d-p)//3+1),2):
                    # d+z+p is 0 mod 2
                    # 2d-2z-2p >= -4 (or z <= d+2-p)
                    for Z in Partitions(2*d-2*z-p,length=z):
                        Q = QuadraticStratum(Z+[-1]*p)
                        if not Q.is_empty():
                            yield Q

class QuadraticStrata_gd(QuadraticStrata_class):
    r"""
    Abelian strata with presrcribed genus and dimension.

    INPUT:

    - ``genus`` - an integer - the genus of the surfaces

    - ``dimension`` - an integer - the dimension of strata

    """
    def __init__(self,genus,dimension,min_nb_poles,max_nb_poles):
        r"""
        TESTS::

            sage: s = AbelianStrata(genus=4,dimension=10)
            sage: loads(dumps(s)) == s
            True
        """
        Parent.__init__(self, category=FiniteEnumeratedSets(), name="AbelianStrata_gd")
        self._genus = genus
        self._dimension = dimension
        self._min_nb_poles = min_nb_poles
        self._max_nb_poles = max_nb_poles

    def __eq__(self, other):
        r"""
        Equality test.
        """
        return (isinstance(other, QuadraticStrata_gd) and
                self._genus == other._genus and
                self._dimension == other._dimension and
                self._min_nb_poles == other._min_nb_poles and
                self._max_nb_poles == other._max_nb_poles)

    def __reduce__(self, other):
        return (QuadraticStrata_gd, (self._genus, self._dimension, self._min_nb_poles, self._max_nb_poles))

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
        if not isinstance(c, QuadraticStratum):
            return False

        return (c.dimension() == self._dimension and
                c.genus() == self._genus and
                c.nb_poles() >= c._min_nb_poles and
                c.nb_poles() <= c._max_nb_poles)

    def __reduce__(self):
        r"""
        Support for pickling.

        EXAMPLES::

            sage: QuadraticStrata(genus=4,dimension=10).__reduce__()
            (<class 'sage.dynamics.flat_surfaces.quadratic_strata.QuadraticStrata_gd'>, (4, 10))
        """
        return (QuadraticStrata_gd, (self._genus,self._dimension))

    def _repr_base_(self):
        r"""
        TESTS::

            sage: QuadraticStrata(genus=2,dimension=4)._repr_base_()
            'Quadratic strata of genus 2 surfaces and dimension 4'
        """
        return "Quadratic strata of genus %d surfaces and dimension %d"  %(self._genus, self._dimension)

    def __iter__(self):
        r"""
        Iterator.

        TESTS::

            sage: for a in QuadraticStrata(genus=1, dimension=6): print a
            Q_1(1^3, -1^3)
            Q_1(3, 1, -1^4)
            Q_1(2^2, -1^4)
            Q_1(5, -1^5)
        """
        from sage.combinat.partition import Partitions
        d = self._dimension
        g = self._genus
        if d == 2:
            if g == 0:
                yield QuadraticStratum(-1,-1,-1,-1)
        else:
            m = max(0,self._min_nb_poles)
            M = min(2*d-2,self._max_nb_poles+1)
            for p in xrange(m,M):
                z = d - p - 2*g + 2
                for Z in Partitions(2*d-2*z-p,length=z):
                    Q = QuadraticStratum(Z+[-1]*p)
                    if not Q.is_empty():
                        yield Q

