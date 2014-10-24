r"""
Strata of differential on Riemann surfaces

The templates above are intended to be unifying for both Abelian and quadratic
differentials. Moreover it provides the general structure for the parentship
relation between stratum and its connected components. It also provides the
necessary __reduce__ method for pickling.
"""
#*****************************************************************************
#       Copyright (C) 2009 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent

def list_to_exp_list(l):
    r"""
    Convert list into exponential notation.

    EXAMPLES::

        sage: from sage.dynamics.flat_surfaces.strata import list_to_exp_list
        sage: l = [0,0,2,2,3,2,0,0,0]
        sage: list_to_exp_list(l)
        [(0, 2), (2, 2), (3, 1), (2, 1), (0, 3)]
    """
    d = []
    i = 0
    while i < len(l):
        j = i
        while j < len(l) and l[j] == l[i]:
            j += 1
        d.append((l[i],j-i))
        i = j
    return d

#
# Stratum, stratum component
#

class Stratum(SageObject):
    r"""
    Generic class for stratum of flat surfaces.

    Assumes there are

    - a method .zeros() which returns the list of all zeros

    - a method .nb_zeros() which returns the number of zeros with an option
      fake_zeros which could be true or false

    - a method .nb_fake_zeros() which returns the number of fake zeros

    - a method .dimension() which returns the dimension of the stratum

    - an attribute ._cc which is a list of classes associated to the
      connected components of self

    There may be

    - an attribute ._name which corresponds to the begining of the string
      representation (default is the empty string)

    - an attribute ._latex_name which corresponds to the begining of the latex
      string representation (uses _name by default)

    """
    _name = ''
    _latex_name = ''

    #
    # String representation
    #

    def _flat_zero_str(self):
        r"""
        String representation of the zeros.

        EXAMPLES::

            sage: a = AbelianStratum({2:3})
            sage: a._flat_zero_str()
            '2, 2, 2'
        """
        return ', '.join(map(str,self.zeros()))

    def _exp_zero_str(self):
        r"""
        String representation with exponential notation

        EXAMPLES::

            sage: a = AbelianStratum(2,2,2)
            sage: a._exp_zero_str()
            '2^3'
        """
        return ', '.join('%d^%d' %(i,e) if e != 1 else '%d' %i for (i,e) in list_to_exp_list(self.zeros()))

    # this attribute can be switched between _flat_zero_str and _exp_zero_str
    _zero_str = _exp_zero_str

    def _repr_(self):
        """
        TESTS::

            sage: repr(AbelianStratum(1,1))   #indirect doctest
            'H_2(1^2)'
        """
        return self._name + "_" + str(self.genus()) + "(" + self._zero_str() + ")"

    def _latex_(self):
        r"""
        Latex string representation

        EXAMPLES::

            sage: AbelianStratum(0)._latex_()
            '\\mathcal{H}_1(0)'
            sage: QuadraticStratum({-1:4})._latex_()
            '\\mathcal{Q}_0(-1^4)'
        """
        return self._latex_name + '_' + str(self.genus()) + "(" + self._zero_str() + ")"

    #
    # Equality and comparisons
    #

    def __hash__(self):
        r"""
        Hash value for self

        EXAMPLES::

            sage: hash(AbelianStratum(0)) == hash(AbelianStratum(0))
            True
            sage: hash(AbelianStratum(2,2)) == hash(AbelianStratum(1,1))
            False
            sage: s = AbelianStratum(1,1,1,1)
            sage: s == loads(dumps(s))
            True

            sage: AbelianStratum('no','way')
            Traceback (most recent call last):
            ...
            ValueError: input must be a list of integers

            sage: AbelianStratum([1,1,1,1], marked_separatrix='full')
            Traceback (most recent call last):
            ...
            ValueError: marked_separatrix must be one of 'no', 'in', 'out'
        """
        if l == ():
            pass

        elif hasattr(l[0], "__iter__") and len(l) == 1:
            l = l[0]

        if not all(isinstance(i, (Integer, int)) for i in l):
            raise ValueError("input must be a list of integers")

        if 'marked_separatrix' in d:
            m = d['marked_separatrix']

            if m is None:
                m = 'no'

            if (m != 'no' and m != 'in' and m != 'out'):
                raise ValueError("marked_separatrix must be one of 'no', "
                                 "'in', 'out'")
            self._marked_separatrix = m

        else:  # default value
            self._marked_separatrix = 'no'

        self._zeroes = list(l)

        if not self._marked_separatrix is 'no':
            self._zeroes[1:] = sorted(self._zeroes[1:], reverse=True)
        else:
            self._zeroes.sort(reverse=True)

        self._genus = sum(l)/2 + 1

        self._genus = Integer(self._genus)

        zeroes = sorted(x for x in self._zeroes if x > 0)

        if self._genus == 1:
            self._cc = (HypCCA,)

        elif self._genus == 2:
            self._cc = (HypCCA,)

        elif self._genus == 3:
            if zeroes == [2, 2] or zeroes == [4]:
                self._cc = (HypCCA, OddCCA)
            else:
                self._cc = (CCA,)

        elif len(zeroes) == 1:
            # just one zeros [2g-2]
            self._cc = (HypCCA, OddCCA, EvenCCA)

        elif zeroes == [self._genus-1, self._genus-1]:
            # two similar zeros [g-1, g-1]
            if self._genus % 2 == 0:
                self._cc = (HypCCA, NonHypCCA)

            else:
                self._cc = (HypCCA, OddCCA, EvenCCA)

        elif len([x for x in zeroes if x % 2]) == 0:
            # even zeroes [2 l_1, 2 l_2, ..., 2 l_n]
            self._cc = (OddCCA, EvenCCA)

        else:
            self._cc = (CCA, )

    def _repr_(self):
        """
        TESTS::

            sage: repr(AbelianStratum(1,1))   #indirect doctest
            'H(1, 1)'
        """
        if self._marked_separatrix == 'no':
            return "H(" + str(self._zeroes)[1:-1] + ")"
        else:
            return ("H" +
                    '^' + self._marked_separatrix +
                    "(" + str(self._zeroes)[1:-1] + ")")

    def __str__(self):
        r"""
        TESTS::

            sage: str(AbelianStratum(1,1))
            'H(1, 1)'
>>>>>>> FETCH_HEAD
        """
        return hash(self._name) + hash(tuple(self.zeros()))

    def __eq__(self, other):
        r"""
        Equality test

        EXAMPLES::

            sage: AbelianStratum(0) == AbelianStratum(0)
            True
            sage: QuadraticStratum(5,-1) == QuadraticStratum(5,-1)
            True
            sage: AbelianStratum(12) == QuadraticStratum(12)
            False
            sage: QuadraticStratum(12,0,0) == QuadraticStratum(12,0)
            False
            sage: AbelianStratum(2,0) == AbelianStratum(2)
            False
        """
        if not isinstance(other, Stratum):
            return False
        return type(self) == type(other) and self.zeros() == other.zeros()

    def __cmp__(self, other):
        r"""
        Comparison

        ALGORITHM:

        First compare the class, then the dimension and then the list of zeros.

        TESTS::

            sage: AbelianStratum(1,1) < AbelianStratum(1,1,0)
            True
            sage: AbelianStratum(1,1,0) < AbelianStratum(1,1)
            False
            sage: AbelianStratum(1,1,0) < AbelianStratum(1,1,0,0)
            True
            sage: AbelianStratum(2) < AbelianStratum(1,1)
            True
            sage: AbelianStratum(4,0) > AbelianStratum(1,1,1,1)
            False
            sage: AbelianStratum(4,0,0,0) > AbelianStratum(1,1,1,1)
            True

        ::

            sage: QuadraticStratum(2,2) < QuadraticStratum(2,2,0)
            True
            sage: QuadraticStratum(2,2,0) < QuadraticStratum(2,2)
            False
            sage: QuadraticStratum(2,2,0) < QuadraticStratum(2,2,0,0)
            True
            sage: QuadraticStratum(4) < QuadraticStratum(2,2)
            True
            sage: QuadraticStratum(4,0) > QuadraticStratum(1,1,1,1)
            False
            sage: QuadraticStratum(4,0,0,0) > QuadraticStratum(1,1,1,1)
            True
        """
        if not isinstance(other,Stratum):
            return NotImplemented

        # compare the type of stratum
        test = cmp(type(self), type(other))
        if test: return test

        # compare the dimension
        test = cmp(self.dimension(), other.dimension())
        if test: return test

        # compare the list of zeros
        test = cmp(self.zeros(),other.zeros())
        if test: return test

        # they are equal
        return 0

    #
    # Connected components
    #

    def is_connected(self):
        r"""
        Test if the strata is connected.

        EXAMPLES:

        ::

            sage: AbelianStratum([2]).is_connected()
            True
            sage: AbelianStratum([2,2]).is_connected()
            False
            sage: QuadraticStratum([-1,-1,-1,-1]).is_connected()
            True
            sage: QuadraticStratum([12]).is_connected()
            False
        """
        return len(self._cc) <= 1

    def permutation_representative(self, *args, **kwds):
        r"""
        Return a permutation of interval exchanges associated to this stratum.

        EXAMPLES:

        Examples from Abelian differentials::

            sage: a = AbelianStratum([3,2,1,0,0])
            sage: p = a.permutation_representative()
            sage: p.stratum()
            H_4(3, 2, 1, 0^2)
            sage: a = AbelianStratum([2, 2, 2])
            sage: p = a.permutation_representative()
            sage: p.stratum()
            H_4(2^3)

        Examples from quadratic differentials::

            sage: a = QuadraticStratum([6,-1,-1])
            sage: p = a.permutation_representative()
            sage: p.stratum()
            Q_2(6, -1^2)
            sage: a = QuadraticStratum([-1,-1,-1,-1,0,0])
            sage: p = a.permutation_representative()
            sage: p.stratum()
            Q_0(0^2, -1^4)
        """
        return self.one_component().permutation_representative(*args,**kwds)

    def is_empty(self):
        r"""
        Return True if the stratum is empty

        EXAMPLES::

            sage: AbelianStratum(2).is_empty()
            False
            sage: QuadraticStratum(1,-1).is_empty()
            True
        """
        return len(self._cc) == 0

    def number_of_components(self):
        r"""
        Returns the number of connected components of self

        EXAMPLES::

            sage: AbelianStratum(2).number_of_components()
            1
            sage: AbelianStratum(4).number_of_components()
            2
            sage: AbelianStratum(3,3).number_of_components()
            2
        """
        return len(self._cc)

    def one_component(self):
        r"""
        Returns a connected component of this stratum.

        EXAMPLES::

            sage: AbelianStratum(2).one_component()
            H_2(2)^hyp
        """
        if self.components():
            return self.components()[-1]
        from sage.categories.sets_cat import EmptySetError
        raise EmptySetError, "The stratum is empty"

    def random_component(self):
        r"""
        Returns a random connected component of this stratum.

        EXAMPLES::

            sage: Q = QuadraticStratum(6,6)
            sage: Q.random_component()
            Q_4(6^2)^hyp
            sage: Q.random_component()
            Q_4(6^2)^reg
        """
        if self.components():
            from sage.misc.prandom import choice
            return choice(self.components())
        from sage.categories.sets_cat import EmptySetError
        raise EmptySetError, "The stratum is empty"

    def components(self):
        """
        Lists the connected components of the Stratum.

        OUTPUT:

        list -- a list of connected components of stratum

        EXAMPLES:

        ::

            sage: AbelianStratum(0).components()
            [H_1(0)^hyp]

        ::

            sage: AbelianStratum(2).components()
            [H_2(2)^hyp]
            sage: AbelianStratum(1,1).components()
            [H_2(1^2)^hyp]

        ::

            sage: AbelianStratum(4).components()
            [H_3(4)^hyp, H_3(4)^odd]
            sage: AbelianStratum(3,1).components()
            [H_3(3, 1)^c]
            sage: AbelianStratum(2,2).components()
            [H_3(2^2)^hyp, H_3(2^2)^odd]
            sage: AbelianStratum(2,1,1).components()
            [H_3(2, 1^2)^c]
            sage: AbelianStratum(1,1,1,1).components()
            [H_3(1^4)^c]
        """
        return map(lambda x: x(self), self._cc)

class StratumComponent(SageObject):
    r"""
    Generic class for connected component of a stratum of flat surfaces.

    Assumes there are implemented

    - a method .permutation_representative()

    There may be

    - an attribute ._name

    - an attribute ._latex_name

    """
    _name = ''
    _latex_name = ''

    def __init__(self, stratum):
        r"""
        TEST::

            sage: a = AbelianStratum(4,4).one_component()
            sage: a == loads(dumps(a))
            True
            sage: q = QuadraticStratum(5,5,-1,-1).one_component()
            sage: q == loads(dumps(q))
            True
        """
        self._stratum = stratum

    def __reduce__(self):
        r"""
        Reduce method for pickling

        TESTS:

        Tests for Abelian strata::

            sage: a = AbelianStratum(2,2)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True
            sage: a = AbelianStratum(3,3)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True
            sage: a = AbelianStratum(6)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True
            sage: a = AbelianStratum(1,1,1,1)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True

        Tests for quadratic strata::

            sage: q = QuadraticStratum(-1,-1,-1,-1)
            sage: all(loads(dumps(cc)) == cc for cc in q.components())
            True
            sage: q = QuadraticStratum(12)
            sage: all(loads(dumps(cc)) == cc for cc in q.components())
            True
        """
        return (self.__class__, (self._stratum,))

    def _repr_(self):
        r"""
        String representation

        EXAMPLES::

            sage: a_hyp = AbelianStratum(4).hyperelliptic_component()
            sage: a_hyp._repr_()
            'H_3(4)^hyp'
            sage: a_odd = AbelianStratum(4).odd_component()
            sage: a_odd._repr_()
            'H_3(4)^odd'
        """
        g = self._parent._genus
        zeroes = [x for x in self._parent._zeroes if x > 0]
        n = self._parent._zeroes.count(0)

        l0 = range(0, 4*g-3)
        l1 = [4, 3, 2]
        for k in range(5, 4*g-6, 4):
            l1 += [k, k+3, k+2, k+1]
        l1 += [1, 0]
        k = 3
        for d in zeroes:
            for i in range(d-1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 2
            k += 2

        if n != 0:
            interval = range(4*g-3, 4*g-3+n)

            if self._parent._zeroes[0] == 0:
                k = l0.index(4)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if self._parent._marked_separatrix == 'in':
            l0.reverse()
            l1.reverse()

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            return ReducedPermutationIET([l0, l1], alphabet=alphabet)

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            return LabelledPermutationIET([l0, l1], alphabet=alphabet)

    def stratum(self):
        r"""
        Return the stratum associated to self

        EXAMPLES::

            sage: a = AbelianStratum(4,4)
            sage: all([c.stratum() == a for c in a.components()])
            True
        """
        return self._stratum

    def __eq__(self,other):
        r"""
        Equality test

        EXAMPLES::

            sage: c_hyp = AbelianStratum(6).hyperelliptic_component()
            sage: c_odd = AbelianStratum(6).odd_component()
            sage: c_hyp == c_hyp
            True
            sage: c_hyp == c_odd
            False
        """
        if not isinstance(other, StratumComponent):
            return NotImplemented

        return (type(self) == type(other) and self._stratum == other._stratum)

    def __cmp__(self, other):
        r"""
        Comparison

        TESTS::

            sage: a1 = AbelianStratum(1,1,1,1)
            sage: c1 = a1.components()[0]
            sage: a2 = AbelianStratum(3,1)
            sage: c2 = a2.components()[0]
            sage: c1 == c1
            True
            sage: c1 == c2
            False
            sage: a1 = AbelianStratum(1,1,1,1)
            sage: c1 = a1.components()[0]
            sage: a2 = AbelianStratum(2, 2)
            sage: c2_hyp, c2_odd = a2.components()
            sage: c1 != c1
            False
            sage: c1 != c2_hyp
            True
            sage: c2_hyp != c2_odd
            True
            sage: c1 == True
            Traceback (most recent call last):
            ...
            TypeError: other must be a connected component
        """
        if not isinstance(other, CCA):
            raise TypeError("other must be a connected component")

        if isinstance(self, type(other)):
            if self._parent._zeroes > other._parent._zeroes:
                return 1
            elif self._parent._zeroes < other._parent._zeroes:
                return -1
            return 0

        return cmp(type(self), type(other))

CCA = ConnectedComponentOfAbelianStratum


class HypConnectedComponentOfAbelianStratum(CCA):
    """
    Hyperelliptic component of Abelian stratum.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'hyp'

    def representative(self, reduced=True, alphabet=None):
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

        EXAMPLES:

        ::

            sage: c = AbelianStratum(0).connected_components()[0]
            sage: c
            H_hyp(0)
            sage: p = c.representative(alphabet="01")
            sage: p
            0 1
            1 0
            sage: p.connected_component()
            H_hyp(0)

        ::

            sage: c = AbelianStratum(0,0).connected_components()[0]
            sage: c
            H_hyp(0, 0)
            sage: p = c.representative(alphabet="abc")
            sage: p
            a b c
            c b a
            sage: p.connected_component()
            H_hyp(0, 0)

        ::

            sage: c = AbelianStratum(2).connected_components()[0]
            sage: c
            H_hyp(2)
            sage: p = c.representative(alphabet="ABCD")
            sage: p
            A B C D
            D C B A
            sage: p.connected_component()
            H_hyp(2)

        ::

            sage: c = AbelianStratum(1,1).connected_components()[0]
            sage: c
            H_hyp(1, 1)
            sage: p = c.representative(alphabet="01234")
            sage: p
            0 1 2 3 4
            4 3 2 1 0
            sage: p.connected_component()
            H_hyp(1, 1)
        """
        g = self._parent._genus
        n = self._parent._zeroes.count(0)
        m = len(self._parent._zeroes) - n

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

            if self._parent._zeroes[0] == 0:
                l0[-1:-1] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1[1:1] = interval

        else:  # H(g-1,g-1,0^n) or H(0,g-1,g-1,0^(n-1))
            l0 = range(1, 2*g+2)
            l1 = range(2*g+1, 0, -1)
            interval = range(2*g+2, 2*g+n+2)

            if self._parent._zeroes[0] == 0:
                l0[-1:-1] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1[1:1] = interval

        if self._parent._marked_separatrix == 'in':
            l0.reverse()
            l1.reverse()

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            return ReducedPermutationIET([l0, l1], alphabet=alphabet)

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            return LabelledPermutationIET([l0, l1], alphabet=alphabet)

HypCCA = HypConnectedComponentOfAbelianStratum


class NonHypConnectedComponentOfAbelianStratum(CCA):
    """
    Non hyperelliptic component of Abelian stratum.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'nonhyp'

NonHypCCA = NonHypConnectedComponentOfAbelianStratum


class EvenConnectedComponentOfAbelianStratum(CCA):
    """
    Connected component of Abelian stratum with even spin structure.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'even'

    def representative(self, reduced=True, alphabet=None):
        r"""
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitely interval exchange
        transformations for each stratum in [Zor08]_.

        EXAMPLES:

        ::

            sage: c = AbelianStratum(6).connected_components()[2]
            sage: c
            H_even(6)
            sage: p = c.representative(alphabet=range(8))
            sage: p
            0 1 2 3 4 5 6 7
            5 4 3 2 7 6 1 0
            sage: p.connected_component()
            H_even(6)

        ::

            sage: c = AbelianStratum(4,4).connected_components()[2]
            sage: c
            H_even(4, 4)
            sage: p = c.representative(alphabet=range(11))
            sage: p
            0 1 2 3 4 5 6 7 8 9 10
            5 4 3 2 6 8 7 10 9 1 0
            sage: p.connected_component()
            H_even(4, 4)
        """
        zeroes = [x for x in self._parent._zeroes if x > 0]
        n = self._parent._zeroes.count(0)
        g = self._parent._genus

        l0 = range(3*g-2)
        l1 = [6, 5, 4, 3, 2, 7, 9, 8]
        for k in range(10, 3*g-4, 3):
            l1 += [k, k+2, k+1]
        l1 += [1, 0]

        k = 4
        for d in zeroes:
            for i in range(d/2-1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 3
            k += 3

        # if there are marked points we transform 0 in [3g-2, 3g-3, ...]
        if n != 0:
            interval = range(3*g-2, 3*g - 2 + n)

            if self._parent._zeroes[0] == 0:
                k = l0.index(6)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if self._parent._marked_separatrix == 'in':
            l0.reverse()
            l1.reverse()

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            return ReducedPermutationIET([l0, l1], alphabet=alphabet)

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            return LabelledPermutationIET([l0, l1], alphabet=alphabet)

EvenCCA = EvenConnectedComponentOfAbelianStratum


class OddConnectedComponentOfAbelianStratum(CCA):
    r"""
    Connected component of an Abelian stratum with odd spin parity.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'odd'

    def representative(self, reduced=True, alphabet=None):
        """
        if not isinstance(other, StratumComponent):
            return NotImplemented

        test = cmp(self._stratum, other._stratum)

        if test: return test

        return cmp(type(self),type(other))

    def genus(self):
        return self._stratum.genus()

#
# Strata (family of strata)
#

           sage: a = AbelianStratum(4,4).connected_components()[1]
            sage: print a.representative(alphabet=range(11))
            0 1 2 3 4 5 6 7 8 9 10
            3 2 5 4 6 8 7 10 9 1 0
        """
        zeroes = [x//2 for x in self._parent._zeroes if x > 0]

        n = self._parent._zeroes.count(0)
        g = self._parent._genus

        l0 = range(3*g-2)
        l1 = [3, 2]
        for k in range(4, 3*g-4, 3):
            l1 += [k, k+2, k+1]
        l1 += [1, 0]

        k = 4
        for d in zeroes:
            for i in range(d-1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 3
            k += 3

        # marked points
        if n != 0:
            interval = range(3*g-2, 3*g-2+n)

            if self._parent._zeroes[0] == 0:
                k = l0.index(3)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if self._parent._marked_separatrix == 'in':
            l0.reverse()
            l1.reverse()

        if reduced:
            from sage.dynamics.interval_exchanges.reduced import ReducedPermutationIET
            return ReducedPermutationIET([l0, l1], alphabet=alphabet)

        else:
            from sage.dynamics.interval_exchanges.labelled import LabelledPermutationIET
            return LabelledPermutationIET([l0, l1], alphabet=alphabet)

OddCCA = OddConnectedComponentOfAbelianStratum
