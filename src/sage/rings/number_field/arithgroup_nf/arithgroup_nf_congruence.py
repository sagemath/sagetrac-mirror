r"""
Congruence subgroups of `{\rm SL}_2(\mathcal{O})` and `{\rm GL}_2(\mathcal{O})`,
where `\mathcal{O}` is a (maximal) order in a number field `K`.

AUTHORS:

- Fredrik Stromberg (2013): initial version based on arithgroup_generic.py
- Marco Streng (2014): more examples in the documentation
"""

################################################################################
#
#       Copyright (C) 2009, The Sage Group -- http://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
#
################################################################################

from sage.rings.number_field.arithgroup_nf.arithgroup_nf_generic import ArithmeticSubgroup_NF_class
from sage.rings.integer import Integer
from sage.all import matrix


class HilbertModularGroup_CongruenceSubgroup_class(ArithmeticSubgroup_NF_class):
    r"""
    Base class for congruence subgroups of `{\rm SL}_{2}(O)`.
    """
    def __init__(self, order, group='SL', name='', ltx=''):
        r"""
        Standard init routine.

        INPUT:

        - `order` -- maximal order in a number field
        - `group` -- string ('SL', 'GL' or 'GL+')
        - `name` -- string.
        - `ltx` -- string

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_class(O,'SL')
            sage: H

        """
        super(HilbertModularGroup_CongruenceSubgroup_class, self).__init__(order, group, name, ltx)
        self._ncusps = None

    def level(self):
        r"""
        Return the level (an ideal) of self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_class(O,'SL')
            sage: H.level()
        """
        if not hasattr(self, '_level'):
            raise NotImplementedError("Level not implemented for %s" % self)
        return self._level

    def is_congruence(self):
        r"""
        Check if self is congruence.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_class(O,'SL')
            sage: H.is_congruence()
            True

        """
        return True

#    def __action__()

    def index(self, G=None):
        r"""
        Return the index of self in SL(2,O).

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_class(O,'SL')
            sage: H.index()
            

        """
        if self.level() == self.base_ring().ideal(1):
            return 1
        raise NotImplementedError("Index not implemented for %s" % self)

    def coset_reps(self, G=None):
        r"""
        Return right coset representatives of self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_class(O,'SL')
            sage: H.coset_reps()
            
        """
        if self.index() != 1 or G is not None:
            raise NotImplementedError
        return matrix(self.number_field(), 2, 2, [1, 0, 0, 1])

    def gens(self):
        r"""
        Return generators of self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_class(O,'SL')
            sage: H.gens()

        """
        if self.index() != 1:
            raise NotImplementedError("Generators not implemented for %s" %
                                      self)
        if self.group() != 'SL':
            raise NotImplementedError("Generators not implemented for %s" %
                                      self)
        gens = [matrix(self.number_field(), 2, 2, [0, 1, -1, 0])]
        for eps in self.base_ring().basis():
            gens.append(matrix(self.number_field(), 2, 2, [1, eps, 0, 1]))
            gens.append(matrix(self.number_field(), 2, 2, [eps, 0,
                                                           0, eps ** -1]))
        return gens

    def ncusps(self):
        r"""
        Return the number of inequivalent cusps of self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_class(O,'SL')
            sage: H.ncusps()

        """
        if self._ncusps is not None:
            return self._ncusps
        if self.index() == 1:
            return self.number_field().class_number()
        raise NotImplementedError("Number of cusps not implemented for %s" %
                                  self)

    def _find_cusps(self):
        r"""
        Find a set of cusp representatives.

        .. NOTE::
        
            Currently only implemented for the full Hilbert modular
            group SL(2,O)

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_class(O,'SL')
            sage: H._find_cusps()

        """
        if self.level() != self.base_ring().ideal(1) or not self.is_special():
            raise NotImplementedError
        K = self.base_field()
        if self.base_ring() != self.base_field().ring_of_integers():
            # Then we have representatives given by ideal classes
            lreps = map(lambda x: x.ideal(), self._class_group.list())
            ncusps = len(lreps)
            for a in lreps:
                if self._verbose > 0:
                    print "Set cusp info for a={0}".format(a)
                if a.is_trivial():
                    ca = NFCusp(K(1), K(0), lreps=lreps)
                else:
                    ag = a.gens()
                    ca = NFCusp(K, ag[0], ag[1], lreps=lreps)
                cusps.append(ca)
        else:
            raise NotImplementedError
        return cusps


class HilbertModularGroup_CongruenceSubgroup_Gamma_class(HilbertModularGroup_CongruenceSubgroup_class):
    r"""
    Principal congruence subgroup `\Gamma(N)` of `{\rm SL}_{2}(O)`.
    where `O` is a maximal order in a number field and `N` an ideal in `O`.
    """
    def __init__(self, order, level, group, name='', ltx=''):
        r"""
        Initialize a principal congruence subgroup of SL(2,O).

        INPUT:

        - `order` -- maximal order in a number field
        - `level` -- ideal in order
        - `group` -- string ('SL', 'GL' or 'GL+')
        - `name` -- string.
        - `ltx` -- string

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma_class(O, 3*O, 'SL')
            sage: H

        """
        assert order.ideal(level) == level
        self._level = level
        super(HilbertModularGroup_CongruenceSubgroup_Gamma_class,
              self).__init__(order, group, name, ltx)

    def __contains__(self, A):
        r"""
        Check if A is in self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma_class(O, 3*O, 'SL')
            sage: Matrix([[1,3],[0,1]]) in H  # indirect
            True
            sage: Matrix([[1,2],[0,1]]) in H  # indirect
            False
            sage: Matrix([[2,3],[3,5]]) in H  # indirect
            False
            sage: Matrix([[1,3],[3,1]]) in H  # indirect
            False

        """
        if A.parent() == self:
            return True
        try:
            a, b, c, d = self(A).matrix().list()   # ?
        except TypeError:
            return False
        if a - 1 not in self.level() or d - 1 not in self.level():
            return False
        if b not in self.level() or c not in self.level():
            return False
        return True

    def __eq__(self, other):
        r"""
        Check if other is equal to self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H1 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma_class(O, 3*O, 'SL')
            sage: H2 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma_class(O, 1*O, 'SL')
            sage: H3 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma_class(O, -3*O, 'SL')
            sage: H1 == H2  # indirect
            False
            sage: H1 == H3  # indirect
            True
            sage: H1 is H3  # indirect
            True

        """
        if not isinstance(other, HilbertModularGroup_CongruenceSubgroup_Gamma_class):
            return False
        return self.level() == other.level()

    def is_subgroup(self, other):
        r"""
        Check if self is a subgroup of other.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: Gamma6 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma_class(O, 6*O, 'SL')
            sage: Gamma2 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma_class(O, 2*O, 'SL')
            sage: Gamma3 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma_class(O, 3*O, 'SL')
            sage: Gamma6.is_subgroup(Gamma2)
            True
            sage: Gamma2.is_subgroup(Gamma3)
            False
            sage: Gamma3.is_subgroup(Gamma6)
            False
            sage: Gamma6.is_subgroup(Gamma3)
            True
            sage: Gamma3.is_subgroup(Gamma2)
            False
            sage: Gamma2.is_subgroup(Gamma6)
            False
            
            sage: Gamma06 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma0_class(O, 6*O, 'SL')
            sage: Gamma06.is_subgroup(Gamma6)
            False
            sage: Gamma6.is_subgroup(Gamma06)
            True
            
            
            
        """
        if isinstance(other, HilbertModularGroup_CongruenceSubgroup_Gamma_class):
            if other.level().divides(self.level()):
                return True
            else:
                return False
        return super(HilbertModularGroup_CongruenceSubgroup_Gamma_class, self).is_subgroup(other)


class HilbertModularGroup_CongruenceSubgroup_Gamma0_class(HilbertModularGroup_CongruenceSubgroup_class):
    r"""
    Congruence subgroup `\Gamma_{0}(N)` of `{\rm SL}_{2}(O)`,
    where `O` is a maximal order in a number field and `N` an ideal in `O`.

    """
    def __init__(self, order, level=None, group='SL', name='', ltx=''):
        r"""
        Congruence subgroup `\Gamma_{0}(N)` of `{\rm SL}_{2}(O)`.

        INPUT:

        - `order` -- maximal order in a number field
        - `level` -- ideal in order
        - `group` -- string ('SL', 'GL' or 'GL+')
        - `name` -- string.
        - `ltx` -- string

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma0_class(O, 3*O, 'SL')
            sage: H
            
        """
        if level is None:
            self._level = order
        else:
            self._level = level
        super(HilbertModularGroup_CongruenceSubgroup_Gamma0_class, self).__init__(order, group, name, ltx)

    def __contains__(self, A):
        r"""
        Check if A is in self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma0_class(O, 3*O, 'SL')
            sage: Matrix([[1,3],[0,1]]) in H  # indirect
            True
            sage: Matrix([[1,2],[0,1]]) in H  # indirect
            False
            sage: Matrix([[2,3],[3,5]]) in H  # indirect
            True
            sage: Matrix([[1,3],[3,1]]) in H  # indirect
            False


        """
        try:
            a, b, c, d = A
        except:
            return False
        if a not in self.base_ring() or b not in self.base_ring() or d not in self.base_ring():
            return False
        if c not in self.level():
            return False
        return True

    def __eq__(self, other):
        r"""
        Check if other is equal to self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: Gamma01 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma0_class(O, 1*O, 'SL')
            sage: Gamma02 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma0_class(O, 2*O, 'SL')
            sage: Gamma01 == Gamma02 # indirect
            False
            
            sage: Gamma1 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma_class(O, 1*O, 'SL')
            sage: Gamma1 == Gamma01  # indirect
            True
            sage: Gamma1 == Gamma02  # indirect
            False

        """
        if not isinstance(other, HilbertModularGroup_CongruenceSubgroup_Gamma0_class):
            return False
        return self.level() == other.level()

    def is_subgroup(self, other):
        r"""
        Check if self is a subgroup of other.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: Gamma06 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma0_class(O, 6*O, 'SL')
            sage: Gamma02 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma0_class(O, 2*O, 'SL')
            sage: Gamma03 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma0_class(O, 3*O, 'SL')
            sage: Gamma06.is_subgroup(Gamma02)
            True
            sage: Gamma02.is_subgroup(Gamma03)
            False
            sage: Gamma03.is_subgroup(Gamma06)
            False
            sage: Gamma06.is_subgroup(Gamma03)
            True
            sage: Gamma03.is_subgroup(Gamma02)
            False
            sage: Gamma02.is_subgroup(Gamma06)
            False
            
            sage: Gamma6 = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma_class(O, 6*O, 'SL')
            sage: Gamma06.is_subgroup(Gamma6)
            False
            sage: Gamma6.is_subgroup(Gamma06)
            True

        """
        if isinstance(other, HilbertModularGroup_CongruenceSubgroup_Gamma0_class):
            if other.level().divides(self.level()):
                return True
            else:
                return False
        return super(HilbertModularGroup_CongruenceSubgroup_Gamma0_class, self).is_subgroup(other)

    def fundamental_domain(self, **kwds):
        r"""
        Return a fundamental domain for self. Not yet implemented.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.arithgroup_nf_congruence.HilbertModularGroup_CongruenceSubgroup_Gamma0_class(O, 3*O, 'SL')
            sage: H.fundamental_domain()
            Traceback (most recent call last):
            ...
            NotImplementedError:

        """
        if self._fundamental_domain is not None:
            return self._fundamental_domain
        #if self.index()<>1:
        raise NotImplementedError
        #self._fundamental_domain =  HilbertFundamentalDomain(self,*kwds)
        #return self._fundamental_domain


class HilbertModularGroup_CongruenceSubgroup_Gamma00_class(
        HilbertModularGroup_CongruenceSubgroup_Gamma0_class):
    r"""
    Base class for congruence subgroups Gamma_{0}^{0} of SL(2,O)
    i.e of matrices [[a,b],[c,d]] with b and c in given ideals.
    """
    def __init__(self, order, b_ideal, c_ideal, group='SL', name='', ltx=''):
        r"""
        Init a congruence group Gamma_0^0

        EXAMPLES::

            sage: "example needed"
            "example needed"

        """
        self._b_ideal = b_ideal
        self._c_ideal = c_ideal
        self._level = b_ideal * c_ideal
        super(CongruenceSubgroup_Gamma00_class, self).__init__(order, c_ideal,
                                                               group, name,
                                                               ltx)

    def b_ideal(self):
        r"""
        Return the ideal bi such that self contains matrices
        congruent to [[*,bi],[ci,*]]

        EXAMPLES::

            sage: "example needed"
            "example needed"

        """
        return self._b_ideal

    def c_ideal(self):
        r"""
        Return the ideal ci such that self contains matrices
        congruent to [[*,bi],[ci,*]]

        EXAMPLES::

            sage: "example needed"
            "example needed"

        """
        return self._c_ideal

    def __contains__(self, A):
        r"""
        Check if A is in self.

        EXAMPLES::

            sage: "example needed"
            "example needed"

        """
        try:
            a, b, c, d = A
        except:
            return False
        if a not in self.base_ring() or d not in self.base_ring():
            return False
        if b not in self.b_ideal():
            return False
        if c not in self.c_ideal():
            return False
        return True

    def __eq__(self, other):
        r"""
        Check if other is equal to self.

        EXAMPLES::

            sage: "example needed"
            "example needed"

        """
        if not isinstance(other, HilbertModularGroup_CongruenceSubgroup_Gamma00_class):
            return False
        return self.level() == other.level()


class HilbertModularGroup_Conjugate_class(ArithmeticSubgroup_NF_class):
    r"""
    Base class for arithmetic subgroups of `SL(2,K)' of the form
    `{\rm SL}_{2}(O+a)` where `O` is a maximal order in the number field `K`
    and `a` is an ideal in `O`.
    """
    def __init__(self, order, a, group='SL', name='', ltx=''):
        r"""
        Init  a conjugate of a Hilbert modular group.

        EXAMPLES::

            sage: "example needed"
            "example needed"

        """
        super(HilbertModularGroup_Conjugate_class, self).__init__(order, group, name, ltx)
        self._a = a

    def plusa(self):
        r"""
        Return the ideal a such that self = SL_2(O+a).

        EXAMPLES::

            sage: "example needed"
            "example needed"

        """
        return self._a

    def _find_cusps(self):
        r"""
        Find a set of cusp representatives.

        .. NOTE::
        
            Currently only implemented for the full Hilbert modular
            group SL(2,O)

        EXAMPLES::

            sage: "example needed"
            "example needed"

        """
        if self.plusa() != self.base_ring().ideal(1) or not self.is_special():
            raise NotImplementedError
        K = self.base_field()
        if self.base_ring() != self.base_field().ring_of_integers():
            # Then we have representatives given by ideal classes
            lreps = map(lambda x: x.ideal(), self._class_group.list())
            ncusps = len(lreps)
            for a in lreps:
                if self._verbose > 0:
                    print "Set cusp info for a={0}".format(a)
                if a.is_trivial():
                    ca = NFCusp(K(1), K(0), lreps=lreps)
                else:
                    ag = a.gens()
                    ca = NFCusp(K, ag[0], ag[1], lreps=lreps)
                cusps.append(ca)
        else:
            raise NotImplementedError
        return cusps


def HilbertModularGroup(O, group='SL', a=None, **kwds):
    r"""
    Returns the Hilbert modular group `SL_{2}(O\oplus \frak{a})`

    Here `O` is a maximal order in a number field and `\frak{a}` an
    ideal in `O`, consisting of matrices of the form `[a b // c d ]`
    with `a,d` in `O`, `c \in \frak{a}` and `d \in \frak{a}^{-1}`.

    INPUT:

    - 'O' -- Order in number field, *or*
          -- Integer, in which case we set O to be the ring of integers in Q[\sqrt(D))]
    - 'a' -- Ideal in O. Default = None (meaning it is O itself)

    EXAMPLES::

        sage: from sage.rings.number_field.arithgroup_nf.all import *
        sage: K = QuadraticField(41)
        sage: O = K.ring_of_integers()
        sage: G = HilbertModularGroup(O); G
        Hilbert modular group `SL_{2}(O)`
    """
    if isinstance(O, (Integer, int)):
        from sage.rings.number_field.number_field import QuadraticField
        O = QuadraticField(O).ring_of_integers()
    if a is None:
        name = 'Hilbert modular group `{group}_{{2}}(O)`'.format(group=group)
        ltx = 'Hilbert modular group `{group}_{{2}}(O)`'.format(group=group)
        return HilbertModularGroup_CongruenceSubgroup_Gamma0_class(O, O.ideal(1), group, name, ltx)
    name = 'Hilbert modular group `{group}_{{2}}(O+a)`'.format(group=group)
    ltx = '{group}_{{2}}(O+a)'.format(group=group)
    return HilbertModularGroup_Conjugate_class(O, a, group, name, ltx)
