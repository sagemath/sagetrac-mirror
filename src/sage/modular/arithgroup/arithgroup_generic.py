r"""
Arithmetic subgroups (finite index subgroups of `{\rm SL}_2(\ZZ)`)
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
################################################################################


import sage.groups.old as group
from sage.rings.all import ZZ
import sage.rings.arith as arith
from sage.misc.cachefunc import cached_method
from copy import copy # for making copies of lists of cusps
from sage.modular.modsym.p1list import lift_to_sl2z
from sage.modular.cusps import Cusp

from sage.misc.lazy_import import lazy_import
lazy_import('sage.modular.arithgroup.congroup_sl2z', 'SL2Z')
from sage.structure.element import parent

from arithgroup_element import ArithmeticSubgroupElement

def is_ArithmeticSubgroup(x):
    r"""
    Return True if x is of type ArithmeticSubgroup.

    EXAMPLE::

        sage: from sage.modular.arithgroup.all import is_ArithmeticSubgroup
        sage: is_ArithmeticSubgroup(GL(2, GF(7)))
        False
        sage: is_ArithmeticSubgroup(Gamma0(4))
        True
    """

    return isinstance(x, ArithmeticSubgroup)


class ArithmeticSubgroup(group.Group):
    r"""
    Base class for arithmetic subgroups of `{\rm SL}_2(\ZZ)`. Not
    intended to be used directly, but still includes quite a few
    general-purpose routines which compute data about an arithmetic subgroup
    assuming that it has a working element testing routine.
    """

    Element = ArithmeticSubgroupElement

    def __init__(self):
        r"""
        Standard init routine.

        EXAMPLE::

            sage: G = Gamma1(7)
            sage: G.category() # indirect doctest
            Category of groups
        """
        group.Group.__init__(self)

    def _repr_(self):
        r"""
        Return the string representation of self.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup()._repr_()
            'Generic arithmetic subgroup of SL2Z'
        """
        return "Generic arithmetic subgroup of SL2Z"

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: Gamma1(7)._repr_option('element_ascii_art')
            True
        """
        if key == 'element_ascii_art':
            return True
        return super(ArithmeticSubgroup, self)._repr_option(key)

    def __reduce__(self):
        r"""
        Used for pickling self.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().__reduce__()
            Traceback (most recent call last):
            ...
            NotImplementedError: all subclasses must define a __reduce__ method
        """
        raise NotImplementedError("all subclasses must define a __reduce__ method")

    def _element_constructor_(self, x, check=True):
        r"""
        Create an element of this congruence subgroup from x.

        If the optional flag check is True (default), check whether
        x actually gives an element of self.

        EXAMPLES::

            sage: G = Gamma(5)
            sage: G([1, 0, -10, 1]) # indirect doctest
            [ 1   0]
            [-10  1]
            sage: G(matrix(ZZ, 2, [26, 5, 5, 1]))
            [26  5]
            [ 5  1]
            sage: G([1, 1, 6, 7])
            Traceback (most recent call last):
            ...
            TypeError: matrix [1 1]
            [6 7] is not an element of Congruence Subgroup Gamma(5)
        """
        # Do not override this function! Derived classes should override
        # _contains_sl2.
        x = SL2Z(x, check)
        if not check or x in self:
            return x
        raise TypeError("matrix %s is not an element of %s" % (x, self))

    def __contains__(self, x):
        r"""
        Test if x is an element of this group. This checks that x defines (is?) a 2x2 integer matrix of determinant 1, and
        then hands over to the routine _contains_sl2, which derived classes should implement.

        EXAMPLES::

            sage: [1,2] in SL2Z # indirect doctest
            False
            sage: [1,2,0,1] in SL2Z # indirect doctest
            True
            sage: SL2Z([1,2,0,1]) in Gamma(3) # indirect doctest
            False
            sage: -1 in SL2Z
            True
            sage: 2 in SL2Z
            False
        """
        # Do not override this function! Derived classes should override
        # _contains_sl2.
        if isinstance(x, type([])) and len(x) == 4:
            if not (x[0] in ZZ and x[1] in ZZ and x[2] in ZZ and x[3] in ZZ):
                return False
            a,b,c,d = map(ZZ, x)
            if a*d - b*c != 1: return False
            return self._contains_sl2(a,b,c,d)
        else:
            if parent(x) is not SL2Z:
                try:
                    y = SL2Z(x)
                except TypeError:
                    return False
                x = y
            return self._contains_sl2(x.a(),x.b(),x.c(),x.d())

    def _contains_sl2(self, a,b,c,d):
        r"""
        Test whether the matrix [a,b;c,d], which may be assumed to have
        determinant 1, is an element of self. This must be overridden by all
        subclasses.

        EXAMPLE::

            sage: G = sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup()
            sage: 1 in G
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
        """
        raise NotImplementedError("Please implement _contains_sl2 for %s" % self.__class__)

    def __hash__(self):
        r"""
        Return a hash of self.

        EXAMPLES::

            sage: Gamma0(11).__hash__()
            -545929996 # 32-bit
            466678398374495476 # 64-bit
            sage: Gamma1(11).__hash__()
            -830809815 # 32-bit
            4909638266971150633 # 64-bit
        """
        return hash(str(self))

    def is_parent_of(self, x):
        r"""
        Check whether this group is a valid parent for the element x. Required
        by Sage's testing framework.

        EXAMPLE::

            sage: Gamma(3).is_parent_of(ZZ(1))
            False
            sage: Gamma(3).is_parent_of([1,0,0,1])
            False
            sage: Gamma(3).is_parent_of(SL2Z([1,1,0,1]))
            False
            sage: Gamma(3).is_parent_of(SL2Z(1))
            True
        """
        return (parent(x) == SL2Z and x in self)

    def coset_reps(self, G=None):
        r"""
        Return right coset representatives for self \\ G, where G is another
        arithmetic subgroup that contains self.  If G = None, default to G =
        SL2Z.

        For generic arithmetic subgroups G this is carried out by Todd-Coxeter
        enumeration; here G is treated as a black box, implementing nothing but
        membership testing.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().coset_reps()
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.coset_reps(Gamma0(3))
            [
            [1 0]  [ 0 -1]  [ 0 -1]  [ 0 -1]
            [0 1], [ 1  0], [ 1  1], [ 1  2]
            ]
        """
        return self.todd_coxeter(G)[0]

    @cached_method
    def todd_coxeter(self, G=None, on_right=True):
        r"""
        Compute coset representatives for self \\ G and action of standard
        generators on them via Todd-Coxeter enumeration.

        If ``G`` is ``None``, default to ``SL2Z``. The method also computes
        generators of the subgroup at same time.

        INPUT:

        - ``G`` - intermediate subgroup (currently not implemented if diffferent
          from SL(2,Z))

        - ``on_right`` - boolean (default: True) - if True return right coset
          enumeration, if False return left one.

        This is *extremely* slow in general.

        OUTPUT:

        - a list of coset representatives

        - a list of generators for the group

        - ``l`` - list of integers that correspond to the action of the
          standard parabolic element [[1,1],[0,1]] of `SL(2,\ZZ)` on the cosets
          of self.

        - ``s`` - list of integers that correspond to the action of the standard
          element of order `2` [[0,-1],[1,0]] on the cosets of self.

        EXAMPLES::

            sage: L = SL2Z([1,1,0,1])
            sage: S = SL2Z([0,-1,1,0])

            sage: G = Gamma(2)
            sage: reps, gens, l, s = G.todd_coxeter()
            sage: len(reps) == G.index()
            True
            sage: all(reps[i] * L * ~reps[l[i]] in G for i in xrange(6))
            True
            sage: all(reps[i] * S * ~reps[s[i]] in G for i in xrange(6))
            True

            sage: G = Gamma0(7)
            sage: reps, gens, l, s = G.todd_coxeter()
            sage: len(reps) == G.index()
            True
            sage: all(reps[i] * L * ~reps[l[i]] in G for i in xrange(8))
            True
            sage: all(reps[i] * S * ~reps[s[i]] in G for i in xrange(8))
            True

            sage: G = Gamma1(3)
            sage: reps, gens, l, s = G.todd_coxeter(on_right=False)
            sage: len(reps) == G.index()
            True
            sage: all(~reps[l[i]] * L * reps[i] in G for i in xrange(8))
            True
            sage: all(~reps[s[i]] * S * reps[i] in G for i in xrange(8))
            True

            sage: G = Gamma0(5)
            sage: reps, gens, l, s = G.todd_coxeter(on_right=False)
            sage: len(reps) == G.index()
            True
            sage: all(~reps[l[i]] * L * reps[i] in G for i in xrange(6))
            True
            sage: all(~reps[s[i]] * S * reps[i] in G for i in xrange(6))
            True
        """
        if G is None:
            G = SL2Z
        if G != SL2Z:
            raise NotImplementedError("Don't know how to compute coset reps for subgroups yet")

        id = SL2Z([1,0,0,1])
        l = SL2Z([1,1,0,1])
        s = SL2Z([0,-1,1,0])

        reps = [id]       # coset representatives
        reps_inv = {id:0} # coset representatives index

        l_wait_back = [id] # rep with no incoming s_edge
        s_wait_back = [id] # rep with no incoming l_edge
        l_wait = [id]      # rep with no outgoing l_edge
        s_wait = [id]      # rep with no outgoing s_edge

        l_edges = [None]    # edges for l
        s_edges = [None]    # edges for s

        gens = []

        while l_wait or s_wait:
            if l_wait:
                x = l_wait.pop(0)
                y = x
                not_end = True
                while not_end:
                    if on_right:
                        y = y*l
                    else:
                        y = l*y
                    for i in xrange(len(l_wait_back)):
                        v = l_wait_back[i]
                        if on_right:
                            yy = y*~v
                        else:
                            yy = ~v*y
                        if yy in self:
                            l_edges[reps_inv[x]] = reps_inv[v]
                            del l_wait_back[i]
                            if yy != id:
                                gens.append(self(yy))
                            not_end = False
                            break
                    else:
                        reps_inv[y] = len(reps)
                        l_edges[reps_inv[x]] = len(reps)
                        reps.append(y)
                        l_edges.append(None)
                        s_edges.append(None)
                        s_wait_back.append(y)
                        s_wait.append(y)
                    x = y

            if s_wait:
                x = s_wait.pop(0)
                y = x
                not_end = True
                while not_end:
                    if on_right:
                        y = y*s
                    else:
                        y = s*y
                    for i in xrange(len(s_wait_back)):
                        v = s_wait_back[i]
                        if on_right:
                            yy = y*~v
                        else:
                            yy = ~v*y
                        if yy in self:
                            s_edges[reps_inv[x]] = reps_inv[v]
                            del s_wait_back[i]
                            if yy != id:
                                gens.append(self(yy))
                            not_end = False
                            break
                    else:
                        reps_inv[y] = len(reps)
                        s_edges[reps_inv[x]] = len(reps)
                        reps.append(y)
                        l_edges.append(None)
                        s_edges.append(None)
                        l_wait_back.append(y)
                        l_wait.append(y)
                    x = y

        return reps, gens, l_edges, s_edges

    def nu2(self):
        r"""
        Return the number of orbits of elliptic points of order 2 for this
        arithmetic subgroup.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().nu2()
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.nu2(Gamma0(1105)) == 8
            True
        """

        # Subgroups not containing -1 have no elliptic points of order 2.

        if not self.is_even():
            return 0

        # Cheap trick: if self is a subgroup of something with no elliptic points,
        # then self has no elliptic points either.

        from all import Gamma0, is_CongruenceSubgroup
        if is_CongruenceSubgroup(self):
            if self.is_subgroup(Gamma0(self.level())) and Gamma0(self.level()).nu2() == 0:
                return 0

        # Otherwise, the number of elliptic points is the number of g in self \
        # SL2Z such that the stabiliser of g * i in self is not trivial. (Note
        # that the points g*i for g in the coset reps are not distinct, but it
        # still works, since the failure of these points to be distinct happens
        # precisely when the preimages are not elliptic.)

        count = 0
        for g in self.coset_reps():
            if g * SL2Z([0,1,-1,0]) * (~g) in self:
                count += 1
        return count

    def nu3(self):
        r"""
        Return the number of orbits of elliptic points of order 3 for this
        arithmetic subgroup.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().nu3()
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.nu3(Gamma0(1729)) == 8
            True

        We test that a bug in handling of subgroups not containing -1 is fixed: ::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.nu3(GammaH(7, [2]))
            2
        """

        # Cheap trick: if self is a subgroup of something with no elliptic points,
        # then self has no elliptic points either.

        from all import Gamma0, is_CongruenceSubgroup
        if is_CongruenceSubgroup(self):
            if self.is_subgroup(Gamma0(self.level())) and Gamma0(self.level()).nu3() == 0:
                return 0

        count = 0
        for g in self.coset_reps():
            if g * SL2Z([0,1,-1,-1]) * (~g) in self:
                count += 1

        if self.is_even():
            return count
        else:
            return count // 2

    def __cmp__(self, other):
        r"""
        Compare self to other.

        NOTE: This function must be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().__cmp__(ZZ)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_abelian(self):
        r"""
        Return True if this arithmetic subgroup is abelian.

        Since arithmetic subgroups are always nonabelian, this always
        returns False.

        EXAMPLES::

            sage: SL2Z.is_abelian()
            False
            sage: Gamma0(3).is_abelian()
            False
            sage: Gamma1(12).is_abelian()
            False
            sage: GammaH(4, [3]).is_abelian()
            False
        """
        return False

    def is_finite(self):
        r"""
        Return True if this arithmetic subgroup is finite.

        Since arithmetic subgroups are always infinite, this always
        returns False.

        EXAMPLES::

            sage: SL2Z.is_finite()
            False
            sage: Gamma0(3).is_finite()
            False
            sage: Gamma1(12).is_finite()
            False
            sage: GammaH(4, [3]).is_finite()
            False
        """
        return False

    def is_subgroup(self, right):
        r"""
        Return True if self is a subgroup of right, and False otherwise. For
        generic arithmetic subgroups this is done by the absurdly slow
        algorithm of checking all of the generators of self to see if they are
        in right.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().is_subgroup(SL2Z)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.is_subgroup(Gamma1(18), Gamma0(6))
            True
        """
        # ridiculously slow generic algorithm

        w = self.gens()
        for g in w:
            if not (g in right):
                return False
        return True

    def is_normal(self):
        r"""
        Return True precisely if this subgroup is a normal subgroup of SL2Z.

        EXAMPLES::

            sage: Gamma(3).is_normal()
            True
            sage: Gamma1(3).is_normal()
            False
        """
        for x in self.gens():
            for y in SL2Z.gens():
                if y*SL2Z(x)*(~y) not in self:
                    return False
        return True

    def is_odd(self):
        r"""
        Return True precisely if this subgroup does not contain the
        matrix -1.

        EXAMPLES::

            sage: SL2Z.is_odd()
            False
            sage: Gamma0(20).is_odd()
            False
            sage: Gamma1(5).is_odd()
            True
            sage: GammaH(11, [3]).is_odd()
            True
        """
        return not self.is_even()

    def is_even(self):
        r"""
        Return True precisely if this subgroup contains the matrix -1.

        EXAMPLES::

            sage: SL2Z.is_even()
            True
            sage: Gamma0(20).is_even()
            True
            sage: Gamma1(5).is_even()
            False
            sage: GammaH(11, [3]).is_even()
            False
        """
        return [-1, 0, 0, -1] in self

    def to_even_subgroup(self):
        r"""
        Return the smallest even subgroup of `SL(2, \ZZ)` containing self.

        EXAMPLE::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().to_even_subgroup()
            Traceback (most recent call last):
            ...
            NotImplementedError: Please implement _contains_sl2 for <class 'sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup_with_category'>
        """
        if self.is_even():
            return self
        else:
            raise NotImplementedError

    def order(self):
        r"""
        Return the number of elements in this arithmetic subgroup.

        Since arithmetic subgroups are always infinite, this always returns
        infinity.

        EXAMPLES::

            sage: SL2Z.order()
            +Infinity
            sage: Gamma0(5).order()
            +Infinity
            sage: Gamma1(2).order()
            +Infinity
            sage: GammaH(12, [5]).order()
            +Infinity
        """
        from sage.rings.infinity import infinity
        return infinity

    def reduce_cusp(self, c):
        r"""
        Given a cusp `c \in \mathbb{P}^1(\QQ)`, return the unique reduced cusp
        equivalent to c under the action of self, where a reduced cusp is an
        element `\tfrac{r}{s}` with r,s coprime non-negative integers, s as
        small as possible, and r as small as possible for that s.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().reduce_cusp(1/4)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def cusps(self, algorithm='default'):
        r"""
        Return a sorted list of inequivalent cusps for self, i.e. a set of
        representatives for the orbits of self on `\mathbb{P}^1(\QQ)`.
        These should be returned in a reduced form where this makes sense.

        INPUTS:

        - ``algorithm`` -- which algorithm to use to compute the cusps of self.
          ``'default'`` finds representatives for a known complete set of
          cusps. ``'modsym'`` computes the boundary map on the space of weight
          two modular symbols associated to self, which finds the cusps for
          self in the process.

        EXAMPLES::

            sage: Gamma0(36).cusps()
            [0, 1/18, 1/12, 1/9, 1/6, 1/4, 1/3, 5/12, 1/2, 2/3, 5/6, Infinity]
            sage: Gamma0(36).cusps(algorithm='modsym') == Gamma0(36).cusps()
            True
            sage: GammaH(36, [19,29]).cusps() == Gamma0(36).cusps()
            True
            sage: Gamma0(1).cusps()
            [Infinity]
        """
        try:
            return copy(self._cusp_list[algorithm])
        except (AttributeError,KeyError):
            self._cusp_list = {}

        from congroup_sl2z import is_SL2Z
        if is_SL2Z(self):
            s = [Cusp(1,0)]

        if algorithm == 'default':
            s = self._find_cusps()
        elif algorithm == 'modsym':
            s = sorted([self.reduce_cusp(c) for c in self.modular_symbols().cusps()])
        else:
            raise ValueError("unknown algorithm: %s"%algorithm)

        self._cusp_list[algorithm] = s
        return copy(s)

    def _find_cusps(self):
        r"""
        Calculate a list of inequivalent cusps.

        EXAMPLES::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5)._find_cusps()
            Traceback (most recent call last):
            ...
            NotImplementedError

        NOTE: There is a generic algorithm implemented at the top level that
        uses the coset representatives of self. This is *very slow* and for all
        the standard congruence subgroups there is a quicker way of doing it,
        so this should usually be overridden in subclasses; but it doesn't have
        to be.
        """
        i = Cusp([1,0])
        L = [i]
        for a in self.coset_reps():
            ai = i.apply([a.a(), a.b(), a.c(), a.d()])
            new = 1
            for v in L:
                if self.are_equivalent(ai, v):
                    new = 0
                    break
            if new == 1:
                L.append(ai)
        return L

    def are_equivalent(self, x, y, trans = False):
        r"""
        Test whether or not cusps x and y are equivalent modulo self.  If self
        has a reduce_cusp() method, use that; otherwise do a slow explicit
        test.

        If trans = False, returns True or False. If trans = True, then return
        either False or an element of self mapping x onto y.

        EXAMPLE::

            sage: Gamma0(7).are_equivalent(Cusp(1/3), Cusp(0), trans=True)
            [  3  -1]
            [-14   5]
            sage: Gamma0(7).are_equivalent(Cusp(1/3), Cusp(1/7))
            False
        """
        x = Cusp(x)
        y = Cusp(y)
        if not trans:
            try:
                xr = self.reduce_cusp(x)
                yr = self.reduce_cusp(y)
                if xr != yr:
                    return False
                if xr == yr:
                    return True
            except NotImplementedError:
                pass

        vx = lift_to_sl2z(x.numerator(),x.denominator(), 0)
        dx = SL2Z([vx[2], -vx[0], vx[3], -vx[1]])
        vy = lift_to_sl2z(y.numerator(),y.denominator(), 0)
        dy = SL2Z([vy[2], -vy[0], vy[3], -vy[1]])

        for i in xrange(self.index()):
            # Note that the width of any cusp is bounded above by the index of self.
            # If self is congruence, then the level of self is a much better bound, but
            # this method is written to work with non-congruence subgroups as well,
            if dy * SL2Z([1,i,0,1])*(~dx) in self:
                if trans:
                    return dy * SL2Z([1,i,0,1]) * ~dx
                else:
                    return True
            elif (self.is_odd() and dy * SL2Z([-1,-i,0,-1]) * ~dx in self):
                if trans:
                    return dy * SL2Z([-1,-i,0,-1]) * ~dx
                else:
                    return True
        return False

    def cusp_data(self, c):
        r"""
        Return a triple (g, w, t) where g is an element of self generating the
        stabiliser of the given cusp, w is the width of the cusp, and t is 1 if
        the cusp is regular and -1 if not.

        EXAMPLES::

            sage: Gamma1(4).cusp_data(Cusps(1/2))
            (
            [ 1 -1]
            [ 4 -3], 1, -1
            )
        """
        c = Cusp(c)

        # first find an element of SL2Z sending infinity to the given cusp
        w = lift_to_sl2z(c.denominator(), c.numerator(), 0)
        g = SL2Z([w[3], w[1], w[2],w[0]])

        for d in xrange(1,1+self.index()):
            if g * SL2Z([1,d,0,1]) * (~g) in self:
                return (g * SL2Z([1,d,0,1]) * (~g), d, 1)
            elif g * SL2Z([-1,-d,0,-1]) * (~g) in self:
                return (g * SL2Z([-1,-d,0,-1]) * (~g), d, -1)
        raise ArithmeticError("Can't get here!")

    def is_regular_cusp(self, c):
        r"""
        Return True if the orbit of the given cusp is a regular cusp for self,
        otherwise False. This is automatically true if -1 is in self.

        EXAMPLES::

            sage: Gamma1(4).is_regular_cusp(Cusps(1/2))
            False
            sage: Gamma1(4).is_regular_cusp(Cusps(oo))
            True
        """
        if self.is_even(): return True
        return (self.cusp_data(c)[2] == 1)

    def cusp_width(self, c):
        r"""
        Return the width of the orbit of cusps represented by c.

        EXAMPLES::

            sage: Gamma0(11).cusp_width(Cusps(oo))
            1
            sage: Gamma0(11).cusp_width(0)
            11
            sage: [Gamma0(100).cusp_width(c) for c in Gamma0(100).cusps()]
            [100, 1, 4, 1, 1, 1, 4, 25, 1, 1, 4, 1, 25, 4, 1, 4, 1, 1]
        """
        return self.cusp_data(c)[1]

    def index(self):
        r"""
        Return the index of self in the full modular group.

        EXAMPLES::

            sage: Gamma0(17).index()
            18
            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5).index()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """

        return len(list(self.coset_reps()))

    def generalised_level(self):
        r"""
        Return the generalised level of self, i.e. the least common multiple of
        the widths of all cusps.

        If self is *even*, Wohlfart's theorem tells us that this is equal to
        the (conventional) level of self when self is a congruence subgroup.
        This can fail if self is odd, but the actual level is at most twice the
        generalised level. See the paper by Kiming, Schuett and Verrill for
        more examples.

        EXAMPLE::

            sage: Gamma0(18).generalised_level()
            18
            sage: sage.modular.arithgroup.arithgroup_perm.HsuExample18().generalised_level()
            24

        In the following example, the actual level is twice the generalised
        level. This is the group `G_2` from Example 17 of K-S-V.

        ::

            sage: G = CongruenceSubgroup(8, [ [1,1,0,1], [3,-1,4,-1] ])
            sage: G.level()
            8
            sage: G.generalised_level()
            4
        """
        return arith.lcm([self.cusp_width(c) for c in self.cusps()])

    def projective_index(self):
        r"""
        Return the index of the image of self in `{\rm PSL}_2(\ZZ)`. This is equal
        to the index of self if self contains -1, and half of this otherwise.

        This is equal to the degree of the natural map from the modular curve
        of self to the `j`-line.

        EXAMPLE::

            sage: Gamma0(5).projective_index()
            6
            sage: Gamma1(5).projective_index()
            12
        """

        if self.is_even():
            return self.index()
        else:
            return self.index() // 2

    def is_congruence(self):
        r"""
        Return True if self is a congruence subgroup.

        EXAMPLE::

            sage: Gamma0(5).is_congruence()
            True
            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup().is_congruence()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """

        raise NotImplementedError

    def genus(self):
        r"""
        Return the genus of the modular curve of self.

        EXAMPLES::

            sage: Gamma1(5).genus()
            0
            sage: Gamma1(31).genus()
            26
            sage: Gamma1(157).genus() == dimension_cusp_forms(Gamma1(157), 2)
            True
            sage: GammaH(7, [2]).genus()
            0
            sage: [Gamma0(n).genus() for n in [1..23]]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 2, 2]
            sage: [n for n in [1..200] if Gamma0(n).genus() == 1]
            [11, 14, 15, 17, 19, 20, 21, 24, 27, 32, 36, 49]


        """

        return ZZ(1 + (self.projective_index()) / ZZ(12)  - (self.nu2())/ZZ(4) - (self.nu3())/ZZ(3) - self.ncusps()/ZZ(2))

    def farey_symbol(self):
        r"""
        Return the Farey symbol associated to this subgroup. See the
        :mod:`~sage.modular.arithgroup.farey_symbol` module for more
        information.

        EXAMPLE::

            sage: Gamma1(4).farey_symbol()
            FareySymbol(Congruence Subgroup Gamma1(4))
        """
        from farey_symbol import Farey
        return Farey(self)

    @cached_method
    def generators(self, algorithm="farey"):
        r"""
        Return a list of generators for this congruence subgroup. The result is cached.

        INPUT:

        - ``algorithm`` (string): either ``farey`` or ``todd-coxeter``.

        If ``algorithm`` is set to ``"farey"``, then the generators will be
        calculated using Farey symbols, which will always return a *minimal*
        generating set. See :mod:`~sage.modular.arithgroup.farey_symbol` for
        more information.

        If ``algorithm`` is set to ``"todd-coxeter"``, a simpler algorithm
        based on Todd-Coxeter enumeration will be used. This is *exceedingly*
        slow for general subgroups, and the list of generators will be far from
        minimal (indeed it may contain repetitions).

        EXAMPLE::

            sage: Gamma(2).generators()
            [
            [1 2]  [ 3 -2]  [-1  0]
            [0 1], [ 2 -1], [ 0 -1]
            ]
            sage: Gamma(2).generators(algorithm="todd-coxeter")
            [
            [1 2]  [-1  0]  [ 1  0]  [-1  0]  [-1  2]  [-1  0]  [1 0]
            [0 1], [ 0 -1], [-2  1], [ 0 -1], [-2  3], [ 2 -1], [2 1]
            ]
        """
        if algorithm=="farey":
            return self.farey_symbol().generators()
        elif algorithm == "todd-coxeter":
            return self.todd_coxeter()[1]
        else:
            raise ValueError("Unknown algorithm '%s' (should be either 'farey' or 'todd-coxeter')" % algorithm)

    def gens(self, *args, **kwds):
        r"""
        Return a tuple of generators for this congruence subgroup.

        The generators need not be minimal. For arguments, see :meth:`~generators`.

        EXAMPLES::

            sage: SL2Z.gens()
            (
            [ 0 -1]  [1 1]
            [ 1  0], [0 1]
            )
        """
        return tuple(self.generators(*args, **kwds))

    def gen(self, i):
        r"""
        Return the i-th generator of self, i.e. the i-th element of the
        tuple self.gens().

        EXAMPLES::

            sage: SL2Z.gen(1)
            [1 1]
            [0 1]
        """
        return self.generators()[i]

    def ngens(self):
        r"""
        Return the size of the minimal generating set of self returned by
        :meth:`generators`.

        EXAMPLES::

            sage: Gamma0(22).ngens()
            8
            sage: Gamma1(14).ngens()
            13
            sage: GammaH(11, [3]).ngens()
            3
            sage: SL2Z.ngens()
            2
        """
        return len(self.generators())

    def ncusps(self):
        r"""
        Return the number of cusps of this arithmetic subgroup. This is
        provided as a separate function since for dimension formulae in even
        weight all we need to know is the number of cusps, and this can be
        calculated very quickly, while enumerating all cusps is much slower.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.ncusps(Gamma0(7))
            2
        """

        return ZZ(len(self.cusps()))

    def nregcusps(self):
        r"""
        Return the number of cusps of self that are "regular", i.e. their
        stabiliser has a generator with both eigenvalues +1 rather than -1. If
        the group contains -1, every cusp is clearly regular.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.nregcusps(Gamma1(4))
            2
        """
        return self.ncusps() - self.nirregcusps()

    def nirregcusps(self):
        r"""
        Return the number of cusps of self that are "irregular", i.e. their
        stabiliser can only be generated by elements with both eigenvalues -1
        rather than +1. If the group contains -1, every cusp is clearly
        regular.

        EXAMPLES::

            sage: sage.modular.arithgroup.arithgroup_generic.ArithmeticSubgroup.nirregcusps(Gamma1(4))
            1
        """
        if self.is_even():
            return 0
        else:
            return ZZ(len([c for c in self.cusps() if not self.is_regular_cusp(c)]))

    def dimension_modular_forms(self, k=2):
        r"""
        Return the dimension of the space of weight k modular forms for this
        group. This is given by a standard formula in terms of k and various
        invariants of the group; see Diamond + Shurman, "A First Course in
        Modular Forms", section 3.5 and 3.6. If k is not given, defaults to k =
        2.

        For dimensions of spaces of modular forms with character for Gamma1, use
        the standalone function dimension_modular_forms().

        For weight 1 modular forms this function only works in cases where one
        can prove solely in terms of Riemann-Roch theory that there aren't any
        cusp forms (i.e. when the number of regular cusps is strictly greater
        than the degree of the canonical divisor). Otherwise a
        NotImplementedError is raised.

        EXAMPLE::

            sage: Gamma1(31).dimension_modular_forms(2)
            55
            sage: Gamma1(3).dimension_modular_forms(1)
            1
            sage: Gamma1(4).dimension_modular_forms(1) # irregular cusp
            1
            sage: Gamma1(31).dimension_modular_forms(1)
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of dimensions of weight 1 cusp forms spaces not implemented in general
        """

        k = ZZ(k)
        if k < 0: return ZZ(0)
        if k == 0: return ZZ(1)

        if not (k % 2):
            # k even

            return (k-1) * (self.genus() - 1) + (k // ZZ(4))*self.nu2() + (k // ZZ(3))*self.nu3() + (k // ZZ(2))*self.ncusps()

        else:
            # k odd
            if self.is_even():
                return ZZ(0)
            else:
                e_reg = self.nregcusps()
                e_irr = self.nirregcusps()

                if k > 1:
                    return (k-1)*(self.genus()-1) + (k // ZZ(3)) * self.nu3() + (k * e_reg)/ZZ(2) + (k-1)/ZZ(2) * e_irr
                else:
                    if e_reg > 2*self.genus() - 2:
                        return e_reg / ZZ(2)
                    else:
                        raise NotImplementedError("Computation of dimensions of weight 1 modular forms spaces not implemented in general")

    def dimension_cusp_forms(self, k=2):
        r"""
        Return the dimension of the space of weight k cusp forms for this
        group. This is given by a standard formula in terms of k and various
        invariants of the group; see Diamond + Shurman, "A First Course in
        Modular Forms", section 3.5 and 3.6. If k is not given, default to k =
        2.

        For dimensions of spaces of cusp forms with character for Gamma1, use
        the standalone function dimension_cusp_forms().

        For weight 1 cusp forms this function only works in cases where one can
        prove solely in terms of Riemann-Roch theory that there aren't any cusp
        forms (i.e. when the number of regular cusps is strictly greater than
        the degree of the canonical divisor). Otherwise a NotImplementedError is
        raised.

        EXAMPLE::

            sage: Gamma1(31).dimension_cusp_forms(2)
            26
            sage: Gamma1(3).dimension_cusp_forms(1)
            0
            sage: Gamma1(4).dimension_cusp_forms(1) # irregular cusp
            0
            sage: Gamma1(31).dimension_cusp_forms(1)
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of dimensions of weight 1 cusp forms spaces not implemented in general
        """
        k = ZZ(k)
        if k <= 0: return ZZ(0)

        if not (k % 2):
            # k even

            if k == 2:
                return self.genus()

            else:
                return (k-1) * (self.genus() - 1) + (k // ZZ(4))*self.nu2() + (k // ZZ(3))*self.nu3() + (k // ZZ(2) - 1)*self.ncusps()

        else:
            # k odd

            if self.is_even():
                return ZZ(0)

            else:
                e_reg = self.nregcusps()
                e_irr = self.nirregcusps()

                if k > 1:
                    return (k-1)*(self.genus()-1) + (k // ZZ(3)) * self.nu3() + (k-2)/ZZ(2) * e_reg + (k-1)/ZZ(2) * e_irr
                else:
                    if e_reg > 2*self.genus() - 2:
                        return ZZ(0)
                    else:
                        raise NotImplementedError("Computation of dimensions of weight 1 cusp forms spaces not implemented in general")

    def dimension_eis(self, k=2):
        r"""
        Return the dimension of the space of weight k Eisenstein series for
        this group, which is a subspace of the space of modular forms
        complementary to the space of cusp forms.

        INPUT:

        - ``k`` - an integer (default 2).

        EXAMPLES::

            sage: GammaH(33,[2]).dimension_eis()
            7
            sage: GammaH(33,[2]).dimension_eis(3)
            0
            sage: GammaH(33, [2,5]).dimension_eis(2)
            3
            sage: GammaH(33, [4]).dimension_eis(1)
            4
        """

        if k < 0: return ZZ(0)
        if k == 0: return ZZ(1)

        if not (k % 2): # k even
            if k > 2:
                return self.ncusps()
            else: # k = 2
                return self.ncusps() - 1

        else: # k odd
            if self.is_even():
                return 0
            if k > 1:
                return self.nregcusps()
            else: # k = 1
                return ZZ(self.nregcusps()/ ZZ(2))

    def as_permutation_group(self):
        r"""
        Return self as an arithmetic subgroup defined in terms of the
        permutation action of `SL(2,\ZZ)` on its right cosets.

        This method uses Todd-Coxeter enumeration (via the method
        :meth:`~todd_coxeter`) which can be extremely slow for arithmetic
        subgroups with relatively large index in `SL(2,\ZZ)`.

        EXAMPLES::

            sage: G = Gamma(3)
            sage: P = G.as_permutation_group(); P
            Arithmetic subgroup of index 24
            sage: G.ncusps() == P.ncusps()
            True
            sage: G.nu2() == P.nu2()
            True
            sage: G.nu3() == P.nu3()
            True
            sage: G.an_element() in P
            True
            sage: P.an_element() in G
            True
        """
        _,_,l_edges,s2_edges=self.todd_coxeter()
        n = len(l_edges)
        s3_edges = [None] * n
        r_edges = [None] * n
        for i in xrange(n):
            ii = s2_edges[l_edges[i]]
            s3_edges[ii] = i
            r_edges[ii] = s2_edges[i]
        if self.is_even():
            from sage.modular.arithgroup.arithgroup_perm import EvenArithmeticSubgroup_Permutation
            g=EvenArithmeticSubgroup_Permutation(S2=s2_edges,S3=s3_edges,L=l_edges,R=r_edges)
        else:
            from sage.modular.arithgroup.arithgroup_perm import OddArithmeticSubgroup_Permutation
            g=OddArithmeticSubgroup_Permutation(S2=s2_edges,S3=s3_edges,L=l_edges,R=r_edges)
        g.relabel()
        return g

    def sturm_bound(self, weight=2):
        r"""
        Returns the Sturm bound for modular forms of the given weight and level
        this subgroup.

        INPUT:

        -  ``weight`` - an integer `\geq 2` (default: 2)

        EXAMPLES::
            sage: Gamma0(11).sturm_bound(2)
            2
            sage: Gamma0(389).sturm_bound(2)
            65
            sage: Gamma0(1).sturm_bound(12)
            1
            sage: Gamma0(100).sturm_bound(2)
            30
            sage: Gamma0(1).sturm_bound(36)
            3
            sage: Gamma0(11).sturm_bound()
            2
            sage: Gamma0(13).sturm_bound()
            3
            sage: Gamma0(16).sturm_bound()
            4
            sage: GammaH(16,[13]).sturm_bound()
            8
            sage: GammaH(16,[15]).sturm_bound()
            16
            sage: Gamma1(16).sturm_bound()
            32
            sage: Gamma1(13).sturm_bound()
            28
            sage: Gamma1(13).sturm_bound(5)
            70

        FURTHER DETAILS: This function returns a positive integer
        `n` such that the Hecke operators
        `T_1,\ldots, T_n` acting on *cusp forms* generate the
        Hecke algebra as a `\ZZ`-module when the character
        is trivial or quadratic. Otherwise, `T_1,\ldots,T_n`
        generate the Hecke algebra at least as a
        `\ZZ[\varepsilon]`-module, where
        `\ZZ[\varepsilon]` is the ring generated by the
        values of the Dirichlet character `\varepsilon`.
        Alternatively, this is a bound such that if two cusp forms
        associated to this space of modular symbols are congruent modulo
        `(\lambda, q^n)`, then they are congruent modulo
        `\lambda`.

        REFERENCES:

        - See the Agashe-Stein appendix to Lario and Schoof,
          *Some computations with Hecke rings and deformation rings*,
          Experimental Math., 11 (2002), no. 2, 303-311.

        - This result originated in the paper Sturm,
          *On the congruence of modular forms*,
          Springer LNM 1240, 275-280, 1987.

        REMARK: Kevin Buzzard pointed out to me (William Stein) in Fall
        2002 that the above bound is fine for `\Gamma_1(N)` with
        character, as one sees by taking a power of `f`. More
        precisely, if `f \cong 0 \pmod{p}` for first
        `s` coefficients, then `f^r \cong 0 \pmod{p}` for
        first `sr` coefficients. Since the weight of `f^r`
        is `r\cdot k(f)`, it follows that if
        `s \geq b`, where `b` is the Sturm bound for
        `\Gamma_0(N)` at weight `k(f)`, then `f^r`
        has valuation large enough to be forced to be `0` at
        `r*k(f)` by Sturm bound (which is valid if we choose
        `r` correctly). Thus `f \cong 0 \pmod{p}`.
        Conclusion: For `\Gamma_1(N)` with fixed character, the
        Sturm bound is *exactly* the same as for `\Gamma_0(N)`.

        A key point is that we are finding
        `\ZZ[\varepsilon]` generators for the Hecke algebra
        here, not `\ZZ`-generators. So if one wants
        generators for the Hecke algebra over `\ZZ`, this
        bound must be suitably modified (and I'm not sure what the
        modification is).

        AUTHORS:

        - William Stein
        """
        return ZZ((self.index() * weight / ZZ(12)).ceil())

# below follow the additions around the integrality algorithm
# Copyright (C) 2011 by    Georg S. Weber
# Distributed under the terms of the GNU General Public License (GPL)
# (either version 2 of the license, or any later version, at your option)

    def z_an_introduction(self):
        r"""
        AUTHOR:

        - Georg S. Weber (February 2011): Initial version

        William Stein wrote in chapter 6 of an early preliminary draft
        (a file named "257.ps", 151 pages, dated "December 16, 2004"
        titled "Algorithms For Computing With Modular Forms") of his
        book "Explicitly Computing Modular Forms" (there this chapter
        became chapter 8. I could not find the following remark in the
        book):

        'Remark 6.6.11. There is rumored to be a "geometric" way to
        compute a presentation for M_2(Gamma0(N)) more directly,
        without resorting to general linear algebra techniques. I am
        unaware of such a method having ever been published, but it
        was sketched to me independently by Georg Weber in 1999 and
        Robert Pollack in 2004. The computations we do after computing
        a presentation for M_2(Gamma0(N)) are usually significantly
        more time consuming than computation of a presentation in the
        first place, so it's unclear how useful this algorithm would
        be in practice. (I have not heard of a method for directly
        obtaining a presentation for M_k(Gamma0(N)).)'

        Robert Pollack's code (written in magma) is available e.g. as
        part of the Stark-Heegner Package (shp) under:
        http://www.math.mcgill.ca/darmon/programs/shp/shp.html and
        there's a years old trac ticket to include these algorithms in
        Sage http://trac.sagemath.org/sage_trac/ticket/812

        The following exposition and code is independent of his
        (although on a technical level, related), and centers around
        the integrality algorithm.

        The only published reference I know of for the integrality
        algorithm is the "nice trick" mentioned under the third point
        of Remark 3.10 in: G. Frey and M. M"uller, Arithmetic of
        modular curves and applications, Algorithmic algebra and
        number theory (Heidelberg, 1997), Springer, Berlin, 1999,
        pp. 11-48 which reads:

        "3. One problem to calculate in the ZZ-module M (see 18) is to find a
            basis for this space. Therefore we have to use Gauss elimination or
            sparse matrix techniques to find a basis which can be rather painful
            since the matrices are very large. But with a nice trick due to
            X. Wang (see [M"ul98]) we have a very efficient implementation,
            which avoids entirely such matrix operations."

        Unfortunately, the reference [M"ul98] cited there is to my knowledge
        still unpublished, and I never could persuade Michael M"uller to share
        some of his code with me. On the other hand, he did explain to me the
        basic idea, how one e.g. can deduce the "binary tree heuristic" from
        that, and did outline the generalizaton to higher weights k > 2..

        The foundation, on which the integrality algorithm relies on, is a
        certain enumeration of the cosets (or a set of coset representatives;
        where we always mean "right" cosets), for arithmetic subgroups G.
        These enumerations have, what one could call a "least index" property
        for the 2- and 3-relations.

        Let us work with the example G = Gamma0(11), an even group of
        index 12 in SL2Z. Here is such a list of a full set of its
        coset representatives::

            sage: reps, gens, basis = Gamma0(11).z_sl2z_todd_coxeter(); reps
            [[1 0]
            [0 1], [-1  1]
            [-1  0], [ 0 -1]
            [ 1 -1], [-1  0]
            [-1 -1], [ 0  1]
            [-1  2], [ 1 -1]
            [ 2 -1], [1 0]
            [2 1], [-1 -1]
            [-1 -2], [ 0 -1]
            [ 1 -3], [-1  1]
            [-3  2], [-1  2]
            [-2  3], [ 2 -1]
            [ 3 -1]]
            sage: len(reps)
            12

        From now, we refer to elements of this list by their index there.
        For Gamma0(11), both nu2 and nu3 vanish, so we have six 2-relations
        (none of which is degenerated)::

            sage: Gamma0(11).z_wellformed_2_relations()
            [[0, 1], [2, 3], [4, 6], [5, 7], [8, 10], [9, 11]]

        and four 3-relations (neither of which is degenerated, too)::

            sage: Gamma0(11).z_wellformed_3_relations()
            [[0, 2, 1], [3, 4, 5], [6, 8, 9], [7, 10, 11]]

        Each index (coset) appears exactly once in some 2-relation, and also
        exactly once in some 3-relation.

        Let RAW(G) be the ZZ-module of rank 12 generated "on" the 12 cosets,
        and let REL(G) by the ZZ-submodule generated by the 2- and 3-relations.
        Obviously "0" + "1" +... + "11" is generated by both 2-relations and
        3-relations alone, so the dimension will be lower than six plus four.
        It is well known, that this is all the cancelling that occurs, so the
        dimension of REL(G) is equal to nine (one less than six plus four).
        The integrality algorithm shows among other things, that REL(G) is
        saturated in RAW(G), so M_2(G, ZZ) = RAW(G)/REL(G) is a free ZZ-module
        of rank three.
        It is always possible to find generators of RAW(G) that become basis
        elements of M_2(G, ZZ), in our case::

            sage: reps, gens, basis = Gamma0(11).z_sl2z_todd_coxeter(); basis
            [1, 10, 11]

        Explicitly, the ZZ-module generated by this basis is saturated (i.e.
        has trivial index one) in RAW(G).

        The integrality algorithm now is a "nice trick" to reduce in a very
        efficient way an arbitrary element of RAW(G) modulo REL(G), without
        relying on general linear algebra techniques, using the following
        property of the coset enumeration used.

        If we throw away one of the relations, the others are generators of
        REL(11), let us throw away the 2-relation "containing" the identity
        2x2 matrix, to be more precise the coset with index 0, i.e. [0, 1]
        (see above).
        Then sort the other relations in ascending order, with respect to the
        cosets of least index they "contain"::

            sage: Gamma0(11).z_wellformed_all_relations_but_one()
            [[0, 2, 1],
             [2, 3],
             [3, 4, 5],
             [4, 6],
             [5, 7],
             [6, 8, 9],
             [7, 10, 11],
             [8, 10],
             [9, 11]]

        The important property of the coset enumeration used is the fact, that
        each index occurs at most once as the "least" index "contained" in some
        relation (this clearly cannot hold if we add [0, 1] back; we must have
        removed one of the two relations with the index 0).

        The generators to become basis elements of M_2(G, ZZ) are simply the
        cosets which do not appear in some relation as a "least index". ::

            sage: reps, gens, basis = Gamma0(11).z_sl2z_todd_coxeter(); basis
            [1, 10, 11]

        The integrality algorithm now simply runs once through this list of
        sorted relations, and applies each relation once. This "cleans" the
        index this relation has as least index, and by construction, this
        index will never get "dirtied" again. The running time is obviously
        O( G.index() ), and by the nature of the relations, the integrality
        algorithm works over ZZ, or arbitrary modules. ::

            sage: Pol12 = PolynomialRing(ZZ, "X", 12); Pol12
            Multivariate Polynomial Ring in X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11 over Integer Ring
            sage: indata = [Pol12("X"+str(i)) for i in range(12)]; indata
            [X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11]
            sage: Gamma0(11).z_integrality_algorithm(indata)
            [0, -X0 + X1, 0, 0, 0, 0, 0, 0, 0, 0, -X4 + X5 + X6 - X7 - X8 + X10, -X4 + X5 + X6 - X7 - X9 + X11]

        Explicitly, one could take the module of homogeneous
        polynomials in two variables in degree k - 1 over the integers
        (or any ring), and use the integrality algorithm to do the
        bulk of the reduction work.  (There are several caveats, of
        course. The considerations for nu2 or nu3 not trivial become
        not trivial, the action of SL2Z on that module is not trivial
        and has to be taken into account, and last but not least we
        must add back resp. take care in a later step for the
        relations coming from [0, 1], which cannot be thrown away
        anymore. I lack data, and code, but maybe for higher weigths
        an index may creep in, in that final step.)

        The claims from above follow by looking at the matrix of the
        basis transformation form the original generators of RAW(G) to
        the twelve generators one gets by replacing the "least index"
        cosets by "their" relation. This matrix is triangular (by
        construction, due to the least index property), and the
        diagonal elements are all ones.  This proves at one and the
        same time that REL(G) is saturated in RAW(G) and thus
        RAW(G)/REL(G) as a ZZ-module is torsion-free, and that the
        "left over" generators which aren't some "least index" form
        a basis of this factor module, i.e. the submodule they generate in the
        factor module is of trivial index one there.

        For how to obtain such enumerations (having the "least index"
        property), a "geometric" point of view is helpful. Either one
        has in mind the graph with the 2- and 3-relations as vertices,
        or else, dually, the graph with the cusps as vertices. In the
        first case, the resulting algorithm is more like the venerable
        group theoretic Todd-Coxeter method for coset enumaration, in
        the second case, more like the Farey symbol driven
        construction of a fundamental domain.  But since these are
        dual notions, the difference does not really matter.  (The
        output and the intermediate data of the method
        "form_list_of_cusps" in the file
        "shp_package/moments/funddomain.magma" from Robert Pollacks
        code is easily recognized to be a Farey symbol, a web search
        brought up
        e.g. 'http://www.cs.st-andrews.ac.uk/~alexk/papers/Aachen06.pdf',
        which is both introductory, and has some references to
        literature on those.  The timings for computing the data of
        Gamma(2^n) for small n on page 17 there seem to be
        outdated. The code below computes coset representatives ::

            sage: reps, gens, basis = Gamma(64).z_sl2z_todd_coxeter()
            sage: ( len(reps), len(gens), len(basis) )
            (196608, 16385, 16385)

            sage: reps, gens, basis = Gamma(128).z_sl2z_todd_coxeter()#long time
            sage: ( len(reps), len(gens), len(basis) )                #long time
            (1572864, 131073, 131073)

        and generators for Gamma(64) in about 17 seconds on my 2.0 GHz
        Core2Duo system, for Gamma(128) in about 478 seconds.)

        One last note about the 'binary tree heuristic'. It essentially says
        that the following matrices have only O( G.index() * log(G.index()) )
        non-zero entries, i.e. are rather sparse, and they also are 'as sparse
        as possible'::

            sage: nrows, ncolumns, matrixdict = Gamma0(20011).z_wf_manin_gens_to_basis_dict()
            sage: nrows, ncolumns
            (20012, 3335)
            sage: len(matrixdict)
            155154
            sage: 1.0 * 155154 / (20012 * 3335)
            0.00232475207529156

        (Even the mere creation of the corresponding matrix::

            sage: Matrix(ZZ, nrows, ncolumns, matrixdict)
            20012 x 3335 sparse matrix over Integer Ring

        is an order of magnitude slower than the computations leading to it.)

        Compare this density with (which is not really bad either!)::

            sage: m = ModularSymbols(20011)                   #long time
            sage: m.manin_gens_to_basis()                     #long time
            20012 x 3335 sparse matrix over Rational Field
            sage: m.manin_gens_to_basis().density() * 1.0     #long time
            0.00558522457739749

        This has consequences for the density of matrices of the Hecke
        operators expressed in the different manin bases::

            sage: G = Gamma0(840)
            sage: MS = G.modular_symbols()
            sage: MS.manin_gens_to_basis().density() * 1.0
            0.0311846139971140
            sage: wf_reps, wf_gens, wf_manin_basis = G.z_sl2z_todd_coxeter()
            sage: G.z_wf_manin_gens_to_basis().density() * 1.0
            0.0146915584415584
            sage: temp = []
            sage: for i in wf_manin_basis:
            ....:   temp.append(MS.manin_symbol( [ wf_reps[i].c(), wf_reps[i].d() ] ).list())
            sage: Mat = matrix(temp)
            sage: Mat.det()
            1
            sage: InvMat = ~Mat
            sage: for p in prime_range(2, 41):
            ....:   [ MS.hecke_matrix(p).density()*1.0, (Mat*MS.hecke_matrix(p)*InvMat).density()*1.0 ]
            [0.0446550851745657, 0.0273840445269017]
            [0.0785764884466183, 0.0549030190588632]
            [0.158401079440040, 0.101055827289594]
            [0.212460785967279, 0.158616967448136]
            [0.393826952268511, 0.278664192949907]
            [0.461008601787823, 0.335800303592511]
            [0.527502108281329, 0.415105414066453]
            [0.552814977230562, 0.448669252825097]
            [0.622310676336650, 0.510035419126328]
            [0.662013830325519, 0.573283858998145]
            [0.687637038286389, 0.590561646146062]
            [0.726017878225670, 0.646071850227694]
        """
        return 42

    def z_integrality_algorithm(self, indata):
        wf_list = self.z_wellformed_coset_reps_and_relations()
        if len(list(indata)) != len(wf_list):
            raise TypeError("indata does not have the correct length")
        wf_relations = self.z_wellformed_all_relations_but_one()
        for relation in wf_relations:
            index = relation[0]
            if index != relation[1]:         # not a degenerated relation
                indata[relation[1]] -= indata[index]
                if len(relation) == 3:
                    indata[relation[2]] -= indata[index]
            indata[index] -= indata[index]   # now clean the data of this index
        return indata

    def z_wellformed_SZ_graph(self):
        from sage.all import Graph
        graph = Graph(multiedges=True)
        #vertices == orbits of i and zeta_6, on the Riemann surface of self
        #edges == orbits of the geodesic path joining i and zeta_6
        #
        # for a 2-relation [A, AS], both A and AS map i to the same point
        # for a 3-relation [A, AZ, AZZ], all three map zeta_6 to the same point
        firstrel = self.z_wellformed_2_relations()[0]
        wf_list = self.z_wellformed_coset_reps_and_relations()
        sublist = wf_list[0]
        graph.add_edge(str(firstrel), str(sublist[3]),
                       label = str(sublist[0]))
        for index in xrange(1, len(wf_list)):
            sublist = wf_list[index]
            thisrel = sublist[3]
            if len(thisrel) == 1:
                thisrel = [ thisrel[0], index ]
            graph.add_edge(str(wf_list[ sublist[4] ][3]), str(thisrel),
                label = sublist[0])
        return graph

    def z_wellformed_cusps_graph(self):
        from sage.all import Graph
        graph = Graph(multiedges=True)
        #vertices == cusps
        #edges == unimodular paths
        return graph

    def z_wellformed_coset_reps_and_relations(self):
        r"""
        If G is odd, we go over to the even subgroup containing it, and with
        the same projective index +-G, generated by G and SL2Z([-1, 0, 0, -1]),
        i.e. +-G == G union -G.

        The return value is a list of quintuples, of length of the index of
        +-G in SL2Z. (One have to bare this in mind working with odd groups.)
        The five entries of the sublists have the following meaning:

        1. a (right) coset representative for +-G in SL2Z, as matrix in SL2Z

        2. the same matrix written as a word in S and Z (for Id the empty word)

        3. the "round" in which this rep was found, and at the same the length
           (with a suitable definition of "length") of a shortest possible
           word of a representative of this coset.
           In most cases this is just the length of the entry before, if not,
           then it's the length of the word of predecessor, with an "S" added,
           and the entry before is of length exactly one longer than that.
           (The "overwriting" case in the code, which is done because it is
           both customary and practical to always have the three elements
           of a 3-relation "together" --- if nu3 is zero, this e.g. gives a
           fundamental domain whose boundary consists of unimodular paths only,
           and which is made up of copies of the "bigger" triangles 0, 1, oo.)
           For the definition of "length":
            - each consecutive sequence of "S" counts as one digit,
              e.g "S", "SS" and "SSS" all have length 1
            - each consecutive sequence of "R" counts as one digit
              e.g."SZ", "SSZ", "SSSZ", "SZZ", "SSZZ" and "SSSZZ" have length 2
            - the empty word has length 1 (!!), and counts before a "Z" at the
              front, i.e. "Z" and "ZZ" both have length 2, "ZS", "ZZS" length 3

        4. If this index is a "least" index, then the relation of which it is
           the first entry.
           Degenerated relations are written as [a, a] resp. [a, a, a].
           If this entry is not a "least" index, i.e. this index belongs to
           the (manin) basis, then this entry is a list of length one,
           containing the index of the entry, which is not the predecessor
           of this coset here, but which also has this index as second
           member of its relation (only the case of 2-relations occur, because
           the 3-relations are always "kept together".)

        5. For the zeroeth coset (with "Id" as representative), "None".
           For all others, the predecessor coset.
           This can be gotten by chopping off the last digit of the word
           (either "S", or "Z" resp. "ZZ", if the word does not end with "S").
           Alternatively, it is the least index having this index as a
           second (or third) entry in its relation.
           (Thus every start of a word of a coset representative is itself the
           the word of another, earlier, coset representative in this list.)

        EXAMPLE::

            sage: Gamma0(11).z_wellformed_coset_reps_and_relations()
            [[[1 0]
            [0 1], '', 1, [0, 2, 1], None], [[-1  1]
            [-1  0], 'ZZ', 1, [0], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 3], 0], [[-1  0]
            [-1 -1], 'ZS', 3, [3, 4, 5], 2], [[ 0  1]
            [-1  2], 'ZSZ', 4, [4, 6], 3], [[ 1 -1]
            [ 2 -1], 'ZSZZ', 4, [5, 7], 3], [[1 0]
            [2 1], 'ZSZS', 5, [6, 8, 9], 4], [[-1 -1]
            [-1 -2], 'ZSZZS', 5, [7, 10, 11], 5], [[ 0 -1]
            [ 1 -3], 'ZSZSZ', 6, [8, 10], 6], [[-1  1]
            [-3  2], 'ZSZSZZ', 6, [9, 11], 6], [[-1  2]
            [-2  3], 'ZSZZSZ', 6, [8], 7], [[ 2 -1]
            [ 3 -1], 'ZSZZSZZ', 6, [9], 7]]
        """
        import sage.modular.arithgroup.all as arithgroup
        if arithgroup.is_Gamma0(self):
            ret_list, _, __ = self.z_wellformed_coset_reps_and_relations_Gamma0()
        elif arithgroup.is_Gamma1(self):
            ret_list = self.z_wellformed_coset_reps_and_relations_Gamma1()
        elif arithgroup.is_Gamma(self):
            ret_list = self.z_wellformed_coset_reps_and_relations_Gamma()
        else:
            ret_list = self.z_wellformed_coset_reps_and_relations_generic()
        return ret_list

    @cached_method
    def z_wellformed_coset_reps_and_relations_generic(self):
        r"""

        EXAMPLE::

            sage: w = SymmetricGroup(2)([2,1])
            sage: G = ArithmeticSubgroup_Permutation(w, w)
            sage: G.z_wellformed_coset_reps_and_relations()
            [[[1 0]
            [0 1], '', 1, [0, 0, 0], None], [[ 0 -1]
            [ 1  0], 'S', 1, [1, 1, 1], 0]]

            sage: a2 = SymmetricGroup(7)([(1,2),(3,4),(5,6)])
            sage: a3 = SymmetricGroup(7)([(1,3,5),(2,6,7)])
            sage: G2 = ArithmeticSubgroup_Permutation(a2*a3, ~a2 * ~a3)
            sage: G2.z_wellformed_coset_reps_and_relations()
            [[[1 0]
            [0 1], '', 1, [0, 2, 3], None], [[ 0 -1]
            [ 1  0], 'S', 1, [1, 4, 5], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 6], 0], [[-1  1]
            [-1  0], 'ZZ', 2, [3, 4], 0], [[-1  1]
            [ 0 -1], 'SZ', 2, [3], 1], [[ 1  0]
            [-1  1], 'SZZ', 2, [5, 5], 1], [[-1  0]
            [-1 -1], 'ZS', 3, [6, 6, 6], 2]]

            sage: import sage.modular.arithgroup.arithgroup_perm as ap
            sage: ap.HsuExample10().z_wellformed_coset_reps_and_relations()
            [[[1 0]
            [0 1], '', 1, [0, 2, 3], None], [[ 0 -1]
            [ 1  0], 'S', 1, [1, 4, 5], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 6], 0], [[-1  1]
            [-1  0], 'ZZ', 2, [3, 4], 0], [[-1  1]
            [ 0 -1], 'SZ', 2, [3], 1], [[ 1  0]
            [-1  1], 'SZZ', 2, [5, 7], 1], [[-1  0]
            [-1 -1], 'ZS', 3, [6, 8, 7], 2], [[ 1 -1]
            [ 2 -1], 'ZSZZ', 3, [5], 6], [[ 0  1]
            [-1  2], 'ZSZ', 4, [8, 9], 6], [[1 0]
            [2 1], 'ZSZS', 5, [9, 9, 9], 8]]
            sage: ap.HsuExample18().z_wellformed_coset_reps_and_relations()
            [[[1 0]
            [0 1], '', 1, [0, 2, 3], None], [[ 0 -1]
            [ 1  0], 'S', 1, [1, 4, 5], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 6], 0], [[-1  1]
            [-1  0], 'ZZ', 2, [3, 4], 0], [[-1  1]
            [ 0 -1], 'SZ', 2, [3], 1], [[ 1  0]
            [-1  1], 'SZZ', 2, [5, 7], 1], [[-1  0]
            [-1 -1], 'ZS', 3, [6, 7, 8], 2], [[ 0  1]
            [-1  2], 'ZSZ', 3, [5], 6], [[ 1 -1]
            [ 2 -1], 'ZSZZ', 4, [8, 9], 6], [[-1 -1]
            [-1 -2], 'ZSZZS', 5, [9, 10, 11], 8], [[-1  2]
            [-2  3], 'ZSZZSZ', 6, [10, 12], 9], [[ 2 -1]
            [ 3 -1], 'ZSZZSZZ', 6, [11, 13], 9], [[2 1]
            [3 2], 'ZSZZSZS', 7, [12, 14, 15], 10], [[-1 -2]
            [-1 -3], 'ZSZZSZZS', 7, [13, 16, 17], 11], [[ 1 -3]
            [ 2 -5], 'ZSZZSZSZ', 8, [14, 17], 12], [[-3  2]
            [-5  3], 'ZSZZSZSZZ', 8, [15, 16], 12], [[-2  3]
            [-3  4], 'ZSZZSZZSZ', 8, [15], 13], [[ 3 -1]
            [ 4 -1], 'ZSZZSZZSZZ', 8, [14], 13]]
        """
        from all import SL2Z

        Id = SL2Z.one()

        minusId = SL2Z([-1, 0,     #   minusId == -Id  and  minusId^2 == Id
                         0,-1])

        S       = SL2Z([ 0,-1,     #   S^2 == minusId  and  S^4 == Id
                         1, 0])    #   S fixes i in the upper half plane

        Z       = SL2Z([ 0,-1,     #   Z^3 == Id
                         1,-1])    #   Z fixes zeta_6 in the upper half plane

        ret_list   = []     #list of sublists (currently empty) to return
        listlength = 0      #helper variable,  listlength == len(ret_list)
        start      = 0      #first of the sublists with relations to be filled
        ret_list.append([Id,       #coset representative as matrix in SL2Z
                         "",       #and as word in S and R (empty word for Id)
                         1,        #Id counts as one "digit", alone or before Z
                         [],       #relation beginning with Id (filled later)
                         None ])   #(index of) predecessor (None for Id)
        listlength += 1            #ret_list now has one more sublist
        if not S in self:
            ret_list.append([S,    #possibly later overwritten by Z or Z*Z
                             "S",  #possibly later overwritten by "Z" or "ZZ"
                             1,    #S also counts as one "digit"
                             [],   #to be filled next round
                             0 ])  #the predecesor of S is Id (and has index 0)
            listlength += 1        #ret_list now has one more sublist
        #  Remark
        #Now, the relations (sublist[3]) still have to be filled in.
        #In the course of doing so, ret_list may grow further.
        #The relations belonging to such new sublists can't be filled in
        #immediately, leading to another round to follow after the current one.
        #The process stops after the first round, in which no new sublist
        #had to be added.
        #The rounds handle alternatingly either a "row" of sublists where
        #the words (sublist[1]) end with either Z or ZZ, and one "S" is to be
        #added; or else a "row" of sublists, where to the words one or two "Z"
        #are to be added.
        while start < listlength:
            cur_wordlength = ret_list[start][2]   #length of the word in digits
            cur_row = range(start, listlength)    #may later become "holes"
            start = listlength                    #(remember for next round)
            if cur_wordlength%2 == 0: #all of the words end with either Z or ZZ
                for a in cur_row:
                    a_sublist = ret_list[a]
                    A = a_sublist[0]
                    AS = A * S
                    if AS * ~A in self:    #no need to test for the negative
                        a_sublist[3] = [a, a]    #degenerated 2-relation
                    else:
                        #  Remark
                        #The following "inner" loop essentially makes this
                        #whole algorithm quadratic in the index of the group
                        #self in SL2Z (i.e. in the number of cosets).
                        #If at the start, there is some list of the coset reps
                        #already available and a kind of lookup function, that
                        #gives the index of the coset an arbitrary matrix
                        #belongs to, then the inner loop(s) can be avoided, and
                        #replaced by just a lookup in a table of coset reps
                        #"already visited" so far, thus making this step O(1).
                        #And making the entire algorithm linear in the index
                        #(which is the best possible).
                        #This is the case having M-symbols at hand, e.g. for
                        #Gamma0(N), Gamma1(N), Gamma(N), see below.
                        for b in cur_row[ cur_row.index(a) + 1 :  ]:
                            b_sublist = ret_list[b]
                            B = b_sublist[0]
                            ASinvB = AS * ~B
                            if (ASinvB in self) or (minusId * ASinvB in self):
                                a_sublist[3] = [a, b]         #2-relation
                                b_sublist[3] = [a]
                                cur_row.remove(b)  #remove from outer for-loop
                                break              #leave inner for-loop
                        if len(a_sublist[3]) != 2: #results in new coset rep
                            a_sublist[3] = [a, listlength]    #2-relation
                            ret_list.append([AS,
                                             a_sublist[1] + "S",
                                             a_sublist[2] + 1,
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            listlength += 1
            else:  #cur_wordlength%2 != 0 , the words end with neither Z nor ZZ
                for a in cur_row:
                    a_sublist = ret_list[a]
                    A = a_sublist[0]
                    AZ = A * Z
                    if AZ * ~A in self:    #no need to test for the negative
                        a_sublist[3] = [a, a, a]    #degenerated 3-relation
                    else:
                        for b in cur_row[ cur_row.index(a) + 1 :  ]:
                            b_sublist = ret_list[b]
                            B = b_sublist[0]
                            AZinvB = AZ * ~B
                            if (AZinvB in self) or (minusId * AZinvB in self):
                                a_sublist[3] = [a, b]   #incomplete 3-relation
                                b_sublist[0] = AZ                   #overwrite
                                b_sublist[1] = a_sublist[1] + "Z"   #overwrite
                                #b_sublist[2]==a_sublist[2]  #don't overwrite!
                                b_sublist[3] = [ b_sublist[4] ]     #old pred.
                                b_sublist[4] = a                    #overwrite
                                cur_row.remove(b)  #remove from outer for-loop
                                break              #leave inner for-loop
                        if len(a_sublist[3]) != 2: #results in new coset rep
                            a_sublist[3] = [a, listlength]   #incomplete 3-rel.
                            ret_list.append([AZ,
                                             a_sublist[1] + "Z",
                                             a_sublist[2] + 1,
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            listlength += 1
                        AZZ = AZ * Z
                        for b in cur_row[ cur_row.index(a) + 1 :  ]:
                            b_sublist = ret_list[b]
                            B = b_sublist[0]
                            AZZinvB = AZZ * ~B
                            if (AZZinvB in self) or (minusId*AZZinvB in self):
                                a_sublist[3] = a_sublist[3] + [b] #comp. 3-rel.
                                b_sublist[0] = AZZ                  #overwrite
                                b_sublist[1] = a_sublist[1] + "ZZ"  #overwrite
                                #b_sublist[2]==a_sublist[2]  #don't overwrite!
                                b_sublist[3] = [ b_sublist[4] ]     #old pred.
                                b_sublist[4] = a                    #overwrite
                                cur_row.remove(b)  #remove from outer for-loop
                                break              #leave inner for-loop
                        if len(a_sublist[3]) != 3: #results in new coset rep
                            a_sublist[3] = a_sublist[3] + [listlength]  #3-rel.
                            ret_list.append([AZZ,
                                             a_sublist[1] + "ZZ", #ZZ counts as
                                             a_sublist[2] + 1, #only one digit!
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            listlength += 1
        return ret_list

    @cached_method
    def z_wellformed_coset_reps_and_relations_Gamma1(self):
        r"""

        EXAMPLE::

            sage: Gamma1(1).z_wellformed_coset_reps_and_relations_Gamma1()
            [[[1 0]
            [0 1], '', 1, [0, 0, 0], None]]

            sage: Gamma1(2).z_wellformed_coset_reps_and_relations_Gamma1()
            [[[1 0]
            [0 1], '', 1, [0, 2, 1], None], [[-1  1]
            [-1  0], 'ZZ', 1, [0], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 2], 0]]

            sage: Gamma1(3).z_wellformed_coset_reps_and_relations_Gamma1()
            [[[1 0]
            [0 1], '', 1, [0, 2, 1], None], [[-1  1]
            [-1  0], 'ZZ', 1, [0], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 3], 0], [[-1  0]
            [-1 -1], 'ZS', 3, [3, 3, 3], 2]]

            sage: Gamma1(4).z_wellformed_coset_reps_and_relations_Gamma1()
            [[[1 0]
            [0 1], '', 1, [0, 2, 1], None], [[-1  1]
            [-1  0], 'ZZ', 1, [0], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 3], 0], [[-1  0]
            [-1 -1], 'ZS', 3, [3, 4, 5], 2], [[ 0  1]
            [-1  2], 'ZSZ', 4, [4, 5], 3], [[ 1 -1]
            [ 2 -1], 'ZSZZ', 4, [4], 3]]

            sage: Gamma1(5).z_wellformed_coset_reps_and_relations_Gamma1()
            [[[1 0]
            [0 1], '', 1, [0, 2, 1], None], [[-1  1]
            [-1  0], 'ZZ', 1, [0], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 3], 0], [[-1  0]
            [-1 -1], 'ZS', 3, [3, 4, 5], 2], [[ 0  1]
            [-1  2], 'ZSZ', 4, [4, 6], 3], [[ 1 -1]
            [ 2 -1], 'ZSZZ', 4, [5, 7], 3], [[1 0]
            [2 1], 'ZSZS', 5, [6, 7, 8], 4], [[ 0 -1]
            [ 1 -3], 'ZSZSZ', 5, [5], 6], [[-1  1]
            [-3  2], 'ZSZSZZ', 6, [8, 9], 6], [[1 1]
            [2 3], 'ZSZSZZS', 7, [9, 10, 11], 8], [[ 1 -2]
            [ 3 -5], 'ZSZSZZSZ', 8, [10, 11], 9], [[-2  1]
            [-5  2], 'ZSZSZZSZZ', 8, [10], 9]]

        """
        import sage.modular.arithgroup.all as arithgroup
        if not arithgroup.is_Gamma1(self):
            raise TypeError, "... only for Gamma1(N) ..."

        N = self.level()
        checklist = [-1 for x in xrange(N*N)] #crude checklist for M-Symbols

        from all import SL2Z

        Id      = SL2Z([ 1, 0,
                         0, 1])

        minusId = SL2Z([-1, 0,     #   minusId == -Id  and  minusId^2 == Id
                         0,-1])

        S       = SL2Z([ 0,-1,     #   S^2 == minusId  and  S^4 == Id
                         1, 0])    #   S fixes i in the upper half plane

        Z       = SL2Z([ 0,-1,     #   Z^3 == Id
                         1,-1])    #   Z fixes zeta_6 in the upper half plane

        ret_list   = []     #list of sublists (currently empty) to return
        listlength = 0      #helper variable,  listlength == len(ret_list)
        start      = 0      #first of the sublists with relations to be filled
        ret_list.append([Id,       #coset representative as matrix in SL2Z
                         "",       #and as word in S and R (empty word for Id)
                         1,        #Id counts as one "digit", alone or before Z
                         [],       #relation beginning with Id (filled later)
                         None ])   #(index of) predecessor (None for Id)
        checklist[N*0 + 1%N] = listlength       #Id.c() == 0    Id.d() == 1
        listlength += 1            #ret_list now has one more sublist
        if not S in self:
            ret_list.append([S,    #possibly later overwritten by Z or Z*Z
                             "S",  #possibly later overwritten by "Z" or "ZZ"
                             1,    #S also counts as one "digit"
                             [],   #to be filled next round
                             0 ])  #the predecesor of S is Id (and has index 0)
            checklist[N*1 + 0] = listlength     #S.c() == 1    S.d() == 0
            listlength += 1        #ret_list now has one more sublist
        while start < listlength:
            cur_wordlength = ret_list[start][2]   #length of the word in digits
            cur_row = range(start, listlength)    #may later become "holes"
            start = listlength                    #(remember for next round)
            if cur_wordlength%2 == 0:
                for a in cur_row:
                    a_sublist = ret_list[a]
                    A = a_sublist[0]
                    AS = A * S
                    c    = AS.c()%N
                    d    = AS.d()%N
                    negc = -c
                    negd = -d
                    if c < 0:
                        c += N
                    if negc < 0:
                        negc += N
                    if d < 0:
                        d += N
                    if negd < 0:
                        negd += N
                    checklist_index = N*c + d
                    b = checklist[ checklist_index ]
                    if b == -1:
                        checklist_index = N*negc + negd
                        b = checklist[ checklist_index ]
                    if a == b:
                        a_sublist[3] = [a, a]    #degenerated 2-relation
                    else:
                        if b != -1:
                            a_sublist[3] = [a, b]         #2-relation
                            b_sublist    = ret_list[b]
                            b_sublist[3] = [a]
                            cur_row.remove(b) #remove from for-loop
                        else: # b == -1  results in new coset rep
                            a_sublist[3] = [a, listlength]    #2-relation
                            ret_list.append([AS,
                                             a_sublist[1] + "S",
                                             a_sublist[2] + 1,
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            checklist[ checklist_index ] = listlength
                            listlength += 1
            else:    #cur_wordlength%2 != 0
                for a in cur_row:
                    a_sublist = ret_list[a]
                    A = a_sublist[0]
                    AZ = A * Z
                    c    = AZ.c()%N
                    d    = AZ.d()%N
                    negc = -c
                    negd = -d
                    if c < 0:
                        c += N
                    if negc < 0:
                        negc += N
                    if d < 0:
                        d += N
                    if negd < 0:
                        negd += N
                    checklist_index = N*c + d
                    b = checklist[ checklist_index ]
                    if b == -1:
                        checklist_index = N*negc + negd
                        b = checklist[ checklist_index ]
                    if a == b:
                        a_sublist[3] = [a, a, a]    #degenerated 3-relation
                    else:
                        if b != -1:
                            a_sublist[3] = [a, b]   #incomplete 3-relation
                            b_sublist    = ret_list[b]
                            b_sublist[0] = AZ                   #overwrite
                            b_sublist[1] = a_sublist[1] + "Z"   #overwrite
                            #b_sublist[2]==a_sublist[2]      #don't overwrite!
                            b_sublist[3] = [ b_sublist[4] ]     #old pred.
                            b_sublist[4] = a                    #overwrite
                            cur_row.remove(b)  #remove from for-loop
                        else: # b == -1  results in new coset rep
                            a_sublist[3] = [a, listlength]   #incomplete 3-rel.
                            ret_list.append([AZ,
                                             a_sublist[1] + "Z",
                                             a_sublist[2] + 1,
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            checklist[ checklist_index ] = listlength
                            listlength += 1
                        AZZ = AZ * Z
                        c    = AZZ.c()%N
                        d    = AZZ.d()%N
                        negc = -c
                        negd = -d
                        if c < 0:
                            c += N
                        if negc < 0:
                            negc += N
                        if d < 0:
                            d += N
                        if negd < 0:
                            negd += N
                        checklist_index = N*c + d
                        b = checklist[ checklist_index ]
                        if b == -1:
                            checklist_index = N*negc + negd
                            b = checklist[ checklist_index ]
                        if b != -1:             #then b != a  by the above
                            a_sublist[3] = a_sublist[3] + [b] #comp. 3-rel.
                            b_sublist    = ret_list[b]
                            b_sublist[0] = AZZ                  #overwrite
                            b_sublist[1] = a_sublist[1] + "ZZ"  #overwrite
                            #b_sublist[2]==a_sublist[2]      #don't overwrite!
                            b_sublist[3] = [ b_sublist[4] ]     #old pred.
                            b_sublist[4] = a                    #overwrite
                            cur_row.remove(b)  #remove from for-loop
                        else: # b == -1  results in new coset rep
                            a_sublist[3] = a_sublist[3] + [listlength]  #3-rel.
                            ret_list.append([AZZ,
                                             a_sublist[1] + "ZZ",
                                             a_sublist[2] + 1, #only one digit!
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            checklist[ checklist_index ] = listlength
                            listlength += 1

        return ret_list

    @cached_method
    def z_wellformed_coset_reps_and_relations_Gamma0(self):
        r"""

        EXAMPLE::

            sage: Gamma0(1).z_wellformed_coset_reps_and_relations_Gamma0()
            ([[[1 0]
            [0 1], '', 1, [0, 0, 0], None]], [0], [0])

            sage: Gamma0(2).z_wellformed_coset_reps_and_relations_Gamma0()
            ([[[1 0]
            [0 1], '', 1, [0, 2, 1], None], [[-1  1]
            [-1  0], 'ZZ', 1, [0], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 2], 0]], [0, 1, 2], [0, 1, 2])

            sage: Gamma0(3).z_wellformed_coset_reps_and_relations_Gamma0()
            ([[[1 0]
            [0 1], '', 1, [0, 2, 1], None], [[-1  1]
            [-1  0], 'ZZ', 1, [0], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 3], 0], [[-1  0]
            [-1 -1], 'ZS', 3, [3, 3, 3], 2]], [0, 1, 3, 2], [0, 1, 3, 2])

            sage: Gamma0(4).z_wellformed_coset_reps_and_relations_Gamma0()
            ([[[1 0]
            [0 1], '', 1, [0, 2, 1], None], [[-1  1]
            [-1  0], 'ZZ', 1, [0], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 3], 0], [[-1  0]
            [-1 -1], 'ZS', 3, [3, 4, 5], 2], [[ 0  1]
            [-1  2], 'ZSZ', 4, [4, 5], 3], [[ 1 -1]
            [ 2 -1], 'ZSZZ', 4, [4], 3]], [0, 1, 3, 4, 2, 5], [0, 1, 4, 2, 3, 5])

            sage: Gamma0(5).z_wellformed_coset_reps_and_relations_Gamma0()
            ([[[1 0]
            [0 1], '', 1, [0, 2, 1], None], [[-1  1]
            [-1  0], 'ZZ', 1, [0], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 3], 0], [[-1  0]
            [-1 -1], 'ZS', 3, [3, 4, 5], 2], [[ 0  1]
            [-1  2], 'ZSZ', 4, [4, 4], 3], [[ 1 -1]
            [ 2 -1], 'ZSZZ', 4, [5, 5], 3]], [0, 1, 3, 5, 4, 2], [0, 1, 5, 2, 4, 3])

            sage: Gamma0(11).z_wellformed_coset_reps_and_relations_Gamma0()
            ([[[1 0]
            [0 1], '', 1, [0, 2, 1], None], [[-1  1]
            [-1  0], 'ZZ', 1, [0], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 3], 0], [[-1  0]
            [-1 -1], 'ZS', 3, [3, 4, 5], 2], [[ 0  1]
            [-1  2], 'ZSZ', 4, [4, 6], 3], [[ 1 -1]
            [ 2 -1], 'ZSZZ', 4, [5, 7], 3], [[1 0]
            [2 1], 'ZSZS', 5, [6, 8, 9], 4], [[-1 -1]
            [-1 -2], 'ZSZZS', 5, [7, 10, 11], 5], [[ 0 -1]
            [ 1 -3], 'ZSZSZ', 6, [8, 10], 6], [[-1  1]
            [-3  2], 'ZSZSZZ', 6, [9, 11], 6], [[-1  2]
            [-2  3], 'ZSZZSZ', 6, [8], 7], [[ 2 -1]
            [ 3 -1], 'ZSZZSZZ', 6, [9], 7]], [0, 1, 3, 7, 9, 10, 5, 6, 11, 8, 4, 2], [0, 1, 11, 2, 10, 6, 7, 3, 9, 4, 5, 8])

        """
        import sage.modular.arithgroup.all as arithgroup
        if not arithgroup.is_Gamma0(self):
            raise TypeError, "... only for Gamma0(N) ..."

        N = self.level()
        import sage.modular.modsym.p1list as p1list
        P1N = p1list.P1List(N)
        lookup_wf  = [-1 for x in xrange(len(P1N))]    #checklist
        lookup_P1N = [ 0 for x in xrange(len(P1N))]

        from all import SL2Z

        Id      = SL2Z([ 1, 0,
                         0, 1])

        minusId = SL2Z([-1, 0,     #   minusId == -Id  and  minusId^2 == Id
                         0,-1])

        S       = SL2Z([ 0,-1,     #   S^2 == minusId  and  S^4 == Id
                         1, 0])    #   S fixes i in the upper half plane

        Z       = SL2Z([ 0,-1,     #   Z^3 == Id
                         1,-1])    #   Z fixes zeta_6 in the upper half plane

        ret_list   = []     #list of sublists (currently empty) to return
        listlength = 0      #helper variable,  listlength == len(ret_list)
        start      = 0      #first of the sublists with relations to be filled
        ret_list.append([Id,       #coset representative as matrix in SL2Z
                         "",       #and as word in S and R (empty word for Id)
                         1,        #Id counts as one "digit", alone or before Z
                         [],       #relation beginning with Id (filled later)
                         None ])   #(index of) predecessor (None for Id)
        checklist_index = P1N.index_of_normalized_pair(0, 1%N) #Id.c==0 Id.d==1
        lookup_wf[checklist_index] = listlength
        lookup_P1N[listlength]     = checklist_index
        listlength += 1            #ret_list now has one more sublist
        if not S in self:
            ret_list.append([S,    #possibly later overwritten by Z or Z*Z
                             "S",  #possibly later overwritten by "Z" or "ZZ"
                             1,    #S also counts as one "digit"
                             [],   #to be filled next round
                             0 ])  #the predecesor of S is Id (and has index 0)
            checklist_index            = P1N.index(1, 0)  #S.c()==1  S.d()==0
            lookup_wf[checklist_index] = listlength
            lookup_P1N[listlength]     = checklist_index
            listlength += 1        #ret_list now has one more sublist
        while start < listlength:
            cur_wordlength = ret_list[start][2]   #length of the word in digits
            cur_row = range(start, listlength)    #may later become "holes"
            start = listlength                    #(remember for next round)
            if cur_wordlength%2 == 0:
                for a in cur_row:
                    a_sublist = ret_list[a]
                    A = a_sublist[0]
                    AS = A * S
                    checklist_index = P1N.index(AS.c(), AS.d())
                    b = lookup_wf[ checklist_index ]
                    if a == b:
                        a_sublist[3] = [a, a]    #degenerated 2-relation
                    else:
                        if b != -1:
                            a_sublist[3] = [a, b]         #2-relation
                            b_sublist    = ret_list[b]
                            b_sublist[3] = [a]
                            cur_row.remove(b) #remove from for-loop
                        else: # b == -1  results in new coset rep
                            a_sublist[3] = [a, listlength]    #2-relation
                            ret_list.append([AS,
                                             a_sublist[1] + "S",
                                             a_sublist[2] + 1,
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            lookup_wf[ checklist_index ] = listlength
                            lookup_P1N[ listlength ]     = checklist_index
                            listlength += 1
            else:    #cur_wordlength%2 != 0
                for a in cur_row:
                    a_sublist = ret_list[a]
                    A = a_sublist[0]
                    AZ = A * Z
                    checklist_index = P1N.index(AZ.c(), AZ.d())
                    b = lookup_wf[ checklist_index ]
                    if a == b:
                        a_sublist[3] = [a, a, a]    #degenerated 3-relation
                    else:
                        if b != -1:
                            a_sublist[3] = [a, b]   #incomplete 3-relation
                            b_sublist    = ret_list[b]
                            b_sublist[0] = AZ                   #overwrite
                            b_sublist[1] = a_sublist[1] + "Z"   #overwrite
                            #b_sublist[2]==a_sublist[2]      #don't overwrite!
                            b_sublist[3] = [ b_sublist[4] ]     #old pred.
                            b_sublist[4] = a                    #overwrite
                            cur_row.remove(b)  #remove from for-loop
                        else: # b == -1  results in new coset rep
                            a_sublist[3] = [a, listlength]   #incomplete 3-rel.
                            ret_list.append([AZ,
                                             a_sublist[1] + "Z",
                                             a_sublist[2] + 1,
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            lookup_wf[ checklist_index ] = listlength
                            lookup_P1N[ listlength ]     = checklist_index
                            listlength += 1
                        AZZ = AZ * Z
                        checklist_index = P1N.index(AZZ.c(), AZZ.d())
                        b = lookup_wf[ checklist_index ]
                        if b != -1:             #then b != a  by the above
                            a_sublist[3] = a_sublist[3] + [b] #comp. 3-rel.
                            b_sublist    = ret_list[b]
                            b_sublist[0] = AZZ                  #overwrite
                            b_sublist[1] = a_sublist[1] + "ZZ"  #overwrite
                            #b_sublist[2]==a_sublist[2]      #don't overwrite!
                            b_sublist[3] = [ b_sublist[4] ]     #old pred.
                            b_sublist[4] = a                    #overwrite
                            cur_row.remove(b)  #remove from for-loop
                        else: # b == -1  results in new coset rep
                            a_sublist[3] = a_sublist[3] + [listlength]  #3-rel.
                            ret_list.append([AZZ,
                                             a_sublist[1] + "ZZ",
                                             a_sublist[2] + 1, #only one digit!
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            lookup_wf[ checklist_index ] = listlength
                            lookup_P1N[ listlength ]     = checklist_index
                            listlength += 1

        return ret_list, lookup_wf, lookup_P1N

    @cached_method
    def z_wellformed_coset_reps_and_relations_Gamma(self):
        r"""

        EXAMPLE::

            sage: Gamma(1).z_wellformed_coset_reps_and_relations_Gamma()
            [[[1 0]
            [0 1], '', 1, [0, 0, 0], None]]

            sage: Gamma(2).z_wellformed_coset_reps_and_relations_Gamma()
            [[[1 0]
            [0 1], '', 1, [0, 2, 3], None], [[ 0 -1]
            [ 1  0], 'S', 1, [1, 4, 5], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 5], 0], [[-1  1]
            [-1  0], 'ZZ', 2, [3, 4], 0], [[-1  1]
            [ 0 -1], 'SZ', 2, [3], 1], [[ 1  0]
            [-1  1], 'SZZ', 2, [2], 1]]

            sage: Gamma(3).z_wellformed_coset_reps_and_relations_Gamma()
            [[[1 0]
            [0 1], '', 1, [0, 2, 3], None], [[ 0 -1]
            [ 1  0], 'S', 1, [1, 4, 5], 0], [[ 0 -1]
            [ 1 -1], 'Z', 2, [2, 6], 0], [[-1  1]
            [-1  0], 'ZZ', 2, [3, 7], 0], [[-1  1]
            [ 0 -1], 'SZ', 2, [4, 8], 1], [[ 1  0]
            [-1  1], 'SZZ', 2, [5, 9], 1], [[-1  0]
            [-1 -1], 'ZS', 3, [6, 9, 10], 2], [[1 1]
            [0 1], 'ZZS', 3, [7, 11, 8], 3], [[-2  1]
            [-1  0], 'ZZSZZ', 3, [4], 7], [[ 0  1]
            [-1  2], 'ZSZ', 3, [5], 6], [[ 1 -1]
            [ 2 -1], 'ZSZZ', 4, [10, 11], 6], [[ 1 -2]
            [ 1 -1], 'ZZSZ', 4, [10], 7]]
        """
        import sage.modular.arithgroup.all as arithgroup
        if not arithgroup.is_Gamma(self):
            if not arithgroup.is_SL2Z(self):     #which isn't is_Gamma ... sigh
                raise TypeError, "... only for Gamma(N) ..."

        N = self.level()
        import sage.modular.modsym.p1list as p1list
        P1N = p1list.P1List(N)
        checklist = [-1 for x in xrange(len(P1N) * N * N)]

        from all import SL2Z

        Id = SL2Z.one()

        minusId = SL2Z([-1, 0,     #   minusId == -Id  and  minusId^2 == Id
                         0,-1])

        S       = SL2Z([ 0,-1,     #   S^2 == minusId  and  S^4 == Id
                         1, 0])    #   S fixes i in the upper half plane

        Z       = SL2Z([ 0,-1,     #   Z^3 == Id
                         1,-1])    #   Z fixes zeta_6 in the upper half plane

        ret_list   = []     #list of sublists (currently empty) to return
        listlength = 0      #helper variable,  listlength == len(ret_list)
        start      = 0      #first of the sublists with relations to be filled
        ret_list.append([Id,       #coset representative as matrix in SL2Z
                         "",       #and as word in S and R (empty word for Id)
                         1,        #Id counts as one "digit", alone or before Z
                         [],       #relation beginning with Id (filled later)
                         None ])   #(index of) predecessor (None for Id)
        checklist[ N*N*P1N.index_of_normalized_pair(1%N, 0) + N*0 + 1%N ] =  \
                                                        listlength #Id.a,b,c,d
        listlength += 1            #ret_list now has one more sublist
        if not S in self:
            ret_list.append([S,    #possibly later overwritten by Z or Z*Z
                             "S",  #possibly later overwritten by "Z" or "ZZ"
                             1,    #S also counts as one "digit"
                             [],   #to be filled next round
                             0 ])  #the predecesor of S is Id (and has index 0)
            checklist[ N*N*P1N.index(0, -1) + N*1 + 0 ] = listlength #S.a,b,c,d
            listlength += 1        #ret_list now has one more sublist
        while start < listlength:
            cur_wordlength = ret_list[start][2]   #length of the word in digits
            cur_row = range(start, listlength)    #may later become "holes"
            start = listlength                    #(remember for next round)
            if cur_wordlength%2 == 0:
                for a in cur_row:
                    a_sublist = ret_list[a]
                    A = a_sublist[0]
                    AS = A * S
                    upper_index = P1N.index(AS.a(), AS.b())
                    c    = AS.c()%N
                    d    = AS.d()%N
                    negc = -c
                    negd = -d
                    if c < 0:
                        c += N
                    if negc < 0:
                        negc += N
                    if d < 0:
                        d += N
                    if negd < 0:
                        negd += N
                    checklist_index = N*N*upper_index + N*c + d
                    b = checklist[ checklist_index ]
                    if b == -1:
                        checklist_index = N*N*upper_index + N*negc + negd
                        b = checklist[ checklist_index ]
                    if a == b:
                        a_sublist[3] = [a, a]    #degenerated 2-relation
                    else:
                        if b != -1:
                            a_sublist[3] = [a, b]         #2-relation
                            b_sublist    = ret_list[b]
                            b_sublist[3] = [a]
                            cur_row.remove(b) #remove from for-loop
                        else: # b == -1  results in new coset rep
                            a_sublist[3] = [a, listlength]    #2-relation
                            ret_list.append([AS,
                                             a_sublist[1] + "S",
                                             a_sublist[2] + 1,
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            checklist[ checklist_index ] = listlength
                            listlength += 1
            else:    #cur_wordlength%2 != 0
                for a in cur_row:
                    a_sublist = ret_list[a]
                    A = a_sublist[0]
                    AZ = A * Z
                    upper_index = P1N.index(AZ.a(), AZ.b())
                    c    = AZ.c()%N
                    d    = AZ.d()%N
                    negc = -c
                    negd = -d
                    if c < 0:
                        c += N
                    if negc < 0:
                        negc += N
                    if d < 0:
                        d += N
                    if negd < 0:
                        negd += N
                    checklist_index = N*N*upper_index + N*c + d
                    b = checklist[ checklist_index ]
                    if b == -1:
                        checklist_index = N*N*upper_index + N*negc + negd
                        b = checklist[ checklist_index ]
                    if a == b:
                        a_sublist[3] = [a, a, a]    #degenerated 3-relation
                    else:
                        if b != -1:
                            a_sublist[3] = [a, b]   #incomplete 3-relation
                            b_sublist    = ret_list[b]
                            b_sublist[0] = AZ                   #overwrite
                            b_sublist[1] = a_sublist[1] + "Z"   #overwrite
                            #b_sublist[2]==a_sublist[2]      #don't overwrite!
                            b_sublist[3] = [ b_sublist[4] ]     #old pred.
                            b_sublist[4] = a                    #overwrite
                            cur_row.remove(b)  #remove from for-loop
                        else: # b == -1  results in new coset rep
                            a_sublist[3] = [a, listlength]   #incomplete 3-rel.
                            ret_list.append([AZ,
                                             a_sublist[1] + "Z",
                                             a_sublist[2] + 1,
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            checklist[ checklist_index ] = listlength
                            listlength += 1
                        AZZ = AZ * Z
                        upper_index = P1N.index(AZZ.a(), AZZ.b())
                        c    = AZZ.c()%N
                        d    = AZZ.d()%N
                        negc = -c
                        negd = -d
                        if c < 0:
                            c += N
                        if negc < 0:
                            negc += N
                        if d < 0:
                            d += N
                        if negd < 0:
                            negd += N
                        checklist_index = N*N*upper_index + N*c + d
                        b = checklist[ checklist_index ]
                        if b == -1:
                            checklist_index = N*N*upper_index + N*negc + negd
                            b = checklist[ checklist_index ]
                        if b != -1:             #then b != a  by the above
                            a_sublist[3] = a_sublist[3] + [b] #comp. 3-rel.
                            b_sublist    = ret_list[b]
                            b_sublist[0] = AZZ                  #overwrite
                            b_sublist[1] = a_sublist[1] + "ZZ"  #overwrite
                            #b_sublist[2]==a_sublist[2]      #don't overwrite!
                            b_sublist[3] = [ b_sublist[4] ]     #old pred.
                            b_sublist[4] = a                    #overwrite
                            cur_row.remove(b)  #remove from for-loop
                        else: # b == -1  results in new coset rep
                            a_sublist[3] = a_sublist[3] + [listlength]  #3-rel.
                            ret_list.append([AZZ,
                                             a_sublist[1] + "ZZ",
                                             a_sublist[2] + 1, #only one digit!
                                             [],   #to be filled next round
                                             a ])  #predecessor
                            checklist[ checklist_index ] = listlength
                            listlength += 1

        return ret_list

    def z_wellformed_2_relations(self):
        r"""

        EXAMPLE::

            sage: Gamma0(11).z_wellformed_2_relations()
            [[0, 1], [2, 3], [4, 6], [5, 7], [8, 10], [9, 11]]
        """
        from all import SL2Z
        S = SL2Z([ 0,-1, 1, 0])
        if S in self:
            Srels = [[0, 0]]
        else:
            Srels = [[0, 1]]
        for sublist in self.z_wellformed_coset_reps_and_relations():
            rel = sublist[3]
            if len(rel) == 2:
                Srels.append(rel)
        return Srels

    def z_wellformed_3_relations(self):
        r"""

        EXAMPLE::

            sage: Gamma0(11).z_wellformed_3_relations()
            [[0, 2, 1], [3, 4, 5], [6, 8, 9], [7, 10, 11]]
        """
        Zrels = []
        for sublist in self.z_wellformed_coset_reps_and_relations():
            rel = sublist[3]
            if len(rel) == 3:
                Zrels.append(rel)
        return Zrels

    def z_wellformed_all_relations_but_one(self):
        r"""

        EXAMPLE::

            sage: Gamma0(11).z_wellformed_all_relations_but_one()
            [[0, 2, 1], [2, 3], [3, 4, 5], [4, 6], [5, 7], [6, 8, 9], [7, 10, 11], [8, 10], [9, 11]]
        """
        rels = []
        for sublist in self.z_wellformed_coset_reps_and_relations():
            rel = sublist[3]
            if len(rel) != 1:
                rels.append(rel)
        return rels

    @cached_method
    def z_sl2z_todd_coxeter(self):
        r"""

        EXAMPLE::

            sage: Gamma0(11).z_sl2z_todd_coxeter()
            ([[1 0]
            [0 1], [-1  1]
            [-1  0], [ 0 -1]
            [ 1 -1], [-1  0]
            [-1 -1], [ 0  1]
            [-1  2], [ 1 -1]
            [ 2 -1], [1 0]
            [2 1], [-1 -1]
            [-1 -2], [ 0 -1]
            [ 1 -3], [-1  1]
            [-3  2], [-1  2]
            [-2  3], [ 2 -1]
            [ 3 -1]], [[-1  1]
            [ 0 -1], [ -3   2]
            [-11   7], [ -4   3]
            [-11   8], [-1  0]
            [ 0 -1]], [1, 10, 11])
        """
        from all import SL2Z
        S = SL2Z([ 0,-1, 1, 0])
        Z = SL2Z([ 0,-1, 1, -1])
        minusId = - SL2Z.one()
        reps = []
        gens = []
        if S in self:
            twogens = [S]
        else:
            twogens = []
        threegens = []
        wf_manin_basis = []
        wf_list = self.z_wellformed_coset_reps_and_relations()
        for index in xrange(len(wf_list)):
            sublist = wf_list[index]
            B = sublist[0]
            reps.append(B)
            rel = sublist[3]
            if len(rel) == 1:
                Asublist = wf_list[ rel[0] ]
                A = Asublist[0]
                gens.append(A * S * ~B)
                wf_manin_basis.append(index)
            else:
                if rel[0] == rel[1]:    #degenerated relation
                    if len(rel) == 2:
                        twogens.append(B * S * ~B)
                    else:
                        threegens.append(B * Z * ~B)
        if self.is_odd():
            negreps = []
            for A in reps:
                negreps.append(minusId * A)
            reps += negreps
            for a in xrange(len(gens)):
                if not gens[a] in self:
                    gens[a] *= minusId
        else:    #maybe we have to care for minusId being added to the gens
            if len(twogens) == 0:    #if A = twogens[0] exists, minusId == A^2
                if len(threegens) != 0:    # all threegens are of order 3
                    threegens[0] *= minusId    #now the first has order 6, and
                else:                          #minusId is its third power
                    gens += [minusId]      #often necessary (maybe not always)
        gens = gens + twogens + threegens
        return reps, gens, wf_manin_basis

    def z_sl2z_SZ_word_problem(self, word):
        """
        #The words of the wellformed coset representatives have the property,
        #that all of their start-words are also coset representatives.
        #So among the coset reps which coincide with the start of the input
        #word, there is exactly one with longest length (in "digits").
        #Splitting off to the left (the word of a) generator of self, this
        #start of the input word can then either replaced by a word with fewer
        #digits, or by a word with the same number of digits but such that the
        #"matching" length of the start with some of the coset reps grows.
        #Since the latter possibility is exhausted after a finite number of
        #times (there are only finitely many coset reps, so their length is
        #bounded), ultimately the first possibility must be the case.
        #So splitting off generatos of self to the left, we could reduce the
        #word size (in digits) of the input word.
        #By induction on the length of the input word, we finally get a
        #representation of the input word as a product of generators of self,
        #and one of the wellformed coset representatives.
        #So one can either use this algorithm to get the representation of
        #a group element of self in the generators of the group, or else
        #one can calculate the coset (representative) with respect to self
        #of an arbitrary element of SL2Z.
        #A bit care has to be taken to carry around whether or not (and how
        #many times) we have had to multiply with minusId.
        #(Easy if self is even, but not really difficult if self is odd.)
        """
        list_of_gens = []
        coset_rep = "TODO"
        return list_of_gens, coset_rep

    def z_sl2z_SZ_word_to_matrix(self, word):
        r"""

        EXAMPLE::

            sage: wf_list = Gamma0(20011).z_wellformed_coset_reps_and_relations()
            sage: wf_list[-1]
            [[ 173 -300]
            [ 440 -763], 'ZSZSZZSZSZZSZZSZZSZZSZZSZSZSZSZZSZSZSZZSZ', 34, [20010], 19990]
            sage: Gamma0(20011).z_sl2z_SZ_word_to_matrix(wf_list[-1][1])
            [ 173 -300]
            [ 440 -763]
        """
        from all import SL2Z
        S = SL2Z([ 0,-1, 1, 0])
        Z = SL2Z([ 0,-1, 1, -1])
        ret_mat = SL2Z([ 1,0, 0, 1])
        for i in word:
            if i == "S":
                ret_mat = ret_mat * S
            if i == "Z":
                ret_mat = ret_mat * Z
        return ret_mat

    def z_sl2z_SZ_matrix_to_word(self, matrix):
        """
        #let A be the matix, then compare where zeta_6 lands by the action
        #under AS, AZ, AZZ --- exactly one of the last alternatives maps zeta_6
        #"nearest to i". To check this in Sage, we may and will work in
        #CyclotomicField(3), looking at the values of real and imaginary parts.
        #Knowing the result, we therefore know the last digit of A,
        #can remove that, and restart from the beginning.
        #(Of course the result of the algorithm in "arithgroup_perm" could also
        #be used, translating from a word in Lm and Rm to a word in S and Z.)
        """
        return "TODO"

    def z_farey_symbol_data(self):
        """
        #can be easily read off the output from
        #self.z_wellformed_coset_reps_and_relations()
        #One can then read off generators of the cuspidal subspace from the
        #result. (Maybe all of them, by a pigeonhole-principle like argument??)
        """
        return "TODO"

    def z_wf_manin_gens_to_basis_slow(self):    #soooooo ssslloooooowwww
        r"""

        EXAMPLE::

            sage: Gamma0(11).z_wf_manin_gens_to_basis_slow()
            [-1  0  0]
            [ 1  0  0]
            [ 0  0  0]
            [ 0  0  0]
            [ 0 -1 -1]
            [ 0  1  1]
            [ 0  1  1]
            [ 0 -1 -1]
            [ 0 -1  0]
            [ 0  0 -1]
            [ 0  1  0]
            [ 0  0  1]

        """
        wf_list = self.z_wellformed_coset_reps_and_relations()
        reps, gens, wf_manin_basis = self.z_sl2z_todd_coxeter()
        from sage.rings.all import ZZ
        import sage.matrix.all as mx
        MaSp = mx.MatrixSpace(ZZ,len(wf_list),len(wf_manin_basis),sparse=True)
        ret_matrix = MaSp(0)
        row_list = range(len(wf_list))
        for el_index in xrange(len(wf_manin_basis)):    #basis elements
            row_index = wf_manin_basis[el_index]
            ret_matrix[ row_index, el_index ] = 1
            row_list.remove(row_index)
            Srow_index = wf_list[row_index][3][0]       #S-related to this one
            ret_matrix[ Srow_index, el_index ] = -1
            row_list.remove(Srow_index)
        for row_index in reversed(row_list):    #build "bottom-up"
            rel = wf_list[row_index][3]
            if rel[0] == rel[1]:    #degenerated relation, nothing to do
                pass
            else:       #not a basis element, and we have  row_index == rel[0]
                ret_matrix.set_row_to_multiple_of_row(rel[0], rel[1], -1)
                if len(rel) == 3:
                    ret_matrix.add_multiple_of_row(rel[0], rel[2], -1)
        return ret_matrix

    @cached_method
    def z_wf_manin_gens_to_basis_dict(self):
        r"""

        EXAMPLE::

            sage: nrows, ncolumns, matrixdict = Gamma0(20011).z_wf_manin_gens_to_basis_dict()
            sage: Gamma0(20011).z_wf_manin_gens_to_basis()
            20012 x 3335 sparse matrix over Integer Ring
            sage: len(matrixdict)
            155154
            sage: 1.0 * 155154 / (20012 * 3335)
            0.00232475207529156
        """
        wf_list = self.z_wellformed_coset_reps_and_relations()
        reps, gens, wf_manin_basis = self.z_sl2z_todd_coxeter()
        mdict = dict()
        #fill columns one after another
        for el_index in xrange(len(wf_manin_basis)):    #basis elements
            gen_index = wf_manin_basis[el_index]
            alt_index = wf_list[gen_index][3][0]        #S-related to this one
            mdict[(gen_index, el_index)] =  1
            mdict[(alt_index, el_index)] = -1
            if alt_index == 0:
                continue    #special case ("S" was overwritten by "Z" or "ZZ")
            if wf_list[gen_index][2] > wf_list[alt_index][2]: #they differ by 2
                gen_index = wf_list[gen_index][4]   #go back once
                mdict[(gen_index, el_index)] = -1
                gen_index = wf_list[gen_index][4]   #go back twice
                mdict[(gen_index, el_index)] =  1
            #  Remark
            #By construction, we have from now on that
            # gen_index[2] == alt_index[2]
            #Going up both ways, we must end either "inside" one and the same
            #Z-relation, or else in the S-relation [Id, S] by coming from left
            #and from right.
            #In each of the two cases, the job is done.
            gen_index = wf_list[gen_index][4]       #chop off last Z resp. ZZ
            alt_index = wf_list[alt_index][4]       #chop off last Z resp. ZZ
            while gen_index != alt_index:       #if not, same Z-relation
                mdict[(gen_index, el_index)] = -1   #alternate!
                mdict[(alt_index, el_index)] =  1
                if gen_index == 0 or alt_index == 0:
                    break   #special case {gen, alt} is {0, 1} ([Id, S]) as set
                gen_index = wf_list[gen_index][4]   #chop off last S
                alt_index = wf_list[alt_index][4]   #chop off last S
                mdict[(gen_index, el_index)] =  1
                mdict[(alt_index, el_index)] = -1
                gen_index = wf_list[gen_index][4]   #chop off last Z resp. ZZ
                alt_index = wf_list[alt_index][4]   #chop off last Z resp. ZZ
        return len(wf_list), len(wf_manin_basis), mdict

    @cached_method
    def z_wf_manin_gens_to_basis(self):
        r"""

        EXAMPLE::

            sage: Gamma0(11).z_wf_manin_gens_to_basis()
            [-1  0  0]
            [ 1  0  0]
            [ 0  0  0]
            [ 0  0  0]
            [ 0 -1 -1]
            [ 0  1  1]
            [ 0  1  1]
            [ 0 -1 -1]
            [ 0 -1  0]
            [ 0  0 -1]
            [ 0  1  0]
            [ 0  0  1]
        """
        nrows, ncolumns, matrixdict = self.z_wf_manin_gens_to_basis_dict()
        from sage.rings.all import ZZ
        import sage.matrix.all as mx
        return mx.Matrix(ZZ, nrows, ncolumns, matrixdict)

    @cached_method
    def z_is_normalized_by_J(self):
        r"""

        EXAMPLE::

            sage: w = SymmetricGroup(2)([2,1])
            sage: G1 = ArithmeticSubgroup_Permutation(w, w)
            sage: G1.z_is_normalized_by_J()
            True

            sage: a2 = SymmetricGroup(7)([(1,2),(3,4),(5,6)])
            sage: a3 = SymmetricGroup(7)([(1,3,5),(2,6,7)])
            sage: G2 = ArithmeticSubgroup_Permutation(a2*a3, ~a2 * ~a3);
            sage: G2.z_is_normalized_by_J()
            False
        """
        import sage.modular.arithgroup.all as arithgroup
        if arithgroup.is_Gamma0(self):
            return True
        elif arithgroup.is_Gamma1(self):
            return True
        elif arithgroup.is_Gamma(self):
            return True
        else:
            #      J = GL2Z([-1, 0,     #   J^2 == Id
            #                 0, 1])    #   J not in SL2Z
            #
            #  J * GL2Z([a, b, c, d]) * J   ==  GL2Z([a, -b, -c, d])
            #
            from all import SL2Z
            ret_val = True
            reps, gens, wf_manin_basis = self.z_sl2z_todd_coxeter()
            for g in gens:
                if not SL2Z([g.a(), - g.b(), - g.c(), g.d()]) in self:
                    ret_val = False
                    break
            return ret_val

    @cached_method
    def z_normalized_by_J_symmetry_data(self):
        r"""
        If self is normalized by J = GL2Z([-1, 0, 0, 1]), then the 2-relations,
        the 3-relations, and the z_wellformed_SZ_graph are all symmetric w.r.t.
        to the following involution on the right cosets:

            +-G * A -> +-G * J * A * J * S

        The basis may not always be chosen in a symmetric way, but nevertheless
        one can read off a basis for both H+ and H- in this case; giving rise
        to integrality algorithms as in the "sign 0" case.
        Explicitly, there never is the necessity of a denominator.

        As Hminus base either the part of the Hplus may be taken which is
        not self-symmetric (i.e. has not a 0 in the symmetry_list), or else
        the "rest" of the wf_manin_basis (all the elements which are not in
        wf_manin_basis_hplus).

        EXAMPLE::

            sage: Gamma0(1).z_normalized_by_J_symmetry_data()
            ([], [], [0])

            sage: Gamma0(11).z_normalized_by_J_symmetry_data()
            ([1, 10], [0, -1], [1, 0, 2, 3, 5, 4, 7, 6, 11, 10, 9, 8])
        """
        if not self.z_is_normalized_by_J():
            raise TypeError, "self is not normalized by J"
        wf_list = self.z_wellformed_coset_reps_and_relations()
        symmetry_list = range(len(wf_list))   #preset as if all were self-symm.
        if len(wf_list) > 1:
            if wf_list[1][2] == 1:
                symmetry_list[0] = 1
                symmetry_list[1] = 0
        max_round  = wf_list[-1][2]    #last entry
        listlength = 0
        cur_round  = 1
        while cur_round < max_round: #drop max_round, it brings no new reps
            cur_row = []
            while wf_list[listlength][2] == cur_round:
                relation = wf_list[listlength][3]
                if len(relation) > 1:              #don't append if basis el.
                    if relation[0] != relation[1]: #don't append if degenerated
                        cur_row.append(listlength)
                listlength += 1
            if cur_round%2 == 0:
                for a in cur_row:
                    b = symmetry_list[a]
                    if b != a:                 #symmetric pair
                        a_newrep  = wf_list[a][3][1]
                        if a_newrep >= listlength:  #really some new coset rep.
                            b_newrep = wf_list[b][3][1]    #also new coset rep.
                            symmetry_list[a_newrep] = b_newrep
                            symmetry_list[b_newrep] = a_newrep
                            cur_row.remove(b)
            else:  #cur_round%2 != 0
                for a in cur_row:
                    az_newrep  = wf_list[a][3][1]
                    azz_newrep = wf_list[a][3][2]
                    b = symmetry_list[a]
                    if b == a:
                        if az_newrep >= listlength:  #azz_newrep >= listlength
                            symmetry_list[az_newrep]  = azz_newrep
                            symmetry_list[azz_newrep] = az_newrep
                    elif (b != az_newrep) and (b != azz_newrep):
                        bz_newrep  = wf_list[b][3][1]
                        bzz_newrep = wf_list[b][3][2]
                        if az_newrep >= listlength:  #bzz_newrep >= listlength
                            symmetry_list[az_newrep]  = bzz_newrep
                            symmetry_list[bzz_newrep] = az_newrep
                        if azz_newrep >= listlength: #bz_newrep >= listlength
                            symmetry_list[azz_newrep] = bz_newrep
                            symmetry_list[bz_newrep]  = azz_newrep
                        cur_row.remove(b)
            cur_round += 1
        wf_manin_basis_hplus = []
        sign_list            = []
        cur_wf_list          = range(len(wf_list))    #will get holes
        for a in cur_wf_list:
            a_rel = wf_list[a][3]
            if len(a_rel) == 1:
                wf_manin_basis_hplus.append(a)
                b = symmetry_list[a]    # b != a   by construction
                if a_rel[0] == b:       #may happen for both S- and Z-relations
                    sign_list.append(0)
                else:
                    b_rel = wf_list[b][3]
                    if len(b_rel) == 1: #also happens with both S- and Z-parity
                        sign_list.append(+1)
                        cur_wf_list.remove(b)         #remove from for-loop
                    else:                       #may happen only when adding S
                        sign_list.append(-1)    #for c will hold: len(c_rel)==1
                        c = b_rel[1]            #c == symmetry_list[a_rel[0]]
                        cur_wf_list.remove(c)         #remove from for-loop
        return wf_manin_basis_hplus, sign_list, symmetry_list
