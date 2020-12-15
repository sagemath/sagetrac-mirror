# -*- coding: utf-8 -*-
r"""
Algebraic Tori
==============


Tori are implemented using the equivalence of categories with  character lattices
with an action of the Galois group (at least as large as the Galois group of
a splitting field). For now this Galois group will just be an abstract group,
either a permutation group or a finite matrix group in ``GL(n,ZZ)``.


To define a torus we use ``AlgebraicTorus(character_lattice)``

::

    sage: L = Lattice_ambient([], 1)
    sage: AlgebraicTorus(L)
    Algebraic Torus of rank 1 defined by the following lattice:
    Ambient free module of rank 1 over the principal ideal domain Integer Ring
    and an action by the Galois group of the form:
    Permutation Group with generators [()]

This is the split torus ``\mathbb{G}_m``, with action of the trivial Galois group.

::

    sage: LL = Lattice_ambient(SymmetricGroup(3), 1)
    sage: AlgebraicTorus(LL)
    Algebraic Torus of rank 1 defined by the following lattice:
    Ambient free module of rank 1 over the principal ideal domain Integer Ring
    and an action by the Galois group of the form:
    Symmetric group of order 3! as a permutation group

This is still ``\mathbb{G}_m``, with trivial action of a Galois group isomorphic to `S_3`. Note that
this Galois group is not necessarily the one of a minimal splitting extension.

::

    sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
    sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
    sage: LLL = Lattice_ambient(SymmetricGroup(3), [act1, act2])
    sage: T3 = AlgebraicTorus(LLL); T3
    Algebraic Torus of rank 3 defined by the following lattice:
    Ambient free module of rank 3 over the principal ideal domain Integer Ring
    and an action by the Galois group of the form:
    Symmetric group of order 3! as a permutation group

This is a non-split anisotropic torus with Galois group of splitting field isomorphic to ``S_3``.

::

    sage: SL = SubLattice(L, [2*L.basis()[0]])
    sage: AlgebraicTorus(SL)
    Algebraic Torus of rank 1 defined by the following lattice:
    Free module of degree 1 and rank 1 over Integer Ring
    Echelon basis matrix:
    [2]
    and an action by the Galois group of the form:
    Permutation Group with generators [()]

This torus is obtained from the sublattice of the first lattice `L`. The torus obtained is isomorphic to `L`.


Attributes of a Torus
=====================

- ``torus._lattice`` -- the character lattice of the torus

- ``torus._base_field`` -- optional, the field over which the torus is defined.

- ``torus._splitting_field`` -- optional, a field over which the torus splits


Methods of a Torus
==================

- :meth:`AlgebraicTorus.rank` -- the rank of the torus.

- :meth:`AlgebraicTorus.galois_group` -- the Galois group (as abstract group)
    of a splitting field of the torus.

- :meth:`AlgebraicTorus.character_lattice` -- the character lattice of the torus

- :meth:`AlgebraicTorus.cocharacter_lattice` -- the cocharacter lattice of the torus

- :meth:`AlgebraicTorus.is_rational` -- tests if a point is rational

- :meth:`AlgebraicTorus.Tate_Cohomology` -- the isomorphism
    type of Tate Cohomology groups of the Torus over a local field

- :meth: AlgebraicTorus.product -- the product of the torus with another specified torus
with same base and splitting field.

- :meth:`AlgebraicTorus.restriction_of_scalars` -- returns the torus obtained by
    restriction of scalars

- :meth:`AlgebraicTorus.norm_one_restriction_of_scalars` -- the torus of norm 1
elements in the restriction of scalars
"""

###########################################################################
#       Copyright (C) 2018-2019 Thomas Rüd <tompa.rud@gmail.com>
#                                David Roe <roed.math@gmail.com>

#
#  Distributed under the terms of the GNU General Public License (GPL)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
###########################################################################


from __future__ import print_function, absolute_import

from sage.schemes.generic.scheme import Scheme
from sage.matrix.constructor import matrix
from sage.modules.glattice import GLattice
from sage.rings.rational_field import QQ




def RestrictionOfScalars(nfield, rk = 1):
    """
    Creates the torus defined by restriction of scalars of a split torus to the field of rational numbers.

    INPUT:

    - ``nfield`` -- Number field corresponding to the restriction of scalars.

    - ``rk`` -- rank of the split torus of which we take the restriction of scalars. The default value is 1.

    EXAMPLES::

        sage: from sage.schemes.group_schemes.tori import RestrictionOfScalars
        sage: F.<a> = QuadraticField([2])
        sage: T1 = RestrictionOfScalars(F); T1
        Algebraic torus of rank 2 over Rational Field split by a degree 2 extension
        sage: K.<t> = NumberField(x^4+x^3+x^2+x+1)
        sage: T2 = RestrictionOfScalars(K); T2
        Algebraic torus of rank 4 over Rational Field split by a degree 4 extension
        sage: T2.galois_group()
        Galois group 4T1 (4) with order 4 of x^4 + x^3 + x^2 + x + 1
        sage: T3 = RestrictionOfScalars(K,3); T3
        Algebraic torus of rank 12 over Rational Field split by a degree 4 extension
    """
    L = GLattice(rk).induced_lattice(nfield.galois_group())
    return AlgebraicTorus(L)

def NormOneRestrictionOfScalars(nfield, rk = 1):
    """
    Creates the torus of norm one elements inside the restriction of scalars of a split torus to the field of rational numbers.

    INPUT:

    - ``nfield`` -- Number field corresponding to the restriction of scalars.

    - ``rk`` -- rank of the split torus of which we take the restriction of scalars. The default value is 1.

    EXAMPLES::

        sage: from sage.schemes.group_schemes.tori import NormOneRestrictionOfScalars
        sage: F.<a> = QuadraticField([2])
        sage: T1 = NormOneRestrictionOfScalars(F); T1
        Algebraic torus of rank 1 over Rational Field split by a degree 2 extension
        sage: T1.character_lattice()._action_matrices
        [[-1]]
        sage: T2 = NormOneRestrictionOfScalars(K); T2
        Algebraic torus of rank 3 over Rational Field split by a degree 4 extension
        sage: T3 = NormOneRestrictionOfScalars(K,3); T3
        Algebraic torus of rank 12 over Rational Field split by a degree 4 extension
    """

    L = GLattice(rk).norm_one_restriction_of_scalars(nfield.galois_group())
    return AlgebraicTorus(L)



class AlgebraicTorus(Scheme):
    r"""
    Creates an algebraic torus through its equivalence of categories with the action
    of a Galois group on an integral lattice.

    INPUT:

    - ``lattice`` -- the character lattice with Galois action
      defining the torus.
    """
    def __init__(self, lattice):
        r"""
        Construct an object of the albegraic torus class.

        EXAMPLES::

            sage: L = Lattice_ambient(DihedralGroup(4), 4)
            sage: AlgebraicTorus(L)
            Algebraic Torus of rank 4 defined by the following lattice:
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
            and an action by the Galois group of the form:
            Dihedral group of order 8 as a permutation group
        """
        Scheme.__init__(self)
        self._galois_group = lattice.group()
        self._lattice = lattice
        self._field = self._galois_group.number_field()
        self._base_field = self._field.base_field()
        self._splitting_field = self._field

    def _repr_(self):
        r"""
        The print representation of an algebraic torus.

        EXAMPLES::

            sage: AlgebraicTorus(Lattice_ambient(CyclicPermutationGroup(3), 2))
            Algebraic Torus of rank 2 defined by the following lattice:
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            and an action by the Galois group of the form:
            Cyclic group of order 3 as a permutation group

        ::

            sage: L = Lattice_ambient(CyclicPermutationGroup(3), 2).norm_one_restriction_of_scalars(SymmetricGroup(3
            ....: ),False)
            sage: AlgebraicTorus(L)
            Algebraic Torus of rank 4 defined by the following lattice:
            Free module of degree 4 and rank 3 over Integer Ring
            Echelon basis matrix:
            [ 1  0  0 -1]
            [ 0  1  0 -1]
            [ 0  0  1 -1]
            and an action by the Galois group of the form:
            Symmetric group of order 3! as a permutation group
        """
        split = self._lattice.action_is_trivial()
        if split:
            s = "Split algebraic"
        else:
            s = "Algebraic"
        s += " torus of rank %s" % self.rank()
        if self._base_field is not None:
            s += " over %s" % (self._base_field)
        if not split:
            s += " split by a degree %s extension" % self.splitting_degree()
        return s

    def rank(self):
        r"""
        The rank of the torus.

        EXAMPLES::

            sage: L = Lattice_ambient(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = Lattice_ambient(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = Lattice_ambient(SymmetricGroup(3), [act1, act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T1.rank()
            1
            sage: T2.rank()
            1
            sage: T3.rank()
            3
        """
        return self._lattice._rank

    def splitting_degree(self):
        r"""
        The degree of the minimal Galois extension trivializing this torus.

        EXAMPLES::

            sage: L = Lattice_ambient(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = Lattice_ambient(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = Lattice_ambient(SymmetricGroup(3), [act1,act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T1.splitting_degree()
            1
            sage: T2.splitting_degree()
            1
            sage: T3.splitting_degree()
            6
        """
        return self._lattice._group.order() // self._lattice.action_kernel().order()

    def galois_group(self):
        r"""
        The abstract Galois group of a splitting field of the torus.

        .. NOTE::

            It doesn't have to be a minimal splitting field, therefore we
            allow trivial Galois action.

        EXAMPLES::

            sage: L = Lattice_ambient(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = Lattice_ambient(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = Lattice_ambient(SymmetricGroup(3), [act1,act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T1.galois_group()
            Permutation Group with generators [()]
            sage: T2.galois_group()
            Symmetric group of order 3! as a permutation group
            sage: T3.galois_group()
            Symmetric group of order 3! as a permutation group
        """
        return self._lattice._group

    def is_rational(self, ruple, subgp=None):
        r"""
        Detect if a point is rational over the base field.

        INPUT:

        - ``ruple`` -- a point of the torus given as an r-tuple of points over the
          splitting field, where r is the rank of the torus

        - ``subgp`` -- optional, subgroup of the Galois group corresponding to an
          intermediate extension; if specified, the algorithm will test the rationality
          over the intermediate extension

        EXAMPLES::

            sage: K.<w> = QuadraticField(5); G = K.galois_group()
            sage: Lat = GLattice(G, 1)
            sage: T = AlgebraicTorus(Lat)
            sage: T.is_rational([w, w])
            False
            sage: T.is_rational([1, 1])
            True
            sage: T.is_rational([1, w])
            False
        """
        if subgp is None:
            elt = matrix(self.rank(), ruple)
            res = elt
            galgen = self.galois_group().gens()
            for g in range(len(galgen)):
                res = self.cocharacter_lattice()._action_matrices[g] * matrix(self.rank(), [galgen[g](x) for x in ruple])
                if not res == elt:
                    return False
            return True

        else:
            elt = matrix(self.rank(), ruple)
            res = elt
            hgen = subgp.gens()
            colat = self.cocharacter_lattice().subgroup_lattice(subgp)
            for g in range(len(hgen)):
                res = colat._action_matrices[g] * matrix(self.rank(), [hgen[g](x) for x in ruple])
                if not res == elt:
                    return False
            return True

    def character_lattice(self):
        r"""
        The character lattice of the torus.

        EXAMPLES::

            sage: L = Lattice_ambient(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = Lattice_ambient(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = Lattice_ambient(SymmetricGroup(3), [act1, act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T1.character_lattice()
            Ambient free module of rank 1 over the principal ideal domain Integer Ring
            sage: T2.character_lattice()
            Ambient free module of rank 1 over the principal ideal domain Integer Ring
            sage: T3.character_lattice()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        return self._lattice

    def cocharacter_lattice(self):
        r""""
        The cocharacter lattice of the torus.

        EXAMPLES::

            sage: L = Lattice_ambient(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = Lattice_ambient(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = Lattice_ambient(SymmetricGroup(3), [act1, act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T3.character_lattice()._action_matrices
            [
            [0 1 0]  [0 1 0]
            [0 0 1]  [1 0 0]
            [1 0 0], [0 0 1]
            ]
            sage: T3.cocharacter_lattice()._action_matrices
            [
            [0 1 0]  [0 1 0]
            [0 0 1]  [1 0 0]
            [1 0 0], [0 0 1]
            ]

        The matrices are all orthogonal so we get the same lattice. The action will be
        different in the next example.

        ::

            sage: Lattice_ambient(PermutationGroup([(1,2,3,4,5,6)]), [matrix(2, [0,1,-1,-1])])
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: T = AlgebraicTorus(_)
            sage: T.character_lattice()._action_matrices
            [
            [ 0  1]
            [-1 -1]
            ]
            sage: T.cocharacter_lattice()._action_matrices
            [
            [-1  1]
            [-1  0]
            ]
        """
        return self._lattice.colattice()

    def Tate_Cohomology(self, n):
        r"""
        Gives the isomorphism type of the nth cohomology group using Tate-Nakayama duality. 
        Only works when the base field is local.

        INPUT:

        - ``n`` -- which cohomology group to compute

        NOTE::

            This currently only works for tori over local p-adic fields. For global fields,
            Tate-Nakayama gives the cohomology of the class group, not the torus itself.

        EXAMPLES::

            sage: L = Lattice_ambient(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = Lattice_ambient(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = Lattice_ambient(SymmetricGroup(3), [act1, act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: for i in range(-5, 6):
            ....:     print("H^"+str(i)+": ") , T1.Tate_Cohomology(i)
            H^-5:  []
            H^-4:  []
            H^-3:  []
            H^-2:  []
            H^-1:  []
            H^0:  []
            H^1:  []
            H^2:  []
            H^3:  []
            H^4:  []

        The Galois group is trivial and has obviously trivial cohomology

        ::

            sage: for i in range(-5, 6):
            ....:     print("H^"+str(i)+": ") , T2.Tate_Cohomology(i)
            H^-5:  []
            H^-4:  [2]
            H^-3:  []
            H^-2:  [6]
            H^-1:  []
            H^0:  [2]
            H^1:  []
            H^2:  [6]
            H^3:  []
            H^4:  [2]

        We can recognize from class field theory that H^2 (the Brauer group of
        our extension) is isomorphic to the cyclic group Cn where n is the order of the
        Galois group. Here the group is S3, which has order 6 so we get C6.

        Another way to see it is seeing this H^2 as H^0 of its character lattice. Since the group
        acts trivially, the fixed elements are the whole lattice, and the trace map is multiplication
        by the order of the group, which is 6, so we get C6^(rank of T1)

        Also, ``H^0`` can be seen as the abelianization of the Galois group (indeed, by Tate-Nakayama it is
        ``H^2(\ZZ)``), which here has order 2 (it is the group of signatures)::

            sage: for i in range(-5, 6):
            ....:     print("H^"+str(i)+": ") , T3.Tate_Cohomology(i)
            H^-5:  []
            H^-4:  [2]
            H^-3:  []
            H^-2:  [2]
            H^-1:  []
            H^0:  [2]
            H^1:  []
            H^2:  [2]
            H^3:  []
            H^4:  [2]

        In this example, we can see the 2-periodicity of the cohomology groups, consequence
        of the group being cyclic.
        """
        return self._lattice.Tate_Cohomology(2 - n)

    def product(self, torus):
        r"""
        Return the product of two tori defined and (for now) also splitting over the same field.

        INPUT:

        - ``torus`` -- another algebraic torus

        EXAMPLES::

            sage: T1 = AlgebraicTorus(Lattice_ambient([], 1));
            sage: T2 = T1.norm_one_restriction(G)
            sage: G = SymmetricGroup(3)
            sage: T3 = T1.restriction_of_scalars(G)
            sage: PT = T2.product(T3); PT; PT.character_lattice()._action_matrices
            Algebraic Torus of rank 11 defined by the following lattice:
            Ambient free module of rank 11 over the principal ideal domain Integer Ring
            and an action by the galois group of the form:
            Symmetric group of order 3! as a permutation group
            [
            [ 0 -1  0  0  1| 0  0  0  0  0  0]  [ 0  0  1 -1  0| 0  0  0  0  0  0]
            [ 0 -1  1  0  0| 0  0  0  0  0  0]  [ 0  0  0 -1  1| 0  0  0  0  0  0]
            [ 0 -1  0  0  0| 0  0  0  0  0  0]  [ 1  0  0 -1  0| 0  0  0  0  0  0]
            [ 1 -1  0  0  0| 0  0  0  0  0  0]  [ 0  0  0 -1  0| 0  0  0  0  0  0]
            [ 0 -1  0  1  0| 0  0  0  0  0  0]  [ 0  1  0 -1  0| 0  0  0  0  0  0]
            [--------------+-----------------]  [--------------+-----------------]
            [ 0  0  0  0  0| 0  0  0  0  1  0]  [ 0  0  0  0  0| 0  0  1  0  0  0]
            [ 0  0  0  0  0| 0  0  1  0  0  0]  [ 0  0  0  0  0| 0  0  0  0  1  0]
            [ 0  0  0  0  0| 0  0  0  0  0  1]  [ 0  0  0  0  0| 1  0  0  0  0  0]
            [ 0  0  0  0  0| 1  0  0  0  0  0]  [ 0  0  0  0  0| 0  0  0  0  0  1]
            [ 0  0  0  0  0| 0  0  0  1  0  0]  [ 0  0  0  0  0| 0  1  0  0  0  0]
            [ 0  0  0  0  0| 0  1  0  0  0  0], [ 0  0  0  0  0| 0  0  0  1  0  0]
            ]
        """
        return AlgebraicTorus(self.character_lattice().sum_lattice(torus.character_lattice()))

        #gives the torus representing the Restriction of scalars.
        #Right now, for a torus defined over K, splitting over L,
        #to compute the restriction of scalars to k inside K,
        #the user has to enter the galois group of the extension L/k
        #In the future, when we will have a better notion for Galois group
        #perhaps we can deal with fields directly.

    def restriction_of_scalars(self, group):
        r"""
        The torus obtained through restriction of scalars.

        INPUT:

        - ``group`` -- the bigger group corresponding the the Galois group
        of the splitting field over the subfield one wishes to restrict scalars.

        .. NOTE::

            The user has to input the (larger) Galois group for this extension.
            More concretely, if the torus is defined over K, splits over L, and
            is defined by the action of Gal(L/K) on its character lattice, then
            if one wants the restriction of scalars to a smaller field k,
            one has to enter Gal(L/k) as argument of this method. This is because 
            galois groups of relative extensions are not supported so far. For 
            now the user has to input the Galois group of the relative extension
            as abstract group.

        EXAMPLES::

            sage: L = Lattice_ambient(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = Lattice_ambient(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = Lattice_ambient(SymmetricGroup(3), [act1, act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T1.restriction_of_scalars(PermutationGroup([(1,2), (3,4), (5,6), (7,8)]))
            Algebraic Torus of rank 16 defined by the following lattice:
            Ambient free module of rank 16 over the principal ideal domain Integer Ring
            and an action by the Galois group of the form:
            Permutation Group with generators [(7,8), (5,6), (3,4), (1,2)]

        ::

            sage: T2.restriction_of_scalars(SymmetricGroup(4))
            Algebraic Torus of rank 4 defined by the following lattice:
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
            and an action by the Galois group of the form:
            Symmetric group of order 4! as a permutation group
            sage: _.character_lattice()._action_matrices
            [
            [0|0|0|1]  [1|0|0|0]
            [-+-+-+-]  [-+-+-+-]
            [1|0|0|0]  [0|1|0|0]
            [-+-+-+-]  [-+-+-+-]
            [0|1|0|0]  [0|0|0|1]
            [-+-+-+-]  [-+-+-+-]
            [0|0|1|0], [0|0|1|0]
            ]

        ::

            sage: T3.restriction_of_scalars(SymmetricGroup(4))
            Algebraic Torus of rank 12 defined by the following lattice:
            Ambient free module of rank 12 over the principal ideal domain Integer Ring
            and an action by the Galois group of the form:
            Symmetric group of order 4! as a permutation group
            sage: _.character_lattice()._action_matrices
            [
            [0 0 0|0 0 0|0 0 0|1 0 0]  [0 1 0|0 0 0|0 0 0|0 0 0]
            [0 0 0|0 0 0|0 0 0|0 1 0]  [1 0 0|0 0 0|0 0 0|0 0 0]
            [0 0 0|0 0 0|0 0 0|0 0 1]  [0 0 1|0 0 0|0 0 0|0 0 0]
            [-----+-----+-----+-----]  [-----+-----+-----+-----]
            [0 1 0|0 0 0|0 0 0|0 0 0]  [0 0 0|0 1 0|0 0 0|0 0 0]
            [0 0 1|0 0 0|0 0 0|0 0 0]  [0 0 0|1 0 0|0 0 0|0 0 0]
            [1 0 0|0 0 0|0 0 0|0 0 0]  [0 0 0|0 0 1|0 0 0|0 0 0]
            [-----+-----+-----+-----]  [-----+-----+-----+-----]
            [0 0 0|0 1 0|0 0 0|0 0 0]  [0 0 0|0 0 0|0 0 0|1 0 0]
            [0 0 0|0 0 1|0 0 0|0 0 0]  [0 0 0|0 0 0|0 0 0|0 1 0]
            [0 0 0|1 0 0|0 0 0|0 0 0]  [0 0 0|0 0 0|0 0 0|0 0 1]
            [-----+-----+-----+-----]  [-----+-----+-----+-----]
            [0 0 0|0 0 0|0 1 0|0 0 0]  [0 0 0|0 0 0|1 0 0|0 0 0]
            [0 0 0|0 0 0|0 0 1|0 0 0]  [0 0 0|0 0 0|0 1 0|0 0 0]
            [0 0 0|0 0 0|1 0 0|0 0 0], [0 0 0|0 0 0|0 0 1|0 0 0]
            ]
        """
        return AlgebraicTorus(self._lattice.induced_lattice(group))

    def norm_one_restriction(self,group):
        r"""
        Torus of norm one elements in the restriction of scalars.

        INPUT:

        - ``group`` -- the abstract Galois group of the restriction of
        scalars.

        EXAMPLES::

            sage: L = Lattice_ambient(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = Lattice_ambient(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = Lattice_ambient(SymmetricGroup(3), [act1,act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T1.norm_one_restriction(PermutationGroup([(1,2), (3,4), (5,6), (7,8)]))
            Algebraic Torus of rank 15 defined by the following lattice:
            Ambient free module of rank 15 over the principal ideal domain Integer Ring
            and an action by the Galois group of the form:
            Permutation Group with generators [(7,8), (5,6), (3,4), (1,2)]
            sage: _.character_lattice()._action_matrices[0]
            [ 0  1  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 1  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  1  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  1  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  1  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  1  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  1  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  1  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  1  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  1  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  1  0  0]
            [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
            sage: T2.norm_one_restriction(SymmetricGroup(4))
            Algebraic Torus of rank 3 defined by the following lattice:
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            and an action by the Galois group of the form:
            Symmetric group of order 4! as a permutation group
            sage: _.character_lattice()._action_matrices
            [
            [-1 -1 -1]  [ 1  0  0]
            [ 1  0  0]  [ 0  1  0]
            [ 0  1  0], [-1 -1 -1]
            ]
            sage: N1 = T3.norm_one_restriction(SymmetricGroup(4)); N1
            Algebraic Torus of rank 11 defined by the following lattice:
            Ambient free module of rank 11 over the principal ideal domain Integer Ring
            and an action by the Galois group of the form:
            Symmetric group of order 4! as a permutation group
            sage: N1.character_lattice()._action_matrices
            [
            [ 0  0  0  0  0  0  0  0  0  1  0]  [ 0  1  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  1]  [ 1  0  0  0  0  0  0  0  0  0  0]
            [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]  [ 0  0  1  0  0  0  0  0  0  0  0]
            [ 0  1  0  0  0  0  0  0  0  0  0]  [ 0  0  0  0  1  0  0  0  0  0  0]
            [ 0  0  1  0  0  0  0  0  0  0  0]  [ 0  0  0  1  0  0  0  0  0  0  0]
            [ 1  0  0  0  0  0  0  0  0  0  0]  [ 0  0  0  0  0  1  0  0  0  0  0]
            [ 0  0  0  0  1  0  0  0  0  0  0]  [ 0  0  0  0  0  0  0  0  0  1  0]
            [ 0  0  0  0  0  1  0  0  0  0  0]  [ 0  0  0  0  0  0  0  0  0  0  1]
            [ 0  0  0  1  0  0  0  0  0  0  0]  [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
            [ 0  0  0  0  0  0  0  1  0  0  0]  [ 0  0  0  0  0  0  1  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0], [ 0  0  0  0  0  0  0  1  0  0  0]
            ]

        We now compute the cohomologies of all those tori.
        For the two latter examples, we check that we get different cohomologies with the
        same group, and the norm one restriction of the split torus has nontrivial cohomology::

            sage: ROS = T1.norm_one_restriction(PermutationGroup([(1,2), (3,4), (5,6), (7,8)]))
            sage: for i in range(-4, 6):
            ....:     print("H^"+str(i)+": ") , ROS.Tate_Cohomology(i)
            H^-4:  [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            H^-3:  [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            H^-2:  [2, 2, 2, 2, 2, 2]
            H^-1:  [2, 2, 2, 2]
            H^0:  []
            H^1:  [16]
            H^2:  []
            H^3:  [16]
            H^4:  [2, 2, 2, 2, 2, 2]
            H^5:  [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

        This torus ``ROS`` is the example of Ono where he applies his formula for the Tamagawa
        number of a Torus. See his paper 'On the Tamagawa Number of Algebraic Tori'.

        ::

            sage: ROS2 = T2.norm_one_restriction(SymmetricGroup(4))
            sage: for i in range(-4, 6):
            ....:     print("H^"+str(i)+": ") , T2.Tate_Cohomology(i)
            H^-4:  [2]
            H^-3:  []
            H^-2:  [6]
            H^-1:  []
            H^0:  [2]
            H^1:  []
            H^2:  [6]
            H^3:  []
            H^4:  [2]
            H^5:  []

            sage: ROS3 = T3.norm_one_restriction(SymmetricGroup(4))
            sage: for i in range(-4, 6):
            ....:     print("H^"+str(i)+": ") , ROS3.Tate_Cohomology(i)
            H^-4:  [2, 2]
            H^-3:  [2, 12]
            H^-2:  [2, 2]
            H^-1:  [2]
            H^0:  [2]
            H^1:  [12]
            H^2:  []
            H^3:  [12]
            H^4:  [2]
            H^5:  [12]
        """
        return AlgebraicTorus(self._lattice.norm_one_restriction_of_scalars(group, True))
