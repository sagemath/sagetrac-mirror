# -*- coding: utf-8 -*-
r"""
Cubic Braid Groups

This module is devoted to factor groups of the Artin braid groups, such
that the images `s_i` of the braid generators have order three:

.. MATH::

    s_i^3 = 1

In general these groups have firstly been investigated by Coxeter, H.S.M.
in: "Factor groups of the braid groups, Proceedings of the Fourth Canadian
Mathematical Congress (Vancover 1957), pp. 95-122".

Coxeter showed, that these groups are finite as long as the number of
strands is less than 6 and infinite elsewise. More explicitely the factor
group on three strand braids is isomorphic to `SL(2,3)`, on four strand
braids to `GU(3,2)` and on five strand braids to `Sp(4,3)  \times C_3`.
Today, these finite groups are known as irreducible complex reflection groups
enumerated in the Shephard-Todd classification as `G_{4}`, `G_{25}` and
`G_{32}`.

Coxeter realized these groups as subgroups of unitary groups with respect
to a certain hermitian form over the complex numbers (in fact over `\QQ`
adjoined with a primitive 12-th root of unity).

In "Einige endliche Faktorgruppen der Zopfgruppen" (Math. Z., 163 (1978),
291-302) J. Assion considered two series `S(m)` and `U(m)` of finite
factors of these groups. The additional relations on the braid group
generators `\{ s_1, \cdot , s_{m-1}\}` are

.. MATH::

    \begin{array}{lll}
    \mbox{S:} & s_3 s_1 t_2 s_1 t_2^{-1} t_3 t_2 s_1 t_2^{-1} t_3^{-1} = 1
              & \mbox{ for } m >= 5 \mbox{ in case } S(m)\\
    \mbox{U:} & t_1 t_3 = 1
              & \mbox{ for } m >= 5 \mbox{ in case } U(m)
    \end{array}

where `t_i = (s_i s_{i+1})^3`. He showed that each series of finite cubic
braid group factors must be an epimorphic image of one of his two series,
as long as the groups with less than 5 strands are the full cubic braid
groups, whereas the group on 5 strands is not. He realized the groups `S(m)`
as symplectic groups over `GF(3)` (resp. subgroups therein) and `U(m)` as
general unitary groups over `GF(4)` (resp. subgroups therein).

This class implements all the groups considered by Coxeter and Assion as
finitely presented groups together with the classical realizations given
by the authors. It also contains the conversion maps between the two ways
of realization. In addition the user can construct other realizations and
maps to matrix groups with help of the burau representation. In case gap3
and CHEVIE are installed the reflection groups (via the gap3 interface)
are availlable, too. The methods for all this functionality are
:meth:`as_classical_group`, :meth:`as_matrix_group`, :meth:`as_permutatopn_group`
and :meth:`as_reflection_group`.


REFERENCES:

- [Cox1957]_
- [Ass1978]_

AUTHORS:

- Sebastian Oehms 2019-02-16, initial version.
"""


# ****************************************************************************
#       Copyright (C) 2019 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************


from sage.rings.integer import Integer
from sage.groups.free_group import FreeGroup
from sage.groups.finitely_presented import FinitelyPresentedGroup, FinitelyPresentedGroupElement
from sage.groups.braid import BraidGroup
from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.finite_rings.finite_field_constructor import GF
from enum import Enum


##############################################################################
#
#                  Enum for the type of the group
#
##############################################################################
class CubicBraidGroupType(Enum):
    r"""
    Enum class to select the type of the group:

    - ``Coxeter`` -- 'C' the full cubic braid group.
    - ``AssionS`` -- 'S' finite factor group of type S considered by Assion.
    - ``AssionU`` -- 'U' finite factor group of type U considered by Assion.
    """
    Coxeter = 'C'
    AssionS = 'S'
    AssionU = 'U'


##############################################################################
#
#   Functions to create Instances of the CubicBraidGroup_class
#
##############################################################################
def CubicBraidGroup(n=None, names='c', cbg_type=CubicBraidGroupType.Coxeter):
    r"""
    Constructs cubic braid groups as instance of :class:`CubicBraidGroup_class`.

    INPUT:

    - ``n`` -- integer or ``None`` (default). The number of
      strands. If not specified the ``names`` are counted and the
      group is assumed to have one more strand than generators.

    - ``names`` -- string or list/tuple/iterable of strings (default:
      ``'c'``). The generator names or name prefix.

    - ``cbg_type`` --  (optional keyword, default = CubicBraidGroupType.Coxeter)
      is passed to the corresponding keyword argument of the constructor of
      :class:`CubicBraidGroup_class`.

    EXAMPLES::

        sage: C3 = CubicBraidGroup(3); C3.generators()
        (c0, c1)
        sage: CubicBraidGroup(3, 'g').generators()
        (g0, g1)
        sage: from sage.groups.cubic_braid import CubicBraidGroupType
        sage: U3.<u1,u2>=CubicBraidGroup(3, cbg_type=CubicBraidGroupType.AssionU); U3.generators()
        (u1, u2)
    """
    # this code is adapted from :func:`BraidGroup`
    # Support Freegroup('a,b') syntax
    if n is not None:
        try:
            n = Integer(n)-1
        except TypeError:
            names = n

            n = None
    # derive n from counting names
    if n is None:
        import six
        if isinstance(names, six.string_types):
            n = len(names.split(','))
        else:
            names = list(names)
            n = len(names)

    from sage.structure.category_object import normalize_names
    names = tuple(normalize_names(n, names))
    return CubicBraidGroup_class(names, cbg_type=cbg_type)


# ----------------------------------------------------------------------------------
# Short-Hands for AssionsGroups:
# ----------------------------------------------------------------------------------
def AssionGroupS(n=None, names='s'):
    r"""
    Constructs cubic braid groups as instance of :class:`CubicBraidGroup_class` which
    have been investigated by J.Assion using the notation S(m). This function is a short
    hand cut for :func:`CubicBraidGroup` setting ``cbg_type=CubicBraidGroupType.AssionS``
    and default ``names='s'``.

    For more information type ``CubicBraidGroup_class?``

    INPUT:

    - ``n`` -- integer or None (default). The number of strands. This argument is passed
      to the corresponding argument of :func:`CubicBraidGroup`.

    - ``names`` -- string or list/tuple/iterable of strings (default:'s'). This argument is
      passed to the corresponding argument of :func:`CubicBraidGroup`.

    EXAMPLES::

        sage: S3 = AssionGroupS(3);  S3
        Assion group on 3 strands of type S
        sage: from sage.groups.cubic_braid import CubicBraidGroupType
        sage: S3x = CubicBraidGroup(3, names='s', cbg_type=CubicBraidGroupType.AssionS); S3x
        Assion group on 3 strands of type S
        sage: S3 == S3x
        True
    """
    return CubicBraidGroup(n = n, names = names, cbg_type=CubicBraidGroupType.AssionS)


def AssionGroupU(n=None, names='u'):
    r"""
    Constructs cubic braid groups as instance of :class:`CubicBraidGroup_class` which
    have been investigated by J.Assion using the notation U(m). This function is a short
    hand cut for :func:`CubicBraidGroup` setting ``cbg_type=CubicBraidGroupType.AssionU``
    and default ``names='u'``.

    For more information type ``CubicBraidGroup_class?``

    INPUT:

    - ``n`` -- integer or None (default). The number of strands. This argument is passed
      to the corresponding argument of :func:`CubicBraidGroup`.

    - ``names`` -- string or list/tuple/iterable of strings (default:'u'). This argument is
      passed to the corresponding argument of :func:`CubicBraidGroup`.

    EXAMPLES::

        sage: U3 = AssionGroupU(3);  U3
        Assion group on 3 strands of type U
        sage: from sage.groups.cubic_braid import CubicBraidGroupType
        sage: U3x = CubicBraidGroup(3, names='u', cbg_type=CubicBraidGroupType.AssionU); U3x
        Assion group on 3 strands of type U
        sage: U3 == U3x
        True

    """
    return CubicBraidGroup(n = n, names = names, cbg_type=CubicBraidGroupType.AssionU)



##############################################################################
#
#                  Class CubicBraidElement (for elements)
#
##############################################################################
class CubicBraidElement(FinitelyPresentedGroupElement):
    r"""
    This class models elements of cubic factor groups of the braid group.
    It is the element class of the CubicBraidGroup_class.

    For more information see the documentation of the parent
    :class:`CubicBraidGroup_class`.

    EXAMPLES::

        sage: C4.<c1,c2,c3> = CubicBraidGroup(4); C4
        Cubic Braid group on 4 strands
        sage: ele1 = c1*c2*c3^-1*c2^-1
        sage: ele2 = C4((1, 2, -3, -2))
        sage: ele1 == ele2
        True
    """

    def braid(self):
        r"""
        Return the canonical braid preimage of ``self`` as Object of the
        class :class:`Braid`.

        OUTPUT:

        The preimage of ``self`` as instance of :class:`Braid`.

        EXAMPLES::

            sage: C3.<c1, c2> = CubicBraidGroup(3)
            sage: c1.parent()
            Cubic Braid group on 3 strands
            sage: c1.braid().parent()
            Braid group on 3 strands
        """

        braid_group = self.parent().braid_group()
        return braid_group(self)


    def burau_matrix(self, root_bur = None, domain = None, characteristic = None, var='t', reduced=False):
        r"""
        Return the Burau matrix of the cubic braid coset.

        This method uses the same method belonging to :class:`Braid`, but
        reduces the indeterminant to a primitive sixth (resp. twelfth in case
        reduced='unitary') root of unity.

        INPUT (all arguments are optional keywords):

        - ``root_bur`` -- six (resp. twelfth) root of unity in some field
          (default root of unity over `\QQ`).
        - ``domain``  -- base_ring for the burau matrix (default is Cyclotomic
          Field of order 3 and degree 2, resp. the domain of `root_bur` if given).
        - ``characteristic`` - integer giving the characteristic of the
          domain (default is 0 or the characteristic of `domain` if given).
        - ``var`` -- string used for the indeterminant name in case root_bur
          must be constructed in a splitting field.
        - ``reduced`` -- boolean (default: ``False``) or string; for more
          information see the documentation of :meth:`burau_matrix` of
          :class:`Braid`.

        OUTPUT:

        The Burau matrix of the cubic braid coset with entries in the
        domain given by the options. In case the option `reduced='unitary'`
        is given a triple consisting of the burau matrix, its adjoined and
        the hermitain form is returned.

        EXAMPLES::

            sage: C3.<c1, c2> = CubicBraidGroup(3)
            sage: ele = c1*c2*c1
            sage: BuMa = ele.burau_matrix(); BuMa
            [  -zeta3         1     zeta3]
            [  -zeta3 zeta3 + 1         0]
            [       1         0         0]
            sage: BuMa.base_ring()
            Cyclotomic Field of order 3 and degree 2
            sage: BuMa == ele.burau_matrix(characteristic = 0)
            True
            sage: BuMa = ele.burau_matrix(domain=QQ); BuMa
            [-t + 1      1  t - 1]
            [-t + 1      t      0]
            [     1      0      0]
            sage: BuMa.base_ring()
            Number Field in t with defining polynomial t^2 - t + 1
            sage: BuMa = ele.burau_matrix(domain = QQ[I, sqrt(3)]); BuMa
            [ 1/2*sqrt3*I + 1/2                  1 -1/2*sqrt3*I - 1/2]
            [ 1/2*sqrt3*I + 1/2 -1/2*sqrt3*I + 1/2                  0]
            [                 1                  0                  0]
            sage: BuMa.base_ring()
            Number Field in I with defining polynomial x^2 + 1 over its base field
            sage: BuMa = ele.burau_matrix(characteristic=7); BuMa
            [3 1 4]
            [3 5 0]
            [1 0 0]
            sage: BuMa.base_ring()
            Finite Field of size 7
            sage: BuMa = ele.burau_matrix(characteristic=2); BuMa
            [t + 1     1 t + 1]
            [t + 1     t     0]
            [    1     0     0]
            sage: BuMa.base_ring()
            Finite Field in t of size 2^2
            sage: F4.<r64> = GF(4)
            sage: BuMa = ele.burau_matrix(root_bur=r64); BuMa
            [r64 + 1       1 r64 + 1]
            [r64 + 1     r64       0]
            [      1       0       0]
            sage: BuMa.base_ring()
            Finite Field in r64 of size 2^2
            sage: BuMa = ele.burau_matrix(domain=GF(5)); BuMa
            [2*t + 2       1 3*t + 3]
            [2*t + 2 3*t + 4       0]
            [      1       0       0]
            sage: BuMa.base_ring()
            Finite Field in t of size 5^2
            sage: BuMa, BuMaAd, H = ele.burau_matrix(reduced='unitary'); BuMa
            [       0 zeta12^3]
            [zeta12^3        0]
            sage: BuMa * H * BuMaAd == H
            True
            sage: BuMa.base_ring()
            Cyclotomic Field of order 12 and degree 4
            sage: BuMa, BuMaAd, H  = ele.burau_matrix(domain = QQ[I, sqrt(3)], reduced='unitary'); BuMa
            [0 I]
            [I 0]
            sage: BuMa.base_ring()
            Number Field in I with defining polynomial x^2 + 1 over its base field
        """
        braid = self.braid()

        unitary = False
        min_pol_root_bur_cf  = [1 , -1 , 1] # minimal polynomial of sixth root of unity
        if type(reduced) == str:
            if reduced == 'unitary':
                unitary = True
                min_pol_root_bur_cf  = [1 , 0 , -1 , 0 , 1] # minimal polynomial of twelfth root of unity

        burau_ori = braid.burau_matrix(reduced=reduced)

        if unitary:
            burau_ori, burau_ori_adj, herm_form_ori = burau_ori

        if domain is not None:
            if isinstance(domain, UniversalCyclotomicField):
                if  root_bur is None:
                    if unitary:
                        root_bur = domain.gen(12)
                    else:
                        root_bur = domain.gen(6)

        if root_bur is None:
            def find_root(domain):
                min_pol_root_bur = domain[var](min_pol_root_bur_cf)
                root_list = min_pol_root_bur.roots()
                if len(root_list) == 0:
                    domain = min_pol_root_bur.splitting_field(min_pol_root_bur.variable_name())
                    min_pol_root_bur = domain[var](min_pol_root_bur_cf)
                else:
                    domain = domain
                root_list = min_pol_root_bur.roots()
                for root in root_list:
                   if root[0] == 0:
                        continue
                   root_bur = root[0]
                   if root[1]  == 1:
                        break
                return root_bur

            if domain is None:
                if (characteristic is None):
                    characteristic = 0
                try:
                    characteristic = Integer(characteristic)
                except ValueError:
                    raise ValueError('characteristic must be in integer')

                if  characteristic <> 0  and not characteristic.is_prime():
                    raise ValueError('characteristic must be a prime')
                if characteristic == 0:
                    if unitary:
                        domain = CyclotomicField(12)
                    else:
                        domain = CyclotomicField(3)
                else:
                    domain = GF(characteristic)
                root_bur = find_root(domain)
                domain = root_bur.parent()

            else: # domain <> None
                if characteristic is None:
                    characteristic = domain.characteristic()
                elif characteristic <> domain.characteristic():
                    raise ValueError('characteristic of domain does not match given characteristic')
                root_bur = find_root(domain)

        else:  # root_bur <> None
            if domain is None:
                domain = root_bur.parent()
            if characteristic is None:
                characteristic = domain.characteristic()
            elif characteristic <> domain.characteristic():
                raise ValueError('characteristic of domain does not match given characteristic')

            if 1 not in domain:
                raise ValueError('root_bur must belong to a domain containing 1')
            if unitary:
                if domain(root_bur**4  -root_bur**2  +1) != 0:
                    raise ValueError('root_bur must vanish on x**4-x**2+1')
            else:
                if domain(root_bur**2  -root_bur +1) != 0:
                    raise ValueError('root_bur must vanish on x**2-x+1')

        def conv2domain (laur_pol):
            l1, l2 = laur_pol.polynomial_construction()
            R = domain['tt']; (tt,) = R._first_ngens(1)
            p1 = R(l1)
            p2 = root_bur**(l2)
            res = p1.substitute(tt=root_bur) * p2
            return res

        from sage.matrix.constructor import matrix

        d1, d2 = burau_ori.dimensions()
        burau_mat = matrix(d1, d2, lambda i,j: conv2domain(burau_ori[i,j]))

        if unitary:
            burau_mat_adj = matrix(d1, d2, lambda i,j: conv2domain(burau_ori_adj[i,j]))
            herm_form = matrix(d1, d2, lambda i,j: conv2domain(herm_form_ori[i,j]))
            return burau_mat, burau_mat_adj, herm_form

        return burau_mat



##############################################################################
#
#                  Class CubicBraidGroup_class
#
##############################################################################
class CubicBraidGroup_class(FinitelyPresentedGroup):
    r"""
    This class implements factor groups of the Artin braid group mapping
    their generators to elements of order 3 (see the module header for more
    informations on these groups).

    These groups are implemented as a particular case of finitely presented
    groups similar to the :class:`BraidGroup_class`.

    A cubic braid group can be created by giving the number of strands, and
    the name of the generators in a similar way as it works for the
    :class:`BraidGroup_class`.

    INPUT (to the constructor):

    - ``names`` -- (see the documentation of :class:`BraidGroup_class`).

    - ``cbg_type`` -- (keyword, explanation see below).


    Setting the keyword ``cbg_type`` to one on the values ``CubicBraidGroupType.AssionS``
    or ``CubicBraidGroupType.AssionU`` the additional relations due to Assion are added:

    .. MATH::

        \begin{array}{lll}
        \mbox{S:} & s_3 s_1 t_2 s_1 t_2^{-1} t_3 t_2 s_1 t_2^{-1} t_3^{-1} = 1
                  & \mbox{ for } m >= 5 \mbox{ in case } S(m)\\
        \mbox{U:} & t_1 t_3 = 1
                  & \mbox{ for } m >= 5 \mbox{ in case } U(m)
        \end{array}

    where `t_i = (s_i s_{i+1})^3`. If ``cbg_type == CubicBraidGroupType.Coxeter`` (default)
    only the cubic relation on the generators is active (Coxeters case of investigation).
    Note that for `n = 2, 3, 4` the groups do not differ between the three possible
    values of cbg_type (as finitely presented groups). But anyway, the instances for
    ``CubicBraidGroupType.Coxeter, CubicBraidGroupType.AssionS`` and ``CubicBraidGroupType.AssionU``
    are different, since they have different classical realizations implemented.

    The creation of instances of this class can also be done more easily by help
    of  :func:`CubicBraidGroup`, :func:`AssionGroupS` and :func:`AssionGroupU`
    (similar to :func:`BraidGroup` with respect to :class:`BraidGroup_class`).


    EXAMPLES::

        sage: from sage.groups.cubic_braid import CubicBraidGroupType
        sage: U3 = CubicBraidGroup(3, cbg_type=CubicBraidGroupType.AssionU); U3
        Assion group on 3 strands of type U
        sage: U3.gens()
        (c0, c1)

    alternative possibilities defining U3::

        sage: U3 = AssionGroupU(3); U3
        Assion group on 3 strands of type U
        sage: U3.gens()
        (u0, u1)
        sage: U3.<u1,u2> = AssionGroupU(3); U3
        Assion group on 3 strands of type U
        sage: U3.gens()
        (u1, u2)

    alternates naming the generators::

        sage: U3 = AssionGroupU(3, 'a, b'); U3
        Assion group on 3 strands of type U
        sage: U3.gens()
        (a, b)
        sage: C3 = CubicBraidGroup(3, 't'); C3
        Cubic Braid group on 3 strands
        sage: C3.gens()
        (t0, t1)
        sage: U3.is_isomorphic(C3)
        #I  Forcing finiteness test
        True
        sage: U3.as_classical_group()
        Subgroup generated by [(1,7,6)(3,19,14)(4,15,10)(5,11,18)(12,16,20), (1,12,13)(2,15,19)(4,9,14)(5,18,8)(6,21,16)] of (The projective general unitary group of degree 3 over Finite Field of size 2)
        sage: C3.as_classical_group()
        Subgroup with 2 generators (
        [  E(3)^2        0]  [       1 -E(12)^7]
        [-E(12)^7        1], [       0   E(3)^2]
        ) of General Unitary Group of degree 2 over Universal Cyclotomic Field with respect to positive definite hermitian form
        [-E(12)^7 + E(12)^11                  -1]
        [                 -1 -E(12)^7 + E(12)^11]

    TESTS:

        sage: C3 = CubicBraidGroup(3)
        sage: TestSuite(C3).run()
        sage: C4 = CubicBraidGroup(4)
        sage: TestSuite(C4).run()        # long time
        sage: C6 = CubicBraidGroup(6)
        sage: TestSuite(C6).run()        # long time
        sage: S3 = AssionGroupS(3)
        sage: TestSuite(S3).run()
        sage: S5 = AssionGroupS(5)
        sage: TestSuite(S5).run()        # long time
        sage: U3 = AssionGroupU(3)
        sage: TestSuite(U3).run()
        sage: U4 = AssionGroupU(4)
        sage: TestSuite(U4).run()        # long time
        sage: U5 = AssionGroupU(5)
        sage: TestSuite(U5).run()        # long time


    REFERENCES:

    - [Cox1957]_
    - [Ass1978]_
    """

    Element = CubicBraidElement

    ###########################################################################################
    # private methods
    ###########################################################################################
    def __init__(self, names, cbg_type=CubicBraidGroupType.Coxeter):
        """
        Python constructor.
        """
        n  = Integer(len(names))
        if n < 1:
            raise ValueError("the number of strands must be an integer larger than one")

        free_group        = FreeGroup(names)
        self._cbg_type    = cbg_type
        self._nstrands    = n+1
        self._ident       = self._cbg_type.value + self._nstrands.str()
        self._braid_group = BraidGroup(names)

        # internal naming of elements for convinience
        b  = [free_group([i])  for i in range(1 , n+1)]
        t  = [free_group([i,  i+1])**3  for i in range(1 , n)]
        ti = [free_group([-i, -i-1])**3  for i in range(1 , n)]

        # first the braid relation
        rels = list(self._braid_group.relations())

        # than the cubic relation
        for i in range(0, n):
            rels.append(b[i]**3)

        # than Assions relation Satz 2.2 for cbg_type=CubicBraidGroupType.AssionS
        # and Satz 2.4 for cbg_type=CubicBraidGroupType.AssionU
        if n > 3:
            for i in range(0, n-3):
                if cbg_type == CubicBraidGroupType.AssionU:
                    rels.append((t[i]*t[i+2])**3)
                elif cbg_type == CubicBraidGroupType.AssionS:
                    rels.append(b[i+2]*b[i]*t[i+1]*b[i]*ti[i+1]*t[i+2]*t[i+1]*b[i]*ti[i+1]*ti[i+2])

        FinitelyPresentedGroup.__init__(self, free_group, tuple(rels))
        self._free_group = free_group

        # ------------------------------------------------------------------------------------------------
        # the following global pointers to classical group realizations will be set in the private method
        # _create_classical_realization
        # ------------------------------------------------------------------------------------------------
        self._classical_group             = None   # This is the classical Group returned by as_classical_group
        self._classical_base_group        = None   # this only differs for special cases for AssionGroup from the former
        self._classical_invariant_form    = None   # invariant form of the classical base group
        self._classical_embedding         = None   # if self._classical_group different from self._classical_base_group
        self._centralizing_matrix         = None   # for AssionGroup: element in classical base group commuting with self
        self._centralizing_element        = None   # image under nat. map of the former one in the proj. classical group
        return


    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:

        String describing ``self``.
        """
        if self._cbg_type == CubicBraidGroupType.Coxeter:
            return "Cubic Braid group on %s strands"%(self.strands())
        else:
            return "Assion group on %s strands of type %s"%(self.strands() ,self._cbg_type.value)


    # -------------------------------------------------------------------------------
    # Methods for test_suite
    # -------------------------------------------------------------------------------
    def _internal_test_attached_group(self, attached_group, tester):
        r"""
        It tests conversion maps from ``self`` to the given attached Group
        which must have been defined using the :meth:`as_classical_group`,
        :meth:`as_matrix_group`, meth:`as_permutation_group` or
        :meth:`as_reflection_group`.

        INPUT:

         - ``attached_group`` -- attached group to be tested as specified above.
        """
        elem = self.an_element()
        att_grp_elem = attached_group(elem)
        if self.is_finite() and self.strands() <= 7: # not realistic for larger number of strands
            att_grp_elem_back= self(att_grp_elem)
            tester.assert_(att_grp_elem_back == elem)
        return


    def _test_classical_group(self, **options):
        r"""
        Method called by TestSuite.

        The following is checked:
           - construction of classical group was faithfull.
           - coercion maps to and from classical group exist and are inverse to each other.
        """
        tester = self._tester(**options)
        classic_grp = self.as_classical_group()
        if self.is_finite():
            self._internal_test_attached_group(classic_grp, tester)
        return


    def _test_permutation_group(self, **options):
        r"""
        Method called by TestSuite.

        The following is checked:
           - construction of permutation group was faithfull.
           - coercion maps to and from permutation group exist and are inverse to each other.
        """
        if self.is_finite():
            tester = self._tester(**options)
            permgrp = self.as_permutation_group()
            self._internal_test_attached_group(permgrp, tester)
        return


    def _test_matrix_group(self, **options):
        r"""
        Method called by TestSuite.

        The following is checked:
           - construction of matrix group was faithfull in several cases.
           - coercion maps to and from matrix group exist.
        """
        tester = self._tester(**options)
        F3 = GF(3); r63 = F3(2); F4 = GF(4); r64 = F4.gen()
        MatDEF = self.as_matrix_group()
        self._internal_test_attached_group(MatDEF, tester)

        if self._cbg_type != CubicBraidGroupType.AssionU or self.strands() < 5 : # not well defined elsewise
            matrix_grpF3 = self.as_matrix_group(root_bur=r63)
            self._internal_test_attached_group(matrix_grpF3, tester)

        if self._cbg_type != CubicBraidGroupType.AssionS or self.strands() < 5 : # not well defined elsewise
            matrix_grpF4 = self.as_matrix_group(root_bur=r64)
            self._internal_test_attached_group(matrix_grpF4, tester)

        if self.strands() < 5  or self._cbg_type == CubicBraidGroupType.Coxeter:
            matrix_grpF5 = self.as_matrix_group(characteristic=5)
            self._internal_test_attached_group(matrix_grpF5, tester)

            matrix_grpF7 = self.as_matrix_group(domain=GF(7))
            self._internal_test_attached_group(matrix_grpF7, tester)
        return


    def _test_reflection_group(self, **options):
        r"""
        Method called by TestSuite.

        The following is checked:
           - construction of reflection group was faithfull.
           - coercion maps to and from reflection group exist and are inverse to each other.
        """
        if self._cbg_type == CubicBraidGroupType.Coxeter and self.is_finite() and  self.strands() > 2:
            from sage.combinat.root_system.reflection_group_real import is_chevie_available
            if is_chevie_available():
                tester = self._tester(**options)
                reflgrp = self.as_reflection_group()
                self._internal_test_attached_group(reflgrp, tester)
        return


    # -------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------
    # local utility-methods
    # -------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------
    def _create_classical_realization(self, just_embedded=False):
        r"""
        Internal method to create the classical groups attached to ``self``.

        This methods sets the following attributes of ``self``:

         - self._classical_group            This is the classical group returned by as_classical_group method.
         - self._classical_base_group       this only differs in special cases for AssionGroup from the former.
         - self._classical_invariant_form   invariant form of the classical base group.
         - self._centralizing_matrix        for AssionGroup: element in classical base group commuting with self.
         - self._centralizing_element       image under natural map of the former one in the projective classical group.
         - self._classical_embedding        as subgroup of classical base group (if different from classical group).
        """

        # -------------------------------------------------------------------------------
        # Set up all incrediences of the classical Assion group (generic part)
        # -------------------------------------------------------------------------------
        def set_classical_realization(self, base_group, proj_group, centralizing_matrix, transvec_matrices):
            r"""
            Internal method to create classical group for Assion groups.

            This is a local function of :meth:`_create_classical_realization`.

            It handels the common part of symplectic and unitary version and creates conversion maps.

            INPUT:

            - ``base_group`` -- The symplectic or unitary groups Sp(m,3) resp. GU(m,2).
            - ``proj_group`` -- The corresponding projeective group of base_group.
            - ``centralizing_matrix`` -- The centralizing matrix according to Assion.
            - ``transvec_matrices`` -- List of transvection matrices according to Assion.

            OUTPUT:

            No output, but the function sets the attributes of ``self`` decribed above.
            """
            centralizing_element = None

            # ------------------------------------------------------------------------------
            # Setting the List of Braid Images
            # ------------------------------------------------------------------------------
            im_gens  = [base_group(m) for m in transvec_matrices]

            # ------------------------------------------------------------------------------
            # By the work of Assion no check on the group homomorphism is needed, at all.
            # But to take care of software bugs they are perfomed in cases where they are
            # not really expansive.
            # ------------------------------------------------------------------------------
            check = False
            if self.strands() < 7:
                check = True

            # ------------------------------------------------------------------------------
            # Do the projective group realization if needed
            # ------------------------------------------------------------------------------
            embedding      = self._classical_embedding
            classical_group = None
            if proj_group is None:
                classical_group = base_group
                hom_to_classic = self.hom(im_gens, check=check)
                classical_group.register_conversion(hom_to_classic)
            else:
                if embedding is None:
                    im_gens.pop()
                    embedding = base_group.subgroup(im_gens, check=check)
                    embedding.register_conversion(self.hom(embedding.gens(), check=check))
                    hom_to_base = self.hom(im_gens, check=check, codomain=base_group)
                    base_group.register_conversion(hom_to_base)
                if not just_embedded:
                    transvec_matrices.pop()
                    nat_hom = base_group.hom(proj_group.gens(), check=check)
                    centralizing_element = nat_hom(centralizing_matrix)
                    classical_group_gens = [nat_hom(m) for m in transvec_matrices]
                    classical_group     = proj_group.subgroup(classical_group_gens, canonicalize=False)
                    hom_to_classic = self.hom(classical_group.gens(), check=check)
                    classical_group.register_conversion(hom_to_classic)

            # ------------------------------------------------------------------------------
            # register constructed items
            # ------------------------------------------------------------------------------
            self._classical_group             = classical_group
            self._classical_base_group        = base_group
            self._classical_invariant_form    = base_group.invariant_form()
            self._centralizing_matrix         = centralizing_matrix
            self._centralizing_element        = centralizing_element
            self._classical_embedding         = embedding
            return

        # -------------------------------------------------------------------------------
        # local methods to set up the classical group (spezific part)
        # -------------------------------------------------------------------------------
        # Case for symlectic groups
        # -------------------------------------------------------------------------------
        def create_sympl_realization(self, m):
            r"""
            Internal method to create classical group for symplectic
            Assion groups (`cbg_type CubicBraidGroupType.AssionS`).

            INPUT:

            - ``m`` --  Integer, the dimension of the classical groups vectorspace of operation.

            The function calculates the centralizing matrix and the transvections as given by Assion
            and then uses set_classical_realization to complete the construction.
            """

            # -----------------------------------------------------------
            # getting the invariant bilinearform of the group
            # and setting constants.
            # -----------------------------------------------------------
            n = self.strands()

            from sage.groups.matrix_gps.symplectic import Sp
            base_group = Sp(m, 3)
            proj_group = None
            if m == n:
                from sage.groups.perm_gps.permgroup_named import PSp
                proj_group = PSp(m, 3)

            bform = base_group.invariant_form()
            bas = bform.column_space().basis()
            F = bform.base_ring()

            if m % 2  == 0:
                mhalf = m/2
            else:
                mhalf = (m-1)/2

            # -----------------------------------------------------------
            # computing a hyperbolic decomposition basis with respect
            # to the invariant bilinear form.
            # -----------------------------------------------------------
            xbas =[bas[mhalf -i -1] for i in range(mhalf)]
            ybas =[bas[mhalf +i]    for i in range(mhalf)]

            # -----------------------------------------------------------
            # computing the List of transvection vectors according to
            # the Assion paper, page 292.
            # -----------------------------------------------------------
            transvections =[xbas[0]]                                   # t_1      = x_1
            for i in range(mhalf-1):
                transvections.append(ybas[i])                          # t_{2i}   = y_i
                transvections.append(xbas[i] + xbas[i+1])              # t_{2i+1} = x_j + x_(j+1)
            transvections.append(ybas[mhalf-1])                        # t_n      = y_m

            # -----------------------------------------------------------
            # Conversion-Map from transvection vector to transvection
            # matrix.
            # -----------------------------------------------------------
            from sage.matrix.constructor import matrix
            def transvec2mat(v, bas=bas, bform=bform, fact=1):
                t = [x + fact*(x * bform * v) * v for x in bas]
                return matrix(bform.base_ring(),  t)

            # ------------------------------------------------------------------------------
            # setting the centralizing matrix for the case of projective group realizaton
            # ------------------------------------------------------------------------------
            centralizing_vector = xbas[mhalf-1]
            centralizing_matrix = base_group(transvec2mat(centralizing_vector, fact=1))
            transvec_matrices   = [transvec2mat(v) for v in transvections]

            set_classical_realization(self, base_group, proj_group, centralizing_matrix, transvec_matrices)
            return

        # -------------------------------------------------------------------------------
        # Case for unitary groups
        # -------------------------------------------------------------------------------
        def create_unitary_realization(self, m):
            """
            Internal method to create classical group for
            unitary Assion groups (`cbg_type CubicBraidGroupType.AssionU`).

            INPUT:

            - ``m`` --  Integer, the dimension of the classical groups vectorspace of operation.

            The function calculates the centralizing_matrix and the transvections as given by Assion
            and then uses set_classical_realization to complete the construction.
            """

            # ---------------------------------------------------------------------
            # getting the invariant bilinearform of the group
            # and setting constants
            # ---------------------------------------------------------------------
            n = self.strands()

            from sage.groups.matrix_gps.unitary import GU
            base_group = GU(m, 2)
            proj_group = None
            if m == n:
                from sage.groups.perm_gps.permgroup_named import PGU
                proj_group = PGU(m, 2)

            bform = base_group.invariant_form()
            bas = bform.column_space().basis()
            F = bform.base_ring()
            a = F.gen()

            if   m % 3  == 1:
                mthird = (m-1)/3
            elif m % 3  == 2:
                mthird = (m-2)/3
            else:
                mthird = m/3

            # -----------------------------------------------------------
            # computing a orthonormal basis with respect
            # to the invariant bilinear form.
            # -----------------------------------------------------------
            xbas =[]
            for i in range(m):
                if 2*i == m-1:
                    xbas.append(bas[i])
                else:
                    xbas.append(a*bas[i] + a.frobenius()*bas[m-1 -i])

            # -----------------------------------------------------------
            # computing the List of transvection vectors according to
            # Assion paper, page 293.
            # -----------------------------------------------------------
            transvections =[xbas[0]]                                          # t_1 = x_1
            if m > 1:
                transvections.append(xbas[0]+xbas[1]+xbas[2])                 # t_2 = x_1 + x_2 + x_3
            for j in range(mthird):
                pos = 3*(j+1)-1
                transvections.append(xbas[pos-1])                             # t_{3i}   = x_{3i-1}
                if  pos +1  < m:
                    transvections.append(xbas[pos-1]+xbas[pos]+xbas[pos+1])   # t_{3i+1} = x_{3i-1} + x_{3i} + x_{3i+1}
                if  pos +3  < m:
                    transvections.append(xbas[pos+1]+xbas[pos+2]+xbas[pos+3]) # t_{3i+2} = x_{3i+1} + x_{3i+2} + x_{3i+3}

            # -----------------------------------------------------------
            # Conversion-Map from transvection vector to transvection
            # matrix.
            # -----------------------------------------------------------
            from sage.matrix.constructor import matrix
            def transvec2mat(v, bas=bas, bform=bform, fact=a):
                # note x does not change under conjugation, since it belongs to standard basis
                t = [x + fact *(x * bform * v.conjugate()) * v for x in bas]
                return matrix(F, t)


            # ------------------------------------------------------------------------------
            # setting the centralizing matrix for the case of projective group realizaton.
            # ------------------------------------------------------------------------------
            centralizing_vector = xbas[m-2]+xbas[m-1]
            centralizing_matrix = base_group(transvec2mat(centralizing_vector, fact=1))
            transvec_matrices   = [transvec2mat(v) for v in transvections]

            set_classical_realization(self, base_group, proj_group, centralizing_matrix, transvec_matrices)
            return

        #----------------------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------------------
        # local functions declaration section finishes here
        #----------------------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------
        # initialization of constants
        # -------------------------------------------------------------------------------
        n = self.strands()

        # -------------------------------------------------------------------------------
        # Setting the Classical group
        # -------------------------------------------------------------------------------
        if   self._cbg_type == CubicBraidGroupType.AssionS:
            dim_sympl_group   = n-1              # S(n-1) = Sp(n-1, 3)
            if n % 2  == 0:
                dim_sympl_group   = n            # S(n-1) = subgroup of PSp(n, 3)
            create_sympl_realization(self,  dim_sympl_group)
        elif self._cbg_type == CubicBraidGroupType.AssionU:
            dim_unitary_group = n-1              # U(n-1) = GU(n-1, 2)
            if n % 3  == 0:
                dim_unitary_group = n            # U(n-1) = subgroup PGU(n, 3)
            create_unitary_realization(self, dim_unitary_group)
        else:
            # -----------------------------------------------------------------------------------------------
            # connection between Coxeter realization and unitary Burau representation according to Squier:
            # -----------------------------------------------------------------------------------------------
            # Notation of Coxeter: p = 3, \theta =\pi/3 = primitive 6th root of unity
            #     i\theta = \pi/3+\pi/2 = 5\pi/6 = 5th power of primitive 12th root of unity
            #     i\theta = z12^5 = - ~z12 where z12 = UCF.gen(12)
            # Let f be the unitary Form of Coxeter and J(s) the one of Squier. Then we have
            #     J(z12) = 2 f
            # Let `buc` be the unitary Burau Matrix of Coxeter for the first braid generator and `bus`
            # the corresponding one according to Squier. Then we have:
            #     buc_[i,i-1] = buc_[i,i+1]= - i\theta, buc_[i,i] = 1 + 2*cos(\pi/6)*i\theta
            #     bus_[i,i-1] = bus_[i,i+1]=   s,       bus_[i,i] = -s^2
            # now 1 + 2*cos(\pi/6)*i\theta = 1 + sqrt(3)*(-sqrt(3)/2 + I/2) = 1- 3/2 + sqrt(3)I/2 = z12^4 = - ~z12^2
            # finally: Coxeters Realization is the unitary Burau representation of Squier for s = ~z12
            # -----------------------------------------------------------------------------------------------
            UCF = UniversalCyclotomicField(); z12 = UCF.gen(12)
            classical_group = self.as_matrix_group(root_bur=~z12, domain=UCF, reduced='unitary')
            self._classical_group            = classical_group
            self._classical_base_group       = classical_group
            self._classical_embedding        = classical_group
            if self._classical_invariant_form is None:
                self._classical_invariant_form  = classical_group.ambient().invariant_form()
        return

    def _element_constructor_(self, x, **kwds):
        r"""
        Extensions to the _element constructor of :class:`FinitelyPresentedGroup`:
        new functionalities are:

            -- constructing element from an element of the attached classical group
               (embedded and not embedded)
            -- constructing element from an element of the attached permutation group
            -- constructing element from an element of the attached refelction group

        INPUT::

        - ``x`` -- can be one of the following:
                -- an instance of the element class of ``self`` (but possible to a different parent).
                -- an instance of the element class of the braid group.
                -- a tuple representing a braid in Tietze form.
                -- an instance of an element class of a parent P such that there is a map from ``self`` to P
                   having :meth:`lift`, for example an element of an alternative realization of ``self``, such
                   as the classical realization.
                -- any other object which works for the element constructor of :class:`FinitelyPresentedGroup`.

        OUTPUT:

        instance of the element class of ``self``

        EXAMPLES::

            sage: S3 = AssionGroupS(3)
            sage: S3Cl = S3.as_classical_group()
            sage: g = mul(S3Cl.gens())
            sage: S3(g)  # indirect doctest
            s0*s1*s0^-1
        """
        if hasattr(x, 'parent'):
            parent = x.parent()
            map_to   = parent.convert_map_from(self)
            if map_to is not None:
                if hasattr(map_to, 'lift'):
                    return map_to.lift(x)
        return super(CubicBraidGroup_class, self)._element_constructor_(x)


    #######################################################################################################################
    # ----------------------------------------------------------------------------------
    # public methods
    # ----------------------------------------------------------------------------------
    #######################################################################################################################

    # ---------------------------------------------------------------------------------------------------------------------
    # strands
    # ---------------------------------------------------------------------------------------------------------------------
    def strands(self):
        r"""
        Return the number of strands of the braid group whose image is ``self``.

        OUTPUT: Integer.

        EXAMPLES::

            sage: C4 = CubicBraidGroup(4)
            sage: C4.strands()
            4
        """
        return self._nstrands


    # ----------------------------------------------------------------------------------
    # braid_group
    # ----------------------------------------------------------------------------------
    def braid_group(self):
        r"""
        Return an Instance of :class:`BraidGroup` with identical generators, such that
        there exists an epimorhism to ``self``.

        OUTPUT:

        Instance of :class:`BraidGroup` having conversion maps to and from ``self``
        (which is just a section in the latter case).

        EXAMPLES::

            sage: U5 = AssionGroupU(5); U5
            Assion group on 5 strands of type U
            sage: B5 = U5.braid_group(); B5
            Braid group on 5 strands
            sage: b = B5([4,3,2,-4,1])
            sage: u = U5([4,3,2,-4,1])
            sage: u == b
            False
            sage: b.burau_matrix()
            [ 1 - t      t      0      0      0]
            [ 1 - t      0      t      0      0]
            [ 1 - t      0      0      0      t]
            [ 1 - t      0      0      1 -1 + t]
            [     1      0      0      0      0]
            sage: u.burau_matrix()
            [  -zeta3 zeta3 + 1         0         0         0]
            [  -zeta3         0 zeta3 + 1         0         0]
            [  -zeta3         0         0         0 zeta3 + 1]
            [  -zeta3         0         0         1     zeta3]
            [       1         0         0         0         0]
            sage: bU = U5(b)
            sage: uB = B5(u)
            sage: bU == u
            True
            sage: uB == b
            True
        """
        return self._braid_group


    # ----------------------------------------------------------------------------------
    # as_matrix_group
    # ----------------------------------------------------------------------------------
    def as_matrix_group(self, root_bur=None, domain=None, characteristic=None, var='t', reduced=False):
        r"""
       Creates an epimorphic image of ``self`` as a matrix group by use of the burau representation.

        INPUT (all arguments are optional by keyword):

        - ``root_bur`` -- six (resp. twelfth) root of unity in some field
          (default root of unity over `\QQ`).
        - ``domain``  -- base_ring for the burau matrix (default is Cyclotomic
          Field of order 3 and degree 2, resp. the domain of `root_bur` if given).
        - ``characteristic`` - integer giving the characteristic of the
          domain (default is 0 or the characteristic of `domain` if given)
          If none of the keywords `root_bur`, `domain` and `characteristic` is
          given the default characteristic is 3 (resp. 2) if ``self`` is of ``cbg_type
          CubicBraidGroupType.AssionS`` (resp. ``CubicBraidGroupType.AssionU``).
        - ``var`` -- string used for the indeterminant name in case `root_bur`
          must be constructed in a splitting field.
        - ``reduced`` -- boolean (default: ``False``); for more information
          see the documentation of :meth:`burau_matrix` of :class:`Braid`.

        OUTPUT:

        An instance of the class :class:`FinitelyGeneratedMatrixGroup_gap` according to the
        input arguments together with a group homomorphism registerd as a conversion
        from ``self`` to it.

        EXAMPLES::

            sage: C5 = CubicBraidGroup(5)
            sage: C5Mch5 = C5.as_matrix_group(characteristic=5); C5Mch5
            Matrix group over Finite Field in t of size 5^2 with 4 generators (
            [2*t + 2 3*t + 4       0       0       0]
            [     1       0       0       0       0]
            [     0       0       1       0       0]
            [     0       0       0       1       0]
            [     0       0       0       0       1],
            <BLANKLINE>
            [     1       0       0       0       0]
            [     0 2*t + 2 3*t + 4       0       0]
            [     0       1       0       0       0]
            [     0       0       0       1       0]
            [     0       0       0       0       1],
            <BLANKLINE>
            [     1       0       0       0       0]
            [     0       1       0       0       0]
            [     0       0 2*t + 2 3*t + 4       0]
            [     0       0       1       0       0]
            [     0       0       0       0       1],
            <BLANKLINE>
            [     1       0       0       0       0]
            [     0       1       0       0       0]
            [     0       0       1       0       0]
            [     0       0       0 2*t + 2 3*t + 4]
            [     0       0       0       1       0]
            )
            sage: c = C5([3,4,-2,-3,1]); c
            c2*c3*c1^-1*c2^-1*c0
            sage: m = C5Mch5(c); m
            [2*t + 2 3*t + 4       0       0       0]
            [     0       0       0       1       0]
            [2*t + 1       0 2*t + 2     3*t 3*t + 3]
            [2*t + 2       0       0 3*t + 4       0]
            [     0       0 2*t + 2 3*t + 4       0]
            sage: m_back = C5(m)
            sage: m_back == c
            True
            sage: U5 = AssionGroupU(5); U5
            Assion group on 5 strands of type U
            sage: U5Mch3 = U5.as_matrix_group(characteristic=3)
            Traceback (most recent call last):
            ...
            ValueError: burau representation does not factor through the relations
        """
        # -------------------------------------------------------------------------------
        # define matrix group by generators using the Burau representation
        # -------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------
        # adaption of default options for Assion groups on more then 4 stands
        # (in order to achieve the maps beeing well defined)
        # -------------------------------------------------------------------------------
        if root_bur is None and domain is None and characteristic is None:
            if   self._cbg_type == CubicBraidGroupType.AssionS:
                characteristic = 3
            elif self._cbg_type == CubicBraidGroupType.AssionU:
                characteristic = 2
            else:
                characteristic = 0

        unitary = False
        if type(reduced) == str:
            if reduced == 'unitary':
                unitary = True

        gen_list =[]
        for braid_gen in self.gens():
            bur_mat = braid_gen.burau_matrix(root_bur=root_bur, domain=domain, characteristic=characteristic,
                     var=var, reduced=reduced)
            if unitary:
                 bur_mat, bur_mat_ad, herm_form = bur_mat

            if domain is None:
                 domain = bur_mat.base_ring()

            gen_list.append(bur_mat)

        if unitary and herm_form.is_singular():
            unitary = False  # since a degenertated hermitian form doesn't define a unitary group
            if self._classical_invariant_form is None:
                self._classical_invariant_form  = herm_form

        if unitary:
            from sage.rings.finite_rings.finite_field_base import is_FiniteField
            from sage.groups.matrix_gps.unitary import GU
            d, d = herm_form.dimensions()
            if is_FiniteField(domain):
                base_group = GU(d, domain, var=domain.gen(), invariant_form=herm_form)
            else:
                base_group = GU(d, domain, invariant_form=herm_form)

            matrix_group = base_group.subgroup(gen_list)
        else:
            from sage.groups.matrix_gps.finitely_generated import MatrixGroup
            matrix_group = MatrixGroup(gen_list)

        # -------------------------------------------------------------------------------
        # check if there is a well defined group homomorphism to matrix_group
        # Register map from ``self`` to matrix_group.
        # Since GAP' check is very expansive (on time and memory), the check is performed
        # here.
        # -------------------------------------------------------------------------------
        hom_to_mat = self.hom(matrix_group.gens(), check=False)
        if not all( hom_to_mat(rel).is_one()  for rel in self.relations()):
            raise ValueError("burau representation does not factor through the relations")
        matrix_group.register_conversion(hom_to_mat)
        return matrix_group


    # ----------------------------------------------------------------------------------
    # Although this method is availlable for finitely presented group
    # we use the classical group implementation (by performace reason) to get
    # the permutation_group.
    # ----------------------------------------------------------------------------------
    def as_permutation_group(self, use_classical=True):
        r"""
        This method returns a permutation group isomporhic to ``self`` together
        with group isomorphism from ``self`` as a conversion.

        INPUT (all arguments are optional by keyword):

        - ``use_classical`` -- (boolean, default True) by default the permutation
          group is calculated via the attached classical matrix group, since this
          results in a smaller degree. If set to False the permutation group will
          be calculated using ``self`` (as finitely presented group).

        OUTPUT:

        An instance of class :class:`PermutationGroup_generic` together with a group homomorphism
        from ``self`` registerd as a conversion.

        EXAMPLES::

            sage: C3 = CubicBraidGroup(3)
            sage: PC3 = C3.as_permutation_group()
            sage: C3.is_isomorphic(PC3)
            #I  Forcing finiteness test
            True
            sage: PC3.degree()
            8
            sage: c = C3([2,1-2])
            sage: C3(PC3(c)) == c
            True
        """
        if not self.is_finite():
            raise ValueError('cubic braid group is infinite!')

        if use_classical:
            CG = self.as_classical_group()
            from sage.groups.perm_gps.permgroup import PermutationGroup_generic
            if isinstance(CG, PermutationGroup_generic):
                return CG
            CGM = CG.as_matrix_group()
            PG = CGM.as_permutation_group()
            img_gens = [PG(CGM(CG(gen))) for gen in self.gens()]
        else:
            PG = super(CubicBraidGroup_class, self).as_permutation_group()
            img_gens = PG.gens()

        img_gens = [PG(gen) for gen in img_gens]
        hom_to_perm = self.hom(img_gens)
        PG.register_conversion(hom_to_perm)
        return PG


    # ----------------------------------------------------------------------------------
    # as_classical_group
    # ----------------------------------------------------------------------------------
    def as_classical_group(self, embedded=False):
        r"""
        Creates an isomorphic image of ``self`` as a classical group according
        to the construction given by Coxeter resp. Assion.

        INPUT (optional keyword):

        - ``embedded`` -- boolean (default = False). This boolean does effect the
          cases of Assion groups when they are realized as projective groups, only.
          More precisely: if ``self`` is of ``cbg_type CubicBraidGroupType.AssionS``
          (for example) and the number of strands ``n`` is even, than its classical group
          is a subgroup of ``PSp(n,3)`` (being centralized by the element
          ``self.centralizing_element(projective=True))``. By default this group will be
          given. Setting ``embedded = True`` the classical realization is given as
          subgroup of its classical enlargement with one more strand (in this
          case as subgroup of ``Sp(n,3))``.

        OUTPUT:

        Depending on the type of ``self`` and the number of strands an instance of ``Sp(n-1,3)``,
        ``GU(n-1,2)``, subgroup of ``PSp(n,3), PGU(n,2)`` or a subgroup of ``GU(n-1, UCF)``
        (``cbg_type CubicBraidGroupType.Coxeter``) with respect to a certain hermitian form
        attached to the burau representation (used by Coxeter and Squier). Here ``UCF`` stands
        for the universal cyclotomic field.

        EXAMPLES::

            sage: U3 = AssionGroupU(3)
            sage: U3Cl = U3.as_classical_group(); U3Cl
            Subgroup generated by [(1,7,6)(3,19,14)(4,15,10)(5,11,18)(12,16,20), (1,12,13)(2,15,19)(4,9,14)(5,18,8)(6,21,16)] of (The projective general unitary group of degree 3 over Finite Field of size 2)
            sage: U3Clemb = U3.as_classical_group(embedded=True); U3Clemb
            Subgroup with 2 generators (
            [0 0 a]  [a + 1     a     a]
            [0 1 0]  [    a a + 1     a]
            [a 0 a], [    a     a a + 1]
            ) of General Unitary Group of degree 3 over Finite Field in a of size 2^2
            sage: u = U3([-2,1,-2,1]); u
            (u1^-1*u0)^2
            sage: uCl = U3Cl(u); uCl
            (1,16)(2,9)(3,10)(4,19)(6,12)(7,20)(13,21)(14,15)
            sage: uCle = U3Clemb(u); uCle
            [a + 1 a + 1     1]
            [a + 1     0     a]
            [   1     a     a]
            sage: U3(uCl) == u
            True
            sage: U3(uCle) == u
            True
            sage: U4 = AssionGroupU(4)
            sage: U4Cl = U4.as_classical_group(); U4Cl
            General Unitary Group of degree 3 over Finite Field in a of size 2^2
            sage: U3Clemb.ambient() == U4Cl
            True
            sage: C4 = CubicBraidGroup(4)
            sage: C4Cl = C4.as_classical_group(); C4Cl
            Subgroup with 3 generators (
            [  E(3)^2        0        0]  [       1 -E(12)^7        0]
            [-E(12)^7        1        0]  [       0   E(3)^2        0]
            [       0        0        1], [       0 -E(12)^7        1],
            <BLANKLINE>
            [       1        0        0]
            [       0        1 -E(12)^7]
            [       0        0   E(3)^2]
            ) of General Unitary Group of degree 3 over Universal Cyclotomic Field with respect to positive definite hermitian form
            [-E(12)^7 + E(12)^11                  -1                   0]
            [                 -1 -E(12)^7 + E(12)^11                  -1]
            [                  0                  -1 -E(12)^7 + E(12)^11]
        """
        # -------------------------------------------------------------------------------
        # create the classical group if not already done
        # -------------------------------------------------------------------------------
        if self._classical_group is None:
            if embedded and self._classical_embedding is None:
                # this is separated to avoid unneccessary (for larger group exhaustive) calculations
                self._create_classical_realization(just_embedded=True)
            else:
                self._create_classical_realization()

        if embedded and self._classical_embedding is not None:
            # ----------------------------------------------------------------------------------------
            # there is a difference between self._classical_group and self._classical_embedding
            # only in the cases where self.strands() divides by 2 (AssionGroupS) resp. 3
            # (AssionGroupU). In this case the embedding is the subgroup of the classical group
            # of one strand more (self.strands() +1) generated by the first self.strands() -1
            # generators
            # ----------------------------------------------------------------------------------------
           return self._classical_embedding
        elif self._classical_group is not None:
           return self._classical_group
        else:
           raise ValueError("no classical embedding defined")


    # ----------------------------------------------------------------------------------
    # as_refection_group
    # ----------------------------------------------------------------------------------
    def as_reflection_group(self):
        r"""
        Creates an isomorphic image of ``self`` as irreducible complex reflection group.
        This is possible only for the finite cubic braid groups of ``cbg_type
        CubicBraidGroupType.Coxeter``.

        This method uses the sage implementation of reflection group via the gap3 CHEVIE
        package. To use this method you must have gap3 together with CHEVIE installed!

        OUTPUT:

        An instance of the class :class:`IrreducibleComplexReflectionGroup` together with
        a group isomorphism from ``self`` registered as a conversion.


        EXAMPLES::

           sage: C3.<c1,c2> = CubicBraidGroup(3)           # optional - gap3
           sage: R3 = C3.as_reflection_group(); R3         # optional - gap3
           Irreducible complex reflection group of rank 2 and type ST4
           sage: R3.cartan_matrix()                        # optional - gap3
           [-2*E(3) - E(3)^2           E(3)^2]
           [        -E(3)^2 -2*E(3) - E(3)^2]
           sage: R3.simple_roots()                         # optional - gap3
           Finite family {1: (0, -2*E(3) - E(3)^2), 2: (2*E(3)^2, E(3)^2)}
           sage: R3.simple_coroots()                       # optional - gap3
           Finite family {1: (0, 1), 2: (1/3*E(3) - 1/3*E(3)^2, 1/3*E(3) - 1/3*E(3)^2)}

       Conversion maps::

           sage: r = R3.an_element()                       # optional - gap3
           sage: cr = C3(r); cr                            # optional - gap3
           c1*c2
           sage: mr = r.matrix(); mr                       # optional - gap3
           [ 1/3*E(3) - 1/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
           [-2/3*E(3) + 2/3*E(3)^2  2/3*E(3) + 1/3*E(3)^2]
           sage: C3Cl = C3.as_classical_group()            # optional - gap3
           sage: C3Cl(cr)                                  # optional - gap3
           [ E(3)^2    -E(4)]
           [-E(12)^7        0]

       The refelection groups can also be viewed as subgroups of unitary groups
       over the universal cyclotomic field. Note that the unitary group corresponding
       to the reflection group is isomorphic but different from the classical group due
       to different hermitian forms for the unitary groups they live in::

           sage: C4 = CubicBraidGroup(4)                   # optional - gap3
           sage: R4 = C4.as_reflection_group()             # optional - gap3
           sage: R4.invariant_form()                       # optional - gap3
           [1 0 0]
           [0 1 0]
           [0 0 1]
           sage: _ == C4.classical_invariant_form()        # optional - gap3
           False
        """

        # -------------------------------------------------------------------------------
        # the reflection groups are called according to the Shephard-Todd classification:
        #    2 strands -> G(2,1,1)
        #    3 strands -> G4
        #    4 strands -> G25
        #    5 strands -> G32
        # -------------------------------------------------------------------------------
        from sage.combinat.root_system.reflection_group_real import is_chevie_available
        if not is_chevie_available():
            raise ImportError("the GAP3 package 'chevie' is needed to obtain the corresponding reflection groups")

        if self._cbg_type  != CubicBraidGroupType.Coxeter or self.strands() > 5  or self.strands() < 2:
            raise ValueError("no reflection group defined")

        # -------------------------------------------------------------------------------
        # define reflection group assiossiated to self
        # -------------------------------------------------------------------------------
        reflection_group = None

        from sage.combinat.root_system.reflection_group_real import ReflectionGroup

        if   self.strands() == 2:
            reflection_group = ReflectionGroup([2 ,1 ,1])
        elif self.strands() == 3:
            reflection_group = ReflectionGroup(4)
        elif self.strands() == 4:
            reflection_group = ReflectionGroup(25)
        elif self.strands() == 5:
            reflection_group = ReflectionGroup(32)

        hom_to_refl = self.hom(reflection_group.gens())
        reflection_group.register_conversion(hom_to_refl)
        return reflection_group


    # ----------------------------------------------------------------------------------
    # classical invariant form returns the invariant form of the classical realization
    # ----------------------------------------------------------------------------------
    def classical_invariant_form(self):
        r"""
        Returns the invariant form of the classical realization of ``self``.

        OUTPUT:

        A square matrix of dimension according to the space the classical realization is
        operating on. In the case of the full cubic braid groups and of the assion groups
        of ``cbg_type CubicBraidGroupType.AssionU`` the matrix is hermitian. In the case of
        the Assion groups of ``cbg_type CubicBraidGroupType.AssionS`` it is alternating.
        Note that the invariant form of the full cubic braid group on more than 5 strands
        is degenerated (causing the group to be infinite).

        In the case of Assion groups having projective classical groups the invariant form
        corresponds to the ambient group of its classical embedding.

        EXAMPLES::

            sage: S3 = AssionGroupS(3)
            sage: S3.classical_invariant_form()
            [0 1]
            [2 0]
            sage: S4 = AssionGroupS(4)
            sage: S4.classical_invariant_form()
            [0 0 0 1]
            [0 0 1 0]
            [0 2 0 0]
            [2 0 0 0]
            sage: S5 = AssionGroupS(5)
            sage: S4.classical_invariant_form() == S5.classical_invariant_form()
            True
            sage: U4 = AssionGroupU(4)
            sage: U4.classical_invariant_form()
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: C5 = CubicBraidGroup(5)
            sage: C5.classical_invariant_form()
            [-E(12)^7 + E(12)^11                  -1                   0                   0]
            [                -1 -E(12)^7 + E(12)^11                  -1                   0]
            [                 0                  -1 -E(12)^7 + E(12)^11                  -1]
            [                 0                   0                  -1 -E(12)^7 + E(12)^11]
            sage: _.is_singular()
            False
            sage: C6 = CubicBraidGroup(6)
            sage: C6.classical_invariant_form().is_singular()
            True
        """
        # -------------------------------------------------------------------------------
        # create the classical_invariant_form if not already done
        # -------------------------------------------------------------------------------
        if self._classical_invariant_form is None:
            self._create_classical_realization()

        if self._classical_invariant_form is None:
            raise ValueError("no classical invariant form defined!")

        return self._classical_invariant_form


    # ----------------------------------------------------------------------------------
    # centralizing element in the classical symplectic resp. unitary group
    # ----------------------------------------------------------------------------------
    def centralizing_element(self, embedded=False):
        r"""
        Returns the centralizing element defined by the work of Assion
        (Hilfssatz 1.1.3 and 1.2.3).

        INPUT (optional):

        - ``embedded`` -- boolean (default = False). This boolean just effects
          the cases of Assion groups when they are realized as projective groups.
          More precisely: if ``self`` is of ``cbg_type CubicBraidGroupType.AssionS``
          (for example) and the number of strands ``n`` is even, than its classical
          group is a subgroup of ``PSp(n,3)`` being centralized by the element return
          for option ``embedded=False``. Otherwise the image of this element inside
          the embedded classical group will be returned (see option embedded of
          :meth:`classical_group`)!

        OUTPUT:

        Depending on the optional keyword a permutation as an element of ``PSp(n,3)``
        (type S) or ``PGU(n,2)`` (type U) for ``n = 0 mod 2`` (type S) reps. ``n = 0 mod 3``
        (type U) is returned. Elsewise, the centralizing element is a matrix
        belonging to ``Sp(n,3)`` reps. ``GU(n,2)``.

        EXAMPLES::

            sage: U3 = AssionGroupU(3);  U3
            Assion group on 3 strands of type U
            sage: U3Cl = U3.as_classical_group(); U3Cl
            Subgroup generated by [(1,7,6)(3,19,14)(4,15,10)(5,11,18)(12,16,20), (1,12,13)(2,15,19)(4,9,14)(5,18,8)(6,21,16)] of (The projective general unitary group of degree 3 over Finite Field of size 2)
            sage: c = U3.centralizing_element(); c
            (1,16)(2,9)(3,10)(4,19)(6,12)(7,20)(13,21)(14,15)
            sage: c in U3Cl
            True
            sage: P = U3Cl.ambient_group()
            sage: P.centralizer(c) == U3Cl
            True

        embedded Version::

            sage: cm = U3.centralizing_element(embedded=True); cm
            [a + 1 a + 1     1]
            [a + 1     0     a]
            [   1     a     a]
            sage: U4 = AssionGroupU(4)
            sage: U4Cl = U4.as_classical_group()
            sage: cm in U4Cl
            True
            sage: [cm * U4Cl(g) == U4Cl(g) * cm for g in U4.gens()]
            [True, True, False]
        """

        # -------------------------------------------------------------------------------
        # create the centralizing elements if not already done
        # -------------------------------------------------------------------------------
        if self._centralizing_matrix is None:
            self._create_classical_realization()

        if self._centralizing_matrix is None:
            raise ValueError("no centralizing element defined")
        else:
            if embedded == True or self._centralizing_element is None:
                return self._centralizing_matrix
            else:
                return self._centralizing_element


    # ----------------------------------------------------------------------------------
    # calculating the order by formula
    # ----------------------------------------------------------------------------------
    def order(self):
        r"""
        To avoid long waittime on calculations the order will be obtained
        using the classical realization.

        OUTPUT:

        Cardinality of the group as Integer or infinity.

        EXAMPLES::

            sage: S15 = AssionGroupS(15)
            sage: S15.order()
            109777561863482259035023554842176139436811616256000
            sage: C6 = CubicBraidGroup(6)
            sage: C6.order()
            +Infinity
        """
        from sage.rings.infinity import infinity
        n = self.strands()

        if self._cbg_type == CubicBraidGroupType.Coxeter and n > 5:
            order = infinity
        else:
            order = self.as_classical_group(embedded=True).order()

        return order


    cardinality = order


    def is_finite(self):
        r"""
        Method from :class:`GroupMixinLibGAP` overwriten because of performance reason.
        """
        from sage.rings.infinity import infinity
        return not self.order() is infinity


    # ----------------------------------------------------------------------------------
    # creating a CubicBraidGroup as subgroup of self on less strands
    # ----------------------------------------------------------------------------------
    def cubic_braid_subgroup(self, nstrands = None):
        r"""
        Creates a cubic braid group as subgroup of ``self`` on the first ``nstrands`` strands.

        INPUT:

        - ``nstrands`` -- integer > 0 and < ``self.strands()`` giving the number of strands
          for the subgroup. The default is one strand less than ``self`` has.

        OUTPUT:

        An instance of this class realizing the subgroup.

        .. NOTE::

            Since ``self`` is inherited from :class:`UniqueRepresentation` the obtained instance
            is identical to other instances created with the same arguments (see example
            below). The ambient group corresponds to the last call of this method.

        EXAMPLES::

            sage: U5 = AssionGroupU(5)
            sage: U3s = U5.cubic_braid_subgroup(3)
            sage: u1, u2 = U3s.gens()
            sage: u1 in U5
            False
            sage: U5(u1) in U5.gens()
            True
            sage: U3s is AssionGroupU(3)
            True
            sage: U3s.ambient() == U5
            True
        """
        if nstrands is None:
            nstrands = self.strands() -1

        n = self.strands()

        nstrands = Integer(nstrands)

        if nstrands >= n or nstrands <= 0:
            raise ValueError("nstrands must be positive and less than %s" %(self.strands()))

        gens = self.gens()
        gens_red = tuple([gens[i] for i in range(nstrands -1)])
        subgrp = CubicBraidGroup(names=gens_red, cbg_type=self._cbg_type)
        subgrp._ambient = self
        return subgrp
