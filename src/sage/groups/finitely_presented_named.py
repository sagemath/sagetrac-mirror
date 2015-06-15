r"""
Named Finitely Presented Groups

Construct groups of small order and "named" groups as quotients of free groups. These
groups are available through tab completion by typing ``groups.presentation.<tab>``
or by importing the required methods. Tab completion is made available through
Sage's :ref:`group catalog <sage.groups.groups_catalog>`. Some examples are engineered
from entries in [THOMAS-WOODS]_.

Groups available as finite presentations:

- Alternating group, `A_n` of order `n!/2` --
  :func:`groups.presentation.Alternating <sage.groups.finitely_presented_named.AlternatingPresentation>`

- Cyclic group, `C_n` of order `n` --
  :func:`groups.presentation.Cyclic <sage.groups.finitely_presented_named.CyclicPresentation>`

- Dicyclic group, nonabelian groups of order `4n` with a unique element of
  order 2 --
  :func:`groups.presentation.DiCyclic <sage.groups.finitely_presented_named.DiCyclicPresentation>`

- Dihedral group, `D_n` of order `2n` --
  :func:`groups.presentation.Dihedral <sage.groups.finitely_presented_named.DihedralPresentation>`

- Finitely generated abelian group, `\ZZ_{n_1} \times \ZZ_{n_2} \times \cdots \times \ZZ_{n_k}` --
  :func:`groups.presentation.FGAbelian <sage.groups.finitely_presented_named.FinitelyGeneratedAbelianPresentation>`

- Klein four group, `C_2 \times C_2` --
  :func:`groups.presentation.KleinFour <sage.groups.finitely_presented_named.KleinFourPresentation>`

- Quaternion group of order 8 --
  :func:`groups.presentation.Quaternion <sage.groups.finitely_presented_named.QuaternionPresentation>`

- Symmetric group, `S_n` of order `n!` --
  :func:`groups.presentation.Symmetric <sage.groups.finitely_presented_named.SymmetricPresentation>`

AUTHORS:

- Davis Shurbert (2013-06-21): initial version

EXAMPLES::

    sage: groups.presentation.Cyclic(4)
    Finitely presented group < a | a^4 >

You can also import the desired functions::

    sage: from sage.groups.finitely_presented_named import CyclicPresentation
    sage: CyclicPresentation(4)
    Finitely presented group < a | a^4 >
"""
#*****************************************************************************
#       Copyright (C) 2013 Davis Shurbert <davis.sprout@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all      import Integer
from sage.groups.free_group import FreeGroup
from sage.groups.finitely_presented import FinitelyPresentedGroup
from sage.libs.gap.libgap import libgap
from sage.matrix.constructor import diagonal_matrix
from sage.modules.fg_pid.fgp_module import FGP_Module
from sage.rings.integer_ring import ZZ
from sage.sets.set import Set

def CyclicPresentation(n):
    r"""
    Build cyclic group of order `n` as a finitely presented group.

    INPUT:

    - ``n`` -- The order of the cyclic presentation to be returned.

    OUTPUT:

    The cyclic group of order `n` as finite presentation.

    EXAMPLES::

        sage: groups.presentation.Cyclic(10)
        Finitely presented group < a | a^10 >
        sage: n = 8; C = groups.presentation.Cyclic(n)
        sage: C.as_permutation_group().is_isomorphic(CyclicPermutationGroup(n))
        True

    TESTS::

        sage: groups.presentation.Cyclic(0)
        Traceback (most recent call last):
        ...
        ValueError: finitely presented group order must be positive
    """
    n = Integer(n)
    if n < 1:
        raise ValueError('finitely presented group order must be positive')
    F = FreeGroup( 'a' )
    rls = F([1])**n,
    return FinitelyPresentedGroup( F, rls )

def FinitelyGeneratedAbelianPresentation(int_list):
    r"""
    Return canonical presentation of finitely generated abelian group.

    INPUT:

    - ``int_list`` -- List of integers defining the group to be returned, the defining list
      is reduced to the invariants of the input list before generating the corresponding
      group.

    OUTPUT:

    Finitely generated abelian group, `\ZZ_{n_1} \times \ZZ_{n_2} \times \cdots \times \ZZ_{n_k}`
    as a finite presentation, where `n_i` forms the invariants of the input list.

    EXAMPLES::

        sage: groups.presentation.FGAbelian([2,2])
        Finitely presented group < a, b | a^2, b^2, a^-1*b^-1*a*b >
        sage: groups.presentation.FGAbelian([2,3])
        Finitely presented group < a | a^6 >
        sage: groups.presentation.FGAbelian([2,4])
        Finitely presented group < a, b | a^2, b^4, a^-1*b^-1*a*b >

    You can create free abelian groups::

        sage: groups.presentation.FGAbelian([0])
        Finitely presented group < a |  >
        sage: groups.presentation.FGAbelian([0,0])
        Finitely presented group < a, b | a^-1*b^-1*a*b >
        sage: groups.presentation.FGAbelian([0,0,0])
        Finitely presented group < a, b, c | a^-1*c^-1*a*c, a^-1*b^-1*a*b, c^-1*b^-1*c*b >

    And various infinite abelian groups::

        sage: groups.presentation.FGAbelian([0,2])
        Finitely presented group < a, b | a^2, a^-1*b^-1*a*b >
        sage: groups.presentation.FGAbelian([0,2,2])
        Finitely presented group < a, b, c | a^2, b^2, a^-1*c^-1*a*c, a^-1*b^-1*a*b, c^-1*b^-1*c*b >

    Outputs are reduced to minimal generators and relations::

        sage: groups.presentation.FGAbelian([3,5,2,7,3])
        Finitely presented group < a, b | a^3, b^210, a^-1*b^-1*a*b >
        sage: groups.presentation.FGAbelian([3,210])
        Finitely presented group < a, b | a^3, b^210, a^-1*b^-1*a*b >

    The trivial group is an acceptable output::

        sage: groups.presentation.FGAbelian([])
        Finitely presented group <  |  >
        sage: groups.presentation.FGAbelian([1])
        Finitely presented group <  |  >
        sage: groups.presentation.FGAbelian([1,1,1,1,1,1,1,1,1,1])
        Finitely presented group <  |  >

    Input list must consist of positive integers::

        sage: groups.presentation.FGAbelian([2,6,3,9,-4])
        Traceback (most recent call last):
        ...
        ValueError: input list must contain nonnegative entries
        sage: groups.presentation.FGAbelian([2,'a',4])
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 'a' to an integer

    TESTS::

        sage: ag = groups.presentation.FGAbelian([2,2])
        sage: ag.as_permutation_group().is_isomorphic(groups.permutation.KleinFour())
        True
        sage: G = groups.presentation.FGAbelian([2,4,8])
        sage: C2 = CyclicPermutationGroup(2)
        sage: C4 = CyclicPermutationGroup(4)
        sage: C8 = CyclicPermutationGroup(8)
        sage: gg = (C2.direct_product(C4)[0]).direct_product(C8)[0]
        sage: gg.is_isomorphic(G.as_permutation_group())
        True
        sage: all([groups.presentation.FGAbelian([i]).as_permutation_group().is_isomorphic(groups.presentation.Cyclic(i).as_permutation_group()) for i in [2..35]])
        True
    """
    from sage.groups.free_group import _lexi_gen
    check_ls = [Integer(x) for x in int_list if Integer(x) >= 0]
    if len(check_ls) != len(int_list):
        raise ValueError('input list must contain nonnegative entries')

    col_sp = diagonal_matrix(int_list).column_space()
    invariants = FGP_Module(ZZ**(len(int_list)), col_sp).invariants()
    name_gen = _lexi_gen()
    F = FreeGroup([next(name_gen) for i in invariants])
    ret_rls = [F([i+1])**invariants[i] for i in range(len(invariants)) if invariants[i]!=0]

    # Build commutator relations
    gen_pairs = list(Set(F.gens()).subsets(2))
    ret_rls = ret_rls + [x[0]**(-1)*x[1]**(-1)*x[0]*x[1] for x in gen_pairs]
    return FinitelyPresentedGroup(F, tuple(ret_rls))

def DihedralPresentation(n):
    r"""
    Build the Dihedral group of order `2n` as a finitely presented group.

    INPUT:

    - ``n`` -- The size of the set that `D_n` is acting on.

    OUTPUT:

    Dihedral group of order `2n`.

    EXAMPLES::

        sage: D = groups.presentation.Dihedral(7); D
        Finitely presented group < a, b | a^7, b^2, (a*b)^2 >
        sage: D.as_permutation_group().is_isomorphic(DihedralGroup(7))
        True

    TESTS::

        sage: n = 9
        sage: D = groups.presentation.Dihedral(n)
        sage: D.ngens() == 2
        True
        sage: groups.presentation.Dihedral(0)
        Traceback (most recent call last):
        ...
        ValueError: finitely presented group order must be positive
    """
    n = Integer( n )
    if n < 1:
        raise ValueError('finitely presented group order must be positive')
    F = FreeGroup([ 'a', 'b' ])
    rls = F([1])**n, F([2])**2, (F([1])*F([2]))**2
    return FinitelyPresentedGroup( F, rls )

def DiCyclicPresentation(n):
    r"""
    Build the dicyclic group of order `4n`, for `n \geq 2`, as a finitely
    presented group.

    INPUT:

    - ``n`` -- positive integer, 2 or greater, determining the order of
      the group (`4n`).

    OUTPUT:

    The dicyclic group of order `4n` is defined by the presentation

    .. MATH::

        \langle a, x \mid a^{2n}=1, x^{2}=a^{n}, x^{-1}ax=a^{-1} \rangle

    .. NOTE::

        This group is also available as a permutation group via
        :class:`groups.permutation.DiCyclic <sage.groups.perm_gps.permgroup_named.DiCyclicGroup>`.

    EXAMPLES::

        sage: D = groups.presentation.DiCyclic(9); D
        Finitely presented group < a, b | a^18, b^2*a^-9, b^-1*a*b*a >
        sage: D.as_permutation_group().is_isomorphic(groups.permutation.DiCyclic(9))
        True

    TESTS::

        sage: Q = groups.presentation.DiCyclic(2)
        sage: Q.as_permutation_group().is_isomorphic(QuaternionGroup())
        True
        sage: all([groups.presentation.DiCyclic(i).as_permutation_group(
        ....: ).is_isomorphic(groups.permutation.DiCyclic(i)) for i in [5,8,12,2^5]])
        True
        sage: groups.presentation.DiCyclic(1)
        Traceback (most recent call last):
        ...
        ValueError: input integer must be greater than 1
    """
    n = Integer(n)
    if n < 2:
        raise ValueError('input integer must be greater than 1')

    F = FreeGroup(['a','b'])
    rls =  F([1])**(2*n), F([2,2])*F([-1])**n, F([-2,1,2,1])
    return FinitelyPresentedGroup(F, rls)

def SymmetricPresentation(n):
    r"""
    Build the Symmetric group of order `n!` as a finitely presented group.

    INPUT:

    - ``n`` -- The size of the underlying set of arbitrary symbols being acted
      on by the Symmetric group of order `n!`.

    OUTPUT:

    Symmetric group as a finite presentation, implementation uses GAP to find an
    isomorphism from a permutation representation to a finitely presented group
    representation. Due to this fact, the exact output presentation may not be
    the same for every method call on a constant ``n``.

    EXAMPLES::

        sage: S4 = groups.presentation.Symmetric(4)
        sage: S4.as_permutation_group().is_isomorphic(SymmetricGroup(4))
        True

    TESTS::

        sage: S = [groups.presentation.Symmetric(i) for i in range(1,4)]; S[0].order()
        1
        sage: S[1].order(), S[2].as_permutation_group().is_isomorphic(DihedralGroup(3))
        (2, True)
        sage: S5 = groups.presentation.Symmetric(5)
        sage: perm_S5 = S5.as_permutation_group(); perm_S5.is_isomorphic(SymmetricGroup(5))
        True
        sage: groups.presentation.Symmetric(8).order()
        40320
    """
    from sage.groups.perm_gps.permgroup_named import SymmetricGroup
    from sage.groups.free_group import _lexi_gen

    n = Integer(n)
    perm_rep = SymmetricGroup(n)
    GAP_fp_rep = libgap.Image(libgap.IsomorphismFpGroupByGenerators(perm_rep, perm_rep.gens()))
    image_gens = GAP_fp_rep.FreeGeneratorsOfFpGroup()
    name_itr = _lexi_gen() # Python generator object for variable names
    F = FreeGroup([next(name_itr) for x in perm_rep.gens()])
    ret_rls = tuple([F(rel_word.TietzeWordAbstractWord(image_gens).sage())
                for rel_word in GAP_fp_rep.RelatorsOfFpGroup()])
    return FinitelyPresentedGroup(F,ret_rls)

def QuaternionPresentation():
    r"""
    Build the Quaternion group of order 8 as a finitely presented group.

    OUTPUT:

    Quaternion group as a finite presentation.

    EXAMPLES::

        sage: Q = groups.presentation.Quaternion(); Q
        Finitely presented group < a, b | a^4, b^2*a^-2, a*b*a*b^-1 >
        sage: Q.as_permutation_group().is_isomorphic(QuaternionGroup())
        True

    TESTS::

        sage: Q = groups.presentation.Quaternion()
        sage: Q.order(), Q.is_abelian()
        (8, False)
        sage: Q.is_isomorphic(groups.presentation.DiCyclic(2))
        True
    """
    F = FreeGroup(['a','b'])
    rls = F([1])**4, F([2,2,-1,-1]), F([1,2,1,-2])
    return FinitelyPresentedGroup(F, rls)

def AlternatingPresentation(n):
    r"""
    Build the Alternating group of order `n!/2` as a finitely presented group.

    INPUT:

    - ``n`` -- The size of the underlying set of arbitrary symbols being acted
      on by the Alternating group of order `n!/2`.

    OUTPUT:

    Alternating group as a finite presentation, implementation uses GAP to find an
    isomorphism from a permutation representation to a finitely presented group
    representation. Due to this fact, the exact output presentation may not be
    the same for every method call on a constant ``n``.

    EXAMPLES::

        sage: A6 = groups.presentation.Alternating(6)
        sage: A6.as_permutation_group().is_isomorphic(AlternatingGroup(6)), A6.order()
        (True, 360)

    TESTS::

        sage: #even permutation test..
        sage: A1 = groups.presentation.Alternating(1); A2 = groups.presentation.Alternating(2)
        sage: A1.is_isomorphic(A2), A1.order()
        (True, 1)
        sage: A3 = groups.presentation.Alternating(3); A3.order(), A3.as_permutation_group().is_cyclic()
        (3, True)
        sage: A8 = groups.presentation.Alternating(8); A8.order()
        20160
    """
    from sage.groups.perm_gps.permgroup_named import AlternatingGroup
    from sage.groups.free_group import _lexi_gen

    n = Integer(n)
    perm_rep = AlternatingGroup(n)
    GAP_fp_rep = libgap.Image(libgap.IsomorphismFpGroupByGenerators(perm_rep, perm_rep.gens()))
    image_gens = GAP_fp_rep.FreeGeneratorsOfFpGroup()
    name_itr = _lexi_gen() # Python generator object for variable names
    F = FreeGroup([next(name_itr) for x in perm_rep.gens()])
    ret_rls = tuple([F(rel_word.TietzeWordAbstractWord(image_gens).sage())
                for rel_word in GAP_fp_rep.RelatorsOfFpGroup()])
    return FinitelyPresentedGroup(F,ret_rls)

def KleinFourPresentation():
    r"""
    Build the Klein group of order `4` as a finitely presented group.

    OUTPUT:

    Klein four group (`C_2 \times C_2`) as a finitely presented group.

    EXAMPLES::

        sage: K = groups.presentation.KleinFour(); K
        Finitely presented group < a, b | a^2, b^2, a^-1*b^-1*a*b >
    """
    F = FreeGroup(['a','b'])
    rls = F([1])**2, F([2])**2, F([-1])*F([-2])*F([1])*F([2])
    return FinitelyPresentedGroup(F, rls)
