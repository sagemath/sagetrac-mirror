# -*- coding: utf-8 -*-
r"""
Creation of modular symbols spaces

EXAMPLES: We create a space and output its category.

::

    sage: C = HeckeModules(RationalField()); C
    Category of Hecke modules over Rational Field
    sage: M = ModularSymbols(11)
    sage: M.category()
    Category of Hecke modules over Rational Field
    sage: M in C
    True

We create a space compute the charpoly, then compute the same but
over a bigger field. In each case we also decompose the space using
`T_2`.

::

    sage: M = ModularSymbols(23,2,base_ring=QQ)
    sage: M.T(2).charpoly('x').factor()
    (x - 3) * (x^2 + x - 1)^2
    sage: M.decomposition(2)
    [
    Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field,
    Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
    ]

::

    sage: M = ModularSymbols(23,2,base_ring=QuadraticField(5, 'sqrt5'))
    sage: M.T(2).charpoly('x').factor()
    (x - 3) * (x - 1/2*sqrt5 + 1/2)^2 * (x + 1/2*sqrt5 + 1/2)^2
    sage: M.decomposition(2)
    [
    Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Number Field in sqrt5 with defining polynomial x^2 - 5,
    Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Number Field in sqrt5 with defining polynomial x^2 - 5,
    Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Number Field in sqrt5 with defining polynomial x^2 - 5
    ]

We compute some Hecke operators and do a consistency check::

    sage: m = ModularSymbols(39, 2)
    sage: t2 = m.T(2); t5 = m.T(5)
    sage: t2*t5 - t5*t2 == 0
    True

This tests the bug reported in :trac:`1220`::

    sage: G = GammaH(36, [13, 19])
    sage: G.modular_symbols()
    Modular Symbols space of dimension 13 for Congruence Subgroup Gamma_H(36) with H generated by [13, 19] of weight 2 with sign 0 and over Rational Field
    sage: G.modular_symbols().cuspidal_subspace()
    Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 13 for Congruence Subgroup Gamma_H(36) with H generated by [13, 19] of weight 2 with sign 0 and over Rational Field

This test catches a tricky corner case for spaces with character::

    sage: ModularSymbols(DirichletGroup(20).1**3, weight=3, sign=1).cuspidal_subspace()
    Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 6 and level 20, weight 3, character [1, -zeta4], sign 1, over Cyclotomic Field of order 4 and degree 2
"""

#*****************************************************************************
#       Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

import weakref

import sage.modular.arithgroup.all as arithgroup
import sage.modular.dirichlet as dirichlet
import sage.rings.rational_field as rational_field
import sage.rings.all as rings


def canonical_parameters(group, weight, sign, base_ring):
    """
    Return the canonically normalized parameters associated to a choice
    of group, weight, sign, and base_ring. That is, normalize each of
    these to be of the correct type, perform all appropriate type
    checking, etc.

    EXAMPLES::

        sage: p1 = sage.modular.modsym.modsym.canonical_parameters(5,int(2),1,QQ) ; p1
        (Congruence Subgroup Gamma0(5), 2, 1, Rational Field)
        sage: p2 = sage.modular.modsym.modsym.canonical_parameters(Gamma0(5),2,1,QQ) ; p2
        (Congruence Subgroup Gamma0(5), 2, 1, Rational Field)
        sage: p1 == p2
        True
        sage: type(p1[1])
        <type 'sage.rings.integer.Integer'>
    """
    sign = rings.Integer(sign)
    if not (sign in [-1,0,1]):
        raise ValueError("sign must be -1, 0, or 1")

    weight = rings.Integer(weight)
    if weight <= 1:
        raise ValueError("the weight must be at least 2")

    if isinstance(group, (int, rings.Integer)):
        group = arithgroup.Gamma0(group)
    elif isinstance(group, dirichlet.DirichletCharacter):
        if group.is_trivial():
            group = arithgroup.Gamma0(group.modulus())
        else:
            eps = group.minimize_base_ring()
            group = (eps, eps.parent())
            if base_ring is None:
                base_ring = eps.base_ring()

    if base_ring is None:
        base_ring = rational_field.RationalField()

    if not isinstance(base_ring, rings.CommutativeRing):
        raise TypeError("base_ring (=%s) must be a commutative ring"%base_ring)

    if not base_ring.is_field():
        raise TypeError("(currently) base_ring (=%s) must be a field"%base_ring)

    return group, weight, sign, base_ring

_cache = {}

def ModularSymbols_clear_cache():
    """
    Clear the global cache of modular symbols spaces.

    EXAMPLES::

        sage: sage.modular.modsym.modsym.ModularSymbols_clear_cache()
        sage: sage.modular.modsym.modsym._cache.keys()
        []
        sage: M = ModularSymbols(6,2)
        sage: sage.modular.modsym.modsym._cache.keys()
        [(Congruence Subgroup Gamma0(6), 2, 0, Rational Field)]
        sage: sage.modular.modsym.modsym.ModularSymbols_clear_cache()
        sage: sage.modular.modsym.modsym._cache.keys()
        []

    TESTS:

    Make sure :trac:`10548` is fixed::

        sage: import gc
        sage: m=ModularSymbols(Gamma1(29))
        sage: m=[]
        sage: ModularSymbols_clear_cache()
        sage: gc.collect() # random
        3422
        sage: a=[x for x in gc.get_objects() if isinstance(x,sage.modular.modsym.ambient.ModularSymbolsAmbient_wtk_g1)]
        sage: a
        []

    """
    global _cache
    _cache = {}

def ModularSymbols(group  = 1,
                   weight = 2,
                   sign   = 0,
                   base_ring = None,
                   use_cache = True,
                   custom_init=None):
    r"""
    Create an ambient space of modular symbols.

    INPUT:

    - ``group`` - A congruence subgroup or a Dirichlet character eps.
    - ``weight`` - int, the weight, which must be = 2.
    - ``sign`` - int, The sign of the involution on modular symbols
      induced by complex conjugation. The default is 0, which means
      "no sign", i.e., take the whole space.
    - ``base_ring`` - the base ring. Defaults to `\QQ` if no character
      is given, or to the minimal extension of `\QQ` containing the
      values of the character.
    - ``custom_init`` - a function that is called with self as input
      before any computations are done using self; this could be used
      to set a custom modular symbols presentation.  If self is
      already in the cache and use_cache=True, then this function is
      not called.

    EXAMPLES: First we create some spaces with trivial character::

        sage: ModularSymbols(Gamma0(11),2).dimension()
        3
        sage: ModularSymbols(Gamma0(1),12).dimension()
        3

    If we give an integer N for the congruence subgroup, it defaults to
    `\Gamma_0(N)`::

        sage: ModularSymbols(1,12,-1).dimension()
        1
        sage: ModularSymbols(11,4, sign=1)
        Modular Symbols space of dimension 4 for Gamma_0(11) of weight 4 with sign 1 over Rational Field

    We create some spaces for `\Gamma_1(N)`.

    ::

        sage: ModularSymbols(Gamma1(13),2)
        Modular Symbols space of dimension 15 for Gamma_1(13) of weight 2 with sign 0 and over Rational Field
        sage: ModularSymbols(Gamma1(13),2, sign=1).dimension()
        13
        sage: ModularSymbols(Gamma1(13),2, sign=-1).dimension()
        2
        sage: [ModularSymbols(Gamma1(7),k).dimension() for k in [2,3,4,5]]
        [5, 8, 12, 16]
        sage: ModularSymbols(Gamma1(5),11).dimension()
        20

    We create a space for `\Gamma_H(N)`::

        sage: G = GammaH(15,[4,13])
        sage: M = ModularSymbols(G,2)
        sage: M.decomposition()
        [
        Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Congruence Subgroup Gamma_H(15) with H generated by [4, 13] of weight 2 with sign 0 and over Rational Field,
        Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 5 for Congruence Subgroup Gamma_H(15) with H generated by [4, 13] of weight 2 with sign 0 and over Rational Field
        ]

    We create a space with character::

        sage: e = (DirichletGroup(13).0)^2
        sage: e.order()
        6
        sage: M = ModularSymbols(e, 2); M
        Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
        sage: f = M.T(2).charpoly('x'); f
        x^4 + (-zeta6 - 1)*x^3 - 8*zeta6*x^2 + (10*zeta6 - 5)*x + 21*zeta6 - 21
        sage: f.factor()
        (x - zeta6 - 2) * (x - 2*zeta6 - 1) * (x + zeta6 + 1)^2

    We create a space with character over a larger base ring than the values of the character::

        sage: ModularSymbols(e, 2, base_ring = CyclotomicField(24))
        Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta24^4], sign 0, over Cyclotomic Field of order 24 and degree 8

    More examples of spaces with character::

        sage: e = DirichletGroup(5, RationalField()).gen(); e
        Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -1

        sage: m = ModularSymbols(e, 2); m
        Modular Symbols space of dimension 2 and level 5, weight 2, character [-1], sign 0, over Rational Field

    ::

        sage: m.T(2).charpoly('x')
        x^2 - 1
        sage: m = ModularSymbols(e, 6); m.dimension()
        6
        sage: m.T(2).charpoly('x')
        x^6 - 873*x^4 - 82632*x^2 - 1860496

    We create a space of modular symbols with nontrivial character in
    characteristic 2.

    ::

        sage: G = DirichletGroup(13,GF(4,'a')); G
        Group of Dirichlet characters modulo 13 with values in Finite Field in a of size 2^2
        sage: e = G.list()[2]; e
        Dirichlet character modulo 13 of conductor 13 mapping 2 |--> a + 1
        sage: M = ModularSymbols(e,4); M
        Modular Symbols space of dimension 8 and level 13, weight 4, character [a + 1], sign 0, over Finite Field in a of size 2^2
        sage: M.basis()
        ([X*Y,(1,0)], [X*Y,(1,5)], [X*Y,(1,10)], [X*Y,(1,11)], [X^2,(0,1)], [X^2,(1,10)], [X^2,(1,11)], [X^2,(1,12)])
        sage: M.T(2).matrix()
        [    0     0     0     0     0     0     1     1]
        [    0     0     0     0     0     0     0     0]
        [    0     0     0     0     0 a + 1     1     a]
        [    0     0     0     0     0     1 a + 1     a]
        [    0     0     0     0 a + 1     0     1     1]
        [    0     0     0     0     0     a     1     a]
        [    0     0     0     0     0     0 a + 1     a]
        [    0     0     0     0     0     0     1     0]

    We illustrate the custom_init function, which can be used to make
    arbitrary changes to the modular symbols object before its
    presentation is computed::

        sage: ModularSymbols_clear_cache()
        sage: def custom_init(M):
        ....:     M.customize='hi'
        sage: M = ModularSymbols(1,12, custom_init=custom_init)
        sage: M.customize
        'hi'

    We illustrate the relation between custom_init and use_cache::

        sage: def custom_init(M):
        ....:     M.customize='hi2'
        sage: M = ModularSymbols(1,12, custom_init=custom_init)
        sage: M.customize
        'hi'
        sage: M = ModularSymbols(1,12, custom_init=custom_init, use_cache=False)
        sage: M.customize
        'hi2'

    TESTS:

    We test use_cache::

        sage: ModularSymbols_clear_cache()
        sage: M = ModularSymbols(11,use_cache=False)
        sage: sage.modular.modsym.modsym._cache
        {}
        sage: M = ModularSymbols(11,use_cache=True)
        sage: sage.modular.modsym.modsym._cache
        {(Congruence Subgroup Gamma0(11), 2, 0, Rational Field): <weakref at ...; to 'ModularSymbolsAmbient_wt2_g0_with_category' at ...>}
        sage: M is ModularSymbols(11,use_cache=True)
        True
        sage: M is ModularSymbols(11,use_cache=False)
        False
    """
    from . import ambient
    key = canonical_parameters(group, weight, sign, base_ring)

    if use_cache and key in _cache:
         M = _cache[key]()
         if not (M is None): return M

    (group, weight, sign, base_ring) = key

    M = None
    if arithgroup.is_Gamma0(group):
            if weight == 2:
                M = ambient.ModularSymbolsAmbient_wt2_g0(
                    group.level(),sign, base_ring, custom_init=custom_init)
            else:
                M = ambient.ModularSymbolsAmbient_wtk_g0(
                    group.level(), weight, sign, base_ring, custom_init=custom_init)

    elif arithgroup.is_Gamma1(group):

        M = ambient.ModularSymbolsAmbient_wtk_g1(group.level(),
                            weight, sign, base_ring, custom_init=custom_init)

    elif arithgroup.is_GammaH(group):

        M = ambient.ModularSymbolsAmbient_wtk_gamma_h(group,
                            weight, sign, base_ring, custom_init=custom_init)

    elif isinstance(group, tuple):
        eps = group[0]
        M = ambient.ModularSymbolsAmbient_wtk_eps(eps,
                            weight, sign, base_ring, custom_init=custom_init)

    if M is None:
        raise NotImplementedError("computation of requested space of modular symbols not defined or implemented")

    if use_cache:
        _cache[key] = weakref.ref(M)
    return M

