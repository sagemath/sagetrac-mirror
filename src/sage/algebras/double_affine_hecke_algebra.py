#*****************************************************************************
#       Copyright (C) 2014 Mark Shimozono <mshimo at vt.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
#
# Double affine Hecke algebras
#
################################

import functools

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.abstract_method import abstract_method
from sage.misc.bindable_class import BindableClass
from sage.misc.cachefunc import cached_method
from sage.misc.functional import is_even
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.category import Category
from sage.categories.homset import Hom
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.morphism import SetMorphism
from sage.categories.realizations import Category_realization_of_parent, Realizations
from sage.categories.tensor import tensor
from sage.sets.family import Family, FiniteFamily
from sage.sets.set import Set
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.algebras.multiparameter_hecke_algebra import MultiParameterHeckeAlgebra, ParameterFamilies
from sage.algebras.smash_product_algebra import SmashProductAlgebra
from sage.algebras.affine_hecke_algebra import ExtendedAffineHeckeAlgebra
from sage.modules.free_module_element import vector
from sage.structure.sage_object import SageObject
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement

class DoubleAffineType(SageObject):
    r"""
    Class which specifies a double affine Hecke algebra
    of not-necessarily-reduced possibly-extended affine type.

    INPUTS::

        - ``cartan_type`` -- An irreducible finite Cartan type

    The following are given as optional keyword parameters::

        - ``untwisted`` -- True or False (default: True); whether to use the untwisted or twisted affinization
        - ``reduced`` -- True or False (default: True); whether the affine root system is reduced
        - ``dual_reduced`` -- True or False (default: True); whether the "dual" affine root system is reduced
        - ``general_linear`` -- True or False (default: None, meaning False); if the root system is
        untwisted affine type A, use the weight lattice of the general linear group
        - ``parameters`` -- a Family, dictionary, tuple, list, or None (default: None).
        If a Family is given it should have keys among the strings::

            ('null_root','short','long','zero','doubled','zero_doubled')

        and values in some parent `K` (usually a field). A dictionary can also be given.
        A list or tuple of pairs is also accepted. If ``parameters`` is None then a base ring is created.
        It is the fraction field of a polynomial ring over the rationals whose generators in order are the null
        root parameter `q` followed by Hecke parameters (a sublist of `[v, vl, v0, v2, vz]`) depending on
        which orbits of roots exist for the given double affine type.
        - ``extra_parameters`` -- a Family, dictionary, tuple, list, True, or None (default: None). If a Family,
        dictionary, tuple, or list then the keys are the same as for ``parameters``. The values should be
        independent elements in the same parent `K`. If for an orbit name `o` we have parameter `v` and
        extra parameter `c` then there is a built-in specialization :meth:`restrict_base_ring_map` which
        when invoked, sends `c` to `-1/v+v`. If ``extra_parameters`` is None then the formulas use
        `-1/v+v` rather than the extra parameter `c`, that is, the extra parameters are not used.
        If ``extra_parameters`` is True and ``parameters`` is None, then the above polynomial ring is constructed
        with an extra generator for every root orbit. The list of extra variables is `[c, cl, c0, c2, cz]`.

    EXAMPLES::

        sage: DoubleAffineType("A1")
        Double Affine Type ['A', 1, 1] reduced dual-reduced
        sage: DoubleAffineType(['C',3], untwisted=True, reduced=False, dual_reduced=True)
        Double Affine Type ['C', 3, 1] reduced dual-reduced
        sage: DoubleAffineType(['F',4], untwisted=False, reduced=True, dual_reduced=False)
        Double Affine Type ['F', 4, 1]^* relabelled by {0: 0, 1: 4, 2: 3, 3: 2, 4: 1} reduced dual-reduced
        sage: K = QQ['q,v,vl,v0,v2,vz'].fraction_field()
        sage: q,v,vl,v0,v2,vz = K.gens()
        sage: DoubleAffineType(['B',2], untwisted=False, reduced=False, dual_reduced=False, parameters=[['null_root',q],['short',v],['long',vl],['zero',v0],['doubled',v2],['zero_doubled',vz]])
        Double Affine Type ['C', 2, 1]^* nonreduced nondual-reduced
        sage: DoubleAffineType(['A',2], general_linear=True)
        Double Affine Type GL(3) reduced dual-reduced

    ..RUBRIC:: Untwisted, dual untwisted, and mixed.

    There is a trichotomy of irreducible affine Cartan types: untwisted, dual untwisted (meaning the dual of an
    untwisted affine type) and mixed. The mixed types, `A_{2n}^{(2)}` and its dual, are the ones with three root lengths.

    ..RUBRIC:: Pair of affine Cartan types `\tilde{X}` and `\tilde{Y}`

    The inputs ``cartan_type`` and ``untwisted`` determine a pair of affine Cartan types `\tilde{X}` and `\tilde{Y}`
    which are both untwisted or both dual untwisted. Define

        - `X` -- The irreducible reduced finite Cartan type given by ``cartan_type``

        - `\tilde{X} -- if ``untwisted``, the untwisted affinization of `X`; otherwise, obtained from `X` by taking the dual, then the untwisted affinization, then the affine dual.

        - `W_a(\tilde{X})` -- The nonextended affine Weyl group of type `\tilde{X}`. It is generated by the simple reflections `s_i` for `i` in the affine Dynkin node set `I^X`.

        - `R_{red}(\tilde{X})` -- The set of affine real roots for the reduced affine root system. It is the set of affine root lattice elements of the form `u \alpha_i` where `u \in W_a(\tilde{X})` and `\alpha_i` is a simple root for `i \in I^X`.

        - `Y` -- Equal to `X` if not ``untwisted``; otherwise the dual of `X`

        - `\tilde{Y}` -- Equal to `\tilde{X}` if not ``untwisted``; otherwise the untwisted affinization of `Y`
        
    `W_a(\tilde{Y})` and `R_{red}(\tilde{Y})` are analogously defined.

    Since `X` and `Y` are either equal or dual, there is a natural bijection between their classical
    Dynkin node sets `I^X_0 = I^X \setminus \{0^X\}` and `I^Y_0`. This induces an isomorphism `W(X)\cong W(Y)` of
    the finite Weyl groups.

    ..RUBRIC:: Extended affine Weyl groups

    See :meth:`sage.combinat.root_system.ExtendedAffineWeylGroup` for definitions for extended affine Weyl groups.
    Define

        - `W_e(\tilde{X})` -- The extended affine Weyl group of type `\tilde{X}`.

        - `F^X` -- The fundamental group of `W_e(\tilde{X})` of length zero elements.

    Using the notation `W` for either of the isomorphic groups `W(X)` or `W(Y)` and letting
    `Q^Y` and `P^Y` be the root and weight lattices of type `Y`, there are isomorphisms

    ..MATH::

        W_a(\tilde{X}) \cong W \ltimes Q^Y \\
        W_e(\tilde{X}) \cong W \ltimes P^Y \cong F^X \ltimes W_a(\tilde{X})

    Note: The usual convention in the untwisted case is to use the coroot and coweight lattices of type `X`.

    The group `F^X` acts on `W_a(\tilde{X})` by the group automorphisms

    ..MATH::

        \pi^X s^X_i (\pi^X)^{-1} = s^X_{\pi^X(i)}

    for all `\pi^X \in F^X` and `i \in I^X` where the affine Dynkin automorphism `\pi^X` is regarded as a
    permutation of `I^X`. Similarly we have

    ..MATH::

        W_a(\tilde{Y}) \cong W \ltimes Q^X \\
        W_e(\tilde{Y}) \cong W \ltimes P^X \cong F^Y \ltimes W_a(\tilde{Y})

    ..RUBRIC:: Affine Hecke algebra `H(W_a(\tilde{X}))`

    Let `K` be a commutative ring. For `i \in I^X` let `v_{\alpha_i^X} \in K` be an invertible element.
    They must satisfy `v_{\alpha_i^X} = v_{\alpha_j^X}` if `\alpha_i^X` and `\alpha_j^X` are in the
    same `W_a(\tilde{X})`-orbit. The mnemonic for these parameters is `v^X`.

    The Hecke algebra `H(W_a(\tilde{X});v^X)` is the `K`-algebra with algebra generators
    `T^X_i` for `i \in I^X` with the `T^X_i` satisfying the same braid relations as the
    `s^X_i` do in `W_a(\tilde{X})`, together with the quadratic relations

    ..MATH::

        (T^X_i - v_{\alpha_i^X})(T^X_i + v_{\alpha_i^X}^{-1}) = 0

    for all `i \in I^X`.

    For `w \in W_a(\tilde{X})` let `w=s^X_{i_1}\dotsm s^X_{i_l}` be a reduced expression
    (one with `l` minimal) where `i_1,\dotsc,i_l \in I^X`. There is a well-defined element

    ..MATH::

        T^X_w = T^X_{i_1} \dotsm T^X_{i_l}.

    These form a `K`-basis of `H(W_a(\tilde{X}))`.

    ..RUBRIC:: Extended affine Hecke algebra of `W_e(\tilde{X})`

    Let `K` and `v_{\alpha_i^X}` be as for `H(W_a(\tilde{X}))` except that we impose the additional
    condition that `v_{\alpha_i^X} = v_{\alpha_j^X}` if `\alpha_i^X` and `\alpha_j^X` are in the
    same `W_e(\tilde{X})`-orbit (as opposed to the same `W_a(\tilde{X})`-orbit).
    Recalling that `W_e(\tilde{X}) \cong F^X \ltimes W_a(\tilde{X})`,
    the extended affine Hecke algebra `H(W_e(\tilde{X}),v^X)` is linearly the tensor product

    ..MATH::

        K[F^X] \otimes H(W_a(\tilde{X}),v^X)

    where `K[F^X]` is the group algebra, and `K[F^X]` acts on `H(W_a(\tilde{X}),v^X)` by

    ..MATH::

        \pi^X T^X_i (\pi^X)^{-1} = T^X_{\pi^X(i)}.

    for all `\pi^X \in F^X` and `i \in I^X`.

    Let `v \in W_e(\tilde{X})` be written `v = \pi^X w` with `\pi^X \in F^X` and `w \in W_a(\tilde{X})`.
    `H(W_e(\tilde{X}))` has a basis given by the elements

    ..MATH::

        T^X_v = \pi^X T^X_w.

    ..RUBRIC:: Not-necessarily-reduced affine root systems

    The pair of boolean inputs ``reduced`` and ``dual_reduced`` affects the choice of root versus weight lattices,
    nonextended and extended affine Weyl groups, and nonreduced versus reduced root systems.

    A node `i \in I^X` is *doubleable* if the evaluation `\alpha_i^\vee(\alpha_j)`
    of the `i`-th simple coroot `\alpha_i^\vee` on the `j`-th simple root, is even for all `j \in I^X`.

    Here are some examples of affine cartan types and their sets of doubleable nodes::

        sage: CartanType(['A',1,1]).doubled_nodes()
        (0, 1)
        sage: CartanType(['A',2,1]).doubled_nodes()
        ()
        sage: CartanType(['B',3,1]).doubled_nodes()
        (3,)
        sage: CartanType(['D',4,2]).doubled_nodes()
        (0, 3)
        sage: CartanType(['C',2,1]).doubled_nodes()
        (1,)
        sage: CartanType(['C',3,1]).doubled_nodes()
        ()

    Fact: The only untwisted or dual untwisted affine root systems
    with doubleable nodes, are `A_1^{(1)}` nodes `0,1`; `C_2^{(1)}` node `1`; `B_n^{(1)}` node `n`;
    and `D_{n+1}^{(2)}`, nodes `0,n`.

    ..RUBRIC: Notation

        - `S^X` -- if the input ``reduced`` is False, the set of doubleable nodes in `I^X`;
        otherwise, the empty set. We call this is the set of doubled nodes.
        - `R(\tilde{X})` -- the union of `R_{red}(\tilde{X})` with the elements `u (2\alpha_i)` 
        for `u \in W_a(\tilde{X})` and `i \in S^X`. This is the not-necessarily reduced
        affine root system.

    We say that `R(\tilde{X})` is reduced if it equals `R_{red}(\tilde{X})` (that is, `S^X` is empty)
    and is nonreduced otherwise. Even if ``reduced`` is set to False, if there are no
    doubleable nodes then `R(\tilde{X})` is still reduced.
     
    `S^Y` and `R(\tilde{Y})` are defined analogously, using `\tilde{Y}` instead of `\tilde{X}` and the input boolean
    ``dual_reduced`` instead of ``reduced``.

    ..RUBRIC:: Lattices `X` and `Y` and possibly extended affine Weyl groups `W(\tilde{X})` and `W(\tilde{Y})`

        - `X` -- By abuse of notation this symbol stands for the weight lattice `P^X` if `R(\tilde{Y})` is reduced (no typo with
          the Y here); otherwise, for the root lattice `Q^X`.

        - `\tilde{X}` -- This stands for the direct sum of the lattice `X` with the integer multiples of the affine null root
          `\delta^X`. It is a sublattice of the level-zero part of the affine weight lattice.

        - `W(\tilde{X})` -- If `R(\tilde{Y})` is reduced, this is the extended affine Weyl group `W_e(\tilde{X})`, and otherwise it is the nonextended affine Weyl group `W_a(\tilde{X})`.

        - `Y`, `\tilde{Y}`, `W(\tilde{Y})` -- Defined analogously with the roles of `X` and `Y` interchanged.

    Observe that the lattice `Y` (resp. `X`) provides the translation elements for the "other" affine Weyl group
    `W(\tilde{X})` (resp. `W(\tilde{Y})`):

    ..MATH::

        W(\tilde{X}) \cong W(Y) \ltimes Y
        W(\tilde{Y}) \cong W(X) \ltimes X

    Let `\Pi^X \subset F^X` be the trivial subgroup if `R(\tilde{Y})` is nonreduced and otherwise let `\Pi^X=F^X`.
    Make similar definitions for `\Pi^Y \subset F^Y`. Then

    ..MATH::

        W(\tilde{X}) \cong \Pi^X \ltimes W_a(\tilde{X})
        W(\tilde{Y}) \cong \Pi^Y \ltimes W_a(\tilde{Y})

    ..RUBRIC:: The base ring `K` of the DAHA and parameters

    We assume `K` is a commutative ring and `QQ`-algebra.

    For every `i \in I^X` there is an invertible element `v_{\alpha_i^X}\in K` and
    for every `i \in S^X` an invertible element `v_{2\alpha_i^X} \in K`.
    For convenience we define `v_{2\alpha_i^X}=v_{\alpha_i^X}` for `i \in I^X \setminus S^X`.
    They must satisfy the compatibility that `v_\alpha=v_\beta` if `\alpha` and `\beta` are in the
    same `W(\tilde{X})`-orbit in `R(\tilde{X})`.

    Finally, there is an invertible element `q \in K`. It represents the exponential of the null root
    `\delta^X` and also the exponential of `-\delta^Y`.
    Under this identification there are group homomorphisms
    `\tilde{X} \to K[X]` and `\tilde{Y} \to K[Y]`. We assume `q` is transcendental over `QQ` so that these maps are injective.
    In particular it makes sense to refer to `X^{\alpha_0^X} \in K[X]`.
    
    There is also an action of `W(\tilde{X})` on `K[X]` by `K`-algebra automorphisms:

    ..MATH::

        \begin{align*}
            s^X_i \cdot X^\lambda &= X^{s^X_i(\lambda)} \\
            \pi^X \cdot X^{\alpha_i^X} &= X^{\alpha_{\pi^X(i)}}
        \end{align*}

    where `i \in I^X`, `\lambda \in \tilde{X}`, and `\pi^X \in F^X`. 
    The action of `s^X_0` involves `q` since `\alpha^X_0` does.
    The action of `\Pi^X` has only been defined on the affine root lattice,
    which is only a sublattice of `\tilde{X}`, but this sublattice spans over
    the rationals. The action of `\Pi^X` can also involve `q`.

    ..RUBRIC:: The possibly-extended affine Hecke algebra `H(W(\tilde{X}))`.

    Define the possibly extended affine Hecke algebra `H(W(\tilde{X}))` to be `H(W_e(\tilde{X}))` or
    `H(W_a(\tilde{X}))` according as `W(\tilde{X})` equals `W_e(\tilde{X})` or `W_a(\tilde{X})`.

    ..RUBRIC:: The DAHA, at last!

    Linearly it is the tensor product

    ..MATH::

        H(X,Y;v^X) = H(W(\tilde{X}),v^X) \otimes K[X]

    To define the product it is enough to specify the commutation of the `T^X_i` with `X^\lambda`
    for `i \in I^X` and `\lambda \in X`. This is achieved using the very interesting Demazure-Lusztig operators.

    For nondoubled nodes `i \in I^X \setminus S^X` we have

    ..MATH::

        T^X_i X^\lambda - X^{s_i^X(\lambda)} T^X_i = (v_{\alpha_i^X} - v_{\alpha_i^X}^{-1})
            \dfrac{X^\lambda - X^{s_i^X(\lambda)}}{1 - X^{\alpha_i^X}}

    For doubled nodes `i \in S^X` we have

    ..MATH::

        T^X_i X^\lambda - X^{s_i^X(\lambda)} T^X_i =
            (v_{\alpha_i^X} - v_{\alpha_i^X}^{-1} + (v_{2\alpha_i^X} - v_{2\alpha_i^X}^{-1})X^{\alpha_i^X}
            \dfrac{X^\lambda - X^{s_i^X(\lambda)}}{1 - X^{2 \alpha_i^X}}

    In particular it follows that the DAHA `H(X,Y;v^X)` has bases of the form

    ..MATH::

        T_v X^\lambda \qquad\text{ and }\qquad X^\lambda T_v

    where `v \in W(\tilde{X})` and `\lambda \in X`. We call these bases the "TX_X" and "X_TX" bases respectively.

    ..RUBRIC:: Dual Hecke parameters

    Recall the parameters `v_{\alpha_i^X}` and `v_{2\alpha_i^X}`. 
    Recall that `v_{2\alpha_i^X} = v_{\alpha_i^X}` if `i \in I^X \setminus S^X`.

    The "dual Hecke parameters" are the following relabelling of the Hecke
    parameters::

        - `v_{\alpha_i^Y} = v_{\alpha_i^X}` -- for `i\ne 0`, 
        - `v_{\alpha_0^Y} = v_{2\alpha_r^X}` -- for `r \in I^X \setminus\{0\}` with `\alpha_r^X` short
        - `v_{2\alpha_i^Y} = v_{\alpha_0^X}` -- if `i \in S^Y\setminus \{0\}`
        - `v_{2\alpha_0^Y} = v_{2\alpha_0^X}` -- if `0 \in S^Y`

    where `v_{2\alpha_i^Y} = v_{\alpha_i^Y}` if `i \in I^Y \setminus S^Y`.
    The parameters `v_{\alpha_i^Y}` and `v_{2\alpha_i^Y}` are collectively labeled `v^Y`.

    ..RUBRIC:: Dual presentation of the DAHA

    The DAHA `H(X,Y;v^X)` has a `K`-algebra isomorphism

    ..MATH::

        H(X,Y;v^X) \cong H(W(\tilde{Y});v^Y) \otimes K[Y]

    We first give the commutation between `T^Y_i` and `Y^\lambda` to define the multiplication on
    the right hand side, and then describe the preimages of the right hand side's generators.

    For `i \in I^Y \setminus S^Y` and `\lambda \in Y` we have

    ..MATH::

        T^Y_i Y^\lambda - Y^{s_i^Y(\lambda)} T^Y_i = (v_{\alpha_i^Y} - v_{\alpha_i^Y}^{-1})
            \dfrac{Y^\lambda - X^{s_i^Y(\lambda)}}{1 - Y^{-\alpha_i^Y}}

    which looks like the `X`-analogue except for the `-\alpha_i^Y`.
    For doubled nodes `i \in S^Y` we have

    ..MATH::

        T^Y_i Y^\lambda - Y^{s_i^Y(\lambda)} T^Y_i =
            (v_{\alpha_i^Y} - v_{\alpha_i^Y}^{-1} + (v_{2\alpha_i^Y} - v_{2\alpha_i^Y}^{-1})Y^{-\alpha_i^Y}
            \dfrac{Y^\lambda - Y^{s_i^Y(\lambda)}}{1 - Y^{- 2 \alpha_i^Y}}

    For nonzero `i` we have `T^Y_i = T^X_i`. We have

    ..MATH::

        T^Y_0 = (X^\varphi T_{s_{\varphi}})^{-1}

    where `\varphi \in X` is the dominant short root (its associated coroot is the highest coroot).
    
    For `\pi_i^Y \in \Pi^Y` where `i` is a special node in `I^Y`, if `i=0^Y` then `\pi_i^Y` is the identity.
    For `i \ne 0^Y` let `\omega_i^X` be the `i`-th fundamental weight and `u_i \in W(X)` the shortest element
    such that `u_i(\omega_i)` is antidominant. Then

    ..MATH::

        \pi^Y_i = X^{\omega_i^X} T_{u_i^{-1}}.

    In particular it follows that the DAHA `H(X,Y;v^X)` has bases of the form

    ..MATH::

        T_v Y^\lambda \qquad\text{ and }\qquad Y^\lambda T_v

    where `v \in W(\tilde{Y})` and `\lambda \in Y`. We call these bases the "TY_Y" and "Y_TY" bases respectively.

    ..RUBRIC:: Bases

    Using the isomorphism `K[X] \otimes H(W) \cong H(W(\tilde{Y}))`
    of the type `\tilde{Y}` affine Hecke algebra, we have
    the following presentations of the DAHA.

        - "XT" -- `K[X] \otimes H(W(\tilde{X}))`
        - "XtyY" -- `K[X] \otimes H(W(X)) \otimes K[Y]`
        - "XtxY" -- `K[X] \otimes H(W(Y)) \otimes K[Y]`
        - "TY" -- `H(W(\tilde{Y})) \otimes K[Y]`

    Note that the classical Weyl groups `W(X)` and `W(Y)` are isomorphic.
    The DAHA has a basis corresponding to each of the permutations of the tensor factors
    of these forms. The various bases of the DAHA are implemented by smash product algebras.

    ..RUBRIC:: Implementation of Hecke parameters

    The default base ring is the field `K` given by

    ::

        sage: K = QQ['q,v,vl,v0,v2,vz'].fraction_field()
        sage: q,v,vl,v0,v2,vz = K.gens()

    This can be customized using the ``parameters`` option. Every parameter must be an invertible element of the
    base ring. The default parameters have the following names.

        - `q`  -- 'null_root'. Recall that `q=X^{\delta^X}= Y^{-\delta^Y}`)
        - `v`  -- 'short' (orbit of short `\alpha_i^X` for some nonzero `i`)
        - `vl` -- 'long' (orbit of long `\alpha_i^X` for some nonzero `i` if not simply-laced)
        - `v0` -- 'zero' (orbit of `\alpha_0^X`)
        - `v2` -- 'doubled' (orbit of `2\alpha_i^X` for nonzero doubled node `i`)
        - `vz` -- 'zero_doubled' (orbit of `2\alpha_0^X`)

    In any given case only some of these parameters are used and often they may have the same value.

    EXAMPLES::

        sage: K = QQ['q,v'].fraction_field()
        sage: dat = DoubleAffineType(['A',2], untwisted=True, reduced=True, dual_reduced=True, parameters=[['null_root',K.gen(0)],['short',K.gen(1)]])
        sage: dat.parameter('short')
        v
        sage: K = QQ['q,v,vl,v0,v2,vz'].fraction_field()
        sage: dat = DoubleAffineType(['B',2], untwisted=False, reduced=False, dual_reduced=False, parameters=[['null_root',K.gen(0)],['short',K.gen(1)],['long',K.gen(2)],['zero',K.gen(3)],['doubled',K.gen(4)],['zero_doubled',K.gen(5)]])
        sage: dat.parameter('null_root')
        q
        sage: dat.parameter('doubled')
        v2
        sage: dat.parameter('zero')
        v0
        sage: dat.parameter('zero_doubled')
        vz

    Setting ``extra_parameters`` to True causes extra variables to be added to the base ring::

        sage: d = DoubleAffineType("B2", untwisted=False, reduced=False, dual_reduced=False, extra_parameters=True)
        sage: d.orbit_names()
        ['short', 'long', 'zero', 'doubled', 'zero_doubled']
        sage: d.hecke_parameter()
        Finite family {'zero': v0, 'short': v, 'long': vl, 'doubled': v2, 'zero_doubled': vz}
        sage: d.base_ring()
        Fraction Field of Multivariate Polynomial Ring in q, v, vl, v0, v2, vz, c, cl, c0, c2, cz over Rational Field
        sage: d.c_parameter()
        Finite family {'zero': c0, 'short': c, 'long': cl, 'doubled': c2, 'zero_doubled': cz}
        sage: d.c_parameter('doubled')
        c2
        sage: d.c_parameter(d.hecke_parameter('doubled'))
        c2

    The meaning of the extra parameter `c` corresponding to the Hecke parameter `v` is
    `-v+1/v`, the negative difference reciprocal. This is handy for the Ram-Yip formula
    for Macdonald-Koornwinder polynomials for non-reduced double affine types,
    where the doubled parameters v2 and vz cannot cancel and
    only appear as functions of c2=-v2+1/v2 and cz=-vz+1/vz, in which case we put c2 and cz
    into the formulas instead::

        sage: [[orbit, d.c_parameter(orbit), d.restrict_base_ring(d.c_parameter(orbit))] for orbit in d.hecke_parameter().keys()]
        [['zero', c0, (-v0^2 + 1)/v0], ['short', c, (-v^2 + 1)/v], ['doubled', c2, (-v2^2 + 1)/v2], ['long', cl, (-vl^2 + 1)/vl], ['zero_doubled', cz, (-vz^2 + 1)/vz]]

    REFERENCES:

    .. [Haiman_ICM] M. Haiman, Cherednik algebras, Macdonald polynomials and combinatorics,
       Proceedings of the International Congress of Mathematicians,
       Madrid 2006, Vol. III, 843-872.
    """

    def __init__(self, cartan_type, **keywords):
        from sage.combinat.root_system.cartan_type import CartanType
        self._cartan_type_classical = CartanType(cartan_type)
        if not self._cartan_type_classical.is_finite():
            raise ValueError("Cartan type is not finite")
        if self._cartan_type_classical.rank() == 1:
            self._cartan_type_classical = CartanType(['A',1])

        self._untwisted = keywords.get('untwisted', True)
        self._reduced = keywords.get('reduced', True)
        self._dual_reduced = keywords.get('dual_reduced', True)
        self._general_linear = keywords.get('general_linear', False)

        def check_parameters(parms):
            if isinstance(parms, (list,tuple)):
                parms = dict(parms)
            if isinstance(parms, dict):
                parms = Family(parms)
            if parms is not None:
                for par_name in parms.keys():
                    if par_name not in ('null_root', 'short', 'long', 'zero', 'doubled', 'zero_doubled'):
                        raise ValueError("Parameter name {} is not recognized".format(par_name))
                    try:
                        (parms[par_name])**(-1)
                    except ZeroDivisionError:
                        raise ValueError("Parameter {} is not invertible".format(parms[par_name]))
            return parms

        parameters = keywords.get('parameters', None)
        parameters = check_parameters(parameters)
        extra_parameters = keywords.get('extra_parameters', None)
        if extra_parameters is not True:
            extra_parameters = check_parameters(extra_parameters)
        
        if self._cartan_type_classical.is_simply_laced():
            # in simply-laced type, override `untwisted`
            self._untwisted = True
        if self._untwisted:
            self._cartan_type = self._cartan_type_classical.affine()
            self._other_affine_type = self._cartan_type_classical.dual().affine()
        else:
            self._cartan_type = self._cartan_type_classical.dual().affine().dual()
            self._other_affine_type = self._cartan_type

        I = self._cartan_type.index_set()
        cartan_matrix = self._cartan_type.cartan_matrix()
        # Compute the set of doubled nodes for \tilde{X}
        if self._reduced:
            self._doubled_nodes = tuple([])
        else:
            # uses that the index set starts with 0
            self._doubled_nodes = self._cartan_type.doubled_nodes()
        self._reduced = (len(self._doubled_nodes) == 0)
        self._node_is_doubled = Family(I, lambda x: x in self._doubled_nodes)

        # Compute whether \tilde{Y} is reduced.
        if not self._dual_reduced:
            self._dual_reduced = (len(self._other_affine_type.doubled_nodes()) == 0)

        # general linear root data
        if self.general_linear():
            if self._untwisted and self._cartan_type.type() == 'A' and self._reduced and self._dual_reduced:
                self._n = self._cartan_type.n + 1
            else:
                raise ValueError("General linear double affine data should be untwisted type A, reduced, and dual-reduced")

        # Make the extended affine Weyl group W_e(\tilde{X})
        from sage.combinat.root_system.extended_affine_weyl_group import ExtendedAffineWeylGroup
        self._E = ExtendedAffineWeylGroup(self._cartan_type, general_linear = self.general_linear(), fundamental="")
        # Specify the subgroup of allowable special affine Dynkin automorphisms `\Pi^X`.
        if self._dual_reduced:
            # Set `\Pi^X = F^X`
            self._special_nodes = self._E.cartan_type().special_nodes()
        else:
            # Set `\Pi^X = \{0\}`
            self._special_nodes = tuple([0])

        # establish the specialization of Hecke parameters by computing
        # W(\tilde{X}) orbits of simple roots and their doubles
        vi = dict()
        v2i = dict()

        I0 = self._cartan_type_classical.index_set()
        # set the v_{\alpha_i} for i nonzero
        if self._cartan_type_classical.is_simply_laced():
            # for one root length we call it short for the purposes of the parameters.
            for i in I0:
                vi[i] = 'short'
        else:
            root_lattice = self._cartan_type_classical.root_system().root_lattice()
            for i in I0:
                vi[i] = 'short' if root_lattice.simple_root(i).is_short_root() else 'long'
        # set v_{\alpha_0}
        if self.general_linear():
            vi[0] = vi[1]
        elif len(self._special_nodes) >= 2:
            # \alpha_0 is in the orbit of another simple root by a Dynkin automorphism in W(\tilde{X})
            vi[0] = vi[self._special_nodes[1]]
        else:
            # check the attachment bond of 0 
            for i in I0:
                a0i = cartan_matrix[0,i]
                if a0i != 0:
                    break
            if a0i * cartan_matrix[i,0] == 1:
                # simple bond, and `\alpha_0` and `\alpha_i` are in the same orbit
                vi[0] = vi[i]
            else:
                # double bond, and `\alpha_0` is not in the orbit of any other simple root
                vi[0] = 'zero'

        # make v_{2\alpha_i}
        for i in I:
            v2i[i] = vi[i]

        # there are at most two doubled nodes and if there are two, one must be 0.
        for i in self._doubled_nodes:
            v2i[i] = 'zero_doubled' if i == 0 else 'doubled'

        # for BC with the nontrivial Dynkin automorphism allowed, 2\alpha_0 is in the orbit of 2\alpha_n.
        if 0 in self._doubled_nodes and self._dual_reduced:
            v2i[0] = 'doubled'

        # dictionaries for the orbits of simple (possibly doubled) roots under W(\tilde{X}).
        # the values are the orbit names like 'short', 'long', etc.
        self._vi = Family(I, lambda i: vi[i])
        self._v2i = Family(I, lambda i: v2i[i])

        # remember the distinct orbit names (in a particular order)
        orbit_names = [x for x in Set([self._vi[i] for i in I]+[self._v2i[i] for i in I])]
        # a dictionary from the orbit names to the generic variable names that will be used
        names_variables = [['null_root','q'],['short','v'],['long','vl'],['zero','v0'],['doubled','v2'],['zero_doubled','vz']]
        names_variables_fam = Family(dict(names_variables))
        names = [x[0] for x in names_variables]
        self._orbit_names = []
        for name in names:
            if name in orbit_names:
                self._orbit_names = self._orbit_names + [name]

        self._extended_base_ring = False
        if parameters is None: # if the parameters are not specified make the generic coefficient ring
            def make_var_string(lis):
                names = ''
                for i in range(len(lis)-1):
                    names = names + lis[i] + ','
                names = names + lis[-1]
                return names
            # make a variable for the null root and a Hecke parameter for every orbit
            variable_name_string = make_var_string([names_variables_fam[o] for o in ['null_root']+self._orbit_names])
            if extra_parameters is True: # make an extra c variable for each Hecke parameter
                self._extended_base_ring = True
                c_names = Family(dict([['short','c'],['long','cl'],['zero','c0'],['doubled','c2'],['zero_doubled','cz']]))
                variable_name_string = variable_name_string + ',' + make_var_string([c_names[o] for o in self._orbit_names])
            K = QQ[variable_name_string].fraction_field()
            self._hecke_parameters = Family(dict([[self._orbit_names[i], K.gen(1+i)] for i in range(len(self._orbit_names))]))
            self._parameters = Family(dict([['null_root',K.gen(0)]]+[[name,self._hecke_parameters[name]] for name in self._orbit_names]))
            if extra_parameters is True:
                self._c_parameters = Family(dict([[self._orbit_names[i],K.gen(1+self._hecke_parameters.cardinality()+i)] for i in range(len(self._orbit_names))]))
        else: # if the parameters are supplied externally (in a Family here)
            if 'null_root' not in parameters.keys():
                raise ValueError("The null root parameter is not specified")
            q = parameters['null_root']
            K = q.parent()
            param_dict = dict({'null_root':q})
            for orbit in self._orbit_names:
                if orbit not in parameters.keys():
                    raise ValueError("There is no parameter supplied for the orbit {}".format(orbit))
                param_dict[orbit] = parameters[orbit]
                if param_dict[orbit] not in K:
                    raise ValueError, "The parameter {} supplied for the orbit {} is not in the parent of {}".format(param_dict[orbit],orbit,q)
            self._parameters = Family(['null_root']+self._orbit_names, lambda x: param_dict[x])
            self._hecke_parameters = Family(self._orbit_names, lambda x: self._parameters[x])
            if isinstance(extra_parameters, FiniteFamily):
                self._extended_base_ring = True
                extra_dict = dict()
                for orbit in self._orbit_names:
                    if orbit not in extra_parameters.keys():
                        raise ValueError("There is no extra parameter supplied for the orbit {}".format(orbit))
                    extra_dict[orbit] = extra_parameters[orbit]
                    if extra_dict[orbit] not in K:
                        raise ValueError, "The extra parameter {} supplied for the orbit {} is not in the parent of {}".format(extra_dict[orbit], orbit, q)
                self._c_parameters = Family(self._orbit_names, lambda x: extra_dict[x])
        self._base_ring = K

        # Families for first and second eigenvalues of the T generators
        self._q1 = Family(dict([[i, self._parameters[self._vi[i]]] for i in I]))
        self._q2 = Family(dict([[i, -1/self._q1[i]] for i in I]))
        if self.general_linear():
            self._m = self._n
        else:
            self._m = len(self._special_nodes)

    def _repr_(self):
        r""" A string representing ``self``.
        """
        def non_string(string, bool):
            if bool:
                return string
            return "non"+string
        return "Double Affine Type {} {} {}{}".format(self.cartan_type() if not self.general_linear() else "GL({})".format(self._n), non_string("reduced", self.reduced()), non_string("dual-reduced", self.dual_reduced()),"" if not self.extended_base_ring() else " with extended base ring")

    def untwisted(self):
        return self._untwisted

    def reduced(self):
        r"""
        Is ``self`` reduced?

        EXAMPLES::

            sage: DoubleAffineType(['B',2], untwisted=False, reduced=True, dual_reduced=False).reduced()
            True
            sage: DoubleAffineType(['B',2], untwisted=False, reduced=False, dual_reduced=False).reduced()
            False
            sage: DoubleAffineType(['A',2], untwisted=False, reduced=False, dual_reduced=False).reduced()
            True
        """
        return self._reduced

    def dual_reduced(self):
        r"""
        Is the DAHA dual of ``self`` reduced?

        EXAMPLES::

            sage: DoubleAffineType(['B',2], untwisted=False, reduced=False, dual_reduced=True).dual_reduced()
            True
            sage: DoubleAffineType(['B',2], untwisted=False, reduced=False, dual_reduced=False).dual_reduced()
            False
            sage: DoubleAffineType(['A',2], untwisted=False, reduced=False, dual_reduced=False).dual_reduced()
            True
        """
        return self._dual_reduced

    def properly_extended(self):
        r"""
        Does the group `W(\tilde{X})` have any nontrivial affine Dynkin automorphisms?

        EXAMPLES::

            sage: DoubleAffineType(['A',2]).properly_extended()
            True
            sage: DoubleAffineType(['A',1], dual_reduced=False).properly_extended()
            False
            sage: DoubleAffineType(['B',2], untwisted=False, reduced=False, dual_reduced=True).properly_extended()
            True
            sage: DoubleAffineType(['B',2], untwisted=False, reduced=False, dual_reduced=False).properly_extended()
            False
        """
        return len(self.special_nodes()) > 1

    def general_linear(self):
        r"""
        Is this the general linear double affine type?

            sage: DoubleAffineType(['A',2]).general_linear()
            False
            sage: DoubleAffineType(['A',2],general_linear=True).general_linear()
            True
        """
        return self._general_linear

    def extended_base_ring(self):
        r"""
        Is the larger base ring being used?

            sage: DoubleAffineType(['B',2],reduced=False,extra_parameters=True).extended_base_ring()
            True
            sage: DoubleAffineType(['B',2],reduced=False).extended_base_ring()
            False
        """
        return self._extended_base_ring

    def base_ring(self):
        return self._base_ring

    def doubled_nodes(self):
        r"""
        The set of doubled nodes of type `X`.

        EXAMPLES::

            sage: DoubleAffineType(['B',2], untwisted=False, reduced=False, dual_reduced=False).doubled_nodes()
            (0, 2)
        """
        return self._doubled_nodes

    def cartan_type(self):
        r"""
        The affine Cartan type of ``self``.

        EXAMPLES::

            sage: DoubleAffineType(['B',2], untwisted=False, reduced=True, dual_reduced=False).cartan_type()
            ['C', 2, 1]^*

        """
        return self._cartan_type

    def cartan_type_classical(self):
        r"""
        The finite Cartan type of ``self``.

        EXAMPLES::

            sage: DoubleAffineType(['B',2], untwisted=False, reduced=True, dual_reduced=False).cartan_type_classical()
            ['B', 2]

        """
        return self._cartan_type_classical

    def special_nodes(self):
        r"""
        Return the set of special nodes indexing length-zero elements of `W(\tilde{X})`.

        EXAMPLES::

            sage: DoubleAffineType(['A',3], untwisted=True, reduced=True, dual_reduced=True).special_nodes()
            (0, 1, 2, 3)
            sage: DoubleAffineType(['B',3], untwisted=True, reduced=True, dual_reduced=True).special_nodes()
            (0, 1)
            sage: DoubleAffineType(['B',3], untwisted=False, reduced=False, dual_reduced=False).special_nodes()
            (0,)
            sage: DoubleAffineType(['B',3], untwisted=False, reduced=False, dual_reduced=True).special_nodes()
            (0, 3)

        """
        return self._special_nodes

    def properly_extended(self):
        r"""
        Is the group `W(\tilde{X})` strictly larger than `W_a(\tilde{X})`?

        EXAMPLES::

            sage: DoubleAffineType(['A',3], untwisted=True, reduced=True, dual_reduced=True).properly_extended()
            True
            sage: DoubleAffineType(['F',4], untwisted=True, reduced=True, dual_reduced=True).properly_extended()
            False
            sage: DoubleAffineType(['B',3], untwisted=False, reduced=True, dual_reduced=False).properly_extended()
            False
            sage: DoubleAffineType(['A',3], general_linear=True).properly_extended()
            True
        """

        return self._general_linear or len(self._special_nodes) >= 2

    def orbit_names(self):
        r"""
        A list of strings representing the orbits of roots in ``self``.

        EXAMPLES::

            sage: DoubleAffineType("A2").orbit_names()
            ['short']
            sage: DoubleAffineType("B2").orbit_names()
            ['short', 'long']
            sage: DoubleAffineType("B2",untwisted=False,reduced=False,dual_reduced=False).orbit_names()
            ['short', 'long', 'zero', 'doubled', 'zero_doubled']
        """
        return self._orbit_names

    def parameter(self, param=None):
        r"""
        A parameter of ``self``.

        INPUTS::

        - ``param`` -- A string (default: None) among 'null_root', 'short', 'long', 'zero', 'doubled', 'zero_doubled'

        If ``param`` is None then the dictionary of all parameters is returned.

        EXAMPLES::

            sage: dat = DoubleAffineType("B2",untwisted=False,reduced=False,dual_reduced=False)
            sage: dat.parameter('null_root')
            q
            sage: dat.parameter('short')
            v
            sage: dat.parameter('long')
            vl
            sage: dat.parameter('doubled')
            v2
            sage: dat.parameter('zero_doubled')
            vz
            sage: dat.parameter()
            Finite family {'zero_doubled': vz, 'short': v, 'long': vl, 'zero': v0, 'doubled': v2, 'null_root': q}
        """
        if param is None:
            return self._parameters
        if param not in self._parameters.keys():
            raise ValueError("{} is not a valid parameter".format(param))
        return self._parameters[param]

    def hecke_parameter(self, param=None):
        r"""
        A Hecke parameter.

        INPUTS::

        - ``param`` -- A string (default: None) one of 'short', 'long', 'zero', 'doubled', and 'zero_doubled'.

        If ``param`` is None then the family of all Hecke parameters is returned.
        The values are base ring variables.

        EXAMPLES::

            sage: d = DoubleAffineType("B2", reduced=False)
            sage: d.hecke_parameter()
            Finite family {'short': v, 'long': vl, 'doubled': v2}
            sage: d.hecke_parameter('doubled')
            v2
        """
        if param is None:
            return self._hecke_parameters
        if not param in self.orbit_names():
            raise ValueError("Invalid Hecke parameter")
        return self._hecke_parameters[param]

    def c_parameter(self, param=None):
        r"""
        The extra variable corresponding to an orbit of roots.

        INPUTS::

        - ``param`` -- A string (default: None) one of 'short', 'long', 'zero', 'doubled', and 'zero_doubled',
        or a variable giving a Hecke parameter.

        If ``param`` is None then the family of all Hecke parameters is returned.
        The values are base ring variables.

        EXAMPLES::

            sage: d = DoubleAffineType("B2", reduced=False, extra_parameters=True)
            sage: d.c_parameter()
            Finite family {'short': c, 'long': cl, 'doubled': c2}
            sage: d.c_parameter('doubled')
            c2
        """
        if not self.extended_base_ring():
            raise ValueError("Only valid with extended base ring")
        hecke_parameters = self.hecke_parameter()
        if param is None:
            return self._c_parameters
        if isinstance(param, str):
            if param not in self.orbit_names():
                raise ValueError("Invalid string for Hecke parameter")
            return self._c_parameters[param]
        elif param in self.base_ring():
            for orbit in self.orbit_names():
                if hecke_parameters[orbit] == param:
                    return self._c_parameters[orbit]
            raise ValueError("Ring element is not a Hecke parameter")
        raise TypeError("Not a valid parameter")

    @cached_method
    def q(self):
        r"""
        The null root parameter.

        EXAMPLES::

            sage: DoubleAffineType("B2").q()
            q
        """
        return self.parameter('null_root')

    @cached_method
    def v(self, i=None):
        r"""
        Hecke parameter for a nondoubled simple root.

        INPUTS::

        - `i` -- A Dynkin node index (default: None)

        If `i` is None the family of all Hecke parameters is returned.

        EXAMPLES::

            sage: DoubleAffineType("B2").v()
            Finite family {0: vl, 1: vl, 2: v}
            sage: DoubleAffineType("B2",dual_reduced=False).v(0)
            v0
        """
        if i is None:
            return Family(dict([[i, self.parameter()[self._vi[i]]] for i in self.cartan_type().index_set()]))
        return self.parameter(self._vi[i])

    @cached_method
    def v2(self, i=None):
        r"""
        Hecke parameter for a doubled simple root.

        INPUTS::

        - `i` -- A Dynkin node index (default: None)

        If `i` is None the family of all doubled Hecke parameters is returned.

        EXAMPLES::

            sage: DoubleAffineType("B2",untwisted=False,reduced=False,dual_reduced=False).v2()
            Finite family {0: vz, 1: vl, 2: v2}
            sage: DoubleAffineType("B2",untwisted=False,dual_reduced=False).v2(0)
            v0
        """
        if i is None:
            return Family(dict([[i, self.parameter(self._v2i[i])] for i in self.cartan_type().index_set()]))
        return self.parameter(self._v2i[i])

    def q1(self, i=None):
        r"""
        The first eigenvalue of the Hecke generator `T_i`.

        INPUTS::

        - `i` -- A Dynkin node index (default: None).

        If `i` is None then the output is the family of all first eigenvalues.

        EXAMPLES::

            sage: d = DoubleAffineType("B2", reduced=False)
            sage: [[i, d.q1(i)] for i in d.cartan_type().index_set()]
            [[0, vl], [1, vl], [2, v]]
            sage: DoubleAffineType("B2").q1()
            Finite family {0: vl, 1: vl, 2: v}
        """
        if i is None:
            return self._q1
        return self._q1[i]

    def q2(self, i=None):
        r"""
        The second eigenvalue of the Hecke generator `T_i`.

        INPUTS::

        - `i` -- A Dynkin node index (default: None).

        If `i` is None then the output is the family of all second eigenvalues.

        EXAMPLES::

            sage: d = DoubleAffineType("B2", reduced=False)
            sage: [[i, d.q2(i)] for i in d.cartan_type().index_set()]
            [[0, (-1)/vl], [1, (-1)/vl], [2, (-1)/v]]
            sage: d = DoubleAffineType("B2")
            sage: d.q2()
            Finite family {0: (-1)/vl, 1: (-1)/vl, 2: (-1)/v}
        """
        if i is None:
            return self._q2
        return self._q2[i]

    def extended_affine_weyl(self):
        r"""
        The extended affine Weyl group of ``self``.

        EXAMPLES::

            sage: DoubleAffineType("B2", reduced=False).extended_affine_weyl()
            Extended affine Weyl group of type ['B', 2, 1]
        """
        return self._E

    @cached_method
    def doubled_parameters(self):
        r"""
        A family on the set of doubled nodes whose values are the extra factor in the nonreduced Demazure-Lusztig operators.

        EXAMPLES::

            sage: DoubleAffineType("B2",untwisted=False,reduced=False,dual_reduced=False).doubled_parameters()
            Finite family {0: (vz^2 - 1)/vz, 2: (v2^2 - 1)/v2}
            sage: DoubleAffineType("B2",untwisted=False,reduced=False,dual_reduced=False,extra_parameters=True).doubled_parameters()
            Finite family {0: -cz, 2: -c2}
        """
        if self.extended_base_ring():
            return Family(dict([[i, -self.c_parameter(self.parameter(self._v2i[i]))] for i in self.doubled_nodes()]))
        return Family(dict([[i, self.diff_reciprocal(self.parameter(self._v2i[i]))] for i in self.doubled_nodes()]))

    def diff_reciprocal(self, x):
        return x - 1/x

    @cached_method
    def eta_base_ring_map(self):
        r"""
        Cherednik's eta automorphism of the base field.

        It inverts all of the parameters, sending `q` and each `v_i` to its reciprocal.

        EXAMPLES::

            sage: d = DoubleAffineType("B2", untwisted=False, reduced=False, dual_reduced=False)
            sage: [[v, d.eta_base_ring_map()(v)] for v in d.parameter().values()]
            [[vz, 1/vz], [v, 1/v], [v2, 1/v2], [v0, 1/v0], [vl, 1/vl], [q, 1/q]]
            sage: d = DoubleAffineType("B2", untwisted=False, reduced=False, dual_reduced=False, extra_parameters=True)
            sage: [[x, d.eta_base_ring_map()(x)] for x in d.base_ring().gens()]
            [[q, 1/q], [v, 1/v], [vl, 1/vl], [v0, 1/v0], [v2, 1/v2], [vz, 1/vz], [c, -c], [cl, -cl], [c0, -c0], [c2, -c2], [cz, -cz]]
        """
        K = self.base_ring()
        images = [1/self.q()] + [1/self.parameter(o) for o in self.orbit_names()]
        if self.extended_base_ring():
            images = images + [-self.c_parameter(o) for o in self.orbit_names()]
        return K.hom(images, codomain=K)

    def dual(self):
        r"""
        The DAHA dual double affine type.

        EXAMPLES::

            sage: dat = DoubleAffineType(['B',2], untwisted=False, reduced=False); dat
            Double Affine Type ['C', 2, 1]^* nonreduced dual-reduced
            sage: dat.dual()
            Double Affine Type ['C', 2, 1]^* reduced nondual-reduced
            sage: dat = DoubleAffineType(['B',3], untwisted=True, reduced=False, dual_reduced=True); dat
            Double Affine Type ['B', 3, 1] nonreduced dual-reduced
            sage: dat.cartan_type().dynkin_diagram()
                O 0
                |
                |
            O---O=>=O
            1   2   3
            B3~
            sage: ddat = dat.dual(); ddat
            Double Affine Type ['C', 3, 1] reduced nondual-reduced
            sage: ddat.cartan_type().dynkin_diagram()
            O=>=O---O=<=O
            0   1   2   3
            C3~
            sage: dat = DoubleAffineType(['B',2], untwisted=True, reduced=False, dual_reduced=True); dat
            Double Affine Type ['B', 2, 1] nonreduced dual-reduced
            sage: dat.cartan_type().dynkin_diagram()
            O=>=O=<=O
            0   2   1
            B2~
            sage: ddat = dat.dual(); ddat
            Double Affine Type ['C', 2, 1] reduced nondual-reduced            
            sage: ddat.cartan_type().dynkin_diagram()
            O=>=O=<=O
            0   1   2
            C2~
            sage: dat = DoubleAffineType(['B',2],untwisted=False,reduced=False,dual_reduced=False,extra_parameters=True); dat
            Double Affine Type ['C', 2, 1]^* nonreduced nondual-reduced with extended base ring
            sage: ddat = dat.dual(); ddat
            Double Affine Type ['C', 2, 1]^* nonreduced nondual-reduced with extended base ring
        """
        if self.untwisted():
            dual_classical_type = self.cartan_type_classical().dual()
        else:
            dual_classical_type = self.cartan_type_classical()
        dual_parameters = dict()
        for x in self.parameter().keys():
            if x == 'doubled':
                dual_parameters['zero'] = self.parameter('doubled')
            elif x == 'zero':
                dual_parameters['doubled'] = self.parameter('zero')
            elif x in ('zero_doubled','null_root'):
                dual_parameters[x] = self.parameter(x)
            elif x == 'short':
                if not self.untwisted() or dual_classical_type.is_simply_laced():
                    dual_parameters[x] = self.parameter(x)
                else:
                    dual_parameters['long'] = self.parameter('short')
            elif x == 'long':
                if not self.untwisted(): # existence of long orbit means not simply laced
                    dual_parameters[x] = self.parameter(x)
                else:
                    dual_parameters['short'] = self.parameter('long')
            else:
                raise ValueError("Illegal parameter")
        dual_parameters = Family(dual_parameters)
        if not self.extended_base_ring():
            extra_dual_parameters = None
        else:
            extra_dual_parameters = dict()
            for x in self._orbit_names:
                if x == 'doubled':
                    extra_dual_parameters['zero'] = self.c_parameter('doubled')
                elif x == 'zero':
                    extra_dual_parameters['doubled'] = self.c_parameter('zero')
                elif x == 'zero_doubled':
                    extra_dual_parameters[x] = self.c_parameter(x)
                elif x == 'short':
                    if not self.untwisted() or dual_classical_type.is_simply_laced():
                        extra_dual_parameters[x] = self.c_parameter(x)
                    else:
                        extra_dual_parameters['long'] = self.c_parameter('short')
                elif x == 'long':
                    if not self.untwisted() or dual_classical_type.is_simply_laced():
                        extra_dual_parameters[x] = self.c_parameter(x)
                    else:
                        extra_dual_parameters['short'] = self.c_parameter('long')
                else:
                    raise ValueError("Illegal parameter")
            extra_dual_parameters = Family(extra_dual_parameters)

        return DoubleAffineType(dual_classical_type, untwisted=self.untwisted(), reduced=self.dual_reduced(), dual_reduced=self.reduced(), general_linear=self.general_linear(), parameters=dual_parameters, extra_parameters=extra_dual_parameters)

    @cached_method
    def restrict_base_ring_map(self):
        r"""
        The endomorphism of the extended base ring that sends each `c` parameter to `-v + 1/v`
        where `v` is the corresponding Hecke parameter.

        Only valid if the base ring is extended.
        """
        if not self.extended_base_ring():
            raise ValueError("Only valid if the base ring is extended")
        K = self.base_ring()
        return K.hom([self.q()] + [self.hecke_parameter(orbit) for orbit in self.orbit_names()] + 
[-self.diff_reciprocal(self.hecke_parameter(orbit)) for orbit in self.orbit_names()], codomain=K)

    def restrict_base_ring(self, f):
        return self.restrict_base_ring_map()(f)

    @cached_method
    def partial_restrict_base_ring_map(self):
        r"""
        The endomorphism of the extended base ring that sends certain `c` parameters to `-v + 1/v`
        where `v` is the corresponding Hecke parameter.

        This is only done for the 'short', 'long', and 'zero' orbits.
        The c parameters for the 'doubled' and 'zero_doubled' orbits are left alone.

        Only valid if the base ring is extended.
        """
        if not self.extended_base_ring():
            raise ValueError("Only valid if the base ring is extended")
        K = self.base_ring()
        return K.hom([self.q()] + [self.hecke_parameter(orbit) for orbit in self.orbit_names()] + 
[-self.diff_reciprocal(self.hecke_parameter(orbit)) for orbit in self.orbit_names() if orbit in ('short','long','zero')] + [self.c_parameter(orbit) for orbit in self.orbit_names() if orbit in ('doubled','zero_doubled')], codomain=K)

    def partial_restrict_base_ring(self, f):
        return self.partial_restrict_base_ring_map()(f)

class DoubleAffineHeckeAlgebraSansDuality(UniqueRepresentation, Parent):
    r"""
    The double affine Hecke algebra, with all its structure except duality.

    The bases that are supported by this class are:

    - "LT" -- `X^\mu \pi^X T_w` where `\pi^X` is in the fundamental group `F^X`, `w` is in the affine Weyl group `W_a(\tilde{X})`,
                   and `\mu` is in the lattice `X`
    - "TL" -- `\pi^X T_w X^\mu` 
    - "LtvLv" -- `X^\mu T_w Y^\nu` where `\mu \in X`, `w` is in the finite Weyl group `W(Y)`, and `\nu \in Y`.
    - "LLvtv" -- `X^\mu Y^\nu T_w`

    EXAMPLES::

        sage: HH = DoubleAffineHeckeAlgebraSansDuality(['A',2])
        sage: LT = HH.LT(); LT
        LT basis of The double affine Hecke algebra of type ['A', 2, 1]
        sage: m = LT.a_monomial(); m
        X[(2, 2, 3)] piX[2] TX[0,1,2]
        sage: TL = HH.TL(); TL
        TL basis of The double affine Hecke algebra of type ['A', 2, 1]
        sage: TL(m)
        ((q*v^2-q)/v)*piX[2] TX[0,1] X[(2, 2, 3)] + ((v^6-3*v^4+3*v^2-1)/(q^2*v^3))*piX[2] X[(2, 2, 3)] + ((-v^4+2*v^2-1)/(q^2*v^2))*piX[2] TX[2] X[(2, 2, 3)] + ((-v^4+2*v^2-1)/(q^2*v^2))*piX[2] TX[1] X[(2, 2, 3)] + ((v^2-1)/(q^2*v))*piX[2] TX[1,2] X[(2, 2, 3)] + q*piX[2] TX[0,1,2] X[(2, 3, 2)]
        sage: LtvLv = HH.LtvLv(); LtvLv
        LtvLv basis of The double affine Hecke algebra of type ['A', 2, 1]
        sage: LtvLv(m)
        X[(2, 2, 3)] Ty[2] Y[(0, 0, -1)]
        sage: LLvtv = HH.LLvtv(); LLvtv
        LLvtv basis of The double affine Hecke algebra of type ['A', 2, 1]
        sage: LLvtv(m)
        ((v^2-1)/v)*X[(2, 2, 3)] Y[(1, 1, 0)] + X[(2, 2, 3)] Y[(1, 0, 1)] Ty[2]        

    The first two have "native" product; the other two have product via coercion.
    This means that serious computations should be done in the "LT" or "TL" bases
    and only at the end should they be coerced into the other bases.

    There is an input option ``dual_side`` which, if True, essentially applies a form of Macdonald duality:
    The interpretation of the bases become:

    - "LT" -- `Y^\mu \pi^Y T_w` where `\pi^Y` is in the fundamental group `F^Y`, `w` is in the affine Weyl group `W_a(\tilde{Y})`,
                   and `\mu` is in the lattice `Y`
    - "TL" -- `\pi^Y T_w Y^\mu` 
    - "LtvLv" -- `Y^\mu T_w X^\nu` where `\mu \in Y`, `w` is in the finite Weyl group `W(X)`, and `\nu \in X`.
    - "LLvtv" -- `Y^\mu X^\nu T_w`

    Other differences: 

    - The Demazure-Lusztig operators use the "dominant" convention unless ``dual_side`` is True, whence they use the
    "antidominant" convention
    - The null root parameter `q` is interpreted as the exponential of the null root `\delta^X` unless ``dual_side`` is True,
    in which case `q` is the exponential of `- \delta^Y`.

    By invoking this class with a double affine type and again with the dual double affine type,
    one obtains many realizations of the same DAHA. However the interaction between these two kinds of realizations
    requires additional machinery, which is included in a separate class.

    ..RUBRIC:: General Linear case

    Here the usual extended affine Hecke algebra of `SL_n` is replaced by that for `GL_n` and the
    group algebra of the weight lattice of `SL_n` is replaced by that of `GL_n`.

    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, untwisted=True, reduced=True, dual_reduced=True, dual_side=False, general_linear=False, parameters=None, extra_parameters=None):
        from sage.combinat.root_system.cartan_type import CartanType
        cartan_type = CartanType(cartan_type)
        def fix_parameters(parameters):
            if isinstance(parameters, list):
                parameters = dict(parameters)
            if isinstance(parameters, dict):
                parameters = Family(parameters)
            return parameters
        parameters = fix_parameters(parameters)
        extra_parameters = fix_parameters(extra_parameters)

        return super(DoubleAffineHeckeAlgebraSansDuality, cls).__classcall__(cls, cartan_type, untwisted, reduced, dual_reduced, dual_side, general_linear, parameters, extra_parameters)

    def __init__(self, cartan_type, untwisted, reduced, dual_reduced, dual_side, general_linear, parameters, extra_parameters):
        self._dat = DoubleAffineType(cartan_type, untwisted=untwisted, reduced=reduced, dual_reduced=dual_reduced, general_linear=general_linear, parameters=parameters, extra_parameters=extra_parameters)
        self._base_ring = self.double_affine_type()._base_ring # this makes .base_ring() work
        self._dual_dat = self.double_affine_type().dual() # set this before any call to `dual_double_affine_type`
        self._m = max(self.double_affine_type()._m, self.dual_double_affine_type()._m)
        assert self.double_affine_type().base_ring() == self.dual_double_affine_type().base_ring(), "The base ring (%s) of the double affine type does not agree with that (%s) of its dual"%(self.double_affine_type().base_ring(), self.dual_double_affine_type().base_ring())

        Parent.__init__(self, category = AlgebrasWithBasis(self.base_ring()).WithRealizations())

        self._dual_side = dual_side
        if self.dual_side():
            prefixL = "Y"
            prefixLv = "X"
            self._convention = "antidominant"
            self._the_q_unscaled = self.double_affine_type().q()**(-1)
        else:
            prefixL = "X"
            prefixLv = "Y"
            self._convention = "dominant"
            self._the_q_unscaled = self.double_affine_type().q()
        self._the_q = self._the_q_unscaled ** self._m

        self._L = self.double_affine_type().cartan_type_classical().root_system().ambient_space()
        self._KL = self._L.algebra(self.base_ring(), prefix=prefixL)
        self._KL_generic = self._L.algebra(self._base_ring)
        self._from_generic_lattice_algebra = self._KL_generic.module_morphism(on_basis = lambda x: self._KL.monomial(x), codomain=self._KL)
        self._W = self._L.weyl_group()
        self._H = ExtendedAffineHeckeAlgebra(self.cartan_type(), self.double_affine_type().q1(), self.double_affine_type().q2(), self.double_affine_type().properly_extended(), self.dual_side(), self.double_affine_type().general_linear())
        self._Lv = self._H.dual_lattice()
        self._KLv = self._Lv.algebra(self.base_ring(), prefix=prefixLv)
        self._Wv = self._H.dual_classical_weyl()

        assert self.base_ring() == self._KL.base_ring()
        if self.base_ring() != self._H.base_ring():
            raise TypeError("parent of parameters:\n {}\n q1: {}\n q2: {}\n base ring of self:\n {}\n base ring of extended affine Hecke:\n {}".format(self.double_affine_type().q1(0).parent(),self.double_affine_type().q1(), self.double_affine_type().q2(), self.base_ring(), self._H.base_ring()))

        # set up eigenvalue computations
        QYvee = self._Lv.cartan_type().root_system().coroot_space()
        parameters = self.double_affine_type().parameter()
        if 'long' in parameters.keys() and parameters['long'] != parameters['short']:
            self._long_orbit_exists = True
            # sum of short coroots of type Y; short here means strictly shorter than some other root
            self._two_rho_short = QYvee.sum([alpha for alpha in QYvee.positive_roots() if alpha.is_short_root()]).to_ambient()
            # sum of long coroots of type Y
            self._two_rho_long = QYvee.sum([alpha for alpha in QYvee.positive_roots() if not alpha.is_short_root()]).to_ambient()
            if not self.double_affine_type().untwisted():
                self._two_rho_short, self._two_rho_long = self._two_rho_long, self._two_rho_short
        else:
            self._long_orbit_exists = False
            self._two_rho_short = QYvee.sum([alpha for alpha in QYvee.positive_roots()]).to_ambient()
            self._two_rho_long = QYvee.zero().to_ambient()
        self._two_rho_zero = QYvee.zero().to_ambient()

        if 'zero' in parameters.keys():
            if parameters['zero'] != parameters['short']:
                cartan_type_classical = self.cartan_type().classical()
                if self._long_orbit_exists and parameters['zero'] != parameters['long'] and not cartan_type_classical.root_system().root_space().simple_root(cartan_type_classical.n).is_short_root():
                    self._two_rho_zero = self._two_rho_long/2
                    self._two_rho_long = self._two_rho_long/2
                else:
                    self._two_rho_zero = self._two_rho_short/2
                    self._two_rho_short = self._two_rho_short/2
            
        # the ingredients of different components of the DAHA
        self._T = self._H.T() # T-basis of extended affine Hecke algebra
        self._E = self._H.extended_affine_weyl()
        self._W0Pv = self._E.W0Pv()
        self._F = self._H.fundamental_group()
        self._tvLv = self._H.tvLv()
        self._Lvtv = self._H.Lvtv()
        self._tv = self._H.dual_classical_hecke()

        cat = ModulesWithBasis(self._base_ring)
        tcat = cat.TensorProducts()
        self._LTmod = tensor([self._KL, self._T], category = tcat)
        self._TLmod= tensor([self._T, self._KL], category = tcat)

        # working here
        LT = self.LT()
        TL = self.TL()
        LT.register_opposite(TL)

        SetMorphism(Hom(LT.factor(0),LT,category=cat),LT.factor_embedding(0)).register_as_coercion()
        SetMorphism(Hom(LT.factor(1),LT,category=cat),LT.factor_embedding(1)).register_as_coercion()
        SetMorphism(Hom(TL.factor(0),TL,category=cat),TL.factor_embedding(0)).register_as_coercion()
        SetMorphism(Hom(TL.factor(1),TL,category=cat),TL.factor_embedding(1)).register_as_coercion()

        self._LtvLvmod = tensor([self._KL, self._tvLv], category = tcat)
        self._tvLvLmod= tensor([self._tvLv, self._KL], category = tcat)

        self._LLvtvmod = tensor([self._KL, self._Lvtv], category = tcat)
        self._LvtvLmod= tensor([self._Lvtv, self._KL], category = tcat)

        LtvLv = self.LtvLv()
        LLvtv = self.LLvtv()

        SetMorphism(Hom(LtvLv.factor(0),LtvLv,category=cat),LtvLv.factor_embedding(0)).register_as_coercion()
        SetMorphism(Hom(LtvLv.factor(1),LtvLv,category=cat),LtvLv.factor_embedding(1)).register_as_coercion()

        SetMorphism(Hom(LLvtv.factor(0),LLvtv,category=cat),LLvtv.factor_embedding(0)).register_as_coercion()
        SetMorphism(Hom(LLvtv.factor(1),LLvtv,category=cat),LLvtv.factor_embedding(1)).register_as_coercion()

        def LtvLv_to_LT_on_basis((mu, w, nu)):
            return LT.from_direct_product((self._KL.monomial(mu),self._T(self._tvLv.monomial((w,nu)))))

        LtvLv_to_LT = LtvLv.module_morphism(on_basis=LtvLv_to_LT_on_basis, codomain=LT)
        LtvLv_to_LT.register_as_coercion()

        def LT_to_LtvLv_on_basis((mu, pi, w)):
            return LtvLv.from_direct_product((self._KL.monomial(mu),self._tvLv(self._T.monomial((pi,w)))))

        LT_to_LtvLv = LT.module_morphism(on_basis=LT_to_LtvLv_on_basis, codomain=LtvLv)
        LT_to_LtvLv.register_as_coercion()

        def LLvtv_to_LT_on_basis((mu, nu, w)):
            return LT.from_direct_product((self._KL.monomial(mu),self._T(self._Lvtv.monomial((nu,w)))))

        LLvtv_to_LT = LLvtv.module_morphism(on_basis=LLvtv_to_LT_on_basis, codomain=LT)
        LLvtv_to_LT.register_as_coercion()

        def LT_to_LLvtv_on_basis((mu, pi, w)):
            return LLvtv.from_direct_product((self._KL.monomial(mu),self._Lvtv(self._T.monomial((pi,w)))))

        LT_to_LLvtv = LT.module_morphism(on_basis=LT_to_LLvtv_on_basis, codomain=LLvtv)
        LT_to_LLvtv.register_as_coercion()

    def double_affine_type(self):
        r"""
        Return the double affine type of ``self``.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B3")
            sage: HH.double_affine_type()
            Double Affine Type ['B', 3, 1] reduced dual-reduced
        """
        return self._dat

    def dual_double_affine_type(self):
        r"""
        Return the dual double affine type of ``self``.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B3")
            sage: HH.dual_double_affine_type()
            Double Affine Type ['C', 3, 1] reduced dual-reduced
        """
        return self._dual_dat

    def dual_side(self):
        return self._dual_side

    def cartan_type(self):
        return self.double_affine_type().cartan_type()

    def base_ring(self):
        return self._base_ring

    @cached_method
    def Y_to_X(self):
        r"""
        The morphism from the ambient space of type `Y` to that of type `X`.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B3")
            sage: dct = HH.cartan_type().classical().dual()
            sage: I0 = dct.index_set()
            sage: [(i, HH._Lv.fundamental_weight(i)) for i in I0]
            [(1, (1, 0, 0)), (2, (1, 1, 0)), (3, (1, 1, 1))]
            sage: [(i, HH.Y_to_X()(HH._Lv.fundamental_weight(i))) for i in I0]
            [(1, (1, 0, 0)), (2, (1, 1, 0)), (3, (1, 1, 1))]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("C3")
            sage: dct = HH.cartan_type().classical().dual()
            sage: I0 = dct.index_set()
            sage: [(i, HH._Lv.fundamental_weight(i)) for i in I0]
            [(1, (1, 0, 0)), (2, (1, 1, 0)), (3, (1/2, 1/2, 1/2))]
            sage: [(i, HH.Y_to_X()(HH._Lv.fundamental_weight(i))) for i in I0]
            [(1, (1, 0, 0)), (2, (1, 1, 0)), (3, (1/2, 1/2, 1/2))]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("F4")
            sage: dct = HH.cartan_type().classical().dual()
            sage: I0 = dct.index_set()
            sage: [(i, HH._Lv.fundamental_weight(i)) for i in I0]
            [(1, (1, 0, 0, 0)), (2, (3/2, 1/2, 1/2, 1/2)), (3, (2, 1, 1, 0)), (4, (1, 1, 0, 0))]
            sage: [(i, HH.Y_to_X()(HH._Lv.fundamental_weight(i))) for i in I0]
            [(1, (1, 1, 0, 0)), (2, (2, 1, 1, 0)), (3, (3, 1, 1, 1)), (4, (2, 0, 0, 0))]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("G2")
            sage: dct = HH.cartan_type().classical().dual()
            sage: I0 = dct.index_set()
            sage: [(i, HH._Lv.fundamental_weight(i)) for i in I0]
            [(1, (2, -1, -1)), (2, (1, 0, -1))]
            sage: [(i, HH.Y_to_X()(HH._Lv.fundamental_weight(i))) for i in I0]
            [(1, (1, 0, -1)), (2, (2/3, -1/3, -1/3))]

        TODO: Move to ambient space code?

        """
        if not self.cartan_type().is_untwisted_affine() or self._dat._general_linear:
            return lambda x: x
        cartan_type = self.cartan_type().classical()
        if cartan_type.is_simply_laced():
            return lambda x: x
        typ = cartan_type.type()
        if typ in ('B','C'):
            return lambda v: self._L.from_vector(v.to_vector())
        if typ == 'F':
            # F4 and G2 have dual implemented as Dynkin reversals of themselves
            def YXF(mu):
                return self._L.from_vector(vector((mu[0]+mu[1],mu[0]-mu[1],mu[2]+mu[3],mu[2]-mu[3]),QQ))
            return YXF
        if typ == 'G':
            def YXG(mu):
                return self._L.from_vector(vector(((2*mu[0]+mu[1])/3,(mu[0]+2*mu[2])/3,(2*mu[1]+mu[2])/3),QQ))
            return YXG
        raise TypeError, "%s should be an irreducible finite Cartan type"%self.cartan_type()

    @cached_method
    def Y_pair_X(self):
        r"""
        Return the pairing between `Y` and `X`.

        This computes the pairing in [Haiman_ICM]_.

        If `Y` equals `X` (dual untwisted case), the pairing is the unique Weyl-invariant one
        with short roots of square length two. This is given by the ambient space pairing except in
        types `B_n` and `F_4`, in which case the ambient space pairing must be doubled.
        If `Y` is the dual of `X` (untwisted case) then the weight lattice of `Y` is identified with the coweight
        lattice of `X` which is then paired with the weight lattice of `X`.

        EXAMPLES::

            sage: E = DoubleAffineHeckeAlgebraSansDuality("A1")
            sage: I0 = E.cartan_type().classical().index_set()
            sage: Matrix([[E.Y_pair_X()(E._Lv.simple_root(i), E._L.simple_root(j)) for j in I0] for i in I0])
            [2]
            sage: E.cartan_type().classical().cartan_matrix()
            [2]
            sage: E = DoubleAffineHeckeAlgebraSansDuality("A2")
            sage: I0 = E.cartan_type().classical().index_set()
            sage: Matrix([[E.Y_pair_X()(E._Lv.simple_root(i), E._L.simple_root(j)) for j in I0] for i in I0])
            [ 2 -1]
            [-1  2]
            sage: E.cartan_type().classical().cartan_matrix()
            [ 2 -1]
            [-1  2]
            sage: E = DoubleAffineHeckeAlgebraSansDuality("B2")
            sage: I0 = E.cartan_type().classical().index_set()
            sage: Matrix([[E.Y_pair_X()(E._Lv.simple_root(i), E._L.simple_root(j)) for j in I0] for i in I0])
            [ 2 -1]
            [-2  2]
            sage: E.cartan_type().classical().cartan_matrix()
            [ 2 -1]
            [-2  2]

            sage: E = DoubleAffineHeckeAlgebraSansDuality("C2")
            sage: I0 = E.cartan_type().classical().index_set()
            sage: Matrix([[E.Y_pair_X()(E._Lv.simple_root(i), E._L.simple_root(j)) for j in I0] for i in I0])
            [ 2 -2]
            [-1  2]
            sage: E.cartan_type().classical().cartan_matrix()
            [ 2 -2]
            [-1  2]

            sage: E = DoubleAffineHeckeAlgebraSansDuality("F4")
            sage: I0 = E.cartan_type().classical().index_set()
            sage: Matrix([[E.Y_pair_X()(E._Lv.simple_root(i), E._L.simple_root(j)) for j in I0] for i in I0])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -2  2 -1]
            [ 0  0 -1  2]
            sage: E.cartan_type().classical().cartan_matrix()
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -2  2 -1]
            [ 0  0 -1  2]

            sage: E = DoubleAffineHeckeAlgebraSansDuality("G2")
            sage: I0 = E.cartan_type().classical().index_set()
            sage: Matrix([[E.Y_pair_X()(E._Lv.simple_root(i), E._L.simple_root(j)) for j in I0] for i in I0])
            [ 2 -3]
            [-1  2]            
            sage: E.cartan_type().classical().cartan_matrix()
            [ 2 -3]
            [-1  2]

            sage: E = DoubleAffineHeckeAlgebraSansDuality("B2", untwisted=False)
            sage: I0 = E.cartan_type().classical().index_set()
            sage: Matrix([[E.Y_pair_X()(E._Lv.simple_root(i), E._L.simple_root(j)) for j in I0] for i in I0])
            [ 4 -2]
            [-2  2]

            sage: E = DoubleAffineHeckeAlgebraSansDuality("C2", untwisted=False)
            sage: I0 = E.cartan_type().classical().index_set()
            sage: Matrix([[E.Y_pair_X()(E._Lv.simple_root(i), E._L.simple_root(j)) for j in I0] for i in I0])
            [ 2 -2]
            [-2  4]

            sage: E = DoubleAffineHeckeAlgebraSansDuality("F4", untwisted=False)
            sage: I0 = E.cartan_type().classical().index_set()
            sage: Matrix([[E.Y_pair_X()(E._Lv.simple_root(i), E._L.simple_root(j)) for j in I0] for i in I0])
            [ 4 -2  0  0]
            [-2  4 -2  0]
            [ 0 -2  2 -1]
            [ 0  0 -1  2]

            sage: E = DoubleAffineHeckeAlgebraSansDuality("G2", untwisted=False)
            sage: I0 = E.cartan_type().classical().index_set()
            sage: Matrix([[E.Y_pair_X()(E._Lv.simple_root(i), E._L.simple_root(j)) for j in I0] for i in I0])
            [ 2 -3]
            [-3  6]

        """
        cartan_type = self.cartan_type().classical()
        typ = cartan_type.type()
        if typ == 'A' and not self._dat._general_linear:
            h = cartan_type.n + 1
            # this normalization is needed for the type A ambient space
            def type_a_pairing(y, x):
                return y.scalar(x) - (QQ.sum(vector(x))*QQ.sum(vector(y)))/h
            return type_a_pairing
        if cartan_type.is_simply_laced():
            return lambda y,x: y.scalar(x)
        if self.cartan_type().is_untwisted_affine():
            return lambda y, x: x.scalar(self.Y_to_X()(y))
        if typ in ('B','F'):
            return lambda y, x: 2 * y.scalar(x)
        return lambda y, x: y.scalar(x)

    def Y_pair_X_m(self, y, x, a=None):
        r"""
        This is a hack to deal with fractional powers of q.

        INPUTS::

            - `y` -- Element of `Y` lattice
            - `x` -- Element of `X` lattice
            - `a` -- optional (default: None, meaning zero) integer giving the coefficient of the null root

        """
        if a is None:
            a = 0
        return (self.Y_pair_X()(y, x)+a) * self._m

    def lattice(self):
        r"""
        The lattice.
        """
        return self._L

    def dual_lattice(self):
        r"""
        The dual lattice.
        """
        return self._Lv

    def Y_eigenvalue(self, mu, la=None,a=None):
        r"""
        The function which, given the pair `(\mu,\lambda)`, returns the eigenvalue of
       `Y^{\mu+a\delta}` on the nonsymmetric Macdonald polynomial `E_\lambda`.

        EXAMPLES::

            sage: HH=DoubleAffineHeckeAlgebraSansDuality("A2")
            sage: mu = HH.dual_lattice().an_element(); mu
            (2, 2, 3)
            sage: la = HH.lattice().an_element(); la
            (2, 2, 3)
            sage: HH.Y_eigenvalue(mu,la)
            1/(q^2*v^2)

            sage: K = QQ['q,v'].fraction_field()
            sage: HH=DoubleAffineHeckeAlgebraSansDuality("B2", parameters = [['null_root',K.gen(0)],['short',K.gen(1)],['long',K.gen(1)]])
            sage: mu = HH.dual_lattice().an_element(); mu
            (2, 2)
            sage: la = HH.lattice().an_element(); la
            (2, 2)
            sage: HH.Y_eigenvalue(mu,la)
            1/(q^16*v^8)            

            sage: HH=DoubleAffineHeckeAlgebraSansDuality("B2", dual_reduced=False, parameters = [['null_root',K.gen(0)],['short',K.gen(1)],['long',K.gen(1)],['zero',K.gen(1)]])
            sage: mu = HH.dual_lattice().an_element(); mu
            (2, 2)
            sage: la = HH.lattice().an_element(); la
            (2, 2)
            sage: HH.Y_eigenvalue(mu,la)
            1/(q^16*v^8)            

            sage: HH=DoubleAffineHeckeAlgebraSansDuality("B3", untwisted=False,reduced=False,dual_reduced=False)
            sage: mu = HH.dual_lattice().from_vector(vector((-3,-2,0))); mu
            (-3, -2, 0)
            sage: la = HH.lattice().from_vector(vector((2,0,-1))); la
            (2, 0, -1)
            sage: HH.Y_eigenvalue(mu,la)
            q^12*v*vl^12*v0
            sage: HH.Y_eigenvalue(mu,la,3)
            q^9*v*vl^12*v0
        """
        if la is None:
            la = self.lattice().zero()
        ulainv = self._Wv.from_reduced_word(la.reduced_word(positive=False))
        ev = self._dat.q() ** (- self.Y_pair_X_m(mu, la, a)) * (self._dat.parameter('short') ** (ulainv.action(self._two_rho_short).scalar(mu)))
        zed = self._two_rho_short.parent().zero()
        if self._two_rho_long != zed:
            ev = ev * (self._dat.parameter('long')  **ulainv.action(self._two_rho_long).scalar(mu))
        if self._two_rho_zero != zed:
            ev = ev * (self._dat.parameter('zero')  **ulainv.action(self._two_rho_zero).scalar(mu))
        return ev

    @cached_method
    def T_eigenvalue(self, u):
        r"""
        The eigenvalue of `T_u` on the generator of the polynomial module, where `u` is an element of the finite Weyl group.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
            sage: u = HH.double_affine_type().extended_affine_weyl().dual_classical_weyl().an_element(); u
            s1*s2
            sage: HH.T_eigenvalue(u)
            v^2

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2", untwisted=False, reduced=False, dual_reduced=False)
            sage: u = HH.double_affine_type().extended_affine_weyl().dual_classical_weyl().an_element(); u
            s1*s2
            sage: HH.T_eigenvalue(u)
            v*vl

        """
        i = u.first_descent()
        if i is None:
            return self.base_ring().one()
        return self._dat.q1(i) * self.T_eigenvalue(u.apply_simple_reflection(i))

    def lattice_algebra(self):
        r"""
        The group algebra of the lattice `L`

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2")
            sage: KL = HH.lattice_algebra(); KL
            Group algebra of the Ambient space of the Root system of type ['B', 2] over Fraction Field of Multivariate Polynomial Ring in q, v, vl over Rational Field
            sage: KL.an_element()
            X[(2, 2)]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2",dual_side=True)
            sage: KL = HH.lattice_algebra(); KL
            Group algebra of the Ambient space of the Root system of type ['B', 2] over Fraction Field of Multivariate Polynomial Ring in q, v, vl over Rational Field
            sage: KL.an_element()
            Y[(2, 2)]
        """
        return self._KL

    def Xu_on_poly_gen(self, u):
        r"""
        The action of the elements `X^u` for `u\in W(\tilde{Y})`
        on the generator of the polynomial module of the DAHA.

        .. RUBRIC:: Definition of `X^u`

        Let `H_Y` be the extended affine Hecke algebra (AHA) generated over the base ring `K`
        by `T_i^Y` for `i\in I^Y` and `\pi^Y\in \Pi^Y`.
        Its quadratic relations involve the dual Hecke parameters `v_{\alpha_i}^Y\in K` for `i\in I^Y`.
        Let `u=\pi^Y s_{i_1}\dotsm s_{i_\ell}` be any factorization of `u\in W(\tilde{Y})`
        with `\pi^Y\in \Pi^Y` and `i_1,\dotsc,i_\ell\in I^Y`. Define the sign `\epsilon_j \in \{ \pm 1\}`
        by the sign in `\pi^Y s_{i_1} \dotsm s_{i_{j-1}} \alpha_{i_j}^Y \in \mathbb{Z} \delta^Y \pm R_+(Y)`,
        where `R_+(Y)` is the set of positive roots for `Y`. Define `X^u \in H_Y` by

        .. MATH::

            X^u = \pi^Y T_{i_1}^{Y \epsilon_1} \dotsm T_{i_\ell}^{Y \epsilon_\ell}.

        The element `X^u` does not depend on the factorization. `H_Y` is a free `K[q^{\pm 1}]`-module with basis
        `X^u` for `u\in W(\tilde{Y}) = X \rtimes W_0`. Writing `u = t_\lambda w` with `\lambda \in X` and
        w\in W_0`, in the polynomial module we have

       .. MATH::

            X^u 1 = v_w X^\lambda

        with `v_w` the eigenvalue of the action of $T_w$ on $1$ in the polynomial module.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['C',2])
            sage: dat = HH.double_affine_type()
            sage: ddat = dat.dual()
            sage: w = ddat.extended_affine_weyl().PvW0().an_element(); w
            t[2*Lambda[1] + 2*Lambda[2]] * s1*s2
            sage: w.parent()
            Extended affine Weyl group of type ['B', 2, 1] realized by Semidirect product of Multiplicative form of Weight lattice of the Root system of type ['C', 2] acted upon by Weyl Group of type ['C', 2] (as a matrix group acting on the weight lattice)
            sage: HH.Xu_on_poly_gen(w)
            v*vl*X[(4, 2)]
            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',2],untwisted=False,reduced=False,dual_reduced=False)
            sage: dat = HH.double_affine_type()
            sage: ddat = dat.dual()
            sage: w = ddat.extended_affine_weyl().PvW0().an_element(); w
            t[2*Lambda[1] + 2*Lambda[2]] * s1*s2
            sage: w.parent()
            Extended affine Weyl group of type ['C', 2, 1]^* realized by Semidirect product of Multiplicative form of Weight lattice of the Root system of type ['B', 2] acted upon by Weyl Group of type ['B', 2] (as a matrix group acting on the weight lattice)
            sage: HH.Xu_on_poly_gen(w)
            v*vl*X[(3, 1)]
        """
        return self.T_eigenvalue(u.to_classical_weyl())*self.lattice_algebra().monomial(self.lattice()(u.to_dual_translation_left()))

    @cached_method
    def lattice_algebra_simple_roots(self):
        r"""
        The family of affine simple roots embedded into the lattice algebra.

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").lattice_algebra_simple_roots()
            Finite family {0: q^3*X[(-1, 0, 1)], 1: X[(1, -1, 0)], 2: X[(0, 1, -1)]}
            sage: DoubleAffineHeckeAlgebraSansDuality("B2",untwisted=False,reduced=False,dual_reduced=False).lattice_algebra_simple_roots()
            Finite family {0: q*X[(-1, 0)], 1: X[(1, -1)], 2: X[(0, 1)]}
        """
        dat = self.double_affine_type()
        X = self.lattice()
        simple_roots = X.simple_roots()   # get the classical simple roots
        I0 = simple_roots.keys()
        affine_cartan_type = dat.cartan_type()
        special_node = affine_cartan_type.special_node()
        assert special_node is 0
        col_ann = affine_cartan_type.col_annihilator()
        cl_alpha_0 = - X.sum([col_ann[i]*simple_roots[i] for i in I0]) # the classical projection of the zeroth simple root
        KL = self.lattice_algebra()
        return Family(dict([[0,dat.q()**dat._m * KL.monomial(cl_alpha_0)]] + [[i, KL.monomial(simple_roots[i])] for i in I0]))

    def dual_lattice_algebra(self):
        r"""
        The group algebra of the dual lattice.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2")
            sage: KLv = HH.dual_lattice_algebra(); KLv
            Group algebra of the Ambient space of the Root system of type ['C', 2] over Fraction Field of Multivariate Polynomial Ring in q, v, vl over Rational Field
            sage: KLv.an_element()
            Y[(2, 2)]
        """
        return self._KLv

    def F(self):
        r"""
        The fundamental group.
        """
        return self._F

    @cached_method
    def _F_to_W0Pv(self, pi):
        r"""
        Given an element `\pi` of the fundamental group, returns `(u, \mu)` where `u` is an element of the
        classical Weyl group and `\mu` is in the ambient space of the dual lattice.
        """
        x = self._W0Pv(pi)
        return (self._W.from_reduced_word(x.to_dual_classical_weyl().reduced_word()), x.to_dual_translation_right().to_ambient())

    @cached_method
    def _F_on_L(self, pi):
        r"""
        Returns a function which sends an element `\mu` of the weight lattice
        to its image under the action of the fundamental group element `\pi`.

        The value lies in the group algebra of the weight lattice. If `\pi=\pi_i` then

        ..MATH::

            \pi_i^X = u_i^{-1} t_{-\omega^Y_{i^*}}

        where ` - \omega^Y_{i^*} = w_0(\omega^Y_i) `

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
            sage: F = HH.double_affine_type().extended_affine_weyl().fundamental_group()
            sage: print "q should be replaced by q**(%s)"%(1/HH._m)
            q should be replaced by q**(1/3)
            sage: for i in F.special_nodes():
            ...       for j in HH.cartan_type().classical().index_set():
            ...           print i, j, HH._F_on_L(F(i))(HH._L.fundamental_weight(j))
            0 1 X[(1, 0, 0)]
            0 2 X[(1, 1, 0)]
            1 1 q*X[(0, 1, 0)]
            1 2 q^2*X[(0, 1, 1)]
            2 1 q^2*X[(0, 0, 1)]
            2 2 q*X[(1, 0, 1)]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2",dual_side=True)
            sage: F = HH.double_affine_type().extended_affine_weyl().fundamental_group()
            sage: print "q should be replaced by q**(%s)"%(1/HH._m)
            q should be replaced by q**(1/3)
            sage: for i in F.special_nodes():
            ...       for j in HH.cartan_type().classical().index_set():
            ...           print i, j, HH._F_on_L(F(i))(HH._L.fundamental_weight(j))
            0 1 Y[(1, 0, 0)]
            0 2 Y[(1, 1, 0)]
            1 1 q*Y[(0, 1, 0)]
            1 2 q^2*Y[(0, 1, 1)]
            2 1 q^2*Y[(0, 0, 1)]
            2 2 q*Y[(1, 0, 1)]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B4",dual_side=True)
            sage: F = HH.double_affine_type().extended_affine_weyl().fundamental_group()
            sage: print "q should be replaced by q**(%s)"%(1/HH._m)
            q should be replaced by q**(1/2)
            sage: for i in F.special_nodes():
            ...       for j in HH.cartan_type().classical().index_set():
            ...           print i, j, HH._F_on_L(F(i))(HH._L.fundamental_weight(j))
            0 1 Y[(1, 0, 0, 0)]
            0 2 Y[(1, 1, 0, 0)]
            0 3 Y[(1, 1, 1, 0)]
            0 4 Y[(1/2, 1/2, 1/2, 1/2)]
            1 1 q^2*Y[(-1, 0, 0, 0)]
            1 2 q^2*Y[(-1, 1, 0, 0)]
            1 3 q^2*Y[(-1, 1, 1, 0)]
            1 4 q*Y[(-1/2, 1/2, 1/2, 1/2)]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2", general_linear=True)
            sage: F = HH.double_affine_type().extended_affine_weyl().fundamental_group()
            sage: print "q should be replaced by q**(%s)"%(1/HH._m)
            q should be replaced by q**(1/3)
            sage: for i in range(5):
            ...       for j in HH.cartan_type().classical().index_set():
            ...           print i, j, HH._F_on_L(F(i))(HH._L.fundamental_weight(j))
            0 1 X[(1, 0, 0)]
            0 2 X[(1, 1, 0)]
            1 1 X[(0, 1, 0)]
            1 2 X[(0, 1, 1)]
            2 1 X[(0, 0, 1)]
            2 2 1/q^3*X[(1, 0, 1)]
            3 1 1/q^3*X[(1, 0, 0)]
            3 2 1/q^6*X[(1, 1, 0)]
            4 1 1/q^3*X[(0, 1, 0)]
            4 2 1/q^6*X[(0, 1, 1)]

            TODO: Fix the powers of q.

        """
        u, nu = self._F_to_W0Pv(pi)
        if not self.dual_side():
            nu = - nu
        if self.double_affine_type().general_linear():
            return lambda mu: self.lattice_algebra().term(pi.act_on_classical_ambient(mu), self._the_q_unscaled**(self.Y_pair_X_m(nu,mu)))
        else:
            return lambda mu: self.lattice_algebra().term(u.action(mu), self._the_q_unscaled**(self.Y_pair_X_m(nu, mu)))

    @cached_method
    def s0_on_L(self, mu):
        r"""
        Act on the element `\mu` of `L` by `s_0`.

        TODO:: CHECK!!!!

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
            sage: print "q should be replaced by q**(1/%s)"%HH._m
            q should be replaced by q**(1/3)
            sage: [(i, HH.s0_on_L(HH._L.fundamental_weight(i))) for i in HH.cartan_type().classical().index_set()]
            [(1, q^3*X[(0, 0, 1)]), (2, q^3*X[(0, 1, 1)])]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2")
            sage: print "q should be replaced by q**(1/%s)"%HH._m
            q should be replaced by q**(1/2)
            sage: [(i, HH.s0_on_L(HH._L.fundamental_weight(i))) for i in HH.cartan_type().classical().index_set()]
            [(1, q^2*X[(0, -1)]), (2, q^2*X[(-1/2, -1/2)])]            

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2", dual_side=True)
            sage: print "q should be replaced by q**(1/%s)"%HH._m
            q should be replaced by q**(1/2)
            sage: [(i, HH.s0_on_L(HH._L.fundamental_weight(i))) for i in HH.cartan_type().classical().index_set()]
            [(1, 1/q^2*Y[(0, -1)]), (2, 1/q^2*Y[(-1/2, -1/2)])]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("C2", dual_side=True)
            sage: print "q should be replaced by q**(1/%s)"%HH._m
            q should be replaced by q**(1/2)
            sage: [(i, HH.s0_on_L(HH._L.fundamental_weight(i))) for i in HH.cartan_type().classical().index_set()]
            [(1, 1/q^2*Y[(-1, 0)]), (2, 1/q^2*Y[(-1, 1)])]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2",untwisted=False)
            sage: print "q should be replaced by q**(1/%s)"%HH._m
            q should be replaced by q**(1/2)
            sage: [(i, HH.s0_on_L(HH._L.fundamental_weight(i))) for i in HH.cartan_type().classical().index_set()]
            [(1, q^4*X[(-1, 0)]), (2, q^2*X[(-1/2, 1/2)])]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2",untwisted=False, dual_side=True)
            sage: print "q should be replaced by q**(1/%s)"%HH._m
            q should be replaced by q**(1/2)
            sage: [(i, HH.s0_on_L(HH._L.fundamental_weight(i))) for i in HH.cartan_type().classical().index_set()]
            [(1, 1/q^4*Y[(-1, 0)]), (2, 1/q^2*Y[(-1/2, 1/2)])]

        """
        dat = self.double_affine_type()
        if dat.cartan_type().is_untwisted_affine():
            phi = dat.cartan_type().classical().root_system().root_space().highest_root().to_ambient()
            phiv = phi.associated_coroot()
        else:
            phiv = dat.cartan_type().classical().root_system().coroot_space().highest_root().to_ambient()
            phi = phiv.associated_coroot()
        sp = mu.scalar(phiv)
        return self._the_q_unscaled**(self._m * sp) * self._KL.monomial(mu - sp * phi)

    @cached_method
    def s0_on_lattice_algebra(self):
        r"""
        The operator that acts on the group algebra of `L` by `s_0`.

        TODO:: FIX!!!!!!!!!!!!

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
            sage: print "q should be replaced by q**(1/%s)"%HH._m
            q should be replaced by q**(1/3)
            sage: [(i, HH.s0_on_lattice_algebra()(HH._KL.monomial(HH._L.fundamental_weight(i)))) for i in HH.cartan_type().classical().index_set()]
            [(1, q^3*X[(0, 0, 1)]), (2, q^3*X[(0, 1, 1)])]
            sage: [(i, HH.s0_on_lattice_algebra()(HH.s0_on_L(HH._L.fundamental_weight(i)))) for i in HH.cartan_type().classical().index_set()]
            [(1, X[(1, 0, 0)]), (2, X[(1, 1, 0)])]
            sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2", dual_side=True)
            sage: print "q should be replaced by q**(1/%s)"%HH._m
            q should be replaced by q**(1/2)
            sage: [(i, HH.s0_on_lattice_algebra()(HH._KL.monomial(HH._L.fundamental_weight(i)))) for i in HH.cartan_type().classical().index_set()]
            [(1, 1/q^2*Y[(0, -1)]), (2, 1/q^2*Y[(-1/2, -1/2)])]
            sage: [(i, HH.s0_on_lattice_algebra()(HH.s0_on_L(HH._L.fundamental_weight(i)))) for i in HH.cartan_type().classical().index_set()]
            [(1, Y[(1, 0)]), (2, Y[(1/2, 1/2)])]
        """
        return self._KL.module_morphism(on_basis=self.s0_on_L, codomain=self._KL, category=ModulesWithBasis(self.base_ring()))

    def T(self):
        r"""
        The extended affine Hecke algebra corresponding to the lattice `L`.

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").T()
            T basis of The affine Hecke algebra of type ['A', 2, 1]
        """
        return self._T

    def LT(self):
        r"""
        The "LT" basis.

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").LT()
            LT basis of The double affine Hecke algebra of type ['A', 2, 1]

        """
        return self.DoubleAffineHeckeAlgebraSansDualityLT()

    def TL(self):
        r"""
        The "TL" basis.

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").TL()
            TL basis of The double affine Hecke algebra of type ['A', 2, 1]

        """
        return self.DoubleAffineHeckeAlgebraSansDualityTL()

    def LtvLv(self):
        r"""
        The LtvLv basis.

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").LtvLv()
            LtvLv basis of The double affine Hecke algebra of type ['A', 2, 1]

        """
        return self.DoubleAffineHeckeAlgebraSansDualityLtvLv()

    def LLvtv(self):
        r"""
        The LLvtv basis.

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").LLvtv()
            LLvtv basis of The double affine Hecke algebra of type ['A', 2, 1]

        """
        return self.DoubleAffineHeckeAlgebraSansDualityLLvtv()

    def tv_Lv(self):
        return self._tvLv

    def Lv_tv(self):
        return self._Lvtv

    def tv(self):
        return self._tv

    @cached_method
    def _F_on_LTmod_left_morphism(self, pi):
        r"""
        Returns the morphism that gives the action of the fundamental group element `pi` on the product form of the "LT" basis.
        """
        cat = ModulesWithBasis(self.base_ring())
        on_basis = lambda (mu,pi0,w): tensor([self._F_on_L(pi)(mu), self._T.product_by_fundamental_group_element_on_basis((pi0,w),pi,side='left')],category=cat.TensorProducts())
        return self._LTmod.module_morphism(on_basis=on_basis, category=cat, codomain = self._LTmod)

    @cached_method
    def _F_on_TLmod_right_morphism(self, pi):
        r"""
        Returns the morphism that gives the right action of the fundamental group element `pi` on the product form of the "TL" basis.
        """
        cat = ModulesWithBasis(self.base_ring())
        on_basis = lambda (pi0,w,mu): tensor([self._T.product_by_fundamental_group_element_on_basis((pi0,w),pi,side='right'),self._F_on_L(pi.inverse())(mu)],category=cat.TensorProducts())
        return self._TLmod.module_morphism(on_basis=on_basis, category=cat, codomain = self._TLmod)

    def a_realization(self):
        r"""
        Returns the default realization.

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").a_realization()
            LT basis of The double affine Hecke algebra of type ['A', 2, 1]

        """
        return self.LT()

    def _repr_(self):
        r"""
        A string representing ``self``.

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2")
            The double affine Hecke algebra of type ['A', 2, 1]
        """
        return "The double affine Hecke algebra of type %s"%self.cartan_type()


    @cached_method
    def is_koorn(self):
        r"""
        Is this a Koornwinder DAHA?
        """
        dat = self.double_affine_type()
        if dat.cartan_type_classical().letter != 'B':
            return False
        if dat.untwisted():
            return False
        if dat.reduced():
            return False
        if dat.dual_reduced():
            return False
        return True

    @cached_method
    def one_parameter_map(self):
        r"""
        The map which specializes all Hecke parameters to a single one.

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality(['B',4],reduced=False).one_parameter_map()
            Ring endomorphism of Fraction Field of Multivariate Polynomial Ring in q, v, vl, v2 over Rational Field Defn: q |--> q v |--> v vl |--> v v2 |--> v
        """
        R = self.base_ring()
        dat = self.double_affine_type()
        v = dat.parameter('short')
        images = [dat.parameter('null_root')] + [v for x in dat.hecke_parameter()]
        if dat.extended_base_ring():
            images = images + [-v+1/v for x in dat.hecke_parameter()]
        return R.hom(images, codomain = R)

    def one_parameter(self, f):
        r"""
        Specialize the coefficients of the element `f` of the polynomial module, to have a single
        Hecke parameter.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',2],reduced=False)
            sage: f = HH.ram_yip_nonsymmetric_macdonald_polynomial((0,1)); f
            ((-q^2*v*vl^2*v2^2+q^2*v*vl^2)/(-q^2*v^2*vl^2*v2+v2))*X[(0, 0)] + ((-vl^2+1)/(-q^2*v^2*vl^2+1))*X[(1, 0)] + X[(0, 1)]
            sage: HH.one_parameter(f)
            ((-q^2*v^4+q^2*v^2)/(-q^2*v^4+1))*X[(0, 0)] + ((-v^2+1)/(-q^2*v^4+1))*X[(1, 0)] + X[(0, 1)]
            sage: dat = HH.double_affine_type()
            sage: E = NonSymmetricMacdonaldPolynomials(['B',2,1],q=dat.parameter('null_root')**2,q1=dat.parameter('short'),q2=-1/dat.parameter('short'))
            sage: E[E.L0()((0,1))]
            ((-q^2*v^4+q^2*v^2)/(-q^2*v^4+1))*B[(0, 0)] + ((-v^2+1)/(-q^2*v^4+1))*B[(1, 0)] + B[(0, 1)]
        """
        return f.map_coefficients(self.one_parameter_map())

    def normalize_ambient_A_poly(self, f):
        r"""
        Normalize weights in a group algebra over the ambient weight lattice.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['A',2],reduced=False)
            sage: R = HH.base_ring()
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(),R.gen(0),R.gen(1),-1/R.gen(1))
            sage: L0 = E.L0()
            sage: E[L0((1,1,1))]
            B[(1, 1, 1)]
            sage: HH.normalize_ambient_A_poly(_)
            B[(0, 0, 0)]
        """
        def normalize_ambient_A_wt(la):
            v = la.to_vector()
            return la - la.parent().from_vector(vector([v[-1]]*len(v)))
        return f.map_support(normalize_ambient_A_wt)

    def koorn_to_D2_red(self, f):
        r"""
        Specialization from Koornwinder polynomial to the Macdonald polynomial of type `D_{n+1}^{(2)}`.
        """
        q,v,vl,v0,v2,vz = self.base_ring().gens()
        return f.map_coefficients(lambda x: x.subs({vl:v,v0:v,v2:v,vz:v}))

    def koorn_to_A2dag_red(self, f):
        r"""
        Specialization from Koornwinder polynomial to the Macdonald polynomial of type `A_{2n}^{(2)\dagger}`.
        """
        R = self.base_ring()
        q,v,vl,v0,v2,vz = R.gens()
        return f.map_coefficients(lambda x: x.subs({vl:v,v0:v,v2:v,vz:R.one()}))

    def koorn_to_C(self, f):
        r"""
        Specialization from Koornwinder polynomial to the Macdonald polynomial of type `C_n^{(1)}`.
        """
        R = self.base_ring()
        q,v,vl,v0,v2,vz = R.gens()
        return self.dualize_support(f.map_coefficients(lambda x: x.subs({vl:v,v0:v,v2:R.one(),vz:R.one()})))

    def koorn_to_A2_red(self, f):
        r"""
        Specialization from Koornwinder polynomial to the Macdonald polynomial of type `A_{2n}^{(2)}`.
        """
        R = f.base_ring()
        q,v,vl,v0,v2,vz = R.gens()
        return self.dualize_support(f.map_coefficients(lambda x: x.subs({vl:v,v0:v,v2:R.one(),vz:v})))

    @cached_method
    def dualize_morphism(self):
        r"""
        The morphism which changes a polynomial to have weights in the ambient lattice of dual type.
        Only designed for `B_n` and `C_n`.
        """
        ct = self.cartan_type().classical()
        dct = ct.dual()
        damb = dct.root_system().ambient_space()
        newalg = damb.algebra(self.base_ring())
        return self.lattice_algebra().module_morphism(on_basis=lambda x: newalg.monomial(damb(tuple(x.to_vector()))), codomain=newalg)

    def dualize_support(self, f):
         r"""
         Changes the polynomial to have weights in the ambient lattice of dual type.
         Only designed for `B_n` and `C_n`.
         """
         return self.dualize_morphism()(f)

    def ram_yip_nonsymmetric_macdonald_polynomial(self, weight, start=None):
        r"""
        Nonsymmetric Macdonald-Koornwinder polynomial via the Ram-Yip formula.

        INPUTS::

            - ``weight`` -- Tuple that will be converted into a weight in the ambient space of `X`.
            - ``start`` -- An element of the (extended) affine Weyl group element of type `\tilde{Y}`.
            This can be an actual element or a reduced word for such an element.

        Let `\lambda \in X` be the weight given by the input tuple ``weight`` and `u \in W(\tilde{Y})` the
        affine Weyl group element given by ``start``. This function returns the value of

        .. MATH::

            X^u E_\lambda

        where `X^u` was defined in :meth:`Xu_on_poly_gen` and `E_\lambda` is the nonsymmetric Macdonald-Koornwinder polynomial.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',2],untwisted=False,reduced=False,dual_reduced=False)
            sage: HH.ram_yip_nonsymmetric_macdonald_polynomial((1,0))
            ((-q^2*v*vl^2*v0^2*v2^2*vz+q^2*v*vl^2*v0^2*vz-q*v0*v2*vz^2+q*v0*v2)/(-q^2*v^2*vl^4*v0^2*v2*vz+v2*vz))*X[(0, 0)] + X[(1, 0)]
            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',2],untwisted=False,reduced=False,dual_reduced=False, extra_parameters=True)
            sage: HH.ram_yip_nonsymmetric_macdonald_polynomial((1,0))
            ((q^2*v*vl^2*v0^2*c2+q*v0*cz)/(-q^2*v^2*vl^4*v0^2+1))*X[(0, 0)] + X[(1, 0)]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['A',1])
            sage: HH.ram_yip_nonsymmetric_macdonald_polynomial((0,2))
            ((-q^2*v^2+q^2-v^2+1)/(-q^4*v^2+1))*X[(0, 0)] + X[(-2, 0)] + ((-v^2+1)/(-q^4*v^2+1))*X[(2, 0)]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['A',1], dual_reduced=False)
            sage: HH.ram_yip_nonsymmetric_macdonald_polynomial((2,0))
            ((-q^4*v^2*v0^2+q^4*v0^2-q^2*v0^2+q^2)/(-q^4*v^2*v0^2+1))*X[(0, 0)] + X[(2, 0)]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['A',1], reduced=False, dual_reduced=False)
            sage: tup = (2,0)
            sage: f = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); f
            ((-q^2*v*v0^2*v2^2*vz+q^2*v*v0^2*vz-q*v0*v2*vz^2+q*v0*v2)/(-q^2*v^2*v0^2*v2*vz+v2*vz))*X[(0, 0)] + X[(2, 0)]
            sage: fo = HH.one_parameter(f); fo
            ((-q*v^2+q)/(-q*v^2+1))*X[(0, 0)] + X[(2, 0)]
            sage: g = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup)
            sage: dat = HH.double_affine_type()
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), dat.q(), dat.parameter('short'), -1/dat.parameter('short'))
            sage: HH.normalize_ambient_A_poly(HH._from_generic_lattice_algebra(E[E.L0()(tup)])-fo)
            0
            sage: tup = (2,0)
            sage: f = HH.one_parameter(HH.ram_yip_nonsymmetric_macdonald_polynomial(tup))
            sage: HH.normalize_ambient_A_poly(f - HH._from_generic_lattice_algebra(E[E.L0()(tup)]))
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['C',2])
            sage: tup = (0,1)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((-v^2+1)/(-q^2*v^2*vl^4+1))*X[(1, 0)] + X[(0, 1)]
            sage: R = HH.lattice_algebra().base_ring()
            sage: q = R.gen(0); v = R.gen(1)
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), q, v, -1/v)
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(E[E.L0()(tup)]).map_coefficients(lambda x: x.subs({q:q*q}))
            0
            sage: tup = (-1,1)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup)
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(E[E.L0()(tup)]).map_coefficients(lambda x: x.subs({q:q*q}))
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['C',3])
            sage: tup = (0,1,0)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((-v^2+1)/(-q^2*v^6*vl^4+1))*X[(1, 0, 0)] + X[(0, 1, 0)]
            sage: f = HH.one_parameter(nsry); f
            ((-v^2+1)/(-q^2*v^10+1))*X[(1, 0, 0)] + X[(0, 1, 0)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), q, v, -1/v)
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-v^2+1)/(-q*v^10+1))*B[(1, 0, 0)] + B[(0, 1, 0)]
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(nsm).map_coefficients(lambda x: x.subs({q:q*q}))
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['F',4])
            sage: tup = (0,1,0,0)
            sage: f = HH.one_parameter(HH.ram_yip_nonsymmetric_macdonald_polynomial(tup)); f
            ((q*v^8-q*v^6)/(q*v^12-1))*X[(0, 0, 0, 0)] + ((-v^2+1)/(-q*v^12+1))*X[(1/2, 1/2, -1/2, -1/2)] + ((-v^2+1)/(-q*v^12+1))*X[(1/2, 1/2, -1/2, 1/2)] + ((-v^2+1)/(-q*v^12+1))*X[(1/2, 1/2, 1/2, -1/2)] + ((-v^2+1)/(-q*v^12+1))*X[(1/2, 1/2, 1/2, 1/2)] + ((v^2-1)/(q*v^12-1))*X[(1, 0, 0, 0)] + X[(0, 1, 0, 0)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), q, v, -1/v)
            sage: f - HH._from_generic_lattice_algebra(E[E.L0()(tup)])
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['G',2])
            sage: tup = (1,0,-1)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((-q*v^2*vl^4+q*vl^4)/(-q*v^4*vl^6+1))*X[(0, 0, 0)] + X[(1, 0, -1)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), q, v, -1/v)
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-q*v^6+q*v^4)/(-q*v^10+1))*B[(0, 0, 0)] + B[(1, 0, -1)]
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(nsm)
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['C',2],untwisted=False)
            sage: tup = (0,2)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((-q^6*v^4*vl^2+q^6*v^4-q^4*v^4+q^4*v^2*vl^2+q^4*v^2-q^2*v^4-q^4+2*q^2*v^2-q^2)/(-q^6*v^4*vl^2+q^4*v^2*vl^2+q^2*v^2-1))*X[(0, 0)] + ((-q^2*v^2+q^2)/(-q^2*v^2+1))*X[(-1, 1)] + ((q^2*v^4-2*q^2*v^2+q^2)/(q^6*v^4*vl^2-q^4*v^2*vl^2-q^2*v^2+1))*X[(1, -1)] + ((-q^4*v^4*vl^2+q^4*v^2*vl^2-q^2*v^4+2*q^2*v^2-q^2+v^2-1)/(-q^6*v^4*vl^2+q^4*v^2*vl^2+q^2*v^2-1))*X[(1, 1)] + ((-v^2+1)/(-q^4*v^2*vl^2+1))*X[(2, 0)] + X[(0, 2)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), q, v, -1/v)
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-q^3*v^6+q^3*v^4-q*v^4+q^2*v^2+2*q*v^2-q^2-q)/(-q^3*v^6+q^2*v^4+q*v^2-1))*B[(0, 0)] + ((-q*v^2+q)/(-q*v^2+1))*B[(-1, 1)] + ((-q*v^4+2*q*v^2-q)/(-q^3*v^6+q^2*v^4+q*v^2-1))*B[(1, -1)] + ((q^2*v^6-q^2*v^4+q*v^4-2*q*v^2-v^2+q+1)/(q^3*v^6-q^2*v^4-q*v^2+1))*B[(1, 1)] + ((-v^2+1)/(-q^2*v^4+1))*B[(2, 0)] + B[(0, 2)]
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(nsm).map_coefficients(lambda y: y.subs({q:q*q}))
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['F',4],untwisted=False)
            sage: tup = (0,1,0,0)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((q^3*v^16*vl^10-q^3*v^14*vl^10+q^2*v^10*vl^4-q^2*v^8*vl^6-q^2*v^8*vl^4+q^2*v^6*vl^6-q*v^2+q)/(q^3*v^16*vl^14-q^2*v^12*vl^10-q*v^4*vl^4+1))*X[(0, 0, 0, 0)] + ((-v^2+1)/(-q*v^4*vl^4+1))*X[(1/2, 1/2, -1/2, -1/2)] + ((-v^2+1)/(-q*v^4*vl^4+1))*X[(1/2, 1/2, -1/2, 1/2)] + ((-v^2+1)/(-q*v^4*vl^4+1))*X[(1/2, 1/2, 1/2, -1/2)] + ((-v^2+1)/(-q*v^4*vl^4+1))*X[(1/2, 1/2, 1/2, 1/2)] + ((q*v^12*vl^6-q*v^10*vl^6-q*v^8*vl^4+q*v^6*vl^6+v^6-v^4*vl^2-v^2+1)/(q^3*v^16*vl^14-q^2*v^12*vl^10-q*v^4*vl^4+1))*X[(1, 0, 0, 0)] + X[(0, 1, 0, 0)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), q, v, -1/v)
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-q^3*v^26+q^3*v^24+q*v^2-q)/(-q^3*v^30+q^2*v^22+q*v^8-1))*B[(0, 0, 0, 0)] + ((-v^2+1)/(-q*v^8+1))*B[(1/2, 1/2, -1/2, -1/2)] + ((-v^2+1)/(-q*v^8+1))*B[(1/2, 1/2, -1/2, 1/2)] + ((-v^2+1)/(-q*v^8+1))*B[(1/2, 1/2, 1/2, -1/2)] + ((v^2-1)/(q*v^8-1))*B[(1/2, 1/2, 1/2, 1/2)] + ((-q*v^18+q*v^16+v^2-1)/(-q^3*v^30+q^2*v^22+q*v^8-1))*B[(1, 0, 0, 0)] + B[(0, 1, 0, 0)]
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(nsm)
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['G',2],untwisted=False)
            sage: tup = (0,1,-1)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((q^3*v^6*vl^4-q^3*v^4*vl^4+q^2*v^4*vl^2-q^2*v^2*vl^2+q*v^2-q)/(q^3*v^6*vl^4-1))*X[(0, 0, 0)] + ((-vl^2+1)/(-q^3*v^6*vl^4+1))*X[(1, -1, 0)] + ((-q^2*v^6*vl^4+q^2*v^4*vl^4-q*v^4*vl^2+q*v^2*vl^2-v^2+1)/(-q^3*v^6*vl^4+1))*X[(1, 0, -1)] + X[(0, 1, -1)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), q, v, -1/v)
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-q^3*v^10+q^3*v^8-q^2*v^6+q^2*v^4-q*v^2+q)/(-q^3*v^10+1))*B[(0, 0, 0)] + ((-v^2+1)/(-q^3*v^10+1))*B[(1, -1, 0)] + ((-q^2*v^10+q^2*v^8-q*v^6+q*v^4-v^2+1)/(-q^3*v^10+1))*B[(1, 0, -1)] + B[(0, 1, -1)]
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(nsm)
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['A',1],reduced=False,dual_reduced=False)
            sage: tup = (1,-1)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((-q^2*v*v0^2*v2^2*vz+q^2*v*v0^2*vz-q*v0*v2*vz^2+q*v0*v2)/(-q^2*v^2*v0^2*v2*vz+v2*vz))*X[(0, 0)] + X[(2, 0)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), dat.q(), dat.parameter('short'), -1/dat.parameter('short'))
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-q*v^2+q)/(-q*v^2+1))*B[(0, 0)] + B[(1, -1)]
            sage: HH.normalize_ambient_A_poly(HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(nsm))
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',3],reduced=False); dat = HH.double_affine_type()
            sage: tup = (0,1,0)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((-q^2*v*vl^4*v2^2+q^2*v*vl^4)/(-q^2*v^2*vl^6*v2+v2))*X[(0, 0, 0)] + ((-vl^2+1)/(-q^2*v^2*vl^6+1))*X[(1, 0, 0)] + X[(0, 1, 0)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), dat.q(), dat.parameter('short'), -1/dat.parameter('short'))
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-q*v^6+q*v^4)/(-q*v^8+1))*B[(0, 0, 0)] + ((-v^2+1)/(-q*v^8+1))*B[(1, 0, 0)] + B[(0, 1, 0)]
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(nsm).map_coefficients(lambda y: y.subs({dat.q():dat.q()**2}))
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['C',3],dual_reduced=False); dat = HH.double_affine_type()
            sage: tup = (0,1,0)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((-v^2+1)/(-q^2*v^6*vl^2*v0^2+1))*X[(1, 0, 0)] + X[(0, 1, 0)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), dat.q(), dat.parameter('short'), -1/dat.parameter('short'))
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-v^2+1)/(-q*v^10+1))*B[(1, 0, 0)] + B[(0, 1, 0)]
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(nsm).map_coefficients(lambda y: y.subs({dat.q():dat.q()**2}))
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['C',2],reduced=False,dual_reduced=False); dat=HH.double_affine_type()
            sage: tup = (0,1)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((-v*v2^2+v)/(-q*v^2*vl^2*v0^2*v2+v2))*X[(1, 0)] + X[(0, 1)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), dat.q(), dat.parameter('short'), -1/dat.parameter('short'))
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-v^2+1)/(-q*v^6+1))*B[(1, 0)] + B[(0, 1)]
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(nsm)
            0

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',2],untwisted=False,reduced=False,dual_reduced=False); dat = HH.double_affine_type()
            sage: tup = (-1,0)
            sage: nsry = HH.ram_yip_nonsymmetric_macdonald_polynomial(tup); nsry
            ((q^5*v^2*vl^6*v0*v2*vz^2+q^4*v*vl^6*v0^2*v2^2*vz-q^5*v^2*vl^6*v0*v2-q^4*v*vl^6*v0^2*vz-q^3*v^2*vl^2*v0*v2*vz^2+q^3*vl^4*v0*v2*vz^2+q^2*v*vl^4*v2^2*vz-q^2*v*vl^2*v0^2*v2^2*vz+q^3*v^2*vl^2*v0*v2-q^3*vl^4*v0*v2-q^2*v*vl^4*vz+q^2*v*vl^2*v0^2*vz-q*v0*v2*vz^2-v*v2^2*vz+q*v0*v2+v*vz)/(q^6*v^2*vl^6*v0^2*v2*vz-q^4*v^2*vl^4*v0^2*v2*vz-q^2*vl^2*v2*vz+v2*vz))*X[(0, 0)] + X[(-1, 0)] + ((q^4*v^2*vl^6*v0^2-q^4*v^2*vl^6+q^2*v^2*vl^6-q^2*v^2*vl^4*v0^2+q^2*v^2*vl^2-q^2*vl^4-v^2*vl^2+vl^4-vl^2+1)/(q^6*v^2*vl^6*v0^2-q^4*v^2*vl^4*v0^2-q^2*vl^2+1))*X[(1, 0)] + ((-vl^2+1)/(-q^2*vl^2+1))*X[(0, -1)] + ((-vl^2+1)/(-q^2*vl^2+1))*X[(0, 1)]
            sage: E = NonSymmetricMacdonaldPolynomials(HH.cartan_type(), dat.q(), dat.parameter('short'), -1/dat.parameter('short'))
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-q*v^2-v^2+q+1)/(-q^2*v^2+1))*B[(0, 0)] + B[(-1, 0)] + ((-v^2+1)/(-q^2*v^2+1))*B[(1, 0)] + ((-v^2+1)/(-q^2*v^2+1))*B[(0, -1)] + ((v^2-1)/(q^2*v^2-1))*B[(0, 1)]
            sage: HH.one_parameter(nsry) - HH._from_generic_lattice_algebra(nsm)
            0

            sage: E = NonSymmetricMacdonaldPolynomials(CartanType(['A',4,2]).dual(), dat.q(), dat.parameter('short'), -1/dat.parameter('short'))
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-v^2+1)/(-q^2*v^2+1))*B[(0, 0)] + B[(-1, 0)] + ((-v^2+1)/(-q^2*v^2+1))*B[(1, 0)] + ((-v^2+1)/(-q^2*v^2+1))*B[(0, -1)] + ((v^2-1)/(q^2*v^2-1))*B[(0, 1)]
            sage: HH.koorn_to_A2dag_red(nsry) - HH._from_generic_lattice_algebra(nsm)
            0

            sage: E = NonSymmetricMacdonaldPolynomials(CartanType(['C',2,1]), dat.q(), dat.parameter('short'), -1/dat.parameter('short'))
            sage: nsm = E[E.L0()(tup)]; nsm
            B[(-1, 0)] + ((-v^2+1)/(-q*v^2+1))*B[(1, 0)] + ((-v^2+1)/(-q*v^2+1))*B[(0, -1)] + ((v^2-1)/(q*v^2-1))*B[(0, 1)]
            sage: HH.koorn_to_C(nsry) - nsm.map_coefficients(lambda x: x.subs({dat.q():dat.q()**2}))
            0

            sage: E = NonSymmetricMacdonaldPolynomials(CartanType(['A',4,2]), dat.q(), dat.parameter('short'), -1/dat.parameter('short'))
            sage: nsm = E[E.L0()(tup)]; nsm
            ((-q*v^2+q)/(-q^2*v^2+1))*B[(0, 0)] + B[(-1, 0)] + ((-v^2+1)/(-q^2*v^2+1))*B[(1, 0)] + ((-v^2+1)/(-q^2*v^2+1))*B[(0, -1)] + ((v^2-1)/(q^2*v^2-1))*B[(0, 1)]
            sage: HH.koorn_to_A2_red(nsry)-nsm
            0
        """
        terms = RamYipFormula(self, weight, start=start).alcove_paths()
        return self.lattice_algebra().sum([x.evaluate() for x in terms])

    @cached_method
    def symmetric_macdonald_koornwinder_polynomial(self, weight):
        r"""
        The symmetric Macdonald-Koornwinder polynomial.

        This version is computed using the Ram-Yip formula.

        INPUTS::

            - ``weight`` -- a tuple

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['A',2])
            sage: HH.symmetric_macdonald_koornwinder_polynomial((1,0,0))
            X[(-1, -1, 0)] + X[(1, 0, 0)] + X[(0, 1, 0)]

        The type A ambient space has many different ways to represent the same SL weight.
        For better or for worse, the above weights have been normalized uniquely to have last part zero.

        Can we normalize terms so that they are in the same root lattice coset as ``weight``?
        This is the best normalization.

        EXAMPLES::

            sage: HH.symmetric_macdonald_koornwinder_polynomial((1,1,0))
            X[(-1, 0, 0)] + X[(1, 1, 0)] + X[(0, -1, 0)]

            To avoid fractional powers of variables in the double affine Hecke algebra,
            the Hecke parameters, often denoted by `t`, are represented here by `v^2`
            and the null root parameter `q`, has been replaced by a certain positive integer power of itself,
            which in this case is 3::

            sage: HH.double_affine_type()._m
            3

        EXAMPLES::

            sage: HH.symmetric_macdonald_koornwinder_polynomial((2,0,0))
            X[(-2, -2, 0)] + ((q^3*v^2-q^3+v^2-1)/(q^3*v^2-1))*X[(-1, 0, 0)] + ((q^3*v^2-q^3+v^2-1)/(q^3*v^2-1))*X[(1, 1, 0)] + X[(2, 0, 0)] + ((q^3*v^2-q^3+v^2-1)/(q^3*v^2-1))*X[(0, -1, 0)] + X[(0, 2, 0)]

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',2],untwisted=False,reduced=False,dual_reduced=False)
            sage: HH.symmetric_macdonald_koornwinder_polynomial((1,0))
            ((q^2*v*vl^4*v0^2*v2^2*vz+q*v^2*vl^4*v0*v2*vz^2-q^2*v*vl^4*v0^2*vz+q^2*v*vl^2*v0^2*v2^2*vz-q*v^2*vl^4*v0*v2+q*v^2*vl^2*v0*v2*vz^2-q^2*v*vl^2*v0^2*vz-q*v^2*vl^2*v0*v2+q*vl^2*v0*v2*vz^2+v*vl^2*v2^2*vz-q*vl^2*v0*v2+q*v0*v2*vz^2-v*vl^2*vz+v*v2^2*vz-q*v0*v2-v*vz)/(q^2*v^2*vl^4*v0^2*v2*vz-v2*vz))*X[(0, 0)] + X[(-1, 0)] + X[(1, 0)] + X[(0, -1)] + X[(0, 1)]
        """
        KX = self.lattice_algebra()
        terms = RamYipFormula(self, weight, symmetric=True).alcove_paths()
        return KX.sum([x.evaluate() for x in terms])

    @cached_method
    def orbit_sum(self, weight):
        r"""
        Return the orbit sum of ``weight``.

        INPUTS::

        - ``weight`` -- an element of the ambient space or a tuple representing such an element

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality(['A',2]).orbit_sum((2,1,0))
            X[(1, 2, 0)] + X[(1, 0, 2)] + X[(2, 1, 0)] + X[(2, 0, 1)] + X[(0, 1, 2)] + X[(0, 2, 1)]
            sage: DoubleAffineHeckeAlgebraSansDuality(['B',2]).orbit_sum((1,1))
            X[(-1, -1)] + X[(-1, 1)] + X[(1, -1)] + X[(1, 1)]
        """
        if isinstance(weight, tuple):
            weight = self.lattice()(weight)
        elif weight.parent() != self.lattice():
            raise TypeError("Must be a tuple or element of the ambient space")
        orbit = weight.orbit()
        return self.lattice_algebra().sum_of_monomials(orbit)

    def expansion_into_symmetric_macdonald_koornwinder_basis(self, f, **keywords):
        r"""
        Return the expansion of the symmetric element `f` into symmetric Macdonald-Koornwinder
        polynomials.

        INPUTS::

        - `f` -- An element to be expanded into symmetric Macdonald Koornwinder polynomials

        There is an optional keyword argument ``verbose`` which, if set to True, prints extra data.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['C',2])
            sage: HH.expansion_into_symmetric_macdonald_koornwinder_basis(HH.orbit_sum((2,0)))
            ((-q^4*v^4*vl^6+2*q^4*v^4*vl^4-q^4*v^2*vl^6-q^2*v^4*vl^6-q^4*v^2*vl^4+2*q^2*v^4*vl^4-q^2*v^2*vl^6+q^2*v^4*vl^2+q^4*vl^4-q^2*v^2*vl^4-q^2*v^2*vl^2+v^4*vl^2+q^2*vl^4-q^2*v^2+2*q^2*vl^2-v^2*vl^2-q^2-v^2+2*vl^2-1)/(q^4*v^4*vl^6-q^2*v^2*vl^4-q^2*v^2*vl^2+1))*X[(0, 0)] + ((-q^2*v^2+q^2-v^2+1)/(q^2*v^2-1))*X[(1, 1)] + X[(2, 0)]
            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',2],untwisted=False,reduced=False,dual_reduced=False)
            sage: HH.expansion_into_symmetric_macdonald_koornwinder_basis(HH.orbit_sum((1,0)))
            ((-q^2*v*vl^4*v0^2*v2^2*vz-q*v^2*vl^4*v0*v2*vz^2+q^2*v*vl^4*v0^2*vz-q^2*v*vl^2*v0^2*v2^2*vz+q*v^2*vl^4*v0*v2-q*v^2*vl^2*v0*v2*vz^2+q^2*v*vl^2*v0^2*vz+q*v^2*vl^2*v0*v2-q*vl^2*v0*v2*vz^2-v*vl^2*v2^2*vz+q*vl^2*v0*v2-q*v0*v2*vz^2+v*vl^2*vz-v*v2^2*vz+q*v0*v2+v*vz)/(q^2*v^2*vl^4*v0^2*v2*vz-v2*vz))*X[(0, 0)] + X[(1, 0)]

        Doesn't work for type A because of the nonuniqueness of representations of weights.
        Fix this!!!

        """
        verbose = keywords.get('verbose',False)
        lattice = self.lattice()
        lattice_algebra = self.lattice_algebra()
        assert f.parent() == lattice_algebra
        expansion = lattice_algebra.zero()
        while not f == lattice_algebra.zero():
            if verbose:
                print "current function: "
                print f
            dominant_weight = None
            for weight in f.support():
                if weight.is_dominant():
                    if verbose:
                        print "current dominant weight: %s"%weight
                    if dominant_weight is None:
                        dominant_weight = weight
                        if verbose:
                            print "first dominant weight: %s"%weight
                    else:
                        diff = weight - dominant_weight
                        if all([diff.scalar(lattice.simple_coroot(i)) >= 0 for i in self.cartan_type().classical().index_set()]):
                            if verbose:
                                print "replacing ", dominant_weight, " with ", weight
                            dominant_weight = weight
            if dominant_weight is None:
                raise(ValueError, "input is not symmetric")
            expansion = expansion + lattice_algebra.term(dominant_weight, f[dominant_weight])
            if verbose:
                print "expansion ", expansion
            f = f - f[dominant_weight] * self.symmetric_macdonald_koornwinder_polynomial(tuple(dominant_weight.to_vector()))
        return expansion

    class _DAHABasesCategory(Category_realization_of_parent):
        r"""
        The category of realizations of a double affine Hecke algebra (without duality)
        """
        def __init__(self, basis, prefix=None):
            r"""
            Initialize a basis.

            INPUT:

            - ``basis`` -- a basis
            - ``prefix`` -- a prefix

            TESTS::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: bases = HH._DAHABasesCategory()
                sage: HH.LT() in bases
                True
            """
            Category_realization_of_parent.__init__(self, basis)
            basis._prefix = prefix

        def super_categories(self):
            r"""
            EXAMPLES::

                sage: DoubleAffineHeckeAlgebraSansDuality("A2")._DAHABasesCategory().super_categories()
                [Category of realizations of The double affine Hecke algebra of type ['A', 2, 1], Category of algebras with basis over Fraction Field of Multivariate Polynomial Ring in q, v over Rational Field, Category of monoids with realizations, Category of additive unital additive magmas with realizations]

            """
            return [Realizations(self.base())] + self.base().category().super_categories()

        def _repr_(self):
            r"""
            Return the representation of ``self``.

            EXAMPLES::

                sage: DoubleAffineHeckeAlgebraSansDuality("A2")._DAHABasesCategory()
                Category of bases of The double affine Hecke algebra of type ['A', 2, 1]
            """
            return "Category of bases of %s" % self.base()

        class ParentMethods:

            def T_generators(self):
                r"""
                The family of elements `T_i` in the given realization.

                It is a family on the affine Dynkin node set `I` with values in ``self``.

                ..warning:: Must be implemented for "LT".

                EXAMPLES::

                    sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                    sage: HH.LT().T_generators()
                    Finite family {0: TX[0], 1: TX[1], 2: TX[2]}
                    sage: HH.TL().T_generators()
                    Finite family {0: TX[0], 1: TX[1], 2: TX[2]}
                    sage: HH.LtvLv().T_generators()
                    Finite family {0: Ty[1,2,1] Y[(-1, 0, 1)] + ((v^2-1)/v), 1: Ty[1], 2: Ty[2]}

                TODO:: other realizations

                """
                HH = self.realization_of()
                return Family(dict([[i, self(HH.LT().T_generators()[i])] for i in HH.double_affine_type().cartan_type().index_set()]))

            def from_fundamental(self, i):
                r"""
                The image of the fundamental group element of `L` indexed by the special node `i`, in ``self``.

                ..warning:: Must be implemented in "LT".

                EXAMPLES::

                    sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                    sage: F = HH._F
                    sage: [(i, HH.LT().from_fundamental(i)) for i in F.special_nodes()]
                    [(0, 1), (1, piX[1]), (2, piX[2])]
                    sage: HH.LT().from_fundamental(1).parent()
                    LT basis of The double affine Hecke algebra of type ['A', 2, 1]
                    sage: [(i, HH.LtvLv().from_fundamental(i)) for i in F.special_nodes()]
                    [(0, 1), (1, Ty[1,2] Y[(-1, -1, 0)]), (2, Ty[2,1] Y[(-1, 0, 0)])]
                """
                return self(self.realization_of().LT().from_fundamental(i))

            def from_lattice_algebra(self, a):
                r"""
                The image of `a` under the morphism from the group algebra of the `L` lattice into ``self``.

                ..warning:: Must be implemented in "LT".

                EXAMPLES::

                    sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                    sage: LT = HH.LT()
                    sage: KL = HH.lattice_algebra()
                    sage: a = KL.monomial(HH.lattice().simple_root(1)); a
                    X[(1, -1, 0)]
                    sage: a.parent()
                    Group algebra of the Ambient space of the Root system of type ['A', 2] over Fraction Field of Multivariate Polynomial Ring in q, v over Rational Field
                    sage: aa = LT.from_lattice_algebra(a); aa
                    X[(1, -1, 0)]
                    sage: aa.parent()
                    LT basis of The double affine Hecke algebra of type ['A', 2, 1]

                """
                HH = self.realization_of()
                return self(HH.LT().factor_embedding(0)(a))

            def from_dual_lattice_algebra(self, b):
                r"""
                The image of `b` under the morphism from the group algebra of the `Lv` lattice into ``self``.

                ..warning:: Must be implemented in basis `LtvLv`.

                EXAMPLES::

                    sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                    sage: LtvLv = HH.LtvLv()
                    sage: KLv = HH.dual_lattice_algebra()
                    sage: Lv = HH.dual_lattice()
                    sage: b = KLv.monomial(Lv.simple_root(1)); b
                    Y[(1, -1, 0)]
                    sage: b.parent()
                    Group algebra of the Ambient space of the Root system of type ['A', 2] over Fraction Field of Multivariate Polynomial Ring in q, v over Rational Field
                    sage: bb = LtvLv.from_dual_lattice_algebra(b); bb
                    Y[(1, -1, 0)]
                    sage: bb.parent()
                    LtvLv basis of The double affine Hecke algebra of type ['A', 2, 1]
                    sage: c = HH.LT().from_dual_lattice_algebra(b); c
                    TX[0,2,0,1] + ((-v^2+1)/v)*TX[0,2,1]
                    sage: LtvLv(c) == b
                    True
                    sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2", dual_side=True)
                    sage: LtvLv = HH.LtvLv()
                    sage: KLv = HH.dual_lattice_algebra()
                    sage: Lv = HH.dual_lattice()
                    sage: b = KLv.monomial(Lv.simple_root(1)); b
                    X[(1, -1, 0)]
                    sage: b.parent()
                    Group algebra of the Ambient space of the Root system of type ['A', 2] over Fraction Field of Multivariate Polynomial Ring in q, v over Rational Field
                    sage: c = HH.TL()(b); c
                    ((v^4-2*v^2+1)/v^2)*TY[2,0] + ((-v^2+1)/v)*TY[2,0,1] + ((-v^2+1)/v)*TY[0,2,0] + ((v^4-2*v^2+1)/v^2) + ((-v^2+1)/v)*TY[1] + TY[0,2,0,1]
                    sage: c.parent()
                    TL basis of The double affine Hecke algebra of type ['A', 2, 1]

                """
                HH = self.realization_of()
                return self(HH.LtvLv().factor_embedding(1)(HH.tv_Lv().factor_embedding(1)(b)))

            def from_extended_affine_hecke(self, a):
                r"""
                Returns the image of `a` from the extended affine Hecke algebra "T" to ``self``.

                ..warning:: Must be implemented in type "LT".

                EXAMPLES::

                    sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                    sage: LT = HH.LT()
                    sage: T = HH.T()
                    sage: a = T.an_element()
                    sage: b = LT.from_extended_affine_hecke(a); b
                    2*TX[0] + 3*TX[0,1] + 1 + TX[0,1,2] + 4*piX[1] TX[0] + 6*piX[1] TX[0,1] + 2*piX[1] + 2*piX[1] TX[0,1,2] + 8*piX[2] TX[0] + 12*piX[2] TX[0,1] + 4*piX[2] + 4*piX[2] TX[0,1,2]
                    sage: b == LT(a)
                    True
                """
                return self(self.realization_of().LT().from_extended_affine_hecke(a))

            def from_dual_classical_hecke(self, a):
                r"""
                Returns the image of `a` from the finite Hecke algebra "tv" of dual type, to ``self``.

                ..warning:: Must be implemented in type "LtvLv".

                EXAMPLES::

                    sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                    sage: LT = HH.LT()
                    sage: tv_Lv = HH.tv_Lv()
                    sage: tv = HH.tv()
                    sage: a = tv.an_element(); a
                    2*Ty[1,2,1] + 4*Ty[1,2] + 1
                    sage: b = LT.from_dual_classical_hecke(a); b
                    1 + 4*TX[1,2] + 2*TX[1,2,1]
                    sage: b == LT(a)
                    True
                """
                return self(self.realization_of().LtvLv().from_dual_classical_hecke(a))

            def signed_generator_product(self, word, signs=None):
                r"""
                Given a reduced word for an element `w` in the affine Weyl group, return the image of the basis element `T_w`
                in ``self``.

                The list ``signs``, if present, is a list of entries `+1` or `-1`. The `-1` indicates that
                the corresponding generator `T_i` should be replaced by its inverse.

                .. warning::

                    Must be implemented in style "LT".

                EXAMPLES::

                    sage: DoubleAffineHeckeAlgebraSansDuality("B2").LT().signed_generator_product([0,2,1])
                    TX[0,2,1]

                    sage: DoubleAffineHeckeAlgebraSansDuality("B2").LT().signed_generator_product([0,2,1], signs=[1,1,-1])
                    TX[0,2,1] + ((-vl^2+1)/vl)*TX[0,2]

                    sage: DoubleAffineHeckeAlgebraSansDuality("B2", dual_side=True).LT().signed_generator_product([0,2,1])
                    TY[0,2,1]
                """
                return self(self.realization_of().LT().signed_generator_product(word,signs=signs))

            def signed_generator_product_dual(self, word, signs=None):
                r"""
                Similar to `meth`:signed_generator_product` except use the dual classical Hecke algebra.

                .. warning::

                    Must be implemented in basis "LtvLv".

                EXAMPLES::

                    sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2")
                    sage: a = HH.LtvLv().signed_generator_product_dual([2,1]); a
                    Ty[2,1]
                    sage: LT = HH.LT()
                    sage: b = LT.signed_generator_product_dual([2,1]); b
                    TX[2,1]
                    sage: b == LT(a)
                    True
                    sage: HH = DoubleAffineHeckeAlgebraSansDuality("B2", dual_side=True)
                    sage: a = HH.LtvLv().signed_generator_product_dual([2,1]); a
                    Tx[2,1]
                    sage: LT = HH.LT()
                    sage: b = LT.signed_generator_product_dual([2,1]); b
                    TY[2,1]
                    sage: b == LT(a)
                    True
                    sage: LT.signed_generator_product_dual([2,1],[1,-1])
                    TY[2,1] + ((-vl^2+1)/vl)*TY[2]

                """
                return self(self.realization_of().LtvLv().signed_generator_product(word,signs=signs))

            def to_polynomial_module(self, x):
                r"""
                Project `x` into the polynomial module.

                ..warning:: Must be implemented by basis "LtvLv".

                EXAMPLES::

                    sage: HH=DoubleAffineHeckeAlgebraSansDuality("A2",general_linear=True); M = HH.LtvLv()
                    sage: x = M.an_element(); x
                    2*X[(2, 2, 3)] Ty[1,2,1] Y[(2, 2, 3)] + 4*X[(2, 2, 3)] Ty[1,2] Y[(2, 2, 3)] + X[(2, 2, 3)] Y[(2, 2, 3)]
                    sage: M.to_polynomial_module(x)
                    ((2*v^3+4*v^2+1)/v^2)*X[(2, 2, 3)]                    
                    sage: LT = HH.LT()
                    sage: y = LT.a_monomial(); y
                    X[(2, 2, 3)] piX[5] TX[0,1,2]
                    sage: M(y)
                    X[(2, 2, 3)] Ty[2] Y[(2, 2, 1)]
                    sage: LT.to_polynomial_module(y)
                    v^3*X[(2, 2, 3)]

                """
                M = self.realization_of().LtvLv()
                return M.to_polynomial_module(M(x))

            def eta_involution(self, x):
                r"""
                Apply Cherednik's eta involution to `x`.

                It is an automorphism over the rationals that is the identity on length zero extended affine Weyl group
                elements, sends `X^{\lambda}` to `X^{-\lambda}`, `T_i` to `T_i^{-1}` for all `i` including `0`, and
                `q` to `q^{-1}` and `v` to `v^{-1}`.

                ..warning:: Must be implemented by basis "LT".

                EXAMPLES::

                    sage: HH=DoubleAffineHeckeAlgebraSansDuality("A2",general_linear=True); LT = HH.LT()
                    sage: x = LT.a_monomial(); x
                    X[(2, 2, 3)] piX[5] TX[0,1,2]
                    sage: LT.eta_involution(x)
                    ((v^4-2*v^2+1)/v^2)*X[(-2, -2, -3)] piX[5] TX[0] + ((-v^2+1)/v)*X[(-2, -2, -3)] piX[5] TX[0,1] + ((-v^2+1)/v)*X[(-2, -2, -3)] piX[5] TX[0,2] + ((-v^6+3*v^4-3*v^2+1)/v^3)*X[(-2, -2, -3)] piX[5] + ((v^4-2*v^2+1)/v^2)*X[(-2, -2, -3)] piX[5] TX[2] + ((v^4-2*v^2+1)/v^2)*X[(-2, -2, -3)] piX[5] TX[1] + ((-v^2+1)/v)*X[(-2, -2, -3)] piX[5] TX[1,2] + X[(-2, -2, -3)] piX[5] TX[0,1,2]
                    sage: M = HH.LtvLv()
                    sage: y = M.a_monomial(); y
                    X[(2, 2, 3)] Ty[1,2] Y[(2, 2, 3)]
                    sage: M.eta_involution(y)
                    X[(-2, -2, -3)] Ty[1,2] Y[(2, 2, 3)]

                """
                LT = self.realization_of().LT()
                return self(LT.eta_involution(LT(x)))

    class _DAHABases(UniqueRepresentation, BindableClass):
        r"""
        The class of realizations of a double affine Hecke algebra without duality.
        """

        def _repr_(self):
            r"""
            EXAMPLES::

                sage: DoubleAffineHeckeAlgebraSansDuality("A2").LT() # indirect doctest
                LT basis of The double affine Hecke algebra of type ['A', 2, 1]

            """
            return "%s basis of the %s"%(self._prefix, self.realization_of())

        def is_parent_of(self, x):
            return x.parent() == self

    class DoubleAffineHeckeAlgebraSansDualityLT(SmashProductAlgebra, _DAHABases):
        r"""
        DAHA basis LT.

        INPUT:

        - `HH` -- DAHA realization parent

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").LT()
            LT basis of The double affine Hecke algebra of type ['A', 2, 1]

        """

        def __init__(self, HH):
            dat = HH.double_affine_type()
            # To define the product on LT we need the left action of left multiplication by T on KL # T.
            # We start with the left action of the Hecke algebra of the nonextended affine Weyl group of `L`
            # on the group algebra KL.
            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            HH._KLaf = dat.cartan_type().root_system().ambient_space().algebra(HH.base_ring())
            HH._ML = HH._KLaf.nonreduced_demazure_lusztig_operators_on_classical(HH._the_q, dat.q1(), dat.q2(), convention=HH._convention, doubled_parameters=dat.doubled_parameters(),domain=HH.lattice_algebra())
            mcat = ModulesWithBasis(HH.base_ring())
            tmcat = mcat.TensorProducts()
            # Now we define the action of the above nonextended affine Hecke algebra on KL # T
            # This is done by first giving the action of the generators `T_i` for `i \in I`:
            # T_i * (e^\mu # pi T_w) = (T_i(e^\mu) - q1[i] e^{s_i(\mu)}) # pi T_w + e^{s_i(\mu)} # T_i pi T_w
            def Ti_on_LTmod_left((mu,pi,w), i):
                if i == 0:
                    smupoly = HH.s0_on_L(mu)
                else:
                    smupoly = HH._KL.monomial(mu.simple_reflection(i))
                return tensor([HH._ML[i](HH._KL.monomial(mu)) - dat.q1(i) * smupoly, HH._T.monomial((pi,w))],category=tmcat) + tensor([smupoly, HH._T.product_by_generator_on_basis((pi,w), i, side='left')],category=tmcat)
            # Use the HeckeAlgebraRepresentation tool to obtain a left action of the nonextended affine Hecke on KL # T.

            HH._MLTmod = HeckeAlgebraRepresentation(HH._LTmod, Ti_on_LTmod_left, dat.cartan_type(), dat._q1, dat.q2(), HH._the_q, side='left', doubled_parameters=dat.doubled_parameters)

            # The left action of the extended affine Hecke algebra on the module LTmod.
            def left_LTmod_on_basis((ppi,ww),mu,(pi,w)):
                return HH._F_on_LTmod_left_morphism(ppi)(HH._MLTmod.Tw(ww)(HH._LTmod.monomial((mu,pi,w))))

            SmashProductAlgebra.__init__(self, HH._KL, HH._T, left_action=left_LTmod_on_basis, category=Category.join((HH._DAHABasesCategory(),AlgebrasWithBasis(HH.base_ring()).TensorProducts())))
            self._style = "LT"

        def _repr_(self):
            HH = self.realization_of()
            return "%s basis of %s"%(self._style, HH._repr_())

        @cached_method
        def T_generators(self):
            r"""
            A family on the affine Dynkin node set with values given by `T_i`.

            EXAMPLES::

                sage: HH=DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: LT = HH.LT()
                sage: LT.T_generators()
                Finite family {0: TX[0], 1: TX[1], 2: TX[2]}

            """
            return Family(dict([[i, self.factor_embedding(1)(self.factor(1).signed_generator_product([i]))] for i in self.realization_of().double_affine_type().cartan_type().index_set()]))

        @cached_method
        def from_fundamental(self, i):
            r"""
            The image of the fundamental group element of `L` indexed by the special node `i`, in ``self``.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: [(i, HH.LT().from_fundamental(i)) for i in HH._F.special_nodes()]
                [(0, 1), (1, piX[1]), (2, piX[2])]

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2",general_linear=True)
                sage: [(i, HH.LT().from_fundamental(i)) for i in range(5)]
                [(0, 1), (1, piX[1]), (2, piX[2]), (3, piX[3]), (4, piX[4])]

            """
            return self.factor_embedding(1)(self.factor(1).from_fundamental(self.realization_of()._F(i)))

        def from_lattice_algebra(self, a):
            r"""
            The image of the element `a` of `KX` in ``self``.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: LT = HH.LT()
                sage: KL = LT.factor(0)
                sage: a = KL.an_element(); a
                X[(2, 2, 3)]
                sage: b = LT.from_lattice_algebra(a); b
                X[(2, 2, 3)]
                sage: b.parent()
                LT basis of The double affine Hecke algebra of type ['A', 2, 1]

            """
            return self.factor_embedding(0)(a)

        def from_extended_affine_hecke(self, a):
            r"""
            The image of the element `a` of `T` in ``self``.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: LT = HH.LT()
                sage: T = LT.factor(1)
                sage: a = T.an_element(); a
                2*TX[0] + 3*TX[0,1] + 1 + TX[0,1,2] + 4*piX[1] TX[0] + 6*piX[1] TX[0,1] + 2*piX[1] + 2*piX[1] TX[0,1,2] + 8*piX[2] TX[0] + 12*piX[2] TX[0,1] + 4*piX[2] + 4*piX[2] TX[0,1,2]
                sage: b = LT.from_extended_affine_hecke(a); b
                2*TX[0] + 3*TX[0,1] + 1 + TX[0,1,2] + 4*piX[1] TX[0] + 6*piX[1] TX[0,1] + 2*piX[1] + 2*piX[1] TX[0,1,2] + 8*piX[2] TX[0] + 12*piX[2] TX[0,1] + 4*piX[2] + 4*piX[2] TX[0,1,2]
                sage: b.parent()
                LT basis of The double affine Hecke algebra of type ['A', 2, 1]                

            """
            return self.factor_embedding(1)(a)

        def signed_generator_product(self, word, signs=None):
            r"""
            Given a reduced word for an element `w` in the affine Weyl group of `L`, return the image of the basis element `T_w`
            in ``self``.

            The list ``signs``, if present, is a list of entries `+1` or `-1`. The `-1` indicates that
            the corresponding generator `T_i` should be replaced by its inverse.

            EXAMPLES::

                sage: DoubleAffineHeckeAlgebraSansDuality("A2").LT().signed_generator_product([0,2,1])
                TX[0,2,1]

            """
            return self.from_extended_affine_hecke(self.factor(1).signed_generator_product(word,signs=signs))

        def _eta_on_basis(self, (mu, pi, w)):
            r"""
            The eta image of the `LT` basis element indexed by a triple.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: LT = HH.LT()
                sage: mu = HH.lattice().an_element(); mu
                (2, 2, 3)
                sage: pi = HH.F().an_element(); pi
                [2]
                sage: w = HH.T().realization_of().affine_weyl().from_reduced_word([1])
                sage: LT._eta_on_basis((mu,pi,w))
                ((-v^2+1)/v)*X[(-2, -2, -3)] piX[2] + X[(-2, -2, -3)] piX[2] TX[1]

            """
            rw = w.reduced_word()
            signs = [-1] * len(rw)
            return self.from_lattice_algebra(self.factor(0).monomial(-mu)) * self.from_fundamental(pi) * self.signed_generator_product(rw, signs)

        @cached_method
        def _eta_map_on_basis_only(self):
            r"""
            The module map which applies eta to the basis (but not the coefficients).
            """
            return self.module_morphism(on_basis = self._eta_on_basis, codomain=self)

        def eta_involution(self, z):
            r"""
            Cherednik's eta involution.

            It sends `X^\lambda` to `X^{-\lambda}`, `T_i` to its inverse, and `q` and `v_i` to their inverses.

            EXAMPLES::

                sage: LT=DoubleAffineHeckeAlgebraSansDuality("A2").LT()
                sage: z = LT.a_monomial(); z
                X[(2, 2, 3)] piX[2] TX[0,1,2]
                sage: LT.eta_involution(z)
                ((v^4-2*v^2+1)/v^2)*X[(-2, -2, -3)] piX[2] TX[0] + ((-v^2+1)/v)*X[(-2, -2, -3)] piX[2] TX[0,1] + ((-v^2+1)/v)*X[(-2, -2, -3)] piX[2] TX[0,2] + ((-v^6+3*v^4-3*v^2+1)/v^3)*X[(-2, -2, -3)] piX[2] + ((v^4-2*v^2+1)/v^2)*X[(-2, -2, -3)] piX[2] TX[2] + ((v^4-2*v^2+1)/v^2)*X[(-2, -2, -3)] piX[2] TX[1] + ((-v^2+1)/v)*X[(-2, -2, -3)] piX[2] TX[1,2] + X[(-2, -2, -3)] piX[2] TX[0,1,2]                
                sage: LT.eta_involution(z*LT.realization_of().double_affine_type().q())
                ((v^4-2*v^2+1)/(q*v^2))*X[(-2, -2, -3)] piX[2] TX[0] + ((-v^2+1)/(q*v))*X[(-2, -2, -3)] piX[2] TX[0,1] + ((-v^2+1)/(q*v))*X[(-2, -2, -3)] piX[2] TX[0,2] + ((-v^6+3*v^4-3*v^2+1)/(q*v^3))*X[(-2, -2, -3)] piX[2] + ((v^4-2*v^2+1)/(q*v^2))*X[(-2, -2, -3)] piX[2] TX[2] + ((v^4-2*v^2+1)/(q*v^2))*X[(-2, -2, -3)] piX[2] TX[1] + ((-v^2+1)/(q*v))*X[(-2, -2, -3)] piX[2] TX[1,2] + 1/q*X[(-2, -2, -3)] piX[2] TX[0,1,2]                

            """
            z = z.map_coefficients(self.realization_of().double_affine_type().eta_base_ring_map()) # apply eta to the coefficients
            return self._eta_map_on_basis_only()(z) # then to the basis elements

    class DoubleAffineHeckeAlgebraSansDualityTL(SmashProductAlgebra, _DAHABases):
        r"""
        DAHA basis TL.

        INPUT:

        - `HH` -- DAHA realization parent

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").TL()
            TL basis of The double affine Hecke algebra of type ['A', 2, 1]

        """

        def __init__(self, HH):
            dat = HH.double_affine_type()
            # Define the left action of the Hecke algebra of the affine Weyl group of `L` on the group algebra KL.
            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            HH._LM = HH._KLaf.nonreduced_demazure_lusztig_operators_on_classical(HH._the_q, dat.q1(), dat.q2(), convention=HH._convention,side='right', doubled_parameters=dat.doubled_parameters(),domain=HH._KL)
            mcat = ModulesWithBasis(HH.base_ring())
            tmcat = mcat.TensorProducts()
            # The above nonextended affine Hecke algebra has a right action on T # KL
            # via the action of the generators `T_i` for `i \in I`:
            # (pi T_w # e^\mu) # T_i = pi T_w # (T_i(e^\mu) - q1[i] e^{s_i(\mu)}) + pi T_w T_i # e^{s_i(\mu)}
            def Ti_on_TLmod_right((pi,w,mu), i):
                if i == 0:
                    smupoly = HH.s0_on_L(mu)
                else:
                    smupoly = HH._KL.monomial(mu.simple_reflection(i))
                return tensor([HH._T.monomial((pi,w)), HH._LM[i](HH._KL.monomial(mu)) - dat.q1(i) * smupoly],category=tmcat) + tensor([HH._T.product_by_generator_on_basis((pi,w), i, side='right'),smupoly],category=tmcat)

            # Use the HeckeAlgebraRepresentation tool to obtain a right action of the nonextended affine Hecke on T # KL.
            HH._TLmodM = HeckeAlgebraRepresentation(HH._TLmod, Ti_on_TLmod_right, dat.cartan_type(), dat.q1(), dat.q2(), HH._the_q, side='right', doubled_parameters=dat.doubled_parameters)

            # Finally we obtain the right action of the extended affine Hecke algebra on the module TLmod.
            def right_TLmod_on_basis((ppi,ww),(pi,w),mu):
                return HH._TLmodM.Tw(ww)(HH._F_on_TLmod_right_morphism(ppi)(HH._TLmod.monomial((pi,w,mu))))

            SmashProductAlgebra.__init__(self, HH._T, HH._KL, right_action=right_TLmod_on_basis, category=Category.join((HH._DAHABasesCategory(),AlgebrasWithBasis(HH.base_ring()).TensorProducts())))

            self._style = "TL"

        def _repr_(self):
            HH = self.realization_of()
            return "%s basis of %s"%(self._style, HH._repr_())

        @cached_method
        def T_generators(self):
            r"""
            A family on the affine Dynkin node set with values given by `T_i`.

            EXAMPLES::

                sage: HH=DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: TL = HH.TL()
                sage: TL.T_generators()
                Finite family {0: TX[0], 1: TX[1], 2: TX[2]}

            """
            return Family(dict([[i, self.factor_embedding(0)(self.factor(0).signed_generator_product([i]))] for i in self.realization_of().double_affine_type().cartan_type().index_set()]))

        @cached_method
        def from_fundamental(self, i):
            r"""
            The image of the fundamental group element of `L` indexed by the special node `i`, in ``self``.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: [(i, HH.TL().from_fundamental(i)) for i in HH._F.special_nodes()]
                [(0, 1), (1, piX[1]), (2, piX[2])]
            """
            return self.factor_embedding(0)(self.factor(0).from_fundamental(self.realization_of()._F(i)))

        def from_lattice_algebra(self, a):
            r"""
            The image of the element `a` of `KL` in ``self``.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: TL = HH.TL()
                sage: KL = TL.factor(1)
                sage: a = KL.an_element(); a
                X[(2, 2, 3)]
                sage: b = TL.from_lattice_algebra(a); b
                X[(2, 2, 3)]
                sage: b.parent()
                TL basis of The double affine Hecke algebra of type ['A', 2, 1]

            """
            return self.factor_embedding(1)(a)

        def from_extended_affine_hecke(self, a):
            r"""
            The image of the element `a` of `T` in ``self``.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: TL = HH.TL()
                sage: T = TL.factor(0)
                sage: a = T.an_element(); a
                2*TX[0] + 3*TX[0,1] + 1 + TX[0,1,2] + 4*piX[1] TX[0] + 6*piX[1] TX[0,1] + 2*piX[1] + 2*piX[1] TX[0,1,2] + 8*piX[2] TX[0] + 12*piX[2] TX[0,1] + 4*piX[2] + 4*piX[2] TX[0,1,2]
                sage: b = TL.from_extended_affine_hecke(a); b
                2*TX[0] + 3*TX[0,1] + 1 + TX[0,1,2] + 4*piX[1] TX[0] + 6*piX[1] TX[0,1] + 2*piX[1] + 2*piX[1] TX[0,1,2] + 8*piX[2] TX[0] + 12*piX[2] TX[0,1] + 4*piX[2] + 4*piX[2] TX[0,1,2]
                sage: b.parent()
                TL basis of The double affine Hecke algebra of type ['A', 2, 1]                

            """
            return self.factor_embedding(0)(a)

        def signed_generator_product(self, word, signs=None):
            r"""
            Given a reduced word for an element `w` in the affine Weyl group of `L`, return the image of the basis element `T_w`
            in ``self``.

            EXAMPLES::

                sage: DoubleAffineHeckeAlgebraSansDuality("A2").TL().signed_generator_product([0,2,1],[1,-1,1])
                ((-v^2+1)/v)*TX[0,1] + TX[0,2,1]

            """
            return self.from_extended_affine_hecke(self.factor(0).signed_generator_product(word,signs=signs))


    class DoubleAffineHeckeAlgebraSansDualityLtvLv(SmashProductAlgebra, _DAHABases):
        r"""
        DAHA basis LtvLv.

        INPUT:

        - `HH` -- DAHA realization parent

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").LtvLv()
            LtvLv basis of The double affine Hecke algebra of type ['A', 2, 1]


        """

        def __init__(self, HH):
            dat = HH.double_affine_type()
            # To define the product on LtvLv we just coerce to LT and back.
            # That is, we define the twist for LtvLv in terms of the twist for LT.
            mcat = ModulesWithBasis(HH.base_ring())
            tmcat = mcat.TensorProducts()
            KL = HH._KL
            coerceoid = tensor([HH.T().coerce_map_from(HH.tv_Lv()), KL._identity_map()], category=tmcat)
            idocoerce = tensor([KL._identity_map(), HH.tv_Lv().coerce_map_from(HH.T())], category=tmcat)

            HH._twist_LtvLv = SetMorphism(Hom(HH._tvLvLmod,HH._LtvLvmod,mcat),idocoerce * HH.LT().twist() * coerceoid)

            SmashProductAlgebra.__init__(self, HH._KL, HH.tv_Lv(), twist_morphism=HH._twist_LtvLv, category=Category.join((HH._DAHABasesCategory(),AlgebrasWithBasis(HH.base_ring()).TensorProducts())))
            self._style = "LtvLv"

            self._to_polynomial_morphism = self.module_morphism(on_basis=self._to_polynomial_module_on_basis, category=ModulesWithBasis(self.base_ring()), codomain = HH._KL)

        def _repr_(self):
            HH = self.realization_of()
            return "%s basis of %s"%(self._style, HH._repr_())

        def from_lattice_algebra(self, a):
            r"""
            The image of the element `a` of `KL` in ``self``.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: LtvLv = HH.LtvLv()
                sage: KL = HH.lattice_algebra()
                sage: a = KL.an_element(); a
                X[(2, 2, 3)]
                sage: b = LtvLv.from_lattice_algebra(a); b
                X[(2, 2, 3)]
                sage: b.parent()
                LtvLv basis of The double affine Hecke algebra of type ['A', 2, 1]
            """
            return self.factor_embedding(0)(a)


        def from_dual_lattice_algebra(self, a):
            r"""
            The embedding of the group algebra of `Lv` into ``self``.

            EXAMPLES::

                sage: HH=DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: KLv = HH.dual_lattice_algebra()
                sage: a = KLv.an_element(); a
                Y[(2, 2, 3)]
                sage: a.parent()
                Group algebra of the Ambient space of the Root system of type ['A', 2] over Fraction Field of Multivariate Polynomial Ring in q, v over Rational Field
                sage: LtvLv = HH.LtvLv()
                sage: b = LtvLv.from_dual_lattice_algebra(a); b
                Y[(2, 2, 3)]
                sage: b.parent()
                LtvLv basis of The double affine Hecke algebra of type ['A', 2, 1]

            """
            return self.factor_embedding(1)(self.factor(1).factor_embedding(1)(a))

        def from_dual_classical_hecke(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra "tv" of dual type, to ``self``.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: tv_Lv = HH.tv_Lv()
                sage: tv = HH.tv()
                sage: a = tv.an_element(); a
                2*Ty[1,2,1] + 4*Ty[1,2] + 1
                sage: LtvLv = HH.LtvLv()
                sage: b = LtvLv.from_dual_classical_hecke(a); b
                2*Ty[1,2,1] + 4*Ty[1,2] + 1
                sage: b == LtvLv(a)
                True
            """
            return self.factor_embedding(1)(self.factor(1).factor_embedding(0)(a))

        def signed_generator_product_dual(self, word, signs=None):
            r"""
            Given a reduced word for an element `w` in the classical Weyl group of `Lv`, return the image of the basis element `T_w`
            in ``self``.

            EXAMPLES::

                sage: DoubleAffineHeckeAlgebraSansDuality("A2").LtvLv().signed_generator_product_dual([2,1],[1,-1])
                Ty[2,1] + ((-v^2+1)/v)*Ty[2]

            """
            return self.factor_embedding(1)(self.factor(1).signed_generator_product(word, signs=signs))

        def _to_polynomial_module_on_basis(self, (mu, u, nu)):
            r"""
            Given a basis triple for ``self``, return the projected element in the polynomial module.

                sage: HH=DoubleAffineHeckeAlgebraSansDuality("A3",general_linear=True); M = HH.LtvLv()
                sage: mu = HH.lattice().an_element(); mu
                (2, 2, 3, 0)
                sage: u = HH.double_affine_type().extended_affine_weyl().dual_classical_weyl().an_element(); u
                s1*s2*s3
                sage: nu = HH.dual_lattice().an_element(); nu
                (2, 2, 3, 0)
                sage: HH.LtvLv()._to_polynomial_module_on_basis((mu,u,nu))
                v^8*X[(2, 2, 3, 0)]

            """
            HH = self.realization_of()
            return HH._KL.term(mu, HH.T_eigenvalue(u) * HH.Y_eigenvalue(nu))

        def to_polynomial_module(self, x):
            r"""
            Project an element into the polynomial module.

            EXAMPLES::

                sage: HH=DoubleAffineHeckeAlgebraSansDuality("A2",general_linear=True); M = HH.LtvLv()
                sage: x = M.an_element(); x
                2*X[(2, 2, 3)] Ty[1,2,1] Y[(2, 2, 3)] + 4*X[(2, 2, 3)] Ty[1,2] Y[(2, 2, 3)] + X[(2, 2, 3)] Y[(2, 2, 3)]
                sage: M.to_polynomial_module(x)
                ((2*v^3+4*v^2+1)/v^2)*X[(2, 2, 3)]
            """
            return self._to_polynomial_morphism(x)


    class DoubleAffineHeckeAlgebraSansDualityLLvtv(SmashProductAlgebra, _DAHABases):
        r"""
        DAHA basis LLvtv.

        INPUT:

        - `HH` -- DAHA realization parent

        EXAMPLES::

            sage: DoubleAffineHeckeAlgebraSansDuality("A2").LLvtv()
            LLvtv basis of The double affine Hecke algebra of type ['A', 2, 1]


        """

        def __init__(self, HH):
            dat = HH.double_affine_type()
            # To define the product on LLvtv we just coerce to LT and back.
            # That is, we define the twist for LLvtv in terms of the twist for LT.
            mcat = ModulesWithBasis(HH.base_ring())
            tmcat = mcat.TensorProducts()
            KL = HH._KL
            coerceoid = tensor([HH.T().coerce_map_from(HH.Lv_tv()), KL._identity_map()], category=tmcat)
            idocoerce = tensor([KL._identity_map(), HH.Lv_tv().coerce_map_from(HH.T())], category=tmcat)

            HH._twist_LLvtv = SetMorphism(Hom(HH._LvtvLmod,HH._LLvtvmod,mcat),idocoerce * HH.LT().twist() * coerceoid)

            SmashProductAlgebra.__init__(self, HH._KL, HH.Lv_tv(), twist_morphism=HH._twist_LLvtv, category=Category.join((HH._DAHABasesCategory(),AlgebrasWithBasis(HH.base_ring()).TensorProducts())))
            self._style = "LLvtv"

            self._to_polynomial_morphism = self.module_morphism(on_basis=self._to_polynomial_module_on_basis, category=ModulesWithBasis(self.base_ring()), codomain = HH._KL)

        def _repr_(self):
            HH = self.realization_of()
            return "%s basis of %s"%(self._style, HH._repr_())

        def from_lattice_algebra(self, a):
            r"""
            The image of the element `a` of `KL` in ``self``.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: LLvtv = HH.LLvtv()
                sage: KL = HH.lattice_algebra()
                sage: a = KL.an_element(); a
                X[(2, 2, 3)]
                sage: b = LLvtv.from_lattice_algebra(a); b
                X[(2, 2, 3)]
                sage: b.parent()
                LLvtv basis of The double affine Hecke algebra of type ['A', 2, 1]
            """
            return self.factor_embedding(0)(a)


        def from_dual_lattice_algebra(self, a):
            r"""
            The embedding of the group algebra of `Lv` into ``self``.

            EXAMPLES::

                sage: HH=DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: KLv = HH.dual_lattice_algebra()
                sage: a = KLv.an_element(); a
                Y[(2, 2, 3)]
                sage: a.parent()
                Group algebra of the Ambient space of the Root system of type ['A', 2] over Fraction Field of Multivariate Polynomial Ring in q, v over Rational Field
                sage: LLvtv = HH.LLvtv()
                sage: b = LLvtv.from_dual_lattice_algebra(a); b
                Y[(2, 2, 3)]
                sage: b.parent()
                LLvtv basis of The double affine Hecke algebra of type ['A', 2, 1]

            """
            return self.factor_embedding(1)(self.factor(1).factor_embedding(0)(a))

        def from_dual_classical_hecke(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra "tv" of dual type, to ``self``.

            EXAMPLES::

                sage: HH = DoubleAffineHeckeAlgebraSansDuality("A2")
                sage: Lv_tv = HH.Lv_tv()
                sage: tv = HH.tv()
                sage: a = tv.an_element(); a
                2*Ty[1,2,1] + 4*Ty[1,2] + 1
                sage: LLvtv = HH.LLvtv()
                sage: b = LLvtv.from_dual_classical_hecke(a); b
                2*Ty[1,2,1] + 4*Ty[1,2] + 1
                sage: b == LLvtv(a)
                True
            """
            return self.factor_embedding(1)(self.factor(1).factor_embedding(1)(a))

        def signed_generator_product_dual(self, word, signs=None):
            r"""
            Given a reduced word for an element `w` in the classical Weyl group of `Lv`, return the image of the basis element `T_w`
            in ``self``.

            EXAMPLES::

                sage: DoubleAffineHeckeAlgebraSansDuality("A2").LLvtv().signed_generator_product_dual([2,1],[1,-1])
                Ty[2,1] + ((-v^2+1)/v)*Ty[2]

            """
            return self.factor_embedding(1)(self.factor(1).signed_generator_product(word, signs=signs))

        def _to_polynomial_module_on_basis(self, (mu, nu, u)):
            r"""
            Given a basis triple for ``self``, return the projected element in the polynomial module.

                sage: HH=DoubleAffineHeckeAlgebraSansDuality("A3",general_linear=True); M = HH.LLvtv()
                sage: mu = HH.lattice().an_element(); mu
                (2, 2, 3, 0)
                sage: u = HH.double_affine_type().extended_affine_weyl().dual_classical_weyl().an_element(); u
                s1*s2*s3
                sage: nu = HH.dual_lattice().an_element(); nu
                (2, 2, 3, 0)
                sage: HH.LLvtv()._to_polynomial_module_on_basis((mu,nu,u))
                v^8*X[(2, 2, 3, 0)]

            """
            HH = self.realization_of()
            return HH._KL.term(mu, HH.T_eigenvalue(u) * HH.Y_eigenvalue(nu))

        def to_polynomial_module(self, x):
            r"""
            Project an element into the polynomial module.

            EXAMPLES::

                sage: HH=DoubleAffineHeckeAlgebraSansDuality("A2",general_linear=True); M = HH.LLvtv()
                sage: x = M.an_element(); x
                2*X[(2, 2, 3)] Y[(2, 2, 3)] Ty[1,2,1] + 4*X[(2, 2, 3)] Y[(2, 2, 3)] Ty[1,2] + X[(2, 2, 3)] Y[(2, 2, 3)]
                sage: M.to_polynomial_module(x)
                ((2*v^3+4*v^2+1)/v^2)*X[(2, 2, 3)]

            """
            return self._to_polynomial_morphism(x)

class RamYipFormula(SageObject):
    r"""
    Container for data governing the Ram-Yip formula for (non)symmetric Macdonald polynomials.

    INPUTS::

    - ``daha`` -- A :class:`DoubleAffineHeckeAlgebraSansDuality` object.
    - ``weight`` -- a weight in the ambient `X` lattice of ``daha``. It can be a weight or tuple.

    ``weight`` is used to make a translation element `t`. The Ram-Yip formula requires a reduced word
    for the minimum length coset representative `m` of this translation with respect to the finite
    Weyl group.

    Optional keyword arguments::

    - ``symmetric`` -- True or False. If True produce the alcove paths for the symmetric
    Macdonald-Koornwinder polynomial. More precisely, starting with the Ram-Yip formula for the
    nonsymmetric MK polynomial of lowest weight ``weight``, apply the Hecke symmetrizer over the finite Weyl group.
    - ``reduced_word`` -- Use the given reduced word (in the form of a tuple) for the above element `m`.
    - ``start`` -- An element of the extended affine Weyl group (or a tuple or list for a reduced word)
    This is ignored if ``symmetric`` is True.
    - ``verbose`` -- True or False (default: False) If True print extra information
    - ``coset_reps`` -- True or False (default: True) This is ignored unless ``symmetric`` is True.
    In that case, if ``coset_reps`` is True then perform the Hecke symmetrization over the coset
    representatives for the stabilizer of the weight. If False then use the entire finite Weyl group.

    EXAMPLES::

        sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',2],untwisted=False,reduced=False,dual_reduced=False)
        sage: from sage.algebras.double_affine_hecke_algebra import RamYipFormula
        sage: RamYipFormula(HH, (1,0), start=[2,1])
        Ram Yip Formula for nonsymmetric Macdonald polynomial of weight (1, 0) and Double Affine Type ['C', 2, 1]^* nonreduced nondual-reduced
        sage: RamYipFormula(HH, (0,-1))
        Ram Yip Formula for nonsymmetric Macdonald polynomial of weight (0, -1) and Double Affine Type ['C', 2, 1]^* nonreduced nondual-reduced
        sage: RamYipFormula(HH, (1,0), symmetric=True)
        Ram Yip Formula for symmetric Macdonald polynomial of weight (1, 0) and Double Affine Type ['C', 2, 1]^* nonreduced nondual-reduced
    """

    def __init__(self, daha, weight, **keywords):
        verbose = keywords.get('verbose',False)
        self._daha = daha
        if isinstance(weight, tuple):
            # weight is a tuple
            self._weight_ambient = daha.lattice()(weight)
        elif weight in self._daha.lattice():
            # weight is in the ambient space
            self._weight_ambient = weight
        elif weight in self._daha.double_affine_type().cartan_type().classical().root_system().weight_lattice():
            # weight is in the weight lattice
            self._weight_ambient = weight.to_ambient()
        else:
            raise TypeError("{} should be a tuple or weight".format(weight))
        self._weight = self._weight_ambient.to_weight_space(ZZ)

        # get the extended affine Weyl group of type Y
        dat = daha.double_affine_type()
        ddat = daha.dual_double_affine_type()
        self._EY = ddat.extended_affine_weyl()
        # a realization of the extended affine Weyl group
        self._WY = self._EY.PvW0()
        # affine Weyl group
        self._Wa = self._EY.affine_weyl()
        # finite Weyl group
        self._W = self._EY.dual_classical_weyl()

        # the element that sends ``weight`` to antidominant
        # right multiplication by this element fixes the directions of alcove paths
        self._u_weight_reduced_word = self._weight_ambient.reduced_word(positive=False)
        self._u_weight = self.classical_weyl_Y().from_reduced_word(self._u_weight_reduced_word)

        # compute a reduced word for the min coset rep of t_{weight}
        if ddat.general_linear():
            translation = self.extended_affine_weyl_Y().from_dual_translation(self._weight_ambient)
        else:
            translation = self.extended_affine_weyl_Y().from_dual_translation(self._weight)
        self._min_rep = translation.to_affine_grassmannian()
        # the projections of the min coset rep of the translation by ``weight``
        # into the fundamental group and the affine Weyl group on the right
        self._fund_trans = self._min_rep.to_fundamental_group()
        self._aff_trans = self._min_rep.to_affine_weyl_right()

        # get the preferred reduced word for the min rep
        reduced_word = keywords.get('reduced_word',None)
        if reduced_word:
            if isinstance(reduced_word, (tuple, list)):
                if verbose:
                    print "reduced word: "
                    print reduced_word
                    print "projection of translation element to affine Weyl:"
                    print self._aff_trans.reduced_word()
                self._the_reduced_word = [i for i in reduced_word]
            else:
                raise TypeError("reduced word should be a tuple or list")
        else:
            self._the_reduced_word = [i for i in self._aff_trans.reduced_word()]
        self._N = len(self._the_reduced_word)

        # check for the symmetric case
        self._symmetric = keywords.get('symmetric',False)

        # set the extended affine Weyl group element ``start``
        if not self._symmetric:
            start = keywords.get('start',None)
            if start:
                if isinstance(start, (tuple,list)):
                    self._start = self.extended_affine_weyl_Y().from_reduced_word(start)
                elif start in self.extended_affine_weyl_Y():
                    self._start = start
                else:
                    raise TypeError("{} must be an extended affine Weyl group element or reduced word".format(start))
            else:
                self._start = self.extended_affine_weyl_Y().one()
        else:
            self._coset_reps = keywords.get('coset_reps', True)

        # compute the list of inversions of the reduced word and associated data
        self._reflections = [self.affine_weyl_Y().from_reduced_word([x for x in reversed(self._the_reduced_word[i:])]+[self._the_reduced_word[i]]+self._the_reduced_word[i:]) for i in range(self._N)] # these are reflections in the affine Weyl group
        self._inversions = [x.reflection_to_root() for x in self._reflections] # the associated affine roots
        # get the multiples of delta^Y in the affine roots
        self._inversion_depths = [x[0] for x in self._inversions]
        # negated classical projections of inversions; these are all classical positive roots expressed in the ambient space
        self._inversions_classical = [-x.to_classical().to_ambient() for x in self._inversions]

        if verbose:
            print ""
            print "%s reflections"%(self._N)
            print "affine Grassmannian element: %s %s"%(self._fund_trans, self._aff_trans)
            print "inversions"
            print ""
            for k in range(self._N):
                print "%s: %s*delta-%s"%(k+1, self._inversion_depths[k], self._inversions_classical[k])

        zerovector = self._daha.lattice().zero()

        if ddat.extended_base_ring():
            self._diff_recip_v = [-ddat.c_parameter(ddat.v(i)) for i in self._the_reduced_word]
            self._diff_recip_v2 = [-ddat.c_parameter(ddat.v2(i)) for i in self._the_reduced_word]
        else:
            self._diff_recip_v = [ddat.diff_reciprocal(ddat.v(i)) for i in self._the_reduced_word]
            self._diff_recip_v2 = [ddat.diff_reciprocal(ddat.v2(i)) for i in self._the_reduced_word]

        self._eigen = [self._daha.Y_eigenvalue(self._inversions_classical[k],zerovector,-self._inversion_depths[k]) for k in range(self._N)]
        self._pos_fold_factor = [-(self._diff_recip_v[k]+self._diff_recip_v2[k]*self._eigen[k])/(1-self._eigen[k]**2) for k in range(self._N)]
        self._neg_fold_factor = [-(self._diff_recip_v[k]*self._eigen[k]**2+self._diff_recip_v2[k]*self._eigen[k])/(1-self._eigen[k]**2) for k in range(self._N)]

        if verbose:
            for k in range(self._N):
                print "reflection {}\n    vfactor: {}\n   v2factor: {}\n    eigenvalue: {}\n   posfactor: {}\n    negfactor: {}".format(k+1, self._diff_recip_v[k], self._diff_recip_v2[k], self._eigen[k], self._pos_fold_factor[k], self._neg_fold_factor[k])

    @cached_method
    def alcove_paths(self):
        if not self._symmetric:
            return self.generate_ram_yip_terms(self._start)
        answer = []
        if not self._coset_reps:
            for u in self._W:
                answer = answer + self.generate_ram_yip_terms(u.reduced_word())
        else:
            orbit = self._weight.orbit()
            for wt in orbit:
                u_red = wt.reduced_word()
                answer = answer + self.generate_ram_yip_terms(u_red)
        return answer

    def _repr_(self):
        r"""
        A string representing ``self``.

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',2],untwisted=False,reduced=False,dual_reduced=False,extra_parameters=True)
            sage: from sage.algebras.double_affine_hecke_algebra import RamYipFormula
            sage: RY = RamYipFormula(HH, (1,0)); RY
            Ram Yip Formula for nonsymmetric Macdonald polynomial of weight (1, 0) and Double Affine Type ['C', 2, 1]^* nonreduced nondual-reduced with extended base ring
        """
        return "Ram Yip Formula for {}symmetric Macdonald polynomial of weight {} and {}".format("" if self._symmetric else "non", self._weight_ambient, self.DAHA().double_affine_type())

    def DAHA(self):
        r"""
        The double affine Hecke algebra of ``self``.

        EXAMPLES::

            sage: HH=DoubleAffineHeckeAlgebraSansDuality(['A',1])
            sage: from sage.algebras.double_affine_hecke_algebra import RamYipFormula
            sage: RY = RamYipFormula(HH, (0,1))
            sage: RY.DAHA()
            The double affine Hecke algebra of type ['A', 1, 1]
        """
        return self._daha

    def extended_affine_weyl_Y(self):
        return self._WY

    def affine_weyl_Y(self):
        return self._Wa

    def classical_weyl_Y(self):
        return self._W

    def generate_ram_yip_terms(self, u=None):
        r"""
        Generate the terms of the Ram-Yip formula.

        More precisely it is for for `T_u E_\lambda` where `\lambda` is the weight
        used in the definition of ``self``.

        INPUTS::

        - `u` -- (default: None) an element of the extended affine Weyl group of type Y
        It may also be a reduced word.

        The output is a list of Ram-Yip terms (alcove paths) for `X^u` and the intertwiner
            for the given translation.

        EXAMPLES::

            sage: HH=DoubleAffineHeckeAlgebraSansDuality(['A',2])
            sage: from sage.algebras.double_affine_hecke_algebra import RamYipFormula
            sage: RY = RamYipFormula(HH, (0,1,0))
            sage: RY.alcove_paths()
            [Ram-Yip term: start: 1 end: t[-Lambda[1] + Lambda[2]] * s2 folds: [0] signs: [1] value: v^-1, Ram-Yip term: start: 1 end: t[Lambda[1]] * s1*s2 folds: [1] signs: [1] value: v^-2 * (v - 1) * (v + 1) * (q^3*v^4 - 1)^-1]

            sage: HH=DoubleAffineHeckeAlgebraSansDuality(['B',2],untwisted=False,reduced=False,dual_reduced=False,extra_parameters=True)
            sage: RY = RamYipFormula(HH, (1,0))
            sage: RY.alcove_paths()
            [Ram-Yip term: start: 1 end: t[Lambda[1]] * s1*s2*s1 folds: [0] signs: [-1] value: vl^-2 * v^-1, Ram-Yip term: start: 1 end: 1 folds: [1] signs: [-1] value: (-1) * v0 * q * (q*v*vl^2*v0 - 1)^-1 * (q*v*vl^2*v0 + 1)^-1 * (q*v*vl^2*v0*c2 + cz)]
        """
        def flat_list(lis):
            new_list = []
            for x in lis:
                new_list = new_list + x
            return new_list
        
        WY = self.extended_affine_weyl_Y()
        if u is None:
            u = WY.one()
        elif isinstance(u, (tuple, list)):
            u = WY.from_reduced_word(u)
        elif u not in WY:
            raise TypeError("Must be an element of the extended affine Weyl group of type Y")
        lis = [RamYipTerm(self, u * self._min_rep, None, None, None, u)]
        # generate the terms
        for k in range(self._N):
            lis = flat_list([x.apply_intertwiner(self._inversions[k], self._reflections[k], self._pos_fold_factor[k], self._neg_fold_factor[k]) for x in lis])
        # normalize 
        if self._symmetric:
            norm_const = self._daha.T_eigenvalue(self._u_weight.inverse())/self._daha.T_eigenvalue(u)
        else:
            norm_const = self._daha.T_eigenvalue(u.to_dual_classical_weyl()*self._u_weight.inverse())
        return [x.normalize_term(norm_const) for x in lis]

    def antidominant_normalization_factor(self, mu):
        r"""
        A factor in the normalization constant for a symmetric Macdonald-Koornwinder polynomial.

        INPUTS::

        - `\mu` -- A dominant weight in the ambient `X` lattice.

        The output is the constant `c(\mu)` which is defined as follows.
        The symmetric MK polynomial `P_\mu` can be computed from the nonsymmetric MK polynomial `E_\mu` by::

            P_\mu = \sum_u v_u T_u E_\mu

        where `u` runs over the representatives of the finite Weyl group `W` with respect to the stabilizer of `\mu`.
        One may also get `P_\mu` by the formula::

            c(\mu) P_\mu = \sum_\mu v_u T_u E_{w_0(\mu)}

        with the sum over the same `u`. Let `(i_1,i_2,\dotsc,i_\ell)` be a reduced word for `u_\mu`, the
        shortest Weyl element taking the dominant weight `\mu` to its antidominant orbit element `w_0(\mu)`.
        Then we have::

            c(\mu) = c(\alpha_{i_1}, \mu) c(\alpha_{i_2}, s_{i_1}(\mu) c(\alpha_{i_3}, s_{i_1}s_{i_2}(\alpha_{i_3})) \dotsm

        where::

            c(\beta,\lambda) = v_\beta^2 - v_\beta \dfrac{(v_\beta-v_\beta^{-1}) + (v_{2\beta}-v_{2\beta}^{-1})\chi_{-\beta,\lambda}}{1 - \chi_{-\beta,\lambda}^2}

        EXAMPLES::

            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',2],untwisted=False,reduced=False,dual_reduced=False)
            sage: from sage.algebras.double_affine_hecke_algebra import RamYipFormula
            sage: RY = RamYipFormula(HH, (1,0), symmetric=True, no_terms=True)
            sage: RY.antidominant_normalization_factor(HH.lattice()((1,0))).factor()
            (q*vl - 1)^-1 * (q*vl + 1)^-1 * (q*vl^2 - 1) * (q*vl^2 + 1) * (q*v*vl^2*v0 - 1) * (q*v*vl^2*v0 + 1) * (q^2*v*vl^2*v0 - 1)^-1 * (q^2*v*vl^2*v0 + 1)^-1 * (q^2*v^2*vl^2 + 1)
            sage: RY.antidominant_normalization_factor(HH.lattice()((1,1))).factor()
            (q*v*v0 - 1) * (q*v*v0 + 1) * (q^2*v*v0 - 1)^-1 * (q^2*v*v0 + 1)^-1 * (q*v*vl*v0 - 1) * (q*v*vl*v0 + 1) * (q^2*v^2 + 1) * (q^2*v*vl*v0 - 1)^-1 * (q^2*v*vl*v0 + 1)^-1 * (q^2*v^2*vl^2 + 1)
            sage: HH = DoubleAffineHeckeAlgebraSansDuality(['B',3],untwisted=False,reduced=False,dual_reduced=False)
            sage: RY = RamYipFormula(HH, (1,0,0), symmetric=True, no_terms=True)
            sage: RY.antidominant_normalization_factor(HH.lattice()((1,1,1))).factor()
            (q*v*v0 - 1) * (q*v*v0 + 1) * (q^2*v*v0 - 1)^-1 * (q^2*v*v0 + 1)^-1 * (q*v*vl*v0 - 1) * (q*v*vl*v0 + 1) * (q^2*v^2 + 1) * (q^2*v*vl*v0 - 1)^-1 * (q^2*v*vl*v0 + 1)^-1 * (q*v*vl^2*v0 - 1) * (q*v*vl^2*v0 + 1) * (q^2*v*vl^2*v0 - 1)^-1 * (q^2*v*vl^2*v0 + 1)^-1 * (q^2*v^2*vl^2 + 1) * (q^2*v^2*vl^4 + 1)
        """
        reduced_word = mu.reduced_word(positive=False)
        HH = self.DAHA()
        ddat = HH.dual_double_affine_type()
        X = HH.lattice()
        K = HH.base_ring()
        c = K.one()
        for i in reversed(reduced_word):
            vi = ddat.v(i)
            v2i = ddat.v2(i)
            chi = HH.Y_eigenvalue(-X.simple_root(i), mu)
            factor = (vi**2 - vi*((vi-1/vi) +(v2i-1/v2i)*chi  )/(K.one()-chi**2)) 
            c = c * factor
            mu = mu.simple_reflection(i)
        return c

class RamYipTerm(SageObject):
    r"""
    Data structure for a term of the Ram-Yip formula.

    INPUTS::

        - ``RY`` -- a `class`:RamYipFormula object
        - ``end`` -- Element of the extended affine Weyl group `W_a(\tilde{Y})`
        - ``folds`` -- binary list (default: None)
        - ``signs`` -- list of plus and minus ones (default: None)
        - ``value`` -- coefficient from base ring (default: None)
        - ``start`` -- Element of the finite Weyl group `W(Y)` (default: None)

        If ``folds`` or ``signs`` is None then the empty list is used.
        If ``value`` is None the unit in the base ring is used.
        If ``start`` is None then the identity element is used.

    """
    def __init__(self, RY, end, folds=None, signs=None, value=None, start=None):
        self._RY = RY
        self._end = end
        self._folds = folds if folds is not None else []
        self._signs = signs if signs is not None else []
        self._value = value if value is not None else RY._daha.base_ring().one()
        self._start = start if start is not None else RY._daha.extended_affine_weyl().classical_weyl().one()

    def _repr_(self):
        return "Ram-Yip term:\n start: %s end: %s \n folds: %s\n signs: %s\n value: %s"%(self._start, self._end, self._folds, self._signs, self._value.factor())

    def evaluate(self, coefficient_only=None):
        if coefficient_only:
            return self._value * (self._RY._daha.Xu_on_poly_gen(self._end).coefficients()[0])
        return self._value * self._RY._daha.Xu_on_poly_gen(self._end)

    def build_term(self, fold, sign, refl=None, eigenfactor=None):
        return RamYipTerm(self._RY, self._end if refl is None else self._end * refl, self._folds+[fold], self._signs+[sign], self._value if eigenfactor is None else self._value * eigenfactor, self._start)

    def apply_intertwiner(self, betak, reflk, eigenp, eigenn):
        sign = 1 if self._end.action_on_affine_roots(-betak).to_classical().is_positive_root() else -1
        eigen = eigenp if sign > 0 else eigenn
        return [self.build_term(0, sign), self.build_term(1, sign, reflk, eigen)]

    def normalize_term(self, coef):
        self._value = self._value/coef
        return self

    def weight(self):
        r"""
        The weight of an alcove path.
        """
        return self._end.to_dual_translation_left()

    def alcove_path(self):
        r"""
        The list of extended affine Weyl group elements corresponding to the given list of folds.
        """
        path = [self._end]
        current = self._end
        for i in range(self._RY._N-1, -1, -1):
            if self._folds[i] == 1:
                current = current * self._RY._reflections[i]
            path = [current]+path
        return path

    def alcove_path_directions(self):
        r"""
        The list of finite Weyl group elements corresponding to the given list of folds.
        """
        path = self.alcove_path()
        return [x.to_dual_classical_weyl().reduced_word() for x in path]
