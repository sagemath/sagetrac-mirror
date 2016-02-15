"""
Extended affine Hecke algebra with several realizations

AUTHOR:

- Mark Shimozono (2014)

"""

#*****************************************************************************
#       Copyright (C) 2014 Mark Shimozono <mshimo at math dot vt dot edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.realizations import Category_realization_of_parent, Realizations
from sage.categories.morphism import SetMorphism
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.category import Category
from sage.categories.homset import Hom
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.tensor import tensor
from sage.sets.family import Family, FiniteFamily
from sage.misc.bindable_class import BindableClass
from sage.misc.abstract_method import abstract_method
from sage.misc.functional import is_even
from sage.algebras.multiparameter_hecke_algebra import MultiParameterHeckeAlgebra, ParameterFamilies
from sage.misc.cachefunc import cached_method
from sage.algebras.smash_product_algebra import SmashProductAlgebra
from sage.rings.integer_ring import ZZ

class ExtendedAffineHeckeAlgebra(UniqueRepresentation, Parent):
    r"""
    The possibly-extended affine not-necessarily-reduced unequal-parameter Hecke algebra with several realizations.

    ..warning: Only implemented for untwisted or dual-of-untwisted affine type

    INPUT:

        - ``cartan_type`` -- An affine or finite Cartan type (a finite Cartan type is an abbreviation for its untwisted affinization)
        - ``q1, q2`` -- parameters (default: None)
        - ``extended`` -- (default: None, which is interpreted as True) whether to use the extended affine Hecke algebra
        - ``dual_side`` -- (default: None, which is interpreted as False) whether to exchange the roles of `X` and `Y`
        - ``general_linear`` -- (default: None, which means False) If True and the root system is extended of untwisted affine type A, \
          use the general linear affine Hecke algebra.

    At first everything is explained with ``dual_side`` False.

    The notation `\tilde{X}`, `X`, `\tilde{Y}`, `Y`, `W_a(\tilde{X})`, `W_e(\tilde{X})`, `W(\tilde{X})` etc.
    is borrowed from :class:`DoubleAffineHeckeAlgebra`.

    The parameters `q1` and `q2` represent the pair of eigenvalues `q1[i]` and `q2[i]` of the algebra generator `T_i`
    for `i` in the affine Dynkin node set `I^X` of type `\tilde{X}`.
    See :func:`sage.algebras.multiparameter_hecke_algebra.ParameterFamilies`.
    To match with the notation of :class:`DoubleAffineHeckeAlgebra`, the parameters `q1[i]` and `q2[i]` are
    `v_{\alpha^X_i}` and `-1/v_{\alpha^X_i}`.

    MNEMONICS:

    - "v" -- Generally represents DAHA-duality
    - "T" -- Hecke algebra of the affine Weyl group `W(\tilde{X})`, which may or may not be extended.
    - "tv" -- Hecke algebra of the finite Weyl group `W(Y)`
    - "Lv" -- Group algebra of the lattice `Y`, which is the weight lattice `P^Y` or the root lattice `Q^Y`
       according as `\tilde{X}` is extended or not.
    - "tvLv" -- tensor product of "tv" and "Lv"
    - "Lvtv" -- tensor product of "Lv" and "tv"

    ..warning:: The implementation of "T" is always the extended algebra, but if the input ``extended`` is False,
    then only the nonextended subalgebra is used. The lattice of the group algebra "Lv"
    is implemented with the ambient space of type `Y`, which is defined over the rationals. Membership
    in the appropriate lattice is only checked at certain points.

    Supported bases:

    - "T" -- `pi_i T_w` for `w \in W_a(\tilde{X})` and `i \in I^X` a special node
    - "tvLv" -- `T_w Y^\mu` for `w \in W(Y)` and `\mu \in Y`
    - "Lvtv" -- `Y^\mu T_w` for `w \in W(Y)` and `\mu \in Y`

    For the multiplication in "T", the fundamental group acts by affine Dynkin automorphisms on the subscripts of the
    `T_i`. The multiplication for the last two bases require the Demazure-Lusztig operators for `i \in I^Y_0`, the
    nonzero Dynkin nodes for type `Y`.

    EXAMPLES::

        sage: H = ExtendedAffineHeckeAlgebra("A2")
        sage: T = H.T(); T
        T basis of The affine Hecke algebra of type ['A', 2, 1]
        sage: a = T.an_element(); a
        2*TX[0] + 3*TX[0,1] + 1 + TX[0,1,2] + 4*piX[1] TX[0] + 6*piX[1] TX[0,1] + 2*piX[1] + 2*piX[1] TX[0,1,2] + 8*piX[2] TX[0] + 12*piX[2] TX[0,1] + 4*piX[2] + 4*piX[2] TX[0,1,2]
        sage: Ty_Y = H.tvLv(); Ty_Y
        tvLv basis of The affine Hecke algebra of type ['A', 2, 1]
        sage: Ty_Y.an_element()
        2*Ty[1,2,1] Y[(2, 2, 3)] + 4*Ty[1,2] Y[(2, 2, 3)] + Y[(2, 2, 3)]
        sage: Y_Ty = H.Lvtv(); Y_Ty
        Lvtv basis of The affine Hecke algebra of type ['A', 2, 1]
        sage: Y_Ty.an_element()
        2*Y[(2, 2, 3)] Ty[1,2,1] + 4*Y[(2, 2, 3)] Ty[1,2] + Y[(2, 2, 3)]

    There are built-in coercions between the bases::

        sage: b = T.monomial((H.fundamental_group()(1),H.affine_weyl().one())); b
        piX[1]
        sage: Y_Ty(b)
        Y[(1, 0, 0)] Ty[1,2] + ((-v^2+1)/v)*Y[(1, 0, 0)] Ty[1] + ((-v^2+1)/v)*Y[(1, 0, 0)] Ty[2] + ((v^4-2*v^2+1)/v^2)*Y[(1, 0, 0)]
        sage: Ty_Y(b)
        Ty[1,2] Y[(-1, -1, 0)]
        sage: T(Ty_Y(b)) == b
        True
        sage: Y_Ty(a)
        1 + ((6*v^2-6)/v)*Y[(1, 0, 0)] Ty[1,2,1] + ((4*v^2+2*v-4)/v)*Y[(1, 0, 0)] Ty[1,2] + ((-6*v^4+2*v^3+12*v^2-2*v-6)/v^2)*Y[(1, 0, 0)] Ty[2,1] + ((-4*v^4-2*v^3+8*v^2+2*v-4)/v^2)*Y[(1, 0, 0)] Ty[1] + ((-6*v^4-2*v^3+12*v^2+2*v-6)/v^2)*Y[(1, 0, 0)] Ty[2] + ((4*v^6-4*v^5-12*v^4+8*v^3+12*v^2-4*v-4)/v^3)*Y[(1, 0, 0)] + ((8*v^2+4*v-8)/v)*Y[(1, 1, 0)] Ty[2,1] + ((-8*v^4-4*v^3+16*v^2+4*v-8)/v^2)*Y[(1, 1, 0)] Ty[1] + ((-8*v^4+8*v^3+16*v^2-8*v-8)/v^2)*Y[(1, 1, 0)] Ty[2] + ((8*v^6-8*v^5-20*v^4+16*v^3+20*v^2-8*v-8)/v^3)*Y[(1, 1, 0)] + 2*Y[(1, 0, -1)] Ty[1,2,1] + ((-2*v^2+3*v+2)/v)*Y[(1, 0, -1)] Ty[1,2] + ((-2*v^2+2)/v)*Y[(1, 0, -1)] Ty[2,1] + ((2*v^4-3*v^3-3*v^2+3*v+2)/v^2)*Y[(1, 0, -1)] Ty[1] + ((2*v^4-3*v^3-4*v^2+3*v+2)/v^2)*Y[(1, 0, -1)] Ty[2] + ((-2*v^6+3*v^5+3*v^4-6*v^3-3*v^2+3*v+2)/v^3)*Y[(1, 0, -1)] + 8*Y[(1, 0, 1)] Ty[1] + 4*Y[(1, 0, 1)] Ty[2] + ((-8*v^2+12*v+8)/v)*Y[(1, 0, 1)] + 2*Y[(0, 1, 0)] Ty[1,2,1] + ((-2*v^2+2)/v)*Y[(0, 1, 0)] Ty[1,2] + 6*Y[(0, 1, 0)] Ty[2,1] + ((-6*v^2+6)/v)*Y[(0, 1, 0)] Ty[1] + 4*Y[(0, 1, 0)] Ty[2] + ((-4*v^2+4)/v)*Y[(0, 1, 0)]
        sage: Ty_Y(a)
        ((6*v^2-6)/v)*Ty[1,2,1] Y[(-1, -1, 0)] + 2*Ty[1,2,1] Y[(-1, 0, -1)] + 2*Ty[1,2,1] Y[(-1, 0, 1)] + ((4*v^2+2*v-4)/v)*Ty[1,2] Y[(-1, -1, 0)] + 3*Ty[1,2] Y[(0, -1, 1)] + ((8*v^2+4*v-8)/v)*Ty[2,1] Y[(-1, 0, 0)] + 6*Ty[2,1] Y[(-1, -1, 0)] + ((3*v^2-3)/v)*Ty[1] + 8*Ty[1] Y[(-1, 0, 0)] + Ty[1] Y[(0, 1, -1)] + 4*Ty[2] Y[(-1, -1, 0)] + ((12*v^2-12)/v)*Ty[2] Y[(0, -1, 0)] + 4*Ty[2] Y[(0, 0, -1)] + ((2*v^2+v-2)/v) + 12*Y[(0, -1, 0)]
        sage: T(Ty_Y(a))==a
        True
        sage: K = QQ['v,vl'].fraction_field()
        sage: H = ExtendedAffineHeckeAlgebra(['C',2,1], q1=Family(dict([[0,K.gen(1)],[1,K.gen(0)],[2,K.gen(1)]])))
        sage: T = H.T()
        sage: Ty_Y = H.tvLv()
        sage: Y_Ty = H.Lvtv()
        sage: a = T.from_fundamental(H.fundamental_group()(2)); a
        piX[2]
        sage: Ty_Y(a)
        Ty[2,1,2] Y[(-1/2, -1/2)]
        sage: Y_Ty(a)
        Y[(1/2, 1/2)] Ty[2,1,2] + ((-vl^2+1)/vl)*Y[(1/2, 1/2)] Ty[1,2] + ((-vl^2+1)/vl)*Y[(1/2, 1/2)] Ty[2,1] + ((vl^4-2*vl^2+1)/vl^2)*Y[(1/2, 1/2)] Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Y[(1/2, 1/2)] Ty[2] + ((-v^2*vl^4+v^2*vl^2+vl^4-v^2-vl^2+1)/(v*vl^2))*Y[(1/2, 1/2)]
        sage: Ty_Y(a) == Ty_Y(Y_Ty(a))
        True
        sage: Ta = H.affine_hecke_algebra()
        sage: Ta[1,0,1,0] == Ta[0,1,0,1]
        True
        sage: Ty_Y(Ta[1,0,1,0]) == Ty_Y(Ta[0,1,0,1])
        True
        sage: Y_Ty(Ta[1,2,1,2]) == Y_Ty(Ta[2,1,2,1])
        True

    Here is an example with a nonreduced root system and unequal parameters.

        sage: K = QQ['v,vl,v0'].fraction_field()
        sage: v,vl,v0=K.gens()
        sage: H = ExtendedAffineHeckeAlgebra(['D',3,2], q1=Family(dict([[0,v0],[1,vl],[2,v]])), extended=False)
        sage: H._doubled_parameters
        Finite family {2: (v0^2 - 1)/v0}
        sage: T = H.T()
        sage: Ty_Y = H.tvLv()
        sage: Y_Ty = H.Lvtv()
        sage: mu = H.dual_lattice().fundamental_weight(2); mu
        (1/2, 1/2)
        sage: id = H.dual_classical_weyl().one()
        sage: a = Ty_Y.monomial((id,mu)); a
        Y[(1/2, 1/2)]
        sage: T(a)
        Traceback (most recent call last):
        ...
        ValueError: (1/2, 1/2) should be in the root lattice
        sage: a = T.from_fundamental(H.fundamental_group()(2)); a
        piX[2]
        sage: Ty_Y(a)
        Traceback (most recent call last):
        ...
        ValueError: Nontrivial fundamental group elements disallowed if the dual affine root system is nonreduced
        sage: mu = 4 * H.dual_lattice().fundamental_weight(2); mu
        (2, 2)
        sage: w = H.dual_classical_weyl().from_reduced_word([2,1,2])
        sage: b = Ty_Y.monomial((w,mu)); b
        Ty[2,1,2] Y[(2, 2)]
        sage: Y_Ty(b)
        ((v^4-2*v^2+1)/v^2)*Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Ty[2] + ((vl^2*v0^4-vl^2*v0^2-v0^4+vl^2+v0^2-1)/(vl*v0^2)) + ((v^2-1)/v)*Y[(-2, 0)] Ty[1,2] + Y[(-2, -2)] Ty[2,1,2] + ((v0^2-1)/v0)*Y[(-2, -1)] Ty[1,2] + ((v0^2-1)/v0)*Y[(-2, 1)] Ty[1,2] + ((v^2-1)/v)*Y[(-2, 2)] Ty[1,2] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(-1, 0)] Ty[1] + ((vl^2*v0^2-vl^2-v0^2+1)/(vl*v0))*Y[(-1, 0)] Ty[2] + ((v0^2-1)/v0)*Y[(-1, -2)] Ty[2,1] + ((v0^4-2*v0^2+1)/v0^2)*Y[(-1, -1)] Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Y[(-1, -1)] Ty[2] + ((vl^2-1)/vl)*Y[(-1, -1)] + ((v0^4-2*v0^2+1)/v0^2)*Y[(-1, 1)] Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Y[(-1, 1)] Ty[2] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(-1, 2)] Ty[1] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(1, 0)] Ty[1] + ((v^2*vl^2*v0^2-v^2*vl^2-v^2*v0^2-vl^2*v0^2+v^2+vl^2+v0^2-1)/(v*vl*v0))*Y[(1, 0)] + ((v0^2-1)/v0)*Y[(1, -2)] Ty[2,1] + ((v0^4-2*v0^2+1)/v0^2)*Y[(1, -1)] Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Y[(1, -1)] Ty[2] + ((v0^4-2*v0^2+1)/v0^2)*Y[(1, 1)] Ty[1] + ((v^4*vl^2-v^4-v^2*vl^2+v^2+vl^2-1)/(v^2*vl))*Y[(1, 1)] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(1, 2)] Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(2, 0)] Ty[1] + ((v^2-1)/v)*Y[(2, -2)] Ty[2,1] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(2, -1)] Ty[1] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(2, 1)] Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(2, 2)] Ty[1] + ((vl^2-1)/vl)*Y[(2, 2)] + ((v^2-1)/v)*Y[(0, -2)] Ty[2,1] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(0, -1)] Ty[1] + ((vl^2*v0^2-vl^2-v0^2+1)/(vl*v0))*Y[(0, -1)] Ty[2] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(0, 1)] Ty[1] + ((v^2*vl^2*v0^2-v^2*vl^2-v^2*v0^2-vl^2*v0^2+v^2+vl^2+v0^2-1)/(v*vl*v0))*Y[(0, 1)] + ((v^4-2*v^2+1)/v^2)*Y[(0, 2)] Ty[1]
        sage: T(b)
        TX[2,1,2,0,1,2,0,1,2,0,1,2,0,1,2]
        sage: Y_Ty(T(b)) == Y_Ty(b)
        True
        sage: Ta = H.affine_hecke_algebra()
        sage: Ta[1,0,1,0] == Ta[0,1,0,1]
        True
        sage: Ty_Y(Ta[1,0,1,0]) == Ty_Y(Ta[0,1,0,1])
        True
        sage: Y_Ty(Ta[1,2,1,2]) == Y_Ty(Ta[2,1,2,1])
        True

    The notation for the code will always match the case that ``dual_side`` is False.

    When ``dual_side`` is True then the roles of `X` and `Y` are exchanged with the following adjustments.
    One must invoke this class using the affine Cartan type `\tilde{Y}` instead of `\tilde{X}`,
    and with the parameters `v^Y_{\alpha_i}` for `i \in I^Y` instead of `v^X_{\alpha_i}`
    for `i \in I^X`. The term ``extended`` now refers to whether `W(\tilde{Y})` will be extended or not
    instead of `W(\tilde{X})`. Instead of changing `Y^\mu` for `\mu` in `Y`, to `T^X_i` for `i \in I^X`,
    we must change `X^\mu` for `\mu` in `X` to `T^Y_i` for `i \in I^Y`, and this entails a change in the
    signs of alcove walks. Finally, the Demazure Lusztig operators that control the commutation of `T^X_i`
    for `i` nonzero and `X^\mu` (which are the dual-side analogues of `T^Y_i` and `Y^\mu` above) must use the
    "dominant" convention (whereas before one uses the "antidominant" convention). Finally, the formulas
    for the fundamental group elements of `F^Y` in terms of `X^\mu` and `T^X_i` are different than those
    for `F^X` in terms of `Y^\mu` and `T^Y_i`.

    ..MATH::

        \begin{align*}
            \pi^X_i &= Y^{\omega_i^Y} T^{-1}_{u_i} = T_{u_i^{-1}} Y^{w_0(\omega_i^Y)}
            \pi^Y_i &= X^{\omega_i^X} T_{{u_i^X}^{-1}} = T^{-1}_{u_i} X^{w_0(\omega_i^X)}
        \end{align*}

        sage: Ht = ExtendedAffineHeckeAlgebra("A2", dual_side=True)
        sage: TY = Ht.T(); TY
        T basis of The affine Hecke algebra of type ['A', 2, 1] dual side
        sage: a = TY.an_element(); a
        2*TY[0] + 3*TY[0,1] + 1 + TY[0,1,2] + 4*piY[1] TY[0] + 6*piY[1] TY[0,1] + 2*piY[1] + 2*piY[1] TY[0,1,2] + 8*piY[2] TY[0] + 12*piY[2] TY[0,1] + 4*piY[2] + 4*piY[2] TY[0,1,2]
        sage: Tx_X = Ht.tvLv()
        sage: TYa = Ht.affine_hecke_algebra()
        sage: Tx_X(TY(TYa[0,1,2]))
        Tx[1] X[(0, 1, -1)] + ((-v^2+1)/v)*X[(1, 0, -1)] + ((-v^2+1)/v)*X[(0, 1, -1)]
        sage: Tx_X(TY(TYa[0])).to_opposite()
        X[(1, 0, -1)] Tx[1,2,1] + ((-v^2+1)/v)*X[(1, 0, -1)] Tx[1,2] + ((-v^2+1)/v)*X[(1, 0, -1)] Tx[2,1] + ((v^4-2*v^2+1)/v^2)*X[(1, 0, -1)] Tx[1] + ((v^4-2*v^2+1)/v^2)*X[(1, 0, -1)] Tx[2] + ((-v^6+2*v^4-2*v^2+1)/v^3)*X[(1, 0, -1)]
        sage: TYa[0,1,0] == TYa[1,0,1]
        True
        sage: TYa[0,2,0] == TYa[2,0,2]
        True
        sage: FY = Ht.fundamental_group()
        sage: FY.special_nodes()
        (0, 1, 2)
        sage: pi = Ht.tvLv().from_fundamental(FY(1)); pi
        Tx[1,2] X[(-1, -1, 0)] + ((-v^2+1)/v)*Tx[1] X[(-1, -1, 0)] + ((-v^2+1)/v)*Tx[2] X[(-1, -1, 0)] + ((v^4-2*v^2+1)/v^2)*X[(-1, -1, 0)]
        sage: pi**3
        X[(-2, -2, -2)]
        sage: pi * Ht.tvLv().from_fundamental(FY(2))
        X[(-1, -1, -1)]

    The "dominant" versus "antidominant" conventions can be seen by comparing the
    following::

        sage: Tx = Ht.dual_classical_hecke()
        sage: KX = Ht.dual_lattice_algebra()
        sage: u = Ht.dual_classical_weyl().from_reduced_word([1])
        sage: mu = Ht.dual_lattice().fundamental_weight(1)
        sage: Tx_X.monomial((u,mu))
        Tx[1] X[(1, 0, 0)]
        sage: Tx_X.monomial((u,mu)).to_opposite()
        X[(0, 1, 0)] Tx[1] + ((-v^2+1)/v)*X[(0, 1, 0)]

        sage: H = ExtendedAffineHeckeAlgebra("A2")
        sage: Ty=H.dual_classical_hecke(); KY=H.dual_lattice_algebra(); Ty_Y=H.tvLv()
        sage: u = H.dual_classical_weyl().from_reduced_word([1])
        sage: mu = -KY.basis().keys().fundamental_weight(1)
        sage: Ty_Y.monomial((u,mu))
        Ty[1] Y[(-1, 0, 0)]
        sage: Ty_Y.monomial((u,mu)).to_opposite()
        Y[(0, -1, 0)] Ty[1] + ((-v^2+1)/v)*Y[(0, -1, 0)]

    If the input `general_linear` is True, the type is untwisted type A, and the root system is extended,
    then we use the general linear extended affine Hecke algebra.

        sage: H = ExtendedAffineHeckeAlgebra(['A',2,1], general_linear = True)
        sage: T = H.T()
        sage: x = T.an_element(); x
        2*piX[1] TX[0] + 3*piX[1] TX[0,1] + piX[1] + piX[1] TX[0,1,2] + 2*piX[5] TX[0] + 3*piX[5] TX[0,1] + piX[5] + piX[5] TX[0,1,2]
        sage: z = H.tvLv()(x); z
        Ty[1,2,1] Y[(0, 1, 0)] + ((3*v^2-3)/v)*Ty[1,2,1] Y[(0, 0, 1)] + ((2*v^2+v-2)/v)*Ty[1,2] Y[(0, 0, 1)] + ((2*v^2+v-2)/v)*Ty[2,1] Y[(1, 2, 2)] + 3*Ty[2,1] Y[(0, 0, 1)] + 2*Ty[1] Y[(1, 2, 2)] + ((3*v^2-3)/v)*Ty[2] Y[(2, 1, 2)] + Ty[2] Y[(2, 2, 1)] + 2*Ty[2] Y[(0, 0, 1)] + 3*Y[(2, 1, 2)]
        sage: H.Lvtv()(x) == H.Lvtv()(z)
        True
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, q1=None, q2=None, extended=True, dual_side=False, general_linear=False):
        from sage.combinat.root_system.cartan_type import CartanType
        cartan_type = CartanType(cartan_type)
        if cartan_type.is_reducible():
            raise ValueError, "Cartan type should be irreducible"
        if cartan_type.is_finite(): # a finite Cartan type is an abbreviation for its untwisted affinization
            cartan_type = cartan_type.affine()
        elif not cartan_type.is_affine():
            raise ValueError, "Cartan type must be finite or affine"
        base_ring, q1, q2 = ParameterFamilies(cartan_type.index_set(), q1, q2)
        return super(ExtendedAffineHeckeAlgebra, cls).__classcall__(cls, cartan_type, base_ring, q1, q2, extended, dual_side, general_linear)

    def __init__(self, cartan_type, base_ring, q1, q2, extended, dual_side, general_linear):
        r"""
        INPUTS::

        - ``cartan_type`` -- A Cartan type, either of finite or affine type.
        If it is of finite type then the untwisted affinization is used.

        The following are optional parameters::

        - `q1`, `q2` -- (default: None) Optional families from the Dynkin node set to a ring.
        Their format is described in meth:`sage.algebras.multiparameter_hecke_algebra.ParameterFamilies`.
        Each of `q1` and `q2` can be a single ring element; this means all of the corresponding Hecke
        eigenvalues get set to this common value. If these are absent then generic parameters are used
        from a constructed function field over the rationals.
        - ``extended`` -- If True, use the extended affine Hecke algebra, that is, the tensor product
        with the group algebra of the special affine Dynkin automorphism group
        - ``dual_side`` -- If True, use a dual realization
        - ``general_linear`` -- If True and the root system is untwisted type A, use the general linear
        root datum.
        """
        self._base_ring = base_ring
        self._q1 = q1
        self._q2 = q2
        Parent.__init__(self, category = AlgebrasWithBasis(self._base_ring).WithRealizations())

        self._cartan_type = cartan_type
        self._extended = extended
        self._dual_side = dual_side
        if general_linear and extended and cartan_type.is_untwisted_affine() and cartan_type.type() == 'A':
            self._general_linear = True
            self._n = cartan_type.n + 1
        else:
            self._general_linear = False
        # build the components of the various realizations
        from sage.combinat.root_system.extended_affine_weyl_group import ExtendedAffineWeylGroup
        self._E = ExtendedAffineWeylGroup(cartan_type, self._general_linear, fundamental="")
        self._FW = self._E.FW()
        self._F = self._E.fundamental_group()
        self._Wa = self._E.affine_weyl()
        if cartan_type.is_untwisted_affine():
            self._cartan_type_v = cartan_type.classical().dual()
            self._cartan_type_vt = self._cartan_type_v.affine()
        else:
            self._cartan_type_v = cartan_type.classical()
            self._cartan_type_vt = self._cartan_type
        self._Lv = self._cartan_type_v.root_system().ambient_space()
        self._Wv = self._Lv.weyl_group(prefix="s")

        if self._extended:
            self._doubled_parameters = Family(dict([]))
        else:
            # find the unique nonzero doubled node of DAHA dual affine type, if it exists
            doubled_nodes = self._cartan_type_vt.doubled_nodes()
            if len(doubled_nodes) == 0 or len(doubled_nodes) == 1 and doubled_nodes[0] == 0:
                self._doubled_parameters = Family(dict([]))
            else:
                # DAHA duality forces this choice. There are at most two doubled nodes
                # and if there are two then one of them is 0.
                di = doubled_nodes[0] if doubled_nodes[0] != 0 else doubled_nodes[1]
                self._doubled_parameters = Family(dict([[di, self._q1[0] + self._q2[0]]]))
        self._dual_reduced = (len(self._doubled_parameters.keys()) == 0)

        # create the tensor factors of the realizations
        # the nonextended affine Hecke algebra with several parameters
        self._Ta = MultiParameterHeckeAlgebra(self.affine_weyl(), self.q1(), self.q2(), prefix="TY" if self._dual_side else "TX", category=AlgebrasWithBasis(self.base_ring()))
        # the group algebra of the fundamental group
        self._KF = self._F.algebra(self.base_ring(), prefix="piY" if self._dual_side else "piX")
        self._KF._print_options['bracket'] = ""
        # The classical Hecke algebra of dual type
        self._tv = MultiParameterHeckeAlgebra(self.dual_classical_weyl(), self._q1, self._q2, prefix="Tx" if self._dual_side else "Ty")
        # The group algebra of the dual translation lattice
        self._KLv = self.dual_lattice().algebra(self.base_ring(), prefix="X" if self._dual_side else "Y")

        module_category = ModulesWithBasis(self._base_ring)
        tensor_category = module_category.TensorProducts()
        # the tensor modules underlying the "T" basis
        self._TaKFmod = tensor([self._Ta, self._KF], category = tensor_category)
        self._KFTamod = tensor([self._KF, self._Ta], category = tensor_category)
        # the tensor modules underlying the "tvLv" and "Lvtv" bases
        self._tvLvmod = tensor([self._tv, self._KLv], category = tensor_category)
        self._Lvtvmod = tensor([self._KLv, self._tv], category = tensor_category)

        # create the realizations
        T = self.T()
        tvLv = self.tvLv()
        Lvtv = self.Lvtv()
        # register coercion between tvLv and Lvtv
        tvLv.register_opposite(Lvtv)

        # coercion of tv into the various bases
        tv = self.dual_classical_hecke()
        tv.module_morphism(on_basis=T.from_dual_classical_hecke_on_basis, category=module_category, codomain=T).register_as_coercion()
        tv.module_morphism(on_basis=tvLv.from_dual_classical_hecke_on_basis, category=module_category, codomain=tvLv).register_as_coercion()
        tv.module_morphism(on_basis=Lvtv.from_dual_classical_hecke_on_basis, category=module_category, codomain=Lvtv).register_as_coercion()
        
        # coercion of the group algebra of Lv into the various bases
        KLv = self.dual_lattice_algebra()
        KLv.module_morphism(on_basis=T.from_dual_lattice_algebra_on_basis, category=module_category, codomain=T).register_as_coercion()
        KLv.module_morphism(on_basis=tvLv.from_dual_lattice_algebra_on_basis, category=module_category, codomain=tvLv).register_as_coercion()
        KLv.module_morphism(on_basis=Lvtv.from_dual_lattice_algebra_on_basis, category=module_category, codomain=Lvtv).register_as_coercion()

        # coercions between T and the other two bases
        tvLv.module_morphism(on_basis = T._from_tvLv_on_basis, category=module_category, codomain=T).register_as_coercion()
        Lvtv.module_morphism(on_basis = T._from_Lvtv_on_basis, category=module_category, codomain=T).register_as_coercion()
        T.module_morphism(on_basis = T._to_tvLv_on_basis, category=module_category, codomain=tvLv).register_as_coercion()
        T.module_morphism(on_basis = T._to_Lvtv_on_basis, category=module_category, codomain=Lvtv).register_as_coercion()

    def T(self):
        r"""
        Realizes ``self`` in "T"-style.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").T()
            T basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        return self.ExtendedAffineHeckeAlgebraT()

    def tvLv(self):
        r"""
        Realizes ``self`` in "tvLv" style.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").tvLv()
            tvLv basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        return self.ExtendedAffineHeckeAlgebratvLv()

    def Lvtv(self):
        r"""
        Realizes ``self`` in "Lvtv" style.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").Lvtv()
            Lvtv basis of The affine Hecke algebra of type ['A', 2, 1]            
        """
        return self.ExtendedAffineHeckeAlgebraLvtv()

    def a_realization(self):
        r"""
        Returns the default realization.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").a_realization()
            T basis of The affine Hecke algebra of type ['A', 2, 1]            

        """
        return self.T()

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2")
            The affine Hecke algebra of type ['A', 2, 1]
            sage: ExtendedAffineHeckeAlgebra("A2", general_linear=True)
            The affine Hecke algebra of type GL(3)
        """
        return "The affine Hecke algebra of type %s"%(self.cartan_type() if not self._general_linear else "GL(%s)"%self._n)

    def base_ring(self):
        return self._base_ring

    def cartan_type(self):
        return self._cartan_type

    @cached_method
    def index_set(self):
        r"""
        Returns the index set of the affine Dynkin diagram.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("B2").index_set()
            (0, 1, 2)

        """
        return self.cartan_type().index_set()

    def extended_affine_weyl(self):
        r"""
        Returns the extended affine Weyl group of ``self``.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra(['A',2,1]).extended_affine_weyl()
            Extended affine Weyl group of type ['A', 2, 1]
            sage: ExtendedAffineHeckeAlgebra(['C',2,1]).extended_affine_weyl()
            Extended affine Weyl group of type ['C', 2, 1]
            sage: ExtendedAffineHeckeAlgebra(['D',3,2]).extended_affine_weyl()
            Extended affine Weyl group of type ['C', 2, 1]^*
        """
        return self._E

    def affine_weyl(self):
        r"""
        Returns the affine Weyl group of ``self``.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra(['A',2,1]).affine_weyl()
            Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
            sage: ExtendedAffineHeckeAlgebra(['A',5,2]).affine_weyl()
            Weyl Group of type ['B', 3, 1]^* (as a matrix group acting on the root lattice)
            sage: ExtendedAffineHeckeAlgebra(['C',2,1]).affine_weyl()
            Weyl Group of type ['C', 2, 1] (as a matrix group acting on the root lattice)

        """
        return self._Wa

    def fundamental_group(self):
        r"""
        The fundamental group.

            sage: ExtendedAffineHeckeAlgebra(['A',2,1]).fundamental_group()
            Fundamental group of type ['A', 2, 1]

        """
        return self._F

    def affine_hecke_algebra(self):
        r"""
        Returns the nonextended affine Hecke algebra.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra(['A',2,1]).affine_hecke_algebra()
            Hecke algebra of type ['A', 2, 1]
        """
        return self._Ta

    def dual_lattice(self):
        r"""
        Returns the "dual" lattice.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra(['A',2,1]).dual_lattice()
            Ambient space of the Root system of type ['A', 2]
            sage: ExtendedAffineHeckeAlgebra(['A',5,2]).dual_lattice()
            Ambient space of the Root system of type ['C', 3]
            sage: ExtendedAffineHeckeAlgebra(['C',2,1]).dual_lattice()
            Ambient space of the Root system of type ['B', 2]

        """
        return self._Lv

    def dual_lattice_algebra(self):
        r"""
        The group algebra of the "dual" lattice.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").dual_lattice_algebra()
            Group algebra of the Ambient space of the Root system of type ['A', 2] over Fraction Field of Univariate Polynomial Ring in v over Rational Field
        """
        if self._dual_side:
            prefix = "X"
        else:
            prefix = "Y"
        return self.dual_lattice().algebra(self.base_ring(), prefix=prefix)

    def dual_classical_weyl(self):
        r"""
        Returns the dual classical Weyl group of ``self``.

        This returns the classical Weyl group of type `Y` if ``self`` is of type `\tilde{X}`.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra(['A',2,1]).dual_classical_weyl()
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
            sage: ExtendedAffineHeckeAlgebra(['A',5,2]).dual_classical_weyl()
            Weyl Group of type ['C', 3] (as a matrix group acting on the ambient space)
            sage: ExtendedAffineHeckeAlgebra(['C',2,1]).dual_classical_weyl()
            Weyl Group of type ['B', 2] (as a matrix group acting on the ambient space)

        """
        return self._Wv

    def dual_classical_hecke(self):
        r"""
        The finite Hecke algebra of "dual" type.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").dual_classical_hecke()
            Hecke algebra of type ['A', 2]

        """
        if self._dual_side:
            prefix = "Tx"
        else:
            prefix = "Ty"
        I0 = self._cartan_type_v.index_set()
        return MultiParameterHeckeAlgebra(self._Wv, self._q1, self._q2, prefix=prefix)

    def dual_classical_weyl_to_affine_morphism(self, w):
        r"""
        The image of `w` under the homomorphism from the dual classical Weyl group into the affine Weyl group.

        EXAMPLES::

            sage: H = ExtendedAffineHeckeAlgebra("A2")
            sage: H.dual_classical_weyl_to_affine_morphism(H.dual_classical_weyl().an_element())
            S1*S2

        """
        return self.affine_weyl().from_reduced_word(w.reduced_word())

    def q1(self):
        r"""
        The Family of first eigenvalues for the `T_i` operators.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").q1()
            Finite family {0: v, 1: v, 2: v}

        """
        return self._q1

    def q2(self):
        r"""
        The Family of first eigenvalues for the `T_i` operators.

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").q2()
            Finite family {0: -1/v, 1: -1/v, 2: -1/v}

        """
        return self._q2

    class _BasesCategory(Category_realization_of_parent):
        r"""
        The category of realizations of an affine Hecke algebra
        """
        def __init__(self, basis, prefix=None):
            r"""
            Initialize a basis.

            INPUT:

            - ``basis`` -- a basis
            - ``prefix`` -- a prefix

            TESTS::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: bases = H._BasesCategory()
                sage: H.T() in bases
                True
            """
            Category_realization_of_parent.__init__(self, basis)
            basis._prefix = prefix

        def super_categories(self):
            r"""
            EXAMPLES::

                sage: ExtendedAffineHeckeAlgebra("A2")._BasesCategory().super_categories()
                [Category of realizations of The affine Hecke algebra of type ['A', 2, 1], Category of algebras with basis over Fraction Field of Univariate Polynomial Ring in v over Rational Field, Category of monoids with realizations, Category of additive unital additive magmas with realizations]

            """
            return [Realizations(self.base())]+self.base().category().super_categories()

        def _repr_(self):
            r"""
            Return the representation of ``self``.

            EXAMPLES::

                sage: ExtendedAffineHeckeAlgebra("A2")._BasesCategory()
                Category of bases of The affine Hecke algebra of type ['A', 2, 1]
            """
            return "Category of bases of %s" % self.base()

        class ParentMethods:

            @abstract_method(optional=False)
            def algebra_generators(self):
                r"""
                The family of generators `T_i` in the realization ``self``.

                EXAMPLES::

                    sage: H = ExtendedAffineHeckeAlgebra("A2")
                    sage: H.T().algebra_generators()
                    Finite family {0: TX[0], 1: TX[1], 2: TX[2]}
                    sage: H.tvLv().algebra_generators()
                    Finite family {0: Ty[1,2,1] Y[(-1, 0, 1)] + ((v^2-1)/v), 1: Ty[1], 2: Ty[2]}
                    sage: K = QQ['v,vl'].fraction_field(); v,vl=K.gens()
                    sage: H = ExtendedAffineHeckeAlgebra(['C',3,1],q1=Family(dict([[0,vl],[1,v],[2,v],[3,vl]])))
                    sage: H.tvLv().algebra_generators()
                    Finite family {0: Ty[1,2,3,2,1] Y[(-1, 0, 0)] + ((vl^2-1)/vl), 1: Ty[1], 2: Ty[2], 3: Ty[3]}
                """
                pass

            @abstract_method(optional=False)
            def from_dual_lattice_algebra_on_basis(self, y):
                r"""
                The image of a lattice element `y` under the morphism from the group algebra of the dual lattice into ``self``.

                EXAMPLES::

                    sage: H = ExtendedAffineHeckeAlgebra("A2")
                    sage: T = H.T()
                    sage: KY = H.dual_lattice_algebra()
                    sage: y = H.dual_lattice().simple_root(1); y
                    (1, -1, 0)
                    sage: z = T.from_dual_lattice_algebra_on_basis(y); z
                    TX[0,2,0,1] + ((-v^2+1)/v)*TX[0,2,1]                   
                    sage: Ty_Y = H.tvLv()
                    sage: Ty_Y(z)
                    Y[(1, -1, 0)]
                    sage: Ty_Y(z) == Ty_Y.from_dual_lattice_algebra_on_basis(y)
                    True
                """
                pass

            @abstract_method(optional=False)
            def from_dual_classical_hecke_on_basis(self, u):
                r"""
                Returns the image of the basis element indexed by the classical Weyl group element `u`,
                in ``self``.

                EXAMPLES::

                    sage: H = ExtendedAffineHeckeAlgebra("A2")
                    sage: u = H.dual_classical_weyl().an_element(); u
                    s1*s2
                    sage: x = H.T().from_dual_classical_hecke_on_basis(u); x
                    TX[1,2]
                    sage: x.parent()
                    T basis of The affine Hecke algebra of type ['A', 2, 1]
                """
                pass

            @abstract_method(optional=False)
            def from_affine_hecke_on_basis(self, w):
                r"""
                Return the image in ``self``, of the affine Hecke algebra element indexed by `w`.

                EXAMPLES::

                    sage: H = ExtendedAffineHeckeAlgebra("A2")
                    sage: w = H.affine_weyl().an_element(); w
                    S0*S1*S2
                    sage: H.T().from_affine_hecke_on_basis(w)
                    TX[0,1,2]
                    sage: H.Lvtv().from_affine_hecke_on_basis(w)
                    Y[(1, 0, -1)] Ty[1] + ((-v^2+1)/v)*Y[(1, 0, -1)]

                """
                pass


            @abstract_method(optional=False)
            def from_fundamental(self, pi):
                r"""
                Return the image in ``self``, of the fundamental group element indexed by `pi`.

                EXAMPLES::

                    sage: H = ExtendedAffineHeckeAlgebra("A2")
                    sage: pi = H.fundamental_group().an_element(); pi
                    [2]
                    sage: H.T().from_fundamental(pi)
                    piX[2]
                    sage: H.Lvtv().from_fundamental(pi)
                    Y[(1, 1, 0)] Ty[2,1] + ((-v^2+1)/v)*Y[(1, 1, 0)] Ty[1] + ((-v^2+1)/v)*Y[(1, 1, 0)] Ty[2] + ((v^4-2*v^2+1)/v^2)*Y[(1, 1, 0)]                    
                """
                pass

            def product_by_generator_on_basis(self, b, i, side = 'right'):
                r"""
                Returns the product of the basis element indexed by `b`, by the generator `T_i`.

                Override if there is a more efficient method for the given basis.

                EXAMPLES::

                    sage: H = ExtendedAffineHeckeAlgebra("A2")
                    sage: T = H.T()
                    sage: pi = T.factor(0).basis().keys().an_element(); pi
                    [2]
                    sage: w = H.affine_hecke_algebra().basis().keys().an_element(); w
                    S0*S1*S2
                    sage: [(i, T.product_by_generator_on_basis((pi,w), i)) for i in H.cartan_type().index_set()]
                    [(0, piX[2] TX[0,1,2,0]), (1, piX[2] TX[0,1,2,1]), (2, piX[2] TX[0,1] + ((v^2-1)/v)*piX[2] TX[0,1,2])]

                """
                if side == 'right':
                    return self.monomial(b) * self.algebra_generators()[i]
                else:
                    return self.algebra_generators()[i] * self.monomial(b)

            def product_by_fundamental_group_element_on_basis(self, b, pi, side='right'):
                r"""
                Returns the product of the basis element indexed by `b`, by the fundamental group element `pi`.

                Override if there is a more efficient method for the given basis.

                EXAMPLES::

                    sage: H = ExtendedAffineHeckeAlgebra("A2")
                    sage: T = H.T()
                    sage: pi0 = T.factor(0).basis().keys().an_element(); pi0
                    [2]
                    sage: w = T.factor(1).basis().keys().an_element(); w
                    S0*S1*S2
                    sage: [(pi, T.product_by_fundamental_group_element_on_basis((pi0,w), pi)) for pi in H.fundamental_group()]
                    [([0], piX[2] TX[0,1,2]), ([1], TX[2,0,1]), ([2], piX[1] TX[1,2,0])]

                """
                if side == 'right':
                    return self.monomial(b) * self.from_fundamental(pi)
                return self.from_fundamental(pi) * self.monomial(b)

            def signed_generator_product(self, word, signs=None):
                r"""
                The product of generators `T_i` with indices coming from the word ``word``.

                If ``signs`` is present, the sign `-1` indicates that the corresponding factor
                `T_i` should be replaced by its inverse.

                .. warning::

                    Must be implemented in style "T".

                EXAMPLES::

                    sage: ExtendedAffineHeckeAlgebra("A2").tvLv().signed_generator_product([0,2,1])
                    Ty[2] Y[(1, -1, 0)]

                """
                return self(self.realization_of().T().signed_generator_product(word,signs=signs))

    class _Bases(UniqueRepresentation, BindableClass):
        r"""
        The class of realizations of an extended affine Hecke algebra.
        """

        def _repr_(self):
            r"""
            EXAMPLES::

                sage: ExtendedAffineHeckeAlgebra("A2").T() # indirect doctest
                T basis of The affine Hecke algebra of type ['A', 2, 1]

            """
            return "%s basis of the %s"%(self._prefix, self.realization_of())

        def is_parent_of(self, x):
            return x.parent() == self

    class ExtendedAffineHeckeAlgebraT(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "T" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").T()
            T basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        def __init__(self, E):

            def ext_twist_on_basis((w, f)):
                return E._TaKFmod.monomial((f, f.inverse().act_on_affine_weyl(w)))
            SmashProductAlgebra.__init__(self, E._KF, E._Ta, twist_on_basis=ext_twist_on_basis, category=Category.join((E._BasesCategory(),AlgebrasWithBasis(E.base_ring()).TensorProducts())))
            self._style = "T"

            cat = ModulesWithBasis(E.base_ring())
            SetMorphism(Hom(E._KF,self, cat),self.factor_embedding(0)).register_as_coercion()
            SetMorphism(Hom(E._Ta,self, cat),self.factor_embedding(1)).register_as_coercion()

        def _repr_(self):
            E = self.realization_of()
            dual_side_string = " dual side" if E._dual_side else ""
            return "%s basis of %s"%(self._style, E._repr_()) + dual_side_string

        @cached_method
        def algebra_generators(self):
            r"""
            The generators `T_i` of the affine Hecke algebra.
            """
            I = self.realization_of().cartan_type().index_set()
            Ta = self.factor(1)
            return Family(dict([[i, self.factor_embedding(1)(Ta.algebra_generators()[i])] for i in I]))

        def from_affine_hecke_on_basis(self, w):
            r"""
            The basis element indexed by an element of the nonextended affine Weyl group.

            EXAMPLES::

                sage: ExtendedAffineHeckeAlgebra("A2").T().signed_generator_product([0,2,1]) # indirect doctest
                TX[0,2,1]                

            """
            H = self.realization_of()
            return self.factor_embedding(1)(self.factor(1).monomial(w))

        def signed_generator_product(self, word, signs=None):
            r"""
            The product of generators (and possibly their inverses) indexed by a word.

            EXAMPLES::

                sage: ExtendedAffineHeckeAlgebra("A2").T().signed_generator_product([0,2,1])
                TX[0,2,1]                
                sage: ExtendedAffineHeckeAlgebra("A2").T().signed_generator_product([0,2,1],[1,1,-1])
                ((-v^2+1)/v)*TX[0,2] + TX[0,2,1]

            """
            if signs is None:
                return self.from_affine_hecke_on_basis(self.realization_of().affine_weyl().from_reduced_word(word))
            return self.factor_embedding(1)(self.factor(1).signed_generator_product(word, signs=signs))

        def from_fundamental(self, f):
            r"""
            The basis element for an element of the fundamental group.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: T = H.T()
                sage: T.from_fundamental(H.fundamental_group()(2))
                piX[2]
            """
            return self.factor_embedding(0)(self.factor(0).monomial(f))

        def from_dual_classical_hecke_on_basis(self, u):
            r"""
            Returns the image of the classical Hecke basis element into ``self``.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: u = H.dual_classical_weyl().an_element(); u
                s1*s2
                sage: H.T().from_dual_classical_hecke_on_basis(u)
                TX[1,2]

            """
            return self.signed_generator_product(u.reduced_word())

        def from_dual_lattice_algebra_on_basis(self, mu):
            r"""
            The image of the `Y`-lattice group algebra element `mu` into T.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra(['A', 2, 1], general_linear=True)
                sage: mu = H.dual_lattice().an_element(); mu
                (2, 2, 3)
                sage: H.T().from_dual_lattice_algebra_on_basis(mu)
                piX[7] TX[1,0] + ((-v^2+1)/v)*piX[7] TX[0] + ((v^4-2*v^2+1)/v^2)*piX[7] + ((-v^2+1)/v)*piX[7] TX[1]            

            """
            H = self.realization_of()
            if H._general_linear:
                mu_wt = mu
            else:
                mu_wt = mu.to_weight_space(ZZ)
            pi, word, signs = H._FW(mu_wt).alcove_walk_signs()
            if not H._dual_reduced and pi != pi.parent().one():
                raise ValueError, "%s should be in the root lattice"%mu
            if H._dual_side:
                signs = tuple(signs)
            else:
                signs = tuple([-x for x in signs])
            Ta = H.affine_hecke_algebra()
            return self.from_direct_product((self.factor(0).monomial(pi), Ta.signed_generator_product(word, signs)))

        def product_by_generator_on_basis(self, b, i, side = 'right'):
            r"""
            The product of the basis element indexed by `b`, by the generator `T_i`.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: T = H.T()
                sage: pi = T.factor(0).basis().keys().an_element(); pi
                [2]
                sage: w = T.factor(1).basis().keys().an_element(); w
                S0*S1*S2
                sage: [(i, T.product_by_generator_on_basis((pi,w), i)) for i in H.cartan_type().index_set()]
                [(0, piX[2] TX[0,1,2,0]), (1, piX[2] TX[0,1,2,1]), (2, piX[2] TX[0,1] + ((v^2-1)/v)*piX[2] TX[0,1,2])]

            """
            if side == 'right':
                return self.from_direct_product((self.factor(0).monomial(b[0]), self.factor(1).product_by_generator_on_basis(b[1], i)))
            else:
                return self.algebra_generators()[i] * self.monomial(b)

        def product_by_fundamental_group_element_on_basis(self, b, pi, side = 'right'):
            r"""
            The product of the basis element indexed by `b`, with the fundamental group element `pi`.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: T = H.T()
                sage: pi0 = T.factor(0).basis().keys().an_element(); pi0
                [2]
                sage: w = T.factor(1).basis().keys().an_element(); w
                S0*S1*S2
                sage: [(pi, T.product_by_fundamental_group_element_on_basis((pi0,w), pi, side='left')) for pi in H.fundamental_group()]
                [([0], piX[2] TX[0,1,2]), ([1], TX[0,1,2]), ([2], piX[1] TX[0,1,2])]
                sage: [(pi, T.product_by_fundamental_group_element_on_basis((pi0,w), pi, side='right')) for pi in H.fundamental_group()]
                [([0], piX[2] TX[0,1,2]), ([1], TX[2,0,1]), ([2], piX[1] TX[1,2,0])]

            """
            if side == 'right':
                return self.monomial(b) * self.factor_embedding(0)(self.factor(0).monomial(pi))
            else:
                return self.monomial((pi*b[0],b[1]))

        def _from_Lvtv_on_basis(self, (mu, u)):
            r"""
            The element of T corresponding to the given basis element of Lvtv.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra(['A',2,1])
                sage: mu = H.dual_lattice().an_element()
                sage: u = H.dual_classical_weyl().an_element()
                sage: H.T()._from_Lvtv_on_basis((mu,u))
                ((-v^2+1)/v)*piX[1] TX[2] + piX[1] TX[0,1,0,2] + ((-v^2+1)/v)*piX[1] TX[0,1,2]                

            ..warning: Used to define coercion.

            """
            return self.from_dual_lattice_algebra_on_basis(mu) * self.from_dual_classical_hecke_on_basis(u)

        def _from_tvLv_on_basis(self, (u, mu)):
            r"""
            The element of T corresponding to the given basis element of Lvtv.

            EXAMPLES::


                sage: H = ExtendedAffineHeckeAlgebra(['A',2,1])
                sage: mu = H.dual_lattice().an_element()
                sage: u = H.dual_classical_weyl().an_element()
                sage: H.T()._from_tvLv_on_basis((u,mu))
                piX[1]

            ..warning: Used to define coercion.

            """
            return self.from_dual_classical_hecke_on_basis(u) * self.from_dual_lattice_algebra_on_basis(mu)

        def _to_Lvtv_on_basis(self, (pi, w)):
            r"""
            The element of Lvtv corresponding to the given basis element.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra(['A',2,1])
                sage: T = H.T()
                sage: pi = H.fundamental_group().an_element(); pi
                [2]
                sage: w = H.affine_weyl().an_element(); w
                S0*S1*S2
                sage: T._to_Lvtv_on_basis((pi,w))
                ((v^2-1)/v)*Y[(1, 1, 0)] + Y[(1, 0, 1)] Ty[2]

            ..warning: Used to define coercion.

            """
            Lvtv = self.realization_of().Lvtv()
            return Lvtv.from_fundamental(pi) * Lvtv.from_affine_hecke_on_basis(w)

        def _to_tvLv_on_basis(self, (pi, w)):
            r"""
            The element of tvLv corresponding to the given basis element.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra(['A',2,1])
                sage: T = H.T()
                sage: pi = H.fundamental_group().an_element(); pi
                [2]
                sage: w = H.affine_weyl().an_element(); w
                S0*S1*S2
                sage: T._to_tvLv_on_basis((pi,w))
                 Ty[2] Y[(0, 0, -1)]

            ..warning: Used to define coercion.

            """
            tvLv = self.realization_of().tvLv()
            return tvLv.from_fundamental(pi) * tvLv.from_affine_hecke_on_basis(w)

    class ExtendedAffineHeckeAlgebratvLv(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "tvLv" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").tvLv()
            tvLv basis of The affine Hecke algebra of type ['A', 2, 1]
        """

        def __init__(self, E):
            Lv = E.dual_lattice_algebra()
            tv = E.dual_classical_hecke()
            module_category = ModulesWithBasis(E.base_ring())
            tensor_category = module_category.TensorProducts()
            if E._dual_side:
                convention = "dominant"
            else:
                convention = "antidominant"
            self._HM = Lv.nonreduced_demazure_lusztig_operators(E.q1(), E.q2(), convention=convention, doubled_parameters=E._doubled_parameters, side="right")
            def right_action_of_Ti_on_tvLv((w,mu), i):
                smu = mu.simple_reflection(i)
                return tensor([tv.monomial(w), self._HM[i](Lv.monomial(mu)) - E.q1()[i]*Lv.monomial(smu)]) + tensor([tv.product_by_generator_on_basis(w,i), Lv.monomial(smu)], category=tensor_category)

            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            self._tvoHM = HeckeAlgebraRepresentation(E._tvLvmod, right_action_of_Ti_on_tvLv, E.cartan_type(), E.q1(), E.q2())

            def right_action_on_tvLv(ww, w, mu):
                return self._tvoHM.Tw(ww)(E._tvLvmod.monomial((w, mu)))

            SmashProductAlgebra.__init__(self, tv, Lv, right_action=right_action_on_tvLv, category=Category.join((E._BasesCategory(), AlgebrasWithBasis(E.base_ring()).TensorProducts())))
            self._style = "tvLv"

        def _repr_(self):
            E = self.realization_of()
            dual_side_string = " dual side" if E._dual_side else ""
            return "%s basis of %s"%(self._style, E._repr_()) + dual_side_string

        def from_dual_lattice_algebra_on_basis(self, y):
            r"""
            The image of `y` under the morphism from the group algebra of `Y` into ``self``.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: Y = H.dual_lattice()
                sage: y = Y.an_element(); y
                (2, 2, 3)
                sage: z = H.tvLv().from_dual_lattice_algebra_on_basis(y); z
                Y[(2, 2, 3)]
                sage: z.parent()
                tvLv basis of The affine Hecke algebra of type ['A', 2, 1]
            """
            return self.factor_embedding(1)(self.factor(1).monomial(y))

        def from_dual_classical_hecke_on_basis(self, u):
            r"""
            Returns the image of the basis element of `u` from the finite Hecke algebra into ``self``.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: Wv = H.dual_classical_weyl()
                sage: u = Wv.an_element(); u
                s1*s2
                sage: H.tvLv().from_dual_classical_hecke_on_basis(u)
                Ty[1,2]                

            """
            return self.factor_embedding(0)(self.factor(0).monomial(u))

        def T0(self):
            r"""
            The operator `T_0`.

            EXAMPLES::

                sage: ExtendedAffineHeckeAlgebra("A2").tvLv().T0()
                Ty[1,2,1] Y[(-1, 0, 1)] + ((v^2-1)/v)                
            """
            return self.realization_of().Lvtv().T0().to_opposite()

        @cached_method
        def algebra_generators(self):
            r"""
            The family of generators `T_i` in the given realization.

            EXAMPLES::

                sage: ExtendedAffineHeckeAlgebra("A2").tvLv().algebra_generators()
                Finite family {0: Ty[1,2,1] Y[(-1, 0, 1)] + ((v^2-1)/v), 1: Ty[1], 2: Ty[2]}

            """
            return Family(dict([[i, self.T0() if i == 0 else self.factor_embedding(0)(self.factor(0).algebra_generators()[i])] for i in self.realization_of().index_set()]))

        @cached_method
        def from_affine_hecke_on_basis(self, w):
            r"""
            The image of the `w`-th affine Hecke algebra basis element in "tvLv".

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: w = H.affine_weyl().an_element(); w
                S0*S1*S2
                sage: H.tvLv().from_affine_hecke_on_basis(w)
                Ty[1] Y[(0, 1, -1)]

            """
            i = w.first_descent(side="left")
            if i is None:
                return self.one()
            return self.algebra_generators()[i] * self.from_affine_hecke_on_basis(w.apply_simple_reflection(i, side="left"))

        @cached_method
        def from_fundamental(self, pi):
            r"""
            The image of the fundamental group element `pi` in "tvLv".

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: tvLv = H.tvLv()
                sage: F = H.fundamental_group()
                sage: tvLv.from_fundamental(F(1))
                Ty[1,2] Y[(-1, -1, 0)]
                sage: tvLv.from_fundamental(F(2))
                Ty[2,1] Y[(-1, 0, 0)]
                sage: tvLv.from_fundamental(F(1)) * tvLv.from_fundamental(F(2))
                Y[(-1, -1, -1)]

            ..warning:: In the ambient space of type "A2", (-1, -1, -1) and (0, 0, 0) both represent
            the zero weight.

            """
            i = pi.value()
            if i == 0:
                return self.one()
            H = self.realization_of()
            if not H._dual_reduced:
                raise ValueError, "Nontrivial fundamental group elements disallowed if the dual affine root system is nonreduced"
            # in the extended affine Weyl group, express pi as w t_mu with w in W(Y) and mu in Y.
            E = H.extended_affine_weyl()
            x = E.PvW0().from_fundamental(pi)
            rw = x.to_dual_classical_weyl().reduced_word()
            mu = x.to_dual_translation_right().to_ambient()
            HY = H.dual_classical_hecke()
            if H._dual_side:
                signs = tuple([-1 for i in range(len(rw))])
            else:
                signs = tuple([1 for i in range(len(rw))])
            return self.from_direct_product((HY.signed_generator_product(rw, signs), H.dual_lattice_algebra().monomial(mu)))

    class ExtendedAffineHeckeAlgebraLvtv(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "Lvtv" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: ExtendedAffineHeckeAlgebra("A2").Lvtv()
            Lvtv basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        def __init__(self, E):
            Lv = E.dual_lattice_algebra()
            tv = E.dual_classical_hecke()
            module_category = ModulesWithBasis(E.base_ring())
            tensor_category = module_category.TensorProducts()
            if E._dual_side:
                convention = "dominant"
            else:
                convention = "antidominant"
            self._HM = Lv.nonreduced_demazure_lusztig_operators(E.q1(), E.q2(), convention=convention, doubled_parameters=E._doubled_parameters, side="left")

            def left_action_of_Ti_on_Lvtv((mu,w), i):
                smu = mu.simple_reflection(i)
                return tensor([self._HM[i](Lv.monomial(mu)) - E.q1()[i]*Lv.monomial(smu), tv.monomial(w)]) + tensor([Lv.monomial(smu),tv.product_by_generator_on_basis(w,i,side="left")], category=tensor_category)

            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            self._LvotvHM = HeckeAlgebraRepresentation(E._Lvtvmod, left_action_of_Ti_on_Lvtv, E.cartan_type(), E.q1(), E.q2(), side="left")

            def left_action_on_Lvtv(ww, mu, w):
                return self._LvotvHM.Tw(ww)(E._Lvtvmod.monomial((mu, w)))

            SmashProductAlgebra.__init__(self, Lv, tv, left_action=left_action_on_Lvtv, category=Category.join((E._BasesCategory(), AlgebrasWithBasis(E.base_ring()).TensorProducts())))

            self._style = "Lvtv"

        def _repr_(self):
            E = self.realization_of()
            dual_side_string = " dual side" if E._dual_side else ""
            return "%s basis of %s"%(self._style, E._repr_()) + dual_side_string

        def from_dual_lattice_algebra_on_basis(self, y):
            r"""
            The image of the lattice element `y` into ``self``.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: y = H.dual_lattice().an_element(); y
                (2, 2, 3)
                sage: z = H.Lvtv().from_dual_lattice_algebra_on_basis(y); z
                Y[(2, 2, 3)]
                sage: z.parent()
                Lvtv basis of The affine Hecke algebra of type ['A', 2, 1]
            """
            return self.factor_embedding(0)(self.factor(0).monomial(y))

        def from_dual_classical_hecke_on_basis(self, u):
            r"""
            Returns the image of the `u`-th finite Hecke algebra basis element into ``self``.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: x = H.Lvtv().from_dual_classical_hecke_on_basis(H.dual_classical_weyl().an_element()); x
                Ty[1,2]
                sage: x.parent()
                Lvtv basis of The affine Hecke algebra of type ['A', 2, 1]
            """
            return self.factor_embedding(1)(self.factor(1).monomial(u))

        def T0(self):
            r"""
            The `T_0` operator.

            Let `\varphi` be the short dominant root of type `Y`. Then

            ..MATH::

                T_0^X = Y^\varphi T^{-1}_{s_\varphi}.

            """
            E = self.realization_of()
            phi = E._cartan_type_v.root_system().coroot_lattice().highest_root().associated_coroot()
            s_phi = phi.associated_reflection()
            tv = self.factor(1)
            return self.from_direct_product((self.factor(0).monomial(phi.to_ambient()), tv.signed_generator_product(s_phi, [-1 for i in range(len(s_phi))])))

        @cached_method
        def algebra_generators(self):
            r"""
            The algebra generators `T_i`.

            EXAMPLES::

                sage: ExtendedAffineHeckeAlgebra("A2").Lvtv().algebra_generators()
                Finite family {0: Y[(1, 0, -1)] Ty[1,2,1] + ((-v^2+1)/v)*Y[(1, 0, -1)] Ty[1,2] + ((-v^2+1)/v)*Y[(1, 0, -1)] Ty[2,1] + ((v^4-2*v^2+1)/v^2)*Y[(1, 0, -1)] Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(1, 0, -1)] Ty[2] + ((-v^6+2*v^4-2*v^2+1)/v^3)*Y[(1, 0, -1)], 1: Ty[1], 2: Ty[2]}

            """
            return Family(dict([[i, self.T0() if i == 0 else self.factor_embedding(1)(self.factor(1).algebra_generators()[i])] for i in self.realization_of().index_set()]))

        @cached_method
        def from_affine_hecke_on_basis(self, w):
            r"""
            The image of the `w`-th basis element of the affine Hecke algebra in Lvtv.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("A2")
                sage: w = H.affine_weyl().an_element(); w
                S0*S1*S2
                sage: H.Lvtv().from_affine_hecke_on_basis(w)
                Y[(1, 0, -1)] Ty[1] + ((-v^2+1)/v)*Y[(1, 0, -1)]            

            """
            i = w.first_descent(side="right")
            if i is None:
                return self.one()
            return self.from_affine_hecke_on_basis(w.apply_simple_reflection(i, side="right")) * self.algebra_generators()[i]

        @cached_method
        def from_fundamental(self, pi):
            r"""
            The image of a fundamental group element in Lvtv.

            EXAMPLES::

                sage: H = ExtendedAffineHeckeAlgebra("C2")
                sage: F = H.fundamental_group()
                sage: Lvtv = H.Lvtv()
                sage: Lvtv.from_fundamental(F(2))
                Y[(1/2, 1/2)] Ty[2,1,2] + ((-v^2+1)/v)*Y[(1/2, 1/2)] Ty[1,2] + ((-v^2+1)/v)*Y[(1/2, 1/2)] Ty[2,1] + ((v^4-2*v^2+1)/v^2)*Y[(1/2, 1/2)] Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(1/2, 1/2)] Ty[2] + ((-v^6+2*v^4-2*v^2+1)/v^3)*Y[(1/2, 1/2)]
                sage: Lvtv.from_fundamental(F(2))**2
                1

            """
            if pi == pi.parent().one():
                return self.one()
            H = self.realization_of()
            if not H._dual_reduced:
                raise ValueError, "Nontrivial fundamental group elements disallowed if the dual affine root system is nonreduced"
            # express pi as t_mu w with w in W(Y) and mu in Y.
            x = H.extended_affine_weyl().PvW0().from_fundamental(pi)
            rw = x.to_dual_classical_weyl().reduced_word()
            mu = x.to_dual_translation_left().to_ambient()
            tv = self.factor(1)
            if H._dual_side:
                signs = tuple([1 for i in range(len(rw))])
            else:
                signs = tuple([-1 for i in range(len(rw))])
            return self.from_direct_product((self.factor(0).monomial(mu),tv.signed_generator_product(rw, signs)))

