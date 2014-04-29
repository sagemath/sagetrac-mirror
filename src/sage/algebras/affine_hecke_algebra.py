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

def DoubledNodes(self):
    r"""
    Given an irreducible affine Cartan type, return the set of doubleable nodes,
    that is, those `i` such that `\alpha_i^\vee(\alpha_j)` is even for all affine Dynkin nodes `j`.
    """
    from sage.combinat.root_system.cartan_type import CartanType
    self = CartanType(self)
    I = self.index_set()
    M = self.cartan_matrix()
    return tuple([i for i in I if all(is_even(M[i,j]) for j in I)])

class AffineHeckeAlgebra(UniqueRepresentation, Parent):
    r"""
    The possibly-extended affine not-necessarily-reduced unequal-parameter Hecke algebra with several realizations.

    ..warning: Only implemented for untwisted or dual-of-untwisted affine type

    INPUT:

        - ``cartan_type`` -- An affine or finite Cartan type (a finite Cartan type is an abbreviation for its untwisted affinization)
        - ``q1, q2`` -- parameters (default: None)
        - ``extended`` -- (default: None, which is interpreted as True) whether to use the extended affine Hecke algebra
        - ``dual_side`` -- (default: None, which is interpreted as False) whether to exchange the roles of `X` and `Y`

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
    - "tv_Lv" -- tensor product of "tv" and "Lv"
    - "Lv_tv" -- tensor product of "Lv" and "tv"

    ..warning:: The implementation of "T" is always the extended algebra, but if the input ``extended`` is False,
    then only the nonextended subalgebra is used. The lattice of the group algebra "Lv"
    is implemented with the ambient space of type `Y`, which is defined over the rationals. Membership
    in the appropriate lattice is only checked at certain points.

    Supported bases:

    - "T" -- `pi_i T_w` for `w \in W_a(\tilde{X})` and `i \in I^X` a special node; the `\pi` term is suppressed in the nonextended case
    - "tv_Lv" -- `T_w Y^\mu` for `w \in W(Y)` and `\mu \in Y`
    - "Lv_tv" -- `Y^\mu T_w` for `w \in W(Y)` and `\mu \in Y`

    For the multiplication in "T", the fundamental group acts by affine Dynkin automorphisms on the subscripts of the
    `T_i`. The multiplication for the last two bases require the Demazure-Lusztig operators for `i \in I^Y_0`, the
    nonzero Dynkin nodes for type `Y`.

    EXAMPLES::

        sage: H = AffineHeckeAlgebra("A2")
        sage: T = H.T(); T
        T basis of The affine Hecke algebra of type ['A', 2, 1]
        sage: a = T.an_element(); a
        2*TX[0] + 3*TX[0,1] + 1 + TX[0,1,2] + 4*piX[1] TX[0] + 6*piX[1] TX[0,1] + 2*piX[1] + 2*piX[1] TX[0,1,2] + 8*piX[2] TX[0] + 12*piX[2] TX[0,1] + 4*piX[2] + 4*piX[2] TX[0,1,2]
        sage: Ty_Y = H.tv_Lv(); Ty_Y
        tv_Lv basis of The affine Hecke algebra of type ['A', 2, 1]
        sage: Ty_Y.an_element()
        Ty[1,2,1] Y[(2, 2, 3)] + 3*Ty[1,2] Y[(2, 2, 3)] + 3*Ty[2,1] Y[(2, 2, 3)]
        sage: Y_Ty = H.Lv_tv(); Y_Ty
        Lv_tv basis of The affine Hecke algebra of type ['A', 2, 1]
        sage: Y_Ty.an_element()
        Y[(2, 2, 3)] Ty[1,2,1] + 3*Y[(2, 2, 3)] Ty[1,2] + 3*Y[(2, 2, 3)] Ty[2,1]

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
        sage: v,vl=K.gens()
        sage: H = AffineHeckeAlgebra(['C',2,1], q1=Family(dict([[0,vl],[1,v],[2,vl]])))
        sage: T = H.T()
        sage: Ty_Y = H.tv_Lv()
        sage: Y_Ty = H.Lv_tv()
        sage: a = T.factor_embedding(0)(T.factor(0).monomial(H.fundamental_group()(2))); a
        piX[2]
        sage: Ty_Y(a)
        Ty[2,1,2] Y[(-1/2, -1/2)]
        sage: Y_Ty(a)
        Y[(1/2, 1/2)] Ty[2,1,2] + ((-vl^2+1)/vl)*Y[(1/2, 1/2)] Ty[1,2] + ((-vl^2+1)/vl)*Y[(1/2, 1/2)] Ty[2,1] + ((vl^4-2*vl^2+1)/vl^2)*Y[(1/2, 1/2)] Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Y[(1/2, 1/2)] Ty[2] + ((-v^2*vl^4+v^2*vl^2+vl^4-v^2-vl^2+1)/(v*vl^2))*Y[(1/2, 1/2)]
        sage: Ty_Y(a) == Ty_Y(Y_Ty(a))
        True
        sage: Ta = T.factor(1)
        sage: Ta[1,0,1,0] == Ta[0,1,0,1]
        True
        sage: Ty_Y(T.factor_embedding(1)(Ta[1,0,1,0])) == Ty_Y(T.factor_embedding(1)(Ta[0,1,0,1]))
        True
        sage: Y_Ty(T.factor_embedding(1)(Ta[1,2,1,2])) == Y_Ty(T.factor_embedding(1)(Ta[2,1,2,1]))
        True

    Here is an example with a nonreduced root system and unequal parameters.

        sage: K = QQ['v,vl,v0'].fraction_field()
        sage: v,vl,v0=K.gens()
        sage: H = AffineHeckeAlgebra(['D',3,2], q1=Family(dict([[0,v0],[1,vl],[2,v]])), extended=False)
        sage: H._doubled_parameters
        Finite family {2: (v0^2 - 1)/v0}
        sage: T = H.T()
        sage: Ty_Y = H.tv_Lv()
        sage: Y_Ty = H.Lv_tv()
        sage: mu = H.Lv().fundamental_weight(2); mu
        (1/2, 1/2)
        sage: id = H.dual_classical_weyl().one()
        sage: a = Ty_Y.monomial((id,mu)); a
        Y[(1/2, 1/2)]
        sage: T(a)
        Traceback (most recent call last):
        ...
        ValueError: (1/2, 1/2) should be in the root lattice
        sage: a = T.factor_embedding(0)(T.factor(0).monomial(H.fundamental_group()(2))); a
        piX[2]
        sage: Ty_Y(a)
        Traceback (most recent call last):
        ...
        ValueError: Nontrivial fundamental group elements disallowed if the dual affine root system is nonreduced
        sage: mu = 4 * H.Lv().fundamental_weight(2); mu
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
        sage: Ta = T.factor(1)
        sage: Ta[1,0,1,0] == Ta[0,1,0,1]
        True
        sage: Ty_Y(T.factor_embedding(1)(Ta[1,0,1,0])) == Ty_Y(T.factor_embedding(1)(Ta[0,1,0,1]))
        True
        sage: Y_Ty(T.factor_embedding(1)(Ta[1,2,1,2])) == Y_Ty(T.factor_embedding(1)(Ta[2,1,2,1]))
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

        sage: Ht = AffineHeckeAlgebra("A2", dual_side=True)
        sage: TY = Ht.T(); TY
        T basis of The affine Hecke algebra of type ['A', 2, 1] dual side
        sage: a = TY.an_element(); a
        6*TY[0] + 9*TY[0,1] + 3 + 3*TY[0,1,2] + 2*piY[1] TY[0] + 3*piY[1] TY[0,1] + piY[1] + piY[1] TY[0,1,2] + 6*piY[2] TY[0] + 9*piY[2] TY[0,1] + 3*piY[2] + 3*piY[2] TY[0,1,2]
        sage: Tx_X = Ht.tv_Lv()
        sage: TYa = TY.factor(1)
        sage: Tx_X(TY.factor_embedding(1)(TYa[0,1,2]))
        Tx[1] X[(0, 1, -1)] + ((-v^2+1)/v)*X[(1, 0, -1)] + ((-v^2+1)/v)*X[(0, 1, -1)]
        sage: Tx_X(TY.factor_embedding(1)(TYa[0])).to_opposite()
        X[(1, 0, -1)] Tx[1,2,1] + ((-v^2+1)/v)*X[(1, 0, -1)] Tx[1,2] + ((-v^2+1)/v)*X[(1, 0, -1)] Tx[2,1] + ((v^4-2*v^2+1)/v^2)*X[(1, 0, -1)] Tx[1] + ((v^4-2*v^2+1)/v^2)*X[(1, 0, -1)] Tx[2] + ((-v^6+2*v^4-2*v^2+1)/v^3)*X[(1, 0, -1)]
        sage: TYa[0,1,0] == TYa[1,0,1]
        True
        sage: TYa[0,2,0] == TYa[2,0,2]
        True
        sage: FY = Ht.fundamental_group()
        sage: FY.special_nodes()
        (0, 1, 2)
        sage: pi = Ht.F_to_tv_Lv_func(FY(1)); pi
        Tx[1,2] X[(-1, -1, 0)] + ((-v^2+1)/v)*Tx[1] X[(-1, -1, 0)] + ((-v^2+1)/v)*Tx[2] X[(-1, -1, 0)] + ((v^4-2*v^2+1)/v^2)*X[(-1, -1, 0)]
        sage: pi**3
        X[(-2, -2, -2)]
        sage: pi * Ht.F_to_tv_Lv_func(FY(2))
        X[(-1, -1, -1)]

    The "dominant" versus "antidominant" conventions can be seen by comparing the
    following::

        sage: Tx = Ht.dual_classical_hecke()
        sage: KX = Ht.Lv_algebra()
        sage: w = Tx.basis().keys().from_reduced_word([1])
        sage: mu = KX.basis().keys().fundamental_weight(1)
        sage: Tx_X.monomial((w,mu))
        Tx[1] X[(1, 0, 0)]
        sage: Tx_X.monomial((w,mu)).to_opposite()
        X[(0, 1, 0)] Tx[1] + ((-v^2+1)/v)*X[(0, 1, 0)]

        sage: H = AffineHeckeAlgebra("A2")
        sage: Ty=H.dual_classical_hecke(); KY=H.Lv_algebra(); Ty_Y=H.tv_Lv()
        sage: w = Ty.basis().keys().from_reduced_word([1])
        sage: mu = -KY.basis().keys().fundamental_weight(1)
        sage: Ty_Y.monomial((w,mu))
        Ty[1] Y[(-1, 0, 0)]
        sage: Ty_Y.monomial((w,mu)).to_opposite()
        Y[(0, -1, 0)] Ty[1] + ((-v^2+1)/v)*Y[(0, -1, 0)]

    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, q1=None, q2=None, extended=None, dual_side=None):
        from sage.combinat.root_system.cartan_type import CartanType
        cartan_type = CartanType(cartan_type)
        if isinstance(q1, dict):
            q1 = Family(q1)
        if isinstance(q2, dict):
            q2 = Family(q2)
        return super(AffineHeckeAlgebra, cls).__classcall__(cls, cartan_type, q1, q2, extended, dual_side)

    def __init__(self, cartan_type, q1, q2, extended, dual_side):
        if cartan_type.is_reducible():
            raise ValueError, "Cartan type should be irreducible"
        if cartan_type.is_finite(): # a finite Cartan type is an abbreviation for its untwisted affinization
            cartan_type = cartan_type.affine()
        elif not cartan_type.is_affine():
            raise ValueError, "Cartan type must be finite or affine"
        self._cartan_type = cartan_type
        if extended is None:
            self._extended = True
        elif extended in [True,False]:
            self._extended = extended
        else:
            raise ValueError, "%s should be a boolean"%extended
        if dual_side is None:
            self._dual_side = False
        elif dual_side not in (True,False):
            raise ValueError, "%s should be a boolean"%dual_side
        else:
            self._dual_side = dual_side
        from sage.combinat.root_system.extended_affine_weyl_group import ExtendedAffineWeylGroup
        self._We = ExtendedAffineWeylGroup(cartan_type, style="PvW0", fundamental = "piX" if not self._dual_side else "piY")
        self._FW = self._We.realization_of().FW()
        self._F = self._We.realization_of().fundamental_group()
        self._Wa = self._We.realization_of().affine_weyl()
        if cartan_type.is_untwisted_affine():
            self._cartan_type_v = cartan_type.classical().dual()
            self._cartan_type_vt = self._cartan_type_v.affine()
        else:
            self._cartan_type_v = cartan_type.classical()
            self._cartan_type_vt = self._cartan_type
        self._Lv = self._cartan_type_v.root_system().ambient_space()
        self._Wv = self._Lv.weyl_group()

        I = cartan_type.index_set()
        self._base_ring, self._q1, self._q2 = ParameterFamilies(I, q1, q2)

        Parent.__init__(self, category = AlgebrasWithBasis(self._base_ring).WithRealizations())

        if self._extended:
            self._doubled_parameters = Family(dict([]))
        else:
            # find the unique nonzero doubled node of DAHA dual affine type, if it exists
            from sage.algebras.affine_hecke_algebra import DoubledNodes
            doubled_nodes = DoubledNodes(self._cartan_type_vt)
            if len(doubled_nodes) == 0 or len(doubled_nodes) == 1 and doubled_nodes[0] == 0:
                self._doubled_parameters = Family(dict([]))
            else:
                # DAHA duality forces this choice. There are at most two doubled nodes
                # and if there are two then one of them is 0.
                di = doubled_nodes[0] if doubled_nodes[0] != 0 else doubled_nodes[1]
                self._doubled_parameters = Family(dict([[di, self._q1[0] + self._q2[0]]]))
        self._dual_reduced = (len(self._doubled_parameters.keys()) == 0)

        # create the realizations (they are cached)
        T = self.T()
        tv_Lv = self.tv_Lv()
        Lv_tv = self.Lv_tv()
        # register coercion between tv_Lv and Lv_tv
        tv_Lv.register_opposite(Lv_tv)

        # coercion of tv into the affine Hecke algebra
        tv = tv_Lv.factor(0)
        def tv_to_T_func(w):
            return T.factor_embedding(1)(T.factors()[1].monomial(self.dual_classical_weyl_to_affine_morphism(w)))
        tv_to_T = tv.module_morphism(on_basis=tv_to_T_func, category=ModulesWithBasis(self._base_ring), codomain=T)
        tv_to_T.register_as_coercion()

        # coercion of group algebra of Lv into affine Hecke algebra
        Lv = tv_Lv.factor(1)

        def T_signs(mu):
            pi, word, signs = self._FW(mu.to_weight_space(ZZ)).alcove_walk_signs()
            if not self._dual_reduced and pi != pi.parent().one():
                raise ValueError, "%s should be in the root lattice"%mu
            if self._dual_side:
                signs = tuple(signs)
            else:
                signs = tuple([-x for x in signs])
            Ta = T.factor(1)
            return T.from_direct_product((T.factor(0).monomial(pi), Ta.product_by_signed_generator_sequence(Ta.one(), word, signs)))

        Lv_to_T = Lv.module_morphism(on_basis=T_signs, category=ModulesWithBasis(self._base_ring), codomain=T)
        Lv_to_T.register_as_coercion()

        tv_Lv_to_T = tv_Lv.module_morphism(on_basis = lambda (w,mu): tv_to_T(tv.monomial(w))*Lv_to_T(Lv.monomial(mu)), category=ModulesWithBasis(self._base_ring), codomain=T)
        tv_Lv_to_T.register_as_coercion()

        Lv_tv_to_T = Lv_tv.module_morphism(on_basis = lambda (mu,w): Lv_to_T(Lv.monomial(mu))*tv_to_T(tv.monomial(w)), category=ModulesWithBasis(self._base_ring), codomain=T)
        Lv_tv_to_T.register_as_coercion()

        def T_to_tv_Lv_func((pi,w)):
            return self.F_to_tv_Lv_func(pi)*self.Ta_to_tv_Lv_func(w)

        T_to_tv_Lv = T.module_morphism(on_basis=T_to_tv_Lv_func, category=ModulesWithBasis(self._base_ring),codomain=tv_Lv)
        T_to_tv_Lv.register_as_coercion()

        def T_to_Lv_tv_func((pi,w)):
            return self.F_to_Lv_tv_func(pi)*self.Ta_to_Lv_tv_func(w)

        T_to_Lv_tv = T.module_morphism(on_basis=T_to_Lv_tv_func, category=ModulesWithBasis(self._base_ring),codomain=Lv_tv)
        T_to_Lv_tv.register_as_coercion()

        Lv_to_Lv_tv = SetMorphism(Hom(Lv,Lv_tv,ModulesWithBasis(self._base_ring)),lambda y: Lv_tv.factor_embedding(0)(y))
        Lv_to_Lv_tv.register_as_coercion()
        Lv_to_tv_Lv = SetMorphism(Hom(Lv,tv_Lv,ModulesWithBasis(self._base_ring)),lambda y: tv_Lv.factor_embedding(1)(y))
        Lv_to_tv_Lv.register_as_coercion()

    @cached_method
    def T(self):
        r"""
        Realizes ``self`` in "T"-style.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").T()
            T basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        return self.AffineHeckeAlgebraT()

    @cached_method
    def tv_Lv(self):
        r"""
        Realizes ``self`` in "tv_Lv" style.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").tv_Lv()
            tv_Lv basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        return self.AffineHeckeAlgebratv_Lv()

    @cached_method
    def Lv_tv(self):
        r"""
        Realizes ``self`` in "Lv_tv" style.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").Lv_tv()
            Lv_tv basis of The affine Hecke algebra of type ['A', 2, 1]            
        """
        return self.AffineHeckeAlgebraLv_tv()

    def a_realization(self):
        r"""
        Returns the default realization.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").a_realization()
            T basis of The affine Hecke algebra of type ['A', 2, 1]            

        """
        return self.T()

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: AffineHeckeAlgebra("A2")
            The affine Hecke algebra of type ['A', 2, 1]

        """
        return "The affine Hecke algebra of type %s"%self.cartan_type()

    def base_ring(self):
        return self._base_ring

    def cartan_type(self):
        return self._cartan_type

    @cached_method
    def index_set(self):
        r"""
        Returns the index set of the affine Dynkin diagram.

        EXAMPLES::

            sage: AffineHeckeAlgebra("B2").index_set()
            (0, 1, 2)

        """
        return self.cartan_type().index_set()

    def Lv(self):
        r"""
        Returns the "dual" lattice.

        EXAMPLES::

            sage: AffineHeckeAlgebra(['A',2,1]).Lv()
            Ambient space of the Root system of type ['A', 2]
            sage: AffineHeckeAlgebra(['A',5,2]).Lv()
            Ambient space of the Root system of type ['C', 3]
            sage: AffineHeckeAlgebra(['C',2,1]).Lv()
            Ambient space of the Root system of type ['B', 2]

        """
        return self._Lv

    @cached_method
    def Lv_algebra(self):
        r"""
        The group algebra of the "dual" lattice.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").Lv_algebra()
            Group algebra of the Ambient space of the Root system of type ['A', 2] over Fraction Field of Univariate Polynomial Ring in v over Rational Field
        """
        if self._dual_side:
            prefix = "X"
        else:
            prefix = "Y"
        return self.Lv().algebra(self.base_ring(), prefix=prefix)

    def dual_classical_weyl(self):
        r"""
        Returns the dual classical Weyl group of ``self``.

        This returns the classical Weyl group of type `Y` if ``self`` is of type `\tilde{X}`.

        EXAMPLES::

            sage: AffineHeckeAlgebra(['A',2,1]).dual_classical_weyl()
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
            sage: AffineHeckeAlgebra(['A',5,2]).dual_classical_weyl()
            Weyl Group of type ['C', 3] (as a matrix group acting on the ambient space)
            sage: AffineHeckeAlgebra(['C',2,1]).dual_classical_weyl()
            Weyl Group of type ['B', 2] (as a matrix group acting on the ambient space)

        """
        return self._Wv

    @cached_method
    def dual_classical_hecke(self):
        r"""
        The finite Hecke algebra of "dual" type.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").dual_classical_hecke()
            Hecke algebra of type ['A', 2]

        """
        if self._dual_side:
            prefix = "Tx"
        else:
            prefix = "Ty"
        I0 = self._cartan_type_v.index_set()
        return MultiParameterHeckeAlgebra(self._Wv, self._q1, self._q2, prefix=prefix)

    @cached_method
    def extended_affine_weyl(self):
        r"""
        Returns the extended affine Weyl group of ``self``.

        EXAMPLES::

            sage: AffineHeckeAlgebra(['A',2,1]).extended_affine_weyl()
            Extended affine Weyl group of type ['A', 2, 1] realized by Semidirect product of Multiplicative form of Weight lattice of the Root system of type ['A', 2] acted upon by Weyl Group of type ['A', 2] (as a matrix group acting on the weight lattice)
            sage: AffineHeckeAlgebra(['C',2,1]).extended_affine_weyl()
            Extended affine Weyl group of type ['C', 2, 1] realized by Semidirect product of Multiplicative form of Weight lattice of the Root system of type ['B', 2] acted upon by Weyl Group of type ['B', 2] (as a matrix group acting on the weight lattice)
            sage: AffineHeckeAlgebra(['D',3,2]).extended_affine_weyl()
            Extended affine Weyl group of type ['C', 2, 1]^* realized by Semidirect product of Multiplicative form of Weight lattice of the Root system of type ['B', 2] acted upon by Weyl Group of type ['B', 2] (as a matrix group acting on the weight lattice)
        """
        return self._We

    def affine_weyl(self):
        r"""
        Returns the affine Weyl group of ``self``.

        EXAMPLES::

            sage: AffineHeckeAlgebra(['A',2,1]).affine_weyl()
            Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
            sage: AffineHeckeAlgebra(['A',5,2]).affine_weyl()
            Weyl Group of type ['B', 3, 1]^* (as a matrix group acting on the root lattice)
            sage: AffineHeckeAlgebra(['C',2,1]).affine_weyl()
            Weyl Group of type ['C', 2, 1] (as a matrix group acting on the root lattice)

        """
        return self._Wa

    def fundamental_group(self):
        r"""
        The fundamental group.

            sage: AffineHeckeAlgebra(['A',2,1]).fundamental_group()
            Fundamental group of type ['A', 2, 1]

        """
        return self._F

    def dual_classical_weyl_to_affine_morphism(self, w):
        r"""
        The image of `w` under the homomorphism from the dual classical Weyl group into the affine Weyl group.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("A2")
            sage: H.dual_classical_weyl_to_affine_morphism(H.dual_classical_weyl().an_element())
            S1*S2

        """
        return self.affine_weyl().from_reduced_word(w.reduced_word())

    def q1(self):
        r"""
        The Family of first eigenvalues for the `T_i` operators.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").q1()
            Finite family {0: v, 1: v, 2: v}

        """
        return self._q1

    def q2(self):
        r"""
        The Family of first eigenvalues for the `T_i` operators.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").q2()
            Finite family {0: -1/v, 1: -1/v, 2: -1/v}

        """
        return self._q2

    @cached_method
    def Ta_to_tv_Lv_func(self, w):
        r"""
        The function from the nonextended affine Hecke algebra "Ta" to tv_Lv.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("A2")
            sage: w = H.affine_weyl().an_element(); w
            S0*S1*S2
            sage: H.Ta_to_tv_Lv_func(w)
            Ty[1] Y[(0, 1, -1)]

        """
        i = w.first_descent(side="left")
        if i is None:
            return self.tv_Lv().one()
        return self.tv_Lv().algebra_generators()[i] * self.Ta_to_tv_Lv_func(w.apply_simple_reflection(i, side="left"))

    @cached_method
    def Ta_to_Lv_tv_func(self, w):
        r"""
        The function from the nonextended affine Hecke algebra "Ta" to Lv_tv.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("A2")
            sage: w = H.affine_weyl().an_element(); w
            S0*S1*S2
            sage: H.Ta_to_Lv_tv_func(w)
            Y[(1, 0, -1)] Ty[1] + ((-v^2+1)/v)*Y[(1, 0, -1)]            

        """
        i = w.first_descent(side="right")
        if i is None:
            return self.Lv_tv().one()
        return self.Ta_to_Lv_tv_func(w.apply_simple_reflection(i, side="right")) * self.Lv_tv().algebra_generators()[i]

    @cached_method
    def F_to_tv_Lv_func(self, pi):
        r"""
        The image of a fundamental group element in tv_Lv.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("A2")
            sage: F = H.fundamental_group()
            sage: H.F_to_tv_Lv_func(F(1))
            Ty[1,2] Y[(-1, -1, 0)]
            sage: H.F_to_tv_Lv_func(F(2))
            Ty[2,1] Y[(-1, 0, 0)]
            sage: H.F_to_tv_Lv_func(F(1)) * H.F_to_tv_Lv_func(F(2))
            Y[(-1, -1, -1)]

        Note that in the crappy ambient space of type "A2", (-1, -1, -1) and (0, 0, 0) both represent
        the zero element.

        """
        i = pi.value()
        if i == 0:
            return self.tv_Lv().one()
        if not self._dual_reduced:
            raise ValueError, "Nontrivial fundamental group elements disallowed if the dual affine root system is nonreduced"
        # in the extended affine Weyl group, express pi as w t_mu with w in W(Y) and mu in Y.
        x = self.extended_affine_weyl().realization_of().W0Pv()(pi)
        rw = x.to_dual_classical_weyl().reduced_word()
        mu = x.to_dual_translation_right().to_ambient()
        tv_Lv = self.tv_Lv()
        HY = tv_Lv.factor(0)
        if self._dual_side:
            signs = tuple([-1 for i in range(len(rw))])
        else:
            signs = tuple([1 for i in range(len(rw))])
        return tv_Lv.from_direct_product((HY.product_by_signed_generator_sequence(HY.one(), rw, signs), tv_Lv.factor(1).monomial(mu)))

    @cached_method
    def F_to_Lv_tv_func(self, pi):
        r"""
        The image of a fundamental group element in tv_Lv.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("C2")
            sage: F = H.fundamental_group()
            sage: H.F_to_Lv_tv_func(F(2))
            Y[(1/2, 1/2)] Ty[2,1,2] + ((-v^2+1)/v)*Y[(1/2, 1/2)] Ty[1,2] + ((-v^2+1)/v)*Y[(1/2, 1/2)] Ty[2,1] + ((v^4-2*v^2+1)/v^2)*Y[(1/2, 1/2)] Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(1/2, 1/2)] Ty[2] + ((-v^6+2*v^4-2*v^2+1)/v^3)*Y[(1/2, 1/2)]
            sage: H.F_to_Lv_tv_func(F(2))**2
            1

        """
        i = pi.value()
        if i == 0:
            return self.Lv_tv().one()
        if not self._dual_reduced:
            raise ValueError, "Nontrivial fundamental group elements disallowed if the dual affine root system is nonreduced"
        # express pi as t_mu w with w in W(Y) and mu in Y.
        x = self._We(pi)
        rw = x.to_dual_classical_weyl().reduced_word()
        mu = x.to_dual_translation_left().to_ambient()
        Lv_tv = self.Lv_tv()
        tv = Lv_tv.factor(1)
        if self._dual_side:
            signs = tuple([1 for i in range(len(rw))])
        else:
            signs = tuple([-1 for i in range(len(rw))])
        return Lv_tv.from_direct_product((Lv_tv.factor(0).monomial(mu),tv.product_by_signed_generator_sequence(tv.one(), rw, signs)))

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

                sage: H = AffineHeckeAlgebra("A2")
                sage: bases = H._BasesCategory()
                sage: H.T() in bases
                True
            """
            Category_realization_of_parent.__init__(self, basis)
            basis._prefix = prefix

        def super_categories(self):
            r"""
            EXAMPLES::

                sage: AffineHeckeAlgebra("A2")._BasesCategory().super_categories()
                [Category of realizations of The affine Hecke algebra of type ['A', 2, 1], Category of algebras with basis over Fraction Field of Univariate Polynomial Ring in v over Rational Field, Category of monoids with realizations, Category of additive unital additive magmas with realizations]

            """
            return [Realizations(self.base())]+self.base().category().super_categories()

        def _repr_(self):
            r"""
            Return the representation of ``self``.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2")._BasesCategory()
                Category of bases of The affine Hecke algebra of type ['A', 2, 1]
            """
            return "Category of bases of %s" % self.base()

        class ParentMethods:

            @abstract_method(optional=False)
            def algebra_generators(self):
                r"""
                The family of generators `T_i` in the given realization.

                It should be a Family with key set the Dynkin node set.

                EXAMPLES::

                    sage: H = AffineHeckeAlgebra("A2")
                    sage: H.T().algebra_generators()
                    Finite family {0: TX[0], 1: TX[1], 2: TX[2]}
                    sage: H.tv_Lv().algebra_generators()
                    Finite family {0: Ty[1,2,1] Y[(-1, 0, 1)] + ((v^2-1)/v), 1: Ty[1], 2: Ty[2]}
                    sage: K = QQ['v,vl'].fraction_field(); v,vl=K.gens()
                    sage: H = AffineHeckeAlgebra(['C',3,1],q1=Family(dict([[0,vl],[1,v],[2,v],[3,vl]])))
                    sage: H.tv_Lv().algebra_generators()
                    Finite family {0: Ty[1,2,3,2,1] Y[(-1, 0, 0)] + ((vl^2-1)/vl), 1: Ty[1], 2: Ty[2], 3: Ty[3]}
                """
                pass

            @abstract_method(optional=False)
            def Lv_morphism(self, y):
                r"""
                The image of `y` under the morphism from the group algebra of the dual lattice into ``self``.

                EXAMPLES::

                    sage: H = AffineHeckeAlgebra("A2")
                    sage: T = H.T()
                    sage: KY = H.Lv_algebra()
                    sage: y = KY.monomial(KY.basis().keys().simple_root(1)); y
                    Y[(1, -1, 0)]
                    sage: z = T.Lv_morphism(y); z
                    TX[0,2,0,1] + ((-v^2+1)/v)*TX[0,2,1]                   
                    sage: Ty_Y = H.tv_Lv()
                    sage: Ty_Y(z)
                    Y[(1, -1, 0)]
                    sage: Ty_Y(z) == Ty_Y.Lv_morphism(y)
                    True
                """
                pass

            @abstract_method(optional=False)
            def dual_classical_hecke_morphism(self, a):
                r"""
                Returns the image of `a` from the finite Hecke algebra into ``self``.

                EXAMPLES::

                    sage: H = AffineHeckeAlgebra("A2")
                    sage: Ty = H.dual_classical_hecke()
                    sage: h = Ty.an_element(); h
                    Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                    sage: H.T().dual_classical_hecke_morphism(h)
                    3*TX[2,1] + 3*TX[1,2] + TX[1,2,1]
                """
                pass

            def product_by_generator_on_basis(self, b, i, side = 'right'):
                r"""
                Returns the product of the basis element indexed by `b`, by the generator `T_i`.

                Override if there is a more efficient method for the given basis.

                EXAMPLES::

                    sage: H = AffineHeckeAlgebra("A2")
                    sage: T = H.T()
                    sage: pi = T.factor(0).basis().keys().an_element(); pi
                    piX[2]
                    sage: w = T.factor(1).basis().keys().an_element(); w
                    S0*S1*S2
                    sage: [(i, T.product_by_generator_on_basis((pi,w), i)) for i in H.cartan_type().index_set()]
                    [(0, piX[2] TX[0,1,2,0]), (1, piX[2] TX[0,1,2,1]), (2, piX[2] TX[0,1] + ((v^2-1)/v)*piX[2] TX[0,1,2])]

                """
                if side == 'right':
                    return self.monomial(b) * self.algebra_generators()[i]
                else:
                    return self.algebra_generators()[i] * self.monomial(b)

            @abstract_method(optional=False)
            def product_by_fundamental_group_element_on_basis(self, b, pi, side='right'):
                r"""
                Returns the product of the basis element indexed by `b`, by the fundamental group element `pi`.

                Override if there is a more efficient method for the given basis.

                EXAMPLES::

                    sage: H = AffineHeckeAlgebra("A2")
                    sage: T = H.T()
                    sage: pi0 = T.factor(0).basis().keys().an_element(); pi0
                    piX[2]
                    sage: w = T.factor(1).basis().keys().an_element(); w
                    S0*S1*S2
                    sage: [(pi, T.product_by_fundamental_group_element_on_basis((pi0,w), pi)) for pi in H.fundamental_group()]
                    [(piX[0], piX[2] TX[0,1,2]), (piX[1], TX[2,0,1]), (piX[2], piX[1] TX[1,2,0])]

                """
                pass

            def from_reduced_word(self, word):
                r"""
                Converts an affine or finite reduced word into a group element.

                .. warning::

                    Must be implemented in style "T".

                EXAMPLES::

                    sage: AffineHeckeAlgebra("A2").tv_Lv().from_reduced_word([0,2,1])
                    Ty[2] Y[(1, -1, 0)]

                """
                return self(self.realization_of().T().from_reduced_word(word))

    class _Bases(UniqueRepresentation, BindableClass):
        r"""
        The class of realizations of an affine Hecke algebra.
        """

        def _repr_(self):
            r"""
            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").T() # indirect doctest
                T basis of The affine Hecke algebra of type ['A', 2, 1]

            """
            return "%s basis of the %s"%(self._prefix, self.realization_of())

        def is_parent_of(self, x):
            return x.parent() == self

    class AffineHeckeAlgebraT(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "T" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").T()
            T basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        def __init__(self, E):
            # the nonextended Hecke algebra of type `\tilde{X}`
            if E._dual_side:
                prefix = "TY"
            else:
                prefix = "TX"
            E._Ta = MultiParameterHeckeAlgebra(E.affine_weyl(), E.q1(), E.q2(), prefix=prefix, category=AlgebrasWithBasis(E.base_ring()))
            # the group algebra of the fundamental group
            E._KF = E._F.algebra(E.base_ring())
            E._KF._print_options['prefix'] = ""
            E._KF._print_options['bracket'] = ""
            cat = ModulesWithBasis(E.base_ring())
            mcat = cat.TensorProducts()
            E._TaoKF = tensor([E._Ta, E._KF], category = mcat)
            E._KFoTa = tensor([E._KF, E._Ta], category = mcat)
            def ext_twist_func((w, f)):
                return E._TaoKF.monomial((f, f.inverse().act_on_affine_weyl(w)))
            SmashProductAlgebra.__init__(self, E._KF, E._Ta, twist_on_basis=ext_twist_func, category=Category.join((E._BasesCategory(),AlgebrasWithBasis(E.base_ring()).TensorProducts())))
            self._style = "T"

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

        def from_reduced_word(self, word):
            r"""
            The basis element for a reduced word of affine type.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").T().from_reduced_word([0,2,1])
                TX[0,2,1]                

            """
            H = self.realization_of()
            return self.factor_embedding(1)(self.factor(1).monomial(H.affine_weyl().from_reduced_word(word)))

        def from_fundamental(self, f):
            r"""
            The basis element for an element of the fundamental group.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: T = H.T()
                sage: T.from_fundamental(H.fundamental_group()(2))
                piX[2]
            """
            return self.factor_embedding(0)(self.factor(0).monomial(f))

        def dual_classical_hecke_morphism(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: Ty = H.dual_classical_hecke()
                sage: h = Ty.an_element(); h
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                sage: H.T().dual_classical_hecke_morphism(h)
                3*TX[2,1] + 3*TX[1,2] + TX[1,2,1]

            """
            return self(a)

        def Lv_morphism(self, y):
            r"""
            The image of the `Y`-lattice group algebra element `y` into ``self``.

            EXAMPLES::

                sage: H=AffineHeckeAlgebra("A2")
                sage: KY = H.Lv_algebra()
                sage: z = H.T().Lv_morphism(KY.monomial(KY.basis().keys().simple_root(2))); z
                TX[0,1,0,2] + ((-v^2+1)/v)*TX[0,1,2]
                sage: H.Lv_tv()(z)
                Y[(0, 1, -1)]
                sage: H.tv_Lv()(z)
                Y[(0, 1, -1)]
                sage: H.Lv_tv()(H.tv_Lv()(z))
                Y[(0, 1, -1)]

            """
            return self(self.realization_of().tv_Lv()(y))

        def product_by_generator_on_basis(self, b, i, side = 'right'):
            r"""
            The product of the basis element indexed by `b`, by the generator `T_i`.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: T = H.T()
                sage: pi = T.factor(0).basis().keys().an_element(); pi
                piX[2]
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

                sage: H = AffineHeckeAlgebra("A2")
                sage: T = H.T()
                sage: pi0 = T.factor(0).basis().keys().an_element(); pi0
                piX[2]
                sage: w = T.factor(1).basis().keys().an_element(); w
                S0*S1*S2
                sage: [(pi, T.product_by_fundamental_group_element_on_basis((pi0,w), pi, side='left')) for pi in H.fundamental_group()]
                [(piX[0], piX[2] TX[0,1,2]), (piX[1], TX[0,1,2]), (piX[2], piX[1] TX[0,1,2])]
                sage: [(pi, T.product_by_fundamental_group_element_on_basis((pi0,w), pi, side='right')) for pi in H.fundamental_group()]
                [(piX[0], piX[2] TX[0,1,2]), (piX[1], TX[2,0,1]), (piX[2], piX[1] TX[1,2,0])]

            """
            if side == 'right':
                return self.monomial(b) * self.factor_embedding(0)(self.factor(0).monomial(pi))
            else:
                return self.monomial((pi*b[0],b[1]))

    class AffineHeckeAlgebratv_Lv(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "tv_Lv" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").tv_Lv()
            tv_Lv basis of The affine Hecke algebra of type ['A', 2, 1]
        """

        def __init__(self, E):
            Lv = E.Lv_algebra()
            tv = E.dual_classical_hecke()
            module_category = ModulesWithBasis(E.base_ring())
            self._Lvotv = tensor([Lv,tv],category=module_category)
            self._tvoLv = tensor([tv,Lv],category=module_category)
            if E._dual_side:
                convention = "dominant"
            else:
                convention = "antidominant"
            self._HM = Lv.nonreduced_demazure_lusztig_operators(E.q1(), E.q2(), convention=convention, doubled_parameters=E._doubled_parameters, side="right")
            def right_action_of_Ti_on_tv_Lv((w,mu), i):
                smu = mu.simple_reflection(i)
                return tensor([tv.monomial(w), self._HM[i](Lv.monomial(mu)) - E.q1()[i]*Lv.monomial(smu)]) + tensor([tv.product_by_generator_on_basis(w,i), Lv.monomial(smu)])

            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            self._tvoHM = HeckeAlgebraRepresentation(self._tvoLv, right_action_of_Ti_on_tv_Lv, E.cartan_type(), E.q1(), E.q2())

            def right_action_on_tv_Lv(ww, w, mu):
                return self._tvoHM.Tw(ww)(self._tvoLv.monomial((w, mu)))

            SmashProductAlgebra.__init__(self, tv, Lv, right_action=right_action_on_tv_Lv, category=Category.join((E._BasesCategory(), AlgebrasWithBasis(E.base_ring()).TensorProducts())))
            self._style = "tv_Lv"

        def _repr_(self):
            E = self.realization_of()
            dual_side_string = " dual side" if E._dual_side else ""
            return "%s basis of %s"%(self._style, E._repr_()) + dual_side_string

        def Lv_morphism(self, y):
            r"""
            The image of `y` under the morphism from the group algebra of `Y` into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: KY = H.Lv_algebra()
                sage: y = KY.an_element(); y
                Y[(2, 2, 3)]
                sage: H.tv_Lv().Lv_morphism(y)
                Y[(2, 2, 3)]

            """
            return self.factor_embedding(1)(y)

        def dual_classical_hecke_morphism(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: Ty = H.dual_classical_hecke()
                sage: h = Ty.an_element(); h
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                sage: H.tv_Lv().dual_classical_hecke_morphism(h)
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]

            """
            return self.factor_embedding(0)(a)

        def T0(self):
            r"""
            The operator `T_0`.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").tv_Lv().T0()
                Ty[1,2,1] Y[(-1, 0, 1)] + ((v^2-1)/v)                
            """
            return self.realization_of().Lv_tv().T0().to_opposite()

        @cached_method
        def algebra_generators(self):
            r"""
            The family of generators `T_i` in the given realization.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").tv_Lv().algebra_generators()
                Finite family {0: Ty[1,2,1] Y[(-1, 0, 1)] + ((v^2-1)/v), 1: Ty[1], 2: Ty[2]}

            """
            return Family(dict([[i, self.T0() if i == 0 else self.factor_embedding(0)(self.factor(0).algebra_generators()[i])] for i in self.realization_of().index_set()]))

        def product_by_fundamental_group_element_on_basis(self, b, pi, side='right'):
            r"""
            The product of the basis element indexed by `b`, with the fundamental group element `pi`.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: tv_Lv = H.tv_Lv()
                sage: v = tv_Lv.factor(0).basis().keys().an_element(); v.reduced_word()
                [1, 2]
                sage: mu = tv_Lv.factor(1).basis().keys().an_element(); mu
                (2, 2, 3)
                sage: [(pi, tv_Lv.product_by_fundamental_group_element_on_basis((v,mu), pi)) for pi in H.fundamental_group()]
                [(piX[0], Ty[1,2] Y[(2, 2, 3)]), (piX[1], Ty[2,1] Y[(1, 2, 2)]), (piX[2], Y[(2, 2, 2)])]

            """
            if side == 'right':
                return self.monomial(b) * self.realization_of().F_to_tv_Lv_func(pi)
            else:
                return self.realization_of().F_to_tv_Lv_func(pi) * self.monomial(b)

    class AffineHeckeAlgebraLv_tv(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "Lv_tv" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").Lv_tv()
            Lv_tv basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        def __init__(self, E):
            Lv = E.Lv_algebra()
            tv = E.dual_classical_hecke()
            module_category = ModulesWithBasis(E.base_ring())
            self._Lvotv = tensor([Lv,tv],category=module_category)
            self._tvoLv = tensor([tv,Lv],category=module_category)
            if E._dual_side:
                convention = "dominant"
            else:
                convention = "antidominant"
            self._HM = Lv.nonreduced_demazure_lusztig_operators(E.q1(), E.q2(), convention=convention, doubled_parameters=E._doubled_parameters, side="left")

            def left_action_of_Ti_on_Lv_tv((mu,w), i):
                smu = mu.simple_reflection(i)
                return tensor([self._HM[i](Lv.monomial(mu)) - E.q1()[i]*Lv.monomial(smu), tv.monomial(w)]) + tensor([Lv.monomial(smu),tv.product_by_generator_on_basis(w,i,side="left")])

            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            self._LvotvHM = HeckeAlgebraRepresentation(self._Lvotv, left_action_of_Ti_on_Lv_tv, E.cartan_type(), E.q1(), E.q2(), side="left")

            def left_action_on_Lv_tv(ww, mu, w):
                return self._LvotvHM.Tw(ww)(self._Lvotv.monomial((mu, w)))

            SmashProductAlgebra.__init__(self, Lv, tv, left_action=left_action_on_Lv_tv, category=Category.join((E._BasesCategory(), AlgebrasWithBasis(E.base_ring()).TensorProducts())))

            self._style = "Lv_tv"

        def _repr_(self):
            E = self.realization_of()
            dual_side_string = " dual side" if E._dual_side else ""
            return "%s basis of %s"%(self._style, E._repr_()) + dual_side_string

        def Lv_morphism(self, y):
            r"""
            The image of `y` under the morphism from the group algebra of the dual lattice into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: KY = H.Lv_algebra()
                sage: y = KY.an_element(); y
                Y[(2, 2, 3)]
                sage: H.Lv_tv().Lv_morphism(y)
                Y[(2, 2, 3)]
            """
            return self.factor_embedding(0)(y)

        def dual_classical_hecke_morphism(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: HY = H.dual_classical_hecke()
                sage: h = HY.an_element(); h
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                sage: H.Lv_tv().dual_classical_hecke_morphism(h)
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
            """
            return self.factor_embedding(1)(a)

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
            return self.from_direct_product((self.factor(0).monomial(phi.to_ambient()), tv.product_by_signed_generator_sequence(tv.one(), s_phi, [-1 for i in range(len(s_phi))])))

        @cached_method
        def algebra_generators(self):
            r"""
            The algebra generators `T_i`.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").Lv_tv().algebra_generators()
                Finite family {0: Y[(1, 0, -1)] Ty[1,2,1] + ((-v^2+1)/v)*Y[(1, 0, -1)] Ty[1,2] + ((-v^2+1)/v)*Y[(1, 0, -1)] Ty[2,1] + ((v^4-2*v^2+1)/v^2)*Y[(1, 0, -1)] Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(1, 0, -1)] Ty[2] + ((-v^6+2*v^4-2*v^2+1)/v^3)*Y[(1, 0, -1)], 1: Ty[1], 2: Ty[2]}

            """
            return Family(dict([[i, self.T0() if i == 0 else self.factor_embedding(1)(self.factor(1).algebra_generators()[i])] for i in self.realization_of().index_set()]))

        def product_by_fundamental_group_element_on_basis(self, b, pi, side='right'):
            r"""
            The product of the basis element indexed by `b`, with the fundamental group element `pi`.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: Lv_tv = H.Lv_tv()
                sage: mu = Lv_tv.factor(0).basis().keys().an_element(); mu
                (2, 2, 3)
                sage: v = Lv_tv.factor(1).basis().keys().an_element(); v.reduced_word()
                [1, 2]
                sage: [(pi, Lv_tv.product_by_fundamental_group_element_on_basis((mu,v), pi)) for pi in H.fundamental_group()]
                [(piX[0], Y[(2, 2, 3)] Ty[1,2]), (piX[1], Y[(2, 3, 3)] Ty[2,1] + ((-v^2+1)/v)*Y[(2, 3, 3)] Ty[1] + ((v^2-1)/v)*Y[(3, 2, 3)] Ty[1,2,1] + ((-v^4+2*v^2-1)/v^2)*Y[(3, 2, 3)] Ty[2,1] + ((-v^4+2*v^2-1)/v^2)*Y[(3, 2, 3)]), (piX[2], Y[(2, 3, 4)] + ((v^2-1)/v)*Y[(3, 2, 4)] Ty[1] + ((-v^4+2*v^2-1)/v^2)*Y[(3, 2, 4)] + ((v^2-1)/v)*Y[(3, 3, 3)] Ty[1,2,1] + ((-v^4+2*v^2-1)/v^2)*Y[(3, 3, 3)] Ty[1,2] + ((-v^4+2*v^2-1)/v^2)*Y[(3, 3, 3)])]
            """
            if side == 'right':
                return self.monomial(b) * self.realization_of().F_to_Lv_tv_func(pi)
            else:
                return self.realization_of().F_to_Lv_tv_func(pi) * self.monomial(b)
