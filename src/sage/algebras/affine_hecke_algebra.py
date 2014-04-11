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

class AffineHeckeAlgebra(UniqueRepresentation, Parent):
    r"""
    The possibly-extended affine not-necessarily-reduced unequal-parameter Hecke algebra with several realizations.

    ..warning: Only implemented for untwisted or dual-of-untwisted affine type

    INPUT:

        - ``cartan_type`` -- An affine or finite Cartan type (a finite Cartan type is an abbreviation for its untwisted affinization)
        - ``q1, q2`` -- parameters (default: None)
        - ``extended`` -- (default None, which is interpreted as True) whether to use the extended affine Hecke algebra

    A lot of notation is borrowed from :class:`DoubleAffineHeckeAlgebra`.

    The parameters `q1` and `q2` represent the pair of eigenvalues `q1[i]` and `q2[i]` of the algebra generator `T_i`
    for `i` in the affine Dynkin node set `I^X` of type `\tilde{X}`.
    See :func:`sage.algebras.multiparameter_hecke_algebra.ParameterFamilies`.
    To match with the notation of :class:`DoubleAffineHeckeAlgebra`, the parameters `q1[i]` and `q2[i]` are
    `v_{\alpha^X_i}` and `-1/v_{\alpha^X_i}`.

    MNEMONICS:

    - "AX" -- Hecke algebra of the affine Weyl group `W(\tilde{X})`, which may or may not be extended.
    - "HY" -- Hecke algebra of the finite Weyl group `W(Y)`
    - "KY" -- Group algebra of the lattice `Y`, which is the weight lattice `P^Y` or the root lattice `Q^Y`
       according as `\tilde{X}` is extended or not.
    - "HY_KY" -- tensor product of "HY" and "KY"
    - "KY_HY" -- tensor product of "KY" and "HY"

    ..warning:: The implementation of "AX" is always the extended algebra, but if the input ``extended`` is False,
    then only the nonextended subalgebra is used. "KY" is implemented with the ambient space of type `Y`,
    which is defined over the rationals. Membership is only checked at certain points.

    Supported bases:

    - "AX" -- `pi_i T_w` for `w \in W_a(\tilde{X})` and `i \in I^X` a special node; the `\pi` term is suppressed in the nonextended case
    - "HY_KY" -- `T_w Y^\mu` for `w \in W(Y)` and `\mu \in Y`
    - "KY_HY" -- `Y^\mu T_w` for `w \in W(Y)` and `\mu \in Y`

    The multiplication for the last two bases require the Demazure-Lusztig operators for `i \in I^Y_0`.

    EXAMPLES::

        sage: H = AffineHeckeAlgebra("A2")
        sage: AX = H.AX(); AX
        AX basis of The affine Hecke algebra of type ['A', 2, 1]
        sage: a = AX.an_element(); a
        2*[piX[2]] TX[0] + 3*[piX[2]] TX[0,1] + [piX[2]] + [piX[2]] TX[0,1,2]
        sage: HY_KY = H.HY_KY(); HY_KY
        HY_KY basis of The affine Hecke algebra of type ['A', 2, 1]
        sage: HY_KY.an_element()
        Ty[1,2,1] Y[(2, 2, 3)] + 3*Ty[1,2] Y[(2, 2, 3)] + 3*Ty[2,1] Y[(2, 2, 3)]
        sage: KY_HY = H.KY_HY(); KY_HY
        KY_HY basis of The affine Hecke algebra of type ['A', 2, 1]
        sage: KY_HY.an_element()
        Y[(2, 2, 3)] Ty[1,2,1] + 3*Y[(2, 2, 3)] Ty[1,2] + 3*Y[(2, 2, 3)] Ty[2,1]

    There are built-in coercions between the bases::

        sage: b = AX.monomial((H.fundamental_group()(1),H.affine_weyl().one())); b
        [piX[1]]
        sage: KY_HY(b)
        Y[(1, 0, 0)] Ty[1,2] + ((-v^2+1)/v)*Y[(1, 0, 0)] Ty[1] + ((-v^2+1)/v)*Y[(1, 0, 0)] Ty[2] + ((v^4-2*v^2+1)/v^2)*Y[(1, 0, 0)]
        sage: HY_KY(b)
        Ty[1,2] Y[(-1, -1, 0)]
        sage: AX(HY_KY(b)) == b
        True

        sage: KY_HY(a)
        ((2*v^2+v-2)/v)*Y[(1, 1, 0)] Ty[2,1] + ((-2*v^4-v^3+4*v^2+v-2)/v^2)*Y[(1, 1, 0)] Ty[1] + ((-2*v^4+2*v^3+4*v^2-2*v-2)/v^2)*Y[(1, 1, 0)] Ty[2] + ((2*v^6-2*v^5-5*v^4+4*v^3+5*v^2-2*v-2)/v^3)*Y[(1, 1, 0)] + 2*Y[(1, 0, 1)] Ty[1] + Y[(1, 0, 1)] Ty[2] + ((-2*v^2+3*v+2)/v)*Y[(1, 0, 1)]
        sage: HY_KY(a)
        ((2*v^2+v-2)/v)*Ty[2,1] Y[(-1, 0, 0)] + 2*Ty[1] Y[(-1, 0, 0)] + ((3*v^2-3)/v)*Ty[2] Y[(0, -1, 0)] + Ty[2] Y[(0, 0, -1)] + 3*Y[(0, -1, 0)]
        sage: AX(HY_KY(a))==a
        True

        sage: K = QQ['v,vl'].fraction_field()
        sage: v,vl=K.gens()
        sage: H = AffineHeckeAlgebra(['C',2,1], q1=Family(dict([[0,vl],[1,v],[2,vl]])))
        sage: AX = H.AX()
        sage: HY_KY = H.HY_KY()
        sage: KY_HY = H.KY_HY()
        sage: a = AX.factor_embedding(0)(AX.factors()[0].monomial(H.fundamental_group()(2))); a
        [piX[2]]
        sage: HY_KY(a)
        Ty[2,1,2] Y[(-1/2, -1/2)]
        sage: KY_HY(a)
        Y[(1/2, 1/2)] Ty[2,1,2] + ((-vl^2+1)/vl)*Y[(1/2, 1/2)] Ty[1,2] + ((-vl^2+1)/vl)*Y[(1/2, 1/2)] Ty[2,1] + ((vl^4-2*vl^2+1)/vl^2)*Y[(1/2, 1/2)] Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Y[(1/2, 1/2)] Ty[2] + ((-v^2*vl^4+v^2*vl^2+vl^4-v^2-vl^2+1)/(v*vl^2))*Y[(1/2, 1/2)]
        sage: HY_KY(a) == HY_KY(KY_HY(a))
        True
        sage: AXa = AX.factors()[1]
        sage: AXa[1,0,1,0] == AXa[0,1,0,1]
        True
        sage: HY_KY(AX.factor_embedding(1)(AXa[1,0,1,0])) == HY_KY(AX.factor_embedding(1)(AXa[0,1,0,1]))
        True
        sage: KY_HY(AX.factor_embedding(1)(AXa[1,2,1,2])) == KY_HY(AX.factor_embedding(1)(AXa[2,1,2,1]))
        True

    Here is an example with a nonreduced root system and unequal parameters.

        sage: K = QQ['v,vl,v0'].fraction_field()
        sage: v,vl,v0=K.gens()
        sage: H = AffineHeckeAlgebra(['D',3,2], q1=Family(dict([[0,v0],[1,vl],[2,v]])), extended=False)
        sage: H._doubled_parameters
        Finite family {2: (v0^2 - 1)/v0}
        sage: AX = H.AX()
        sage: HY_KY = H.HY_KY()
        sage: KY_HY = H.KY_HY()
        sage: mu = H.Y_lattice().fundamental_weight(2); mu
        (1/2, 1/2)
        sage: id = H.classical_weyl().one()
        sage: a = HY_KY.monomial((id,mu)); a
        Y[(1/2, 1/2)]
        sage: AX(a)
        Traceback (most recent call last):
        ...
        ValueError: (1/2, 1/2) should be in the root lattice
        sage: a = AX.factor_embedding(0)(AX.factors()[0].monomial(H.fundamental_group()(2))); a
        [piX[2]]
        sage: HY_KY(a)
        Traceback (most recent call last):
        ...
        ValueError: Nontrivial fundamental group elements disallowed if the dual affine root system is nonreduced
        sage: mu = 4 * H.Y_lattice().fundamental_weight(2); mu
        (2, 2)
        sage: w = H.classical_weyl().from_reduced_word([2,1,2])
        sage: b = HY_KY.monomial((w,mu)); b
        Ty[2,1,2] Y[(2, 2)]
        sage: KY_HY(b)
        ((v^4-2*v^2+1)/v^2)*Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Ty[2] + ((vl^2*v0^4-vl^2*v0^2-v0^4+vl^2+v0^2-1)/(vl*v0^2)) + ((v^2-1)/v)*Y[(-2, 0)] Ty[1,2] + Y[(-2, -2)] Ty[2,1,2] + ((v0^2-1)/v0)*Y[(-2, -1)] Ty[1,2] + ((v0^2-1)/v0)*Y[(-2, 1)] Ty[1,2] + ((v^2-1)/v)*Y[(-2, 2)] Ty[1,2] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(-1, 0)] Ty[1] + ((vl^2*v0^2-vl^2-v0^2+1)/(vl*v0))*Y[(-1, 0)] Ty[2] + ((v0^2-1)/v0)*Y[(-1, -2)] Ty[2,1] + ((v0^4-2*v0^2+1)/v0^2)*Y[(-1, -1)] Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Y[(-1, -1)] Ty[2] + ((vl^2-1)/vl)*Y[(-1, -1)] + ((v0^4-2*v0^2+1)/v0^2)*Y[(-1, 1)] Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Y[(-1, 1)] Ty[2] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(-1, 2)] Ty[1] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(1, 0)] Ty[1] + ((v^2*vl^2*v0^2-v^2*vl^2-v^2*v0^2-vl^2*v0^2+v^2+vl^2+v0^2-1)/(v*vl*v0))*Y[(1, 0)] + ((v0^2-1)/v0)*Y[(1, -2)] Ty[2,1] + ((v0^4-2*v0^2+1)/v0^2)*Y[(1, -1)] Ty[1] + ((v^2*vl^2-v^2-vl^2+1)/(v*vl))*Y[(1, -1)] Ty[2] + ((v0^4-2*v0^2+1)/v0^2)*Y[(1, 1)] Ty[1] + ((v^4*vl^2-v^4-v^2*vl^2+v^2+vl^2-1)/(v^2*vl))*Y[(1, 1)] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(1, 2)] Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(2, 0)] Ty[1] + ((v^2-1)/v)*Y[(2, -2)] Ty[2,1] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(2, -1)] Ty[1] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(2, 1)] Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(2, 2)] Ty[1] + ((vl^2-1)/vl)*Y[(2, 2)] + ((v^2-1)/v)*Y[(0, -2)] Ty[2,1] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(0, -1)] Ty[1] + ((vl^2*v0^2-vl^2-v0^2+1)/(vl*v0))*Y[(0, -1)] Ty[2] + ((v^2*v0^2-v^2-v0^2+1)/(v*v0))*Y[(0, 1)] Ty[1] + ((v^2*vl^2*v0^2-v^2*vl^2-v^2*v0^2-vl^2*v0^2+v^2+vl^2+v0^2-1)/(v*vl*v0))*Y[(0, 1)] + ((v^4-2*v^2+1)/v^2)*Y[(0, 2)] Ty[1]
        sage: AX(b)
        TX[2,1,2,0,1,2,0,1,2,0,1,2,0,1,2]
        sage: KY_HY(AX(b)) == KY_HY(b)
        True
        sage: AXa = AX.factors()[1]
        sage: AXa[1,0,1,0] == AXa[0,1,0,1]
        True
        sage: HY_KY(AX.factor_embedding(1)(AXa[1,0,1,0])) == HY_KY(AX.factor_embedding(1)(AXa[0,1,0,1]))
        True
        sage: KY_HY(AX.factor_embedding(1)(AXa[1,2,1,2])) == KY_HY(AX.factor_embedding(1)(AXa[2,1,2,1]))
        True

    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, q1=None, q2=None, extended=None, prefix=None):
        from sage.combinat.root_system.cartan_type import CartanType
        cartan_type = CartanType(cartan_type)
        if isinstance(q1, dict):
            q1 = Family(q1)
        if isinstance(q2, dict):
            q2 = Family(q2)
        return super(AffineHeckeAlgebra, cls).__classcall__(cls, cartan_type, q1, q2, extended, prefix)

    def __init__(self, cartan_type, q1, q2, extended, prefix):
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
        if prefix is None:
            prefix = "piX"
        from sage.combinat.root_system.extended_affine_weyl_group import ExtendedAffineWeylGroup
        self._WeX = ExtendedAffineWeylGroup(cartan_type, style="PvW0", fundamental = prefix)
        self._FW = self._WeX.realization_of().FW()
        self._FX = self._WeX.realization_of().fundamental_group()
        self._WaX = self._WeX.realization_of().affine_weyl()
        if cartan_type.is_untwisted_affine():
            self._cartan_type_Y = cartan_type.classical().dual()
            self._cartan_type_Yt = self._cartan_type_Y.affine()
        else:
            self._cartan_type_Y = cartan_type.classical()
            self._cartan_type_Yt = self._cartan_type
        self._LY = self._cartan_type_Y.root_system().ambient_space()
        self._WY = self._LY.weyl_group()

        I = cartan_type.index_set()
        self._base_ring, self._q1, self._q2 = ParameterFamilies(I, q1, q2)

        if self._extended:
            self._doubled_parameters = Family(dict([]))
        else:
            # find the unique finite doubled node of type `\tilde{Y}`, if it exists
            di = None
            I0 = self._cartan_type_Y.index_set()
            cartan_matrix = self._cartan_type_Yt.cartan_matrix()
            for i in I0:
                if all(is_even(cartan_matrix[i,j]) for j in I0):
                    di = i
                    break
            if di is None:
                self._doubled_parameters = Family(dict([]))
            else:
                # DAHA duality forces this choice
                self._doubled_parameters = Family(dict([[di, self._q1[0] + self._q2[0]]]))
        self._dual_reduced = (len(self._doubled_parameters.keys()) == 0)

        Parent.__init__(self, category = AlgebrasWithBasis(self._base_ring).WithRealizations())

        # create the realizations (they are cached)
        AX = self.AX()
        HY_KY = self.HY_KY()
        KY_HY = self.KY_HY()
        # register coercion between HY_KY and KY_HY
        HY_KY.register_opposite(KY_HY)

        # coercion of HY into the affine Hecke algebra
        HY = HY_KY.factors()[0]
        def HY_to_AX_func(w):
            return AX.factor_embedding(1)(AX.factors()[1].monomial(self.classical_weyl_to_affine_morphism(w)))
        HY_to_AX = HY.module_morphism(on_basis=HY_to_AX_func, category=ModulesWithBasis(self._base_ring), codomain=AX)
        HY_to_AX.register_as_coercion()

        # coercion of group algebra of Y into affine Hecke algebra
        KY = HY_KY.factors()[1]

        def T_signs(mu):
            pi, word, signs = self._FW(mu.to_weight_space(ZZ)).alcove_walk_signs()
            if not self._dual_reduced and pi != pi.parent().one():
                raise ValueError, "%s should be in the root lattice"%mu
            AXa = AX.factors()[1]
            return AX.from_direct_product(AX.factors()[0].monomial(pi), AXa.product_by_signed_generator_sequence(AXa.one(), word, tuple([-x for x in signs])))

        KY_to_AX = KY.module_morphism(on_basis=T_signs, category=ModulesWithBasis(self._base_ring), codomain=AX)
        KY_to_AX.register_as_coercion()

        HY_KY_to_AX = HY_KY.module_morphism(on_basis = lambda (w,mu): HY_to_AX(HY.monomial(w))*KY_to_AX(KY.monomial(mu)), category=ModulesWithBasis(self._base_ring), codomain=AX)
        HY_KY_to_AX.register_as_coercion()

        KY_HY_to_AX = KY_HY.module_morphism(on_basis = lambda (mu,w): KY_to_AX(KY.monomial(mu))*HY_to_AX(HY.monomial(w)), category=ModulesWithBasis(self._base_ring), codomain=AX)
        KY_HY_to_AX.register_as_coercion()

        def AX_to_HY_KY_func((pi,w)):
            return self.FX_to_HY_KY_func(pi)*self.AXa_to_HY_KY_func(w)

        AX_to_HY_KY = AX.module_morphism(on_basis=AX_to_HY_KY_func, category=ModulesWithBasis(self._base_ring),codomain=HY_KY)
        AX_to_HY_KY.register_as_coercion()

        def AX_to_KY_HY_func((pi,w)):
            return self.FX_to_KY_HY_func(pi)*self.AXa_to_KY_HY_func(w)

        AX_to_KY_HY = AX.module_morphism(on_basis=AX_to_KY_HY_func, category=ModulesWithBasis(self._base_ring),codomain=KY_HY)
        AX_to_KY_HY.register_as_coercion()

        KY_to_KY_HY = SetMorphism(Hom(KY,KY_HY,ModulesWithBasis(self._base_ring)),lambda y: KY_HY.factor_embedding(0)(y))
        KY_to_KY_HY.register_as_coercion()
        KY_to_HY_KY = SetMorphism(Hom(KY,HY_KY,ModulesWithBasis(self._base_ring)),lambda y: HY_KY.factor_embedding(1)(y))
        KY_to_HY_KY.register_as_coercion()

    @cached_method
    def AX(self):
        r"""
        Realizes ``self`` in "AX"-style.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").AX()
            AX basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        return self.AffineHeckeAlgebraAX()

    @cached_method
    def HY_KY(self):
        r"""
        Realizes ``self`` in "HY_KY" style.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").HY_KY()
            HY_KY basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        return self.AffineHeckeAlgebraHY_KY()

    @cached_method
    def KY_HY(self):
        r"""
        Realizes ``self`` in "KY_HY" style.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").KY_HY()
            KY_HY basis of The affine Hecke algebra of type ['A', 2, 1]            
        """
        return self.AffineHeckeAlgebraKY_HY()

    def a_realization(self):
        r"""
        Returns the default realization.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").a_realization()
            AX basis of The affine Hecke algebra of type ['A', 2, 1]            

        """
        return self.AX()

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

    def Y_lattice(self):
        r"""
        Returns the `Y` lattice.

        EXAMPLES::

            sage: AffineHeckeAlgebra(['A',2,1]).Y_lattice()
            Ambient space of the Root system of type ['A', 2]
            sage: AffineHeckeAlgebra(['A',5,2]).Y_lattice()
            Ambient space of the Root system of type ['C', 3]
            sage: AffineHeckeAlgebra(['C',2,1]).Y_lattice()
            Ambient space of the Root system of type ['B', 2]

        """
        return self._LY

    @cached_method
    def KY(self, prefix=None):
        r"""
        The group algebra of the `Y` lattice.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").KY()
            Group algebra of the Ambient space of the Root system of type ['A', 2] over Fraction Field of Univariate Polynomial Ring in v over Rational Field
        """
        if prefix is None:
            prefix = "Y"
        return self.Y_lattice().algebra(self.base_ring(), prefix=prefix)

    def classical_weyl(self):
        r"""
        Returns the classical Weyl group of ``self``.

        This returns the classical Weyl group `W(Y)`.

        EXAMPLES::

            sage: AffineHeckeAlgebra(['A',2,1]).classical_weyl()
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
            sage: AffineHeckeAlgebra(['A',5,2]).classical_weyl()
            Weyl Group of type ['C', 3] (as a matrix group acting on the ambient space)
            sage: AffineHeckeAlgebra(['C',2,1]).classical_weyl()
            Weyl Group of type ['B', 2] (as a matrix group acting on the ambient space)

        """
        return self._WY

    @cached_method
    def HY(self, prefix=None):
        r"""
        The finite Hecke algebra of type `Y`.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").HY()
            Hecke algebra of type ['A', 2]

        """
        if prefix is None:
            prefix = "Ty"
        I0 = self._cartan_type_Y.index_set()
        return MultiParameterHeckeAlgebra(self._WY, self._q1, self._q2, prefix=prefix)

    @cached_method
    def extended_affine_weyl(self, prefix=None):
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
        return self._WeX

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
        return self._WaX

    def fundamental_group(self):
        r"""
        The fundamental group.

            sage: AffineHeckeAlgebra(['A',2,1]).fundamental_group()
            Fundamental group of type ['A', 2, 1]

        """
        return self._FX

    def classical_weyl_to_affine_morphism(self, w):
        r"""
        The image of `w` under the homomorphism from the classical Weyl group into the affine Weyl group.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("A2")
            sage: H.classical_weyl_to_affine_morphism(H.classical_weyl().an_element())
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
    def AXa_to_HY_KY_func(self, w):
        r"""
        The function from the nonextended affine Hecke algebra "AXa" to HY_KY.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("A2")
            sage: w = H.affine_weyl().an_element(); w
            S0*S1*S2
            sage: H.AXa_to_HY_KY_func(w)
            Ty[1] Y[(0, 1, -1)]

        """
        i = w.first_descent(side="left")
        if i is None:
            return self.HY_KY().one()
        return self.HY_KY().algebra_generators()[i] * self.AXa_to_HY_KY_func(w.apply_simple_reflection(i, side="left"))

    @cached_method
    def AXa_to_KY_HY_func(self, w):
        r"""
        The function from the nonextended affine Hecke algebra "AXa" to KY_HY.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("A2")
            sage: w = H.affine_weyl().an_element(); w
            S0*S1*S2
            sage: H.AXa_to_KY_HY_func(w)
            Y[(1, 0, -1)] Ty[1] + ((-v^2+1)/v)*Y[(1, 0, -1)]            

        """
        i = w.first_descent(side="right")
        if i is None:
            return self.KY_HY().one()
        return self.AXa_to_KY_HY_func(w.apply_simple_reflection(i, side="right")) * self.KY_HY().algebra_generators()[i]

    @cached_method
    def FX_to_HY_KY_func(self, pi):
        r"""
        The image of a fundamental group element in HY_KY.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("A2")
            sage: F = H.fundamental_group()
            sage: H.FX_to_HY_KY_func(F(1))
            Ty[1,2] Y[(-1, -1, 0)]
            sage: H.FX_to_HY_KY_func(F(2))
            Ty[2,1] Y[(-1, 0, 0)]
            sage: H.FX_to_HY_KY_func(F(1)) * H.FX_to_HY_KY_func(F(2))
            Y[(-1, -1, -1)]

        Note that in the crappy ambient space of type "A2", (-1, -1, -1) and (0, 0, 0) both represent
        the zero element.

        """
        i = pi.value()
        if i == 0:
            return self.HY_KY().one()
        if not self._dual_reduced:
            raise ValueError, "Nontrivial fundamental group elements disallowed if the dual affine root system is nonreduced"
        # in the extended affine Weyl group, express pi as w t_mu with w in W(Y) and mu in Y.
        x = self.extended_affine_weyl().realization_of().W0Pv()(pi)
        w = x.to_dual_classical_weyl().reduced_word()
        mu = x.to_dual_translation_right().to_ambient()
        HY_KY = self.HY_KY()
        HY = HY_KY.factors()[0]
        return HY_KY.from_direct_product(HY.product_by_signed_generator_sequence(HY.one(), w), HY_KY.factors()[1].monomial(mu))

    @cached_method
    def FX_to_KY_HY_func(self, pi):
        r"""
        The image of a fundamental group element in HY_KY.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("C2")
            sage: F = H.fundamental_group()
            sage: H.FX_to_KY_HY_func(F(2))
            Y[(1/2, 1/2)] Ty[2,1,2] + ((-v^2+1)/v)*Y[(1/2, 1/2)] Ty[1,2] + ((-v^2+1)/v)*Y[(1/2, 1/2)] Ty[2,1] + ((v^4-2*v^2+1)/v^2)*Y[(1/2, 1/2)] Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(1/2, 1/2)] Ty[2] + ((-v^6+2*v^4-2*v^2+1)/v^3)*Y[(1/2, 1/2)]
            sage: H.FX_to_KY_HY_func(F(2))**2
            1

        """
        i = pi.value()
        if i == 0:
            return self.KY_HY().one()
        if not self._dual_reduced:
            raise ValueError, "Nontrivial fundamental group elements disallowed if the dual affine root system is nonreduced"
        # express pi as t_mu w with w in W(Y) and mu in Y.
        x = self._WeX(pi)
        w = x.to_dual_classical_weyl().reduced_word()
        mu = x.to_dual_translation_left().to_ambient()
        signs = tuple([-1 for i in range(len(w))])
        KY_HY = self.KY_HY()
        HY = KY_HY.factors()[1]
        return KY_HY.from_direct_product(KY_HY.factors()[0].monomial(mu),HY.product_by_signed_generator_sequence(HY.one(), w, signs))

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
                sage: H.AX() in bases
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
                    sage: H.AX().algebra_generators()
                    Finite family {0: TX[0], 1: TX[1], 2: TX[2]}
                    sage: H.HY_KY().algebra_generators()
                    Finite family {0: Ty[1,2,1] Y[(-1, 0, 1)] + ((v^2-1)/v), 1: Ty[1], 2: Ty[2]}
                    sage: K = QQ['v,vl'].fraction_field(); v,vl=K.gens()
                    sage: H = AffineHeckeAlgebra(['C',3,1],q1=Family(dict([[0,vl],[1,v],[2,v],[3,vl]])))
                    sage: H.HY_KY().algebra_generators()
                    Finite family {0: Ty[1,2,3,2,1] Y[(-1, 0, 0)] + ((vl^2-1)/vl), 1: Ty[1], 2: Ty[2], 3: Ty[3]}
                """
                pass

            @abstract_method(optional=False)
            def Y_morphism(self, y):
                r"""
                The image of `y` under the morphism from the group algebra of `Y` into ``self``.

                EXAMPLES::

                    sage: H = AffineHeckeAlgebra("A2")
                    sage: AX = H.AX()
                    sage: KY = H.KY()
                    sage: y = KY.monomial(KY.basis().keys().simple_root(1)); y
                    Y[(1, -1, 0)]
                    sage: z = AX.Y_morphism(y); z
                    TX[0,2,0,1] + ((-v^2+1)/v)*TX[0,2,1]                   
                    sage: HY_KY = H.HY_KY()
                    sage: HY_KY(z)
                    Y[(1, -1, 0)]
                    sage: HY_KY(z) == HY_KY.Y_morphism(y)
                    True
                """
                pass

            @abstract_method(optional=False)
            def classical_hecke_morphism(self, a):
                r"""
                Returns the image of `a` from the finite Hecke algebra into ``self``.

                EXAMPLES::

                    sage: H = AffineHeckeAlgebra("A2")
                    sage: HY = H.HY()
                    sage: h = HY.an_element(); h
                    Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                    sage: H.AX().classical_hecke_morphism(h)
                    3*TX[2,1] + 3*TX[1,2] + TX[1,2,1]
                """
                pass

            def from_reduced_word(self, word):
                r"""
                Converts an affine or finite reduced word into a group element.

                .. warning::

                    Must be implemented in style "AX".

                EXAMPLES::

                    sage: AffineHeckeAlgebra("A2").HY_KY().from_reduced_word([0,2,1])
                    Ty[2] Y[(1, -1, 0)]

                """
                return self(self.realization_of().AX().from_reduced_word(word))

    class _Bases(UniqueRepresentation, BindableClass):
        r"""
        The class of realizations of an affine Hecke algebra.
        """

        def _repr_(self):
            r"""
            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").AX() # indirect doctest
                AX basis of The affine Hecke algebra of type ['A', 2, 1]

            """
            return "%s basis of the %s"%(self._prefix, self.realization_of())

        def is_parent_of(self, x):
            return x.parent() == self

    class AffineHeckeAlgebraAX(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "AX" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").AX()
            AX basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        def __init__(self, E):
            # the nonextended Hecke algebra of type `\tilde{X}`
            E._HXa = MultiParameterHeckeAlgebra(E.affine_weyl(), E.q1(), E.q2(), prefix="TX", category=AlgebrasWithBasis(E.base_ring()))
            # the group algebra of the fundamental group of type `\tilde{X}`
            E._KFX = E.fundamental_group().algebra(E.base_ring())
            E._KFX._print_options['prefix'] = ""
            mcat = ModulesWithBasis(E.base_ring()).TensorProducts()
            E._HXaoKFX = tensor([E._HXa, E._KFX], category = mcat)
            E._KFXoHXa = tensor([E._KFX, E._HXa], category = mcat)
            def ext_twist_func((w, f)):
                return (f, f.inverse().act_on_affine_weyl(w))
            E._ext_twist = E._HXaoKFX.module_morphism(on_basis=E._HXaoKFX.monomial*ext_twist_func, codomain=E._KFXoHXa, category=mcat)
            SmashProductAlgebra.__init__(self, E._KFX, E._HXa, E._ext_twist, category=Category.join((E._BasesCategory(),AlgebrasWithBasis(E.base_ring()).TensorProducts())))
            self._style = "AX"

        def _repr_(self):
            E = self.realization_of()
            return "%s basis of %s"%(self._style, E._repr_())

        @cached_method
        def algebra_generators(self):
            r"""
            The generators `T_i` of the affine Hecke algebra.
            """
            I = self.realization_of().cartan_type().index_set()
            AXa = self.factors()[1]
            return Family(dict([[i, self.factor_embedding(1)(AXa.algebra_generators()[i])] for i in I]))

        def from_reduced_word(self, word):
            r"""
            The basis element for a reduced word of affine type.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").AX().from_reduced_word([0,2,1])
                TX[0,2,1]                

            """
            H = self.realization_of()
            return self.factor_embedding(1)(self.factors()[1].monomial(H.affine_weyl().from_reduced_word(word)))

        def from_fundamental(self, f):
            r"""
            The basis element for an element of the fundamental group.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: AX = H.AX()
                sage: AX.from_fundamental(H.fundamental_group()(2))
                [piX[2]]
            """
            return self.factor_embedding(0)(self.factors()[0].monomial(f))

        def classical_hecke_morphism(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: HY = H.HY()
                sage: h = HY.an_element(); h
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                sage: H.AX().classical_hecke_morphism(h)
                3*TX[2,1] + 3*TX[1,2] + TX[1,2,1]

            """
            return self(a)

        def Y_morphism(self, y):
            r"""
            The image of the `Y`-lattice group algebra element `y` into ``self``.

            EXAMPLES::

                sage: H=AffineHeckeAlgebra("A2")
                sage: KY = H.KY()
                sage: z = H.AX().Y_morphism(KY.monomial(KY.basis().keys().simple_root(2))); z
                TX[0,1,0,2] + ((-v^2+1)/v)*TX[0,1,2]
                sage: H.KY_HY()(z)
                Y[(0, 1, -1)]
                sage: H.HY_KY()(z)
                Y[(0, 1, -1)]
                sage: H.KY_HY()(H.HY_KY()(z))
                Y[(0, 1, -1)]

            """
            return self(self.realization_of().HY_KY()(y))

    class AffineHeckeAlgebraHY_KY(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "HY_KY" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").HY_KY()
            HY_KY basis of The affine Hecke algebra of type ['A', 2, 1]
        """

        def __init__(self, E):
            KY = E.KY()
            H = E.HY()
            module_category = ModulesWithBasis(E.base_ring())
            self._YoHY = tensor([KY,H],category=module_category)
            self._HYoY = tensor([H,KY],category=module_category)
            self._YM = KY.nonreduced_demazure_lusztig_operators(E.q1(), E.q2(), convention="antidominant", doubled_parameters=E._doubled_parameters, side="right")
            def right_action_on_HY_KY((w,mu), i):
                smu = mu.simple_reflection(i)
                return tensor([H.monomial(w), self._YM[i](KY.monomial(mu)) - E.q1()[i]*KY.monomial(smu)]) + tensor([H.product_by_generator_on_basis(w,i), KY.monomial(smu)])

            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            self._HYoYM = HeckeAlgebraRepresentation(self._HYoY, right_action_on_HY_KY, E.cartan_type(), E.q1(), E.q2())

            def twist_func((mu,w)):
                return self._HYoYM.Tw(w)(self._HYoY.monomial((H.one_basis(), mu)))

            twist = self._YoHY.module_morphism(on_basis=twist_func,codomain=self._HYoY,category=module_category)
            SmashProductAlgebra.__init__(self, H, KY, twist, category=Category.join((E._BasesCategory(), AlgebrasWithBasis(E.base_ring()).TensorProducts())))
            self._style = "HY_KY"

        def _repr_(self):
            E = self.realization_of()
            return "%s basis of %s"%(self._style, E._repr_())

        def Y_morphism(self, y):
            r"""
            The image of `y` under the morphism from the group algebra of `Y` into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: KY = H.KY()
                sage: y = KY.an_element(); y
                Y[(2, 2, 3)]
                sage: H.HY_KY().Y_morphism(y)
                Y[(2, 2, 3)]

            """
            return self.factor_embedding(1)(y)

        def classical_hecke_morphism(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: HY = H.HY()
                sage: h = HY.an_element(); h
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                sage: H.HY_KY().classical_hecke_morphism(h)
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]

            """
            return self.factor_embedding(0)(a)

        def T0(self):
            r"""
            The operator `T_0`.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").HY_KY().T0()
                Ty[1,2,1] Y[(-1, 0, 1)] + ((v^2-1)/v)                
            """
            return self.realization_of().KY_HY().T0().to_opposite()

        @cached_method
        def algebra_generators(self):
            r"""
            The family of generators `T_i` in the given realization.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").HY_KY().algebra_generators()
                Finite family {0: Ty[1,2,1] Y[(-1, 0, 1)] + ((v^2-1)/v), 1: Ty[1], 2: Ty[2]}

            """
            return Family(dict([[i, self.T0() if i == 0 else self.factor_embedding(0)(self.factors()[0].algebra_generators()[i])] for i in self.realization_of().index_set()]))


    class AffineHeckeAlgebraKY_HY(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "KY_HY" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").KY_HY()
            KY_HY basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        def __init__(self, E):
            KY = E.KY()
            H = E.HY()
            module_category = ModulesWithBasis(E.base_ring())
            self._YoHY = tensor([KY,H],category=module_category)
            self._HYoY = tensor([H,KY],category=module_category)
            self._YM = KY.nonreduced_demazure_lusztig_operators(E.q1(), E.q2(), convention="antidominant", doubled_parameters=E._doubled_parameters, side="left")

            def left_action_on_KY_HY((mu,w), i):
                smu = mu.simple_reflection(i)
                return tensor([self._YM[i](KY.monomial(mu)) - E.q1()[i]*KY.monomial(smu), H.monomial(w)]) + tensor([KY.monomial(smu),H.product_by_generator_on_basis(w,i,side="left")])

            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            self._YoHYM = HeckeAlgebraRepresentation(self._YoHY, left_action_on_KY_HY, E.cartan_type(), E.q1(), E.q2(), side="left")

            def untwist_func((w,mu)):
                return self._YoHYM.Tw(w)(self._YoHY.monomial((mu, H.one_basis())))

            untwist = self._HYoY.module_morphism(on_basis=untwist_func,codomain=self._YoHY,category=module_category)
            SmashProductAlgebra.__init__(self, KY, H, untwist, category=Category.join((E._BasesCategory(), AlgebrasWithBasis(E.base_ring()).TensorProducts())))

            self._style = "KY_HY"

        def _repr_(self):
            E = self.realization_of()
            return "%s basis of %s"%(self._style, E._repr_())

        def Y_morphism(self, y):
            r"""
            The image of `y` under the morphism from the group algebra of `Y` into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: KY = H.KY()
                sage: y = KY.an_element(); y
                Y[(2, 2, 3)]
                sage: H.KY_HY().Y_morphism(y)
                Y[(2, 2, 3)]
            """
            return self.factor_embedding(0)(y)

        def classical_hecke_morphism(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: HY = H.HY()
                sage: h = HY.an_element(); h
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                sage: H.KY_HY().classical_hecke_morphism(h)
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
            phi = E._cartan_type_Y.root_system().coroot_lattice().highest_root().associated_coroot()
            s_phi = phi.associated_reflection()
            HY = self.factors()[1]
            return self.from_direct_product(self.factors()[0].monomial(phi.to_ambient()), HY.product_by_signed_generator_sequence(HY.one(), s_phi, [-1 for i in range(len(s_phi))]))

        @cached_method
        def algebra_generators(self):
            r"""
            The algebra generators `T_i`.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").KY_HY().algebra_generators()
                Finite family {0: Y[(1, 0, -1)] Ty[1,2,1] + ((-v^2+1)/v)*Y[(1, 0, -1)] Ty[1,2] + ((-v^2+1)/v)*Y[(1, 0, -1)] Ty[2,1] + ((v^4-2*v^2+1)/v^2)*Y[(1, 0, -1)] Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(1, 0, -1)] Ty[2] + ((-v^6+2*v^4-2*v^2+1)/v^3)*Y[(1, 0, -1)], 1: Ty[1], 2: Ty[2]}

            """
            return Family(dict([[i, self.T0() if i == 0 else self.factor_embedding(1)(self.factors()[1].algebra_generators()[i])] for i in self.realization_of().index_set()]))
