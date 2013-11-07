r"""
Coinvariants of the Symmetric Group : quotient of multivariate polynomials by the symmetric polynomials
"""
#*****************************************************************************
#  Copyright (C) 2013 Nicolas Borie <nicolas.borie at math.u-psud dot fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.misc.abstract_method import abstract_method
import sage.libs.symmetrica.all as symmetrica

from sage.categories.algebras import Algebras
from sage.categories.modules import Modules
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.realizations import Realizations

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.combinat.free_module import CombinatorialFreeModule
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.combinat.permutation import Permutation
from sage.combinat.sf.sf import SymmetricFunctions

class Abstract_Quotient_MulPol_as_sym_module(UniqueRepresentation, Parent):
    r"""
    Returns the abstract quotient of multivariate polynomials by
    the ideal of constant free symmetric polynomials.

    EXAMPLES::


    TESTS::

        sage: A = AbstractPolynomialRingAsSymModule(QQ, 3)
        sage: S = A.schur_schubert()
        sage: H = A.Harmonic()
        sage: D = A.Descent()
        sage: M = A.staircase()
        sage: p = S.an_element(); print p; print p.to_quotient_by_symmetric_polynomials()
        s[]*Y_() + 2*s[]*Y_(2,3) + 3*s[]*Y_(1,2) + s[]*Y_(1,2,3)
        Y_() + 2*Y_(2,3) + 3*Y_(1,2) + Y_(1,2,3)
        sage: p = H.an_element(); print p; print p.to_quotient_by_symmetric_polynomials()
        s[]*H_() + 2*s[]*H_(2,3) + 3*s[]*H_(1,2) + s[]*H_(1,2,3)
        H_() + 2*H_(2,3) + 3*H_(1,2) + H_(1,2,3)
        sage: p = D.an_element(); print p; print p.to_quotient_by_symmetric_polynomials()
        s[]*(1) + 2*s[]*(x1*x3) + 3*s[]*(x2) + s[]*(x2*x3)
        (1) + 2*(x1*x3) + 3*(x2) + (x2*x3)
        sage: p = M.an_element(); print p; print p.to_quotient_by_symmetric_polynomials()
        s[]*(1) + 2*s[]*(x2) + 3*s[]*(x1) + s[]*(x1*x2)
        (1) + 2*(x2) + 3*(x1) + (x1*x2)
    """
    def __init__(self, AbstractModule):
        r"""
        TESTS::


        """
        R = AbstractModule.base_ring().base_ring()
        Parent.__init__(self, base=R, category = (Modules(R).WithRealizations(), Algebras(R)))
        self._abstract_module = AbstractModule
        self._main_repr_var = AbstractModule._main_repr_var

    def _repr_(self):
        r"""
        TESTS::

        """
        return "Abstract quotient of multivariate Polynomial in %s variables over %s by symmetric functions"%(self.nb_variables(), self.base_ring())

    @cached_method
    def nb_variables(self):
        r"""
        Return the number of variables of ``self``.

        EXAMPLES::

        """
        return self.ambient().nb_variables()

    @cached_method
    def ambient(self):
        r"""
        Return the ambient space from which the quotient ``self`` is
        defined.

        EXAMPLES::

        """
        return self._abstract_module

    def a_realization(self):
        r"""
        Return a realization of ``self``.

        EXAMPLES::

        """
        return self.Schubert()

    @cached_method
    def polynomial_algebra(self):
        r"""
        Returns the usual polynomial algebra of multivariate polynomials.

        EXAMPLES::

        """
        return self.ambient().polynomial_algebra()

    @abstract_method
    def lift(self, elt):
        r"""
        Return a lifted element in the ambient space built from ``elt``.

        All concrete realizations of this abstract quotient must
        overwrite this method.

        EXAMPLES::

        """
        pass

    @abstract_method
    def retract(self, elt):
        r"""
        Return a retracted element in ``self`` from an element ``elt``
        living in the ambient space of this quotient.

        All concrete realizations of this abstract quotient must
        overwrite this method.

        EXAMPLES::

        """
        pass

    class Bases(Category_realization_of_parent):
        r"""
        The category of bases as free module over symmetric
        polynomials of the multivariate polynomials.

        TESTS::


        """
        def super_categories(self):
            r"""
            """
            R = self.base().base_ring()
            return [ModulesWithBasis(R), Algebras(R), Realizations(self.base())]

        class ParentMethods:
            @cached_method
            def one_basis(self):
                r"""
                returns the key indexing the one for the
                multiplication in each realization.

                EXAMPLES::

                """
                return self.basis().keys().one()

            @cached_method
            def nb_variables(self):
                r"""
                Returns the number of variables of ``self``.

                EXAMPLES::

                """
                return self.category().base().nb_variables()

            @cached_method
            def polynomial_algebra(self):
                r"""
                Returns the usual polynomial algebra of multivariate polynomials.

                EXAMPLES::

                """
                return self.category().base().polynomial_algebra()

            def lift(self, elt):
                r"""
                Return a lifted element in the ambient space built from ``elt``.

                EXAMPLES::

                """
                SF = SymmetricFunctions(self.base_ring()).schur()
                d = dict()
                for key,coef in elt:
                    d[key] = SF(coef)
                return self.ambient()._from_dict(d)

            def retract(self, elt):
                r"""
                Returns a retracted element in ``self`` from an element ``elt``
                living in the ambient space of this quotient.

                EXAMPLES::

                """
                return elt.to_quotient_by_symmetric_polynomials()

        class ElementMethods:
            def to_multivariate_polynomial(self):
                r"""
                Returns ``self`` expressed in the polynomial algebra
                modelized by the parent of ``self``.

                EXAMPLES::


                """
                R = self.parent().category().base().polynomial_algebra()
                return R.sum([R(coef)*self.parent().basis_element_to_multivariate_polynomial(key) for key,coef in self])

            def lift_to_sym_module(self):
                r"""
                Returns the lifted element from ``self`` living in the
                corresponding realization of the abstract module.

                EXAMPLES::


                """
                return self.parent().lift(self)

    class Schubert(CombinatorialFreeModule, BindableClass):
        r"""
        The Schubert basis of the coinvariant of the symmetric group.

        EXAMPLES::


        TESTS::

        """
        def __init__(self, AbsQuo):
            r"""
            TESTS::

                sage: S = AbstractPolynomialRingAsSymModule(QQ, 3)
                sage: from sage.combinat.multivariate_polynomials.coinvariants_symmetric_group import Abstract_Quotient_MulPol_as_sym_module
                sage: T = Abstract_Quotient_MulPol_as_sym_module(S)
            """
            CombinatorialFreeModule.__init__(self, AbsQuo.base_ring(),
               SymmetricGroup(AbsQuo.nb_variables()), prefix='Y_', latex_prefix='Y',
               latex_bracket=False, bracket=False, category=AbsQuo.Bases())
            self._main_repr_var = AbsQuo._main_repr_var

        def _repr_(self):
            r"""
            TESTS::

                sage: S = AbstractPolynomialRingAsSymModule(QQ, 3).Schubert()
                sage: T = S.quotient_by_symmetric_polynomials(); T
                Coinvariants of the symmetric group of degree 3 in Schubert basis
            """
            return "Coinvariants of the symmetric group of degree %s in Schubert basis"%(self.nb_variables())

        def basis_element_to_multivariate_polynomial(self, key):
            r"""
            Returns the element basis of ``self`` indexed by the key
            ``key`` expressed in the polynomial algebra modelized by
            ``self``.

            EXAMPLES::

                sage: S = AbstractPolynomialRingAsSymModule(QQ, 3).Schubert()
                sage: T = S.quotient_by_symmetric_polynomials()
                sage: for k in T.basis().keys(): T.basis_element_to_multivariate_polynomial(k)
                1
                x1 + x2
                x1
                x1*x2
                x1^2
                x1^2*x2
            """
            return self.polynomial_algebra()(symmetrica.t_SCHUBERT_POLYNOM(list(Permutation(key))))

        def ambient(self):
            r"""
            Returns the ambient space from which ``self`` is a quotient.

            EXAMPLES::

                sage: S = AbstractPolynomialRingAsSymModule(QQ, 3).Schubert()
                sage: T = S.quotient_by_symmetric_polynomials()
                sage: T.ambient()
                Multivariate Polynomial in x1, x2, x3 over Symmetric Functions over Rational Field in the Schur basis in Schur Schubert decomposition
            """
            return self.realization_of().ambient().Schubert()

    S = schur_schubert = Schubert

    class Harmonic(CombinatorialFreeModule, BindableClass):
        r"""

        EXAMPLES::


        TESTS::


        """
        def __init__(self, AbsQuo):
            r"""
            TESTS::

            """
            CombinatorialFreeModule.__init__(self, AbsQuo.base_ring(),
               SymmetricGroup(AbsQuo.nb_variables()), prefix='H_', latex_prefix='H',
               latex_bracket=False, bracket=False, category=AbsQuo.Bases())
            self._main_repr_var = AbsQuo._main_repr_var

        def _repr_(self):
            r"""
            TESTS::

            """
            return "Coinvariants of the symmetric group of degree %s in Harmonic basis"%(self.nb_variables())

        def basis_element_to_multivariate_polynomial(self, key):
            r"""
            Returns the element basis of ``self`` indexed by the key
            ``key`` expressed in the polynomial algebra modelized by
            ``self``.

            EXAMPLES::

            """
            from sage.combinat.multivariate_polynomials.abstract_sym_module import harmonic_polynomial_for_permutation_group_element
            R = self.polynomial_algebra()
            return harmonic_polynomial_for_permutation_group_element(R, key)

        def _repr_term(self, m):
            r"""
            Returns a string representation of an Harmonic monomial
            indexed by a permutation group element.

            TESTS::


            """
            return 'H_' + str(m)

        def ambient(self):
            r"""
            Returns the ambient space from which ``self`` is a quotient.

            EXAMPLES::

            """
            return self.realization_of().ambient().Harmonic()

    H = harmonic = Harmonic

    class Descent(CombinatorialFreeModule, BindableClass):
        r"""


        EXAMPLES::


        TESTS::


        """
        def __init__(self, AbsQuo):
            r"""
            TESTS::

            """
            CombinatorialFreeModule.__init__(self, AbsQuo.base_ring(),
               SymmetricGroup(AbsQuo.nb_variables()), prefix='D_', latex_prefix='D',
               latex_bracket=False, bracket=False, category=AbsQuo.Bases())
            self._main_repr_var = AbsQuo._main_repr_var

        def _repr_(self):
            r"""
            TESTS::


            """
            return "Coinvariants of the symmetric group of degree %s in descent monomial basis"%(self.nb_variables())

        def basis_element_to_multivariate_polynomial(self, key):
            r"""
            Returns the element basis of ``self`` indexed by the key
            ``key`` expressed in the polynomial algebra modelized by
            ``self``.

            EXAMPLES::


            """
            R = self.polynomial_algebra()
            mon = Permutation(key).descent_polynomial()
            return R(mon)

        def _repr_term(self, m):
            r"""
            Returns a string representation of a descent monomial
            indexed by a permutation group element.

            TESTS::


            """
            return '(' + str(self.basis_element_to_multivariate_polynomial(m)) + ')'

        def ambient(self):
            r"""
            Returns the ambient space from which ``self`` is a quotient.

            EXAMPLES::

            """
            return self.realization_of().ambient().Descent()

    D = descent = Descent

    class Staircase(CombinatorialFreeModule, BindableClass):
        r"""

        EXAMPLES::


        TESTS::


        """
        def __init__(self, AbsQuo):
            r"""
            TESTS::

            """
            CombinatorialFreeModule.__init__(self, AbsQuo.base_ring(),
               SymmetricGroup(AbsQuo.nb_variables()), prefix='', latex_prefix='',
               latex_bracket=False, bracket=False, category=AbsQuo.Bases())
            self._main_repr_var = AbsQuo._main_repr_var

        def _repr_(self):
            r"""
            TESTS::

            """
            return "Coinvariants of the symmetric group of degree %s in monomial under the staircase basis"%(self.nb_variables())

        def basis_element_to_multivariate_polynomial(self, key):
            r"""
            Returns the element basis of ``self`` indexed by the key
            ``key`` expressed in the polynomial algebra modelized by
            ``self``.

            EXAMPLES::


            """
            R = self.polynomial_algebra()
            exp = Permutation(key).to_lehmer_code()
            return R.prod([R.gens()[i]**exp[i] for i in range(self.nb_variables())])

        def _repr_term(self, m):
            r"""
            Returns a string representation of a monomial under the
            staircase indexed by a permutation group element.

            TESTS::

            """
            return '(' + str(self.basis_element_to_multivariate_polynomial(m)) + ')'

        def ambient(self):
            r"""
            Returns the ambient space from which ``self`` is a quotient.

            EXAMPLES::

            """
            return self.realization_of().ambient().Staircase()

    M = staircase = Staircase
