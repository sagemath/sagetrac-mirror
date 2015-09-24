# -*- coding: utf-8 -*-
r"""
Tools for expanding methods like product and coproduct on categories


AUTHOR:

- Jean-Baptiste Priez (first version)

"""
#*****************************************************************************
#  Copyright (C) 2015      Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.realizations import RealizationsCategory
from sage.categories.with_realizations import WithRealizationsCategory
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
import itertools


def expand_unary_operator(category, operator):
    """
    This function expands the unary operator ``operator`` along the category ``category`` (of at least modules).

    When we define an (unary) operator ``foo`` in a Sage category ``C``, we generally have to define
    the methods several times:

    * one method ``foo`` in ``ParentMethod`` (with the definition). By Sage convention, it is a good praxis
      to define methods in parent (and an element don't know how to do... Hey Daddy, could you do that for me??)
    * an other one in ``ElementMethod`` (with the same (boring) definition). The user should define methods in parent
      but at the end, the method is use by an element. The good praxis consists in ask to the parent.

    Most of the time a math structure (as an algebra) in Sage is defined as a structure *with realizations*
    and *with basis*, so our method also should be expand in the categories ``WithRealizations``,
    ``Realizations.ParentMethods``.

    In ``WithRealizations.ParentMethods`` the method ``foo` is in charge of find the *good* realization where
    the operator is defined. By default the method ``foo`` (of a realization) is not implemented.
    So by default the realization should ask (by the method ``foo_by_coercion``) if there exists an other
    realizations which implements the method.

    Most of the time, the operator ``foo`` is defined *on basis*, so the operation should be extended by linearity,
    one has:

    * a lazy method ``foo`` in the category ``WithBasis.ParentMethods`` which determine if the method ``foo_on_basis``
      is implemented, if it is then the method ``foo_on_basis`` is used and extended by linearity, else it uses
      the method by coercion,
    * a default not implemented method ``foo_on_basis`` (the user is free to implement in some bases),
    * a generic method ``foo_by_linearity`` to extend the operation by linearity.

    To resume, the minimalist code:

    .. code-block:: python

        class CategoryOfSomeWonderfulStructures(Category):

            def super_categories(self):
                return [...]

            class ParentMethods:

                def foo(self, I):
                    '''
                    My amazing unary operator with definition
                    '''

        expand_unary_operator(CategoryOfSomeWonderfulStructures, CategoryOfSomeWonderfulStructures.ParentMethods.foo)

    provides by the function ``expand_unary_operator`` the following code:

    .. code-block:: python

        class CategoryOfSomeWonderfulStructures(Category):

            def super_categories(self):
                return [...]

            class ParentMethods:

                def foo(self, I):
                    '''
                    My amazing unary operator with definition
                    '''

            class ElementMethods:

                def foo(self):
                    '''
                    My amazing unary operator applied
                    '''
                    return self.parent().foo(self)

            class WithRealizations(WithRealizationsCategory):
                class ParentMethods:

                    def foo(self, I):
                        return self.a_realization().foo(I)

            class Realizations(RealizationsCategory):
                class ParentMethods:

                    def foo_by_coercion(self, I):
                        ... looking for the good realizations which knows how to do


            class WithBasis(CategoryWithAxiom_over_base_ring):
                class ParentMethods:

                    @abstract_method
                    def foo_on_basis(self, i):
                        pass

                    def foo_by_linearity(self, I):
                        return self.linear_combination(...)

    TESTS::

        sage: from sage.categories.misc.operators_tools import expand_unary_operator
        sage: class A: pass
        sage: setattr(Algebras, "ParentMethods", A)
        sage: setattr(Algebras.ParentMethods, "mon_coproduct", lambda self, x: None)
        sage: expand_unary_operator(Algebras, "mon_coproduct")
        sage: hasattr(Algebras.ElementMethods, "mon_coproduct")
        True
        sage: hasattr(Algebras.WithRealizations.ParentMethods, "mon_coproduct")
        True
        sage: hasattr(Algebras.Realizations.ParentMethods, "mon_coproduct_by_coercion")
        True
        sage: hasattr(Algebras.WithBasis.ParentMethods, "mon_coproduct")
        True
        sage: hasattr(Algebras.WithBasis.ParentMethods, "mon_coproduct_on_basis")
        True
        sage: hasattr(Algebras.WithBasis.ParentMethods, "mon_coproduct_by_linearity")
        True

        sage: class PM:
        ....:    def foo(self, a):
        ....:       '''documentation'''
        sage: setattr(HopfAlgebras, "ParentMethods", PM)
        sage: expand_unary_operator(HopfAlgebras, "foo")

        sage: S = NonCommutativeSymmetricFunctions(QQ).S()
        sage: print S.foo.__doc__
        documentation
        sage: print S[3].foo.__doc__
        documentation
        sage: print S.realization_of().foo.__doc__
        documentation

        sage: S.foo_on_basis = lambda I: tensor([S.monomial(I), S.monomial(I)])
        sage: S[3,1,2].foo()
        S[3, 1, 2] # S[3, 1, 2]

        sage: R = NonCommutativeSymmetricFunctions(QQ).R()
        sage: R[3,2].foo()
        R[3, 2] # R[3, 2] + R[3, 2] # R[5] + R[5] # R[3, 2]

        sage: print R.foo.__doc__
        documentation

    """
    doc = getattr(category.ParentMethods, operator).__doc__
    # the method has to be defined on ParentMethods for documentation
    ##########################################################################
    # ElementMethods #########################################################
    ##########################################################################
    # class ElementMethods:
    #     def xxx_coproduct(self, x, y):
    #         return self.parent().xxx_coproduct(x, y)
    ################
    if not hasattr(category, "ElementMethods"):
        class EM:
            pass
        setattr(category, "ElementMethods", EM)
    EM = category.ElementMethods

    if not hasattr(EM, operator):
        setattr(EM, operator,
            lambda self: getattr(self.parent(), operator)(self)
        )
        setattr(getattr(EM, operator).im_func, "__doc__", doc)

    ##########################################################################
    # WithRealizations #######################################################
    ##########################################################################
    # class WithRealizations:
    #     class ParentMethods:
    #         def xxx_coproduct(self, x, y):
    #             return self.a_realization().xxx_coproduct(x)
    if not hasattr(category, "WithRealizations"):
        class WR(WithRealizationsCategory):
            class ParentMethods:
                pass
        setattr(category, "WithRealizations", WR)
    elif not hasattr(category.WithRealizations, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(category.WithRealizations, "ParentMethods", PM)
        except:
            category.WithRealizations.__dict__["ParentMethods"] = PM
    WRPM = category.WithRealizations.ParentMethods

    if not hasattr(WRPM, operator):
        setattr(WRPM, operator,
            lambda self, x: getattr(self.a_realization(), operator)(x)
        )
        setattr(getattr(WRPM, operator).im_func, "__doc__", doc)

    ##########################################################################
    # Realizations ###########################################################
    ##########################################################################
    # class Realizations:
    #     class ParentMethods:
    #         def xxx_coproduct_by_coercion(self, x):
    #             looking for the good realizations for
    #             the xxx_product
    if not hasattr(category, "Realizations"):
        class R(RealizationsCategory):
            pass
        setattr(category, "Realizations", R)
    if not hasattr(category.Realizations, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(category.Realizations, "ParentMethods", PM)
        except:
            category.Realizations.__dict__["ParentMethods"] = PM
    RPM = category.Realizations.ParentMethods

    if not hasattr(RPM, operator + "_by_coercion"):
        def meth_by_coercion(self, x):
            from sage.categories.tensor import tensor
            for R in self.realization_of().realizations():
                self_to_R = R.coerce_map_from(self)
                R_to_self = self.coerce_map_from(R)
                coprod_meth = getattr(R, operator)
                if self_to_R is not None and R_to_self is not None and \
                   coprod_meth != getattr(R, operator + "_by_coercion"):
                    SxS = tensor((self, self))
                    return SxS(coprod_meth(R(self(x))))
            return NotImplementedError

        setattr(RPM, operator + "_by_coercion", meth_by_coercion)
        setattr(getattr(RPM, operator + "_by_coercion").im_func, "__doc__", doc)

    ##########################################################################
    # WithBasis ##############################################################
    ##########################################################################
    # class WithBasis:
    #     class ParentMethods:
    #         @abstract_method
    #         def xxx_coproduct_on_basis(self, x):
    #             pass
    #
    #         def xxx_coproduct_by_linearity(self, X):
    #             return self.linear_combination(...)
    if not hasattr(category, "WithBasis"):
        class WB(CategoryWithAxiom_over_base_ring):
            pass
        setattr(category, "WithBasis", LazyImport(WB, ))
    if not hasattr(category.WithBasis, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(category.WithBasis, "ParentMethods", PM)
        except:
            category.WithBasis.__dict__["ParentMethods"] = PM
    WPM = category.WithBasis.ParentMethods
    ###############
    if not hasattr(WPM, operator + "_on_basis"):
        setattr(
            WPM, operator + "_on_basis",
            abstract_method(lambda self, x: NotImplemented, optional=True)
        )
    ################
    if not hasattr(WPM, operator + "_by_linearity"):
        def meth_by_linearity(self, x):
            return self.tensor_square().linear_combination(itertools.imap(
                # coproduct
                lambda (mon, coeff): \
                    (getattr(self, operator + "_on_basis")(mon), coeff),
                # iteration on the element
                x.monomial_coefficients().iteritems()
            ))
        setattr(WPM, operator + "_by_linearity", meth_by_linearity)
        setattr(getattr(WPM, operator + "_by_linearity").im_func, "__doc__", doc)

    ################
    # lazy attribute of the product
    if not hasattr(WPM, operator):
        def lazy_prod(self):
            if getattr(self, operator + "_on_basis") is NotImplemented:
                if hasattr(self, operator + "_by_coercion"):
                    return getattr(self, operator + "_by_coercion")
                else:
                    return NotImplemented
            return getattr(self, operator + "_by_linearity")

        setattr(WPM, operator, lazy_attribute(lazy_prod))


def expand_binary_operator(category, operator):
    """
    This function expand the binary operator ``operator`` along the category ``category``.

    For the documentation .. see:  :meth:``sage.categories.misc.operators_tools.expand_unary_operator``

    The concept of this function is similar to that of ``expand_unary_operator`` applied
    to binary operators.

    TESTS::

        sage: from sage.categories.misc.operators_tools import expand_binary_operator
        sage: class A: pass
        sage: setattr(Algebras, "ParentMethods", A)
        sage: setattr(Algebras.ParentMethods, "mon_product", lambda self, x, y: None)
        sage: expand_binary_operator(Algebras, "mon_product")
        sage: hasattr(Algebras.ElementMethods, "mon_product")
        True
        sage: hasattr(Algebras.WithRealizations.ParentMethods, "mon_product")
        True
        sage: hasattr(Algebras.Realizations.ParentMethods, "mon_product_by_coercion")
        True
        sage: hasattr(Algebras.WithBasis.ParentMethods, "mon_product")
        True
        sage: hasattr(Algebras.WithBasis.ParentMethods, "mon_product_on_basis")
        True
        sage: hasattr(Algebras.WithBasis.ParentMethods, "mon_product_by_linearity")
        True

        sage: class PM:
        ....:    def bar(self, a, b):
        ....:       '''documentation'''
        sage: setattr(HopfAlgebras, "ParentMethods", PM)
        sage: expand_binary_operator(HopfAlgebras, "bar")

        sage: S = NonCommutativeSymmetricFunctions(QQ).S()
        sage: print S.bar.__doc__
        documentation
        sage: print S[3].bar.__doc__
        documentation
        sage: print S.realization_of().bar.__doc__
        documentation

        sage: S.bar_on_basis = lambda I, J: S.monomial(I)
        sage: S[3,1,2].bar(S[3])
        S[3, 1, 2]
        sage: S[3,1,2].bar(S[3] + S[2,3,1] - 3*S[5])
        -S[3, 1, 2]

        sage: R = NonCommutativeSymmetricFunctions(QQ).R()
        sage: R[3].bar(R[3,1,2])
        0
        sage: R[3,1].bar(R[3])
        R[3, 1]

    """
    doc = getattr(category.ParentMethods, operator).__doc__
    # the method has to be defined on ParentMethods for documentation
    ##########################################################################
    # ElementMethods #########################################################
    ##########################################################################
    # class ElementMethods:
    #     def xxx_product(self, x, y):
    #         return self.parent().xxx_product(x, y)
    ################
    if not hasattr(category, "ElementMethods"):
        class EM:
            pass
        setattr(category, "ElementMethods", EM)
    EM = category.ElementMethods

    if not hasattr(EM, operator):
        setattr(EM, operator,
            lambda self, x: getattr(self.parent(), operator)(self, x)
        )
        setattr(getattr(EM, operator).im_func, "__doc__", doc)

    ##########################################################################
    # WithRealizations #######################################################
    ##########################################################################
    # class WithRealizations:
    #     class ParentMethods:
    #         def xxx_product(self, x, y):
    #             return self.a_realization().xxx_product(x, y)
    if not hasattr(category, "WithRealizations"):
        class WR(WithRealizationsCategory):
            class ParentMethods:
                pass
        setattr(category, "WithRealizations", WR)
    elif not hasattr(category.WithRealizations, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(category.WithRealizations, "ParentMethods", PM)
        except:
            category.WithRealizations.__dict__["ParentMethods"] = PM
    WRPM = category.WithRealizations.ParentMethods

    if not hasattr(WRPM, operator):
        setattr(WRPM, operator,
            lambda self, x, y: getattr(self.a_realization(), operator)(x, y)
        )
        setattr(getattr(WRPM, operator).im_func, "__doc__", doc)

    ##########################################################################
    # Realizations ###########################################################
    ##########################################################################
    # class Realizations:
    #     class ParentMethods:
    #         def xxx_product_by_coercion(self, x, y):
    #             looking for the good realizations for
    #             the xxx_product
    if not hasattr(category, "Realizations"):
        class R(RealizationsCategory):
            pass
        setattr(category, "Realizations", R)
    if not hasattr(category.Realizations, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(category.Realizations, "ParentMethods", PM)
        except:
            category.Realizations.__dict__["ParentMethods"] = PM
    RPM = category.Realizations.ParentMethods

    if not hasattr(RPM, operator + "_by_coercion"):
        def meth_by_coercion(self, x, y):
            for R in self.realization_of().realizations():
                self_to_R = R.coerce_map_from(self)
                R_to_self = self.coerce_map_from(R)
                prod_meth = getattr(R, operator)
                if self_to_R is not None and R_to_self is not None and \
                   prod_meth != getattr(R, operator + "_by_coercion"):
                    return self(prod_meth(R(self(x)), R(self(y))))
            return NotImplementedError

        setattr(RPM, operator + "_by_coercion", meth_by_coercion)
        setattr(getattr(RPM, operator + "_by_coercion").im_func, "__doc__", doc)

    ##########################################################################
    # WithBasis ##############################################################
    ##########################################################################
    # class WithBasis:
    #     class ParentMethods:
    #         @abstract_method
    #         def xxx_product_on_basis(self, x, y):
    #             pass
    #
    #         def xxx_product_by_linearity(self, X, Y):
    #             return self.linear_combination(...)
    if not hasattr(category, "WithBasis"):
        class WB(CategoryWithAxiom_over_base_ring):
            pass
        setattr(category, "WithBasis", WB)
    if not hasattr(category.WithBasis, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(category.WithBasis, "ParentMethods", PM)
        except:
            category.WithBasis.__dict__["ParentMethods"] = PM
    WPM = category.WithBasis.ParentMethods
    ###############
    if not hasattr(RPM, operator + "_on_basis"):
        setattr(
            WPM, operator + "_on_basis",
            abstract_method(lambda self, x, y: NotImplemented, optional=True)
        )
    ################
    if not hasattr(WPM, operator + "_by_linearity"):
        def meth_by_linearity(self, x, y):
            if x.parent() != y.parent():
                return getattr(self, operator + "_by_coercion")(x, y)
            return self.linear_combination(itertools.imap(
                # product
                lambda ((mon_l, coeff_l), (mon_r, coeff_r)): \
                    (getattr(self, operator + "_on_basis")(mon_l, mon_r),
                     coeff_l * coeff_r),
                # iteration on left and right
                itertools.product(
                    x.monomial_coefficients().iteritems(),
                    y.monomial_coefficients().iteritems()))
            )

        setattr(WPM, operator + "_by_linearity", meth_by_linearity)
        setattr(getattr(WPM, operator + "_by_linearity").im_func, "__doc__", doc)

    ################
    # lazy attribute of the product
    if not hasattr(WPM, operator):
        def lazy_prod(self):
            if getattr(self, operator + "_on_basis") is NotImplemented:
                return getattr(self, operator + "_by_coercion")
            else:
                return getattr(self, operator + "_by_linearity")

        setattr(WPM, operator, lazy_attribute(lazy_prod))
