# -*- coding: utf-8 -*-
r"""
Tools for expanding methods like product, coproduct on categories

AUTHORS:

 - Jean-Baptiste Priez (first version)
"""
#*****************************************************************************
#  Copyright (C) 2013      Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.misc.abstract_method import abstract_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.realizations import RealizationsCategory
from sage.categories.with_realizations import WithRealizationsCategory
from sage.categories.category_with_axiom import \
    CategoryWithAxiom_over_base_ring
import itertools

class OperationExpander:

    def __new__(cls, R):
        n = super(OperationExpander, cls).__new__(cls)
        print "****"
        return n

    _list_of_expandable_unary_operations = None
    _list_of_expandable_binary_operations = None

    class ElementMethods:
        pass

    class WithRealizations:
        class ParentMethods:
            pass

    class Realizations:
        class ParentMethods:
            pass

    class WithBasis:
        class ParentMethods:
            pass


def expand_unary_operator(cls, meth):
    """
    TESTS::

        sage: from sage.categories.operators_tools import expand_unary_operator
        sage: class A: pass
        sage: setattr(Algebras, "ParentMethods", A)
        sage: setattr(Algebras, "mon_coproduct", lambda self, x: None)
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
    """
    # the method has to be defined on ParentMethods for documentation
    ##########################################################################
    # ElementMethods #########################################################
    ##########################################################################
    # class ElementMethods:
    #     def xxx_coproduct(self, x, y):
    #         return self.parent().xxx_coproduct(x, y)
    ################
    if not hasattr(cls, "ElementMethods"):
        class EM:
            pass
        setattr(cls, "ElementMethods", EM)
    EM = cls.ElementMethods

    if not hasattr(EM, meth):
        setattr(EM, meth,
            lambda self: getattr(self.parent(), meth)(self)
        )

    ##########################################################################
    # WithRealizations #######################################################
    ##########################################################################
    # class WithRealizations:
    #     class ParentMethods:
    #         def xxx_coproduct(self, x, y):
    #             return self.a_realization().xxx_coproduct(x)
    if not hasattr(cls, "WithRealizations"):
        class WR(WithRealizationsCategory):
            class ParentMethods:
                pass
        setattr(cls, "WithRealizations", WR)
    elif not hasattr(cls.WithRealizations, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(cls.WithRealizations, "ParentMethods", PM)
        except:
            cls.WithRealizations.__dict__["ParentMethods"] = PM
    WRPM = cls.WithRealizations.ParentMethods

    if not hasattr(WRPM, meth):
        setattr(WRPM, meth,
            lambda self, x: getattr(self.a_realization(), meth)(x)
        )

    ##########################################################################
    # Realizations ###########################################################
    ##########################################################################
    # class Realizations:
    #     class ParentMethods:
    #         def xxx_coproduct_by_coercion(self, x):
    #             looking for the good realizations for
    #             the xxx_product
    if not hasattr(cls, "Realizations"):
        class R(RealizationsCategory):
            pass
        setattr(cls, "Realizations", R)
    if not hasattr(cls.Realizations, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(cls.Realizations, "ParentMethods", PM)
        except:
            cls.Realizations.__dict__["ParentMethods"] = PM
    RPM = cls.Realizations.ParentMethods

    if not hasattr(RPM, meth + "_by_coercion"):
        def meth_by_coercion(self, x):
            from sage.categories.tensor import tensor
            for R in self.realization_of().realizations():
                self_to_R = R.coerce_map_from(self)
                R_to_self = self.coerce_map_from(R)
                coprod_meth = getattr(R, meth)
                if self_to_R is not None and R_to_self is not None and \
                   coprod_meth != getattr(R, meth + "_by_coercion"):
                    SxS = tensor((self, self))
                    return SxS(coprod_meth(R(self(x))))
            return NotImplementedError

        setattr(RPM, meth + "_by_coercion", meth_by_coercion)

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
    if not hasattr(cls, "WithBasis"):
        class WB(CategoryWithAxiom_over_base_ring):
            pass
        setattr(cls, "WithBasis", LazyImport(WB, ))
    if not hasattr(cls.WithBasis, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(cls.WithBasis, "ParentMethods", PM)
        except:
            cls.WithBasis.__dict__["ParentMethods"] = PM
    WPM = cls.WithBasis.ParentMethods
    ###############
    if not hasattr(WPM, meth + "_on_basis"):
        setattr(
            WPM, meth + "_on_basis",
            abstract_method(lambda self, x: NotImplemented, optional=True)
        )
    ################
    if not hasattr(WPM, meth + "_by_linearity"):
        def meth_by_linearity(self, x):
            return self.tensor_square().linear_combination(itertools.imap(
                # coproduct
                lambda (mon, coeff): \
                    (getattr(self, meth + "_on_basis")(mon), coeff),
                # iteration on the element
                x.monomial_coefficients().iteritems()
            ))

        setattr(WPM, meth + "_by_linearity", meth_by_linearity)
    ################
    # lazy attribute of the product
    if not hasattr(WPM, meth):
        def lazy_prod(self):
            if getattr(self, meth + "_on_basis") is NotImplemented:
                if hasattr(self, meth + "_by_coercion"):
                    return getattr(self, meth + "_by_coercion")
                else:
                    return NotImplemented
            return getattr(self, meth + "_by_linearity")

        setattr(WPM, meth, lazy_attribute(lazy_prod))


def expand_binary_operator(cls, meth):
    """
    TESTS::

        sage: from sage.combinat.hopf_algebras.categories.tools import expand_binary_operator
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
    """
    # the method has to be defined on ParentMethods for documentation

    ##########################################################################
    # ElementMethods #########################################################
    ##########################################################################
    # class ElementMethods:
    #     def xxx_product(self, x, y):
    #         return self.parent().xxx_product(x, y)
    ################
    if not hasattr(cls, "ElementMethods"):
        class EM:
            pass
        setattr(cls, "ElementMethods", EM)
    EM = cls.ElementMethods

    if not hasattr(EM, meth):
        setattr(EM, meth,
            lambda self, x: getattr(self.parent(), meth)(self, x)
        )

    ##########################################################################
    # WithRealizations #######################################################
    ##########################################################################
    # class WithRealizations:
    #     class ParentMethods:
    #         def xxx_product(self, x, y):
    #             return self.a_realization().xxx_product(x, y)
    if not hasattr(cls, "WithRealizations"):
        class WR(WithRealizationsCategory):
            class ParentMethods:
                pass
        setattr(cls, "WithRealizations", WR)
    elif not hasattr(cls.WithRealizations, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(cls.WithRealizations, "ParentMethods", PM)
        except:
            cls.WithRealizations.__dict__["ParentMethods"] = PM
    WRPM = cls.WithRealizations.ParentMethods

    if not hasattr(WRPM, meth):
        setattr(WRPM, meth,
            lambda self, x, y: getattr(self.a_realization(), meth)(x, y)
        )

    ##########################################################################
    # Realizations ###########################################################
    ##########################################################################
    # class Realizations:
    #     class ParentMethods:
    #         def xxx_product_by_coercion(self, x, y):
    #             looking for the good realizations for
    #             the xxx_product
    if not hasattr(cls, "Realizations"):
        class R(RealizationsCategory):
            pass
        setattr(cls, "Realizations", R)
    if not hasattr(cls.Realizations, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(cls.Realizations, "ParentMethods", PM)
        except:
            cls.Realizations.__dict__["ParentMethods"] = PM
    RPM = cls.Realizations.ParentMethods

    if not hasattr(RPM, meth + "_by_coercion"):
        def meth_by_coercion(self, x, y):
            for R in self.realization_of().realizations():
                self_to_R = R.coerce_map_from(self)
                R_to_self = self.coerce_map_from(R)
                prod_meth = getattr(R, meth)
                if self_to_R is not None and R_to_self is not None and \
                   prod_meth != getattr(R, meth + "_by_coercion"):
                    return self(prod_meth(R(self(x)), R(self(y))))
            return NotImplementedError

        setattr(RPM, meth + "_by_coercion", meth_by_coercion)

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
    if not hasattr(cls, "WithBasis"):
        class WB(CategoryWithAxiom_over_base_ring):
            pass
        setattr(cls, "WithBasis", WB)
    if not hasattr(cls.WithBasis, "ParentMethods"):
        class PM:
            pass
        try:
            setattr(cls.WithBasis, "ParentMethods", PM)
        except:
            cls.WithBasis.__dict__["ParentMethods"] = PM
    WPM = cls.WithBasis.ParentMethods
    ###############
    if not hasattr(RPM, meth + "_on_basis"):
        setattr(
            WPM, meth + "_on_basis",
            abstract_method(lambda self, x, y: NotImplemented, optional=True)
        )
    ################
    if not hasattr(WPM, meth + "_by_linearity"):
        def meth_by_linearity(self, x, y):
            if x.parent() != y.parent():
                return getattr(self, meth + "_by_coercion")(x, y)
            return self.linear_combination(itertools.imap(
                # product
                lambda ((mon_l, coeff_l), (mon_r, coeff_r)): \
                    (getattr(self, meth + "_on_basis")(mon_l, mon_r),
                     coeff_l * coeff_r),
                # iteration on left and right
                itertools.product(
                    x.monomial_coefficients().iteritems(),
                    y.monomial_coefficients().iteritems()))
            )

        setattr(WPM, meth + "_by_linearity", meth_by_linearity)
    ################
    # lazy attribute of the product
    if not hasattr(WPM, meth):
        def lazy_prod(self):
            if getattr(self, meth + "_on_basis") is NotImplemented:
                return getattr(self, meth + "_by_coercion")
            else:
                return getattr(self, meth + "_by_linearity")

        setattr(WPM, meth, lazy_attribute(lazy_prod))
