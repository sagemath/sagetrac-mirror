r"""
Bialgebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_types import Category_over_base_ring
from sage.categories.all import Algebras, Coalgebras
from sage.categories.realizations import RealizationsCategory

class Bialgebras(Category_over_base_ring):
    """
    The category of bialgebras

    EXAMPLES::

        sage: Bialgebras(ZZ)
        Category of bialgebras over Integer Ring
        sage: Bialgebras(ZZ).super_categories()
        [Category of algebras over Integer Ring, Category of coalgebras over Integer Ring]

    TESTS::

        sage: TestSuite(Bialgebras(ZZ)).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: Bialgebras(QQ).super_categories()
            [Category of algebras over Rational Field, Category of coalgebras over Rational Field]
        """
        R = self.base_ring()
        return [Algebras(R), Coalgebras(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        pass

    class Realizations(RealizationsCategory):

        class ParentMethods:

            def coproduct_by_coercion(self, x):
                r"""
                Try to returns the coproduct by coercion if
                `coproduct_by_basis` is not implemented.

                TESTS::

                    sage: Sym = SymmetricFunctions(QQ)
                    sage: m = Sym.monomial()
                    sage: f = m[2,1]
                    sage: m.coproduct_on_basis
                    NotImplemented
                    sage: m.coproduct == m.coproduct_by_coercion
                    True
                    sage: f.coproduct()
                    m[] # m[2, 1] + m[1] # m[2] + m[2] # m[1] + m[2, 1] # m[]

                TESTS::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = N.ribbon()
                    sage: R.coproduct_by_coercion.__module__
                    'sage.categories.bialgebras'
                    sage: R.coproduct_on_basis
                    NotImplemented
                    sage: R.coproduct == R.coproduct_by_coercion
                    True
                    sage: R[1].coproduct()
                    R[] # R[1] + R[1] # R[]
                """
                from sage.categories.tensor import tensor
                for R in self.realization_of().realizations():
                    if R.coerce_map_from(self) is not None and \
                                    self.coerce_map_from(R) is not None and \
                                    R.coproduct != R.coproduct_by_coercion:
                    #                     return self.tensor_square()(
                    #                         R(self(x)).coproduct()
                    #                     )

                        return self.tensor_square().sum(map(
                            lambda ((I, J), coeff): coeff * tensor((
                                self(R(I)),
                                self(R(J))
                            )), R(self(x)).coproduct(). \
                                monomial_coefficients().iteritems()
                        ))
                return NotImplementedError

            def product_by_coercion(self, left, right):
                r"""
                This method try to find a realization to coerce and compute the
                product and coerce back.
                The realizations which are try are specified by
                ``self.realization_of().realizations()``.

                TESTS::

                    sage: nM = NonCommutativeSymmetricFunctions(QQ).monomial()
                    sage: nM.product_on_basis
                    NotImplemented
                    sage: nM.product == nM.product_by_coercion
                    True
                    sage: nM([1,1])*nM([2])
                    3*nM[1, 1, 2] + nM[1, 3] + nM[2, 2]
                """
                for R in self.realization_of().realizations():
                    self_to_R = R.coerce_map_from(self)
                    R_to_self = self.coerce_map_from(R)
                    if self_to_R is not None and R_to_self is not None and \
                                    R.product != R.product_by_coercion:
                        return self(R(self(left)) * R(self(right)))
                return NotImplementedError
