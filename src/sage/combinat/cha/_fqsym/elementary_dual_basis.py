# -*- coding: utf-8 -*-
r"""
The elementary dual basis of FQSym Hopf algebra.
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>,
#                          Rémi Maurice <maurice@univ-mlv.fr>.
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
from sage.combinat.cha.fqsym import FreeQuasiSymmetricFunctions


class ElementaryDual(FreeQuasiSymmetricFunctions.Bases.Base):
    '''
    This basis `(n_\sigma)_{\mathfrak S}` defines an analogy to the basis of
    monomial quasi-symmetric functions (see
    :mod:`sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions`) for ``FQSym``.

    EXAMPLES::

        sage: n = FQSym(QQ).n(); n
        The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the ElementaryDual basis

    .. MATH::

        \mathbb{F}_{\sigma} = \sum_{\mu \\preceq_l \sigma} n_{\mu}

    where `\\preceq_l` is the weak order in the left permutohedron
    (see [AguSot]_). So by a Möbius inversion:

    EXAMPLES::

        sage: F = FQSym(QQ).F()
        sage: F(n[4,1,2,3])
        -F[3, 1, 2, 4] + F[4, 1, 2, 3]
        sage: F(n[2,3,1])
        -F[1, 3, 2] + F[2, 3, 1]

    TESTS::

        sage: n = FQSym(QQ).n(); n
        The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the ElementaryDual basis
        sage: TestSuite(n).run()
    '''
    _prefix = "n"

    def dual_basis(self):
        return self.realization_of().E()

    def build_morphisms(self):
        '''
        TESTS::

            sage: FQS = FQSym(QQ); F = FQS.F(); n = FQS.n()
            sage: F(n[2,1,3])
            -F[1, 2, 3] + F[2, 1, 3]
            sage: n[2,1,3] * n[1]
            n[1, 3, 2, 4] + 2*n[2, 1, 3, 4] + n[2, 1, 4, 3] + 2*n[2, 3, 1, 4] + n[2, 4, 1, 3] + n[3, 2, 1, 4] + n[4, 2, 1, 3]
        '''
        morph = lambda F, T, func, tri = None, comp = None: (
            F._module_morphism(func, codomain=T, triangular=tri, cmp=comp)
            if comp is not None else
            F._module_morphism(func, codomain=T, triangular=tri))

        F = self.realization_of().F()
        n = self

        # n to F and back
        F_to_n = morph(F, n,
            lambda sig: n.sum_of_monomials(
                sig.permutohedron_smaller(side='left')),
            tri="upper")
        F_to_n.register_as_coercion()
        (~F_to_n).register_as_coercion()
