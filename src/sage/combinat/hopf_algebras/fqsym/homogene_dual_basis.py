# -*- coding: utf-8 -*-
r"""
The homogeneous dual basis of FQSym Hopf algebra.

m-basis of FQSym
"""
#*****************************************************************************
#       Copyright (C) 2013 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.fqsym import FreeQuasiSymmetricFunctions
from sage.combinat.hopf_algebras import register_as_realization


class HomogeneDual(FreeQuasiSymmetricFunctions.Bases.Base):
    '''
    This basis `(m_\sigma)_{\mathfrak S}` defines an analogy to the basis of
    monomial quasi-symmetric functions (see
    :mod:`sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions`) for ``FQSym``.

    EXAMPLES::

        sage: m = FQSym(QQ).m(); m
        The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the HomogeneDual basis

    .. MATH::

        \mathbb{F}_{\sigma} = \sum_{\mu \\succeq_l \sigma} m_{\mu}

    where `\\succeq_l` is the weak order in the left permutohedron
    (see [AguSot]_). So by a Möbius inversion:

    EXAMPLES::

        sage: F = FQSym(QQ).F()
        sage: F(m[4,1,2,3])
        F[4, 1, 2, 3] - F[4, 1, 3, 2] - F[4, 2, 1, 3] + F[4, 3, 2, 1]
        sage: F(m[2,3,1])
        F[2, 3, 1] - F[3, 2, 1]

    TESTS::

        sage: m = FQSym(QQ).m(); m
        The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the HomogeneDual basis
        sage: TestSuite(m).run()
    '''
    _prefix = "m"

    def dual_basis(self):
        return self.realization_of().H()

    def build_morphisms(self):
        '''
        TESTS::

            sage: FQS = FQSym(QQ); F = FQS.F(); m = FQS.m()
            sage: F(m[2,1,3])
            F[2, 1, 3] - F[3, 1, 2]
            sage: m[1,2] * m[2,1]
            m[1, 2, 4, 3] + m[1, 3, 4, 2] + m[1, 4, 2, 3] + 3*m[1, 4, 3, 2] + m[2, 3, 4, 1] + 2*m[2, 4, 3, 1] + m[3, 4, 2, 1] + m[4, 1, 2, 3] + 2*m[4, 1, 3, 2] + m[4, 2, 3, 1] + m[4, 3, 1, 2]
        '''
        morph = lambda F, T, func, tri = None, comp = None: (
            F._module_morphism(func, codomain=T, triangular=tri, cmp=comp)
            if comp is not None else
            F._module_morphism(func, codomain=T, triangular=tri))

        F = self.realization_of().F()
        m = self

        # F to m and back
        F_to_m = morph(F, m,
            lambda sig: m.sum_of_monomials(
                sig.permutohedron_greater(side='left')),
            tri="lower")
        F_to_m.register_as_coercion()
        (~F_to_m).register_as_coercion()

register_as_realization(FreeQuasiSymmetricFunctions, HomogeneDual, "m")