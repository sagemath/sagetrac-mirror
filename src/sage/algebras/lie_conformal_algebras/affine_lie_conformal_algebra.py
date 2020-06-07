"""
Affine Lie Conformal Algebra

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer import Integer
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra
from .graded_lie_conformal_algebra import GradedLieConformalAlgebra

class AffineLieConformalAlgebra(GradedLieConformalAlgebra):
    def __init__(self,R,ct,names=None,prefix='',bracket=''):
        r"""
        The current or affine Kac-Moody Lie conformal algebra.

        INPUT:

        - ``R`` -- a commutative Ring; the base ring for this Lie
          conformal algebra.
        - ``ct`` -- a ``str`` or a ``CartanType``; the Cartan Type for
          the corresponding finite dimensional Lie algebra. It must
          correspond to a simple finite dimensional Lie algebra.

        EXAMPLES::

                sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                sage: V
                The affine Lie conformal algebra of type ['A', 1] over Rational Field

                sage: V = AffineLieConformalAlgebra(QQ, CartanType(["A",2,1]))
                Traceback (most recent call last):
                ...
                ValueError: Only affine algebras of simple finite dimensionalLie algebras are implemented
        """
        if type(ct) is str:
            from sage.combinat.root_system.cartan_type import CartanType
            ct = CartanType(ct)
        if not ( ct.is_finite() and ct.is_irreducible ):
            raise ValueError("Only affine algebras of simple finite dimensional"
                "Lie algebras are implemented")
        hv = Integer(ct.dual_coxeter_number())
        g = LieAlgebra(R, cartan_type=ct)
        B = g.basis()
        S = B.keys()
        gdict = {}
        for k1 in S:
            for k2 in S:
                if S.rank(k2) <= S.rank(k1):
                    myb = B[k1].bracket(B[k2]).monomial_coefficients()
                    myf = R(2).inverse_of_unit()*R(hv).inverse_of_unit()\
                          *g.killing_form(B[k1],B[k2])
                    if myb or myf:
                        gdict[(k1,k2)] = {}
                        if myb:
                            gdict[(k1,k2)][0] = {(nk,0):myb[nk] for nk in
                                                 myb.keys()}
                        if myf:
                            gdict[(k1,k2)][1] = {('K',0):myf}

        weights = (1,)*B.cardinality()
        self._ct = ct
        GradedLieConformalAlgebra.__init__(self,
                    R, gdict, index_set=S,
                    central_elements=('K',), weights=weights,
                    names=names, prefix=prefix,bracket=bracket)

    def cartan_type(self):
        """
        The Cartan type of this Lie conformal algbera.

        EXAMPLES::

            sage: R = AffineLieConformalAlgebra(QQ, 'B3')
            sage: R
            The affine Lie conformal algebra of type ['B', 3] over Rational Field
            sage: R.cartan_type()
            ['B', 3]
        """
        return self._ct

    def _repr_(self):
        return "The affine Lie conformal algebra of "\
            "type {} over {}".format(self._ct,self.base_ring())


