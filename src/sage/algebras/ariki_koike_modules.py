r"""
Ariki-Koike Modules

AUTHORS:

- Travis Scrimshaw (2018-05-31): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2018 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, division
from six.moves import range

from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.misc.latex import latex
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.modules import Modules
from sage.categories.rings import Rings
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.integer_ring import ZZ
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.partition_tuple import PartitionTuples
from sage.combinat.permutation import Permutations
from sage.combinat.tableau_tuple import StandardTableauTuples

class ArikiKoikeModule(CombinatorialFreeModule):
    """
    Irreducible representation of the Ariki-Koike algebra.

    REFERENCES:

    - [Mathas2002]_
    """
    @staticmethod
    def __classcall_private__(cls, la, q=None, u=None, R=None):
        """
        Normalize input to ensure a unique representation.
        """
        la = PartitionTuples()(la)
        # TODO: Once Ariki-Koike algebras are merged, pass the construction
        #   parameters off to that.
        n = la.size()
        r = la.level()
        if u is None:
            if q is not None:
                R = q.parent()
            if R is None:
                R = PolynomialRing(ZZ, 'u', r)
                u = R.gens()
                if q is None:
                    R = LaurentPolynomialRing(R, 'q')
                    q = R.gen()
            else:
                u = PolynomialRing(ZZ, 'u', r).gens()
                if q is None:
                    q = 'q'
        else:
            if not isinstance(u, (list,tuple)):
                u = [u]*r
            if R is None:
                from sage.structure.element import get_coercion_model
                cm = get_coercion_model()
                if q is None:
                    R = cm.common_parent(*[val.parent() for val in u])
                    R = LaurentPolynomialRing(R, 'q')
                    q = R.gen()
                else:
                    R = cm.common_parent(q.parent(), *[val.parent() for val in u])
            elif q is None:
                q = 'q'
            u = [R(val) for val in u]
        if R not in Rings().Commutative():
            raise TypeError("base ring must be a commutative ring")
        q = R(q)
        u = tuple(u)
        return super(ArikiKoikeModule, cls).__classcall__(cls, la, q, u, R)

    def __init__(self, la, q, u, R):
        """
        Initialize ``self``.
        """
        self._shape = la
        self._q = q
        self._u = u
        self._Pn = Permutations(la.size())
        indices = StandardTableauTuples(la)
        cat = Modules(R).WithBasis()
        CombinatorialFreeModule.__init__(self, R, indices, category=cat,
                                         bracket=False, prefix='S')

    def _repr_(self):
        return "Ariki-Koike module of shape {} over {}".format(self._shape, self.base_ring())

    def _latex_(self):
        return "S^{{{}}}_{{{}}}".format(latex(self._shape), latex(self.base_ring()))

    def _test_relations(self, **options):
        """
        Test that the relations of the Ariki-Koike algebra are statisfied.

        EXAMPLES::

            sage: from sage.algebras.ariki_koike_modules import ArikiKoikeModule
            sage: q = ZZ['q'].fraction_field().gen()
            sage: S = ArikiKoikeModule([[2,1],[1]], q, [q^2+1,q-3], q.parent())
            sage: S._test_relations(elements=S.basis())
        """
        tester = self._tester(**options)
        r = self._shape.level()
        n = self._shape.size()
        q = self._q

        # Build the polynomial for testing the T0 action
        z = PolynomialRing(self.base_ring(), 'DUMMY').gen()
        T0_poly = -prod(z - val for val in self._u)
        def apply_T0_power(b, exp):
            for i in range(exp):
                b = b.T(0)
            return b

        for b in tester.some_elements():
            t0 = self.linear_combination( (apply_T0_power(b, exp), c)
                                          for exp,c in enumerate(T0_poly) )
            tester.assertEquals(t0, self.zero())

            tester.assertEquals(b.T(0).T(1).T(0).T(1), b.T(1).T(0).T(1).T(0))
            for i in range(1, n):
                tester.assertEquals(b.T(i).T(i), (q-1)*b.T(i) + q*b)
                tester.assertEquals(b.T(i).T(0), b.T(0).T(i))
                if i < n - 1:
                    tester.assertEquals(b.T(i).T(i+1).T(i), b.T(i+1).T(i).T(i+1))
                    for j in range(i+2, n):
                        tester.assertEquals(b.T(i).T(j), b.T(j).T(i))

    def _T_on_basis(self, i, t):
        """
        Return the action of `T_i` on the basis element indexed by
        the standard tableau tuple ``t``.
        """
        R = self.base_ring()
        if i == 0:
            c = t.cells_containing(1)[0]
            if len(c) == 2: # It is of level 1 and a regular tableau
                c = (0,) + c
            coeff = self._q**(c[2] - c[1]) * self._u[c[0]]
            return self._from_dict({t: R(coeff)},
                                   remove_zeros=False, coerce=False)

        c = t.cells_containing(i)[0]
        cp = t.cells_containing(i+1)[0]
        if len(c) == 2: # It is of level 1 and a regular tableau
            c = (0,) + c
            cp = (0,) + cp
        if c[0] == cp[0] and c[2] == cp[2]: # same column
            return self._from_dict({t: -R.one()},
                                   remove_zeros=False, coerce=False)
        elif c[0] == cp[0] and c[1] == cp[1]: # same row
            return self._from_dict({t: self._q},
                                   remove_zeros=False, coerce=False)
        else: # result is standard
            s = t.symmetric_group_action_on_entries(self._Pn.simple_reflection(i))
            assert s.parent() is t.parent()
            def res(cell):
                return self._q**(cell[2] - cell[1]) * self._u[cell[0]]
            # Note that the residue of i in t is given by the cell c
            #   and of i in s corresponds to cell cp because the
            #   corresponding action of the permutation on t.
            ct = (self._q - 1) * res(c) / (res(c) - res(cp))
            cs = (self._q * res(c) - res(cp)) / (res(c) - res(cp))
            return self._from_dict({t: R(ct), s: R(cs)},
                                   remove_zeros=False, coerce=False)

    class Element(CombinatorialFreeModule.Element):
        def T(self, i):
            """
            Return the action of `T_i` on ``self``.
            """
            if not self._monomial_coefficients: # action on 0 is 0
                return self
            P = self.parent()
            return P.linear_combination((P._T_on_basis(i, t), c) for t,c in self)

