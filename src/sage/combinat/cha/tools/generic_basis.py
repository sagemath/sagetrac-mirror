# -*- coding: utf-8 -*-
r"""
Generic basis of a Combinatorial Hopf Algebra

This module implements a generic basis of combinatorial Hopf algebra.
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
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.bindable_class import BindableClass


class GenericBasis(CombinatorialFreeModule, BindableClass):
    """
    Generic basis of a Combinatorial Hopf Algebra
    """
    _prefix = "** TO DEFINE **"
    _basis_indices = "** TO DEFINE **"

    def __init__(self, CHA):
        CombinatorialFreeModule.__init__(
            self,
            CHA.base_ring(),
            self._basis_indices,
            prefix=self._prefix,
            latex_prefix="\\mathbf{%s}" % self._prefix,
            bracket=False,
            category=CHA.Bases())

    def one_basis(self):
        """
        By default the `one` of a combinatorial Hopf algebra is the unique
        element of the homogeneous component of size `0`.
        """
        return self.basis().keys()([])

    def counit_on_basis(self, sigma):
        """
        By default the `counit` is given on basis by `1` (of the field) if
        `\sigma` is the `one` (of the module) and `0` (of the field) in
        otherwise.
        """
        return self.base_ring().one() if len(sigma) == 0 \
            else self.base_ring().zero()

    def _repr_(self):
        """
        TESTS::

            sage: FQSym(QQ).F()
            The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the Fundamental basis 
        """
        return repr(self.realization_of()) + " on the " + \
            self._realization_name() + " basis"

    def __getitem__(self, c, *rest):
        """
        This method implements the abuses of notations::

            F[1,3,2]
            F[[1,3,2]]
            F[FQSym.indices()([2,1])]

        .. todo::

            This should call ``super.monomial`` if the input can't
            be made into a composition so as not to interfere with
            the standard notation ``Psi['x,y,z']``.

            This could possibly be shared with Sym, FQSym, and
            other algebras with bases indexed by list-like objects

            TODO:: généraliser la méthode...
        """
        from sage.rings.integer import Integer
        Keys = self.basis().keys()
        if c in Keys:
            assert len(rest) == 0
            c = Keys(c)
        else:
            if len(rest) > 0 or isinstance(c, (int, Integer)):
                c = Keys([c] + list(rest))
            else:
                c = Keys(list(c))
        return self.monomial(c)
