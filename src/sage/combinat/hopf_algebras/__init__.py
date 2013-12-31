# -*- coding: utf-8 -*-
r"""
Pattern of graded connexe Hopf algebras

Author:
-------

    - Jean-Baptiste Priez
"""
#*****************************************************************************
#       Copyright (C) 2013 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.categories.realizations import Category_realization_of_parent
from sage.misc.abstract_method import abstract_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.rings import Rings


class GenericGradedConnexeHopfAlgebras(UniqueRepresentation, Parent):

    the_category = lambda self, R: HopfAlgebras(R).Graded().Connected().WithRealizations()

    def __init__(self, R):
        """
        """
        # TODO:: tester
        __import__(self.__module__ + ".bases")

        assert(R in Rings()), '%s must be a ring' % R
        Parent.__init__(self, base=R, category=self.the_category(R))

    _external_realizations = []

    def __init_extra__(self):
        for realization_name in self._external_realizations:
            getattr(self, realization_name)

    class Bases(Category_realization_of_parent):

        def super_categories(self):
            R = self.base().base_ring()
            return [GradedHopfAlgebrasWithBasis(R).Realizations()]

def register_as_realization(GCHopfAlgebra_class, realization_class, short_name=None):
    setattr(GCHopfAlgebra_class, realization_class.__name__, realization_class)
    GCHopfAlgebra_class._external_realizations.append(realization_class.__name__)
    if short_name:
        setattr(GCHopfAlgebra_class, short_name, realization_class)

from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.bindable_class import BindableClass

class GenericBasisOfGCHopfAlgebra(CombinatorialFreeModule, BindableClass):
    """
    Generic basis of a Combinatorial Hopf Algebra
    """

    _prefix = "** TO DEFINE **"
    _basis_indices = "** TO DEFINE **"

    def build_morphisms(self):
        """
        Has to be overload with all morphisms registrations
        """

    def __init__(self, CHA):
        #setattr(CHA, "Bases", GenericGradedConnexeHopfAlgebras.Bases)
        CombinatorialFreeModule.__init__(
            self,
            CHA.base_ring(),
            self._basis_indices,
            prefix=self._prefix,
            bracket=False,
            category=CHA.Bases())

        self.build_morphisms()

    def _repr_(self):
        return repr(self.realization_of()) + " on the " + \
            self._realization_name() + " basis"

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

def words_like_getitem(self, c, *rest):
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