# -*- coding: utf-8 -*-
r"""
The combinatorial Hopf algebra of Parking function
"""

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.lazy_import import LazyImport
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.combinat.parking_functions import ParkingFunctions

from sage.combinat.cha.tools.generic_basis import GenericBasis


class ParkingQuasiSymmetricFunctions(UniqueRepresentation, Parent):
    r'''

    '''
    def __init__(self,R):
        '''
        '''
        from sage.categories.rings import Rings
        assert(R in Rings()), '%s must be a ring' % R
        Parent.__init__(
            self, base=R,
            category=GradedHopfAlgebras(R).WithRealizations()
        )

    M = Monomial = LazyImport(
        "monomial_basis", "Monomial")#sage.combinat.cha._my_pqsym.
    S = MonomialDual = LazyImport(
        "monomial_dual_basis", "MonomialDual")#sage.combinat.cha._my_pqsym.
    L = LeftWeakOrder = LazyImport(
        "left_weak_order_basis","LeftWeakOrder")
    R = RightWeakOrder = LazyImport(
        "right_weak_order_basis","RightWeakOrder")
    GL = GreaterLeft = LazyImport(
        "greater_left_basis", "GreaterLeft")
    SR = SmallerRight = LazyImport(
        "smaller_right_basis", "SmallerRight")
    GR = GreaterRight = LazyImport(
        "greater_right_basis", "GreaterRight")
    SL = SmallerLeft = LazyImport(
        "smaller_left_basis", "SmallerLeft")
    
    def _repr_(self):
        return "The combinatorial Hopf algebra of Parking Quasi-Symmetric" + \
            " Functions over the %s" % self.base_ring()

    class Bases(Category_realization_of_parent):

        def super_categories(self):
            R = self.base().base_ring()
            return [GradedHopfAlgebras(R).Realizations(),
                    GradedHopfAlgebrasWithBasis(R)]

        class ParentMethods:

            def build_morphisms(self):
                '''
                Define morphisms associated to the current basis
                '''

            def __init_extra__(self):
                self.build_morphisms()

        class Base(GenericBasis):
            _prefix = "** TO DEFINE **"
            _basis_indices = ParkingFunctions()
