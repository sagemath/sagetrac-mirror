r"""
Examples of twisted group algebras
"""
#*****************************************************************************
#       Copyright (C) 2015 Mark Shimozono <mshimo@math.vt.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __builtin__ import staticmethod
from sage.algebras.twisted_group_algebra import TwistedGroupAlgebra
from sage.combinat.root_system.cartan_type import CartanType
from sage.rings.integer_ring import ZZ
from sage.sets.family import Family

class MyTwistedGroupAlgebra(TwistedGroupAlgebra):
    r"""
    An example of a twisted group algebra.

    In this example the group is the Weyl group of a root system
    and the ring is the fraction field of the integer group algebra
    of the weight lattice.
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type):
        cartan_type = CartanType(cartan_type)
        return super(MyTwistedGroupAlgebra, cls).__classcall__(cls, cartan_type)

    def __init__(self, cartan_type):
        r"""
        EXAMPLES::

            sage: from sage.categories.examples.twisted_group_algebra import MyTwistedGroupAlgebra
            sage: MyTwistedGroupAlgebra("A2")
            Twisted group ring of Fraction Field of Group algebra of the Ambient space of the Root system of type ['A', 2] over Integer Ring and Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
        """
        Phi = cartan_type.root_system()
        if not hasattr(Phi, "ambient_space"):
            raise TypeError("Root system of type %s does not have an ambient space realization"%cartan_type)
        X = Phi.ambient_space()
        W = X.weyl_group(prefix="s")
        ZX = X.algebra(ZZ)
        QX = ZX.rational_function_field(prefix="x")

        self._my_weyl_autos = Family(W, lambda w: QX.induced_endomorphism([w.action(X.basis()[i]) for i in X.basis().keys()]))
        self._my_action = lambda w, a: self._my_weyl_autos[w](a)
        TwistedGroupAlgebra.__init__(self, QX, W, self._my_action)

