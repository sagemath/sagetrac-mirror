"""
Root system data for hyperbolic types
"""
#*****************************************************************************
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def construct_hyperbolic(rank, index):
    """
    Construct the hyperbolic type of the given rank from the table
    in [CCCMNNP]_.

    REFERENCES:

    .. [CCCMNNP] Lisa Carbone, Sjuvon Chung, Leigh Cobbs, Robert McRae,
       Debajyoti Nandi, Yursa Naqvi, and Diego Penta. *Classification of
       hyperbolic Dynkin diagrams, root lengths, and Weyl orbits*. (2009).
    """

from cartan_type import CartanType_crystallographic, CartanType_hyperbolic
from cartan_type import CartanType_lorentzian as TypeLorentzian
from dynkin_diagram import DynkinDiagram

class CartanType(CartanType_hyperbolic, CartanType_crystallographic):
    r"""
    Hyperbolic Cartan types in rank at least 3.

    A Cartan type is hyperbolic if it is indefinate type and every proper
    subdiagram of the Dynkin diagram is finite or affine type.

    We note that our labeling is consistant with Kac and our Dynkin diagrams
    use the opposite convention for arrows in [CCCMNNF]_.
    """
    def __init__(self, rank, index):
        """
        Initialize ``self``.
        """
        self.n = rank
        self.index = index

    def _repr_(self, compact=False):
        r"""
        Return a string representation of ``self``.
        """
        if compact:
            return "H{}({})".format(self.index, self.n)
        return "['H', {}, {}]".format(self.index, self.n)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.
        """
        return "H_{{{}}}^{{({})}}".format(self.index, self.n)

class CartanType_lorentzian(CartanType, TypeLorentzian):
    """
    Concrete base class for hyperbolic Lorentzian Cartan types.
    """
    def __init__(self, affine, index):
        """
        Initialize ``self``.
        """
        TypeLorentzian.__init__(self, affine)
        CartanType.__init__(self, self._affine.rank()+1, index)

class CartanType_Rank3Cycle(CartanType):
    """
    Hyperbolic Cartan types of rank 3 that are a cycle.
    """
    def __init__(self, index):
        """
        Initialize ``self``.
        """
        if index == 1:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,1,2))
        elif index == 2:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,1,-22))
        elif index == 3:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,1,3))
        elif index == 4:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,1,4))
        elif index == 5:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((2,-2,1))
        elif index == 6:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,2,2))
        elif index == 7:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((-2,2,1))
        elif index == 8:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,2,-22))
        elif index == 9:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,2,3))
        elif index == 10:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,2,-3))
        elif index == 11:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,2,4))
        elif index == 12:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,2,-4))
        elif index == 13:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-2,-22))
        elif index == 14:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-2,3))
        elif index == 15:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-2,-3))
        elif index == 16:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-2,4))
        elif index == 17:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-2,-4))
        elif index == 18:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-22,-22))
        elif index == 19:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-22,3))
        elif index == 20:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-22,-3))
        elif index == 21:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-22,4))
        elif index == 22:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-22,-4))
        elif index == 23:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((3,-3,1))
        elif index == 24:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,3,3))
        elif index == 25:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((-3,3,1))
        elif index == 26:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,3,4))
        elif index == 27:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,3,-4))
        elif index == 28:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-3,4))
        elif index == 29:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-3,-4))
        elif index == 30:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,4,-4))
        elif index == 31:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,4,4))
        elif index == 32:
            dd = CartanType_Rank3Cycle._construct_dynkin_diagram((1,-4,4))
        else:
            raise ValueError("invalid index")
        self._dynkin_diagram = dd

    @staticmethod
    def _construct_dynkin_diagram(labels):
        """
        Construct the Dynkin diagram with edges labeled by ``labels``.

        For positive labels, the arrow points from `a` to `a+1` (mod 3),
        for negative labels, the arrow points from `a+1` to `a` (mod 3),
        and if the label is `-22`, there are arrows pointing in each
        direction (corresponding to type `A_1^{(1)}`).
        """
        g = DynkinDiagram()
        for i,l in enumerate(labels):
            i += 1 # +1 for indexing
            if l == -22:
                g.add_edge(i, (i+1)%3, 2)
                g.add_edge((i+1)%3, i, 2)
            elif l > 0:
                g.add_edge(i, (i+1)%3, l)
            elif l < 0:
                g.add_edge((i+1)%3, i, l)
        return g

    def dynkin_diagram(self):
        """
        Return the Dynkin diagram of ``self``.
        """
        return self._dynkin_diagram

    def index_set(self):
        """
        Return the index set of ``self``.
        """
        return (1, 2, 3)

