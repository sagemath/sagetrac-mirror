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

from cartan_type import CartanType_simple, CartanType_crystallographic, CartanType_hyperbolic
class CartanType(CartanType_simple, CartanType_crystallographic, CartanType_hyperbolic):
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

