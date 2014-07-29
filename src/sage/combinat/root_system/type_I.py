"""
Root system data for type I
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cartan_type import CartanType_standard_finite, CartanType_simple
from sage.rings.universal_cyclotomic_field.all import UniversalCyclotomicField

class CartanType(CartanType_standard_finite, CartanType_simple):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['I',5])
            sage: ct
            ['I', 5]
            sage: ct._repr_(compact = True)
            'I2(5)'
            sage: ct.rank()
            2
            sage: ct.index_set()
            [1, 2]

            sage: ct.is_irreducible()
            True
            sage: ct.is_finite()
            True
            sage: ct.is_affine()
            False
            sage: ct.is_crystallographic()
            False
            sage: ct.is_simply_laced()
            False

        TESTS::

            sage: TestSuite(ct).run()
        """
        assert n >= 1
        CartanType_standard_finite.__init__(self, "I", n)

    def _repr_(self, compact = False):
        """
        TESTS::

            sage: ct = CartanType(['I',5])
            sage: repr(ct)
            "['I', 5]"
            sage: ct._repr_(compact=True)
            'I2(5)'
        """
        format = 'I2(%s)' if compact else "['I', %s]"
        return format%self.n

    def dynkin_diagram(self):
        """
        Returns a Dynkin diagram for type G.

        EXAMPLES::

            sage: g = CartanType(['G',2]).dynkin_diagram()
            sage: g
              3
            O=<=O
            1   2
            G2
            sage: sorted(g.edges())
            [(1, 2, 1), (2, 1, 3)]
        """
        from dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        bond = self.n
        UCF = UniversalCyclotomicField()
        #E = UCF.gens()
        if bond % 2 == 0:
            label1 = 1
        else:
            label1 = E(2*bond)+E(2*bond)**-1
        g.add_edge(1,2,label1)
        label2 = -(2+E(bond)+E(bond)**-1)/label1
        g.add_edge(2,1,-label2)
        return g

    def ascii_art(self, label = lambda x: x):
        """
        Returns an ascii art representation of the Dynkin diagram

        EXAMPLES::

            sage: print CartanType(['I',5]).ascii_art(label = lambda x: x+2)
              X
            O=<=O
            3   4
        """
        return "  X\nO=<=O\n%s   %s"%tuple(label(i) for i in (1,2))

    def rank(self):
        """
        Type `I_p` is of rank 2

        EXAMPLES::
            sage: CartanType(['I', 5]).rank()
            2
        """
        return 2

    def index_set(self):
        """
        Type `I_p` is of rank 2

        EXAMPLES::
            sage: CartanType(['I', 5]).index_set()
            [1, 2]
        """
        return [1, 2]

    def coxeter_diagram(self):
        """
        Returns the Coxeter matrix for this type.

        EXAMPLES::

            sage: ct = CartanType(['I', 4])
            sage: ct.coxeter_diagram()
            Graph on 2 vertices
            sage: ct.coxeter_diagram().edges()
            [(1, 2, 4)]
            sage: ct.coxeter_matrix()
            [1 4]
            [4 1]
        """
        from sage.graphs.graph import Graph
        return Graph([[1,2,self.n]], multiedges=False)

    def coxeter_number(self):
        """
        Return the Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['I',3]).coxeter_number()
            3
            sage: CartanType(['I',12]).coxeter_number()
            12
        """
        return self.n
