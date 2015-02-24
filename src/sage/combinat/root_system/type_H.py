"""
Root system data for type H
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from cartan_type import CartanType_standard_finite, CartanType_simple
from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import UniversalCyclotomicField


class CartanType(CartanType_standard_finite, CartanType_simple):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: ct = CartanType(['H',3])
            sage: ct
            ['H', 3]
            sage: ct._repr_(compact = True)
            'H3'
            sage: ct.rank()
            3

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
        assert n in [3, 4]
        CartanType_standard_finite.__init__(self, "H", n)

    def coxeter_diagram(self):
        """
        Return a Coxeter diagram for type H.

        EXAMPLES::

             sage: ct = CartanType(['H',3])
             sage: ct.coxeter_diagram()
             Graph on 3 vertices
             sage: sorted(ct.coxeter_diagram().edges())
             [(1, 2, 3), (2, 3, 5)]
             sage: ct.coxeter_matrix()
             [1 3 2]
             [3 1 5]
             [2 5 1]

             sage: ct = CartanType(['H',4])
             sage: ct.coxeter_diagram()
             Graph on 4 vertices
             sage: sorted(ct.coxeter_diagram().edges())
             [(1, 2, 3), (2, 3, 3), (3, 4, 5)]
             sage: ct.coxeter_matrix()
             [1 3 2 2]
             [3 1 3 2]
             [2 3 1 5]
             [2 2 5 1]
        """
        from sage.graphs.graph import Graph
        n = self.n
        g = Graph(multiedges=False)
        for i in range(1, n):
            g.add_edge(i, i + 1, 3)
        g.set_edge_label(n - 1, n, 5)
        return g

    def coxeter_number(self):
        """
        Return the Coxeter number associated with ``self``.

        EXAMPLES::

            sage: CartanType(['H',3]).coxeter_number()
            10
            sage: CartanType(['H',4]).coxeter_number()
            30
        """
        if self.n == 3:
            return 10
        return 30

    def dynkin_diagram(self):
        """
        Return the Dynkin diagram of type `H`.

        EXAMPLES::

            sage: a = CartanType(['H',3]).dynkin_diagram()
            sage: a
            O---O=5=O
            1   2   3
            H3
            sage: sorted(a.edges())
            [(1, 2, 1), (2, 1, 1), (2, 3, -E(5)^2 - E(5)^3),
            (3, 2, -E(5)^2 - E(5)^3)]
        """
        E = UniversalCyclotomicField().gen
        from dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)
        g.add_edge(1, 2)
        g.add_edge(2, 3)
        g.add_edge(n - 1, n, -E(5) ** 2 - E(5) ** 3)
        g.add_edge(n, n - 1, -E(5) ** 2 - E(5) ** 3)
        return g

    def ascii_art(self, label=None):
        """
        Return an ascii art representation of the Dynkin diagram.

        EXAMPLES::

            sage: print CartanType(['H',3]).ascii_art()
            O---O=5=O
            1   2   3
        """
        if label is None:
            label = lambda x: x
        n = self.n
        ret = (n - 2) * "O---" + "O=5=O\n"
        ret += "   ".join("{}".format(label(i)) for i in self.index_set())
        return ret
