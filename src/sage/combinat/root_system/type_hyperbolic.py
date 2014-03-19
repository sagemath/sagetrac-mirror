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

from cartan_type import CartanType_crystallographic
class CartanType(CartanType_crystallographic):
    r"""
    Hyperbolic Cartan types.

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
        return "['H', {}, {}".format(self.index, self.n)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.
        """
        return "H_{{{}}}^{{({})}}".format(self.index, self.n)

    def index_set(self):
        """
        Return the index set of ``self``.
        """
        return tuple(range(1, self.n+1))

    def is_finite(self):
        r"""
        Return if ``self`` is finite type, which is ``False``.
        """
        return False

    def is_affine(self):
        r"""
        Return if ``self`` is affine type, which is ``False``.
        """
        return False

#####################################################################
## XE_n

class CartanType_XEn(CartanType):
    r"""
    The Cartan type `XE_n` where `X = B,C,D`.

    We have the following indices:

    Rank 7:

    - `DE_7 = H_1^{(7)}`
    - `BE_7 = H_2^{(7)}`
    - `CE_7 = H_3^{(7)}`

    Rank 8:

    - `DE_8 = H_1^{(8)}`
    - `BE_8 = H_2^{(8)}`
    - `CE_8 = H_3^{(8)}`

    Rank 9:

    - `DE_9 = H_1^{(9)}`
    - `BE_9 = H_2^{(9)}`
    - `CE_9 = H_3^{(9)}`

    Rank 10:

    - `DE_{10} = H_1^{(10)}`
    - `BE_{10} = H_2^{(10)}`
    - `CE_{10} = H_3^{(10)}`
    """
    def __init__(self, letter, rank):
        """
        Initialize ``self``.
        """
        index = 1
        if letter == 'B':
            index = 2
        elif letter == 'C':
            index = 3
        self.letter = letter
        CartanType.__init__(self, rank, index)

    def dynkin_diagram(self):
        r"""
        Returns a Dynkin diagram for type `E_{10}`.

        EXAMPLES::
        """
        from dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        g.add_edge(1, 3)
        g.add_edge(2, 4)
        for i in range(3, 8):
            g.add_edge(i, i+1)
        if self.letter == 'D':
            g.add_edge(9, 8)
            g.add_edge(10, 8)
        elif self.letter == 'B':
            g.add_edge(8, 9)
            g.add_edge(9, 10, 2)
        elif self.letter == 'C':
            g.add_edge(8, 9)
            g.add_edge(10, 9, 2)
        return g

    def _latex_dynkin_diagram(self, label = lambda x: x, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::
        """
        ret = "\\draw (0 cm,0) -- ({} cm,0);\n".format(9*node_dist)
        ret += "\\draw ({} cm, 0 cm) -- +(0,{} cm);\n".format(2*node_dist, node_dist)
        ret += "\\draw[fill=white] (0, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(label(1))
        for i in range(2, 9):
            ret += "\\draw[fill=white] ({} cm, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(i*node_dist, label(i+2))
        ret += "\\draw[fill=white] ({} cm, {} cm) circle (.25cm) node[right=3pt]{{${}$}};".format(2*node_dist, node_dist, label(2))
        return ret

    def ascii_art(self, label = lambda x: x):
        r"""
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::
        """
        if self.letter == 'D':
            ret =  "        O {}" + " "*8*(self.n-4) + "O {}\n"
            ret += ("        |   " + " "*8*(self.n-4) + "|")*2
            labels = [2,self.n-1,3,4] + range(5, self.n-1) + [self.n]
            return (ret + "\nO" + "---O"*(self.n-2) + "\n"
                    "{}" + "   {}"*(self.n-3)).format(*map(label, labels))

        ret = "        O {}\n        |\n        |\nO" + "---O"*(self.n-3)
        if self.letter == 'B':
            ret += "=>=O"
        elif self.letter == 'C':
            ret += "=<=O"
        return (ret + "\n{}" + "   {}"*7).format(*map(label, [2,1,3,4] + range(5, self.n+1)))

#####################################################################
## Rank 10

class CartanType_E10(CartanType):
    r"""
    The Cartan type `E_{10}` which is `H_4^{(10)}`.
    """
    def __init__(self):
        """
        Initialize ``self``.
        """
        CartanType.__init__(self, 10, 4)

    def dynkin_diagram(self):
        r"""
        Returns a Dynkin diagram for type `E_{10}`.

        EXAMPLES::
        """
        from dynkin_diagram import DynkinDiagram_class
        g = DynkinDiagram_class(self)
        g.add_edge(1, 3)
        g.add_edge(2, 4)
        for i in range(3, 10):
            g.add_edge(i, i+1)
        return g

    def _latex_dynkin_diagram(self, label = lambda x: x, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::
        """
        ret = "\\draw (0 cm,0) -- ({} cm,0);\n".format(9*node_dist)
        ret += "\\draw ({} cm, 0 cm) -- +(0,{} cm);\n".format(2*node_dist, node_dist)
        ret += "\\draw[fill=white] (0, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(label(1))
        for i in range(2, 9):
            ret += "\\draw[fill=white] ({} cm, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(i*node_dist, label(i+2))
        ret += "\\draw[fill=white] ({} cm, {} cm) circle (.25cm) node[right=3pt]{{${}$}};".format(2*node_dist, node_dist, label(2))
        return ret

    def ascii_art(self, label = lambda x: x):
        r"""
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::
        """
        return ("        O {}\n        |\n        |\nO" + "---O"*8 + "\n"
               "{}" + "   {}"*8).format(*map(label, (2,1,3,4,5,6,7,8,9,10)))

# For unpickling backward compatibility (Sage <= 4.1)

# For unpickling backward compatibility (Sage <= 4.1)
from sage.structure.sage_object import register_unpickle_override
