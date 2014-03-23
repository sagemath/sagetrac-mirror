r"""
Root system data for Lorentzian type E
"""

from sage.combinat.root_system.cartan_type import CartanType_simply_laced
from sage.combinat.root_system.type_hyperbolic import CartanType_lorentzian

class CartanType_E6Lorentzian(CartanType_simply_laced, CartanType_lorentzian):
    r"""
    The Cartan type `E_6^{(1) \wedge}` which is `T_{4,3,3}` or `H_5^{(8)}`.
    """
    def __init__(self):
        """
        Initialize ``self``.
        """
        CartanType_lorentzian.__init__(self, ['E',6,1], 4)

    def ascii_art(self, label=lambda x: x):
        r"""
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::
        """
        return ("        O {}\n        |\n        |\n"*3 +
                "O---O---O---O---O\n{}   {}   {}   {}   {}").format(
                *map(label, [-1,0,2,1,3,4,5,6]))

class CartanType_E7Lorentzian(CartanType_simply_laced, CartanType_lorentzian):
    r"""
    The Cartan type `E_7^{(1) \wedge}` which is `T_{5,4,2}` or `H_5^{(9)}`.
    """
    def __init__(self):
        """
        Initialize ``self``.
        """
        CartanType_lorentzian.__init__(self, ['E',7,1], 4)

    def ascii_art(self, label=lambda x: x):
        r"""
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::
        """
        return ("            O {}\n            |\n            |\n" +
                " O" + "---O"*7 + "\n{}" + "   {}"*7).format(
                *map(label, [2,-1,0,1,3,4,5,6,7]))

class CartanType_E10(CartanType_simply_laced, CartanType_lorentzian):
    r"""
    The Cartan type `E_{10}` which is `E_8^{(1)\wedge}` or `H_4^{(10)}`.
    """
    def __init__(self):
        """
        Initialize ``self``.
        """
        CartanType_lorentzian.__init__(self, ['E',8,1], 4)

    def _latex_dynkin_diagram(self, label=lambda x: x, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::
        """
        ret = "\\draw (0 cm,0) -- ({} cm,0);\n".format(8*node_dist)
        ret += "\\draw ({} cm, 0 cm) -- +(0,{} cm);\n".format(2*node_dist, node_dist)
        ret += "\\draw[fill=white] (0, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(label(1))
        for i in range(1, 7):
            ret += "\\draw[fill=white] ({} cm, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(i*node_dist, label(i+2))
        ret += "\\draw[fill=white] ({} cm, {} cm) circle (.25cm) node[right=3pt]{{${}$}};".format(2*node_dist, node_dist, label(2))
        ret += "\\draw[fill=white] ({} cm, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(7*node_dist, label(0))
        ret += "\\draw[fill=white] ({} cm, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(8*node_dist, label(-1))
        return ret

    def ascii_art(self, label=lambda x: x):
        r"""
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::
        """
        return ("        O {}\n        |\n        |\nO" + "---O"*8 + "\n"
               "{}" + "   {}"*8).format(*map(label, (2,1,3,4,5,6,7,8,0,-1)))

