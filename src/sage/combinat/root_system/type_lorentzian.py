"""
Root system data for Lorentzian types
"""
#*****************************************************************************
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.root_system.cartan_type import CartanType_simply_laced, CartanType_lorentzian

#####################################################################
## XE_n

class CartanType_DEn(CartanType_simply_laced, CartanType_lorentzian):
    r"""
    The Cartan type `D_n^{(1)\wedge}`.

    These are also denoted by `DE_{n+2}` and `H_1^{(n+2)}`.
    """
    def _latex_dynkin_diagram(self, label = lambda x: x, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::
        """
        n = self.rank()
        ret = "\\draw (0 cm,0) -- ({} cm,0);\n".format((n-3)*node_dist)
        ret += "\\draw ({} cm, 0 cm) -- +(0,{} cm);\n".format(2*node_dist, node_dist)
        ret += "\\draw ({} cm, 0 cm) -- +(0,{} cm);\n".format((n-4)*node_dist, node_dist)

        ret += "\\draw[fill=white] (0, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(label(1))
        for i in range(1, n-3):
            ret += "\\draw[fill=white] ({} cm, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(i*node_dist, label(i+2))
        ret += "\\draw[fill=white] ({} cm, {} cm) circle (.25cm) node[right=3pt]{{${}$}};".format(2*node_dist, node_dist, label(2))
        ret += "\\draw[fill=white] ({} cm, {} cm) circle (.25cm) node[right=3pt]{{${}$}};".format((n-4)*node_dist, node_dist, label(n-1))
        ret += "\\draw[fill=white] ({} cm, 0) circle (.25cm) node[below=4pt]{{${}$}};".format((n-3)*node_dist, label(n))
        return ret

    def ascii_art(self, label = lambda x: x):
        r"""
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::
        """
        n = self.rank()
        if n == 6:
            return "    " + special_str + " %s\n    |\n    |\nO---O---O\n%s   |%s  %s\n    |\n    O %s"%tuple(label(i) for i in (4,1,2,3,0,-1))
        ret =  "        O {} "   + "    "*(n-7) + "O {}"
        ret += ("\n        |   " + "    "*(n-7) + "|")*2
        labels = [2,n-1,1,3,4] + range(5, n-1) + [n]
        return (ret + "\nO" + "---O"*(n-3) + "\n"
                "{}" + "   {}"*(n-3)).format(*map(label, labels))

class CartanType_XEn(CartanType_lorentzian):
    r"""
    The Cartan type `XE_n` for `X = B,C` and `n = 7,8,9,10`.

    These are also denoted by `H_i^{(n)}` where `i = 2,3` for `X = B,C`
    respectively.
    """
    def dual(self):
        """
        Types `EB_n` and `EC_n` are in duality:

        EXAMPLES::

            sage: CartanType(["C", 3]).dual()
            ['B', 3]
        """
        if self.letter == 'B':
            return CartanType_XEn('C', self.n)
        # otherwise self.letter == 'C'
        return CartanType_XEn('B', self.n)

    def _latex_dynkin_diagram(self, label = lambda x: x, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::
        """
        n = self.n
        ret = "\\draw (0 cm,0) -- ({} cm,0);\n".format((n-3)*node_dist)
        ret += "\\draw ({} cm, 0 cm) -- +(0,{} cm);\n".format(2*node_dist, node_dist)
        ret += "\\draw ({} cm, 0.1 cm) -- +({} cm,0);\n".format((n-3)*node_dist, node_dist)
        ret += "\\draw ({} cm, -0.1 cm) -- +({} cm,0);\n".format((n-3)*node_dist, node_dist)
        if self.letter == 'B':
            ret += self._latex_draw_arrow_tip((n-2.5)*node_dist+0.2, 0, 0)
        elif self.letter == 'C':
            ret += self._latex_draw_arrow_tip((n-2.5)*node_dist-0.2, 0, 180)

        ret += "\\draw[fill=white] (0, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(label(1))
        for i in range(1, n-1):
            ret += "\\draw[fill=white] ({} cm, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(i*node_dist, label(i+2))
        ret += "\\draw[fill=white] ({} cm, {} cm) circle (.25cm) node[right=3pt]{{${}$}};".format(2*node_dist, node_dist, label(2))
        return ret

    def ascii_art(self, label=lambda x: x):
        r"""
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::
        """
        ret = "        O {}\n        |\n        |\nO" + "---O"*(self.n-3)
        if self.letter == 'B':
            ret += "=>=O"
        elif self.letter == 'C':
            ret += "=<=O"
        return (ret + "\n{}" + "   {}"*(self.n-2)).format(*map(label, [2,1,3,4] + range(5, self.n+1)))

#####################################################################
## Type E

class CartanType_E6Lorentzian(CartanType_simply_laced, CartanType_lorentzian):
    r"""
    The Cartan type `E_7^{(1) \wedge}` which is `T_{4,3,3}` or `H_5^{(8)}`.
    """
    def __init__(self):
        """
        Initialize ``self``.
        """
        CartanType_lorentzian.__init__(self, ['E',6,1])

class CartanType_E7Lorentzian(CartanType_simply_laced, CartanType_lorentzian):
    r"""
    The Cartan type `E_7^{(1) \wedge}` which is `T_{5,4,2}` or `H_5^{(9)}`.
    """
    def __init__(self):
        """
        Initialize ``self``.
        """
        CartanType_lorentzian.__init__(self, ['E',7,1])

class CartanType_E10(CartanType_simply_laced, CartanType_lorentzian):
    r"""
    The Cartan type `E_{10}` which is `E_8^{(1)\wedge}` or `H_4^{(10)}`.
    """
    def __init__(self):
        """
        Initialize ``self``.
        """
        CartanType_lorentzian.__init__(self, ['E',8,1])

    def _latex_dynkin_diagram(self, label = lambda x: x, node_dist=2):
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

    def ascii_art(self, label = lambda x: x):
        r"""
        Return a ascii art representation of the extended Dynkin diagram.

        EXAMPLES::
        """
        return ("        O {}\n        |\n        |\nO" + "---O"*8 + "\n"
               "{}" + "   {}"*8).format(*map(label, (2,1,3,4,5,6,7,8,0,-1)))

