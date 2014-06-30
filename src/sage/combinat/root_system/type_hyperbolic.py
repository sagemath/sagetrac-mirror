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

    .. [Messaoud] Hechmi Ben Massaoud. *Almost split real forms for hyperbolic
       Kac-Moody Lie algebras*.
       http://hal.archives-ouvertes.fr/docs/00/06/10/26/PDF/asrfhkmla_1_.pdf
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

#####################################################################

class CartanType_rank3_two_types(CartanType):
    r"""
    Concrete base class for hyperbolic Cartan types for rank 3 determined
    by 2 types.

    The valid types are (specified in the following order):

    - ``'A'`` -- a single edge ``---`` with letter `A`
    - ``'B'`` -- a type `B_2` edge ``=>=`` with letter `B`
    - ``'C'`` -- a type `C_2` edge ``=<=`` with letter `C`
    - ``'G'`` -- a type `G_2` edge ``=<=`` with label 3 and letter `G`
    - ``'G*'`` -- a type `G_2^{\vee}` edge ``=>=`` with label 3 and
      letter `G^{\vee}`
    - ``'BC'`` -- a type `BC_2^{(1)}` edge ``=<=`` with label 4 and
      letter `A^{(2)}`
    - ``'BC*'`` -- a type `BC_2^{(1)\dagger}` edge ``=>=`` with label 4 and
      letter `A^{(2)\dagger}`
    - ``'A~'`` -- a type `A^{(1)}_1` edge ``<=>`` and letter `\widetilde{A}`
    """
    def __init__(self, X, Y):
        """
        Initialize ``self``.
        """
        self._X = X
        self._Y = Y
        index = 1 # TODO: placeholder
        CartanType.__init__(self, 3, index)

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: ct = CartanType()
            sage: latex(ct)
        """
        ret = ""
        if self._X == 'G*':
            ret += 'G^{{{}}}'.format(self.global_options('dual_latex'))
        elif self._X == 'A~':
            ret += '\\widetilde{A}'
        elif self._X == 'BC':
            ret += 'A^{(2)}'
        elif self._X == 'BC*':
            ret += 'A^{(2)\\dagger}'
        else:
            ret += self._X

        if self._Y == 'G*':
            ret += 'G^{{{}}}'.format(self.global_options('dual_latex'))
        elif self._Y == 'A~':
            ret += '\\widetilde{A}'
        elif self._Y == 'BC':
            ret += 'A^{(2)}'
        elif self._Y == 'BC*':
            ret += 'A^{(2)\\dagger}'
        else:
            ret += self._Y

        return ret

    def dual(self):
        """
        EXAMPLES::

           sage: ct = CartanType().dual()
           sage: ct.dual()
        """
        if self._X == 'A':
            if self._Y == 'A~':
                return self
            return self.relabel({1:1, 2:3, 3:2})

        if self._X == self._Y:
            if self._X == 'A~':
                return self
            return self.relabel({1:3, 2:2, 3:1})

        if self._X == 'A~':
            X = 'A~'
        elif self._X == 'B':
            X = 'C'
        elif self._X == 'C':
            X = 'B'
        elif self._X[-1] == '*':
            X = self._X[:-1]
        else:
            X = self._X + '*'

        if self._Y == 'A~':
            Y = 'A~'
        elif self._Y == 'B':
            Y = 'C'
        elif self._Y == 'C':
            X = 'B'
        elif self._Y[-1] == '*':
            Y = self._Y[:-1]
        else:
            Y = self._Y + '*'

        return self.__class__(X, Y)

    def index_set(self):
        """
        Return the index set of ``self``.
        """
        return (1, 2, 3)

    def dynkin_diagram(self):
        """
        Return the hyperbolic Dynkin diagram.

        EXAMPLES::

            sage: a = CartanType().dynkin_diagram()
            sage: a
             +-------+
             |       |
             |       |
             O---O---O
             1   2   3
            sage: sorted(a.edges())
            [(0, 1, 1),
             (0, 3, 1),
             (1, 0, 1),
             (1, 2, 1),
             (2, 1, 1),
             (2, 3, 1),
             (3, 0, 1),
             (3, 2, 1)]
        """
        from dynkin_diagram import DynkinDiagram_class
        n = self.n
        g = DynkinDiagram_class(self)
        X = self._X
        Y = self._Y
        if X == 'A':
            g.add_edge(0, 1, 1)
        if X == 'B':
            g.add_edge(0, 1, 2)
        elif X == 'C':
            g.add_edge(1, 0, 2)
        elif X == 'G':
            g.add_edge(0, 1, 3)
        elif X == 'G*':
            g.add_edge(1, 0, 3)
        elif X == 'BC':
            g.add_edge(0, 1, 4)
        elif X == 'BC*':
            g.add_edge(1, 0, 4)
        elif X == 'A~':
            g.add_edge(1, 0, 2)
            g.add_edge(0, 1, 2)

        if Y == 'B':
            g.add_edge(1, 2, 2)
        elif Y == 'C':
            g.add_edge(2, 1, 2)
        elif Y == 'G':
            g.add_edge(1, 2, 3)
        elif Y == 'G*':
            g.add_edge(2, 1, 3)
        elif Y == 'BC':
            g.add_edge(1, 2, 4)
        elif Y == 'BC*':
            g.add_edge(2, 1, 4)
        elif Y == 'A~':
            g.add_edge(2, 1, 2)
            g.add_edge(1, 2, 2)

        return g

    def _latex_dynkin_diagram(self, label = lambda x: x, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print CartanType()._latex_dynkin_diagram()
        """
        X = self._X
        Y = self._Y

        if X == 'B':
            ret += "\\draw (0 cm, 0.1 cm) -- +({} cm,0);\n".format(node_dist)
            ret += "\\draw (0 cm, -0.1 cm) -- +({} cm,0);\n".format(node_dist)
            ret += self._latex_draw_arrow_tip(0.5*node_dist+0.2, 0, 0)
        elif X == 'C':
            ret += "\\draw (0 cm, 0.1 cm) -- +({} cm,0);\n".format(node_dist)
            ret += "\\draw (0 cm, -0.1 cm) -- +({} cm,0);\n".format(node_dist)
            ret += self._latex_draw_arrow_tip(0.5*node_dist-0.2, 0, 180)
        elif X == 'A~':
            ret += "\\draw (0 cm, 0.1 cm) -- +({} cm,0);\n".format(node_dist)
            ret += "\\draw (0 cm, -0.1 cm) -- +({} cm,0);\n".format(node_dist)
            ret += self._latex_draw_arrow_tip(0.33*node_dist-0.2, 0, 180)
            ret += self._latex_draw_arrow_tip(0.66*node_dist+0.2, 0, 0)
        elif X[0] == 'G':
            ret = "\\draw (0,0) -- ({} cm,0);\n".format(node_dist)
            ret += "\\draw (0, 0.15 cm) -- +({} cm,0);\n".format(node_dist)
            ret += "\\draw (0, -0.15 cm) -- +({} cm,0);\n".format(node_dist)
            if X[-1] == '*': # G*
                ret += self._latex_draw_arrow_tip(0.5*node_dist+0.2, 0, 0)
            else: # G
                ret += self._latex_draw_arrow_tip(0.5*node_dist-0.2, 0, 180)
        elif X.startswith('BC'):
            ret = "\\draw (0, 0.05 cm) -- +({} cm,0);\n".format(node_dist)
            ret += "\\draw (0, -0.05 cm) -- +({} cm,0);\n".format(node_dist)
            ret += "\\draw (0, 0.15 cm) -- +({} cm,0);\n".format(node_dist)
            ret += "\\draw (0, -0.15 cm) -- +({} cm,0);\n".format(node_dist)
            if X[-1] == '*': # BC*
                ret += self._latex_draw_arrow_tip(0.5*node_dist+0.2, 0, 0)
            else: # BC
                ret += self._latex_draw_arrow_tip(0.5*node_dist-0.2, 0, 180)

        if Y == 'B':
            ret += "\\draw ({0} cm, 0.1 cm) -- +({0} cm,0);\n".format(node_dist)
            ret += "\\draw ({0} cm, -0.1 cm) -- +({0} cm,0);\n".format(node_dist)
            ret += self._latex_draw_arrow_tip(1.5*node_dist+0.2, 0, 0)
        elif Y == 'C':
            ret += "\\draw ({0} cm, 0.1 cm) -- +({0} cm,0);\n".format(node_dist)
            ret += "\\draw ({0} cm, -0.1 cm) -- +({0} cm,0);\n".format(node_dist)
            ret += self._latex_draw_arrow_tip(1.5*node_dist-0.2, 0, 180)
        elif Y == 'A~':
            ret += "\\draw ({0} cm, 0.1 cm) -- +({0} cm,0);\n".format(node_dist)
            ret += "\\draw ({0} cm, -0.1 cm) -- +({0} cm,0);\n".format(node_dist)
            ret += self._latex_draw_arrow_tip(1.33*node_dist-0.2, 0, 180)
            ret += self._latex_draw_arrow_tip(1.66*node_dist+0.2, 0, 0)
        elif Y[0] == 'G':
            ret = "\\draw ({0} cm,0) -- ({0} cm,0);\n".format(node_dist)
            ret += "\\draw ({0} cm, 0.15 cm) -- +({0} cm,0);\n".format(node_dist)
            ret += "\\draw ({0} cm, -0.15 cm) -- +({0} cm,0);\n".format(node_dist)
            if Y[-1] == '*': # G*
                ret += self._latex_draw_arrow_tip(1.5*node_dist+0.2, 0, 0)
            else: # G
                ret += self._latex_draw_arrow_tip(1.5*node_dist-0.2, 0, 180)
        elif Y.startswith('BC'):
            ret = "\\draw ({0} cm, 0.05 cm) -- +({0} cm,0);\n".format(node_dist)
            ret += "\\draw ({0} cm, -0.05 cm) -- +({0} cm,0);\n".format(node_dist)
            ret += "\\draw ({0} cm, 0.15 cm) -- +({0} cm,0);\n".format(node_dist)
            ret += "\\draw ({0} cm, -0.15 cm) -- +({0} cm,0);\n".format(node_dist)
            if Y[-1] == '*': # BC*
                ret += self._latex_draw_arrow_tip(1.5*node_dist+0.2, 0, 0)
            else: # BC
                ret += self._latex_draw_arrow_tip(1.5*node_dist-0.2, 0, 180)

        for i in range(3):
            ret += "\\draw[fill=white] ({} cm, 0) circle (.25cm) node[below=4pt]{{${}$}};\n".format(
                    i*node_dist, label(i+1))
        return ret

    def ascii_art(self, label=lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        EXAMPLES::

            sage: print CartanType().ascii_art()
            +-------+
            |       |
            |       |
            O---O---O
            1   2   3

            sage: print CartanType().ascii_art(label = lambda x: x+2)
            +-------+
            |       |
            |       |
            O---O---O
            1   2   3
        """
        X = self._X
        Y = self._Y
        top_row = False
        ret = ''

        if X[0] == 'G':
            ret += '  3  '
            top_row = True
        elif X.startswith('BC'):
            ret += '  4  '
            top_row = True
        else:
            ret += '     '

        if Y[0] == 'G':
            ret += ' 3'
            top_row = True
        elif Y.startswith('BC'):
            ret += ' 4'
            top_row = True

        if not top_row:
            ret = 'O'
        else:
            ret += '\nO'

        if X == 'A':
            ret += "---O"
        elif X in ['C', 'BC', 'G']:
            ret += "=<=O"
        elif X == 'B' or Y[-1] == '*': # G* or BC*
            ret += "=>=O"
        elif X == 'A~':
            ret += "<=>O"

        if Y in ['C', 'BC', 'G']:
            ret += "=<=O"
        elif Y == 'B' or Y[-1] == '*': # G* or BC*
            ret += "=>=O"
        elif Y == 'A~':
            ret += "<=>O"

        ret += "\n{}   {}   {}".format(*map(label, self.index_set()))
        return ret

    def is_compact(self):
        """
        Return if ``self`` is a compact hyperbolic Cartan type.
        """
        return self._X not in ['BC', 'A~'] and self._Y not in ['BC', 'A~']

# These are all (non necessarily standard) Lorentzian extensions of affine types
class CartanType_rank3_XY(CartanType_rank3_two_types):
    """
    Hyperbolic Cartan types for rank 3 that are of the form `XY`.
    """
    def dynkin_diagram(self):
        """
        Return the hyperbolic Dynkin diagram.

        EXAMPLES::

            sage: a = CartanType().dynkin_diagram()
            sage: a
             +-------+
             |       |
             |       |
             O---O---O
             1   2   3
            sage: sorted(a.edges())
            [(0, 1, 1),
             (0, 3, 1),
             (1, 0, 1),
             (1, 2, 1),
             (2, 1, 1),
             (2, 3, 1),
             (3, 0, 1),
             (3, 2, 1)]
        """
        return CartanType_rank3_two_types.dynkin_diagram(self).copy(immutable=True)

class CartanType_rank3_AXY(CartanType_rank3_two_types):
    """
    Hyperbolic Cartan types for rank 3 that are of the form `AXY`.
    """
    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: ct = CartanType()
            sage: latex(ct)
        """
        return "A" + CartanType_rank3_two_types._latex_(self)

    def dynkin_diagram(self):
        """
        Return the hyperbolic Dynkin diagram.

        EXAMPLES::

            sage: a = CartanType().dynkin_diagram()
            sage: a
             +-------+
             |       |
             |       |
             O---O---O
             1   2   3
            sage: sorted(a.edges())
            [(0, 1, 1),
             (0, 3, 1),
             (1, 0, 1),
             (1, 2, 1),
             (2, 1, 1),
             (2, 3, 1),
             (3, 0, 1),
             (3, 2, 1)]
        """
        g = CartanType_rank3_two_types.dynkin_diagram(self)
        g.add_edge(0, 2, 1)
        return g.copy(immutable=True)

    def _latex_dynkin_diagram(self, label=lambda x: x, node_dist=2):
        r"""
        Return a latex representation of the Dynkin diagram.

        EXAMPLES::

            sage: print CartanType()._latex_dynkin_diagram()
        """
        X = self._X
        Y = self._Y

        # Connecting top line
        ret = "\\draw (0 cm,0 cm) -- (0 cm, 1 cm) -- ({0} cm,1 cm) -- ({0} cm, 0 cm);\n".format(2*node_dist)
        ret += "\\draw (0 cm, 1cm) -- ({} cm,1 cm);\n".format(2*node_dist)
        ret += CartanType_rank3_two_types._latex_dynkin_diagram(label, node_dist)
        return ret

    def ascii_art(self, label=lambda x: x):
        """
        Returns a ascii art representation of the extended Dynkin diagram

        EXAMPLES::

            sage: print CartanType().ascii_art()
            +-------+
            |       |
            |       |
            O---O---O
            1   2   3

            sage: print CartanType().ascii_art(label = lambda x: x+2)
            +-------+
            |       |
            |       |
            O---O---O
            1   2   3
        """
        X = self._X
        Y = self._Y
        ret = "+-------+\n|       |\n"
        if X[0] == 'G':
            ret += '| 3  '
        elif X.startswith('BC'):
            ret += '| 4  '
        else:
            ret += '|    '  
        if Y[0] == 'G':
            ret += ' 3 |'
        elif Y.startswith('BC'):
            ret += ' 4 |'
        else:
            ret += '   |'

        # d for diagram, l for labels
        d,l = CartanType_rank3_two_types.ascii_art(self).splitlines()[-2:]
        return ret + '\n' + d + '\n' + l

class CartanType_rank3_cycle(CartanType):
    """
    Hyperbolic Cartan types of rank 3 that are a cycle.
    """
    def __init__(self, index):
        """
        Initialize ``self``.
        """
        if index == 1:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,1,2))
        elif index == 2:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,1,-22))
        elif index == 3:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,1,3))
        elif index == 4:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,1,4))
        elif index == 5:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((2,-2,1))
        elif index == 6:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,2,2))
        elif index == 7:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((-2,2,1))
        elif index == 8:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,2,-22))
        elif index == 9:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,2,3))
        elif index == 10:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,2,-3))
        elif index == 11:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,2,4))
        elif index == 12:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,2,-4))
        elif index == 13:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-2,-22))
        elif index == 14:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-2,3))
        elif index == 15:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-2,-3))
        elif index == 16:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-2,4))
        elif index == 17:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-2,-4))
        elif index == 18:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-22,-22))
        elif index == 19:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-22,3))
        elif index == 20:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-22,-3))
        elif index == 21:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-22,4))
        elif index == 22:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-22,-4))
        elif index == 23:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((3,-3,1))
        elif index == 24:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,3,3))
        elif index == 25:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((-3,3,1))
        elif index == 26:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,3,4))
        elif index == 27:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,3,-4))
        elif index == 28:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-3,4))
        elif index == 29:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-3,-4))
        elif index == 30:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,4,-4))
        elif index == 31:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,4,4))
        elif index == 32:
            dd = CartanType_rank3_cycle._construct_dynkin_diagram((1,-4,4))
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

