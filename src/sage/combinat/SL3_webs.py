import array

from copy import copy

from sage.functions.trig import cos, sin
from sage.graphs.graph import Graph
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField
from sage.symbolic.constants import pi

from sage.combinat.SL3_webs_coloring import Color


def AreCrossing(Pos, P1, P2, P3, P4):
    """
    Check whether the line from P1 to P2 crosses the line from P3 to P4.

    INPUT:

    Pos -- positions ?

    P1, P2, P3, P4 are points expressed as x,y coordinates.

    EXAMPLES:: TODO
    """
    R = RealField(32)
    Crossing = False
    if P1 != P3 and P1 != P4 and P2 != P3 and P2 != P4:
        # Coordonnees des sommets
        x1 = R(Pos[P1][0])
        x2 = R(Pos[P2][0])
        x3 = R(Pos[P3][0])
        x4 = R(Pos[P4][0])
        y1 = R(Pos[P1][1])
        y2 = R(Pos[P2][1])
        y3 = R(Pos[P3][1])
        y4 = R(Pos[P4][1])
        # Parametre des droites
        a12 = y2 - y1  # Droite entre P1 et P2
        b12 = x1 - x2
        c12 = x2 * y1 - x1 * y2
        a13 = y3 - y1
        b13 = x1 - x3
        c13 = x3 * y1 - x1 * y3
        a23 = y3 - y2
        b23 = x2 - x3
        c23 = x3 * y2 - x2 * y3
        if ((a12 * x3 + b12 * y3 + c12) * (a12 * x4 + b12 * y4 + c12) < 0 and
            (a13 * x2 + b13 * y2 + c13) * (a13 * x4 + b13 * y4 + c13) > 0 and
            (a23 * x1 + b23 * y1 + c23) * (a23 * x4 + b23 * y4 + c23) > 0):
            # Determiner si P4 et P3 ne sont pas du meme cote de la droite passant par P1 et P2, determiner si P2 et P4 sont du meme cote de la droite passsant par P1 et P3 puis determiner si P1 et P4 sont du meme cote de la droite passant par P2 et P3
            return True
    return Crossing


class Web(object):
    """
    EXAMPLES::

        sage: from sage.combinat.SL3_webs import Web
        sage: W = Web([1,2,3,4],[1,1,-1,-1]); W
        Web with 4 vertices (4 on the boundary)
    """
    def __init__(self, InitStates, InitColor, **kwargs):
        # (Graph = , Pos = (En forme KK), BlackVertices =, WhiteVertices =, n =)

        self._InitStates = copy(InitStates)
        self._InitColor = copy(InitColor)

        self._n = len(self._InitStates)
        self._Graph = Graph(multiedges=True)
        for i in range(self._n):
            self._Graph.add_vertex(i + 1)

        self._Pos = {i + 1: [i + 1, 0] for i in range(self._n)}
        self._BlackVertices = [i + 1 for i in range(self._n)
                               if InitColor[i] == 1]
        self._WhiteVertices = [i + 1 for i in range(self._n)
                               if InitColor[i] == -1]

        def SetHeight(G, i, n, Pos, HeightNewVertices, ThisStepVertices,
                      NewVertex1='type D2', NewVertex2='type S'):
            # positions i and i+1 are the positions where something is being attached
            if NewVertex1 == 'type D2':  # Move type D2
                # this is a move where two already-existing vertices are being joined by a line.
                # it is helpful for the layout to adjust the x-coordinate of the lower of the two vertices to bring it
                # closer to the other one.
                if Pos[ThisStepVertices[i]][1] < Pos[ThisStepVertices[i + 1]][1]:
                    # Le premier est plus bas
                    Pos[ThisStepVertices[i]][0] = (Pos[ThisStepVertices[i]][0]) * 2 / 3 + (Pos[ThisStepVertices[i + 1]][0]) / 3

                    for Edge in G.edges(labels=False):
                        if Edge[0] > n and Edge[1] > n:
                            for Neighbor in G.neighbors(ThisStepVertices[i]):
                                while AreCrossing(Pos, ThisStepVertices[i], Neighbor, Edge[0], Edge[1]):
                                    Pos[ThisStepVertices[i]][1] -= 1
                elif Pos[ThisStepVertices[i]][1] > Pos[ThisStepVertices[i + 1]][1]:  # Le deuxieme est plus bas
                    Pos[ThisStepVertices[i]][0] = (Pos[ThisStepVertices[i]][0]) * 1 / 3 + (Pos[ThisStepVertices[i + 1]][0]) * 2 / 3
                    for Edge in G.edges(labels=False):
                        if Edge[0] > n and Edge[1] > n:
                            while AreCrossing(Pos, ThisStepVertices[i], ThisStepVertices[i + 1], Edge[0], Edge[1]):
                                Pos[ThisStepVertices[i]][1] -= 1
                else:  # Hauteurs egales
                    Pos[ThisStepVertices[i]][0] = (Pos[ThisStepVertices[i]][0]) * 2 / 3 + (Pos[ThisStepVertices[i + 1]][0]) * 1 / 3
                    for Edge in G.edges(labels=False):
                        if Edge[0] > n and Edge[1] > n:
                            while AreCrossing(Pos, ThisStepVertices[i], ThisStepVertices[i + 1], Edge[0], Edge[1]):
                                Pos[ThisStepVertices[i]][1] -= 1
            elif NewVertex2 == 'type S':  # Move type S
                Pos[NewVertex1] = [(Pos[ThisStepVertices[i]][0] + Pos[ThisStepVertices[i + 1]][0]) / 2, HeightNewVertices]  # Position si rien ne croise
                for Edge in G.edges(labels=False):
                    if Edge[0] > n and Edge[1] > n:
                        while AreCrossing(Pos, NewVertex1, ThisStepVertices[i], Edge[0], Edge[1]) or AreCrossing(Pos, NewVertex1, ThisStepVertices[i + 1], Edge[0], Edge[1]):
                            Pos[NewVertex1][1] -= 1
            else:  # Move type D1
                Pos[NewVertex1] = [(2 * Pos[ThisStepVertices[i]][0] + Pos[ThisStepVertices[i + 1]][0]) / 3, HeightNewVertices]  # Position si rien ne croise
                Pos[NewVertex2] = [(Pos[ThisStepVertices[i]][0] + 2 * Pos[ThisStepVertices[i + 1]][0]) / 3, HeightNewVertices]
                for Edge in G.edges(labels=False):
                    if Edge[0] > n and Edge[1] > n:
                        while AreCrossing(Pos, NewVertex1, ThisStepVertices[i], Edge[0], Edge[1]) or AreCrossing(Pos, NewVertex2, ThisStepVertices[i + 1], Edge[0], Edge[1]):
                            Pos[NewVertex1][1] -= 1
                            Pos[NewVertex2][1] -= 1
            return Pos

        # building the web
        ThisStepVertices = self._Graph.vertices()
        ThisStepStates = list(self._InitStates)
        ThisStepColors = InitColor
        Step = 0
        NewVertex = self._n + 1
        HeightNewVertices = 0
        LastStepVertices = ["Not a number, not empty"]

        while LastStepVertices != ThisStepVertices:
            LastStepVertices = list(ThisStepVertices)
            i = 0
            Step += 1
            HeightNewVertices -= 1
            while i < len(ThisStepVertices) - 1:
                if ThisStepColors[i] == ThisStepColors[i + 1]:
                    # Same Colors
                    if ThisStepStates[i] > ThisStepStates[i + 1]:
                        # States 1&0,0&-1 and 1&-1: move (type S) on ThisStepVertices[i] and ThisStepVertices[i+1].
                        if ThisStepColors[i] == -1:
                            self._BlackVertices.append(NewVertex)
                        else:
                            self._WhiteVertices.append(NewVertex)
                        self._Graph.add_edges([[ThisStepVertices[i], NewVertex], [ThisStepVertices[i + 1], NewVertex]])
                        self._Graph.set_edge_label(ThisStepVertices[i], NewVertex, ThisStepStates[i] * ThisStepColors[i])
                        self._Graph.set_edge_label(ThisStepVertices[i + 1], NewVertex, ThisStepStates[i + 1] * ThisStepColors[i + 1])
                        self._Pos = SetHeight(self._Graph, i, self._n, self._Pos, HeightNewVertices, ThisStepVertices, NewVertex)
                        ThisStepVertices[i: i + 2] = [NewVertex]
                        ThisStepStates[i: i + 2] = [ThisStepStates[i] + ThisStepStates[i + 1]]
                        ThisStepColors[i: i + 2] = [-ThisStepColors[i]]
                        NewVertex += 1
                        i += 1
                    else:
                        # Other states : no move (type S) on ThisStepVertices[i] and ThisStepVertices[i+1] performed.
                        i += 1
                else:
                    # Different Colors
                    if ThisStepStates[i] < ThisStepStates[i + 1] or ThisStepStates[i] == ThisStepStates[i + 1] != 0:
                        # Other states, nothing to do.
                        i += 1
                    elif ThisStepStates[i] == 1 == -ThisStepStates[i + 1]:
                        # Move (type D2)on ThisStepVertices[i] and ThisStepVertices[i+1].
                        self._Graph.add_edge(ThisStepVertices[i], ThisStepVertices[i + 1])
                        self._Graph.set_edge_label(ThisStepVertices[i], ThisStepVertices[i + 1], ThisStepColors[i])
                        self._Pos = SetHeight(self._Graph, i, self._n, self._Pos, HeightNewVertices, ThisStepVertices)
                        ThisStepVertices[i: i + 2] = []
                        ThisStepStates[i: i + 2] = []
                        ThisStepColors[i: i + 2] = []
                    else:
                        # Move (Type D1) on ThisStepVertices[i] and ThisStepVertices[i+1].
                        if ThisStepColors[i] == -1:
                            self._BlackVertices.append(NewVertex)
                            self._WhiteVertices.append(NewVertex + 1)
                        else:
                            self._WhiteVertices.append(NewVertex)
                            self._BlackVertices.append(NewVertex + 1)
                        self._Graph.add_edges([[ThisStepVertices[i], NewVertex], [ThisStepVertices[i + 1], NewVertex + 1], [NewVertex, NewVertex + 1]])
                        self._Graph.set_edge_label(ThisStepVertices[i], NewVertex, ThisStepStates[i] * ThisStepColors[i])
                        self._Graph.set_edge_label(ThisStepVertices[i + 1], NewVertex + 1, ThisStepStates[i + 1] * ThisStepColors[i + 1])
                        self._Graph.set_edge_label(NewVertex, NewVertex + 1, -ThisStepColors[i])
                        self._Pos = SetHeight(self._Graph, i, self._n, self._Pos, HeightNewVertices, ThisStepVertices, NewVertex, NewVertex + 1)
                        ThisStepVertices[i: i + 2] = [NewVertex, NewVertex + 1]
                        ThisStepStates[i: i + 2] = [ThisStepStates[i] - 1, ThisStepStates[i + 1] + 1]
                        ThisStepColors[i: i + 2] = [-ThisStepColors[i], -ThisStepColors[i + 1]]
                        NewVertex += 2
                        i += 2

        # Calcul des positions FP
        self._PosFP = copy(self._Pos)
        PlusBas = 0  # we find the lowest y-coordinate
        for V in self._PosFP:
            if self._PosFP[V][1] < PlusBas:
                PlusBas = self._PosFP[V][1]

        def PlaneToCircle(L, a, n):
            return [-1 / a * (L[1] - a) * cos(2 / n * pi * L[0]),
                    -1 / a * (L[1] - a) * sin(2 / n * pi * L[0])]

        for i in self._PosFP:
            self._PosFP[i] = PlaneToCircle(self._PosFP[i], PlusBas - 0.5, self._n)

    def __repr__(self):
        """
        Return a string representation.
        """
        return "Web with %s vertices (%s on the boundary)" % (self._Graph.order(), self._n)

    def __eq__(self, other):  # Teste les etiquettes!
        """
        Compare ``self`` with ``other``.
        """
        return (self._InitStates == other._InitStates and
                self._InitClasp == other._InitClasp and
                self._Graph.edges() == other._Graph.edges())

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):  # Teste les etiquettes?
        return hash((tuple(self._InitStates),
                     tuple(self._InitClasp),
                     self._Graph.copy(immutable=True),
                     tuple(self._Graph.edges())))

    def plot(self, **kwargs):
        """
        Return a plot of the web.

        INPUT:

        Layout defaults to circular, but if Layout='KK', then it is laid out the way it is built from the KK encoding.

            FigSize=[x,y], EdgeLabels = True/False

        EXAMPLES::

            sage: from sage.combinat.SL3_webs import Web
            sage: W = Web([1,2,3,4],[1,1,-1,-1]); W
            sage: W.plot()
            Graphics object consisting of 5 graphics primitives
        """
        FigSize = kwargs.get("FigSize")
        if FigSize is None:
            FigSize = [5, 5]

        vcolors = {"grey": self._BlackVertices, "white": self._WhiteVertices}

        if kwargs.get("Layout") == 'KK':
            self._Graph.plot(pos=self._Pos, vertex_colors=vcolors,
                             edge_labels=kwargs.get("EdgeLabels"),
                             figsize=FigSize, aspect_ratio='automatic')
        else:
            self._Graph.plot(pos=self._PosFP, vertex_colors=vcolors,
                             edge_labels=kwargs.get("EdgeLabels"),
                             figsize=FigSize, aspect_ratio='automatic')

    @cached_method
    def theR(self):
        """
        Return the polynomial ring in which the web invariant lives.

        EXAMPLES::

            sage: from sage.combinat.SL3_webs import Web
            sage: W = Web([1,2,3,4],[1,1,-1,-1])
            sage: W.theR()
            Multivariate Polynomial Ring in x11, x12, x13, x21, x22, x23, x31, x32, x33, x41, x42, x43 over Rational Field
        """
        n = self._n
        vars = [f"x{i}{j}" for i in range(1, n + 1) for j in range(1, 4)]
        return PolynomialRing(QQ, vars)

    @cached_method
    def invariant(self):
        """
        Compute the invariant associated to the web.

        EXAMPLES:: ?
        """
        n = self._n
        H = self._Graph
        He = H.edges()
        Hl = array.array('b')
        Hn = []
        Sum = {}
        lie = self.layout()
        # Retire les etiquettes
        for E in He:
            Hl.append(-1)
            Ne = []
            # this is going to be the list of all edges that have a
            # vertex in common with edge E.
            Ne = [He.index(x) for x in H.edges_incident(E[0])]
            Ne += [He.index(x) for x in H.edges_incident(E[1])]
            Hn.append(copy(Ne))

        Color(Hl, Hn, He, 0, Sum, lie, n)
        return self.theR()(Sum)

    @cached_method
    def layout(self):
        """
        Doc ?
        """
        G = self._Graph
        n = self._n
        lie = []
        for I in range(n + 1, len(G.vertices()) + 1):
            Neighbors = G.neighbors(I)
            A = Neighbors[0]
            B = Neighbors[1]
            ea = G.edges(labels=False).index(tuple(sorted([A, I])))
            eb = G.edges(labels=False).index(tuple(sorted([B, I])))
            lie.append([ea, eb])
        return lie
