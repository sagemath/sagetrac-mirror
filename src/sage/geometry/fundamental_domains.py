r"""
Fundamental domains for subgroups of PSL(2,ZZ)

- fundamental domains on the hyperbolic half plane

- congruence test (i.e. does there exists n such that a conjugate of the group
  contains Gamma(n))

- normalisator and other standard group operation

REFERENCES:

- Nielsen, Schreier
- Newman
- Stothers had written some non available stuff about that
- Kulkarni, Verrill
"""
#################################
# DATA STRUCTURE FOR COSET GRAPH

from sage.structure.sage_object import SageObject
from sage.plot.misc import options
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.matrix.constructor import matrix

from sage.graphs.triangle_graph import TriangleGraph_2_3_infinity
from sage.geometry.hyperbolic_space.all import HyperbolicPlane
from sage.all import CC
from sage.functions.log import exp
from sage.symbolic.constants import pi
from sage.symbolic.all import I
from sage.rings.infinity import Infinity


class FundamentalDomain(SageObject):
    r"""
    Fundamental domain for Fuchsian groups (actually subgroups of PSL(2,ZZ))
    """
    def __init__(self, graph, name=None):
        self._graph = graph
        if name is None:
            self._name = 'a subgroup of index %d' % graph.num_verts()
        gens, tree, lifts = graph.spanning_tree()
        self._gens = gens
        self._tree = tree
        self._lifts = lifts

    def _repr_(self):
        """
        Return the string representation of self.
        """
        return "Fundamental domain of %s" % self._name

    @options(color='blue')
    def plot(self, **options):
        HH = HyperbolicPlane().UHP()
        pol = HH.polygon([Infinity, CC(0, 1), CC(0, 0), CC(exp(I * pi / 3))])
        G = pol.plot(face_color=options['color'], alpha='0.5')
        for g in self._lifts[1:]:
            G += (g * pol).plot(face_color=options['color'], alpha=0.2)
        return G

    @options(aspect_ratio=1, ymax=2)
    def show(self, **options):
        self.plot().show(**options)

#############################################################
# TRICK FOR PSL REPRESENTATIVE (to be moved in matrix groups)


def PSL_rep(m, n):
    r"""
    Choose the representative for PSL(2, Z/nZ)

    INPUT:

    - ``m`` -- a matrix

    - ``n`` -- positive integer
    """
    Zn = Zmod(n)
    #print "SPL computation... from\n", m

    G = MatrixSpace(Zn, 2)
    m_mod = G(m)
    #print "reduced mod %d\n" %n, m_mod
    inv = G(matrix([[-1, 0], [0, -1]]))
    a, b, c, d = G(m_mod).list()

    if n % 2:
        limit = Zmod(n)((n + 1) // 2)
        if a != 0:
            if a < limit:
                res = m_mod
            else:
                res = inv * m_mod

        else:
            if b < limit:
                res = m_mod
            else:
                res = inv * m_mod

    else:
        half = Zn(n // 2)
        if a != 0 and a != half:
            if a < half:
                res = m_mod
            else:
                res = inv * m_mod

        else:
            if b < half:
                res = m_mod
            else:
                res = inv * m_mod

    #print "get\n", res
    res.set_immutable()
    return res

########################################################
# COSETS GRAPHS FOR GAMMA'S (to be moved in sage.modular)


def gamma_triangle_graph(n, return_mapping=False):
    r"""
    Coset graph for the principal congruence subgroup.
    """
    g2 = matrix([[0, -1], [1, 0]])
    g2.set_immutable()
    g3 = matrix([[1, -1], [1, 0]])
    g3.set_immutable()

    id = PSL_rep(matrix([[1, 0], [0, 1]]), n)
    g2_mod = PSL_rep(g2, n)
    #print "g2\n",g2_mod
    g3_mod = PSL_rep(g3, n)
    #print "g3\n", g3_mod

    nb = 1
    edges = [[0], [0]]
    matrices = {id: 0}

    l = set([id])

    while l:
        x0 = l.pop()
        i0 = matrices[x0]

        #print "pop\n", x0

        # comuting the g2 neighbor
        #print x0*g2_mod
        x1 = PSL_rep(x0 * g2_mod, n)
        #print "g2 neighbor\n", x1

        if not x1 in matrices:
            matrices[x1] = nb
            nb += 1
            edges[0].append(None)
            edges[1].append(None)
            l.add(x1)

        i1 = matrices[x1]
        edges[0][i0] = i1

        # computing the g3 neighbor
        x1 = PSL_rep(x0 * g3_mod, n)
        #print "g3 neighbor\n", x1

        if not x1 in matrices:
            matrices[x1] = nb
            nb += 1
            edges[0].append(None)
            edges[1].append(None)
            l.add(x1)

        i1 = matrices[x1]
        edges[1][i0] = i1

    if return_mapping:
        return TriangleGraph_2_3_infinity(edges), matrices
    else:
        return TriangleGraph_2_3_infinity(edges)


def gamma0_triangle_graph(n, return_mapping=False):
    r"""
    Return the coset graph for gamma0.
    """
    T = matrix([[1, 1], [0, 1]])
    g2 = matrix([[0, -1], [1, 0]])
    g2.set_immutable()
    g3 = matrix([[1, -1], [1, 0]])
    g3.set_immutable()

    id = PSL_rep(matrix([[1, 0], [0, 1]]), n)
    g2_mod = PSL_rep(g2, n)
    #print "g2\n",g2_mod
    g3_mod = PSL_rep(g3, n)
    #print "g3\n", g3_mod

    nb = 1
    edges = [[0], [0]]
    matrices = {id: 0}
    # adding the coset as known matrices
    for i in xrange(1, n):
        matrices[PSL_rep(T ** i, n)] = 0

    l = set([id])

    while l:
        x0 = l.pop()
        i0 = matrices[x0]

        # computing the g2 neighbor
        x1 = PSL_rep(x0 * g2_mod, n)

        if not x1 in matrices:
            for i in xrange(n):
                matrices[PSL_rep(T ** i * x1, n)] = nb
            nb += 1
            edges[0].append(None)
            edges[1].append(None)
            l.add(x1)

        i1 = matrices[x1]
        edges[0][i0] = i1

        # computing the g3 neighbor
        x1 = PSL_rep(x0 * g3_mod, n)
        #print "g3 neighbor\n", x1

        if not x1 in matrices:
            for i in xrange(n):
                matrices[PSL_rep(T ** i * x1, n)] = nb
            nb += 1
            edges[0].append(None)
            edges[1].append(None)
            l.add(x1)

        i1 = matrices[x1]
        edges[1][i0] = i1

    if return_mapping:
        return TriangleGraph_2_3_infinity(edges), matrices
    else:
        return TriangleGraph_2_3_infinity(edges)
