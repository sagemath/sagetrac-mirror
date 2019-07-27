from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.graphs.graph import DiGraph
#from sage.rings.integer import sign
from six.moves import range
from sage.structure.sage_object import SageObject
from copy import copy
from sage.rings.all import ZZ, CC, infinity
from sage.graphs.all import Graph, DiGraph
from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
from sage.combinat.cluster_algebra_quiver.quiver_mutation_type import QuiverMutationType, QuiverMutationType_Irreducible, QuiverMutationType_Reducible, _edge_list_to_matrix
from sage.combinat.cluster_algebra_quiver.mutation_class import _principal_part, _digraph_mutate, _matrix_to_digraph, _dg_canonical_form, _mutation_class_iter, _digraph_to_dig6, _dig6_to_matrix
from sage.combinat.cluster_algebra_quiver.mutation_type import _connected_mutation_type, _mutation_type_from_data, is_mutation_finite

from sage.plot.colors import rainbow
from sage.graphs.graph_generators import GraphGenerators
from sage.all import e, pi, I

from sage.combinat.cluster_algebra_quiver.interact import cluster_interact


def DoubleBruhatDigraph(CartanType, u, v, word = False):
    '''Returns a quiver from the Weyl group elements u and v, along with lists of exchangeable and frozen vertices, using the algorithm outlined in Cluster Algebras III.

    INPUT:

    - ``u`` and ``v`` -- can be any of the following::
        * WeylGroupElement
        * Permutation
        * str - a string representation of a permutation
        * PermutationGroupElement
        * Matrix
    - CartanType -- the Cartan type of the underlying group, i.e. ['A', 5], ['C',6], etc.
    - ``word`` -- (default: False) a chosen reduced word for (u,v) as an element of W x W, using the convention that simple reflections in v are associated to negative integers.

    Examples::

        sage: from sage.combinat.cluster_algebra_quiver.double_bruhat_digraph import DoubleBruhatDigraph
        sage: W = WeylGroup(['C',4])
        sage: s1,s2,s3,s4 = W.simple_reflections()
        sage: D,F,S = DoubleBruhatDigraph(['C',4],s1*s2*s3*s4,s4*s3*s2*s1)
        sage: D
        Digraph on 12 vertices
        sage: F
        [-4, -3, -2, -1, 5, 6, 7, 8]
        sage: S
        [[-1, 1, 8], [-2, 2, 7], [-3, 3, 6], [-4, 4, 5]]


'''

    typeChar = CartanType[0]
    r = CartanType[1]

    lu = u.length()
    lv = v.length()

    # Constructs a reduced word if none is specified
    if not word:
        word = u.reduced_word()
        for n in v.reduced_word():
            word.append(-n)

    # Pads the list with frozen variables (Note: there will be other frozen variable in the final quiver)
    word = list(range(-r,0)) + word
    indices = list(range(-r,0)) + list(range(1,lu+lv+1))
    M = CartanMatrix(CartanType)

    # Determines the exchangeable vertices
    exchangeables = []
    frozen = []
    for character in indices:
        if iExchangeable(character, word, r):
            exchangeables.append(character)
        else:
            frozen.append(character)

    # Constructs the digraph described by BFZ in Cluster Algebras III
    dg = DiGraph()
    dg.add_vertices(indices)

    iks=[]

    for k0 in range(len(word)):
        iks.append(word[k0])
        for l0 in range(k0+1,len(word)):

            # There are some unfortunate indexing gymnastics to avoid the '0' problem
            ik = word[k0]
            il = word[l0]
            k = indices[k0]
            l = indices[l0]


            kplus = plus(k,word,r)
            lplus = plus(l,word,r)
            if kplus in indices:
                ikplus = word[indices.index(kplus)]
            if lplus in indices:
                ilplus = word[indices.index(lplus)]

            # This runs through BFZ's three conditions under which there can be an edge between k and l
            if k in exchangeables and l in exchangeables:

                # horizontal edges
                if l == kplus:
                    if il>0:
                        dg.add_edge(k,l)
                    else:
                        dg.add_edge(l,k)

                # inclined edges
                elif l<kplus and kplus<lplus and M[abs(ik)-1,abs(il)-1]<0 and il/abs(il) == ikplus/abs(ikplus):
                    if il<0:
                        dg.add_edge(k,l,(-M[abs(ik)-1,abs(il)-1],M[abs(il)-1,abs(ik)-1]))
                    else:
                        dg.add_edge(l,k,(-M[abs(il)-1,abs(ik)-1],M[abs(ik)-1,abs(il)-1]))

                # inclined edges
                elif l<lplus and lplus<kplus and M[abs(ik)-1,abs(il)-1]<0 and il/abs(il)==- ilplus/abs(ilplus):
                    if il<0:
                        dg.add_edge((k,l,(-M[abs(ik)-1,abs(il)-1],M[abs(il)-1,abs(ik)-1])))
                    else:
                        dg.add_edge(l,k,(-M[abs(il)-1,abs(ik)-1],M[abs(ik)-1,abs(il)-1]))

            elif (k in exchangeables).__xor__(l in exchangeables):
                if k in exchangeables:
                    # horizontal edges
                    if l == kplus:
                        if il>0:
                            dg.add_edge(k,l)
                        else:
                            dg.add_edge(l,k)

                    # inclined edges
                    elif l<kplus and kplus<lplus and M[abs(ik)-1,abs(il)-1]<0 and il/abs(il) == ikplus/abs(ikplus):
                        if il<0:
                            dg.add_edge(k,l,(-M[abs(ik)-1,abs(il)-1],M[abs(ik)-1,abs(il)-1]))
                        else:
                            dg.add_edge(l,k,(-M[abs(ik)-1,abs(il)-1],M[abs(ik)-1,abs(il)-1]))

                    # inclined edges
                    elif l<lplus and lplus<kplus and M[abs(ik)-1,abs(il)-1]<0 and il/abs(il)==- ilplus/abs(ilplus):
                        if il<0:
                            dg.add_edge((k,l,(-M[abs(ik)-1,abs(il)-1],M[abs(ik)-1,abs(il)-1])))
                        else:
                            dg.add_edge(l,k,(-M[abs(ik)-1,abs(il)-1],M[abs(ik)-1,abs(il)-1]))

                elif l in exchangeables:
                    # horizontal edges
                    if l == kplus:
                        if il>0:
                            dg.add_edge(k,l)
                        else:
                            dg.add_edge(l,k)

                    # inclined edges
                    elif l<kplus and kplus<lplus and M[abs(ik)-1,abs(il)-1]<0 and il/abs(il) == ikplus/abs(ikplus):
                        if il<0:
                            dg.add_edge(k,l,(-M[abs(il)-1,abs(ik)-1],M[abs(il)-1,abs(ik)-1]))
                        else:
                            dg.add_edge(l,k,(-M[abs(il)-1,abs(ik)-1],M[abs(il)-1,abs(ik)-1]))

                    # inclined edges
                    elif l<lplus and lplus<kplus and M[abs(ik)-1,abs(il)-1]<0 and il/abs(il)==- ilplus/abs(ilplus):
                        if il<0:
                            dg.add_edge((k,l,(-M[abs(il)-1,abs(ik)-1],M[abs(il)-1,abs(ik)-1])))
                        else:
                            dg.add_edge(l,k,(-M[abs(il)-1,abs(ik)-1],M[abs(il)-1,abs(ik)-1]))


    return dg, frozen, strings(indices, iks, r)

def plus(k,word,r):
    '''Returns the index 'k+', being the smallest index l in 'word' such that k<l and the index k and l entries of 'word' agree.  NOTE: the 'index' is assumed to be following the indexing conventions of DoubleBruhatDigraph, not the standard python conventions.

    INPUT:
    - ``k`` -- An index for `word` in the convention of DoubleBruhatQuiver.
    - ``word`` -- The reduced word constructed in DoubleBruhatQuiver.
    - ``r`` -- The number of additional frozen variables in the beginning of `word` (which is the rank of the underlying Cartan matrix).'''

    # fix for index 0
    if k < 0:
        ik = abs(word[k+r])
    else:
        ik = abs(word[k+r - 1])

    for l in range(k + 1,len(word) - r + 1):
        if l < 0 and abs(word[l + r]) == ik:
            return l
        elif l > 0 and abs(word[l+r - 1]) == ik:
            return l

    return len(word) - r + 1

def iExchangeable(k, word, r):
    '''Returns True if index k is exchangeable, and False otherwise.

    INPUT:

    - ``k`` -- An index for `word` in the convention of DoubleBruhatDigraph.
    - ``word`` -- The reduced word constructed in DoubleBruhatDigraph.
    - ``r`` -- The number of additional frozen variables in the beginning of `word` (which is the rank of the underlying Cartan matrix).'''

    if (k >= 1 and k <= (len(word) - r)) and (plus(k, word, r) >= 1 and plus(k, word, r) <= (len(word) - r)):
        return True

    return False

def strings(indices, iks, r):
    listk=[]
    for l in range(0, r):
        listl=[]
        for k in range(0,len(iks)):
            if abs(iks[l]) == abs(iks[k]):
                listl.append(indices[k])

        listk.append(listl)
    listk.reverse()
    return listk


class DoubleBruhatClusterQuiver(SageObject):
    """
    Documentation goes here
    """

    def __init__(self, CartanType, u, v, word = False):
        """
        TESTS:

        EXAMPLES:
        """
        #fix the following "from"
        from sage.combinat.cluster_algebra_quiver.double_bruhat_digraph import *
        D, F, S = DoubleBruhatDigraph(CartanType, u, v, word = False)
        self._quiver = ClusterQuiver(D, frozen = F)
        self._D = D
        self._F = F
        self._S = S
        self._cartan_type = CartanType

    def _repr_(self):
        """
        Return the description of ``self``.

        EXAMPLES:
        """
        self._quiver._repr_()

    def quiver(self):
        """
        change!Return the description of ``self``.

        EXAMPLES:
        """
        return self._quiver

    def show_sheets(self):
        n, m = self._quiver._n, self._quiver._m
        from sage.plot.plot import Graphics
        fig_size = 1.2
        if m == 0:
            width_factor = 1
            fig_size = fig_size*2*n/3
        else:
            width_factor = 2
            fig_size = fig_size*4*n/3

        M = CartanMatrix(self._cartan_type)
        listk = self._S
        plot_sequence = []
        number_of_sheets = 0
        for i in range(0, M.nrows()):
            for j in range(0, i+1):
                if M[i,j]<0:
                    L = self._quiver._digraph.subgraph(listk[i]+listk[j])
                    subFrozens = [value for value in self._F if value in listk[i]+listk[j]]
                    A = ClusterQuiver(L, frozen=subFrozens)
                    A.show()
                    #L
                    #L.show()
                    #figure out what subgraph is?

                    #plot_sequence.append(self._quiver.plot(subgraph = [self._quiver._digraph.subgraph(listk[i]+listk[j]),listk[i],listk[j]], circular=False, linear=True, center=(number_of_sheets*width_factor,0)))
                    number_of_sheets += 1
        #for p in plot_sequence:
        #    p.show()

        #All of this is for printing on one file!
        # sequence = []
        # for p in plot_sequence:
        #     sequence.append( p )
        # plot_obj = Graphics()
        # for elem in sequence:
        #     plot_obj += elem
        # plot_obj.show(axes=False, figsize=[fig_size*number_of_sheets,fig_size])

############saved############

    # def show_sheets(self):
    #     n, m = self._quiver._n, self._quiver._m
    #     from sage.plot.plot import Graphics
    #     fig_size = 1.2
    #     if m == 0:
    #         width_factor = 1
    #         fig_size = fig_size*2*n/3
    #     else:
    #         width_factor = 2
    #         fig_size = fig_size*4*n/3
    #     # if self._construction_type._description != 'DB':
    #     #     raise ValueError("Sheets are only valid for double Bruhat cells")
    #
    #     M = CartanMatrix(self._cartan_type)
    #     listk = self._S
    #     plot_sequence = []
    #     number_of_sheets = 0
    #     for i in range(0, M.nrows()):
    #         for j in range(0, i+1):
    #             if M[i,j]<0:
    #                 #figure out what subgraph is?
    #                 plot_sequence.append(self._quiver.plot(self._quiver._digraph.subgraph(listk[i]+listk[j]), center=(number_of_sheets*width_factor,0)))
    #                 #plot_sequence.append(self._quiver.plot(subgraph = [self._quiver._digraph.subgraph(listk[i]+listk[j]),listk[i],listk[j]], circular=False, linear=True, center=(number_of_sheets*width_factor,0)))
    #                 number_of_sheets += 1
    #     sequence = []
    #     for p in plot_sequence:
    #         sequence.append( p )
    #     plot_obj = Graphics()
    #     for elem in sequence:
    #         plot_obj += elem
    #     plot_obj.show(axes=False, figsize=[fig_size*number_of_sheets,fig_size])



################trying sheets################

# from sage.plot.plot import Graphics
# fig_size = 1.2
# if m == 0:
#     width_factor = 1
#     fig_size = fig_size*2*n/3
# else:
#     width_factor = 2
#     fig_size = fig_size*4*n/3
# if self._construction_type._description != 'DB':
#     raise ValueError("Sheets are only valid for double Bruhat cells")
#
# M = CartanMatrix(self._construction_type._cartan_type)
# listk = self._construction_type._strings
# plot_sequence = []
# number_of_sheets = 0
# for i in range(0, M.nrows()):
#     for j in range(0, i+1):
#         if M[i,j]<0:
#             plot_sequence.append(self.plot(subgraph = [self._digraph.subgraph(listk[i]+listk[j]),listk[i],listk[j]], circular=False, linear=True, center=(number_of_sheets*width_factor,0)))
#             number_of_sheets += 1
# sequence = []
# for p in plot_sequence:
#     sequence.append( p )
# plot_obj = Graphics()
# for elem in sequence:
#     plot_obj += elem
# plot_obj.show(axes=False, figsize=[fig_size*number_of_sheets,fig_size])
# else:
#
# plot = self.plot( circular=circular, directed=directed, mark=mark, save_pos=save_pos, greens=greens)
# if circular:
#     plot.show( figsize=[fig_size*3*(n+m)/4+1,fig_size*3*(n+m)/4+1] )
# else:
#     plot.show( figsize=[fig_size*n+1,fig_size*n+1] )
