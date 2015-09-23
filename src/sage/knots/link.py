r"""
Links

A knot is defined as embedding of the circle `\mathbb{S}^1` in the 3-dimensional
sphere `\mathbb{S}^3`, considered up to ambient isotopy. They represent the physical
idea of a knotted rope, but with the particularity that the rope is closed. That
is, the ends of the rope are joined.

A link is an embedding of one or more copies of `\mathbb{S}^1` in `\mathbb{S}^3`,
considered up to ambient isotopy. That is, a link represents the idea of one or more
tied ropes. Every knot is a link, but not every link is a knot.

Generically, the projection of a link on `\RR^2`
is a curve with crossings. The crossings are represented to show which strand goes
over the other. This curve is called a planar diagram of the link. If we remove the
crossings, the resulting connected components are segments. These segments are
called the edges of the diagram.

REFERENCES:

.. [Wikipedia] Wikipedia article :wikipedia:`Knot_(mathematics)`

AUTHORS:

- Miguel Angel Marco Buzunariz
- Amit Jamadagni
"""

##############################################################################
#       Copyright (C) 2014  Miguel Angel Marco Buzunariz
#                           Amit Jamadagni
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################


from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.symbolic.ring import SR
from sage.rings.integer import Integer
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.functions.generalized import sign
from sage.plot.line import line
from sage.plot.bezier_path import bezier_path
from sage.misc.flatten import flatten
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from copy import deepcopy, copy


class Link(object):
    r"""
    A link.

    A link can be created by using one of the conventions mentioned below:

    Braid:

    - The closure of a braid is a link::

        sage: B = BraidGroup(8)
        sage: L = Link(B([-1, -1, -1, -2,1, -2,3,-2,3]))
        sage: L
        Link with 1 components represented by 9 crossings
        sage: L = Link(B([1, 2, 1, -2,-1]))
        sage: L
        Link with 2 components represented by 5 crossings

      Note that the strands of the braid that have no crossings at all
      are removed.

    - Oriented Gauss Code:

      Label the crossings from `1` to `n` (where `n` is the number of
      crossings) and start moving along the link. Trace every component of
      the link, by starting at a particular point on one component of the
      link and writing down each of the crossings that you encounter until
      returning to the starting point. The crossings are written with sign
      depending on whether we cross them as over or undercrossing. Each
      component is then represented as a list whose elements are the
      crossing numbers. A second list of `+1` and `-1`'s keeps track of
      the orientation of each crossing::

        sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1,-1,-1,-1,+1,+1,-1,+1]])
        sage: L
        Link with 1 components represented by 8 crossings

      For links there may be more than one component and the input is
      as follows::

        sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
        sage: L
        Link with 3 components represented by 4 crossings

    - Planar Diagram (PD) Code:

      The diagram of the link is formed by segments that are adjacent to
      the crossings. Label each one of this segments with a positive number,
      and for each crossing, write down the four incident segments. The
      order of these segments is clockwise, starting with the incoming
      undercrossing.

      There is no particular distinction between knots and links for
      this input.

    EXAMPLES:

    One of the representations of the trefoil knot::

        sage: L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])
        sage: L
        Link with 1 components represented by 3 crossings

    One of the representations of the Hopf link::

        sage: L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
        sage: L
        Link with 2 components represented by 2 crossings

    We can construct links from from the braid group::

        sage: B = BraidGroup(8)
        sage: L = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
        sage: L
        Link with 2 components represented by 8 crossings
        sage: L = Link(B([1, 2, 1]))
        sage: L
        Link with 2 components represented by 3 crossings

    .. WARNING::

        Equality of knots is done by comparing the corresponding braids,
        which may give false negatives.

    .. TODO::

        Implement methods to creating new links from previously created links.
    """
    def __init__(self, data):
        """
        Initialize ``self``.

        TESTS::

            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, -1, -1, -2,1, -2, 3, -2]))
            sage: TestSuite(L).run()
            sage: L = Link(B([1, 2, 1]))
            sage: TestSuite(L).run()
            sage: L = Link([[1, 1, 2, 2]])
            sage: TestSuite(L).run()
        """
        if isinstance(data, list):
            if not data:
                raise ValueError("does not accept empty list as argument")

            if len(data) != 2 or not all(isinstance(i, list) for i in data[0]):
                for i in data:
                    if len(i) != 4:
                        raise ValueError("invalid PD code: crossings must be represented by four segments")
                    else:
                        flat = flatten(data)
                        if any(flat.count(i) != 2 for i in set(flat)):
                            raise ValueError("invalid PD code: each segment must appear twice")
                self._pd_code = data
                self._oriented_gauss_code = None
                self._braid = None

            else:
                flat = flatten(data[0])
                a, b = max(flat), min(flat)
                if 2 * len(data[1]) != len(flat) or set(range(b, a + 1)) - set([0]) != set(flat):
                    raise ValueError("invalid input: data is not a valid oriented Gauss code")
                self._oriented_gauss_code = data
                self._pd_code = None
                self._braid = None

        else:
            from sage.groups.braid import Braid
            if isinstance(data, Braid):
                self._braid = data
                self._oriented_gauss_code = None
                self._pd_code = None

            else:
                raise ValueError("invalid input: data must be either a list or a braid")

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT: string representation

        EXAMPLES::

            sage: B = BraidGroup(8)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L
            Link with 1 components represented by 4 crossings
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L
            Link with 3 components represented by 4 crossings
        """
        number_of_components = self.number_of_components()
        pd_len = len(self.pd_code())
        return 'Link with {} components represented by {} crossings'.format(number_of_components, pd_len)

    def __eq__(self, other):
        """
        Check equality.

        TESTS::

            sage: B = BraidGroup(8)
            sage: L1 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L2 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L1 == L2
            True
            sage: L3 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
            sage: L1 == L3
            False
        """
        if not isinstance(other, self.__class__):
            return False
        if self._pd_code is not None:
            if self.pd_code() == other.pd_code():
                return True
        if self._oriented_gauss_code is not None:
            if self.oriented_gauss_code() == other.oriented_gauss_code():
                return True
        return self.braid() == other.braid()

    def __ne__(self, other):
        """
        Check inequality.

        TESTS::

            sage: B = BraidGroup(8)
            sage: L1 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L2 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L1 != L2
            False
            sage: L3 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
            sage: L1 != L3
            True
        """
        return not self.__eq__(other)

    def braid(self):
        """
        Return a braid representation of ``self``.

        OUTPUT: an element in the braid group

        EXAMPLES::

            sage: L = Link([[2,3,1,4],[4,1,3,2]])
            sage: L.braid()
            s^2
            sage: L = Link([[[-1, 2, -3, 1, -2, 3]], [-1, -1, -1]])
            sage: L.braid()
            s^-3
            sage: L = Link([[1,8,2,7],[8,4,9,5],[3,9,4,10],[10,1,7,6],[5,3,6,2]])
            sage: L.braid()
            (s0*s1^-1)^2*s1^-1
        """
        if self._braid is not None:
            return self._braid

        from sage.groups.braid import BraidGroup
        comp = self._isolated_components_()
        if len(comp) > 1:
            L1 = Link(comp[0])
            L2 = Link(flatten(comp[1:], max_level=1))
            b1 = L1.braid()
            b2 = L2.braid()
            n1 = b1.parent().strands()
            n2 = b2.parent().strands()
            t1 = list(b1.Tietze())
            t2 = [sign(x)*(abs(x) + n1) for x in b2.Tietze()]
            B = BraidGroup(n1 + n2)
            return B(t1 + t2)

        # look for possible Vogel moves, perform them and call recursively to the modified link
        pd_code = self.pd_code()
        newedge = max(flatten(pd_code)) + 1
        for region in self.regions():
            n = len(region)
            for i in range(n-1):
                a = region[i]
                seifcirca = [x for x in self.seifert_circles() if abs(a) in x]
                for j in range(i+1,n):
                    b = region[j]
                    seifcircb = [x for x in self.seifert_circles() if abs(b) in x]
                    if seifcirca != seifcircb and sign(a) == sign(b):
                        tails, heads = self._directions_of_edges_()

                        newPD = deepcopy(pd_code)
                        if sign(a) == 1:
                            C1 = newPD[newPD.index(heads[a])]
                            C1[C1.index(a)] = newedge + 1
                            C2 = newPD[newPD.index(tails[b])]
                            C2[C2.index(b)] = newedge + 2
                            newPD.append([newedge + 3, a, b, newedge])
                            newPD.append([newedge + 2, newedge + 1, newedge + 3, newedge])
                            self._braid = Link(newPD).braid()
                            return self._braid
                        else:
                            C1 = newPD[newPD.index(heads[-a])]
                            C1[C1.index(-a)] = newedge + 1
                            C2 = newPD[newPD.index(tails[-b])]
                            C2[C2.index(-b)] = newedge + 2
                            newPD.append([newedge + 2, newedge, newedge + 3, newedge + 1])
                            newPD.append([newedge + 3, newedge, -b , -a])
                            self._braid = Link(newPD).braid()
                            return self._braid

        # We are in the case where no Vogel moves are necessary.
        G = DiGraph()
        G.add_vertices([tuple(c) for c in self.seifert_circles()])
        for i,c in enumerate(pd_code):
            if self.orientation()[i] == 1:
                a  = [x for x in self.seifert_circles() if c[1] in x][0]
                b  = [x for x in self.seifert_circles() if c[0] in x][0]
                G.add_edge(tuple(a), tuple(b))
            else:
                a  = [x for x in self.seifert_circles() if c[0] in x][0]
                b  = [x for x in self.seifert_circles() if c[3] in x][0]
                G.add_edge(tuple(a), tuple(b))
        ordered_cycles = G.all_simple_paths(starting_vertices=G.sources(), ending_vertices=G.sinks())[0]
        B = BraidGroup(len(ordered_cycles))
        available_crossings = copy(pd_code)
        crossing = [x for x in pd_code if set(ordered_cycles[0]).intersection(set(x))][0]
        available_crossings.remove(crossing)
        status = [None for i in ordered_cycles]
        if self.orientation()[pd_code.index(crossing)] == 1:
            b = B([1])
            status[0] = crossing[2]
            status[1] = crossing[3]
        else:
            b = B([-1])
            status[0] = crossing[1]
            status[1] = crossing[2]
        counter = 0
        while available_crossings:
            possibles = [x for x in available_crossings if status[counter] in x]
            if len(status) < counter + 2 or status[counter + 1] != None:
                possibles = [x for x in possibles if status[counter + 1] in x]
            if possibles:
                added = possibles[0]
                if self.orientation()[self.pd_code().index(added)] == 1:
                    b *= B([counter + 1])
                    status[counter] = added[2]
                    status[counter + 1] = added[3]
                else:
                    b *= B([-counter - 1])
                    status[counter] = added[1]
                    status[counter + 1] = added[2]
                if counter > 0:
                    counter -= 1
                available_crossings.remove(added)
            else:
                counter += 1
        self._braid = b
        return b

    def _directions_of_edges_(self):
        r"""
        Return the directions of the edges given by the PD code of ``self``.

        OUTPUT:

        A tuple of two dictionaries. The first one assigns
        each edge of the PD code to the crossing where it starts.
        The second dictionary assigns it to where it ends.

        EXAMPLES::

            sage: L = Link([[1, 3, 2, 4], [2, 3, 1, 4]])
            sage: L._directions_of_edges_()
            ({1: [2, 3, 1, 4], 2: [1, 3, 2, 4], 3: [1, 3, 2, 4], 4: [2, 3, 1, 4]},
             {1: [1, 3, 2, 4], 2: [2, 3, 1, 4], 3: [2, 3, 1, 4], 4: [1, 3, 2, 4]})

        ::

            sage: L = Link([[1,5,2,4],[5,3,6,2],[3,1,4,6]])
            sage: L._directions_of_edges_()
            ({1: [3, 1, 4, 6],
              2: [1, 5, 2, 4],
              3: [5, 3, 6, 2],
              4: [3, 1, 4, 6],
              5: [1, 5, 2, 4],
              6: [5, 3, 6, 2]},
             {1: [1, 5, 2, 4],
              2: [5, 3, 6, 2],
              3: [3, 1, 4, 6],
              4: [1, 5, 2, 4],
              5: [5, 3, 6, 2],
              6: [3, 1, 4, 6]})

        ::

            sage: L = Link([[1,2,3,3], [2,4,5,5], [4,1,7,7]])
            sage: L._directions_of_edges_()
            ({1: [4, 1, 7, 7],
              2: [1, 2, 3, 3],
              3: [1, 2, 3, 3],
              4: [2, 4, 5, 5],
              5: [2, 4, 5, 5],
              7: [4, 1, 7, 7]},
             {1: [1, 2, 3, 3],
              2: [2, 4, 5, 5],
              3: [1, 2, 3, 3],
              4: [4, 1, 7, 7],
              5: [2, 4, 5, 5],
              7: [4, 1, 7, 7]})
        """
        tails = {}
        heads = {}
        pd_code = self.pd_code()
        for C in pd_code:
            tails[C[2]] = C
            a = C[2]
            D = C
            while not a in heads:
                next_crossing = [x for x in pd_code if a in x and x != D]
                if not next_crossing:
                    heads[a] = D
                    tails[a] = D
                    if D[0] == a:
                        a = D[2]
                    elif D[1] == a:
                        a = D[3]
                    else:
                        a = D[1]
                else:
                    heads[a] = next_crossing[0]
                    tails[a] = D
                    D = next_crossing[0]
                    a = D[(D.index(a)+2) % 4]

        unassigned = set(flatten(pd_code)).difference(set(tails.keys()))
        while unassigned:
            a = unassigned.pop()
            for x in pd_code:
                if a in x:
                    D = x
                    break
            while not a in heads:
                tails[a] = D
                for x in pd_code:
                    if a in x and x != D:
                        next_crossing = x
                        break
                heads[a] = next_crossing
                D = next_crossing
                a = D[(D.index(a)+2) % 4]
                if a in unassigned:
                    unassigned.remove(a)
        return tails, heads

    def oriented_gauss_code(self):
        """
        Return the oriented Gauss code of ``self``.

        The oriented Gauss code has two parts:

        a. the Gauss code

        b. the orientation of each crossing

        The following orientation was taken into consideration for
        construction of knots:

        From the outgoing of the overcrossing if we move in the clockwise
        direction to reach the outgoing of the undercrossing then we label
        that crossing as `-1`.

        From the outgoing of the overcrossing if we move in the anticlockwise
        direction to reach the outgoing of the undercrossing then we label
        that crossing as `+1`.

        One more consideration we take in while constructing the orientation
        is the order of the orientation is same as the ordering of the
        crossings in the Gauss code.

        .. NOTE::

            Convention: under is denoted by `-1`, and over by `+1` in the
            crossing info.

        OUTPUT: oriented Gauss code

        EXAMPLES::

            sage: L = Link([[1, 11, 2, 10], [6, 2, 7, 3], [3, 12, 4, 9], [9, 5, 10, 6], [8, 1, 5, 4], [11, 8, 12, 7]])
            sage: L.oriented_gauss_code()
            [[[-1, 2, -3, 5], [4, -2, 6, -5], [-4, 1, -6, 3]], [-1, 1, 1, 1, -1, -1]]
            sage: L = Link([[1, 4, 2, 3], [6, 1, 3, 2], [7, 4, 8, 5], [5, 8, 6, 7]])
            sage: L.oriented_gauss_code()
            [[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]]
            sage: B = BraidGroup(8)
            sage: b = B([1, 1, 1, 1, 1])
            sage: L = Link(b)
            sage: L.oriented_gauss_code()
            [[[1, -2, 3, -4, 5, -1, 2, -3, 4, -5]], [1, 1, 1, 1, 1]]
        """
        if self._oriented_gauss_code is not None:
            return self._oriented_gauss_code

        pd = self.pd_code()
        orient = self.orientation()
        crossing_info = {}
        for i, j in enumerate(pd):
            if orient[i] == -1:
                crossing_info[(j[0], -1, i + 1)] = j[2]
                crossing_info[(j[3], 1, i + 1)] = j[1]
            elif orient[i] == 1:
                crossing_info[(j[0], -1, i + 1)] = j[2]
                crossing_info[(j[1], 1, i + 1)] = j[3]
        edges = {}
        cross_number = {}
        for i, j in crossing_info.items():
            edges[i[0]] = [j]
            if i[1] == 1:
                cross_number[i[0]] = i[2]
            elif i[1] == -1:
                cross_number[i[0]] = -i[2]
        edges_graph = DiGraph(edges)
        d = edges_graph.all_simple_cycles()
        code = []
        for i in d:
            l = []
            for j in i:
                l.append(cross_number[j])
            del l[-1]
            code.append(l)
        oriented_code = [code, orient]
        self._oriented_gauss_code = oriented_code
        return self._oriented_gauss_code

    def pd_code(self):
        """
        Return the planar diagram code of ``self``.

        The planar diagram is returned in the following format.

        We construct the crossing by starting with the entering component
        of the undercrossing, move in the clockwise direction and then
        generate the list. If the crossing is given by `[a, b, c, d]`,
        then we interpret this information as:

        1. `a` is the entering component of the undercrossing;
        2. `b, d` are the components of the overcrossing;
        3. `c` is the leaving component of the undercrossing.

        OUTPUT: planar diagram representation

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]], [1, 1, -1, -1]])
            sage: L.pd_code()
            [[6, 1, 7, 2], [2, 5, 3, 6], [8, 4, 1, 3], [4, 8, 5, 7]]
            sage: B = BraidGroup(2)
            sage: b = B([1, 1, 1, 1, 1])
            sage: L = Link(b)
            sage: L.pd_code()
            [[2, 1, 3, 4], [4, 3, 5, 6], [6, 5, 7, 8], [8, 7, 9, 10], [10, 9, 1, 2]]
            sage: L = Link([[[2, -1], [1, -2]], [1, 1]])
            sage: L.pd_code()
            [[2, 3, 1, 4], [4, 1, 3, 2]]
            sage: L = Link([[1, 2, 3, 3], [2, 4, 5, 5], [4, 1, 7, 7]])
            sage: L.pd_code()
            [[1, 2, 3, 3], [2, 4, 5, 5], [4, 1, 7, 7]]
        """
        if self._pd_code is not None:
            return self._pd_code

        if self._oriented_gauss_code is not None:
            oriented_gauss_code = self._oriented_gauss_code
            d_dic = {}
            if len(oriented_gauss_code[0]) > 1:
                d = flatten(oriented_gauss_code[0])
                for i, j in enumerate(d):
                    d_dic[j] = [i + 1, i + 2]
                # here we collect the final component in each Gauss code
                last_component = [i[-1] for i in oriented_gauss_code[0]]
                first_component = [i[0] for i in oriented_gauss_code[0]]
                # here we correct the last_component
                for i, j in zip(last_component, first_component):
                    d_dic[i][1] = d_dic[j][0]
                crossing_dic = {}
                for i,x in enumerate(oriented_gauss_code[1]):
                    if x == -1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][1],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][0]]
                    elif x == 1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][0],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][1]]
            elif len(oriented_gauss_code[0]) == 1:
                for i, j in enumerate(oriented_gauss_code[0][0]):
                    d_dic[j] = [i + 1, i + 2]
                d_dic[oriented_gauss_code[0][0][-1]][1] = 1
                crossing_dic = {}
                for i, x in enumerate(oriented_gauss_code[1]):
                    if x == -1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][1],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][0]]
                    elif x == 1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][0],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][1]]

            pd = crossing_dic.values()
            self._pd_code = pd
            return self._pd_code

        if self._braid is not None:
            strings = range(1, self._braid.strands() + 1)
            b = list(self._braid.Tietze())
            pd = []
            strings_max = strings[-1]
            for i in b:
                if i > 0:
                    pd.append(
                        [strings[i], strings[i - 1], strings_max + 1, strings_max + 2])
                else:
                    pd.append(
                        [strings[abs(i) - 1], strings_max + 1, strings_max + 2, strings[abs(i)]])
                strings[abs(i) - 1] = strings_max + 1
                strings[abs(i)] = strings_max + 2
                strings_max = strings_max + 2
            for i in pd:
                for j in range(4):
                    if i[j] in strings:
                        i[j] = strings.index(i[j]) + 1
            self._pd_code = pd
            return pd

        raise TypeError("BUG: invalid state")

    def gauss_code(self):
        """
        Return the Gauss code of ``self``.

        The Gauss code is generated by the following procedure:

        a. Number the crossings from `1` to `n`.
        b. Select a point on the knot and start moving along the component.
        c. At each crossing, take the number of the crossing, along with
           sign, which is `-` if it is an undercrossing and `+` if it is a
           overcrossing.

        OUTPUT: Gauss code representation

        EXAMPLES::

            sage: L = Link([[1,4,2,3],[4,1,3,2]])
            sage: L.gauss_code()
            [[-1, 2], [1, -2]]
            sage: B = BraidGroup(8)
            sage: L = Link(B([1, -2, 1, -2, -2]))
            sage: L.gauss_code()
            [[-1, 3, -4, 5], [1, -2, 4, -5, 2, -3]]
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L.gauss_code()
            [[-1, 2], [-3, 4], [1, 3, -4, -2]]
        """
        return self.oriented_gauss_code()[0]

    def _dowker_notation_(self):
        """
        Return the Dowker notation of ``self``.

        Similar to the PD code we number the components, so every crossing
        is represented by four numbers. We focus on the incoming entities
        of the under and the overcrossing. It is the pair of incoming
        undercrossing and the incoming overcrossing. This information at
        every crossing gives the Dowker notation.

        OUTPUT:

        A list containing the pair of incoming under cross and the incoming
        over cross.

        EXAMPLES::

            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L._dowker_notation_()
            [(1, 6), (7, 2), (3, 10), (11, 4), (14, 5), (13, 8), (12, 9)]
            sage: B = BraidGroup(4)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L._dowker_notation_()
            [(2, 1), (3, 6), (7, 5), (8, 10)]
            sage: L = Link([[1,4,2,3],[4,1,3,2]])
            sage: L._dowker_notation_()
            [(1, 3), (4, 2)]
        """
        pd = self.pd_code()
        orient = self.orientation()
        dn = [(i[0], i[3]) if orient[j] == -1 else (i[0], i[1])
              for j, i in enumerate(pd)]
        return dn

    def _braidwordcomponents_(self):
        """
        Return the disjoint braid components, if any, else return the braid
        of ``self``.

        For example consider the braid `[-1, 3, 1, 3]` this can be viewed
        as a braid with components as `[-1, 1]` and `[3, 3]`. There is no
        common crossing to these two (in sense there is a crossing between
        strand `1` and `2`, crossing between `3` and `4` but no crossing
        between strand `2` and `3`, so these can be viewed as independent
        components in the braid).

        OUTPUT: list containing the components is returned

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._braidwordcomponents_()
            ([-1, 1], [3, 3])
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braidwordcomponents_()
            ([-1, 1, 1, 1], [3], [5, 7, 6])
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braidwordcomponents_()
            ([-2, 1, 1], [4, 4], [6])
        """
        ml = list(self.braid().Tietze())
        if not ml:
            raise ValueError("the braid remains the same with no components")

        l = set(abs(k) for k in ml)
        missing1 = set(range(min(l), max(l) + 1)) - l
        if not missing1:
            return tuple([ml])

        missing = sorted(missing1)
        x = [[] for i in range(len(missing) + 1)]
        for i,a in enumerate(missing):
            for j, mlj in enumerate(ml):
                if mlj != 0 and abs(mlj) < a:
                    x[i].append(mlj)
                    ml[j] = 0
                elif mlj != 0 and abs(mlj) > missing[-1]:
                    x[-1].append(mlj)
                    ml[j] = 0
        return tuple([a for a in x if a])

    def _braidwordcomponentsvector_(self):
        """
        The list from the :meth:`_braidwordcomponents_` is flattened to
        give out the vector form.

        OUTPUT: vector containing braidwordcomponents

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._braidwordcomponentsvector_()
            [-1, 1, 3, 3]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braidwordcomponentsvector_()
            [-1, 1, 1, 1, 3, 5, 7, 6]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braidwordcomponentsvector_()
            [-2, 1, 1, 4, 4, 6]
        """
        bc = self._braidwordcomponents_()
        return flatten(bc)

    def _homology_generators_(self):
        """
        The set of generators for the first homology group of the connected
        Seifert surface of the given link.

        This method uses the :meth:`_braidwordcomponentsvector_` to generate
        the homology generators. The position of the repeated element w.r.t.
        the braidwordcomponentvector list is compiled into a list.

        OUTPUT:

        - The homology generators relating to the braid word representation

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._homology_generators_()
            [1, 0, 3]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._homology_generators_()
            [1, 2, 3, 0, 0, 0, 0]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._homology_generators_()
            [0, 2, 0, 4, 0]
        """
        x4 = self._braidwordcomponentsvector_()
        hom_gen = []
        for j in range(len(x4) - 1):
            a = abs(x4[j])
            for i in range(j + 1, len(x4)):
                if a == abs(x4[i]):
                    hom_gen.append(i)
                    break
            else:
                hom_gen.append(0)
        return hom_gen

    @cached_method
    def seifert_matrix(self):
        """
        Return the Seifert matrix associated with ``self``.

        OUTPUT:

        The intersection matrix of a (not necessarily minimal) Seifert surface.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.seifert_matrix()
            [ 0  0]
            [ 0 -1]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L.seifert_matrix()
            [ 0  0  0]
            [ 1 -1  0]
            [ 0  1 -1]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.seifert_matrix()
            [-1  0]
            [ 0 -1]
        """
        x5 = self._braidwordcomponentsvector_()
        h = self._homology_generators_()
        hl = len(h)
        A = matrix(ZZ, hl, hl)
        for i in range(hl):
            if h[i] != 0:
                for j in range(i, hl):
                    if i == j:
                        A[i, j] = -cmp((x5[i] + x5[h[i]]), 0)
                    elif (h[i] > h[j]):
                        A[i, j] = 0
                        A[j, i] = 0
                    elif (h[i] < j):
                        A[i, j] = 0
                        A[j, i] = 0
                    elif (h[i] == j):
                        if(x5[j] > 0):
                            A[i, j] = 0
                            A[j, i] = 1
                        else:
                            A[i, j] = -1
                            A[j, i] = 0
                    elif abs(abs(x5[i]) - abs(x5[j])) > 1:
                        A[i, j] = 0
                    elif (abs(x5[i]) - abs(x5[j]) == 1):
                        A[i, j] = 0
                        A[j, i] = -1
                    elif (abs(x5[j]) - abs(x5[i]) == 1):
                        A[i, j] = 1
                        A[j, i] = 0
                    else:  # for debugging
                        A[i, j] = 2
                        A[j, i] = 2
            else:
                for k in range(hl):
                    A[k, i] = 0
                    A[i, k] = 0
        k = []
        for i in range(hl):
            if h[i] == 0:
                k.append(i)
        for i in reversed(k):
            A = A.delete_rows([i])
            A = A.delete_columns([i])
        A.set_immutable()
        return A

    @cached_method
    def number_of_components(self):
        """
        Return the number of connected components of ``self``.

        OUTPUT: number of connected components

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.number_of_components()
            4
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.number_of_components()
            5
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.number_of_components()
            1
        """
        G = Graph()
        pd = self.pd_code()
        G.add_vertices(set(flatten(pd)))
        for c in pd:
            G.add_edge(c[0], c[2])
            G.add_edge(c[1], c[3])
        return G.connected_components_number()

    def is_knot(self):
        """
        Return ``True`` if ``self`` is a knot.

        Every knot is a link but the converse is not true.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([1,3,1,-3]))
            sage: L.is_knot()
            False
            sage: B = BraidGroup(8)
            sage: L = Link(B([1, 2, 3, 4, 5, 6]))
            sage: L.is_knot()
            True
        """
        return self.number_of_components() == 1

    def genus(self):
        """
        Return the genus of ``self``.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.genus()
            0
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.genus()
            0
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.genus()
            1
        """
        b = self.braid().Tietze()
        if not b:
            return 0

        B = self.braid().parent()
        x = self._braidwordcomponents_()
        q = []
        genus = 0
        s_tmp = []
        for xi in x:
            tmp = []
            b1 = min(abs(k) for k in xi)
            for xij in xi:
                if xij > 0:
                    xij = xij - b1 + 1
                else:
                    xij = xij + b1 - 1
                tmp.append(xij)
            s_tmp.append(B(tmp))
        s = []
        for i in s_tmp:
            b = i.Tietze()
            s.append(list(b))
        t = [Link(B(s[i])).number_of_components() for i in range(len(s))]
        for i, j in enumerate(s):
            if not j:
                s[i].append(-2)
        for i in s:
            q2 = max(abs(k) + 1 for k in i)
            q.append(q2)
        g = [((2 - t[i]) + len(x[i]) - q[i]) / 2 for i in range(len(x))]
        for i in range(len(g)):
            genus = genus + g[i]
        return genus

    def signature(self):
        """
        Return the signature of ``self``.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.signature()
            -1
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.signature()
            -2
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.signature()
            -2
        """
        m = 2 * (self.seifert_matrix() + self.seifert_matrix().transpose())
        e = m.eigenvalues()
        sum = 0
        s = []
        for i, j in enumerate(e):
            s.append(cmp(j, 0))
            sum = sum + s[i]
        return sum

    def alexander_polynomial(self, var='t'):
        """
        Return the Alexander polynomial of ``self``.

        INPUT:

        - ``var`` -- (default: ``'t'``) the variable in the polynomial

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.alexander_polynomial()
            0
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.alexander_polynomial()
            t^-1 - 2 + t
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.alexander_polynomial()
            t^-1 - 1 + t
        """
        R = LaurentPolynomialRing(ZZ, var)
        t = R.gen()
        f = (self.seifert_matrix() - t *
             (self.seifert_matrix().transpose())).determinant()
        if f != 0:
            exp = f.exponents()
            return t ** ((-max(exp) - min(exp)) / 2) * f
        return f

    def determinant(self):
        """
        Return the determinant of the knot.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 2, 1, 2]))
            sage: L.determinant()
            1
            sage: B = BraidGroup(8)
            sage: L = Link(B([2, 4, 2, 3, 1, 2]))
            sage: L.determinant()
            3
            sage: L = Link(B([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.determinant()
            65
        """
        if self.is_knot():
            a = self.alexander_polynomial()
            return Integer(abs(a(-1)))

        raise NotImplementedError("determinant implmented only for knots")

    def is_alternating(self):
        """
        Return ``True`` if the given knot diagram is alternating else
        returns ``False``.

        Alternating diagram implies every overcross is followed by an
        undercross or the vice-versa.

        We look at the Gauss code if the sign is alternating, ``True``
        is returned else the knot is not alternating ``False`` is returned.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, -1, -1, -1]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([1, -2, -1, 2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([-1, 3, 1,3, 2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([-1,2,-1,2]))
            sage: L.is_alternating()
            True
        """
        if self.is_knot():
            x = self.gauss_code()
            s = [cmp(i, 0) for i in x[0]]
            return (s == [(-1) ** (i + 1) for i in range(len(x[0]))]
                    or s == [(-1) ** i for i in range(len(x[0]))])
        else:
            return False

    def orientation(self):
        r"""
        Return the orientation of the crossings of the link diagram
        of ``self``.

        EXAMPLES::

            sage: L = Link([[1, 4, 5, 2], [3, 5, 6, 7], [4, 8, 9, 6], [7, 9, 10, 11], [8, 1, 13, 10], [11, 13, 2, 3]])
            sage: L.orientation()
            [-1, 1, -1, 1, -1, 1]
            sage: L = Link([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
            sage: L.orientation()
            [-1, -1, -1, -1, 1, -1, 1]
            sage: L = Link([[1, 2, 3, 3], [2, 4, 5, 5], [4, 1, 7, 7]])
            sage: L.orientation()
            [-1, -1, -1]
        """
        directions = self._directions_of_edges_()[0]
        orientation = []
        for C in self.pd_code():
            if C[0] == C[1] or C[2] == C[3]:
                orientation.append(-1)
            elif C[1] == C[2] or C[0] == C[3]:
                orientation.append(1)
            elif directions[C[1]] == C:
                orientation.append(-1)
            else:
                orientation.append(1)
        return orientation

    def seifert_circles(self):
        """
        Return the Seifert circles from the link diagram of ``self``.

        Seifert circles are the circles obtained by smoothing all crossings
        respecting the orientation of the segments.

        Each Seifert circle is represented as a list of the segments
        that form it.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.seifert_circles()
            [[1, 7, 5, 3], [2, 6], [4, 8]]
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.seifert_circles()
            [[1, 13, 9, 3, 15, 5, 11, 7], [2, 10, 6, 12], [4, 16, 8, 14]]
            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1,-1,-1,-1,1,-1,1]])
            sage: L.seifert_circles()
            [[1, 7, 3, 11, 5], [2, 8, 14, 6], [4, 12, 10], [9, 13]]
            sage: L = Link([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
            sage: L.seifert_circles()
            [[1, 7, 3, 11, 5], [2, 8, 14, 6], [4, 12, 10], [9, 13]]
            sage: L = Link([[[-1, 2, -3, 5], [4, -2, 6, -5], [-4, 1, -6, 3]], [-1, 1, 1, 1, -1, -1]])
            sage: L.seifert_circles()
            [[1, 11, 8], [2, 7, 12, 4, 5, 10], [3, 9, 6]]
            sage: B = BraidGroup(2)
            sage: L = Link(B([1, 1, 1]))
            sage: L.seifert_circles()
            [[1, 3, 5], [2, 4, 6]]
        """
        available_segments = set(flatten(self.pd_code()))
        result = []
        tails, heads = self._directions_of_edges_()
        while available_segments:
            a = available_segments.pop()
            if heads[a] == tails[a]:
                result.append([a])
            else:
                C = heads[a]
                par = []
                while not a in par:
                    par.append(a)
                    if tails[C[(C.index(a) + 1) % 4]] == C:
                        a = C[(C.index(a) + 1) % 4]
                    else:
                        a = C[(C.index(a) - 1) % 4]
                    if a in available_segments:
                        available_segments.remove(a)
                    C = heads[a]
                result.append(par)
        return result

    def regions(self):
        """
        Return the regions from the link diagram of ``self``.

        Regions are obtained always turning left at each crossing.

        Then the regions are represented as a list with the segments that form
        its boundary, with a sign deppending on the orientation of the segment
        as part of the boundary.

        EXAMPLES::

            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L.regions()
            [[1, 7, 3, 11, 5], [2, -7], [4, -11], [6, -1], [8, -13, 10, -3], [9, 13], [12, -9, 14, -5], [-14, -8, -2, -6], [-12, -4, -10]]
            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.regions()
            [[1, 7, -4], [2, -5, -7], [3, -8, 5], [4, 8], [6, -1, -3], [-2, -6]]
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.regions()
            [[1, 13, -8], [2, -9, -13], [3, -14, 9], [4, 16, 8, 14], [5, 11, 7, -16], [6, -11], [10, -5, -15, -3], [12, -1, -7], [15, -4], [-12, -6, -10, -2]]
            sage: B = BraidGroup(2)
            sage: L = Link(B([-1, -1, -1]))
            sage: L.regions()
            [[1, 3, 5], [2, -1], [4, -3], [6, -5], [-2, -6, -4]]
            sage: L = Link([[[1, -2, 3, -4], [-1, 5, -3, 2, -5, 4]], [-1, 1, 1, -1, -1]])
            sage: L.regions()
            [[1, -5], [2, -8, 4, 5], [3, 8], [6, -9, -2], [7, -3, 9], [10, -4, -7], [-10, -6, -1]]
            sage: L = Link([[1, 2, 3, 3], [2, 5, 4, 4], [5, 7, 6, 6], [7, 1, 8, 8]])
            sage: L.regions()
            [[-3], [-4], [-6], [-8], [1, 2, 5, 7], [-2, 3, -1, 8, -7, 6, -5, 4]]

        .. NOTE::

            The link diagram is assumed to have only one completely isolated
            component. This is because otherwise some regions would be have
            disconnected boundary.
        """
        pd = self.pd_code()
        tails, heads = self._directions_of_edges_()
        available_edges = set(flatten(pd))
        if len(pd) == 1:
            if pd[0][0] == pd[0][1]:
                return [[-pd[0][2]], [pd[0][0]], [pd[0][2], -pd[0][0]]]
            else:
                return [[pd[0][2]], [-pd[0][0]], [-pd[0][2], pd[0][0]]]

        loops = [i for i in available_edges if heads[i] == tails[i]]
        available_edges = available_edges.union({-i for i in available_edges})
        regions = []

        for edge in loops:
            cros = heads[edge]
            if cros[1] == edge:
                regions.append([edge])
            else:
                regions.append([-edge])
            available_edges.remove(edge)
            available_edges.remove(-edge)

        while available_edges:
            edge = available_edges.pop()
            region = []
            while not edge in region:
                region.append(edge)
                if edge > 0 :
                    cros = heads[edge]
                    ind = cros.index(edge)
                else:
                    cros = tails[-edge]
                    ind = cros.index(-edge)
                next_edge = cros[(ind + 1) % 4]
                if [next_edge] in regions:
                    region.append(-next_edge)
                    next_edge = cros[(ind - 1) % 4]
                elif [-next_edge] in regions:
                    region.append(next_edge)
                    next_edge = cros[(ind - 1) % 4]
                if tails[next_edge] == cros:
                    edge = next_edge
                else:
                    edge = -next_edge
                if edge in available_edges:
                    available_edges.remove(edge)
            regions.append(region)
        return regions

    def writhe(self):
        """
        Return the writhe of ``self``.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.writhe()
            0
            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L.writhe()
            -3
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.writhe()
            -2
        """
        x = self.oriented_gauss_code()
        pos = x[1].count(1)
        neg = (-1) * x[1].count(-1)
        return pos + neg

    @lazy_attribute
    def _jones_polynomial(self):
        """
        Cached version of the Jones polynomial of the trace closure of the
        braid representation of ``self`` in a generic variable with the skein
        normalization.

        The computation of the Jones polynomial uses the representation
        of the braid group on the Temperley--Lieb algebra. We cache the
        part of the calculation which does not depend on the choices of
        variables or normalizations.

        .. SEEALSO::

            :meth:`jones_polynomial`

        TESTS::

            sage: B = BraidGroup(9)
            sage: b = B([1, 2, 3, 4, 5, 6, 7, 8])
            sage: Link(b).jones_polynomial()
            1

            sage: B = BraidGroup(2)
            sage: b = B([])
            sage: Link(b)._jones_polynomial
            -A^-2 - A^2
            sage: b = B([-1, -1, -1])
            sage: Link(b)._jones_polynomial
            -A^-16 + A^-12 + A^-4
        """
        braid = self.braid()
        trace = braid.markov_trace(normalized=False)
        A = trace.parent().gens()[0]
        D = A**2 + A**(-2)
        exp_sum = braid.exponent_sum()
        num_comp = braid.components_in_closure()
        return (-1)**(num_comp-1) * A**(2*exp_sum) * trace // D

    def jones_polynomial(self, variab=None, skein_normalization=False, algorithm='jonesrep'):
        """
        Return the Jones polynomial of ``self``.

        The normalization is so that the unknot has Jones polynomial `1`. If
        ``skein_normalization`` is ``True``, the variable of the result is
        replaced by a itself to the power of `4`, so that the result
        agrees with the conventions of [Lic]_ (which in particular differs
        slightly from the conventions used otherwise in this class), had
        one used the conventional Kauffman bracket variable notation directly.

        If ``variab`` is ``None`` return a polynomial in the variable `A`
        or `t`, depending on the value ``skein_normalization``. In
        particular, if ``skein_normalization`` is ``False``, return the
        result in terms of the variable `t`, also used in [Lic]_.

        The calculation goes through one of two possible algorithms, depending
        on the value of ``algorithm``. Possible values are ``'jonesrep'`` which
        uses the Jones representation of a braid representation of ``self`` to
        compute the polynomial of the trace closure of the braid, and
        ``statesum`` which recursively computes the Kauffman bracket of
        ``self``. Depending on how the link is given, there might be
        significant time gains in using one over the other. When the trace
        closure of the braid is ``self``, the algorithms give the same result.

        INPUT:

        - ``variab`` -- variable (default: ``None``); the variable in the
          resulting polynomial; if unspecified, use either a default variable
          in `ZZ[A,A^{-1}]` or the variable `t` in the symbolic ring

        - ``skein_normalization`` -- boolean (default: ``False``); determines
          the variable of the resulting polynomial

        - ``algorithm`` -- string (default: ``'jonesrep'``); algorithm to use.

        OUTPUT:

        If ``skein_normalization`` if ``False``, this returns an element
        in the symbolic ring as the Jones polynomial of the link might
        have fractional powers when the link is not a knot. Otherwise the
        result is a Laurant polynomial in ``variab``.

        EXAMPLES:

        The unknot::

            sage: B = BraidGroup(9)
            sage: b = B([1, 2, 3, 4, 5, 6, 7, 8])
            sage: Link(b).jones_polynomial()
            1

        Two unlinked unknots::

            sage: B = BraidGroup(4)
            sage: b = B([1, 3])
            sage: Link(b).jones_polynomial()
            -sqrt(t) - 1/sqrt(t)

        The Hopf link::

            sage: B = BraidGroup(2)
            sage: b = B([-1,-1])
            sage: Link(b).jones_polynomial()
            -1/sqrt(t) - 1/t^(5/2)

        Different representations of the trefoil and one of its mirror::

            sage: B = BraidGroup(2)
            sage: b = B([-1, -1, -1])
            sage: Link(b).jones_polynomial(skein_normalization=True)
            -A^-16 + A^-12 + A^-4
            sage: Link(b).jones_polynomial()
            1/t + 1/t^3 - 1/t^4
            sage: B = BraidGroup(3)
            sage: b = B([-1, -2, -1, -2])
            sage: Link(b).jones_polynomial(skein_normalization=True)
            -A^-16 + A^-12 + A^-4
            sage: R.<x> = LaurentPolynomialRing(GF(2))
            sage: Link(b).jones_polynomial(skein_normalization=True, variab=x)
            x^-16 + x^-12 + x^-4
            sage: B = BraidGroup(3)
            sage: b = B([1, 2, 1, 2])
            sage: Link(b).jones_polynomial(skein_normalization=True)
            A^4 + A^12 - A^16

        K11n42 (the mirror of the "Kinoshita-Terasaka" knot) and K11n34 (the
        mirror of the "Conway" knot)::

            sage: B = BraidGroup(4)
            sage: K11n42 = Link(B([1, -2, 3, -2, 3, -2, -2, -1, 2, -3, -3, 2, 2]))
            sage: K11n34 = Link(B([1, 1, 2, -3, 2, -3, 1, -2, -2, -3, -3]))
            sage: cmp(K11n42.jones_polynomial(), K11n34.jones_polynomial())
            0

        The two algorithms for computation give the same result when the trace
        closure of the braid representation is the link itself::

            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6, -5]], [-1, -1, -1, -1, 1, -1, 1]])
            sage: jonesrep = L.jones_polynomial(algorithm='jonesrep')
            sage: statesum = L.jones_polynomial(algorithm='statesum')
            sage: cmp(jonesrep, statesum)
            0

        But when we have thrown away unknots so that the trace closure of the
        braid is not necessarily the link itself, this is only true up to a
        power of the Jones polynomial of the unknot::

            sage: B = BraidGroup(3)
            sage: b = B([1])
            sage: L = Link(b)
            sage: b.components_in_closure()
            2
            sage: L.number_of_components()
            1
            sage: L.jones_polynomial(algorithm='jonesrep')
            -sqrt(t) - 1/sqrt(t)
            sage: L.jones_polynomial(algorithm='statesum')
            1

        REFERENCES:

        .. [Lic] William B. Raymond Lickorish. An Introduction to Knot Theory,
           volume 175 of Graduate Texts in Mathematics. Springer-Verlag,
           New York, 1997. ISBN 0-387-98254-X
        """
        if algorithm == 'statesum':
            poly = self._bracket_()
            t = poly.parent().gens()[0]
            writhe = self.writhe()
            jones = (poly * (-t)**(-3 * writhe))
            # Switch to the variable A to have the result agree with the output
            # of the jonesrep algorithm
            A = LaurentPolynomialRing(ZZ, 'A').gen()
            jones = jones(A**-1)
        elif algorithm == 'jonesrep':
            jones = self._jones_polynomial
        else:
            raise ValueError("bad value of algorithm")

        if skein_normalization:
            if variab is None:
                return jones
            else:
                return jones(variab)
        else:
            if variab is None:
                variab = 't'
            # We force the result to be in the symbolic ring because of the expand
            return jones(SR(variab)**(ZZ(1)/ZZ(4))).expand()

    @cached_method
    def _bracket_(self):
        r"""
        Return the Kauffman bracket polynomial of the diagram.

        Note that this is not an invariant of the link, but of the diagram.
        In particular, it is not invariant under Reidemeister I moves.

        EXAMPLES::

            sage: L = Link([[[-1, 2, 3, -4, 5, -6, 7, 8, -2, -5, 6, 1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L._bracket_()
            -t^-10 + 2*t^-6 - t^-2 + 2*t^2 - t^6 + t^10 - t^14
            sage: L = Link([[2, 1, 3, 4], [4, 3, 1, 2]])
            sage: L._bracket_()
            -t^-4 - t^4
        """
        t = LaurentPolynomialRing(ZZ, 't').gen()
        pd_code = self.pd_code()
        if len(pd_code) == 1:
            if pd_code[0][0] == pd_code[0][1]:
                return -t**(-3)
            else:
                return -t**3

        cross = pd_code[0]
        rest = deepcopy(pd_code[1:])
        [a, b, c, d] = cross
        if a == b and c == d and len(rest) > 0:
            return (~t + t**(-5)) * Link(rest)._bracket_()
        elif a == d and c == b and len(rest) > 0:
            return (t + t**5) * Link(rest)._bracket_()
        elif a == b:
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = c
            return -t**(-3) * Link(rest)._bracket_()
        elif a == d:
            for cross in rest:
                if c in cross:
                    cross[cross.index(c)] = b
            return -t**3 * Link(rest)._bracket_()
        elif c == b:
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = a
            return -t**3 * Link(rest)._bracket_()
        elif c == d:
            for cross in rest:
                if b in cross:
                    cross[cross.index(b)] = a
            return -t**(-3) * Link(rest)._bracket_()
        else:
            rest_2 = deepcopy(rest)
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = a
                if c in cross:
                    cross[cross.index(c)] = b
            for cross in rest_2:
                if d in cross:
                    cross[cross.index(d)] = c
                if b in cross:
                    cross[cross.index(b)] = a
            return t * Link(rest)._bracket_() + ~t * Link(rest_2)._bracket_()

    def _isolated_components_(self):
        r"""
        Return the PD codes of the isolated components of ``self``.

        Isolated components are links corresponding to subdiagrams that don't
        have any common crossing.

        EXAMPLES::

            sage: L = Link([[1, 1, 2, 2], [3, 3, 4, 4]])
            sage: L._isolated_components_()
            [[[1, 1, 2, 2]], [[3, 3, 4, 4]]]
        """
        G = Graph()
        for c in self.pd_code():
            G.add_vertex(tuple(c))
        for i in range(G.num_verts()-1):
            for j in range(i, G.num_verts()):
                if len(set(G.vertices()[i]).intersection(G.vertices()[j])) > 0:
                    G.add_edge(G.vertices()[i], G.vertices()[j])
        return [[list(i) for i in j] for j in G.connected_components()]

    def plot(self, gap=0.1, **kwargs):
        r"""
        Plot ``self``.

        INPUT:

        - ``gap`` -- (default: 0.1) the size of the blank gap left for
          the crossings

        The usual keywords for plots can be used here too.

        EXAMPLES::

            sage: L = Link([[[-1, 2, -3, 1, -2, 3], [4, -5, 6, -4, 5, -6]], [1, 1, 1, 1, 1, 1]])
            sage: L.plot()
            Graphics object consisting of 28 graphics primitives
        """
        comp = self._isolated_components_()
        if len(comp) > 1:
            L1 = Link(comp[0])
            L2 = Link(flatten(comp[1:], max_level=1))
            P1 = L1.plot(gap, **kwargs)
            P2 = L2.plot(gap, **kwargs)
            xtra = P1.get_minmax_data()['xmax'] + P2.get_minmax_data()['xmin'] + 2
            for P in P2:
                if hasattr(P, 'path'):
                    for p in P.path[0]:
                        p[0] += xtra
                    for p in P.vertices:
                        p[0] += xtra
                else:
                    P.xdata = [p + xtra for p in P.xdata]
            return P1 + P2
        if not 'color' in kwargs:
            kwargs['color'] = 'blue'
        if not 'axes' in kwargs:
            kwargs['axes'] = False
        if not 'aspect_ratio' in kwargs:
            kwargs['aspect_ratio'] = 1
        # The idea is the same followed in spherogram, but using MLP instead of
        # network flows.
        # We start by computing a way to bend the edges left or right
        # such that the resulting regions are in fact closed regions
        # with straight angles, and using the minimal number of bends.
        regions = sorted(self.regions(), key=len)
        regions = regions[:-1]
        edges = list(set(flatten(self.pd_code())))
        edges.sort()
        MLP = MixedIntegerLinearProgram(maximization = True)
        # v will be the list of variables in the MLP problem. There will be
        # two variables for each edge: number of right bendings and number of
        # left bendings (at the end, since we are minimizing the total, only one
        # of each will be nonzero
        v = MLP.new_variable(nonnegative=True)
        for i in range(2*len(edges)):
            MLP.set_min(v[i], 0)
        # one condition for each region
        for i in range(len(regions)):
            cond = 0
            r = regions[i]
            es = 4 - len(r)
            for e in r:
                if e > 0:
                    cond = cond + v[2*edges.index(e)] - v[2*edges.index(e) + 1]
                else:
                    cond = cond - v[2*edges.index(-e)] + v[2*edges.index(-e) + 1]
            MLP.add_constraint(cond, min=es, max=es)
        MLP.set_objective(-sum(v.values()))
        MLP.solve()
        # we store the result in a vector s packing right bends as negative left ones
        s = range(len(edges))
        values = MLP.get_values(v)
        for i in range(len(edges)):
            s[i] = int(values[2*i] - values[2*i + 1])
        # segments represents the different parts of the previos edges after bending
        segments = {e: [(e,i) for i in range(abs(s[edges.index(e)])+1)] for e in edges}
        pieces = {tuple(i): [i] for j in segments.values() for i in j}
        nregions = []
        for r in regions:
            nregion = []
            for e in r:
                if e > 0:
                    rev = segments[e][:-1]
                    sig = sign(s[edges.index(e)])
                    nregion += [[a, sig] for a in rev]
                    nregion.append([segments[e][-1], 1])
                else:
                    rev = segments[-e][1:]
                    rev.reverse()
                    sig = sign(s[edges.index(-e)])
                    nregion+=[[a, -sig] for a in rev]
                    nregion.append([segments[-e][0], 1])
            nregions.append(nregion)
        N = max(segments.keys()) + 1
        segments = [i for j in segments.values() for i in j]
        badregions = [nr for nr in nregions if any(-1 == x[1] for x in nr)]
        while len(badregions)>0:
            badregion = badregions[0]
            badturns = []
            a = 0
            while badregion[a][1] != -1:
                a += 1
            c = -1
            b = a
            while c != 2:
                if b == len(badregion)-1:
                    b = 0
                else:
                    b += 1
                c += badregion[b][1]
            otherregion = [nr for nr in nregions
                           if any(badregion[b][0] == x[0] for x in nr)]
            if len(otherregion) == 1:
                otherregion = None
            elif otherregion[0] == badregion:
                otherregion = otherregion[1]
            else:
                otherregion = otherregion[0]
            N1 = N
            N = N + 2
            N2 = N1 + 1
            segments.append(N1)
            segments.append(N2)
            if type(badregion[b][0]) in (int, Integer):
                segmenttoadd = [x for x in pieces.keys()
                                if badregion[b][0] in pieces[x]]
                if len(segmenttoadd) > 0:
                    pieces[segmenttoadd[0]].append(N2)
            else:
                pieces[tuple(badregion[b][0])].append(N2)

            if a < b:
                r1 = badregion[:a] + [[badregion[a][0],0], [N1,1]] + badregion[b:]
                r2 = badregion[a+1:b] + [[N2,1],[N1,1]]
            else:
                r1 = badregion[b:a] + [[badregion[a][0],0], [N1,1]]
                r2 = badregion[:b] + [[N2,1],[N1,1]] + badregion[a+1:]

            if otherregion:
                c = [x for x in otherregion if badregion[b][0] == x[0]]
                c = otherregion.index(c[0])
                otherregion.insert(c+1, [N2,otherregion[c][1]])
                otherregion[c][1] = 0
            nregions.remove(badregion)
            nregions.append(r1)
            nregions.append(r2)
            badregions = [nr for nr in nregions if any(x[1] == -1 for x in nr)]
        MLP = MixedIntegerLinearProgram(maximization = True)
        variables = {}
        for e in segments:
            variables[e] = MLP.new_variable(nonnegative=True)
            MLP.set_min(variables[e][0], 1)
        for r in nregions:
            horp = []
            horm = []
            verp = []
            verm = []
            direction = 0
            for se in r:
                if direction % 4 == 0:
                    horp.append(variables[se[0]][0])
                elif direction == 1:
                    verp.append(variables[se[0]][0])
                elif direction == 2:
                    horm.append(variables[se[0]][0])
                elif direction == 3:
                    verm.append(variables[se[0]][0])
                if se[1] == 1:
                    direction += 1
            MLP.add_constraint(sum(horp)-sum(horm), min=0, max=0)
            MLP.add_constraint(sum(verp)-sum(verm), min=0, max=0)
        MLP.set_objective(-sum([x[0] for x in variables.values()]))
        solved = MLP.solve()
        lengths = {piece: sum(MLP.get_values(variables[a])[0] for a in pieces[piece])
                   for piece in pieces}
        image = line([], **kwargs)
        crossings = {tuple(self.pd_code()[0]): (0,0,0)}
        availables = self.pd_code()[1:]
        used_edges = []
        horizontal_eq = 0
        vertical_eq = 0
        ims = line([], **kwargs)
        while len(used_edges) < len(edges):
            i = 0
            j = 0
            while crossings.keys()[i][j] in used_edges:
                if j < 3:
                    j += 1
                else:
                    j = 0
                    i+=1
            c = crossings.keys()[i]
            e = c[j]
            used_edges.append(e)
            direction = (crossings[c][2] - c.index(e)) % 4
            orien = self.orientation()[self.pd_code().index(list(c))]
            if s[edges.index(e)] < 0:
                turn = -1
            else:
                turn = 1
            lengthse = [lengths[(e,i)] for i in range(abs(s[edges.index(e)])+1)]
            if c.index(e) == 0 or (c.index(e) == 1 and orien == 1) or (c.index(e) == 3 and orien == -1):
                turn = -turn
                lengthse.reverse()
            tailshort = (c.index(e) % 2 == 0)
            x0 = crossings[c][0]
            y0 = crossings[c][1]
            im = []
            for l in lengthse:
                if direction == 0:
                    x1 = x0 + l
                    y1 = y0
                elif direction == 1:
                    x1 = x0
                    y1 = y0 + l
                elif direction == 2:
                    x1 = x0 - l
                    y1 = y0
                elif direction == 3:
                    x1 = x0
                    y1 = y0 -l
                im.append(([[x0,y0],[x1,y1]], l, direction))
                direction = (direction + turn) % 4
                x0 = x1
                y0 = y1
            direction = (direction - turn) % 4
            c2 = [ee for ee in availables if e in ee]
            if len(c2) == 1:
                availables.remove(c2[0])
                crossings[tuple(c2[0])] = (x1, y1, (direction + c2[0].index(e) + 2) % 4)
            c2 = [ee for ee in self.pd_code() if e in ee and ee != list(c)]
            if not c2:
                headshort = not tailshort
            else:
                headshort = (c2[0].index(e) % 2 == 0)
            a = deepcopy(im[0][0])
            b = deepcopy(im[-1][0])
            if tailshort:
                im[0][0][0][0] += cmp(a[1][0], im[0][0][0][0]) * gap
                im[0][0][0][1] += cmp(a[1][1], im[0][0][0][1]) * gap
            if headshort:
                im[-1][0][1][0] -= cmp(b[1][0], im[-1][0][0][0]) * gap
                im[-1][0][1][1] -= cmp(b[1][1], im[-1][0][0][1]) * gap
            l = line([], **kwargs)
            c = 0
            p = im[0][0][0]
            if len(im) == 4 and max([x[1] for x in im]) == 1:
                l = bezier_path([[im[0][0][0], im[0][0][1], im[-1][0][0], im[-1][0][1]]], **kwargs)
                p = im[-1][0][1]
            else:
                while c < len(im)-1:
                    if im[c][1] > 1:
                        (a, b) = im[c][0]
                        if b[0] > a[0]:
                            e = [b[0] - 1, b[1]]
                        elif b[0] < a[0]:
                            e = [b[0] + 1, b[1]]
                        elif b[1] > a[1]:
                            e = [b[0], b[1] - 1]
                        elif b[1] < a[1]:
                            e = [b[0] , b[1] + 1]
                        l += line((p, e), **kwargs)
                        p = e
                    if im[c+1][1] == 1 and c < len(im) - 2:
                        xr = round(im[c+2][0][1][0])
                        yr = round(im[c+2][0][1][1])
                        xp = xr - im[c+2][0][1][0]
                        yp = yr - im[c+2][0][1][1]
                        q = [p[0] + im[c+1][0][1][0] - im[c+1][0][0][0] - xp,
                             p[1] + im[c+1][0][1][1] - im[c+1][0][0][1] - yp]
                        l += bezier_path([[p, im[c+1][0][0], im[c+1][0][1], q]], **kwargs)
                        c += 2
                        p = q
                    else:
                        if im[c+1][1] == 1:
                            q = im[c+1][0][1]
                        else:
                            q = [im[c+1][0][0][0] + sign(im[c+1][0][1][0] - im[c+1][0][0][0]),
                                 im[c+1][0][0][1] + sign(im[c+1][0][1][1] - im[c+1][0][0][1])]
                        l += bezier_path([[p, im[c+1][0][0], q]], **kwargs)
                        p = q
                        c += 1
            l += line([p, im[-1][0][1]], **kwargs)
            image += l
            ims += sum(line(a[0], **kwargs) for a in im)
        return image

