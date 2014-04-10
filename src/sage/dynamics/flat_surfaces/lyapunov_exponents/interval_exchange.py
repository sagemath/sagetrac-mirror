r"""
This file define operations on interval exchange, and associated translation surfaces.
In particular, we give a way to manipulate finite covers of theses translations surfaces.
For example, to compute their lyapunov exponents, and its decomposition in isotropic subspaces.

AUTHORS:

- Charles Fougeron (2014-04-08) <charles.fougeron@ens.fr> : initial version

.. TODO::

    - Merge to the existing code concerning strata and interval echanges
    - Fix a bug for zero lyapunov exponents in isotropic decomposition
"""



#Here we give one function to get a random length vector for the IE
from sage.misc.cachefunc import cached_method
from path_vector import *

def UCF_mat_to_cyclotomic_field(M):
    from sage.rings.number_field.number_field import CyclotomicField
    from sage.rings.arith import lcm
    n = lcm([x.field_order() for x in M.list()])
    return M.change_ring(CyclotomicField(n))

class Interval(object) :
    def __init__(self, label, orientation):
        assert orientation == 1 or orientation == -1
        self.label = label
        self.orientation = orientation
    def __repr__(self):
        if self.orientation == 1:
            return " %s"%(self.label)
        else:
            return "-%s"%(self.label)
    def invert(self):
        if self.orientation == 1:
            self.orientation = -1
        else:
            self.orientation = 1

def index_cycle(i, cycles):
    """
    return index of the cycle in the array cycles where i appears
    """
    d = len(cycles)
    cursor = 0
    while cycles[cursor].count(t) == 0 and cursor <= d:
        cursor += 1
    if cursor == d:
        raise NameError('Index not in the cycle')
    return cursor

class CoveringStratum(FlippedLabelledPermutationLI):
    r"""
    Test
    """
    def __init__(self, intervals=None, alphabet=None, flips=None, permutations=None, lengths=None):
        if intervals is None:
            intervals = [[], []]
        if flips is None: flips = []

        super(FlippedLabelledPermutationLI, self).__init__(intervals, alphabet)
        self._init_flips(intervals, flips)

        self._permutations= deepcopy(permutations)
        self.degree = len(permutations[permutations.keys()[0]])
        self.number_of_intervals = [len(intervals[0]), len(intervals[1])]
        self._lengths = deepcopy(lengths)

        self.total_number_of_intervals_base = len(intervals[0]) + len(intervals[1])
        self.nb_labels = self.total_number_of_intervals_base / 2
        self.total_number_of_intervals_cover = (len(intervals[0]) + len(intervals[1]))*self.degree

class IntExchange(object):
    r"""
    Create an interval exchange with possibility of considering a finite cover by giving
    the permutations associated to the natural generator of the  homotopy group of the suspension surface.

    We have a convention for labelling the saddle connections of the cover of the suspension surface :
    first we define a directed path between the two saddle connections, which will be a cycle in the ground surface,
    the direction will be choosen such that the intersection number with both saddle connections, directed according
    to the orientation is positive (egal to 1).
    The label bearing the number of the copy of the cover will be the one to which the path arrives.


    INPUT:

    - ``intervals`` -- a list of two lists of Interval class elements

    - ``lengths`` -- a dictionnary of positive floats for each label appearing in the intervals;
    the length of intervals of the corresponding label

    - ``permutations`` -- dictionnary of permutation for each label appearing in the intervals;
    The permutation associated to the label `a` gives for the saddle connection of label `a`
    bearing the number `i`, the number of the other copy on which there is the saddle connection
    to which it is identified.


    EXAMPLES::

        sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
        ...   [Interval('b',-1), Interval('c', 1), Interval('c',-1)]]

        sage: sigma_a = Permutation('(1,2)(3,4)(5,6)')
        sage: sigma_b = Permutation('(1,3)(2,4)(5)(6)')
        sage: sigma_c = Permutation('(1,6)(2,5)(3,4)')

        sage: permutations = {'a': sigma_a, 'b': sigma_b, 'c': sigma_c}

        sage: ie = IntExchange(intervals, None, permutations).rand_lg()

        sage: ie
        Intervals    : [ a, -a,  b]
                       [ -b,  c, -c]
        Lengths      : {'a': 0.111439293741037, 'c': 0.111439293741037, 'b': 0.5143475134191677}
        Permutations : {'a': [2, 1, 4, 3, 6, 5], 'c': [6, 5, 4, 3, 2, 1], 'b': [3, 4, 1, 2, 5, 6]}

    """
    def __init__(self, intervals, permutations = None, lengths = None):
        from copy import deepcopy
        self._intervals = deepcopy(intervals)
        self._permutations= deepcopy(permutations)
        self.degree = len(permutations[permutations.keys()[0]])
        self.number_of_intervals = [len(intervals[0]), len(intervals[1])]
        self._lengths = deepcopy(lengths)

        self._twins = [[self.twin(i, j) for j in range(self.number_of_intervals[i])] for i in range(2)]

        self.total_number_of_intervals_base = len(intervals[0]) + len(intervals[1])
        self.nb_labels = self.total_number_of_intervals_base / 2
        self.total_number_of_intervals_cover = (len(intervals[0]) + len(intervals[1]))*self.degree
        self.cycles = self._edge_cover_permutation().to_cycles()
        self.n_cycles = len(self.cycles)
        #self.signature_intersection = self.signature(self.intersection_matrix())

    def __repr__(self):
        return ("Intervals    : %s\n               %s\nPermutations : %s\n"
                %(self._intervals[0], self._intervals[1], self._permutations))

    def twin(self, i, j):
        r"""
        Give the position of the twin interval associated to the one in line i and place j
        
        EXAMPLES::
        
        """
        for n in range(self.number_of_intervals[i]):
            if n <> j and self._intervals[i][n].label == self._intervals[i][j].label:
                return(i,n)
        for n in range(self.number_of_intervals[(i+1)%2]):
            if self._intervals[(i+1)%2][n].label == self._intervals[i][j].label:
                return((i+1)%2,n)

    def degrees(self):
        r"""
        Give the list of numbers given to each copy of the suface in the cover.
        It goes from 1 to the degree of the cover.

        
        """
        return(range(1, self.degree + 1))

    def iter_lab(self):
        r"""
        Gives the list of the labels of intervals in order to iterate
        """
        return(self._permutations.iterkeys())

    def labels(self):
        r"""
        Gives the list of the labels of intervals
        """
        return(self._permutations.keys())

    def canonical_VectPaths(self):
        return(VectPaths(self.labels(), self.degrees()))

    def cover_generators(self):
        return([(d,a) for d in self.degrees() for a in self.labels()])

    @cached_method
    def _free_basis(self):
        r"""
        Return matrix from free basis of H1(X, \Sigma) to R^(d*lab)
        and its pseudo inverse, from R^(d*lab) to H1(X, \Sigma)
        """
        bound_equation = {d: self.canonical_VectPaths() for d in self.degrees()}
        for d in self.degrees():
            for i in xrange(2):
                for j in xrange(self.number_of_intervals[i]):
                    bound_equation[d]._vect[self.name(i,j)(d)][self.label(i,j)] += self._intervals[i][j].orientation
        s = []
        def bind(ind,deg,lab):
            from copy import deepcopy
            line = deepcopy(bound_equation[ind])
            c = line.val(deg, lab)
            for i in self.degrees():
                coeff = bound_equation[i].val(deg,lab)
                for d,a in self.cover_generators():
                    bound_equation[i]._vect[d][a] -= line.val(d,a)*coeff/c
                    if projections[d][a].val(deg, lab) <> 0:
                        coeff_aux = projections[d][a].val(deg, lab)
                        for d_aux, a_aux in self.cover_generators():
                            projections[d][a]._vect[d_aux][a_aux] -= line.val(d_aux,a_aux)*coeff_aux/c
        def non_zero(deg, lab):
            for i in self.degrees():
                if bound_equation[i].val(deg,lab) <> 0:
                    return i
            return None
        projections = {d: {a: self.canonical_VectPaths().id(d,a) for a in self.labels()} for d in self.degrees()}
        for d,a in self.cover_generators():
            i = non_zero(d,a)
            if i == None:
                s.append((d,a))
            else:
                bind(i,d,a)
        def vect_paths_id((d,a)):
            return self.canonical_VectPaths().id(d,a).to_list()
        def eval(vect, (d,a)):
            return vect.val(d,a)
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ
        return matrix(ZZ, [vect_paths_id(s[i]) for i in xrange(len(s))]), matrix(ZZ, [[eval(projections[d][a],s[i]) for i in xrange(len(s))]
                                                                             for d, a in self.cover_generators()])
    @cached_method
    def galois_group(self):
        r"""
        Return the galois group of the cover
        """
        n = self.degree
        gap_list = ""
        for lab in self.iter_lab():
            perm = self._permutations[lab]
            perm_str = perm.cycle_string()
            gap_list += perm_str + ", "
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        from sage.interfaces.gap import gap
        G = gap("Centralizer(SymmetricGroup(%s), "%n + "Group([" + gap_list + "]))")
        return G

    @cached_method
    def _characters(self):
        r"""
        We want to decompose homology space with the observation that it gives a representation of the galois group.

        OUTPUT:
            - character: table of character, character[i][g] give the value of the i-th character on g
            g is given by a number
            - character_degree: list of degree of every character
            - g : Order of the group
            - perm : table s.t. perm[g] give the permutation associated to the group element g on the cover
            - n_cha : number of characters
        """
        from sage.interfaces.gap import gap
        G = self.galois_group()
        G_order, T = gap.Order(G)._sage_(), gap.CharacterTable(G)
        irr_characters = gap.Irr(T)
        n_cha = len(irr_characters)
        from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import UniversalCyclotomicField
        UCF = UniversalCyclotomicField()
        character = [[UCF(irr_characters[i][j]._sage_()) for j in range(1, G_order + 1)] for i in range(1, n_cha + 1)]
        character_degree = [gap.Degree(irr_characters[i])._sage_() for i in range(1, n_cha + 1)]
        gap_size_centralizers = gap.SizesCentralizers(T)
        gap_orders = gap.OrdersClassRepresentatives(T)
        def find_group_element(t):
            cursor = 1
            while gap.Order( gap.Centralizer(gap.SymmetricGroup(n), t) ) <> gap_size_centralizers[cursor] and gap.Order(t) <> gap_orders[cursor]:
                cursor += 1
            return cursor
        elements_group = gap.Elements(G)
        from sage.combinat.permutation import Permutation, Permutations
        perm = [Permutation(str(elements_group[i]))*Permutations(self.degree).identity() for i in range(1, G_order + 1)]
        #identity assures that the permutation has the right size
        return(character, character_degree, G_order, perm, n_cha)

    def character_table(self):
        return self._characters()[0]

    def character_degree(self):
        return self._characters()[1]

    def galois_group_order(self):
        return self._characters()[2]

    def galois_group_permutation(self):
        return self._characters()[3]

    def n_characters(self):
        return self._characters()[4]

    def _which_cycle(self, i_edge):
        c = self.cycles
        for k in range(len(c)):
            c_k = c[k]
            for s in range(len(c_k)):
                if c_k[s] == i_edge:
                    return k
        return None

    def _index_edge_cycle_orientated(self, d, a):
        r"""
        return index of the left edge, the right edge, depending on the orientation
        the departing and arriving edge
        """
        i, j = self.label_position(a)
        i_edge = self._index_edge(i, j, d) + 1
        if i == 1 and j == self.number_of_intervals[1] - 1:                                     #the edge with the biggest number
            i_next = self._index_edge(i, self.number_of_intervals[0], d) + 1
        else:
            i_next = i_edge + 1
        if (i == 0 and self._intervals[i][j].orientation == 1) or (i == 1 and self._intervals[i][j].orientation == -1) :
            return self._which_cycle(i_edge), self._which_cycle(i_next)
        else :
            return self._which_cycle(i_next), self._which_cycle(i_edge)

    def _index_edge_cycle_non_orientated(self, d, a):
        i, j = self.label_position(a)
        i_left = self._index_edge(i, j, d) + 1
        if i == 1 and j == self.number_of_intervals[1] - 1:                                     #the edge with the biggest number
            i_right = self._index_edge(i, self.number_of_intervals[0], d) + 1
        else:
            i_right = i_left + 1
        return self._which_cycle(i_left), self._which_cycle(i_right)

    def _h_zero_sigma(self):
        r"""
        Return matrix of the border application
        """
        from sage.modules.free_module_element import vector
        from sage.rings.integer_ring import ZZ
        projected_vectors = self.canonical_VectPaths().copy(vector(ZZ, self.n_cycles))
        for d, a in self.cover_generators():
            i_depart_cycle, i_arrive_cycle = self._index_edge_cycle_orientated(d, a)
            projected_vectors.val(d,a)[i_arrive_cycle] += 1
            projected_vectors.val(d,a)[i_depart_cycle] -= 1
        from sage.matrix.constructor import matrix
        return self.h_one_rel_to_generator()*matrix(ZZ, [projected_vectors.val(d,a) for d, a in self.cover_generators()])

    @cached_method
    def h_one_projection(self):               #projection in h_one_rel on the kernel of the border application
        ker = self._h_zero_sigma().kernel()
        d = ker.degree()
        r = ker.rank()
        from sage.matrix.constructor import matrix
        def complete(basis):
            from copy import copy
            from sage.modules.free_module_element import vector
            from sage.rings.integer_ring import ZZ
            B = copy(basis)
            k = len(B)
            c = 0
            def id(i):
                res = vector(ZZ, d)
                res[i] = 1
                return res
            while k < d:
                if matrix(B + [id(c)]).rank() == k + 1:
                    B += [id(c)]
                    k += 1
                c += 1
            return matrix(B)
        B = complete(self._h_zero_sigma().kernel().basis())
        P = matrix(d)
        for i in range(r):
            P[i,i] = 1
        return B.inverse()*P*B

    @cached_method
    def h_one_rel_to_h_one(self):
        ker = self._h_zero_sigma().kernel()
        d = ker.degree()
        r = ker.rank()
        from sage.matrix.constructor import matrix
        def complete(basis):
            from copy import copy
            from sage.modules.free_module_element import vector
            from sage.rings.integer_ring import ZZ
            B = copy(basis)
            k = len(B)
            c = 0
            def id(i):
                res = vector(ZZ, d)
                res[i] = 1
                return res
            while k < d:
                if matrix(B + [id(c)]).rank() == k + 1:
                    B += [id(c)]
                    k += 1
                c += 1
            return matrix(B)
        B_c = complete(self._h_zero_sigma().kernel().basis())
        P = matrix(d,r)
        for i in range(r):
            P[i,i] = 1
        return B_c.inverse()*P

    @cached_method
    def h_one_to_h_one_rel(self):               #projection in h_one_rel on the kernel of the border application
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ
        ker = self._h_zero_sigma().kernel()
        return matrix(ZZ, ker.dimension(), ker.degree(), ker.basis())

    @cached_method
    def h_one_to_generator(self):
        return self.h_one_to_h_one_rel()*self.h_one_rel_to_generator()

    @cached_method
    def h_one_rel_to_generator(self):
        return self._free_basis()[0]

    @cached_method
    def generator_to_h_one_rel(self):
        return self._free_basis()[1]

    @cached_method
    def generator_to_h_one(self):
        return self.generator_to_h_one_rel()*self.h_one_rel_to_h_one()

    @cached_method
    def representation_matrix(self, t):
        r"""
        Matrix from H1(X, Sigma) to H1(X, Sigma)
        action of $t \in G$ galois group
        """
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ
        from sage.rings.integer import Integer
        return self.h_one_to_generator()*matrix(ZZ, [[Integer(self.galois_group_permutation()[t](n) == m and a == b) for n in self.degrees()
                                                      for a in self.labels()]
                            for m in self.degrees() for b in self.labels()])*self.generator_to_h_one()

    @cached_method
    def number_of_character_appeareance(self, ind_character):
        s = 0
        for t in range(self.galois_group_order()):
            s += (self.character_table()[ind_character][t])*(self.representation_matrix(t).trace())
        return int(RR(s/self.galois_group_order()))

    def dimension(self, i):
        s = 0
        for t in range(self.galois_group_order()):
            s += self.character_table()[i][t]*self.galois_group_permutation()[t].number_of_fixed_points()*self.degree
        return int(RR(s/self.galois_group_order()))

    def _line_double(self,k):                                     #find a label in the line k which appears twice, return None if there is none
        line_k = self._intervals[k]
        n = len(line_k)
        for i in range(n):
            for j in range(i+1, n) :
                if line_k[i].label == line_k[j].label :
                    return j
        return(None)

    def _double(self):                                            #find a label which appears twice on the same line
        j_0 = self._line_double(0)
        if j_0 <> None:
            return (0, j_0)
        j_1 = self._line_double(1)
        if j_1 <> None:
            return (1, j_1)
        return (None, None)

    def _is_double(self, i, j):
        i_p, j_p = self._twins[i][j]
        return (i_p == i)

    def label(self, i, j):           # label of the i, j interval
        return(self._intervals[i][j].label)

    def label_position(self, lab):
        r"""
        Give position in polygone of the interval which will carry the name corresponding to the label and degree.
        """
        for j in range(self.number_of_intervals[0]):
            if self.label(0, j) == lab and self.give_name(0, j):
                return (0, j)
        for j in range(self.number_of_intervals[1]):
            if self.label(1, j) == lab  and self.give_name(1, j):
                return (1, j)
        return None

    def length(self, lab):
        return(self._lengths(lab))

    def permutation(self,i,j):
        r"""
        Returns the permutation associated to the label of the interval in line `i` and position `j`
        """
        return(self._permutations[self.label(i,j)])

    def min_length(self):
        r"""
        Gives the minimal length of all intervals
        """
        return(min(self._lengths.values()))

    def max_length(self):
        r"""
        Gives the maximal length of all intervals
        """
        return(max(self._lengths.values()))

    def _sum(self,i,label):                                       #return(sum of lengths of intervals on a line, number of times label appears)
        lg = self._lengths
        line = self._intervals[i]
        n = len(line)
        s, compt = 0, 0
        for k in range(n) :
            if line[k].label == label :
                compt += 1
            else :
                s += lg[line[k].label]
        return (s, compt)

    def rand_lg(self):
        r"""
        To get random length associated to your interval exchange. It modify directly your interval exchange
        and returns it. The total length of the two intervals is normalised to `1`.
        """
        lg = {self._intervals[i][j].label : R.random_element(0,1) for i in range(2) for j in range(self.number_of_intervals[i])}
        self._lengths = lg
        lab = lg.keys()
        i_double, j_double = self._double()
        if i_double == None:
            return(self)
        lab_double = self._intervals[i_double][j_double].label
        (s_0, l_0), (s_1, l_1) = self._sum(0,lab_double), self._sum(1,lab_double)
        l = (s_0 - s_1)/(l_1 - l_0)
        if l > 0:
            lg[lab_double] = l
        else:
            j_double = self._line_double((i_double + 1) % 2)
            if j_double == None:
                print self
                raise NameError('No double')
            lab_double = self._intervals[(i_double + 1) % 2][j_double].label
            (s_0, l_0), (s_1, l_1) = self._sum(0,lab_double), self._sum(1,lab_double)
            l = (s_0 - s_1)/(l_1 - l_0)
            if l <= 0:
                print self
                raise NameError('Negative Value in IE')
            lg[lab_double] = l

    def diff_sum(self):
        l_0, l_1 = self._intervals[0], self._intervals[1]
        lg = self._lengths
        s0 = sum([lg[l_0[i].label] for i in range(len(l_0))])
        s1 = sum([lg[l_1[i].label] for i in range(len(l_1))])
        return(abs(s0-s1))

    def normalise(self):
        r"""
        As the length goes to zero with the Rauzy transform, you sometime have to renormalise your total length to 1 in
        order to have non zero length (there is some troubles with computer approximation).
        This normalise the length of both intervals to one.
        """
        l_0, l_1 = self._intervals[0], self._intervals[1]
        lg = self._lengths
        s0 = sum([lg[l_0[i].label] for i in range(len(l_0))])
        s1 = sum([lg[l_1[i].label] for i in range(len(l_1))])
        for a in lg.iterkeys():
            lg[a] = 2*lg[a]/(s0+s1)
        if abs(s0 - s1) > 2**(-50):
            (i,j) = self._double()
            label = self._intervals[i][j].label
            (s_0, l_0) = self._sum(0,label)
            (s_1, l_1) = self._sum(1,label)
            if l_0 <> l_1 :
                val = (s_0 - s_1)/(l_1 - l_0)
                if val > 0 :
                    lg[label] = val
                else :
                    print (self, val)
                    raise NameError('Valeur neg')
        return (s0)

    def rauzy_type(self) :
        r"""
        To apply Rauzy transform, you need to know which one of the two last intervals in the lines is the shortest.
        This returns the couple (lign where the last interval is the shortest, lign where it is the longest)

        EXAMPLES::
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
            ...   [Interval('b',1), Interval('c', 1), Interval('c',-1)]]

            sage: lengths = {'a': 2, 'b': 3, 'c': 2}

            sage: permutations = {lab: Permutations(1).identity() for lab in lengths.iterkeys()}

            sage: IntExchange(intervals, lengths, permutations).rauzy_type()
            (1, 0)

            sage: lengths['b'] = 1
            sage: IntExchange(intervals, lengths, permutations).rauzy_type()
            (0, 1)

        """
        line_0, line_1 = self._intervals[0], self._intervals[1]
        a, b = line_0[len(line_0)-1].label, line_1[len(line_1)-1].label
        if self._lengths[a] <= self._lengths[b] :
            return (0,1)
        else :
            return (1,0)

    def give_name(self, i, j) :
        r"""
        Answers if the interval correspond to the given convention and take the number of the corresponding copy.
        See IntExchange for more information about the convention.
        """
        return (self._intervals[i][j].orientation == 1)

    def move(self, i_init, j_init, i_dest, j_dest):
        r"""
        Move the interval from position `i_{init}, j_{init}` to `i_{dest}, j_{dest}` with new_orientation as orientation

        EXAMPLES::
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
            ...   [Interval('b',1), Interval('c', 1), Interval('c',-1)]]

            sage: lengths = {'a': 2, 'b': 3, 'c': 2}

            sage: permutations = {lab: Permutations(1).identity() for lab in lengths.iterkeys()}

            sage: ie = IntExchange(intervals, lengths, permutations)
            sage: ie
            Intervals    : [ a, -a,  b]
                           [ b,  c, -c]
            Lengths      : {'a': 2, 'c': 2, 'b': 3}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>

            sage: ie.move(0,1,0,0,-1)
            sage: ie
            Intervals    : [-a,  a,  b]
                           [ b,  c, -c]
            Lengths      : {'a': 2, 'c': 2, 'b': 3}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>

        """
        if i_init == i_dest:
            if j_init >= j_dest:
                self._intervals[i_dest].insert(j_dest, self._intervals[i_init][j_init])
                del self._intervals[i_init][j_init + 1]
            else:
                self._intervals[i_dest].insert(j_dest + 1, self._intervals[i_init][j_init])
                del self._intervals[i_init][j_init]
        else:
            self._intervals[i_dest].insert(j_dest, self._intervals[i_init][j_init])
            del self._intervals[i_init][j_init]
            self.number_of_intervals[i_init] -= 1
            self.number_of_intervals[i_dest] += 1
        self._twins = [[self.twin(i, j) for j in range(self.number_of_intervals[i])] for i in range(2)]

    def name(self,i, j) :
        r"""
        Give permutation which associates number of the copy from which interval takes his name.
        """
        from sage.combinat.permutation import Permutations
        if self.give_name(i,j) :
            return Permutations(self.degree).identity()
        else :
            return self.permutation(i,j).inverse()

    def ident_rev(self,i,j):
        r"""
        Return the permutation for the given saddle connection of identifiction according to the level of the cover.
        """
        if self.give_name(i,j):
            return self.permutation(i,j)
        else :
            return self.permutation(i,j).inverse()


    def rauzy(self, debug = 0):

        r"""
        Apply the rauzy transform to the interval exchange.

        EXAMPLES::

            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
            ...   [Interval('b',1), Interval('c', 1), Interval('c',-1)]]

            sage: lengths = {'a': 2, 'b': 3, 'c': 2}

            sage: permutations = {lab: Permutations(1).identity() for lab in lengths.iterkeys()}

            sage: ie = IntExchange(intervals, lengths, permutations)
            sage: ie
            Intervals    : [ a, -a,  b]
                           [ b,  c, -c]
            Lengths      : {'a': 2, 'c': 2, 'b': 3}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>

            sage: ie.rauzy()
            (-1, 'c', 'b', [1], [1], [1], True, [1])

            sage: ie
            Intervals    : [ a, -a,  b]
                           [ b, -c,  c]
            Lengths      : {'a': 2, 'c': 2, 'b': 1}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>
            sage: ie.rauzy()
            (1, 'b', 'c', [1], [1], [1], True, [1])

            sage: ie
            Intervals    : [ a, -a]
                           [ b, -b, -c,  c]
            Lengths      : {'a': 2, 'c': 1, 'b': 1}
            Permutations : {'a': [1], 'c': [1], 'b': [1]}
            <BLANKLINE>


        Non trivial cover case
        ::

            sage: sigma_a = Permutation('(1,2)(3,4)(5,6)')
            sage: sigma_b = Permutation('(1,3)(2,4)(5)(6)')
            sage: sigma_c = Permutation('(1,6)(2,5)(3,4)')

            sage: permutations = {'a': sigma_a, 'b': sigma_b, 'c': sigma_c}
            sage: ie = IntExchange(intervals, lengths, permutations)

            sage: ie.rauzy_rev()
            (-1, 'c', 'b', [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6], [6, 5, 4, 3, 2, 1], True, [3, 4, 1, 2, 5, 6])

            sage: ie
            Intervals    : [ a, -a,  b]
                           [ b, -c,  c]
            Lengths      : {'a': 2, 'c': 2, 'b': 1}
            Permutations : {'a': [2, 1, 4, 3, 6, 5], 'c': [4, 3, 6, 5, 2, 1], 'b': [3, 4, 1, 2, 5, 6]}

        """
        from sage.combinat.permutation import Permutation, Permutations
        (i_0, i_1) = self.rauzy_type()
        d_0, d_1 = self.number_of_intervals[i_0] - 1, self.number_of_intervals[i_1] - 1
        A, B = self._intervals[i_0][d_0], self._intervals[i_1][d_1]        #A is the shortest of both final intervals B the longest
        self._lengths[B.label] = self._lengths[B.label] - self._lengths[A.label]                       #CHANGE LENGTH
        (i_1_twin, d_1_twin) = self._twins[i_1][d_1]                        #CHANGE PERMUTATION #compute image of the long interval to see first return
        (i_0_twin, d_0_twin) = self._twins[i_0][d_0]                        #compute image of the short interval
        if debug :
            print self
            print "win :  ",
            print self.label(i_1, d_1)
            print "los :  ",
            print self.label(i_0, d_0)
        A_twin_lab = self.label(i_0_twin, d_0_twin)
        name_A, name_A_twin, name_B = self.name(i_0, d_0), self.name(i_0_twin,d_0_twin), self.name(i_1,d_1)
        ident_A, ident_A_twin, ident_B = self.ident_rev(i_0,d_0), self.ident_rev(i_0_twin,d_0_twin), self.ident_rev(i_1,d_1)
        if A.orientation == 1 and B.orientation == 1:
            perm_one = Permutations(self.degree).identity()
            perm_two = self._permutations[B.label].inverse()
        if A.orientation == 1 and B.orientation == -1:
            perm_one = self._permutations[B.label]
            perm_two = self._permutations[B.label]
        if A.orientation == -1 and B.orientation == 1:
            perm_one = self._permutations[A.label].inverse()
            perm_two = Permutations(self.degree).identity()
        if A.orientation == -1 and B.orientation == -1:
            perm_one = self._permutations[B.label]*self._permutations[A.label].inverse()
            perm_two = Permutations(self.degree).identity()
        if i_1 <> i_1_twin:
            dj = 1
        else:
            dj = 0                                                          #the way we insert depends on the change of orientation
        self.move(i_0, d_0, i_1_twin, d_1_twin + dj)
        if self.give_name(i_1_twin, d_1_twin+dj) :                          #CHANGE PERMUTATIONS, if the image of A is the one who gives the name
            self._permutations[A.label] = ident_B.inverse()*ident_A         #beware, the sense of composition is inversed !!
        else :
            self._permutations[A.label] = ident_A_twin*ident_B
        return (A.label, B.label, A.orientation * B.orientation, perm_one, perm_two)

    def _type(self, i, j):
        i_twin, _ = self._twins[i][j]
        if i_twin == i:
            return 0
        else:
            return 1

    def edge_permutation(self):
        r"""
        We number edges of the upper part of the polygon, and return the permutation of the edge you
        get while winding around it. It is numbered from 0 to n, the n-1 first edges numbered i correspond to the left edge of the
        i-th interval, the n-th is the rightmost one.
        """
        from sage.combinat.permutation import Permutation
        perm = range(n+1)
        for j in range (n+1):
            if j < self.number_of_intervals[0] :
                i_twin, j_twin = self._twins[0][j]
            else:
                i_twin, j_twin = self._twins[1][len(self._intervals[1])-1]
            if i_twin == 0:
                perm[j] = j_twin + 1
            else:
                while i_twin <> 0 and j_twin <> 0:
                    i_twin, j_twin = self._twins[i_twin][j_twin-1]
                if j_twin == 0 and i_twin == 1:
                    perm[j] = 0
                else:
                    perm[j] = j_twin + 1
        return (Permutation([perm[i] + 1 for i in range(n+1)]))

    def twin_cover(self, i, j, d):
        if j < self.number_of_intervals[i]:
            (i_t, j_t), s = self._twins[i][j], self.ident_rev(i,j)(d)
        if j == self.number_of_intervals[i]:
            i_t, j_t = self._twins[i][j-1]
            j_t += 1
            s = self.ident_rev(i, j-1)(d)
        return(i_t, j_t, s)

    def _index_edge(self, i, j, d):
        """
        Give the number given to the edge at the left of the interval at i, j in level c
        """
        n = self.total_number_of_intervals_base
        if j == 0:
            return((d-1)*self.total_number_of_intervals_base)
        if j == self.number_of_intervals[i]:
            return((d-1)*self.total_number_of_intervals_base + self.number_of_intervals[0])
        return( (d-1)*self.total_number_of_intervals_base + i*self.number_of_intervals[0] + j)

    def _next_edge(self, i, j, d):
        r"""
        return next edge if you go counterclockwise around the point at the left of i,j,c interval
        """
        if i == 0:
            i_t, j_t, d_t = self.twin_cover(0, j, d)
            j_t += 1 - i_t
        if i == 1:
            i_t, j_t, d_t = self.twin_cover(1, j-1, d)
            j_t += 1 - i_t
        return self._index_edge(i_t, j_t, d_t)

    def _next_path(self, i, j, d, side):
        r"""
        return next intersected path if you wind counterclockwise
        """
        if i == 0 and side == "Right":
            if j < self.number_of_intervals[0] - 1:
                return 0, j + 1, d, "Left"
            if j == self.number_of_intervals[0] - 1:
                return 1, self.number_of_intervals[1] - 1, d, "Right"
        if i == 0 and side == "Left":
            i_t, j_t, d_t = self.twin_cover(i, j, d)
            if i_t == 0:
                return i_t, j_t, d_t, "Right"
            if i_t == 1:
                return i_t, j_t, d_t, "Left"
        if i == 1 and side == "Right":
            i_t, j_t, d_t = self.twin_cover(i, j, d)
            if i_t == 0:
                return i_t, j_t, d_t, "Right"
            if i_t == 1:
                return i_t, j_t, d_t, "Left"
        if i == 1 and side == "Left":
            if j == 0:
                return 0, 0, d, "Left"
            if j > 0:
                return 1, j - 1, d, "Right"

    def _next_path_generator(self, d, a, side):
        r"""
        Return next path with format degree and label, st the end point of orientated path is the edge
        """
        i, j = self.label_position(a)
        i_aux, j_aux, d_aux, side_aux = self._next_path(i, j, d, side)
        while not self.give_name(i_aux, j_aux):
            i_aux, j_aux, d_aux, side_aux = self._next_path(i_aux, j_aux, d_aux, side_aux)
        return d_aux, self.label(i_aux, j_aux), side_aux

    def _edge_cover_permutation(self):
        from sage.combinat.permutation import Permutation
        perm = range(self.total_number_of_intervals_cover)
        for c in self.degrees():
            for j in range(self.number_of_intervals[0]):
                perm[self._index_edge(0, j, c)] = self._next_edge(0, j, c) + 1    #index must be greater than 1 for permutations
            for j in range(1, self.number_of_intervals[1] + 1):
                perm[self._index_edge(1, j, c)] = self._next_edge(1, j, c) + 1
        return (Permutation(perm))

    def _is_extremal(self, i):
        l = [len(self._intervals[0]), len(self._intervals[1])]
        n = l[0] + l[1]
        j_base = (i-1) % n
        return( j_base == 0 or j_base == l[0] )

    def singularities_order(self):
        r"""
        ..NOTE::
        To compute the angle around a singularity, we use the following simple argument.
        If you apply the Teichmueller flow, the singularity order won't change, and the angle betwin two intervals
        will tend to ``\pi`` except at extremal point where it tends to 0.
        We first count the angle and then substract it by 2 to get the quadratic singularity order.
 
        EXAMPLES::
        
        Quadratic[{-1, -1, -1, -1}] ::
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)], [Interval('b',-1), Interval('c', 1), Interval('c',-1)]]
            sage: R = orientable_double_cover_IntExchange(intervals)
            sage: R.singularities_order()
            [-1, -1, -1, -1]
        """
        singularities = [0 for k in range(self.n_cycles)]
        for k in range(self.n_cycles):
            c = self.cycles[k]
            for i in range(len(c)):
                if not self._is_extremal(c[i]):
                    singularities[k] += 1
            singularities[k] -= 2
        return(singularities)

    def _occurences(self, l):
        l.sort()
        elements = []
        nb_of_occurence = []
        i = 0
        while i < len(l):
            if i == 0 or l[i] <> l[i-1]:
                elements.append(l[i])
                nb_of_occurence.append(0)
            nb_of_occurence[len(elements) - 1] += 1
            i += 1
        return(elements, nb_of_occurence)

    def stratum(self):
        r"""
        EXAMPLES::
        
        Quadratic[{-1, -1, -1, -1}] ::
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)], [Interval('b',-1), Interval('c', 1), Interval('c',-1)]]
            sage: R = orientable_double_cover_IntExchange(intervals)
            sage: R.stratum()
            'Q(-1^4)'

        Quadratic[{5, -1}] ::
            sage: intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('b',-1), Interval('d',1)], 
            ....:              [Interval('d',1), Interval('e',1), Interval('e',-1), Interval('c',1), Interval('a',1)]]
            sage: R = trivial_IntExchange(intervals)
            sage: R.stratum()
            'Q(-1^1, 5^1)'
        """
        order_of_singularity, nb_of_occurence = self._occurences(self.singularities_order())
        n = len(order_of_singularity)
        if n == 1:
            return "Q(" + str(order_of_singularity[0]) + "^" + str(nb_of_occurence[0]) + ")"
        s = "Q("
        s += str(order_of_singularity[0]) + "^" + str(nb_of_occurence[0]) + ","
        for i in range(1, n - 1):
            s += " " + str(order_of_singularity[i]) + "^" + str(nb_of_occurence[i]) + ","
        s += " " + str(order_of_singularity[n-1]) + "^" + str(nb_of_occurence[n-1]) + ")"
        return(s)

    def genus(self):
        r"""
        EXAMPLES::
        
        Quadratic[{-1, -1, -1, -1}] ::
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)], [Interval('b',-1), Interval('c', 1), Interval('c',-1)]]
            sage: R = orientable_double_cover_IntExchange(intervals)
            sage: R.genus()
            0

        Quadratic[{5, -1}] ::
            sage: intervals = [[Interval('a',1), Interval('b',1), Interval('c',1), Interval('b',-1), Interval('d',1)], 
            ....:              [Interval('d',1), Interval('e',1), Interval('e',-1), Interval('c',1), Interval('a',1)]]
            sage: R = trivial_IntExchange(intervals)
            sage: R.genus()
            2
        """
        sum_d = sum(self.singularities_order())
        return((sum_d + 4)/4)

    def couple_paths(self, v, i_edge_cycle_ref):
        from sage.rings.integer import Integer
        arrive, depart = [], []
        for d, a in self.cover_generators():
            if v.val(d,a) <> 0:
                i_left_cycle, i_right_cycle = self._index_edge_cycle_non_orientated(d, a)
                i_depart_cycle, i_arrive_cycle = self._index_edge_cycle_orientated(d, a)

                if v.val(d,a) < 0:
                    i_depart_cycle, i_arrive_cycle = i_arrive_cycle, i_depart_cycle

                if i_depart_cycle == i_edge_cycle_ref:
                    if i_left_cycle == i_right_cycle:
                        depart += [(d, a, "Left")]*Integer(abs(v.val(d,a)))
                    else:
                        if i_left_cycle == i_edge_cycle_ref:
                            depart += [(d, a, "Left")]*Integer(abs(v.val(d,a)))
                        if i_right_cycle == i_edge_cycle_ref:
                            depart += [(d, a, "Right")]*Integer(abs(v.val(d,a)))

                if i_arrive_cycle == i_edge_cycle_ref:
                    if i_left_cycle == i_right_cycle:
                        arrive += [(d, a, "Right")]*Integer(abs(v.val(d,a)))
                    else:
                        if i_left_cycle == i_edge_cycle_ref:
                            arrive += [(d, a, "Left")]*Integer(abs(v.val(d,a)))
                        if i_right_cycle == i_edge_cycle_ref:
                            arrive += [(d, a, "Right")]*Integer(abs(v.val(d,a)))

        return arrive, depart

    def interval_to_name(self, i, j, d):
        lab = self._intervals[i][j].label
        deg = self.name(i,j)(d)
        return(deg, lab)

    def label_orientation(self, a):
        i, j = self.label_position(a)
        if i == 0 :
            return 1
        else :
            return -1

    def is_cycle(self, d, a):
        left, right = self._index_edge_cycle_non_orientated(d, a)
        return left == right

    def intersection(self, (d_arr, a_arr, side_arr), (d_dep, a_dep, side_dep), v):
        r"""
        return the number of intersection of the intervals arriving to an edge and departing from it, with the intervals
        given by v.

        EXAMPLES::
            sage: intervals = [[Interval('a',1), Interval('a',-1), Interval('b',1)],
            ...   [Interval('b',1), Interval('c', 1), Interval('c',-1)]]
            sage: ie = orientable_double_cover_IntExchange(intervals)

            sage: v = ie.canonical_VectPaths()
            sage: v._vect[1]['a'] = 1
            sage: v._vect[2]['a'] = -1
            sage: ie.intersection((1, 'b', 'Right'), (2, 'b', 'Right'), v)
            0
        """
        d_aux, a_aux, side_aux = self._next_path_generator(d_arr, a_arr, side_arr)
        c = 0
        s = 0
        while (d_aux, a_aux, side_aux) <> (d_dep, a_dep, side_dep):
            if side_aux == "Right":
                c += v.val(d_aux, a_aux)*self.label_orientation(a_aux)
            if side_aux == "Left":
                c -= v.val(d_aux, a_aux)*self.label_orientation(a_aux)
            d_aux, a_aux, side_aux = self._next_path_generator(d_aux, a_aux, side_aux)
            s += 1
            if s > 2^15 :
                print self
                print a_arr, a_dep, a_aux
                print "d_arr :  " + str(d_arr) + "  a_arr :  " + str(a_arr) + "side_arr :  " + str(side_arr)
                print "d_dep :  " + str(d_dep) + "  a_dep :  " + str(a_dep) + "side_dep :  " + str(side_dep)
                print "d_aux :  " + str(d_aux) + "  a_aux :  " + str(a_aux) + "side_aux :  " + str(side_aux) +"\n"
                raise NameError('Problem with intersection : not finite')
        return c

    def global_intersection(self, v, w):
        r"""
        give intersection of VectPath from the canonical generator family
        """
        res = 0
        for i in range(self.n_cycles):
            arrive, depart = self.couple_paths(v, i)
            if len(arrive) <> len(depart):
                print "i : %s"%i
                print "arrive :  " + str(arrive)
                print "depart : " + str(depart) + "\n"
                print "v :  " + str(v),
                print "w :  " + str(w)
                raise NameError('Not an element from kernel of border application')
            for j in range(len(arrive)):
                res += self.intersection(arrive[j], depart[j], w)
        return res

    @cached_method
    def intersection_matrix(self):
        dim = self.h_one_to_generator().nrows()
        def inter(i, j):
            v_i = self.canonical_VectPaths().copy_vector(self.h_one_to_generator()[i])
            v_j = self.canonical_VectPaths().copy_vector(self.h_one_to_generator()[j])
            return self.global_intersection(v_i, v_j)
        from sage.matrix.constructor import matrix
        from sage.rings.universal_cyclotomic_field.universal_cyclotomic_field import E, UniversalCyclotomicField
        UCF = UniversalCyclotomicField()
        return UCF(E(4))*matrix(dim, inter)


    def isotopic_projection(self, i_character, (d_j, a_j)):
        r"""
        Return a new VectPaths corresponding to the projection of self to the isotropic space of the \chi_i character
        recall the formula p_i_character (d, a) = \sum_{t \in G} char_i(t).conjugate (galois_group_permutation(t)(d), a)
        """
        res = self.canonical_VectPaths()
        char = self.character_table()[i_character]
        coeff = self.character_degree()[i_character]/self.galois_group_order()
        somme = 0
        for t in range(self.galois_group_order()):
            res._vect[self.galois_group_permutation()[t](d_j)][a_j] += coeff*char[t]
        return(res)

    @cached_method
    def isotopic_projection_matrix_before(self, i):
        from sage.matrix.constructor import matrix
        return matrix([self.isotopic_projection(i, (d, a)).to_list() for d, a in self.cover_generators()])

    def vector_isotopic_projection(self, i_character, v):
        r"""
        Return a new VectPaths corresponding to the projection of self to the isotopic space of the \chi_i character
        """
        from sage.rings.all import CC
        res = self.canonical_VectPaths()
        coeff = CC(self.character_degree()[i_character]/self.galois_group_order())
        for d, a in self.cover_generators():
            somme = 0
            for t in range(self.galois_group_order()):
                somme += CC(self.character_table()[i_character][t])*v.val(self.galois_group_permutation()[t].inverse()(d), a)
            res._vect[d][a] = coeff*somme
        return(res)

    def isotopic_projection_matrix(self, i_character):
        r"""
        Return the projection matrix to the isotopic space of the \chi_i character
        recall the formula :
        p_i_character (d, a) = \sum_{t \in G} char_i(t).conjugate (galois_group_permutation(t)(d), a)
        """
        res = [[0 for _ in self.cover_generators()] for _ in self.cover_generators()]
        char = self.character_table()[i_character]
        coeff = self.character_degree()[i_character]/self.galois_group_order()
        to_number = {}
        nb = 0
        for a in self.labels() :
            to_number[a] = nb
            nb += 1
        for d, a in self.cover_generators():
            for t in range(self.galois_group_order()):
                res[(d - 1) * self.nb_labels + to_number[a]][(self.galois_group_permutation()[t](d) - 1)* self.nb_labels + to_number[a]] += coeff*char[t]
        from sage.matrix.constructor import matrix
        return matrix(res)

    @cached_method
    def intersection_form_isotopic(self, i):
        from sage.functions.other import conjugate
        M = self.h_one_to_generator()*self.isotopic_projection_matrix(i)*self.generator_to_h_one()
        return M*self.intersection_matrix()*M.transpose().conjugate()

    def signature(self, M):
        from sage.rings.all import CC
        P = UCF_mat_to_cyclotomic_field(M).characteristic_polynomial()
        roots = P.change_ring(CC).roots(CC)
        p, q = 0, 0
        for i in range(len(roots)) :
            if abs(roots[i][0].imag_part()) > 10**(-10) :
                print "Problem"
                print roots[i]
                raise NameError('Valeur propre')
            from sage.functions.other import real_part
            if roots[i][0].real_part() < -10**(-10) :
                p += roots[i][1]
            if roots[i][0].real_part() > 10**(-10) :
                q += roots[i][1]
        return(p,q)

    @cached_method
    def signatures_isotopic(self, i):
        r"""
        Returns the signature of the intersection form on the isotropic subspace for the ``i``-th character.

        EXAMPLES::

            sage: R = cyclic_cover_iet(6, [1, 1, 1, 3])
            sage: [R.signatures_isotopic(i) for i in xrange(6)]
            [(0, 0), (1, 1), (0, 2), (2, 0), (0, 1), (1, 0)]

        """
        return self.signature(self.intersection_form_isotopic(i))

    @cached_method
    def dimensions(self):
        r"""
        Return a list ``res`` of dimensions such that ``res`` correspond to the ``i``-th character of the corresponding isotropic subspace in homology.
        """
        from math import ceil
        return [ceil(self.character_degree()[k]*self.number_of_character_appeareance(k)) for k in xrange(self.n_characters())]

    def lyapunov_exponents_H_plus(self, nb_vectors=None, nb_experiments=10,
                                  nb_iterations=32768, verbose=False, output_file=None):
        r"""
        Compute the H^+ Lyapunov exponents for this covering locus.

        It calls the C-library lyap_exp interfaced with Cython. The computation
        might be significantly faster if ``nb_vectors=1`` (or if it is not
        provided but genus is 1).

        INPUT:

        - ``nb_vectors`` -- the number of exponents to compute. The number of
          vectors must not exceed the dimension of the space!

         - ``nb_experiments`` -- the number of experiments to perform. It might
           be around 100 (default value) in order that the estimation of
           confidence interval is accurate enough.

         - ``nb_iterations`` -- the number of iteration of the Rauzy-Zorich
           algorithm to perform for each experiments. The default is 2^15=32768
           which is rather small but provide a good compromise between speed and
           quality of approximation.

        - ``verbose`` -- if ``True`` provide additional informations rather than
          returning only the Lyapunov exponents (i.e. ellapsed time, confidence
          intervals, ...)

        - ``output_file`` -- if provided (as a file object or a string) output
          the additional information in the given file rather than on the
          standard output.

        EXAMPLES::
            sage: R = cyclic_cover_iet(4, [1, 1, 1, 1])
            sage: R.lyapunov_exponents_H_plus()
            [0.9996553085103, 0.0007776980910571506, 0.00022201024035355403]

        """
        import time
        import sage.dynamics.flat_surfaces.lekz as lekz   # the cython bindings
        if nb_vectors is None:
            nb_vectors = self.genus()

        if output_file is None:
            from sys import stdout
            output_file = stdout
        elif isinstance(output_file, str):
            output_file = open(output_file, "w")

        nb_vectors = int(nb_vectors)
        nb_experiments = int(nb_experiments)
        nb_iterations = int(nb_iterations)

        if verbose:
            output_file.write("Stratum : " + str(self.stratum()))
            output_file.write("\n")

        if nb_vectors <= 0:
            raise ValueError("the number of vectors must be positive")
        if nb_experiments <= 0:
            raise ValueError("the number of experiments must be positive")
        if nb_iterations <= 0:
            raise ValueError("the number of iterations must be positive")

        #Translate our structure to the C structure"
        seen = []
        associate = {}
        count = 0
        for i in range(2):
            for j in range(len(self._intervals[i])):
                lab = self._intervals[i][j].label
                if seen.count(lab) == 0:
                    seen.append(lab)
                    associate[lab] = int(count)
                    count += 1
        k = len(self._intervals[0])
        gp = range(self.total_number_of_intervals_base)
        twin = range(self.total_number_of_intervals_base)
        for i in range(2):
            for j in range(self.number_of_intervals[0]):
                gp[j + i*k] = associate[self._intervals[i][j].label]
        n = self.total_number_of_intervals_base // 2
        for i in range(2):
            for j in range(self.number_of_intervals[i]):
                i_twin, j_twin = self._twins[i][j]
                twin[i * k + j] = int(i_twin * k + j_twin)

        sigma = range(n)
        for lab in self.iter_lab():
            sigma[associate[lab]] = map(int,list(self._permutations[lab]))
        sigma = reduce(lambda x, y: x+y, sigma)
        sigma = map(lambda x : x-1, sigma)

        if self._lengths <> None :
            lengths = range(n)
            for lab in self.iter_lab():
                lengths[associate[lab]] = self._lengths[lab]
        else:
            lengths = None

        t0 = time.time()
        res = lekz.lyapunov_exponents_H_plus_cyclic_cover(
                   gp, twin, int(k), int(n), sigma, self.degree,
                   nb_vectors, nb_experiments, nb_iterations, lengths)
        t1 = time.time()

        res_final = []

        m,d = mean_and_std_dev(res[0])
        if verbose:
            from math import log, floor, sqrt
            output_file.write("sample of %d experiments\n"%nb_experiments)
            output_file.write("%d iterations (~2^%d)\n"%(
                    nb_iterations,
                    floor(log(nb_iterations) / log(2))))
            output_file.write("ellapsed time %s\n"%time.strftime("%H:%M:%S",time.gmtime(t1-t0)))
            output_file.write("Lexp Rauzy-Zorich: %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    m,d, 2.576*d/sqrt(nb_experiments)))
        for i in xrange(1,nb_vectors+1):
            m,d = mean_and_std_dev(res[i])
            if verbose:
                output_file.write("theta%d           : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    i,m,d, 2.576*d/sqrt(nb_experiments)))
            res_final.append(m)

        return res_final


    def lyapunov_exponents_H_plus_isotopic(self, nb_vectors=None, nb_experiments=100,
                                                       nb_iterations=32768, verbose=False, output_file=None):
        r"""
        Compute the Lyapunov exponents for this covering locus.
        And make the isotopic decomposition for the Hodge bundle.

        We use an isotopic decomposition since the galois group of the cover acts on the homology group.
        The finite representation of group theory gives a decomposition by characters of the Hodge bundle,
        invariant by Teichm\"uller flow.

        It calls the C-library lyap_exp interfaced with Cython. The computation
        might be significantly faster if ``nb_vectors=1`` (or if it is not
        provided but genus is 1).

        INPUT:

        - ``nb_vectors`` -- the number of exponents to compute. The number of
          vectors must not exceed the dimension of the space!

         - ``nb_experiments`` -- the number of experiments to perform. It might
           be around 100 (default value) in order that the estimation of
           confidence interval is accurate enough.

         - ``nb_iterations`` -- the number of iteration of the Rauzy-Zorich
           algorithm to perform for each experiments. The default is 2^15=32768
           which is rather small but provide a good compromise between speed and
           quality of approximation.

        - ``verbose`` -- if ``True`` provide additional informations rather than
          returning only the Lyapunov exponents (i.e. ellapsed time, confidence
          intervals, ...)

        - ``output_file`` -- if provided (as a file object or a string) output
          the additional information in the given file rather than on the
          standard output.

        OUTPUT:

        Return a list ``res`` of lists of lyapunov exponents for a given character.
        To have access to the corresponding table of character corresponding to
        ``res[i]`` do ``self.character_table[i]``.
        """
        import time
        import sage.dynamics.flat_surfaces.lekz as lekz    # the cython bindings

        nb_vectors = sum([d//2 for d in self.dimensions()])
        if self.genus() <> nb_vectors:
            print "Warning: Sum of dimensions of isotopic subspaces are not equal to the genus"

        if output_file is None:
            from sys import stdout
            output_file = stdout
        elif isinstance(output_file, str):
            output_file = open(output_file, "w")

        nb_vectors = int(nb_vectors)
        nb_experiments = int(nb_experiments)
        nb_iterations = int(nb_iterations)

        if verbose:
            output_file.write("Stratum : " + str(self.stratum()))
            output_file.write("\n\n")

        if nb_vectors <= 0:
            raise ValueError("the number of vectors must be positive")
        if nb_experiments <= 0:
            raise ValueError("the number of experiments must be positive")
        if nb_iterations <= 0:
            raise ValueError("the number of iterations must be positive")

        #Translate our structure to the C structure"
        seen = []
        associate = {}
        count = 0
        for i in range(2):
            for j in range(len(self._intervals[i])):
                lab = self._intervals[i][j].label
                if seen.count(lab) == 0:
                    seen.append(lab)
                    associate[lab] = int(count)
                    count += 1
        k = len(self._intervals[0])
        gp = range(self.total_number_of_intervals_base)
        twin = range(self.total_number_of_intervals_base)
        for i in xrange(2):
            for j in xrange(self.number_of_intervals[0]):
                gp[j + i*k] = associate[self._intervals[i][j].label]
        n = self.total_number_of_intervals_base // 2
        for i in xrange(2):
            for j in xrange(self.number_of_intervals[i]):
                i_twin, j_twin = self._twins[i][j]
                twin[i * k + j] = int(i_twin * k + j_twin)

        size_of_matrix = self.nb_labels * self.degree
        projections = range(size_of_matrix**2 * self.n_characters())
        from sage.functions.other import real_part
        from sage.rings.all import CC
        for i_char in xrange(self.n_characters()):
            for i in xrange(size_of_matrix):
                for j in xrange(size_of_matrix):
                    projections [i_char * (size_of_matrix**2) + i * size_of_matrix + j] = float(real_part(CC(self.isotopic_projection_matrix(i_char)[j][i])))

        sigma = range(n)
        for lab in self.iter_lab():
            sigma[associate[lab]] = map(int,list(self._permutations[lab]))
        sigma = reduce(lambda x, y: x+y, sigma)
        sigma = map(lambda x : x-1, sigma)

        if self._lengths <> None :
            lengths = range(n)
            for lab in self.iter_lab():
                lengths[associate[lab]] = self._lengths[lab]
        else:
            lengths = None

        #Launch the computation
        t0 = time.time()
        from sage.misc.prandom import random
        from sage.rings.integer import Integer
        v = VectPaths(self.labels(), self.degrees()).random()
        res = lekz.lyapunov_exponents_H_plus_cyclic_cover(
                   gp, twin, int(k), int(n), sigma, self.degree,
                   nb_vectors, nb_experiments, nb_iterations, lengths,
                   self.n_characters(),
                   [Integer(d//2) for d in self.dimensions()],
                   projections)
        t1 = time.time()

        #Write the result
        res_final = []
        m,d = mean_and_std_dev(res[0])
        if verbose:
            from math import log, floor, sqrt
            output_file.write("sample of %d experiments\n"%nb_experiments)
            output_file.write("%d iterations (~2^%d)\n"%(
                    nb_iterations,
                    floor(log(nb_iterations) / log(2))))
            output_file.write("ellapsed time %s\n"%time.strftime("%H:%M:%S",time.gmtime(t1-t0)))
            output_file.write("Lexp Rauzy-Zorich: %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    m,d, 2.576*d/sqrt(nb_experiments)))

        i_0 = 1
        from sage.rings.integer import Integer
        for i_char in xrange(self.n_characters()):
            res_rel = []
            if verbose:
                output_file.write("##### char_%d #####\n"%(i_char))
            for i in xrange(i_0, i_0 + Integer(self.dimensions()[i_char]//2)) :
                m,d = mean_and_std_dev(res[i])
                if verbose:
                    output_file.write("theta%d           : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                            i,m,d, 2.576*d/sqrt(nb_experiments)))
                res_rel.append(m)
            i_0 += Integer(self.dimensions()[i_char]//2)
            res_final.append(res_rel)
        return res_final


def mean_and_std_dev(l):
    r"""
    Return the mean and standard deviation of the floatting point numbers in
    the list l.

    The implementation is very naive and should not be used for large list
    (>1000) of numbers.

    .. NOTE::
    mean and std are implemented in Sage but are quite buggy!
    """
    from math import sqrt
    m = sum(l) / len(l)
    if len(l) == 1:
        d = 0
    else:
        d = sum((x-m)**2 for x in l) / (len(l)-1)
    return m,sqrt(d)
