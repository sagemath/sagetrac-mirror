#*****************************************************************************
#       Copyright (C) 2016 Mark Shimozono <mshimo at vt.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
#
# Compression of the Ram Yip formula for symmetric Koornwinder polynomials
#
############################################

from sage.misc.cachefunc import cached_method
from sage.sets.set import Set, Set_object
from sage.sets.family import Family
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ

from sage.algebras.double_affine_hecke_algebra import DoubleAffineHeckeAlgebraSansDuality, RamYipFormula

class CompressionDomain(RamYipFormula):
    r"""
    A container class for computations involving the Ram-Yip formula
    for symmetric Koornwinder polynomials.

    INPUTS::

    - n -- The rank of the underlying finite type B (or BC) root system.
    - ``highest_weight_coefficients`` -- a list of `n` nonnegative integers `[c_1, c_2,\dotsc, c_n]`.
    The integer `c_i` specifies the number of columns of height `i` in a partition of height at most `n`.

    The corresponding type `B_n` dominant weight is::

        \mu = \sum_{i=1}^{n-1} c_i \omega_i + 2 c_n \omega_n

    where $\omega_i$ is the `i`-th fundamental weight. Note that $\omega_n$ is a spin weight (a half column of
    height $n$) so $2\omega_n$ corresponds to a full-width column column of height $n$.
    The resulting weight lies in the root lattice (which is the correct lattice for the Koornwinder double affine
    type) and its translation element lies in the nonextended affine Weyl group.

    EXAMPLES::

        sage: from sage.algebras.compression_domain import CompressionDomain
        sage: d = CompressionDomain(5, [1,2,0,1,0], no_paths=True)

    Up to normalization, we will compute the Ram-Yip formula for the symmetric Koornwinder polynomial using
    the Hecke symmetrization of the nonsymmetric Koornwinder polynomial `E_\mu` over the finite Weyl group::

        \sum_{u \in W} v_w T_w E_{\mu}

    Here `v_w` is the eigenvalue of `T_w` on the generator of the polynomial module.
    At the cost of a different normalizing constant, one may instead sum over coset representatives
    involving the stabilizer of `\mu` in `W`. We will also consider a variant (via the ``antidominant`` option)
    which uses `E_{w_0(\mu)}` instead of `E_\mu` where `w_0` is the long element of `W`.
    For the Koornwinder DAHA we have `w_0(\mu) = - \mu`.

    The Ram-Yip formula uses a reduced word ``rw`` for the minimum coset representative `m_\mu`
    (resp. `m_{w_0(\mu)}`) for the translation element `t_\mu` (resp. `t_{w_0(\mu)}`) with respect
    to the finite Weyl group $W$ acting on the right. 

    We specify the following reduced word for `m_\mu`.
    In general for a weight `\lambda` let `u_\lambda` be the shortest finite Weyl group element
    sending `\lambda` to antidominant. Let `u_i` be this element for the `i`-th fundamental weight.


    of height `i`, with respect to the subgroup that stabilizes the `j`-th fundamental weight.
    Explicitly we have for `i \le j`

    MATH::

        u_{i,j} &= s_i s_{i+1} \dotsm s_j 
        a_{i,j} &= s_i s_{i+1} \dotsm s_{n-1} s_n s_{n-1} \dotsm s_j
        m_{i,j} &= a_{i, i} a_{i-1,i} \dotsm a_{i-j+1,i} u_{i-j, i-1}, u_{i-j-1,i-2} \dotsm u_{1,j} \\
        z_i &= s_0 (s_1s_0)(s_2s_1s_0) \dotsm (s_{i-1} \dotsm s_1s_0).

    The product of `b`s is the minimum coset representative for the column translation
    for the finite Weyl group. The product of the `a`s is the minimum coset representative of
    the element `u_i` (the shortest element taking the single column weight to antidominant)
    with respect to the stabilizer of the `j`-th column weight.
    If this column was the last, then only the `b`s are used.

    The column is given as the two factors, the product of the `b`s and the product of the `a`s.

    An alcove path is a sequence of finite Weyl group elements. To get a filling we
    apply certain of these elements to certain weights.    


    If ``antidominant`` is True then we use the translation element $t_{w_0(\mu)}$.
    Note that $m_{w_0(\mu)} = t_{w_0(\mu)}. Again we factorize according to the list of columns.
    Each column of height `i` uses the reduced word (with `2i` factors)::

        a_{i,i} \dotsm a_{2,i} a_{1,i} b_1 b_2 \dotsm b_i.

    Note that the `b` part is the same as in the dominant case and the `a` part is similar, except that for our
    antidominant situation we pretend that the next column also has the same height `i`. In particular in the antidominant
    case, for each column of height `i` the reduced words are the same (independent of the height of the next column).

    ``sampling`` works similarly as before. If it is 'standard' then for a column of size `i` we get `2i` factors.
    If it is 'fine' then each `a` word gets split in two as before. If it is 'coarse' then as before there
    are two subwords for each column, here the `a` part followed by the `b` part.
    """

    def __init__(self, n, highest_weight_coefficients, **keywords):
        r"""
        The following optional keyword arguments are valid and have False as default value.

        INPUTS::

        - ``no_paths`` -- Do not generate the alcove paths
        - ``verbose`` -- Print additional information
        - ``n_numerators`` -- When printing a rational function, use abbreviations for numerators for `s_n` reflections
        - ``antidominant`` -- Use the translation by the antidominant element in the Weyl orbit of the highest weight

        There is another optional keyword::

        - ``sampling`` -- 'coarse', 'standard', 'fine' (default: 'coarse')

        CHANGE OF PLAN: only 'coarse' is being maintained. Others are not maintained.

        This option determines the map from alcove paths to fillings, as follows.

        An alcove path (for the purposes of the filling) is a sequence of finite Weyl group elements
        whose positions are labeled `0` through `N` where `N` is the length of a given reduced word for the
        minimum coset representative of the translation element. The reflections live in the gaps between these group elements.

        This reduced word is defined by giving a subword, one for each column.
        The longer columns will go before (to the left of) shorter columns.

        Let `t_i` be the translation element by the dominant weight corresponding to a column of height `i`.
        This weight is the fundamental weight except when `i=n` in which case it is twice the `n`-th
        fundamental weight.

        Let `m_i` be the minimum coset representative of `t_i` with respect to the finite Weyl group `W`.
        We have::

            m_i = (s_0) (s_1s_0) (s_2s_1s_0) \dotsm (s_{i-1}\dotsm s_1s_0).

        For a general weight `\lambda` let `u_\lambda` be the element of `W` that sends `\lambda` to
        the antidominant chamber. Let `u_i` be this element for the `i`-th fundamental weight. We have::

            t_i = m_i u_i.

        In the Koornwinder case we have::

            u_i^{-1} &= u_i \\
            m_i^{-1} &= m_i \\
            t_i^{-1} &= u_i m_i

        Finally for `n \ge i \ge j \ge 1` let `m_{i,j}` be the minimum length coset representative of
        `m_i` with respect to the stabilizer of the `j`-th fundamental weight. Letting `a_{j,i}` be the
        product of reflections of ascending index,::

            a_{j,i} = s_j s_{j+1} \dotsm s_i

        `b_{i,j}` the product of reflections from `s_i` up to `s_n` and back down to `s_j`,::

            b_{i,j} = s_i s_{i+1} \dotsm s_{n-1} s_n s_{n-1} \dotsm s_{j+1} s_j

        we have::

            m_{i,j} = b_{i,i} b_{i-1,i} \dotsm b_{i-j+1,i} a_{i-j,i-1} a_{i-j-1,i-2} \dotsm a_{1,j}.

        Finally, the required reduced word for `m_\mu` is obtained by the concatenation of
        reduced words of elements `c_1 c_2\dotsm` where `c_k` represents the `k`-th column of the partition of
        `\mu`. If this column has height `i` and the column to the right has height `j` with `i \ge j` then::

            c_k =  m_i m_{i,j}.

        For the last column (say of height `i`) we just use `m_i`.

        The reduced word for `m_\mu` is broken into factors, 2 for each column, namely, the
        `m_i` part and the `m_{i,j}` part, except for the last column, which only has the `m_i` part.

        Each alcove path defines a sequence of finite Weyl group elements, one more than the length of
        `m_\mu`. In each gap between reflections in the reduced word for `m_\mu` there is a corresponding
        Weyl group element. We pick the subsequence of Weyl elements corresponding to the gaps specified
        by our factorization of `m_\mu`. Before and after each `z_i` we take the corresponding Weyl elements
        and apply them to the weight for the column of height `i` and remember these weights.
        This sequence of weights (that come in pairs) is called a filling.
        It is a theorem that the weight of the alcove path equals half the sum of the weights in the filling
        associated to the alcove path.

        The 'antidominant' option changes around the conventions. Instead of using a reduced word for
        `m_\mu` we use a reduced word for `t_{w_0(\mu)}`, which equals `t_{-\mu}` due to the fact that
        we are in the Koornwinder DAHA. We have `t_i^{-1} = u_i m_i` and::

            u_i &= m_{i,i} \\
                &= b_{i,i} b_{i-1,i} \dotsm b_{1,i}.

        In this case we obtain the filling by applying to the weight of the column of height `i`,
        to the Weyl elements corresponding to the gaps before and after the `m_i` factors.
        """

        self._verbose = keywords.get('verbose',False)
        self._sampling = keywords.get('sampling',None)
        if self._sampling is None:
            self._sampling = 'coarse'
        else:
            assert self._sampling in ('coarse', 'standard', 'fine')

        self._n_numerators = keywords.get('n_numerators',False)

        if n != len(highest_weight_coefficients):
            raise ValueError("The length of the list of coefficients does not equal the rank")
        self._n = n # the rank

        if not all([c in ZZ and c >= 0 for c in highest_weight_coefficients]):
            raise ValueError("Highest weight coefficients must be nonnegative integers")

        if sum(highest_weight_coefficients) == 0:
            raise ValueError("Highest weight must be nonzero")

        # the Koornwinder double affine Hecke algebra.
        HH = DoubleAffineHeckeAlgebraSansDuality(['B',n],untwisted=False,reduced=False,dual_reduced=False, extra_parameters=True)

        self._XL = HH.lattice() # type X ambient lattice; the dual translation lattice for the
        # extended affine Weyl group of type Y.
        # use this to make signed permutations from finite Weyl elements
        self._I0 = range(1,n+1) # the classical index set
        self._standard_basis = Family(self._I0, lambda i: self._XL.basis()[i-1])

        # remember the column multiplicities in a Family
        self._column_multiplicities = Family(self.index_set(), lambda i: highest_weight_coefficients[i-1])
        self._column_weight = Family(self.index_set(), lambda i: (1 if i!=self._n else 2)*self._XL.fundamental_weight(i))
        # the sequence of columns
        self._column = self._listlist_to_list([[i]*self._column_multiplicities[i] for i in self.index_set()])
        self._column = [x for x in reversed(self._column)]
        self._next_column = self._column[1:] + [None]
        self._last_column = self._column[-1]

        # the highest weight
        weight = self._XL.sum(self._column_multiplicities[i]*self._column_weight[i] for i in self.index_set())
        self._the_dominant_weight = weight

        self._antidominant = keywords.get('antidominant',False) # are we in the dominant or antidominant case?

        if self._antidominant:
            self._my_column_word_factored = [self._column_word_factored(c,c) for c in self._column]
            self._the_reduced_word_factored = self._listlist_to_list(self._my_column_word_factored)
            self._sampling_position = self._partial_sums([len(u) for u in self._the_reduced_word_factored])
        else:
            self._my_column_word_factored = [self._column_word_factored(self._column[i],self._next_column[i]) for i in range(len(self._column))]
            self._the_reduced_word_factored = self._listlist_to_list(self._my_column_word_factored)
            self._sampling_position = [0] + self._partial_sums([len(u) for u in self._the_reduced_word_factored])
        self._column_word = [self._listlist_to_list(lis) for lis in self._my_column_word_factored]
        reduced_word = self._listlist_to_list(self._the_reduced_word_factored)

        # all sampling pairs
        self._pair = []
        base = 0
        for i in range(len(self._column)):
            k = self._column[i]
            self._pair  = self._pair + [(base+offset, self._column_weight[k]) for offset in range(len(self._column_word[i]))]
            base = base + len(self._column_word[i])

        # the sampling pairs for creating fillings
        assert len(self._sampling_position) == 2 * len(self._column)

        self._sampling_pair = []
        for i in range(len(self._column)):
            k = self._column[i]
            self._sampling_pair = self._sampling_pair + [[self._sampling_position[2*i], self._column_weight[k]],[self._sampling_position[2*i+1], self._column_weight[k]]]

        if self._verbose:
            print "factorized reduced word"
            for i in range(len(self._column)):
                print "column {}: {}".format(i+1, self._my_column_word_factored[i])
            print "all sampling pairs: "
            print self._pair
            print "sampling pairs: "
            print self._sampling_pair

        if self._antidominant:
            # move the weight to antidominant
            weight = weight.to_dominant_chamber(positive=False)

        RamYipFormula.__init__(self, HH, weight, symmetric=True, reduced_word = reduced_word, verbose=self._verbose)

        # denominators associated with the inversions of m_\lambda; Koornwinder double affine datum
        self._reflection_denom = [(self.DAHA().base_ring().one() - self._eigen[i]**(2 if self._the_reduced_word[i] in [0,n] else 1)).numerator() for i in range(self._N)]

        if self._verbose:
            for i in range(self._N):
                print "inversion {}".format(i+1)
                print "{} \delta - {}".format(self._inversion_depths[i], self._inversions_classical[i])
                print "eigenvalue: {}".format(self._eigen[i])
                print "denom: {}\n".format(self._reflection_denom[i])

      # special Ram-Yip numerator factors for Koornwinders
        ddat = HH.dual_double_affine_type()
        self._base_ring = ddat.base_ring()
        cz = ddat.c_parameter('zero_doubled')
        c2 = ddat.c_parameter('zero') # this is due to using Y parameters and DAHA duality
        v = ddat.parameter('short')
        v0 = ddat.parameter('doubled') # this is due to using Y parameters and DAHA duality

        self._pos_fold_numerator = dict()
        self._neg_fold_numerator = dict()

        for i in range(self._N):
            if self._the_reduced_word[i] == 0:
                # numerator is called to set the parent to the polynomial ring
                # instead of the fraction field, for divisibility tests
                self._pos_fold_numerator[i] = (self._eigen[i]*cz+c2).numerator()   # \chi_i^z
                self._neg_fold_numerator[i] = (self._eigen[i]*c2+cz).numerator()   # \chi_i^2
            # should we disable special factors other than those involving c2 and cz???
            elif self._the_reduced_word[i] == self._n and self._n_numerators:
                self._pos_fold_numerator[i] = (v*self._eigen[i]*(1/v0 - v0)+(1-v**2)).numerator()  # v\chi_i^0
                self._neg_fold_numerator[i] = (v0*self._eigen[i]*(1/v - v)+(1-v0**2)).numerator() # v_0\chi_i^s


        no_paths = keywords.get('no_paths', False)
        if no_paths:
            return
        if self._verbose:
            print "generating database of alcove paths"
        self.generate_database()
        if self._verbose:
            print "done"
    
    def _partial_sums(self, lis):
        r"""
        The sequence of partial sums of a nonempty list of integers.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: CompressionDomain(2,[1,0],no_paths=True)._partial_sums([1,3,2])
            [1, 4, 6]
        """
        if len(lis) == 1:
            return lis
        ps = self._partial_sums(lis[:-1])
        return ps + [ps[-1]+lis[-1]]

    def generate_database(self):
        r"""
        Generate the alcove paths and store them by weight and filling.
        """
        # Create a dictionary with dominant weights as keys
        # values are dictionaries with fillings as keys and values given by lists of alcove paths
        terms = self.alcove_paths()
        self._dt = dict()
        for term in terms:
            filling = self.alcove_path_to_filling(term.alcove_path_directions(), self._sampling_pair, self._u_weight_reduced_word)
            weight = self.weight_of_filling(filling)
            if weight.is_dominant_weight():
                if weight in self._dt.keys():
                    if filling in self._dt[weight].keys():
                        self._dt[weight][filling] = self._dt[weight][filling] + [term]
                    else:
                        self._dt[weight][filling] = [term]
                else:
                    self._dt[weight] = dict()
                    self._dt[weight][filling] = [term]

        assert self._the_dominant_weight in self._dt.keys()
        fill_dict = self._dt[self._the_dominant_weight]
        assert len(fill_dict.keys()) == 1
        key = fill_dict.keys()[0]
        self._the_leading_coef = self.add_terms(key)

    def _listlist_to_list(self, listlist):
        return [inner for outer in listlist for inner in outer]

    def rank(self):
        r"""
        The rank of the Koornwinder DAHA.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: CompressionDomain(3, (1,0,0), no_paths=True).rank()
            3
        """
        return self._n

    @cached_method
    def index_set(self):
        r"""
        The Dynkin node set for the classical Dynkin diagram of type `B_n`.
        This is because we are in the rank `n` Koornwinder DAHA.
        """
        return [i for i in range(1,self.rank()+1)]

    def DAHA(self):
        return self._daha

    def _b_word(self, i, j):
        r"""
        A reduced word that goes from `i` to `n` and down to `j`.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: bword = CompressionDomain(5, (1,0,0,0,0), no_paths=True)._b_word
            sage: bword(2,4)
            [2, 3, 4, 5, 4]
            sage: bword(3,5)
            [3, 4, 5]
            sage: bword(1,1)
            [1, 2, 3, 4, 5, 4, 3, 2, 1]
        """
        return [x for x in range(i, self.rank())] + [x for x in range(self.rank(),j-1,-1)]

    def _a_word(self, i, j):
        r"""
        A reduced word going from `i` up to `j`.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: CompressionDomain(5,(1,0,0,0,0), no_paths=True)._a_word(2,4)
            [2, 3, 4]
        """
        return [x for x in range(i,j+1)]

    def _min_rep_column(self, k):
        r"""
        The affine Grassmannian part of the canonical reduced word for the translation
        by the weight given by a column of height `k`.

        It is denoted `m_i` above.
        Since the Koornwinder DAHA uses type `B_n`, the column of height `k` is the
        `k`-th fundamental weight for `1 \le k \le n-1` and twice the 
        `n`-th fundamental weight if `k=n`.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(3, (1,0,0), no_paths=True)
            sage: [d._min_rep_column(i) for i in d.index_set()]
            [[0], [0, 1, 0], [0, 1, 0, 2, 1, 0]]
        """
        return self._listlist_to_list([[x for x in reversed(range(i))] for i in range(1,k+1)])

    def _min_min(self, i, j):
        r"""
        The element `m_{i,j}`.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(5, (1,0,0,0,0), no_paths=True)
            sage: d._min_min(4, 2)
            [4, 5, 4, 3, 4, 5, 4, 2, 3, 1, 2]
            sage: d._min_min(4, 4)
            [4, 5, 4, 3, 4, 5, 4, 2, 3, 4, 5, 4, 1, 2, 3, 4, 5, 4]
            sage: d._min_min(4, None)
            []
        """
        if j is None:
            return []
        return self._listlist_to_list([self._b_word(x, i) for x in range(i,i-j,-1)] + [self._a_word(k,k+j-1) for k in range(i-j,0,-1)])

    def _column_word_factored(self, i, j):
        if j is None:
            return [self._min_rep_column(i)]
        if self._antidominant:
            return [self._min_min(i, i),self._min_rep_column(i)]
        return [self._min_rep_column(i), self._min_min(i, j)]

    def _repr_(self, short=False):
        r"""
        A string representing ``self``.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: CompressionDomain(2,(2,0),no_paths=True)
            Compression domain with columns [1, 1] dominant no extra numerators for Ram Yip Formula for symmetric Macdonald polynomial of weight (2, 0) and Double Affine Type ['C', 2, 1]^* nonreduced nondual-reduced with extended base ring
        """
        if short:
            return "Compression domain of rank {} with columns {} {}dominant {} extra numerators".format(self._n, self._column, "anti" if self._antidominant else "", "" if self._n_numerators else "no", super(CompressionDomain, self)._repr_())

        return "Compression domain with columns {} {}dominant {} extra numerators for {}".format(self._column, "anti" if self._antidominant else "", "" if self._n_numerators else "no", super(CompressionDomain, self)._repr_())
        

    def weights(self):
        r"""
        The set of weights that occur in the symmetric Koornwinder polynomial.

        The highest weight of the Koornwinder polynomial was specified during the creation of ``self``.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: CompressionDomain(2,(1,0)).weights()
            [(1, 0), (0, 0)]
        """
        return [wt for wt in reversed(sorted(self._dt.keys()))]

    def alcove_path_to_filling(self, path, pairs, move=None):
        r"""
        Convert an alcove path to a Lenart filling.
        This needs to get fixed.

        INPUTS::

        - ``path`` -- alcove path as a list of reduced words of finite Weyl elements
        - ``pairs`` -- a list of pairs of the form ``(i, weight)`` where `i` is an index into ``path``
        and ``weight`` is a weight.
        - ``move`` -- (default: None) If not None, a reduced word for a Weyl group element to apply
        before applying ``path[i]``.

        OUTPUT::

        The tuple of weights obtained by converting each pair
        ``(i, weight)`` into ``weight`` acted on by ``path[i]``.
        """
        weights = []
        if not move:
            move = []
        for (i,weight) in pairs:
            weights.append(weight.weyl_action(path[i]+move))
        return tuple(weights)

    def weight_of_filling(self, filling):
        r"""
        The weight of a filling.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2,(1,0))
            sage: wt = d.weights()[-1]; wt
            (0, 0)
            sage: f = d.fillings(wt)[-1]; f
            ((-1, 0), (1, 0))
            sage: d.weight_of_filling(f)
            (0, 0)
        """
        return QQ(1)/2 * sum(x for x in filling)

    def fillings(self, weight):
        r"""
        The fillings of a given weight.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2,(1,0))
            sage: wt = d.weights()[-1]; wt
            (0, 0)
            sage: d.fillings(wt)
            [((0, -1), (0, 1)), ((0, 1), (0, -1)), ((1, 0), (-1, 0)), ((-1, 0), (1, 0))]
        """
        if not weight in self._dt.keys():
            return []
        return self._dt[weight].keys()

    def n_fillings(self, weight):
        r"""
        The number of fillings of a given weight.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2,(1,0))
            sage: wt = d._dt.keys()[-1]; wt
            (1, 0)
            sage: d.n_fillings(wt)
            1
        """
        if weight in self.weights():
            return len(self.fillings(weight))
        return 0

    def show_fillings(self, weight):
        r"""
        A latex string for the fillings of the given weight, together will all alcove paths and Ram-Yip sums.
        """
        if self.n_fillings(weight) == 0:
            return "No fillings of weight {}".format(weight)
        string = "\\subsection{Weight $%s$}\n"%(weight)
        i = 0
        for filling in self._dt[weight].keys():
            string = string + "\\noindent\\begin{flushleft}\n filling %s: %s paths \\\\\n %s  \\end{flushleft}\n"%(i, self.n_paths(filling), self.la_filling(filling))
            string = string + self.display_filling(filling,None,None,self._the_leading_coef)
            i = i + 1
        return string

    def res(self, f):
        r"""
        In the rational function f, send c to 1/v - v, c0 to 1/v0 - v0, and cl to 1/vl - vl
        but leave cz and c2 alone.
        """
        return self.DAHA().double_affine_type().partial_restrict_base_ring_map()(f)

    def weyl_to_signed_permutation(self, w):
        r"""
        Convert a finite type `B_n` Weyl group element to a signed permutation.

        A signed permutation is a list of integers whose absolute values form a permutation.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(3, [2,1,0], no_paths=True)
            sage: w = d.DAHA().dual_double_affine_type().extended_affine_weyl().dual_classical_weyl().an_element(); w
            s1*s2*s3
            sage: d.weyl_to_signed_permutation(w)
            [2, 3, -1]
        """
        images = [self._standard_basis[i].weyl_action(w) for i in self._standard_basis.keys()]
        return [self._single_column_weight_to_word(v)[0] for v in images]

    def extract_factor(self, poly, factor, string, factor_string, power=None):
        r"""
        If factor divides poly, return the pair (poly/factor, string + factor_string);
        otherwise return (poly,string). If power is not None, divide poly by the 
        maximum power (say p) of factor that divides it, and return the pair
        (poly/factor^p, string + factor_string^p).
        """
        factor = factor.numerator()
        if factor.divides(poly):
            poly = (poly/factor).numerator() # the .numerator forces the element into the polynomial ring
            p = 1
            string = string + factor_string
            if power:
                while factor.divides(poly):
                    poly = (poly/factor).numerator()
                    p = p + 1
            if power and p > 1:
                string = string + "^{%s}"%p
        return (poly, string)

    def latex_simple_monomial(self, m):
        r"""
        Latex string for a monomial only involving q, v, vl, v0.

        Based on variable order in polynomial ring.
        """
        def latex_variable_power(var_string, power):
            if power == 0:
                return ""
            if power == 1:
                return var_string
            return var_string + "^{%s}"%(power)

        degrees = m.degrees()
        if sum(degrees) == 0:
            return "1"
        str = latex_variable_power("q", degrees[0]) + latex_variable_power("v", degrees[1]) + latex_variable_power("\\"+"vl ", degrees[2]) + latex_variable_power("v_0", degrees[3])
        if len(degrees) >= 9:
            str = str + latex_variable_power("c_2", degrees[9])
        if len(degrees) >= 10:
            str = str + latex_variable_power("c_z", degrees[10])
        return str

    def latex_term(self, coef, monom, first_term=None):
        str_monom = self.latex_simple_monomial(monom)
        str_coef = "%s"%(coef)
        if coef > 0 and not first_term: # when to use a + sign
            str = "+"
        else:
            str = ""
        if coef == 1:
            return str + str_monom
        elif coef > 0:
            return str + str_coef + str_monom
        elif coef == -1:
            return "-" + str_monom
        else:
            return str_coef + str_monom

    def latex_simple_polynomial(self, f):
        r"""
        Latex string for a polynomial.

        Assumes the polynomial only has the variables q, v, vl, v0, c2, cz, and only signs for coefficients.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(3, [2,1,0], no_paths=True)
            sage: ddat = d.DAHA().dual_double_affine_type()
            sage: q = ddat.q(); v = ddat.parameter('short'); vl = ddat.parameter('long'); v0 = ddat.parameter('doubled');
            sage: print d.latex_simple_polynomial((1 - q**2 * v * vl**2 * v0).numerator())
            -q^{2}v\vl ^{2}v_0+1
            sage: print d.latex_simple_polynomial(1)
            1
        """
        if f in QQ:
            return "{}".format(f)
        monomials = f.monomials()
        coefficients = f.coefficients()
        str = self.latex_term(coefficients[0], monomials[0], True)
        for i in range(1, len(monomials)):
            str = str + self.latex_term(coefficients[i], monomials[i])
        return str

    def latex_factorization(self, fact):
        r"""
        A latex string for a factorized polynomial.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2,(1,0))
            sage: ddat = d.DAHA().dual_double_affine_type()
            sage: q = ddat.q(); v = ddat.parameter('short'); vl = ddat.parameter('long'); v0 = ddat.parameter('doubled');
            sage: f = (1 - vl**2)*(1 - q**2*v**2*v0**2).numerator()
            sage: print d.latex_factorization(f.factor())
            (\vl -1)(\vl +1)(qvv_0-1)(qvv_0+1)
        """
        unit = fact.unit()
        if len(fact) == 1:
            # there is a tendency to factor out a -1. If that happens
            # and there is only one factor with no power, multiply the
            # single factor by the sign.
            factor, power = fact[0]
            if power == 1 and unit == -1:
                return "(" + self.latex_simple_polynomial(fact.expand()) + ")"
        one = self.DAHA().base_ring().one().numerator()
        if unit == one:
            string = ""
        elif unit > 0:
            string = "%s"%(unit)
        else:
            string = "(%s)"%(unit)
        for f in fact:
            factor, power = f
            string = string + "(" + self.latex_simple_polynomial(factor) + ")"
            if power > 1:
                string = string + "^{%s}"%(power)
        return string

    def latex_rational_function(self, f):
        r"""
        Convert a rational function to a latex string.

        This recognizes and factors out the denominators that occur in the Ram-Yip formula
        and also certain numerators which will never combine (those with 'c2' and 'cz').

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2,(2,0), n_numerators=True)
            sage: wt = d.weights()[-1]
            sage: f = d.fillings(wt)[0]
            sage: p = d.paths(f)[0]
            sage: ddat = d.DAHA().dual_double_affine_type()
            sage: print d.latex_rational_function(d.res(p.evaluate(True)))
            q^{2}v_0\dfrac{\chi^s_{3}}{e_{3}}
            sage: d = CompressionDomain(2,(2,0), n_numerators=False)
            sage: print d.latex_rational_function(d.res(p.evaluate(True)))
            q^{2}\dfrac{1}{e_{3}}(-q^{2}v^{2}\vl ^{2}v_0^{2}+q^{2}\vl ^{2}v_0^{2}-v_0^{2}+1)
        """
        num = f.numerator()
        den = f.denominator()

        # find the factors of the denominator
        den_string = ""
        den_pos = []
        for i in range(self._N):
            if self._reflection_denom[i].divides(den):
                (den, den_string) = self.extract_factor(den, self._reflection_denom[i], den_string, "e_{%s}"%(i+1))
                den_pos = den_pos + [i]

        if den_string == "":
            den_string = "1"

        if den == self._base_ring.one() or den == - self._base_ring.one():
            num = num * den # transfer the sign to the numerator
        else: # we have some random junk in the denominator
            den_string = den_string + self.latex_factorization(den.factor())

        ddat = self.DAHA().dual_double_affine_type()
        q = ddat.q()
        v = ddat.parameter('short')
        vl = ddat.parameter('long')
        v0 = ddat.parameter('doubled')

        # find the binomial factors of the numerator coming from Ram-Yip for short roots
        num_string_short = ""
        #for i in den_pos:
        for i in range(self._N):
            if self._the_reduced_word[i] == 0:
                if (self._pos_fold_numerator[i]).divides(num):
                    (num, num_string_short) = self.extract_factor(num, self._pos_fold_numerator[i], num_string_short, "\chi^z_{%s}"%(i+1))
                if (self._neg_fold_numerator[i]).divides(num):
                    (num, num_string_short) = self.extract_factor(num, self._neg_fold_numerator[i], num_string_short, "\chi^2_{%s}"%(i+1))
            elif self._the_reduced_word[i] == self._n and self._n_numerators:
                if (self._pos_fold_numerator[i]).divides(num):
                    (num, num_string_short) = self.extract_factor(num, self._pos_fold_numerator[i], num_string_short, "\chi^0_{%s}"%(i+1))
                    num = num * v.numerator()  # we really divided out by v*\chi_i^0
                if (self._neg_fold_numerator[i]).divides(num):
                    (num, num_string_short) = self.extract_factor(num, self._neg_fold_numerator[i], num_string_short, "\chi^s_{%s}"%(i+1))
                    num = num * v0.numerator() # we really divided out by v0*\chi_i^s

        if num_string_short == "":
            num_string_short = "1"

        # now extract all powers of q, v, vl, v0, and (1 - vl)^2
        # which will be written in front
        num_string = ""
        (num, num_string) = self.extract_factor(num, q, num_string, 'q', power=True)
        (num, num_string) = self.extract_factor(num, v, num_string, 'v', power=True)
        (num, num_string) = self.extract_factor(num, vl, num_string, 'v_\ell', power=True)
        (num, num_string) = self.extract_factor(num, v0, num_string, 'v_0', power=True)
        (num, num_string) = self.extract_factor(num, (1-vl**2), num_string, '(1-v_\ell^2)', power=True)
        string = num_string + "\dfrac{%s}{%s}"%(num_string_short,den_string)
        if num != self._base_ring.one().numerator():
            string = string + "%s"%(self.latex_factorization(num.factor()))
        return string

    def paths(self, filling):
        r"""
        A list of alcove paths for the given filling.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2,(2,0))
            sage: wt = d.weights()[1]; wt
            (1, 1)
            sage: f = d.fillings(wt)[-1]; f
            ((1, 0), (1, 0), (0, 1), (0, 1))
            sage: d.paths(f)
            [Ram-Yip term: start: 1 end: t[2*Lambda[2]] * s2*s1 folds: [0, 0, 0, 1, 0] signs: [-1, -1, -1, -1, -1] value: (-1) * cl * v * v0^2 * q^2 * (q*v*vl*v0 - 1)^-1 * (q*v*vl*v0 + 1)^-1, Ram-Yip term: start: 1 end: t[2*Lambda[2]] * s2*s1*s2 folds: [0, 0, 1, 1, 0] signs: [-1, -1, -1, 1, -1] value: cl * v0 * q^2 * (q*v*vl*v0 - 1)^-1 * (q*v*vl*v0 + 1)^-1 * (q^2*v*vl^2*v0 - 1)^-1 * (q^2*v*vl^2*v0 + 1)^-1 * (q^2*v*vl^2*v0*c + c0), Ram-Yip term: start: 1 end: t[2*Lambda[2]] * s2*s1*s2 folds: [0, 1, 1, 0, 0] signs: [-1, -1, 1, -1, -1] value: v^-1 * cl * q^2 * (q*vl - 1)^-1 * (q*vl + 1)^-1 * (q^2*v*vl^2*v0 - 1)^-1 * (q^2*v*vl^2*v0 + 1)^-1 * (q^2*v*vl^2*v0*c0 + c)]
        """
        weight = self.weight_of_filling(filling)
        if weight not in self.weights():
            return []
        if filling not in self.fillings(weight):
            return []
        return self._dt[weight][filling]

    def n_paths(self, filling):
        """
        The number of Ram-Yip paths for the given filling.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2,(2,0))
            sage: wt = d.weights()[-1]; wt
            (0, 0)
            sage: f = d.fillings(wt)[-1]; f
            ((1, 0), (-1, 0), (0, 1), (0, -1))
            sage: d.n_paths(f)
            3
        """
        return len(self.paths(filling))

    def add_terms(self, filling, indices=None, normalizer=None):
        """
        Add the Ram-Yip terms for a filling.

        INPUTS::

        - ``filling`` -- a filling
        - ``indices`` -- (default: None) A list of indices used to select a subset
        of alcove paths for the given filling. If None, sum over all alcove paths.
        - ``normalizer`` -- (default: None) a coefficient to divide by.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2,(2,0))
            sage: wt = d.weights()[-1]
            sage: f = d.fillings(wt)[-1]
            sage: d.add_terms(f).factor()
            (vl - 1) * (vl + 1) * v0^2 * vl^2 * q^4 * (q*vl - 1)^-1 * (q*vl + 1)^-1 * (q*v*vl^2*v0 - 1)^-1 * (q*v*vl^2*v0 + 1)^-1 * (q*v*vl^2*v0*c2 + cz) * (q^3*v*vl^2*v0 - 1)^-1 * (q^3*v*vl^2*v0 + 1)^-1 * (q^3*v*vl^2*v0*c2 + cz)
        """
        wt = self.weight_of_filling(filling)
        if wt not in self._dt.keys():
            raise ValueError("{} does not appear as a weight".format(wt))
        if filling not in self._dt[wt].keys():
            raise ValueError("{} does not appear as a filling".format(filling))
        if indices == None:
            indices=range(len(self._dt[wt][filling]))
        ans = self.res(self._base_ring.sum([self._dt[wt][filling][i].evaluate(True) for i in indices]))
        if normalizer is not None:
            return ans/normalizer
        else:
            return ans

    def sum_str(self, filling, indices=None,normalizer=None):
        """
        A latex string for the sum of Ram-Yip terms.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2,(2,0))
            sage: wt = d.weights()[-1]
            sage: f = d.fillings(wt)[-1]
            sage: print d.sum_str(f)
            q^{4}v_\ell^{2}v_0^{2}(1-v_\ell^2)\dfrac{\chi^2_{1}\chi^2_{5}}{e_{1}e_{2}e_{5}}
        """
        return self.latex_rational_function(self.add_terms(filling,indices,normalizer))

    def path_group_display(self, filling, indices=None, normalizer=None):
        r"""
        Latex string to display a subset of alcove paths belonging to the given filling
        and the sum of their terms in the Ram-Yip formula.
        """
        if indices is None:
            indices = [i for i in range(self.n_paths(filling))]
        # latex the RY paths
        string = "\\begin{flushleft}\n"
        wt = self.weight_of_filling(filling)
        for i in indices:
            term = self._dt[wt][filling][i]
            string = string + "folds: %s \\quad index: %s start: $%s$ end: $%s$ \\\\\n %s\\\\\n"%(self.la_RY(term), i, self.la_weyl(term._start), self.la_weyl(term._end*self._u_weight), self.la_path(term))
        string = string + "\\end{flushleft}\n"
        # latex the rational function given by the sum of the terms
        string = string + "$$\n%s\n$$\n"%self.sum_str(filling, indices, normalizer)
        return string
        
    def display_filling(self, filling, groups=None, group_headers=None,normalizer=None):
        r"""
        A latex string for the grouped paths belonging to a fixed filling.

        INPUTS::

        - ``filling`` -- a filling
        - ``groups`` -- (default: None) a list of lists. Each sublist is a list of indices into the
        paths for the given filling. Each such group of alcove paths is printed, together with the
        sum of corresponding Ram-Yip terms. If None, the entire set of indices is used as one group.
        - ``group_headers`` -- (default: None) a list of strings to be included at the beginning of the
        corresponding groups
        """
        string = ""
        if groups is None:
            groups = [range(self.n_paths(filling))]
        elif len(groups) > 1:
            string = "\\noindent SPLIT\n"
        for i in range(len(groups)):
            if group_headers is not None:
                string = string + group_headers[i]
            string = string + self.path_group_display(filling, groups[i],normalizer)
        return string

    def extended_filling_from_path(self, term):
        """
        Ram Yip term converted into a filling but with every step remembered.
        """
        return self.alcove_path_to_filling(term.alcove_path_directions(), self._pair, self._u_weight_reduced_word)

    def latex_header(self):
        r"""
        Latex header for the compression data.
        """
        str = "\\section{$n=%s$, column sequence $%s$}\n"%(self._n, self._column)
        str = str + "$$\\begin{array}{|r|c|l|l|l|} \\hline \n \\text{index} & &  \\text{root} & \\text{eigenvalue} & \\text{denominator} \\\\\n \\hline \\hline"
        for i in range(self._N):
            str = str + "%s &s_{%s} & %s \delta - %s & %s & %s \\\\\n"%(i+1, self._the_reduced_word[i], self._inversion_depths[i], self._inversions_classical[i], self.latex_simple_polynomial(self._eigen[i].numerator()), self.latex_simple_polynomial(self._reflection_denom[i].numerator()))
        str = str + "\\hline \\end{array}\n$$\n"
        str = str + "$$\\begin{array}{|r|c|l|l|} \\hline \n \\text{index} & &  \\text{$\chi^z$ or $\chi^0$} & \\text{$\chi^2$ or $\chi^s$} \\\\\n \\hline \\hline"
        for i in range(self._N):
            k = self._the_reduced_word[i]
            if k == 0:
                str = str + "%s &s_{%s} & %s & %s \\\\\n"%(i+1, k, self.latex_simple_polynomial(self._pos_fold_numerator[i]), self.latex_simple_polynomial(self._neg_fold_numerator[i]))
            elif k == self._n and self._n_numerators:
                str = str + "%s &s_{%s} & %s & %s \\\\\n"%(i+1, k, self.latex_simple_polynomial(self._pos_fold_numerator[i]), self.latex_simple_polynomial(self._neg_fold_numerator[i]))
        str = str + "\\hline \\end{array}\n$$\n"

        str = str + "sampling pairs: %s\n"%(self._sampling_pair)

        return str

    def latex_start_file(self):
        r"""
        The start of an output latex file for compression data.
        """
        return "\\documentclass{amsart}\n\\usepackage{tabmac}\n\\newcommand{\\vl}{v_\ell}\n\\begin{document}\n" + self._repr_(True) + "\n"

    def latex_end_file(self):
        r"""
        The start of an output latex file for compression data.
        """
        return "\\end{document}\n"

    def latex_all(self, filename):
        r"""
        A latex file for the compression data.
        """
        with open(filename, 'w') as fp:
            print >>fp, self.latex_start_file() + self.latex_header()
            for w in self.weights():
                print >>fp, self.show_fillings(w)
            print >>fp, self.latex_end_file()

    def latex_try_pairs(self, filling, filename):
        """
        Latex file for trying to combine all pairs of paths associated with a given filling.
        """
        p = self.paths(filling)
        npaths = len(p)
        with open(filename, 'w') as fp:
            print >>fp, self.latex_start_file()
            for j in range(npaths):
                for k in range(j+1, npaths):
                    print >>fp, self.display_filling(filling, [[j,k]])
            print >>fp, self.latex_end_file()
        
    def la_fold(self, fold, sign):
        """latex output for a fold with sign"""
        if sign == -1:
            return "\underline{%s}"%fold
        else:
            return "%s"%fold

    def la_weyl(self, w):
        """
        Input is extended affine Weyl group element but it actually lives in the finite
        Weyl group.
        """
        rw = w.to_classical_weyl().reduced_word()
        if len(rw) == 0:
            return "1"
        str = ""
        for i in rw:
            str = str + "s_%s"%i
        return str

    def la_RY(self, ry):
        """
        Latex string for a Ram-Yip path.

        Folds are indicated by a 0 or 1 where 1 means fold = omit reflection.
        Underline means the fold is negative.

        """
        str = ""
        for i in range(len(ry._folds)):
            str = str + self.la_fold(ry._folds[i], ry._signs[i])
        return str

    def _single_column_weight_to_word(self, wt):
        r"""
        The word associated to a single column weight.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2, [0,1])
            sage: w = d.weights()[1]
            sage: f = d.fillings(w)[1]
            sage: [d._single_column_weight_to_word(wt) for wt in f]
            [[1, -2], [1, 2]]
        """
        return [i+1 if wt[i]>0 else -(i+1) for i in wt.support()]

    def filling_to_tableau(self, wts):
        r"""
        A rectangular tableau from a filling.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2, [1,1])
            sage: w = d.weights()[2]
            sage: f = d.fillings(w)[1]; f
            ((-1, 1), (1, 1), (1, 0), (1, 0))
            sage: d.filling_to_tableau(f)
            [[-1, 1, 1, 1], [2, 2]]
        """
        tab = dict()
        for r in range(self._column[0]):
            tab[r] = []
        for wt in wts:
            col = self._single_column_weight_to_word(wt)
            for r in range(len(col)):
                tab[r] = tab[r] + [col[r]]
        return [tab[r] for r in range(self._column[0])]

    def la_filling(self, wts):
        r"""
        Latex string for a filling.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2, [1,1])
            sage: w = d.weights()[2]
            sage: f = d.fillings(w)[1]; f
            ((-1, 1), (1, 1), (1, 0), (1, 0))
            sage: d.filling_to_tableau(f)
            [[-1, 1, 1, 1], [2, 2]]
            sage: print d.la_filling(f)
            $\tableau[s]{\overline{1}&1&1&1\\ 2&2\\ }$
        """
        return self.la_tableau(self.filling_to_tableau(wts))

    def la_tableau(self, tab):
        r"""
        Latex string for a tableau.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2, [0,1], no_paths=True)
            sage: print d.la_tableau([[1,1,2],[3,4,4]])
            $\tableau[s]{1&1&2\\ 3&4&4\\ }$
        """
        str = "$\\tableau[s]{"
        for r in range(len(tab)):
            for i in range(len(tab[r])):
                if i > 0:
                    str = str + "&"
                str = str + self.la_entry(tab[r][i])
            str = str + "\\\\\n"
        str = str + "}$"
        return str

    def la_entry(self, x):
        if x == 0:
            return "\\bl"
        if x > 0:
            return "%s"%x
        if x < 0:
            return "\overline{%s}"%(-x)

    def la_path(self, path):
        r"""
        Latex string for a path.

        The path is given as a sequence of signed permutations.

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2, [0,1])
            sage: w = d.weights()[1]
            sage: f = d.fillings(w)[1]
            sage: print d.la_path(d.paths(f)[1])
            $\tableau[s]{2&2&1&\overline{1}\\ \overline{1}&\overline{1}&\overline{2}&\overline{2}\\ }$
        """
        return self.la_tableau(self.path_to_signed_perm_tableau(path))

    def path_to_signed_perm_tableau(self, path):
        r"""
        A rectangular tableau from a path.

        The columns are signed permutations computed from the elements of the
        alcove path.

        EXAMPLES::

            sage: from sage.algebras.compression_domain import CompressionDomain
            sage: d = CompressionDomain(2, [0,1])
            sage: w = d.weights()[1]
            sage: f = d.fillings(w)[1]; f
            ((1, -1), (1, 1))
            sage: [d.path_to_signed_perm_tableau(p) for p in d.paths(f)]
            [[[2, 2, 2, -2], [-1, -1, -1, -1]], [[2, 2, 1, -1], [-1, -1, -2, -2]], [[2, 2, -1, -1], [-1, 1, -2, -2]]]
        """
        tab = dict()
        for r in range(self.rank()):
            tab[r] = []
        W = self.classical_weyl_Y()
        signed_perms = [self.weyl_to_signed_permutation(W.from_reduced_word(w)) for w in path.alcove_path_directions()]
        for p in signed_perms:
            for r in range(self.rank()):
                tab[r] = tab[r] + [p[r]]
        return [tab[r] for r in range(self.rank())]

    def split(self, filling, positions):
        r"""
        Partition the alcove paths for a given filling by extra sampling.

        The alcove paths for the given filling are sampled at specified positions (the weights
        are predetermined by the positions) producing a map from paths to tuples of weights
        we will call mini-fillings. This data is returned in the form of a :class:`Family`
        whose keys are mini-fillings and whose values are the list of indices of the paths
        which map to the given mini-filling.
        """
        paths = self.paths(filling)
        pairs = [self._pair[i] for i in positions]
        path_dict = dict()
        for i in range(len(paths)):
            val = self.alcove_path_to_filling(paths[i].alcove_path_directions(), pairs)
            if val in path_dict.keys():
                path_dict[val] = path_dict[val] + [i]
            else:
                path_dict[val] = [i]
        return Family(path_dict)

    def merge(self, fam, part):
        r"""
        Merge groups of alcove paths given by :meth:`split`, 

        INPUTS::

        - ``fam`` -- a :class:`Family` provided by :meth:`split`.
        - ``part`` -- a listlist of integers

        :meth:`split` defines a partition of the set of paths with a given filling
        whose blocks are labeled by mini-fillings. Its values are lists of integers that specify
        subsets of alcove paths. Say this partition has `N` blocks.
        Then ``part`` should define a partition of the set ``range(N)``. It specifies which blocks
        to merge. The resulting object is a :class:`Family` whose keys are tuples of mini-fillings
        and whose values are lists of integers which define the merged blocks.
        """
        def my_union(listlist):
            s = Set()
            for lis in listlist:
                s = s.union(Set(lis))
            return s
        keys = fam.keys()
        return Family(dict([[Set([keys[i] for i in block]), my_union([fam[keys[i]] for i in block])] for block in part]))

    def latex_minifillings(self, fillings):
        r"""
        Latex string for a list of fillings.
        """
        if not isinstance(fillings,(Set_object,list)):
            fillings = [fillings]
        string = "\\noindent\\begin{flushleft}\n "
        for i in range(len(fillings)):
            string = string + self.la_filling(fillings[i]) + "\\qquad %s"%i + "\\\\\n"
        string = string + "\\end{flushleft}\n "
        return string

    def latex_family(self, filling, fam):
        r"""
        Latex string for the above family.
        """
        keys = fam.keys()
        group_headers = [self.latex_minifillings(k) for k in keys]
        return self.display_filling(filling, [fam[k] for k in keys], group_headers, self._the_leading_coef)

