"""
Generalized Hidden Markov Models.

Generalized Hidden Markov Models (GHMMs) are generalizations of ordinary Hidden Markov Models
of the output-from-transition type.
"""

from sage.structure.element import is_Matrix
from sage.matrix.constructor import matrix, Matrix
from sage.matrix.special import diagonal_matrix, identity_matrix, block_matrix
from sage.modules.free_module_element import vector
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.stats.hmm import hmm
from sage.misc.prandom import random
from sage.graphs.digraph import DiGraph
from copy import copy

def _initials_equivalent(Ts, init1, init2):
    """
    Return ``True`` if the two initial distributions on states 
    give the same process.

    INPUT:

    - ``Ts`` -- a list of matrices, each of the same dimension,
       and summing to a stochastic matrix.

    - ``init1`` -- a list or vector of probabilities, with as many elements as ``Ts``.
    - ``init2`` -- a list or vector of probabilities, with as many elements as ``Ts``.

    ALGORITHM: Adapted from [Bal1993]_.
    """
    n = len(init1)
    V = [matrix(ZZ, 0, n, [])] # array so it can be updated by a subroutine; ugly Python trick
    rows = []
    def new_vec(v):
        # Consider a new vector which must lie in N
        # Return False if the new vector or some descendent
        # is not in N, else return True
        if v in V[0].row_space(): # was v in span(V[0]):
            return True
        if not sum(v).is_zero():
            return False
        if len(rows) >= n - 1:
            return True # Can this happen?  Everything in V is orthogonal to [1,1,1 ...].
        rows.append(v)
        V[0] = matrix(rows)
        for T in Ts:
            if not new_vec(v*T):
                return False
        return True
    return new_vec(vector(init2) - vector(init1))

def _hmm_to_ghmm_data(discrete_hmm):
    """
    Return the transition matrices and initial distribution for
    a GeneralizedHiddenMarkovModel equivalent to the given
    DiscreteHiddenMarkovModel, in the sense that every output sequence
    is an equally likely output of both.  This is an internal function
    used in the constructor for :class:`GeneralizedHiddenMarkovModel`.

    Given a DiscreteHiddenMarkovModel with transition matrix `A`, output matrix
    `b`, and initial state distribution `pi`, we could almost make do with
    the same set of states, same initial distribution, and transition matrices 
    ``[A * diagonal_matrix(c) for c in B.columns()]``.  However, this fails to take
    into account that the HMM emits its first output symbol before transitioning
    to a new state, while the GHMM transitions and emits together.  So, we add an
    extra 'startup' state to the GHMM, which outputs with the correct distribution
    and transitions appropriately to the regular states, which are as described
    above.

    EXAMPLES::

        sage: M = hmm.DiscreteHiddenMarkovModel([[0.6,0.2,0.2],[0.4,0.1,.5],[.6,.2,.2]], [[0.5,.5],[.1,.9],[0.9,0.1]], [0.3,0.7,0])
        sage: G = hmm.GeneralizedHiddenMarkovModel(M) # indirect doctest
        sage: G.transition_matrix()
        [0.6 0.2 0.2 0.0]
        [0.4 0.1 0.5 0.0]
        [0.6 0.2 0.2 0.0]
        [0.3 0.7 0.0 0.0]
        sage: G
        Generalized Hidden Markov Model with 4 states and 2 output symbols
        sage: G.canonical_form()
        Generalized Hidden Markov Model with 3 states and 2 output symbols
    """
    if not isinstance(discrete_hmm, hmm.DiscreteHiddenMarkovModel):
        raise TypeError
    pi = vector(discrete_hmm.initial_probabilities())
    A = discrete_hmm.transition_matrix()
    B = discrete_hmm.emission_matrix()
    z = vector([0]*(1+A.nrows()))
    Ts = [(A * Bi).stack(pi*Bi).augment(z) for c in B.columns() for Bi in [diagonal_matrix(c)]]
    return Ts, vector([0]*A.nrows() + [1])


def _digraph_to_ghmm_data(g):
    r"""
    Return a list of transition matrices for a GHMM corresponding to digraph g.

    ``g``'s vertices are states, and there must be some integer ``m`` so that
    each edges is labelled by a list of ``m`` numbers, giving the pseudo-probability
    of emitting outputs 0, 1, ..., ``m-1`` when traversing that edge.

    INPUT:

    - ``g`` -- a DiGraph

    OUTPUT:

    A list ``L`` of transition matrices, and a dictionary mapping the graph vertices to
    ``range(0, len(L))`` (which can be used to transform data associated with the graph
    into information about the states of the GHMM associated with ``L``).

    EXAMPLES::

        sage: D = DiGraph([(0, 1, [1/2, 1/3, 1/6]), (1, 0, [1/2, 0, 1/2])])
        sage: hmm.GeneralizedHiddenMarkovModel(D).Ts[1] # indirect doctest
        [ 0  1/3]
        [ 0    0]
        sage: D = DiGraph([('x', 'x', [1/2, 1/3, 1/6])], loops=True)
        sage: hmm.GeneralizedHiddenMarkovModel(D).Ts # indirect doctest
        [[1/2], [1/3], [1/6]]

    TESTS::

        sage: hmm.GeneralizedHiddenMarkovModel(DiGraph()) # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: must have at least one state
        sage: D = DiGraph([('a', 'b', [1]), ('c', 'd', [0,1])])
        sage: hmm.GeneralizedHiddenMarkovModel(D)
        Traceback (most recent call last):
        ...
        ValueError: all labels must have the same length
        sage: D = DiGraph(); D.add_vertex(0)
        sage: hmm.GeneralizedHiddenMarkovModel(D)
        Traceback (most recent call last):
        ...
        ValueError: no edges, cannot derive a GHMM
    """
    G = copy(g)
    f = G.relabel(return_map=True)
    n = G.num_verts()
    if n == 0:
        raise ValueError('must have at least one state')
    k = -1
    for (i,j,lab) in G.edge_iterator():
        if k == -1:
            k = len(lab)
        elif k != len(lab):
            raise ValueError('all labels must have the same length')
    if k == -1:
        raise ValueError('no edges, cannot derive a GHMM')
    Ts = [ [[0]*n for _ in range(n) ] for _ in range(k)]
    for (i,j,lab) in G.edge_iterator():
        for o,p in enumerate(lab):
            Ts[o][i][j] = p
    return [matrix(m) for m in Ts], f

class GeneralizedHiddenMarkovModel:
    r"""
    A GeneralizedHiddenMarkovModel (GHMM) represents a Generalized
    Hidden Markov Model of the output-from-transition type. 
    If there are `n` states and `k` outputs, such a model is characterized by
    a list of `k` `n \times n` matricies `T_i` and an initial distribution
    (stochastic vector of length `k`).

    Each `T_i` is an `n \times n` matrix, where `T_i[j,k]` is the (pseudo)-probability of
    transitioning from state `j` to state `k` and emitting symbol `i`.

    The overall transition matrix `T` is the sum of the `T_i`.  Then `T[j,k]`
    is the pseudo-probability of transitioning from state `j` to state `k`.

    `T` must be stochastic; that is, every row should sum to 1.

    If the initial distribution is not provided, it is computed as the left
    eigenvector of `T` with eigenvalue 1.

    The `T_i` may have negative entries, which clearly cannot be interpreted as
    probabilities.  A GHMM is proper if for every output sequence, the calculated
    likelihood of this sequence (see :meth:`likelihood`) is nonnegative.
    In this case, we do indeed get a probability distribution over outputs 
    of any fixed length.

    EXAMPLES::

        sage: T0 = Matrix([[  0, 1/2],[1/2,   0]])
        sage: T1 = Matrix([[1/2,   0],[1/8, 3/8]])
        sage: G = hmm.GeneralizedHiddenMarkovModel([T0, T1]); G
        Generalized Hidden Markov Model with 2 states and 2 output symbols
        sage: G.transition_matrix()
        [1/2 1/2]
        [5/8 3/8]
        sage: G.initial_probabilities()
        (5/9, 4/9)
        sage: G.initial_probabilities() == G.initial_probabilities() * G.transition_matrix()
        True

    The "Golden Mean Process" (from [Upp1997]_).  This process always emits
    a 1 following a 0, and after a 1 emits 0 or 1 each with probability 1/2::

        sage: T0 = Matrix([[0, 1/2], [0, 0]])
        sage: T1 = Matrix([[1/2, 0], [1, 0]])
        sage: GMP = hmm.GeneralizedHiddenMarkovModel([T0, T1])
        sage: GMP.likelihood([0,0]), GMP.likelihood([1,1])
        (0, 1/3)

    An equivalent presentation with negative entries, also from [Upp1997]_::

        sage: S0 = Matrix([[1, -1/2], [2, -1]])
        sage: S1 = Matrix([[1/2,  0], [0,  0]])
        sage: GMP2 = hmm.GeneralizedHiddenMarkovModel([S0, S1], [4/3, -1/3])
        sage: GMP2.likelihood([0,0]), GMP2.likelihood([1,1])
        (0, 1/3)
        sage: GMP.is_same_process_as(GMP2)
        True

    A GHMM can be derived from a DiscreteHiddenMarkovModel::

        sage: M = hmm.DiscreteHiddenMarkovModel([[0.6,0.4],[0.9,0.1]], [[0.5,.5],[.1,.9]], [0.3,0.7])
        sage: G = hmm.GeneralizedHiddenMarkovModel(M)
        sage: G.transition_matrix()
        [0.6 0.4 0.0]
        [0.9 0.1 0.0]
        [0.3 0.7 0.0]
        sage: exp(M.log_likelihood([0,1,0,1]))
        0.03503412
        sage: G.likelihood([0,1,0,1])
        0.03503412
    """
    def __init__(self, data, initial_distribution = None):
        r"""
        Create a Generalized Hidden Markov Model of the output-from-transition type.

        INPUT:

        - ``data`` -- Data defining the states and transition matricies of the GHMM.
           This can be a list of transition matrices, a DiscreteHiddenMarkovModel,
           a GeneralizedHiddenMarkovModel, or a DiGraph of the sort produced by
           :meth:`graph`.

        - ``initial_distribution`` -- None (the default), or a list or vector of
          preudo-probabilities.

        EXAMPLES::

            sage: T0 = Matrix([[  0, 1/2],[1/2,   0]])
            sage: T1 = Matrix([[1/2,   0],[1/8, 3/8]])
            sage: G = hmm.GeneralizedHiddenMarkovModel([T0, T1], [1/2,1/2]); G
            Generalized Hidden Markov Model with 2 states and 2 output symbols

        When the defining data is a list of transition matricies, the initial
        distribution can be omitted, in which case it will be the left egenvector
        for eigenvalue 1 for the transition matrix (i.e. sum of the listed matrices).

        ::

            sage: G = hmm.GeneralizedHiddenMarkovModel([T0, T1])
            sage: G.initial_probabilities()
            (5/9, 4/9)

        If defining data is an instance of :class:GeneralizedHiddenMarkovModel,
        the transition matricies are copied. If ``initial_distribution`` is given
        then it is used, otherwise the initial distribution of the data is used.

        ::

            sage: hmm.GeneralizedHiddenMarkovModel(G).initial_probabilities()
            (5/9, 4/9)
            sage: hmm.GeneralizedHiddenMarkovModel(G,[1/3,2/3]).initial_probabilities()
            (1/3, 2/3)

        If the defining data is an instance of :class:hmm.DiscreteHiddenMarkovModel,
        an equivalent GHMM is determined.  This has 1 more state than the HMM, as may be
        needed in the general case (but the canonical form of this GHMM may have
        fewer states).  In this case, ``initial_distribution`` must be None and an
        initial distribution corresponding to the DiscreteHiddenMarkovModel's initial
        distribution is computed.

        ::

            sage: M = hmm.DiscreteHiddenMarkovModel([[0.5,0.5],[0.25,0.75]], [[0.5,0.5],[1,0]], [0,1])
            sage: G = hmm.GeneralizedHiddenMarkovModel(M); G
            Generalized Hidden Markov Model with 3 states and 2 output symbols
 
        If the defining data is a DiGraph ``D``, then ``D``'s vertices are
        states, there must be some integer ``m`` so that each edges is labelled
        by a list of ``m`` numbers, giving the probability of emitting
        outputs when traversing that edge.  ``D``'s states must be [0, 1, ..., k]
        for some ``k``.  The initial distribution can be omitted, in which case it
        will be computed as the left egenvector for eigenvalue 1 for the transition matrix.
        Otherwise, it should be a dictionary mapping vertices to pseudo-probabilities.

        ::

            sage: D = DiGraph([('a', 'b', [1/2, 1/3, 1/6]), ('b', 'a', [1/2, 0, 1/2])])
            sage: G = hmm.GeneralizedHiddenMarkovModel(D, {'a': 1}); G
            Generalized Hidden Markov Model with 2 states and 3 output symbols
            sage: G.initial_probabilities()
            (1, 0)
            sage: hmm.GeneralizedHiddenMarkovModel(D)
            Generalized Hidden Markov Model with 2 states and 3 output symbols
            sage: D = DiGraph([(0, 0, [1/2, 1/3, 1/6])], loops=True)
            sage: hmm.GeneralizedHiddenMarkovModel(D)
            Generalized Hidden Markov Model with 1 state and 3 output symbols
        """
        if isinstance(data, hmm.DiscreteHiddenMarkovModel):
            if initial_distribution is not None:
                raise ValueError('must not provide an initial distribution when giving a DiscreteHiddenMarkovModel')
            # Convert a state-output HMM into an edge-output GHMM
            T_matrices, initial_distribution = _hmm_to_ghmm_data(data)
        elif isinstance(data, GeneralizedHiddenMarkovModel):
            # copy data from the input
            T_matrices = data.Ts
            if initial_distribution is None:
                initial_distribution = data.initial_probabilities()
        elif isinstance(data, DiGraph):
            T_matrices, f = _digraph_to_ghmm_data(data)
            if initial_distribution is not None:
                # f gives the relabelling of vertices to states
                if not isinstance(initial_distribution, dict):
                    raise ValueError('initial distribution for a digraph must be a dictionary')
                d = [0]* (T_matrices[0].nrows()) # number of states
                for s, p in initial_distribution.items():
                    d[f[s]] = p
                initial_distribution = d
        else:
            T_matrices = data
        self.n_outputs = len(T_matrices)
        self.n_states = T_matrices[0].nrows()
        self._identity_matrix = identity_matrix(self.n_states)
        self.Ts = T_matrices
        for T in T_matrices:
            if T.dimensions() != (self.n_states, self.n_states):
                raise ValueError('T matrices must all be square and of the same size')
        T = sum(T_matrices) # overall transition matrix
        self._transition_matrix = T
        if T.base_ring().is_exact():
            v = vector([1] * self.n_states)
            if T * v != v:
                raise ValueError('transition matrix must be stochastic')
        if initial_distribution:
            self._initial_distribution = vector(initial_distribution)
            if self._initial_distribution.base_ring().is_exact() and \
               sum(initial_distribution) != 1:
                raise ValueError('initial distribution must be stochastic')
        else:
            # left eigenvector for 1
            def steady_state_probs(transition_matrix):
                M = (transition_matrix - 1).left_kernel().basis_matrix()
                if M.nrows() < 1:
                    raise ValueError('not a stochastic matrix')
                if M.nrows() > 1:
                    raise ValueError('disconnected transition matrix, no unique initial probability distribution')
                r = M.row(0)
                return r / sum(r) # normalize to sum 1
            self._initial_distribution = steady_state_probs(T)

    def __repr__(self):
        """
        Return string representation of this generalized hidden Markov model.

        EXAMPLES::

            sage: T = identity_matrix(3)/2
            sage: m = hmm.GeneralizedHiddenMarkovModel([T,T], [0.1, 0.2, 0.7])
            sage: m.__repr__()
            'Generalized Hidden Markov Model with 3 states and 2 output symbols'
            sage: m = hmm.GeneralizedHiddenMarkovModel([2*T], [0.1, 0.4, 0.4])
            sage: m.__repr__()
            'Generalized Hidden Markov Model with 3 states and 1 output symbol'
        """
        def plural(n, str):
            return ('%d %s' if n==1 else '%d %ss')%(n,str)
        return "Generalized Hidden Markov Model with %s and %s"%(
            plural(self.n_states, 'state'), plural(self.n_outputs, 'output symbol'))

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: m = hmm.GeneralizedHiddenMarkovModel([identity_matrix(2)], [1/3,2/3])
            sage: loads(dumps(m)) == m
            True
        """
        return unpickle_ghmm_v0, (self.Ts, self._initial_distribution)

    def __eq__(self, other):
        """
        Return `True` if `self` equals `other`.

        EXAMPLES::

           sage: T0 = Matrix([[1,4],[3,2]])/10
           sage: T1 = Matrix([[3,2],[1,4]])/10
           sage: G0 = hmm.GeneralizedHiddenMarkovModel([T0, T1], [1/2,1/2])
           sage: G1 = hmm.GeneralizedHiddenMarkovModel([T1, T0], [1/2,1/2])
           sage: G0 == G0 # indirect doctest
           True
           sage: G0.__eq__(T0)
           False
           sage: G0 == G1 # indirect doctest
           False
        """
        return isinstance(other, GeneralizedHiddenMarkovModel) \
            and self.Ts == other.Ts \
            and self._initial_distribution == other._initial_distribution

    def __ne__(self, other):
        """
        Return `True` if `self` equals `other`.

        EXAMPLES::

           sage: T0 = Matrix([[1,4],[3,2]])/10
           sage: T1 = Matrix([[3,2],[1,4]])/10
           sage: G0 = hmm.GeneralizedHiddenMarkovModel([T0, T1], [1/2,1/2])
           sage: G1 = hmm.GeneralizedHiddenMarkovModel([T1, T0], [1/2,1/2])
           sage: G0 != G0 # indirect doctest
           False
           sage: G0.__ne__(T0)
           True
           sage: G0.__ne__(G1)
           True
        """
        return not (self == other)

    def transition_matrix(self):
        """
        Return the transition matrix of this GHMM.

        EXAMPLES::

            sage: from sage.stats.hmm.ghmm import GeneralizedHiddenMarkovModel
            sage: T0 = Matrix([[  0, 1/2],[1/2,   0]])
            sage: T1 = Matrix([[1/2,   0],[1/8, 3/8]])
            sage: G = GeneralizedHiddenMarkovModel([T0, T1])
            sage: G.transition_matrix()
            [1/2 1/2]
            [5/8 3/8]
        """
        return matrix(self._transition_matrix.rows()) # make a copy

    def initial_probabilities(self):
        """
        Return the initial probability distribution on the states of this GHMM.

        EXAMPLES::

            sage: T0 = Matrix([[  0, 1/2],[1/2,   0]])
            sage: T1 = Matrix([[1/2,   0],[1/8, 3/8]])
            sage: G = hmm.GeneralizedHiddenMarkovModel([T0, T1])
            sage: G.initial_probabilities()
            (5/9, 4/9)
            sage: G.initial_probabilities() == G.initial_probabilities() * G.transition_matrix()
            True
        """
        return self._initial_distribution

    def T(self, word):
        r"""
        Return an `n \times n` matrix `T`, where `T[j,k]` is the
        pseudo-probability of transitioning from state `j` to state
        `k` in ``len(word)`` steps, and emitting symbols ``word[0]``,
        ``word[1]`` ... ``word[-1]``.

        EXAMPLES::

            sage: from sage.stats.hmm.all import GeneralizedHiddenMarkovModel
            sage: T0 = Matrix([[  0, 1/2],[1/2,   0]])
            sage: T1 = Matrix([[1/2,   0],[1/8, 3/8]])
            sage: G = GeneralizedHiddenMarkovModel([T0, T1])
            sage: G.T([])
            [1 0]
            [0 1]
            sage: G.T([0,1])
            [1/16 3/16]
            [ 1/4    0]
        """
        return prod([self.Ts[i] for i in word], self._identity_matrix)

    def likelihood(self, word, after=None):
        r"""
        Return the probability of this GHMM emitting the given word.  If ``after`` is
        supplied, return the conditional probability of emitting the word immediately
        following the emission of ``after``; this is the same as

            self.likelihood(after+word) / self.likelihood(after)

        INPUT:

        - ``word`` -- a sequence of outputs.
        - ``after`` -- None (the default) or a sequence of outputs.

        OUTPUT:

        A probability.  Note that these "probabilities" add to 1 over all outputs
        of a given length, not over all outputs

        EXAMPLES::

            sage: from sage.stats.hmm.all import GeneralizedHiddenMarkovModel
            sage: T0 = Matrix([[1/4, 1/2],[1/2, 1/4]])
            sage: T1 = Matrix([[1/4,   0],[1/8, 1/8]])
            sage: G = GeneralizedHiddenMarkovModel([T0, T1], [1, 0])
            sage: G.likelihood([0,0,0])
            27/64
            sage: G.likelihood([1,0,1])
            3/64
            sage: G.likelihood([1,0,1,0]) == sum(G.initial_probabilities() * G.T([1,0,1,0]))
            True
            sage: G.likelihood([0], after=[0,0])
            3/4
            sage: G.likelihood([1], after=[0,0])
            1/4
        """
        v = self._initial_distribution
        d = 1
        if after is not None:
            for s in after:
                v = v * self.Ts[s]
            d = sum(v)
        for s in word:
            v = v * self.Ts[s]
        return sum(v) / d

    def initials_equivalent(self, init1, init2 = None):
        """
        Return True if the two initial distributions on states give the same process.

        INPUT:

        - ``init1`` -- a list of probabilities
        - ``init2`` -- a second list of probabilities.  If ``init2`` is
          omitted or is None, then ``self``'s initial probability distribution
          is used.

        OUTPUT:

        A Boolean.

        EXAMPLES::

            sage: from sage.stats.hmm.all import GeneralizedHiddenMarkovModel
            sage: T0 = Matrix([[1,2],[2,1]])
            sage: T1 = Matrix([[1,1],[0,2]])
            sage: G = GeneralizedHiddenMarkovModel([T0/5, T1/5])
            sage: G.initials_equivalent([1,0], [0,1])
            True
            sage: G.initials_equivalent([1,0])
            True
            sage: G2 = GeneralizedHiddenMarkovModel([T0/8, T0/8, T1/8])
            sage: G2.initials_equivalent([1/3,2/3], [1,0])
            True
            sage: T2 = Matrix([[1,1],[0,3]])
            sage: T3 = Matrix([[3,0],[1,1]])
            sage: GeneralizedHiddenMarkovModel([T2/5,T3/5]).initials_equivalent([1/3,2/3], [1,0])
            False

        ALGORITHM: Adapted from [Bal1993]_.
        """
        if init2 is None:
            init2 = self._initial_distribution
        return _initials_equivalent(self.Ts, init1, init2)

    def is_same_process_as(self, other):
        """
        Return ``True`` is ``other`` defines the same process as ``self``.
        This is true if and only if ``self.likelihood(w) == other.likelihood(w)``
        for all possible output words ``w``.

        INPUT:

        - ``other`` -- An instance of GeneralizedHiddenMarkovModel

        EXAMPLES::

            sage: from sage.stats.hmm.all import GeneralizedHiddenMarkovModel
            sage: T0 = Matrix([[1/2,0],[0,1/2]])
            sage: T1 = Matrix([[0,1/2],[1/2,0]])
            sage: G = GeneralizedHiddenMarkovModel([T0, T1])
            sage: T3 = Matrix([[1/2]])
            sage: G1 = GeneralizedHiddenMarkovModel([T3,T3]) # fair coin
            sage: G.is_same_process_as(G1)
            True
            sage: G2 = GeneralizedHiddenMarkovModel([4/3*T3,2/3*T3]) # unfair coin
            sage: G.is_same_process_as(G2)
            False
            sage: G2.is_same_process_as(G1)
            False

        ALGORITHM:

        Adapted from the algorithm in [Bal1993]_: construct the disjoint sum GHMM,
        then see whether two different initial probabilities are equivalent.
        """
        n1 = self.n_states
        n2 = other.n_states
        Ts1 = self.Ts
        Ts2 = other.Ts
        if len(Ts1) != len(Ts2):
            return False
        z0 = ZZ.zero()
        return _initials_equivalent([block_matrix(2,2,[T1, z0, z0, T2]) for T1,T2 in zip(Ts1, Ts2)],
                                   list(self.initial_probabilities()) + [z0]*n2,
                                   [z0]*n1 + list(other.initial_probabilities()))


    # Functions for minimization/canonicalization.  See [Upp1997] for an
    # explanation of the algorithm.

    def _future_distribution(self, word):
        """
        Return ``self.initial_distribution * self.T(word)``, scaled to sum to 1.

        This is the vector of pseudo-probabilities of lying in each state after the GHMM
        emits ``word``.

        INPUT:

        - ``w`` -- an iterable (e.g., list) of outputs.

        EXAMPLES:

            sage: G = hmm.GeneralizedHiddenMarkovModel([matrix.circulant([1/2,1/4,1/4])], [1,0,0])
            sage: G._future_distribution([0])
            (1/2, 1/4, 1/4)
            sage: G._future_distribution([0,0])
            (3/8, 5/16, 5/16)
        """
        d = self._initial_distribution
        for i in word:
            d = d * self.Ts[i]
        if sum(d) > 0:
            return d / sum(d)
        return d

        
    def _history_sufficient_wordlist(self):
        """
        Return a "sufficient wordlist" for ``self`` together with a (row) basis matrix
        for the set of "mixed states" these words span.

        This basis spans the space of reachable mixed states, that is, the set of vectors
        ``self._future_distribution(w)`` as ``w`` ranges over all possible outputs.

        EXAMPLES:

            sage: T0 = Matrix(3,3,[4,4,0,0,6,0,2,5,0])/12
            sage: T1 = Matrix(3,3,[0,0,4,6,0,0,3,0,2])/12
            sage: H = hmm.GeneralizedHiddenMarkovModel([T0, T1], [0,1,0])
            sage: H._history_sufficient_wordlist()
            (
                               [0 1 0]
                               [1 0 0]
            [[], [1], [1, 1]], [0 0 1]
            )
            sage: H2 = hmm.GeneralizedHiddenMarkovModel([T0, T1], [1,0,0])
            sage: H2._history_sufficient_wordlist()
            (
                            [  1   0   0]
                            [1/2 1/2   0]
            [[], [0], [1]], [  0   0   1]
            )
        """
        words = [[]]
        Q = [[i] for i in range(self.n_outputs)]
        # find a suitable ring, holds initial distribution plus transitions...
        R = self._transition_matrix.stack(self._initial_distribution).base_ring()
        distributions = matrix(R, [self.initial_probabilities()])
        while len(Q) > 0:
            w = Q.pop(0)
            v = self._future_distribution(w)
            if v not in distributions.row_space():
                words.append(w)
                distributions = distributions.stack(v, subdivide=False)
                Q.extend(w+[i] for i in range(self.n_outputs))
        return words, distributions

    def _past_distribution(self, word):
        """
        Return ``self.T(word) * vector([1,1...])``, scaled to sum to 1.

        This is a sort of backwards-looking version of ``self._future_distribution(word)``.

        INPUT:

        - ``w`` -- an iterable (e.g., list) of outputs.

        EXAMPLES:

            sage: T0 = Matrix([[1/3, 0], [0, 2/3]]); T1 = Matrix([[0, 2/3], [1/3, 0]])
            sage: G = hmm.GeneralizedHiddenMarkovModel([T0, T1], [1/4,3/4])
            sage: G._past_distribution([])
            (1/2, 1/2)
            sage: G._past_distribution([1])
            (2/3, 1/3)
            sage: G._past_distribution([0,0])
            (1/5, 4/5)
        """
        d = vector([1]*self.n_states)
        for i in reversed(word):
            d = self.Ts[i] * d
        if sum(d) > 0:
            return d / sum(d)
        return d

    def _future_sufficient_wordlist(self):
        """
        Return a "future sufficient wordlist" for ``self`` together with a row
        basis matrix for the space the future distributions of these words span.

        This basis spans the set of vectors ``self._future_distribution(w)``
        as ``w`` ranges over all possible outputs.

        EXAMPLES:

            sage: T0 = Matrix(3,3,[4,4,0,0,6,0,2,5,0])/12
            sage: T1 = Matrix(3,3,[0,0,4,6,0,0,3,0,2])/12
            sage: H = hmm.GeneralizedHiddenMarkovModel([T0, T1], [0,1,0])
            sage: H._future_sufficient_wordlist()
            (
                       [ 1/3  1/3  1/3]
            [[], [0]], [8/21  2/7  1/3]
            )

        """
        words = [[]]
        Q = [[i] for i in range(self.n_outputs)]
        R = self._transition_matrix.base_ring()
        distributions = matrix(R, [1]*self.n_states)/self.n_states
        while len(Q) > 0:
            w = Q.pop(0)
            v = self._past_distribution(w)
            if v not in distributions.row_space():
                words.append(w)
                distributions = distributions.stack(v, subdivide=False)
                Q.extend([i]+w for i in range(self.n_outputs))
        return words, distributions


    def canonical_form(self):
        """Return Upper's canonical form for this GHMM, as defined in [Upp1997]_.

        This canonical form is a GHMM that has the smallest number of
        states possible for any GHMM that is equivalent to ``self``.

        EXAMPLES::

            sage: T0 = Matrix([[0,1,0],[0,0,0],[0,1,0]])/2
            sage: T1 = Matrix([[1,0,0],[0,0,2],[1,0,0]])/2
            sage: gmp1 = hmm.GeneralizedHiddenMarkovModel([T0, T1], [1/3,1/3,1/3])
            sage: gmp1
            Generalized Hidden Markov Model with 3 states and 2 output symbols
            sage: gmp1.canonical_form()
            Generalized Hidden Markov Model with 2 states and 2 output symbols
            sage: M = Matrix.circulant([1/2, 1/4, 1/6, 1/12])
            sage: G = hmm.GeneralizedHiddenMarkovModel([M/2, M/3, M/6]); G
            Generalized Hidden Markov Model with 4 states and 3 output symbols
            sage: G.canonical_form()
            Generalized Hidden Markov Model with 1 state and 3 output symbols
            sage: G.canonical_form().Ts
            [[1/2], [1/3], [1/6]]
            sage: Ts = [Matrix([[1/2,0],[0,1/3]]), Matrix([[1/2, 0], [0, 2/3]])]
            sage: hmm.GeneralizedHiddenMarkovModel(Ts, [1,0]).canonical_form().Ts
            [[1/2], [1/2]]
            sage: hmm.GeneralizedHiddenMarkovModel(Ts, [0,1]).canonical_form().Ts
            [[1/3], [2/3]]
            sage: hmm.GeneralizedHiddenMarkovModel(Ts, [1/2,1/2]).canonical_form()
            Generalized Hidden Markov Model with 2 states and 2 output symbols

        ALGORITHM: From [Upp1997]_.
        """
        # Compute the H and F matrices
        _, H1 = self._history_sufficient_wordlist()
        _, F1 = self._future_sufficient_wordlist()
        F1t = F1.transpose()
        M = H1 * F1t
        E = M.echelon_form()
        ok_rows = [i for i,er in enumerate(E.rows()) if not er.is_zero()]
        E = M.transpose().echelon_form()
        ok_cols = [i for i,er in enumerate(E.rows()) if not er.is_zero()]
        H = H1.matrix_from_rows(ok_rows)
        F = F1t.matrix_from_columns(ok_cols)
        # Hi is a pseudo-inverse for H
        Hi = F * (H * F).inverse()
        # Now conjugate the Ts and initial distribution
        Ts = [H * T * Hi for T in self.Ts]
        init = self._initial_distribution * Hi
        return GeneralizedHiddenMarkovModel(Ts, init)

    def sample(self, n_samples):
        r"""
        Generate ``n_samples`` outputs from the process represented by ``self``.

        The outputs are yielded one by one.  If ``w = list(self.samples(n)``,
        then ``w`` is generated with probability ``self.likelihood(w)``.

        EXAMPLES::

            sage: T0 = Matrix([[0,1/2],[1,0]])
            sage: T1 = Matrix([[0,1/2],[0,0]])
            sage: G = hmm.GeneralizedHiddenMarkovModel([T0, T1], [0,1])
            sage: set_random_seed(5)
            sage: list(G.sample(15))
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]

        We can collect statistics on the outputs. For ``G``, every odd output is 0
        and the rest are equally likely to be 0 or 1.  So if we sum up outputs of
        length 8, we should see sums of 0,1,2,3, or 4, with probility 1/16, 4/16,
        6/16, 4/16, and 1/16 respectively.  With 128 trials the expected number
        of occurrences of each sum is 8, 32, 48, 32, and 8, respectively.  The
        result we see below is typical, with a p-value of 0.74 for the
        chi-squared test.

        ::

            sage: set_random_seed(101)
            sage: d = [0]*5
            sage: for _ in range(128): d[sum(G.sample(8))] += 1
            sage: d
            [6, 37, 44, 34, 7]

        Systems with negative transition pseudo-probabilities are correctly handled::

            sage: S0 = Matrix([[1, -1/2], [2, -1]])
            sage: S1 = Matrix([[1/2,  0], [0,  0]])
            sage: GMP = hmm.GeneralizedHiddenMarkovModel([S0, S1], [4/3, -1/3])
            sage: set_random_seed(7)
            sage: list(GMP.sample(20))
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1]
        """
        # Method: We know Prob(w) = i * T(w) * one, with i the initial distribution
        # So after emitting w, the prob of output o is
        # Prob(o|w) = Prob(w+[o])/Prob(w) = (i*T(w)*T(a)*one) / (i*T(w)*one)
        # We will keep track of i*T(w), and emit siccessive outputs o with the right
        # probability, then set w = w+[o] and emit the next symbol.
        #
        # Note we cannot do the same as for HMMs, just evolving a single state,
        # because of the possibility of negative "probabilities" in the transition
        # matrices.
        one = vector(self.Ts[0].base_ring(), [1] * self.n_states)
        v = self._initial_distribution # v maintains i * T(w)
        Va = [Ta * one for Ta in self.Ts]
        for _ in range(n_samples):
            p = 0
            n = sum(v) # = v * 1
            r = n * random()
            for i, va in enumerate(Va):
                p += v * va
                if p >= r:
                    yield(i)
                    v = v * self.Ts[i] # w = w + [i]
                    break
            else:
                raise Error('should not happen')

    def graph(self):
        r"""
        Create a labelled directed graph from the Generalized Hidden Markov Model.

        Each vertex corresponds to a state, and each edge is labelled
        with a list giving the pseudo-probability of each output when
        taking the transition from its source to its destination.

        EXAMPLES::

            sage: T0 = Matrix([[0,1/2],[1,0]])
            sage: T1 = Matrix([[0,1/2],[0,0]])
            sage: g = hmm.GeneralizedHiddenMarkovModel([T0, T1], [0,1]).graph()
            sage: g.edges()
            [(0, 1, [1/2, 1/2]), (1, 0, [1, 0])]

        The type of graph output from this function can be used to construct a GHMM::

            sage: hmm.GeneralizedHiddenMarkovModel(g).Ts == [T0, T1]
            True
        """
        G = DiGraph(loops=True)
        for i in range(self.n_states):
            for j in range(self.n_states):
                label = [T[i,j] for T in self.Ts]
                if not all(x.is_zero() for x in label):
                    G.add_edge(i,j,label)
        return G

def unpickle_ghmm_v0(Ts, initial_probs):
    """
    Internal function used for unpickling.  Just calls the GHMM constructor.

    TESTS::

        sage: T0 = Matrix([[1/3,  0], [1/6,1/6]])
        sage: T1 = Matrix([[  0,1/3], [1/3,  0]])
        sage: T2 = Matrix([[1/6,1/6], [  0,1/3]])
        sage: G = hmm.GeneralizedHiddenMarkovModel([T0, T1, T2])
        sage: G1 = sage.stats.hmm.ghmm.unpickle_ghmm_v0(G.Ts, G.initial_probabilities())
        sage: G1
        Generalized Hidden Markov Model with 2 states and 3 output symbols
    """
    return GeneralizedHiddenMarkovModel(Ts, initial_probs)
