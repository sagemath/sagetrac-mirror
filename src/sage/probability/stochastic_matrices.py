r"""
Stochastic Matrices
"""
#*****************************************************************************
#  Copyright (C) 2014 Arvind Ayyer <arvind at math.iisc.ernet.in>
#                 and Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

#from sage.structure.element import Element
from sage.structure.element_wrapper import ElementWrapper
#from sage.matrix.matrix_rational_dense import Matrix_rational_dense

class StochasticMatrix(ElementWrapper):
    r"""

    EXAMPLES::

     """

    def __init__(self,  M, side='row',normalized = False):
        r"""
        Converts a matrix of nonnegative real numbers to the class of (row/ column) stochastic matrices.
        A matrix is said to be row (resp. column) stochastic if the sum of the
        entries of each row (resp. column) is equal to 1.

        INPUT:

        - ``side`` -- (default: 'row') if set to ``row`` checks whether row sums
          satisfy the condition; if set to ``column`` checks whether column sums
          satisfy the condition.

        - ``normalized`` -- if set to ``True``, checks that
          the sums are equal to 1. When set to ``False`` (default), checks that
          the row (resp. column) sums are all equal to some constant
          possibly different from 1.

        EXAMPLES::

            sage: M = Matrix([[1/2,1/2],[1,0]])
            sage: StochasticMatrix(M)
            [1/2 1/2]
            [  1   0]

        TESTS::

        sage: M = matrix([[1, 1],[2,2]])
        sage: StochasticMatrix(M)
        Traceback (most recent call last):
        ...
        ValueError: Input matrix is not stochastic
        sage: StochasticMatrix(M,side='column',normalized = True)
        Traceback (most recent call last):
        ...
        ValueError: Input matrix is not stochastic
        sage: StochasticMatrix(M,side='column')
        [1 1]
        [2 2]
        """
        if not M.is_stochastic(side = side,normalized = normalized):
            raise ValueError("Input matrix is not stochastic")
        self.side = side
        ElementWrapper.__init__(self,  M.parent(), M)

    def __getattr__(self,name):
        try:
            return getattr(self.value,name)
        except AttributeError:
            return super(StochasticMatrix, self).__getattr__(name)

    def stationary_distribution(self, irreducible = True):
        r"""
        Returns the stationary distribution of the stochastic matrix

        EXAMPLES::

            sage: M = Matrix([[1/2,1/2],[1,0]])
            sage: StochasticMatrix(M)
            [1/2 1/2]
            [  1   0]

        TESTS::

        """

        if self.side == 'row':
            sums = map(sum, self.rows())[0]
            K = (self.value-sums).kernel()
        else:
            sums = map(sum, self.columns())[0]
            K = (self.value-sums).transpose().kernel()

        if not irreducible:
            return K
        else:
            if K.dimension() > 1:
                raise ValueError("The stochastic matrix is not irreducible")
            else:
                V = K.basis()[0]
                return V/sum(V)

    def to_digraph(self):
        r"""
        Returns the graph of the stochastic matrix. The Markov chain 
        of the stochastic matrix corresponds
        to a random walk on this directed graph.

        EXAMPLES::

            sage: M = Matrix([[1/2,1/2],[1,0]])
            sage: StochasticMatrix(M)
            [1/2 1/2]
            [  1   0]

        TESTS::

        """
        return DiGraph(self, format = 'weighted_adjacency_matrix', loops = 'true')


