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
        """
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
        sage: StochasticMatrix(M,side='column')
        [1 1]
        [2 2]

        """
        if not M.is_stochastic(side = side,normalized = normalized):
            raise ValueError("Input matrix is not stochastic")

        ElementWrapper.__init__(self,  M.parent(), M)


#    def _repr_(self):
#        """
#        EXAMPLES::
#
#            ERROR here
#            
#        """
#        return "A stochastic matrix: %s"%(self,)




