r"""
Stochastic Matrices
"""
#*****************************************************************************
#  Copyright (C) 2014 Anne Schilling <anne at math.ucdavis.edu>
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

    def __init__(self,  M, side='row'):
        """
        EXAMPLES::

            sage: C = sage.categories.examples.crystals.HighestWeightCrystalOfTypeA(n=4)
            sage: C == Crystals().example(n=4)
            True
        """
        assert M.is_square()

        ElementWrapper.__init__(self,  M.parent(), M)


#    def _repr_(self):
#        """
#        EXAMPLES::
#
#            ERROR here
#            
#        """
#        return "A stochastic matrix: %s"%(self,)




