"""
Neveu-Schwarz Super Lie Conformal Algebra.

AUTHORS:

- Reimundo Heluani (06-03-2020): Initial implementation
"""
#******************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .graded_lie_conformal_algebra import GradedLieConformalAlgebra

class NeveuSchwarzLieConformalAlgebra(GradedLieConformalAlgebra):

    def __init__(self,R,names=('L','G')):
        """
        The Neveu-Schwarz super Lie conformal algebra.

        EXAMPLES::

            sage: R = NeveuSchwarzLieConformalAlgebra(AA); R
            The Neveu-Schwarz super Lie conformal algebra over Algebraic Real Field
            sage: R.structure_coefficients()
            Finite family {('G', 'G'): ((0, L), (2, 1/3*C)),  ('G', 'L'): ((0, 1/2*TG), (1, 3/2*G)),  ('L', 'G'): ((0, TG), (1, 3/2*G)),  ('L', 'L'): ((0, TL), (1, 2*L), (3, 1/2*C))}
            sage: R.inject_variables()
            Defining L, G, C
            sage: G.nproduct(G,0)
            L
            sage: G.degree()
            3/2
        """
        nsdict =  {('L','L'):{0:{('L',1):1}, 1:{('L',0): 2}, 
                    3:{('C', 0):R(2).inverse_of_unit()}},
                   ('L','G'):{0:{('G',1):1}, 1:{('G',0):3*R(2).\
                    inverse_of_unit()}}, ('G','G'): {0:{('L',0):1},
                    2:{('C',0):R(3).inverse_of_unit()}}}
        from sage.rings.rational_field import QQ
        weights = (2,QQ(3/2))
        parity = (0,1)
        GradedLieConformalAlgebra.__init__(self,R,nsdict,names=names,
            central_elements=('C',), weights=weights, parity=parity)

    def _repr_(self):
        return "The Neveu-Schwarz super Lie conformal algebra over {}".\
                format(self.base_ring())


