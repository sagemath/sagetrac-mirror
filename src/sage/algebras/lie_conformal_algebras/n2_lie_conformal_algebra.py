"""
N=2 Super Lie Conformal Algebra.

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

class N2LieConformalAlgebra(GradedLieConformalAlgebra):
    def __init__(self,R, names=('L','J','G1','G2')):
        """
        The N=2 super Lie conformal algebra.

        EXAMPLES::

            sage: R = N2LieConformalAlgebra(QQ); R
            The N=2 super Lie conformal algebra over Rational Field
            sage: R.inject_variables()
            Defining L, J, G1, G2, C
            sage: G1.bracket(G2)
            {0: L + 1/2*TJ, 1: J, 2: 1/3*C}
            sage: G2.bracket(G1)
            {0: L - 1/2*TJ, 1: -J, 2: 1/3*C}
            sage: G1.degree()
            3/2
            sage: J.degree()
            1
        """
        n2dict =\
        {('L','L'):{0:{('L',1):1}, 1:{('L',0): 2}, 
        3:{('C', 0):R(2).inverse_of_unit()}},
        ('L','G1'):{0:{('G1',1):1}, 1:{('G1',0):3*R(2).\
        inverse_of_unit()}}, 
        ('L','G2'):{0:{('G2',1):1}, 1:{('G2',0):3*R(2).\
        inverse_of_unit()}}, 
        ('G1','G2'): {0:{('L',0):1,('J',1):R(2).inverse_of_unit()},
                   1:{('J',0):1}, 2:{('C',0):R(3).inverse_of_unit()}},
        ('L','J'): {0:{('J',1):1},1:{('J',0):1}},
        ('J','J'): {1:{('C',0):R(3).inverse_of_unit()}},
        ('J','G1'): {0:{('G1',0):1}},
        ('J','G2'): {0:{('G2',0):-1}}}
        from sage.rings.rational_field import QQ
        weights = (2,1,QQ(3/2),QQ(3/2))
        parity = (0,0,1,1)
        GradedLieConformalAlgebra.__init__(self,R,n2dict,
                                           names=('L', 'J','G1','G2'),
                                           central_elements=('C',), 
                                           weights=weights, parity=parity)

    def _repr_(self):
        return "The N=2 super Lie conformal algebra over {}".\
                format(self.base_ring())


