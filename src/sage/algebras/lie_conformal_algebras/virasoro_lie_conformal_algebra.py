"""
Virasoro Lie Conformal Algebra.

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .graded_lie_conformal_algebra import GradedLieConformalAlgebra
class VirasoroLieConformalAlgebra(GradedLieConformalAlgebra):
    def __init__(self, R):
        """
        The Virasoro Lie Conformal algebra over `R`

        INPUT:

        - ``R``: a commutative Ring. Behaviour is undefined if `R` is 
          not a Field of characteristic zero. 

        EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ)
                sage: Vir.category()
                Category of finitely generated H-graded Lie conformal algebras with basis over Rational Field
                sage: Vir.gens()
                (L, C)
                sage: L = Vir.0
                sage: sorted(L.bracket(L).items())
                [(0, TL), (1, 2*L), (3, 1/2*C)]
        """
        virdict =  {('L','L'):{0:{('L',1):1}, 1:{('L',0): 2}, 
                    3:{('C', 0):R(2).inverse_of_unit()}}}
        GradedLieConformalAlgebra.__init__(self,R, virdict, 
            names = ('L',), central_elements = ('C',), weights = (2,))

    def _repr_(self):
        return "The Virasoro Lie conformal algebra over {}".format(
                                                            self.base_ring())


