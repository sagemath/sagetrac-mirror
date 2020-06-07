r"""
Virasoro Vertex Algebra

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation.
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


from .universal_enveloping_vertex_algebra import \
                                                UniversalEnvelopingVertexAlgebra
from sage.sets.family import Family

class VirasoroVertexAlgebra(UniversalEnvelopingVertexAlgebra):
    def __init__(self, R, c, arg0 = None, names="L"):
        r"""
        The universal Virasoro vertex algebra

        INPUT:

        - ``R`` a ring, the base ring of this vertex algebra. Undefined behaviour if
          this is not a Field of characteristic zero.

        - ``c`` The central charge of this vertex algebra if ``arg0`` is not
          specified.

        - ``arg0`` a positive integer or ``None`` (Default: ``None``). If specified
          a positive integer `q` then the parameter `c` has to be a positive integer
          `p` coprime with `q` and this vertex algebra returns the irreducible
          quotient of the Virasoro vertex algebra at central charge

        .. MATH::

            c = 1 - 6 \frac{(p-q)^2}{pq}


        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); V
            The Virasoro vertex algebra at central charge 1/2

        """
        from sage.algebras.lie_conformal_algebras.virasoro_lie_conformal_algebra import\
                                                    VirasoroLieConformalAlgebra
        ML = VirasoroLieConformalAlgebra(R)
        if arg0 is not None:
            c = 1  - 6*(c-arg0)**2/(c*arg0)
        cp = Family({ML.gen(1):c})
        super(VirasoroVertexAlgebra,self).__init__(R, ML,
                 central_parameters=cp, names = names)
        self._c = c

    def _repr_(self):
        return "The Virasoro vertex algebra at central charge {}".format(self.central_charge())



