r"""
Examples of Vertex Algebras

The following examples are implemented:

- :mod:`Abelian vertex algebra<.abelian_vertex_algebra>`
- :mod:`Affine vertex algebra<.affine_vertex_algebra>`
- :mod:`Free Bosons<.free_bosons_vertex_algebra>`
- :mod:`Free Fermions<.free_fermions_vertex_algebra>`
- :mod:`N=2 super vertex algebra<.n2_vertex_algebra>`
- :mod:`Neveu-Schwarz super vertex algebra<.neveu_schwarz_vertex_algebra>`
- :mod:`Virasoro vertex algebra<.virasoro_vertex_algebra>`
- :mod:`Weyl vertex algebra<.weyl_vertex_algebra>`

AUTHORS:

- Reimundo Heluani (2020-06-15): Initial implementation.
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

from .abelian_vertex_algebra import AbelianVertexAlgebra as Abelian
from .affine_vertex_algebra import AffineVertexAlgebra as Affine
from .free_bosons_vertex_algebra import FreeBosonsVertexAlgebra as FreeBosons
from .free_fermions_vertex_algebra import FreeFermionsVertexAlgebra as FreeFermions
from .n2_vertex_algebra import N2VertexAlgebra as N2
from .neveu_schwarz_vertex_algebra import NeveuSchwarzVertexAlgebra as NeveuSchwarz
from .virasoro_vertex_algebra import VirasoroVertexAlgebra as Virasoro
from .weyl_vertex_algebra import WeylVertexAlgebra as Weyl

assert Abelian
assert Affine
assert FreeBosons
assert FreeFermions
assert N2
assert NeveuSchwarz
assert Virasoro
assert Weyl
