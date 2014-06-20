r"""
Lattice Polytopes in the Graded Ring Database

The Graded Ring Database at http://grdb.lboro.ac.uk
contains databases of Fano polytopes.

Fano polytopes are the polytopes that contain the origin in their
strict interior such that all their vertices are primitive. Given a
lattice polytope `P` with the origin in its strict interior,
one can construct the spanning fan of `P` and hence a toric
variety `X_P`. The toric variety `X_P` is Fano if and only if `P`
is a Fano polytope.

This module contains functions that will return lattice polytope
objects from a selection of GRDB databases,
which can then be manipulated by the tools in
:mod:`lattice_polytope<sage.geometry.lattice_polytope>`.

The databases are located in the grdb_polytopes package, which
will need to be installed in order to access these functions
successfully.

For documentation on the individual functions see:

    * :func:`lattice_polytopes.CanonicalFano<sage.geometry.lattice_polytopes_backend.CanonicalFano>`
    * :func:`lattice_polytopes.ReflexiveFano<sage.geometry.lattice_polytopes_backend.ReflexiveFano>`
    * :func:`lattice_polytopes.SmallPolygon<sage.geometry.lattice_polytopes_backend.SmallPolygon>`
    * :func:`lattice_polytopes.SmoothFano<sage.geometry.lattice_polytopes_backend.SmoothFano>`
    * :func:`lattice_polytopes.TerminalFano<sage.geometry.lattice_polytopes_backend.TerminalFano>`
    * :func:`lattice_polytopes.lReflexive<sage.geometry.lattice_polytopes_backend.lReflexive>`
    * :func:`lattice_polytopes.LDP<sage.geometry.lattice_polytopes_backend.LDP>`

EXAMPLES::

    sage: C = lattice_polytopes.CanonicalFano(3,34); C #optional - GRDB_polytopes
    A lattice polytope: 3-dimensional, 14 vertices.

For more examples see the individual functions above.

AUTHORS:

- Samuel Gonshaw (2012-08-06): initial version
- Tom Coates (2014-02-24): minor changes to docstring
"""
from lattice_polytopes_backend import (CanonicalFano, ReflexiveFano,
                                        SmallPolygon, SmoothFano,
                                        TerminalFano, lReflexive, LDP)
