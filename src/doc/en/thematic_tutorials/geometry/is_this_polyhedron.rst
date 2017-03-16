.. -*- coding: utf-8 -*-

.. linkall

.. _is_this_polyhedron:

==============================================================
Is this polyhedron ... the one your looking for?
==============================================================

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>

Once a polyhedron object is constructed, it is often useful to check if it has certain 
properties.

Here is a list of properties that Sage can check.

.. note::

    The following list may note be complete due to recent additions of features. You can 
    check by making a :code:`tab` completion after typing :code:`is_` to see which methods
    are available.

Combinatorial isomorphism
==============================================================

Two polyhedra are *combinatorially isomorphic* if their face lattices are isomorphic.

The verification of this is done using the vertices to facets adjacency oriented graph.

::

    sage: P1 = Polyhedron(vertices = [[1, 0], [0, 1]], rays = [[1, 1]])
    sage: P7 = Polyhedron(vertices = [[3, 0], [4, 1]], rays = [[-1, 1]])
    sage: P1_and_P7 = P1 & P7
    sage: Square = Polyhedron(vertices = [[1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1]])
    sage: Square.is_combinatorially_isomorphic(P1_and_P7)
    True

.. end of output

Compactness or is it a polytope
==============================================================

Emptyness
==============================================================

Full-dimension
==============================================================

Lattice polyhedron
==============================================================

Inscribed on a sphere
==============================================================

Minkowski summand
==============================================================

Neighborlyness
==============================================================

Reflexiveness
==============================================================

Simplicity
==============================================================

Simpliciality
==============================================================

Is it the simplex?
==============================================================

Is it the whole space?
==============================================================

 'is_universe'
