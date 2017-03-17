.. -*- coding: utf-8 -*-

.. linkall

.. _is_this_polyhedron:

==============================================================
Is this polyhedron ...
==============================================================

.. MODULEAUTHOR:: Jean-Philippe Labbé <labbe@math.fu-berlin.de>

Once a polyhedron object is constructed, it is often useful to check if it has certain 
properties.

Here is a list of properties that Sage can check.

.. note::

    The following list may note be complete due to recent additions of features. You can 
    check by making a :code:`tab` completion after typing :code:`is_` to see which methods
    are available.

... combinatorially isomorphic to this one?
==============================================================

Two polyhedra are *combinatorially isomorphic* if their face lattices are isomorphic.


Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_combinatorially_isomorphic`

... compact/a polytope?
==============================================================

A compact polyhedron is also called a *polytope*.

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_compact`

... empty?
==============================================================

This function although sounding trivial is very important!

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_empty`

... full-dimensional?
==============================================================

A polyhedron is full dimensional when it does not have any equations in its
:math:`H`-representation.

Reference manual: :meth:`sage.geometry.polyhedron.base.Polyhedron_base.is_full_dimensional`

... a lattice polyhedron?
==============================================================

BUG HERE!



... inscribed on a sphere?
==============================================================


... a Minkowski sum of this other one?
==============================================================

... neighborly?
==============================================================

... reflexive?
==============================================================

... simple? (Aren't they all?)
==============================================================

... simplicial?
==============================================================

... the simplex?
==============================================================

... the whole space?
==============================================================

 'is_universe'
