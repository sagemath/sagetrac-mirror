.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Semantics of Cones and Polytopes
--------------------------------

General Remarks
~~~~~~~~~~~~~~~

The general semantics of a `big
object <https://polymake.org/doku.php/user_guide/lingo#big_object>`__ in
polymake is as follows: a list of properties describes an equivalence
class of mathematical objects. Often this equivalence class consists of
a single element, but this is not necessary.

As an example an object of class
`Polytope <https://polymake.org/release_docs/latest/polytope.html#polytope__Polytope__9>`__
defined by ``VERTICES`` gives such a single element equivalence class. A
typical example of a class with several elements is a polytope given
combinatorially, in terms of ``VERTICES_IN_FACETS``. An extreme case
would be a Polytope object defined by ``VOLUME`` only, defining the set
of polytopes of all possible dimensions which happen to have that
volume. While this is not very useful, a similar example would be a
Polytope object defined by ``F_VECTOR`` only. From this it makes sense
to derive, e.g., ``N_VERTICES`` or ``H_VECTOR``.

All big objects are immutable as mathematical objects. This means it is
possible to add more properties, but only consistent ones. Ideally,
these properties pre-exist (since they are logically derived from the
input description of the object), and the rules only make them explicit.
If a user asks for a property which cannot be derived, this property is
set to ``undef``. This occurs, e.g., if one asks for the ``VERTICES`` of
a combinatorially defined polytope.

To view the list properties that currently constitute your object, you
can use the ``properties`` method.


.. link

.. CODE-BLOCK:: perl

    polymake> $p = new Polytope(POINTS=>[[1,2],[1,3]]);




.. link

.. CODE-BLOCK:: perl

    polymake> $p->properties;
    name: p
    type: Polytope<Rational>
    
    POINTS
    1 2
    1 3
    
    
    CONE_AMBIENT_DIM
    2





::

   POINTS
   1 2
   1 3


   CONE_AMBIENT_DIM
   2

Objects of type ``Polytope``
----------------------------

Polytope theory is nice because this is where combinatorics meets metric
geometry. For mathematical software dealing with such objects it is
necessary to get the semantics straight. Below we describe some
pitfalls.

With coordinates: Geometry
~~~~~~~~~~~~~~~~~~~~~~~~~~

Being non-empty is recorded in the property ``FEASIBLE``. This is
``true`` if and only if the polytope is not empty.


.. link

.. CODE-BLOCK:: perl

    polymake> print cube(3)->FEASIBLE;
    true




A non-empty polytope in R^n is encoded as its homogenization in R^{n+1}.
Hence, any non-empty polytope has at least one facet (which may be the
far hyperplane [1,0,0,…,0]) and one vertex.

Without coordinates: Combinatorics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``VERTICES_IN_FACETS`` always describes the combinatorics of a bounded
polytope: this is any polytope which is projectively equivalent to the
polyhedron defined by ``VERTICES`` or ``POINTS`` or dually modulo its
``LINEALITY_SPACE``.

Each property must clearly specify if it depends on the geometry or only
on the combinatorics.

Special Cases
~~~~~~~~~~~~~

Most of what comes below is a consequence of the design decisions
explained above.

Empty polytopes
^^^^^^^^^^^^^^^

With the introduction of the ``Cone`` class and redefining ``Polytope``
as a derived class (in version 2.9.10) this raises the question of how
to deal with empty polytopes. This is a bit subtle as the cone over an
empty polytope does not have a canonical definition. Most text books
hence exclude this case. For them a polytope is never empty. There was a
time when this was also polymake’s point of view (until version 2.3).

However, this was changed for the reason that often people generate
systems of inequalities and then look at the feasible region. Most of
the time they obtain a polytope and proceed, but sometimes it fails, and
the region is empty. It is therefore necessary to give a definition of
the empty polytope (geometrically) which is consistent:

An empty polytope is recognized by ``FEASIBLE == false``. Such a
polytope is required to have ``VERTICES`` and ``FACETS`` empty.


.. link

.. CODE-BLOCK:: perl

    polymake> $e = new Polytope(POINTS=>[]);
    polymake> print $e->FEASIBLE;
    false







.. link

.. CODE-BLOCK:: perl

    polymake> print $e->FACETS;

This is totally different from having ``VERTICES`` or ``FACETS``
undefined (see above).


.. link

.. CODE-BLOCK:: perl

    polymake> $nc = new Polytope(VERTICES_IN_FACETS => cube(2)->VERTICES_IN_FACETS);

Zero-dimensional polytopes
^^^^^^^^^^^^^^^^^^^^^^^^^^

A zero-dimensional polytope is a single point. In our model it has one
vertex and one facet (the far hyperplane).


.. link

.. CODE-BLOCK:: perl

    polymake> $z = new Polytope(POINTS=>[[1,2,3]]);




.. link

.. CODE-BLOCK:: perl

    polymake> print $z->FACETS;
    1 0 0





``VERTICES_IN_FACETS`` is a 1-by-1 matrix with a zero entry. This means
that the single vertex does *not* lie on the single facet.


.. link

.. CODE-BLOCK:: perl

    polymake> print $z->VERTICES_IN_FACETS;
    {}





Such a polytope is both simple and simplicial, i.e. it is a simplex.


.. link

.. CODE-BLOCK:: perl

    polymake> print $z->SIMPLICIAL,",",$z->SIMPLE;
    true,true




Zero-dimensional fans
^^^^^^^^^^^^^^^^^^^^^

A zero-dimensional fan can e.g. be defined via


.. link

.. CODE-BLOCK:: perl

    polymake> $f = new fan::PolyhedralFan(RAYS=>[], MAXIMAL_CONES=>[[]]);

Summing Up
~~~~~~~~~~

For instance we have four possibilities which can occur for
``VERTICES``. The property

-  does not exist (it is not listed in ``properties``): This basically
   means that the property is not derived/calculated, yet.

-  exists and is set to ``undef``: Polymake is not able to derive this
   property with the given properties. The polytope may be empty or not.

-  exists and is empty: So the polytope is empty.

-  exists and is neither set to ``undef`` nor is empty: Our polytope is
   not empty and the property returns what you expect.
