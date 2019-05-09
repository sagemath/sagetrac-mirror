.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Introduction to topaz
---------------------

This tutorial tries to give the user a first idea about the features of
the ``topaz`` application of ``polymake``. We take a look at a variety
of small examples.

First, we have to make ``topaz`` the current application. For this you
can either start ``polymake`` with the option ``-A topaz``,

::

   polymake -A topaz

or, if you’ve already started ``polymake``, type


::

    polymake> application 'topaz';

in the ``polymake`` shell.

Simplicial complexes
~~~~~~~~~~~~~~~~~~~~

The most important object of the ``topaz`` application is the simplicial
complex. There are several ways of obtaining one.

From faces
^^^^^^^^^^

For example, you can specify some faces of the complex. You can pass
them as an ``Array< Set<Int> >``, or ``Array< Array<Int> >``:


::

    polymake> # $s = new SimplicialComplex(INPUT_FACES=>[new Set(0), new Set(0,1), new Set(1,2,3)]);
    ........> $s = new SimplicialComplex(INPUT_FACES=>[[0],[0,1],[1,2,3]]);

|{{ :tutorial:small_complex.png?400|}}| As you can see, redundancies are
allowed – ``[0]`` is not a facet of the complex, and thus not necessary
for encoding ``$s``. You can compute the inclusion maximal faces like
this:

.. |{{ :tutorial:small_complex.png?400|}}| image:: attachment:small_complex.png


::

    polymake> print $s->FACETS;
    {0 1}
    {1 2 3}
    





You can also pass the ``FACETS`` to the constructor, but be aware that
in that case the vertices must be numbered increasingly starting with
``0`` and redundancies are prohibited.

Take a look at your complex using


::

    polymake> $s->VISUAL;

For more information on visualizing simplicial complex, see the section
below.

|{{:tutorial:face_lattice.png?200 \|}}|\ ``polymake`` can compute the
Hasse diagram of a simplicial complex (watch out, this gets really large
for large complexes!). To print all the faces of the complex together
with their rank in the face lattice, do this:

.. |{{:tutorial:face_lattice.png?200 \|}}| image:: attachment:face_lattice.png


::

    polymake> print $s->HASSE_DIAGRAM->DECORATION;
    ({-1} 4)
    ({0 1} 2)
    ({1 2 3} 3)
    ({0} 1)
    ({1} 1)
    ({1 2} 2)
    ({1 3} 2)
    ({2 3} 2)
    ({} 0)
    ({2} 1)
    ({3} 1)
    





The first entry of each pair denotes the face, the second is the rank.
The ``{-1}``-node is a dummy representing the whole complex. the
``{}``-node is the empty face. If you want to look at a pretty graph
representation, try the visualization:


::

    polymake> $s->VISUAL_FACE_LATTICE;

Using clients
^^^^^^^^^^^^^

There are several clients that construct common simplicial complexes
(for a comprehensive list, see the `topaz
documentation <https://polymake.org/release_docs/latest/topaz.html>`__).
An example is the torus client:


::

    polymake> $t = torus();

Of course, ``polymake`` can compute the reduced integer homology groups
of a simplicial complex, so we can convice ourselves this is a torus:


::

    polymake> print $t->MANIFOLD;
    1





::

    polymake> print $t->HOMOLOGY;
    ({} 0)
    ({} 2)
    ({} 1)
    





The ``i``-th line represents the `i`-th homology module. The curly
braces contain torsion coefficients with multiplicity, the second pair
entry denotes the Betti number. The empty curly braces indicate that
``$t`` is torsion-free. You can see a non-empty torsion group here
(using the ``rows_numbered`` client for a pretty print with the
corresponding dimensions):


::

    polymake> print rows_numbered( real_projective_plane()->HOMOLOGY );
    0:{} 0
    1:{(2 1)} 0
    2:{} 0
    





As expected, the first homology group has torsion coefficient ``2`` with
multiplicity ``1`` and all Betti numbers are zero.

As boundary complex
^^^^^^^^^^^^^^^^^^^

If your complex is a pseudo-manifold, you can obtain a new complex from
its boundary. For example, this produces a triangulation of the
`2`-sphere:


::

    polymake> $bs = simplex(3)->BOUNDARY;
    ........> print $bs->SPHERE;
    1
    





Triangulating polytopes
^^^^^^^^^^^^^^^^^^^^^^^

The triangulation of a polytope is a simplicial complex, too. The
``TRIANGULATION`` gets stored in a property of the polytope. We use the
``cube`` client from the ``polytope`` application to demonstrate:


::

    polymake> $c = polytope::cube(3);
    ........> $tc = $c->TRIANGULATION;
    ........> print $tc->FACETS;
    {0 1 2 4}
    {1 2 3 4}
    {1 3 4 5}
    {2 3 4 6}
    {3 4 5 6}
    {3 5 6 7}
    





Geometric realizations
~~~~~~~~~~~~~~~~~~~~~~

The ``topaz`` application is primarily designed to deal with abstract
simplicial complexes that do not come with coordinates for an embedding
in euclidean space. There is a special object subtype named
``GeometricSimplicialComplex`` that has extra properties for dealing
with coodinates.

You can pass the coordinates to the constructor. Take care to choose an
embedding without crossings!


::

    polymake> $s = new GeometricSimplicialComplex(INPUT_FACES=>[[0],[0,1],[1,2,3]], COORDINATES=>[[1,0],[1,1],[0,2],[2,2]]);

Some clients produce complexes with geometric realization…


::

    polymake> $b = ball(3);
    ........> # print a dense representation of the sparse matrix
    ........> print dense( $b->COORDINATES );
    0 0 0
    1 0 0
    0 1 0
    0 0 1
    





…some others provide the option ``geometric_realization`` so you can
decide whether to invest the extra computing time.


::

    polymake> $bs = barycentric_subdivision($b,geometric_realization=>1);

Again, see the `topaz
documentation <https://polymake.org/release_docs/latest/topaz.html>`__
for a comprehensive list.

Visualization
~~~~~~~~~~~~~

Visualization of simplicial complexes uses the ``VISUAL`` property.
Check out


::

    polymake> help 'objects/SimplicialComplex/methods/Visualization/VISUAL';

|{{ :tutorial:ball_triang.png?300|}}| for a list of available options
and this `tutorial <visual_tutorial>`__ for a general intro to
visualization in polymake.

If your complex is of dimension three or lower, you can visualize a
geometric realization together with the ``GRAPH`` of the complex using
the ``VISUAL`` property. Note that if your complex is not a
``GeometricSimplicialComplex``, ``polymake`` will use the spring
embedder to find an embedding of the graph of the complex, which is not
guaranteed to result in an intersection-free visualization.

.. |{{ :tutorial:ball_triang.png?300|}}| image:: attachment:ball_triang.png


::

    polymake> $bs->VISUAL;

You should give the ``explode`` feature of jReality a try – it gives a
good (and pretty!) overview of the object. You can find it in the left
slot of the jReality interface.

|{{:tutorial:ball_triang_pink.png?250 \|}}| ``topaz`` may also visualize
distinguished subcomplexes or just sets of faces with different
decorations (colors, styles, etc.). For example, to highlight the fourth
facet of ``$bs`` in pink, do this:

.. |{{:tutorial:ball_triang_pink.png?250 \|}}| image:: attachment:ball_triang_pink.png


::

    polymake> $a = new Array<Set<Int>>(1); $a->[0] = $bs->FACETS->[4];
    ........> $bs->VISUAL->FACES($a, FacetColor => 'pink');

The same can be used for the visualization of the face lattice. As an
example, we have a look at a ``morse matching`` of the Klein bottle with
its associated critical faces. In order to see the arrowheads in the
picture clearly, you ought to use graphviz or svg to vizualize it.


::

    polymake> $k =  klein_bottle();
    ........> graphviz($k->VISUAL_FACE_LATTICE->MORSE_MATCHING->FACES($k->MORSE_MATCHING->CRITICAL_FACES));

|{{ :tutorial:kb_mm_faces.gif?400|}}| Here the matching of faces is
denoted by reversed red arrows and the critical faces are marked red.
Check that the graph remains acyclic.

For higher dimensional complexes that cannot be visualized in 3D, you
can still have a look at the graphs while ignoring any specified
coordinates by using ``VISUAL_GRAPH``, ``VISUAL_DUAL_GRAPH``, or
``VISUAL_MIXED_GRAPH``. An easy example:

.. |{{ :tutorial:kb_mm_faces.gif?400|}}| image:: attachment:kb_mm_faces.gif


::

    polymake> polytope::cube(3)->TRIANGULATION->VISUAL_MIXED_GRAPH;

shows the primal and dual graph of the polytope together with an edge
between a primal and a dual node iff the primal node represents a vertex
of the corresponding facet of the dual node.

.. figure:: attachment:cube_graph.png
   :alt: {{ :tutorial:cube_graph.png?600 \|}}

   {{ :tutorial:cube_graph.png?600 \|}}

Visualization of the ``HASSE_DIAGRAM`` is possible via
``VISUAL_FACE_LATTICE``. It renders the graph in a .pdf file. You can
even pipe the tikz code to whatever location using the ``tikz`` client:

::

   tikz($s->VISUAL_FACE_LATTICE, File=>"/path/to/file.tikz");
