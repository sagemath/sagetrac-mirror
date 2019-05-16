.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Objects, Properties and Rules
=============================

Objects
~~~~~~~

In polymake, there is two kinds of objects. A *Big Object* models a
complex mathematical concept, like a Polytope or a SimplicialComplex,
while a *small object* is an instance of one of the many data types
commonly used in computer science, like Integers, Matrices, Sets or
Maps. A big object consists of a collection of other objects (big or
small) describing it, called *properties*, and functions to compute more
properties from the ones already known, called *production rules*.

To find out the type of an object ``$c``, enter

::

   print $c->type->full_name;

To get a more detailed explanation of the ``polymake`` object model and
properties, check out the `scripting
guide <:user_guide:howto:scripting#most_important_interfaces>`__.

You can save polymake objects to disc, as explained
`here <:user_guide:tutorials:data>`__.

Properties
~~~~~~~~~~

Each (big) object has a list of properties of various types. When an
object is ‘born’ it comes with an initial list of properties, and all
other properties will be derived from those. Let’s look at example from
the ``polytope`` application. The following creates a 3-dimensional
cube:


.. link

.. CODE-BLOCK:: perl

    polymake> $c=cube(3);

To find out what the initial set of properties is, use the
``list_properties`` method. It returns an array of strings. The extra
code is just there to print this list nicely.


.. link

.. CODE-BLOCK:: perl

    polymake> print join(", ", $c->list_properties);
    CONE_AMBIENT_DIM, CONE_DIM, FACETS, AFFINE_HULL, VERTICES_IN_FACETS, BOUNDED




To see what a property contains, use the ``->`` syntax:


.. link

.. CODE-BLOCK:: perl

    polymake> print $c->FACETS;
    1 1 0 0
    1 -1 0 0
    1 0 1 0
    1 0 -1 0
    1 0 0 1
    1 0 0 -1





::

   1 1 0 0
   1 -1 0 0
   1 0 1 0
   1 0 -1 0
   1 0 0 1
   1 0 0 -1

You can also get the content of all properties using the ``properties``
method:


.. link

.. CODE-BLOCK:: perl

    polymake> $c->properties;
    name: c
    type: Polytope<Rational>
    description: cube of dimension 3
    
    
    CONE_AMBIENT_DIM
    4
    
    CONE_DIM
    4
    
    FACETS
    1 1 0 0
    1 -1 0 0
    1 0 1 0
    1 0 -1 0
    1 0 0 1
    1 0 0 -1
    
    
    AFFINE_HULL
    
    
    VERTICES_IN_FACETS
    {0 2 4 6}
    {1 3 5 7}
    {0 1 4 5}
    {2 3 6 7}
    {0 1 2 3}
    {4 5 6 7}
    
    
    BOUNDED
    true





::

   CONE_AMBIENT_DIM
   4

   CONE_DIM
   4

   FACETS
   1 1 0 0
   1 -1 0 0
   1 0 1 0
   1 0 -1 0
   1 0 0 1
   1 0 0 -1


   AFFINE_HULL


   VERTICES_IN_FACETS
   {0 2 4 6}
   {1 3 5 7}
   {0 1 4 5}
   {2 3 6 7}
   {0 1 2 3}
   {4 5 6 7}


   BOUNDED
   1

Production Rules
~~~~~~~~~~~~~~~~

The object is changed if we ask for a property which has not been
computed before.


.. link

.. CODE-BLOCK:: perl

    polymake> print $c->VERTICES;
    1 -1 -1 -1
    1 1 -1 -1
    1 -1 1 -1
    1 1 1 -1
    1 -1 -1 1
    1 1 -1 1
    1 -1 1 1
    1 1 1 1





.. link

.. CODE-BLOCK:: perl

    polymake> print join(", ", $c->list_properties);
    CONE_AMBIENT_DIM, CONE_DIM, FACETS, AFFINE_HULL, VERTICES_IN_FACETS, BOUNDED, FEASIBLE, POINTED, N_VERTICES, N_FACETS, VERTICES, LINEALITY_SPACE




The property ``VERTICES`` was added, but a few others were computed on
the way, too. ``polymake`` applied a sequence of *production rules* that
add new properties to the object that can be computed from the
properties the object already posesses.

What properties *can* be computed for a given object depends on the set
of rules defined for it. Here is a short sequence of commands which lets
you find out.


.. link

.. CODE-BLOCK:: perl

    polymake> $t=$c->type;
    polymake> print join(", ", sorted_uniq(sort { $a cmp $b } map { keys %{$_->properties} } $t, @{$t->super}));
    AFFINE_HULL, ALTSHULER_DET, BALANCE, BALANCED, BOUNDARY_LATTICE_POINTS, BOUNDED, CANONICAL, CD_INDEX_COEFFICIENTS, CENTERED, CENTERED_ZONOTOPE, CENTRALLY_SYMMETRIC, CENTROID, CHIROTOPE, CIRCUITS, COCIRCUITS, COCIRCUIT_EQUATIONS, COCUBICAL, COCUBICALITY, COMBINATORIAL_DIM, COMPLEXITY, COMPRESSED, CONE_AMBIENT_DIM, CONE_DIM, CS_PERMUTATION, CUBICAL, CUBICALITY, CUBICAL_H_VECTOR, DEGREE_ONE_GENERATORS, DUAL_BOUNDED_H_VECTOR, DUAL_GRAPH, DUAL_H_VECTOR, EDGE_ORIENTABLE, EDGE_ORIENTATION, EHRHART_POLYNOMIAL_COEFF, EQUATIONS, EXCESS_RAY_DEGREE, EXCESS_VERTEX_DEGREE, F2_VECTOR, FACETS, FACETS_THRU_INPUT_RAYS, FACETS_THRU_POINTS, FACETS_THRU_RAYS, FACETS_THRU_VERTICES, FACET_SIZES, FACET_VERTEX_LATTICE_DISTANCES, FACET_WIDTH, FACET_WIDTHS, FACE_SIMPLICITY, FAR_FACE, FAR_HYPERPLANE, FATNESS, FEASIBLE, FLAG_VECTOR, FOLDABLE_COCIRCUIT_EQUATIONS, FOLDABLE_MAX_SIGNATURE_UPPER_BOUND, FTR_CYCLIC_NORMAL, FTV_CYCLIC_NORMAL, FULL_DIM, F_VECTOR, FacetPerm, FacetPerm.pure, GALE_TRANSFORM, GALE_VERTICES, GORENSTEIN, GORENSTEIN_CONE, GORENSTEIN_INDEX, GORENSTEIN_VECTOR, GRAPH, GROEBNER_BASIS, GROUP, G_VECTOR, HASSE_DIAGRAM, HILBERT_BASIS_GENERATORS, HILBERT_SERIES, HOMOGENEOUS, H_STAR_VECTOR, H_VECTOR, INEQUALITIES, INEQUALITIES_THRU_RAYS, INEQUALITIES_THRU_VERTICES, INPUT_LINEALITY, INPUT_RAYS, INPUT_RAYS_IN_FACETS, INPUT_RAY_LABELS, INTERIOR_LATTICE_POINTS, INTERIOR_RIDGE_SIMPLICES, LATTICE, LATTICE_BASIS, LATTICE_CODEGREE, LATTICE_DEGREE, LATTICE_EMPTY, LATTICE_POINTS_GENERATORS, LATTICE_VOLUME, LATTICE_WIDTH, LATTICE_WIDTH_DIRECTION, LINEALITY_DIM, LINEALITY_SPACE, LINEAR_SPAN, LP, MAHLER_VOLUME, MAX_INTERIOR_SIMPLICES, MILP, MINIMAL_NON_FACES, MINIMAL_VERTEX_ANGLE, MINKOWSKI_CONE, MOEBIUS_STRIP_EDGES, MOEBIUS_STRIP_QUADS, MONOID_GRADING, NEIGHBORLINESS, NEIGHBORLY, NEIGHBOR_RAYS_CYCLIC_NORMAL, NEIGHBOR_VERTICES_CYCLIC_NORMAL, NORMAL, N_01POINTS, N_BOUNDARY_LATTICE_POINTS, N_BOUNDED_VERTICES, N_FACETS, N_HILBERT_BASIS, N_INPUT_RAYS, N_INTERIOR_LATTICE_POINTS, N_LATTICE_POINTS, N_POINTS, N_RAYS, N_RAY_FACET_INC, N_VERTEX_FACET_INC, N_VERTICES, ONE_RAY, ONE_VERTEX, POINTED, POINTS, POINTS_IN_FACETS, POINT_LABELS, POLAR_SMOOTH, QUOTIENT_SPACE, Q_GORENSTEIN_CONE, Q_GORENSTEIN_CONE_INDEX, RAYS, RAYS_IN_FACETS, RAYS_IN_INEQUALITIES, RAYS_IN_RIDGES, RAY_LABELS, RAY_SEPARATORS, RAY_SIZES, REFLEXIVE, RELATIVE_VOLUME, REL_INT_POINT, RIF_CYCLIC_NORMAL, SCHLEGEL_DIAGRAM, SIMPLE, SIMPLEXITY_LOWER_BOUND, SIMPLE_POLYHEDRON, SIMPLICIAL, SIMPLICIALITY, SIMPLICIAL_CONE, SIMPLICITY, SMOOTH, SMOOTH_CONE, SPECIAL_FACETS, SPLITS, SPLIT_COMPATIBILITY_GRAPH, SQUARED_RELATIVE_VOLUMES, STEINER_POINT, STEINER_POINTS, SUBRIDGE_SIZES, TERMINAL, TILING_LATTICE, TOWARDS_FAR_FACE, TRIANGULATION, TRIANGULATION_INT, TWO_FACE_SIZES, UNBOUNDED_FACETS, VALID_POINT, VERTEX_BARYCENTER, VERTEX_LABELS, VERTEX_NORMALS, VERTEX_SIZES, VERTICES, VERTICES_IN_FACETS, VERTICES_IN_INEQUALITIES, VERTICES_IN_RIDGES, VERY_AMPLE, VIF_CYCLIC_NORMAL, VOLUME, VertexPerm, VertexPerm.pure, WEAKLY_CENTERED, ZONOTOPE_INPUT_POINTS




Instead of showing the (lengthy) enumeration have a look at the
`documentation <https://polymake.org/release_docs/latest/polytope.html>`__
for a complete list of properties known for objects of the application
``polytope``.

Schedules
^^^^^^^^^

You may wonder what sequence of rules led to the computation of a
property you request. There usually are several mathematical ways to
compute a property. ``polymake`` uses a nice scheduling algorithm to
find the most efficient procedure, and you can look at what it returns.

Suppose we want to see which sequence of rules leads to the computation
of the F_VECTOR.


.. link

.. CODE-BLOCK:: perl

    polymake> $schedule=$c->get_schedule("F_VECTOR");
    polymake> print join("\n", $schedule->list);
    LINEALITY_DIM : LINEALITY_SPACE
    COMBINATORIAL_DIM : CONE_DIM, LINEALITY_DIM
    precondition : COMBINATORIAL_DIM ( F_VECTOR : N_FACETS, N_RAYS, COMBINATORIAL_DIM )
    F_VECTOR : N_FACETS, N_RAYS, COMBINATORIAL_DIM




So if you ask for the f-vector, ``polymake`` will first compute the
dimension of the lineality space from the basis of the lineality space,
then compute the combinatorial dimension from the lineality and cone
dimensions, and then compute the f-vector from the number of facets,
number of rays, and combinatorial dimension of the polytope. Applying
the schedule to the object yields the same as asking for the property
right away:


.. link

.. CODE-BLOCK:: perl

    polymake> $schedule->apply($c);
    polymake> print join(", ", $c->list_properties);
    CONE_AMBIENT_DIM, CONE_DIM, FACETS, AFFINE_HULL, VERTICES_IN_FACETS, BOUNDED, FEASIBLE, POINTED, N_VERTICES, N_FACETS, VERTICES, LINEALITY_SPACE, LINEALITY_DIM, COMBINATORIAL_DIM, F_VECTOR




As you can see, the things ``polymake`` needed to compute in order to
get to the f-vector are stored in the object as well, so you don’t have
to recompute them later.

If you’re interested, read more about rule scheduling in the `scripting
guide <:user_guide:howto:scripting#rule_planning>`__ and the article on
`writing rules yourself <.extend/rules>`__.
