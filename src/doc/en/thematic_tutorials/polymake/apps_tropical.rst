.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Tutorial: Tropical arithmetics and tropical geometry
====================================================

This tutorial showcases the main features of application tropical, such
as

-  Tropical arithmetics
-  Tropical convex hull computations
-  Tropical cycles and hypersurfaces.
-  Tropical intersection theory

To use the full palette of tools for tropical geometry, switch to the
corresponding application by typing the following in the ``polymake``
shell:


::

    polymake> application 'tropical';

Disclaimer: Min or Max - you have to choose!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most objects and data types related to tropical computations have a
template parameter which tells it whether Min or Max is used as tropical
addition. There is **no default** for this, so you have to choose!

Disclaimer 2: Newest version required
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most of the features described here only work in polymake version 3.0 or
newer.


a-tint
^^^^^^

As of version 2.15-beta3, polymake comes bundled with the extension
`a-tint <https://github.com/simonhampe/atint>`__ by Simon Hampe, which
specializes in (but is not limited to) tropical intersection theory, and
in fact supplies a big part of polymakes functionality related to
tropical geometry. You can find a non-comprehensive list of features
`here <https://github.com/simonhampe/atint/wiki/Feature-list>`__.


Tropical arithmetics
--------------------

You can create an element of the tropical semiring (over the rationals)
simply by writing something like this:


::

    polymake> $a = new TropicalNumber<Max>(4);
    ........> $b = new TropicalNumber<Min>(4);
    ........> $c = new TropicalNumber<Min>("inf");

You can now do basic arithmetic - that is **tropical** addition and
multiplication with these. Note that tropical numbers with different
tropical additions don’t mix!


::

    polymake> print $a * $a;

::

    polymake> print $b + $c*$b;

::

    polymake> # print $a + $b; # this won't work!

Tropical vector/matrix arithmetics also work - you can even ask for the
tropical determinant!


::

    polymake> $m = new Matrix<TropicalNumber<Max>>([[0,1,2],[0,"-inf",3],[0,0,"-inf"]]);
    ........> $v = new Vector<TropicalNumber<Max>>(1,1,2);
    ........> print $m + $m;

::

    polymake> print $m * $v;

::

    polymake> print tdet($m);

Finally, you can also create tropical polynomials. This can be done with
the special toTropicalPolynomial parser:


::

    polymake> $q = toTropicalPolynomial("min(2a,b+c)");
    ........> print $q;

Homogeneous and extended coordinates
------------------------------------

All higher tropical objects (trop. polytopes, cycles, …) live in the
*tropical projective torus*
`\mathbb R^n/\mathbb R\mathbf 1\cong \mathbb R^{n-1}`, i.e. they
are defined in (tropically) homogeneous coorditates, and points are
given via their unique representative whose first entry is `0`
(e.g. `(1,2,3)` becomes `(0,1,2)` etc.).

At the same time many properties (e.g. ``VERTICES``) expect coordinates
to be “extended” with an additional leading entry, which signals whether
the vector should be interpreted as a conventional point in the torus
(leading entry `\mathbf 1`), or as a ray/direction/“point at
infinity” (leading entry `\mathbf 0`). You can read
`here <coordinates.ipynb>`__ for how this makes sense, or just think of
the leading entry as purely symbolic.

Note that this contrasts to conventional real/complex/… projective
coordinates, where this extension is not needed, since there the leading
entry naturally indicates either an affine point (`1`) or point at
infinity (`0`) - again, see `here <coordinates.ipynb>`__.

**Example**: the standard tropical max-line in the trop. 2-torus can be
defined by 4 points in extended homogeneous coordinates: a point
`(\mathbf 1,0,0,0)` and 3 “rays” `(\mathbf 0,0,-1,0)`,
`(\mathbf 0,0,0,-1)`, and `(\mathbf 0,0,1,1)`.


Helper functions
^^^^^^^^^^^^^^^^

Since it can be tedious to write everything in tropically homogeneous
coordinates, especially when you are already working with affine
coordinates, Polymake provides the functions ``thomog`` and ``tdehomog``
to homogenize or dehomogenize (extended) coordinates, respectively.

Both have signature ``(matrix, chart=0, has_leading_coordinate=1)``,
where ``matrix`` is a list of vectors that you want to (de)homogenize,
``has_leading coordinate`` is a boolean indicating whether we are
dealing with extended coordinates (i.e. vectors have an additional
leading coordinate), and ``chart`` is the index of the coordinate that
is set to `0` when identifying
`\mathbb R^n/\mathbb R\mathbf 1` with `\mathbb R^{n-1}`
(i.e. homogenizing works by simply adding a `0` entry at index
``chart``, dehomogenizing works by taking the representative whose entry
at position ``chart`` is `0`, and then deleting that entry), all
while ignoring the first entry if ``has_leading_coordinate`` is ``1``.

Some examples:


::

    polymake> print thomog([[1,3,4],[0,5,6]]);

::

    polymake> print thomog([[2,3,4]], 1, 0);

::

    polymake> print tdehomog([[1,3,4,5]]);

::

    polymake> print tdehomog([[2,3,4,5]], 1, 0);

Tropical convex hull computations
---------------------------------

The basic object for tropical convex hull computations is ``Polytope``
(**Careful:** If you’re not in application tropical, be sure to use the
namespace identifier ``tropical::Polytope`` to distinguish it from the
``polytope::Polytope``).

A tropical polytope should always be created via ``POINTS`` (i.e. not
``VERTICES``), since they determine the combinatorial structure. The
following creates a tropical line segment in the tropical projective
plane. Note that the point (0,1,1) is not a vertex, as it is in the
tropical convex hull of the other two points. However, it does play a
role when computing the corresponding subdivision of the tropical
projective torus into covector cells (see the
`note <apps_tropical#A%20note%20on%20coordinates>`__ below to understand
the different coordinates):


::

    polymake> $c = new Polytope<Min>(POINTS=>[[0,0,0],[0,1,1],[0,2,1]]);
    ........> print $c->VERTICES;

::

    polymake> print rows_labeled($c->PSEUDOVERTICES);

::

    polymake> print $c->MAXIMAL_COVECTOR_CELLS; #Sets of PSEUDOVERTICES. They are maximal cells of the induced subdivision of the torus.

::

    polymake> print $c->POLYTOPE_MAXIMAL_COVECTOR_CELLS;

::

    polymake> $c->VISUAL_SUBDIVISION;

In case you’re just interested in either the subdivision of the full
torus, or the polyhedral structure of the tropical polytope, the
following will give you those structures as ``fan::PolyhedralComplex``
objects in *affine* coordinates:


::

    polymake> $t = $c->torus_subdivision_as_complex;
    ........> $p = $c->polytope_subdivision_as_complex;
    ........> print $p->VERTICES;

::

    polymake> print $p->MAXIMAL_POLYTOPES;

Note that by default, the affine chart is {x_0 = 0}. You can choose any
chart {x_i = 0} by passing i as an argument to
``.._subdivision_as_complex``.

Polymake computes the full subdivision of both the torus and the
polytope as a ``CovectorLattice``, which is just a ``FaceLattice`` with
an additional map that attaches to each cell in the subdivision its
covector. For more details on this data structure see the `reference
documentation <http://polymake.org/release_docs/snapshot/tropical.html>`__.
You can visualize the covector lattice with


::

    polymake> $c->TORUS_COVECTOR_DECOMPOSITION->VISUAL;
    ........> $c->POLYTOPE_COVECTOR_DECOMPOSITION->VISUAL;

Each node in the lattice is a cell of the subdivision. The top row
describes the vertices and rays of the subdivision. The bottom row is
the covector of that cell with respect to the ``POINTS``.


Tropical cycles
---------------


The main object here is ``Cycle``, which represents a weighted and
balanced, rational pure polyhedral complex in the tropical projective
torus (see the `note <apps_tropical#A%20note%20on%20coordinates>`__
below, if you’re confused by coordinates in the following examples).

A tropical cycle can be created, like a ``PolyhedralComplex``, by
specifying its vertices and maximal cells (and possibly a lineality
space). The only additional data are the weights on the maximal cells.


::

    polymake> $x = new Cycle<Max>(PROJECTIVE_VERTICES=>[[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]],MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3]],WEIGHTS=>[1,1,1]);

This creates the standard tropical (max-)line in the plane. There are
two caveats to observe here: 1. The use of ``POINTS`` and
``INPUT_POLYTOPES`` is strongly discouraged. ``WEIGHTS`` always refer to
``MAXIMAL_POLYTOPES`` and the order of the latter can be different from
the order in ``INPUT_POLYTOPES``. 2. You can also define a cycle using
``VERTICES`` instead of ``PROJECTIVE_VERTICES``. However, in that case
all vertices have to be normalized such that the second coordinate
(i.e. the one after the leading 0/1, see
`note <apps_tropical#A%20note%20on%20coordinates>`__) is 0. I.e. in the
above example, the point (0,-1,0,0) would have to be replaced by
(0,0,1,1).

Entering projective coordinates can be a little tedious, since it
usually just means adding a zero in front of your affine coordinates.
There is a convenience function that does this for you. The following
creates the excact same cycle as above:


::

    polymake> $x = new Cycle<Max>(VERTICES=>thomog([[1,0,0],[0,1,1],[0,-1,0],[0,0,-1]]),MAXIMAL_POLYTOPES=>[[0,1],[0,2],[0,3]],WEIGHTS=>[1,1,1]);

One can now ask for basic properties of the cycle, e.g., if it’s
balanced:


::

    polymake> print is_balanced($x);

Hypersurfaces
^^^^^^^^^^^^^

Most of the time you probably won’t want to input your tropical cycle
directly as above. Polymake has a special data type ``Hypersurface`` for
hypersurfaces of *homogeneous* tropical polynomials. The following
creates the standard tropical min-line in the plane:


::

    polymake> $H = new Hypersurface<Min>(POLYNOMIAL=>toTropicalPolynomial("min(a,b,c)"));
    ........> print $H->VERTICES;

::

    polymake> print $H->MAXIMAL_POLYTOPES;

::

    polymake> print $H->WEIGHTS;

Tropical intersection theory
----------------------------

Functionality related to tropical intersection theory is provided by the
bundled extension atint by Simon Hampe, see
`here <https://github.com/simonhampe/atint/wiki/User-Manual>`__ for a
more extensive, but also slightly outdated documentation.


Basic examples
^^^^^^^^^^^^^^


Computing the divisor of a tropical polynomial in `\mathbb R^n`
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


::

    polymake> $f = toTropicalPolynomial("max(0,x,y,z)");
    ........> $div = divisor(projective_torus<Max>(3), rational_fct_from_affine_numerator($f));
    ........> $div->VISUAL;

Here, ``projective_torus`` creates the tropical projective 3-torus (aka
`\mathbb R^3`) as a tropical fan with weight 1.


Visualizing a curve in a tropical surface
'''''''''''''''''''''''''''''''''''''''''

Let’s create the standard tropical hyperplane in 3-space:


::

    polymake> $l = uniform_linear_space<Min>(3,2);

Furthermore, we compute a curve as the divisor of a rational function on
``$l``:


::

    polymake> $div = divisor($l, rational_fct_from_affine_numerator(toTropicalPolynomial("min(3x+4,x-y-z,y+z+3)")));

We now want to visualize these together and since we want to put the
resulting picture in a slide show, we want the picture to look a bit
nicer, e.g. we want to intersect both the surface and the curve with the
same bounding box (we don’t want the curve to continue “beyond the
surface”).


::

    polymake> print $l->bounding_box(1);

::

    polymake> print $div->bounding_box(1);

The command boundingBox tells us, what a “good” bounding box for a given
polyhedral complex should be. Obviously the bounding box for the curve
is larger, so we want to take this one. ``compose`` is used to visualize
several objects at the same time:


::

    polymake> compose($l->VISUAL(VertexStyle=>"hidden",BoundingMode=>"absolute",BoundingBox=>$div->bounding_box(1)), $div->VISUAL(VertexStyle=>"hidden",EdgeColor=>"red", BoundingMode=>"absolute",BoundingBox=>$div->bounding_box(1)));

Alternatively, we could simply specify an explicit bounding box to make
the surface look more symmetric:


::

    polymake> $m = new Matrix<Rational>([[-5,-5,-5],[5,5,5]]);
    ........> compose($l->VISUAL(VertexStyle=>"hidden",BoundingMode=>"absolute",BoundingBox=>$m), $div->VISUAL(VertexStyle=>"hidden",EdgeColor=>"red",BoundingMode=>"absolute",BoundingBox=>$m));

Creation functions for commonly used tropical cycles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Uniform linear spaces
'''''''''''''''''''''

As a special case of tropical linear spaces, one can create the space
associated to the uniform matroid (with trivial valuation) as
``uniform_linear_space<Addition>(n,k; w=1)`` where \* ``Addition`` is
either ``Min`` or ``Max`` \* ``n`` is the dimension of the ambient
tropical projective torus. \* ``k`` is the (projective) dimension of the
linear space. \* ``w`` is the (global) weight of the space. By default,
this is 1. In particular, this gives the matroid fan of the uniform
matroid `U_{n+1,k+1}`.


True linear spaces
''''''''''''''''''

There are creation function for cycles that are supported on an affine
linear space. In the special case, that a cycle is the whole projective
torus, this can be created via ``projective_torus<Addition>(n;w=1)``
where n is the (projective) dimension and w is an integer weight.

If the cycle is given by a basis of the lineality space and a
translation vector, one can use
``affine_linear_space<Addition>(matrix; vector, w =1)``.

Note that both the generators of the lineality space in the matrix and
the translation vector should be given in tropical homogeneous
coordinates, but WITHOUT a leading coordinate. If the translation vector
is not given, it is equal to 0. The weight w is equal to 1 by default.


Halfspace subdivision
'''''''''''''''''''''

``halfspace_subdivision<Addition>(a,g;w=1)`` creates the subdivision of
the projective torus along the (true) hyperplane defined by the equation
`g\cdot x = a`, where `a` is a Rational and `g` is a
Vector. Note that the sum over all entries in `g` must be 0 for
this to be well-defined over tropical homogeneous coordinates.


Matroid fans
''''''''''''

A-tint uses the algorithms developped by Felipe Rincón to compute
matroidal fans. The original algorithm was developed for matrix matroids
but has been adapted here to also work on general matroids. You compute
a matroidal fan via ``$m = matroid_fan<Addition>(m)``, where ``m`` is
either a Matroid object or a matrix, whose column matroid is then
considered.


Operations on tropical varieties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Computing the cartesian product
'''''''''''''''''''''''''''''''

The function ``cartesian_product`` computes the cartesian product of an
arbitrary, comma-separated list of tropical cycles using the same
tropical addition, e.g. the product of two uniform linear spaces could
be computed via:


::

    polymake> $p = cartesian_product(uniform_linear_space<Max>(3,2), uniform_linear_space<Max>(3,1));

Computing the skeleton
''''''''''''''''''''''

This is actually an operation on polyhedral complexes, since taking the
skeleton forgets the weights. The function
``skeleton_complex(cycle,k,preservesRays)`` computes the k-dimensional
skeleton of fan. The last argument is optional, by default false and
should be set to true, if you are certain that all your rays remain in
the skeleton (directional rays may disappear when taking the skeleton of
a polyhedral complex). In this case, ray indices will be preserved.


Comparing complexes
'''''''''''''''''''

You can check if two polyhedral complexes are exactly the same,
i.e. have the same rays, cones and weights up to reordering and
equivalence of tropical coordinates, by calling the function
``check_cycle_equality``. Note that this function will not recognize if
two complex are equal up to subdivision. The function has an optional
parameter, which is TRUE by default, that determines whether weights
ought to be checked for equality as well. This will be ignored if any of
the complexes does not have any weights.


Refining
''''''''

You can refine a weighted complex along any other complex containing it
(as a set), using the function ``intersect_container(X,Y)``.


Local computations
^^^^^^^^^^^^^^^^^^

a-tint has a mechanism to allow computations locally around a given set
of cones: Each Cycle has a property LOCAL_RESTRICTION. This is an
IncidenceMatrix (or list of sets of ray indices). Each of these sets of
ray indices is supposed to describe a cone of the complex - though not
necessarily a maximal one.

If you define a complex with this property, then all computations will
only be done “around” these cones. E.g. a-tint will only recognize those
codimension one faces that contain one of the cones described in
LOCAL_RESTRICTION.

**Note:** It is important, that all maximal cones of a restricted
complex contain one of the LOCAL_RESTRICTION cones.

This way a bounded complex can also be made “balanced”:


::

    polymake> $w = new Cycle<Max>(VERTICES=>thomog([[1,0],[1,-1],[1,1]]),MAXIMAL_POLYTOPES=>[[0,1],[0,2]],WEIGHTS=>[1,1],LOCAL_RESTRICTION=>[[0]]);
    ........> print $w->IS_BALANCED;

The above example describes the bounded line segment
`[-1,1] \in \mathbb R`, subdivided at 0. Since we restrict
ourselves locally to the vertex 0, the outer codimension one faces are
not detected by a-tint and the complex is balanced.

This can improve computation speed, e.g. when one is only interested in
the weight of a divisor at a certain cone.

There are some convenience functions to create a local variety from a
global one: \* ``local_restrict(X,C)``: Takes a variety and localizes at
a set of cones C. It removes all maximal cones that do not contain one
of these cones. Here, C is an IncidenceMatrix, i.e. it is of the same
form as e.g. MAXIMAL_POLYTOPES. \* ``local_vertex(X,r)``: Takes a
variety and localizes at a ray/vertex, given by its row index r in
VERTICES. \* ``local_codim_one(X,c)``: Takes a variety and localizes at
a codimension one face, given by its row index in
CODIMENSION_ONE_POLYTOPES. \* ``local_point(X,v)``: Takes a variety and
refines it such that a given point v (in homog. coordinates, with
leading coordinate) is in the polyhedral structure. Then it localizes at
this point. This returns an error message if v is not actually in the
complex.

Some examples:


::

    polymake> $x = uniform_linear_space<Max>(3, 2);
    ........> $x1 = local_restrict($x, new IncidenceMatrix([[0],[2,3]])); 
    ........>      #The tropical hyperplane, locally around the 0-th ray and the maximal cone spanned by rays 2 and 3.
    ........>      #As a set this can be interpreted as the union of an open neighborhood of ray 0 and the interior of the maximal       
    ........>      #cone <2,3>
    ........> $x2 = local_vertex($x, 0); #An open neighborhood of ray 0.
    ........> $x3 = local_codim_one($x, 0); #An open neighborhood of the codimension one face no. 0 (which is a ray)
    ........> $x4 = local_point($x, new Vector<Rational>([1,0,1,1,0]));
    ........>      #Refines the surface such that it contains the point (0,1,1,0), then takes an open neighborhood of that

For details see the internal documentation of these functions.


Delocalizing
^^^^^^^^^^^^

If you have a locally restricted complex and you would like to obtain
the same complex without the LOCAL_RESTRICTION, you can call its user
method delocalize(), e.g.


::

    polymake> $y = $x->delocalize();