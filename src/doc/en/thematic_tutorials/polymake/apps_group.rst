.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Application ``group``: Groups in polymake
=========================================

``polymake`` can deal with symmetry groups acting on polytopes,
point/vector configurations and simplicial complexes. This can take the
guise of permutation or matrix groups acting on the set of vertices
(points), facets, or coordinates of realizations of polytopes or
point/vector configurations, or abstract permutation groups acting on
the set of vertices or facets of a simplicial complex.

Some functionality builds on ``PermLib``, a ``C++``-library for
permutation computations written by `Thomas
Rehn <http://www.math.uni-rostock.de/~rehn/index.html>`__, but much
functionality is built in natively into ``polymake``. For example, you
can natively calculate conjugation classes of permutation or matrix
groups, projectors to the isotypic components of representations, or the
invariant polynomials of a matrix representation.

General properties of groups
----------------------------

We start with the description of permutation groups in ``polymake``. An
object of type
`Group <https://polymake.org/release_docs/latest/group.html#group__Group__5>`__
records the abstract properties of the groups that do not depend on any
particular representation, which essentially are just the ``ORDER``,
``CHARACTER_TABLE``, and ``CONJUGACY_CLASS_SIZES``. Moreover, a
``Group`` object can contain several subobjects that encode actions
(representations) of the group, most notably a ``PERMUTATION_ACTION``
that encodes permutations of indices. If the Group object is contained
inside a ``Cone``, ``Polytope``, ``PointConfiguration``, or
``VectorConfiguration``, it may be encoded more specifically as a
``RAY_ACTION``, ``FACET_ACTION``, etc. See the
`documentation <https://polymake.org/release_docs/latest/group.html>`__
for more information on the action types.

In order to access the complete set of functions dealing with groups,
you should switch to the corresponding
`application <:user_guide:lingo#%20application>`__.


.. link

.. CODE-BLOCK:: perl

    polymake> application "group";

Conjugacy classes and character tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As we mentioned before, the only properties that are stored directly in
the ``Group`` object are the order, the sizes of the conjugacy classes,
and sometimes the character table. Other properties, such as
representatives for the conjugacy classes themselves, depend on the
action chosen, and are thus stored inside the corresponding ``*_ACTION``
object.

For the symmetric groups, we currently include the ``CHARACTER_TABLE``
up to degree 7; for the ``group::dihedral_group`` we provide the
conjugacy classes in full generality, and the ``CHARACTER_TABLE`` in an
exact form for the dihedral groups of order 10, 16, 20, and 24 (where it
can be expressed using a single quadratic extension of the rationals).
For the other dihedral groups, we provide ``AccurateFloat``
representations of the entries of the character table, as ``polymake``
currently cannot work with arbitrary cyclotomic integers.


.. link

.. CODE-BLOCK:: perl

    polymake> print dihedral_group(20)->CHARACTER_TABLE;
    1 1 1 1 1 1 1 1
    1 1 1 1 1 1 -1 -1
    1 -1 1 -1 1 -1 1 -1
    1 -1 1 -1 1 -1 -1 1
    2 1/2+1/2r5 -1/2+1/2r5 1/2-1/2r5 -1/2-1/2r5 -2 0 0
    2 -1/2+1/2r5 -1/2-1/2r5 -1/2-1/2r5 -1/2+1/2r5 2 0 0
    2 1/2-1/2r5 -1/2-1/2r5 1/2+1/2r5 -1/2+1/2r5 -2 0 0
    2 -1/2-1/2r5 -1/2+1/2r5 -1/2+1/2r5 -1/2-1/2r5 2 0 0
        





.. link

.. CODE-BLOCK:: perl

    polymake> print dihedral_group(22)->CHARACTER_TABLE;
    1 1 1 1 1 1 1
    1 1 1 1 1 1 -1
    2 3788669096982621/2251799813685248 3741725795518811/4503599627370496 -2563716210467435/9007199254740992 -2949230557375555/2251799813685248 -8642344396869719/4503599627370496 0
    2 3741725795518811/4503599627370496 -2949230557375555/2251799813685248 -1080293049608715/562949953421312 -2563716210467439/9007199254740992 3788669096982621/2251799813685248 0
    2 -2563716210467435/9007199254740992 -1080293049608715/562949953421312 7483451591037615/9007199254740992 3788669096982621/2251799813685248 -2949230557375555/2251799813685248 0
    2 -2949230557375555/2251799813685248 -2563716210467439/9007199254740992 3788669096982621/2251799813685248 -8642344396869719/4503599627370496 7483451591037615/9007199254740992 0
    2 -8642344396869719/4503599627370496 3788669096982621/2251799813685248 -2949230557375555/2251799813685248 7483451591037615/9007199254740992 -2563716210467435/9007199254740992 0
    





Important notice
~~~~~~~~~~~~~~~~

For internal consistency, it is crucial that the
``GROUP->*_ACTION->CONJUGACY_CLASS_REPRESENTATIVES`` (and thus the
``->CONJUGACY_CLASSES`` themselves) be ordered in accordance with the
columns of the ``GROUP->CHARACTER_TABLE``. This is guaranteed to be the
case for the character tables provided natively by ``polymake``, but if
you import character tables from GAP or other sources, correctly
ordering the conjugacy classes is up to you.

Permutation groups
------------------

An instance of a group action is created by specifying its property
``GENERATORS``, a set of permutations or matrices that generates the
group. The Group object itself is then constructed by passing the
action. In the following example we create a symmetric group of degree
3, and then compute its order.


.. link

.. CODE-BLOCK:: perl

    polymake> $p = new PermutationAction(GENERATORS => [[1,0,2],[0,2,1]]);
    polymake> $g = new Group(PERMUTATION_ACTION => $p);
    polymake> print $g->ORDER;
    polymake: used package permlib
      A callable C++ library for permutation computations. 
      Written by Thomas Rehn.
      http://www.math.uni-rostock.de/~rehn/software/permlib.html 
        





::

   6

Of course, there is a user function for creating symmetric groups given
the degree, as well as for several other standard constructions. See the
`docs <https://polymake.org/release_docs/latest/group.html#group__Producing_a_group__15>`__
for a comprehensive list.


.. link

.. CODE-BLOCK:: perl

    polymake> $h = symmetric_group(3);
    polymake> print $h->PERMUTATION_ACTION->GENERATORS;
    1 0 2
    0 2 1
    





Properties of permutation actions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can compute some interesting properties of a PermutationAction:


.. link

.. CODE-BLOCK:: perl

    polymake> $p = new PermutationAction(GENERATORS => [[1,0,2],[0,2,1]]);
    polymake> print all_group_elements($p);
    0 1 2
    0 2 1
    1 0 2
    1 2 0
    2 0 1
    2 1 0
    





There also exist basic functions to compute orbits and stabilizers, for
instance:


.. link

.. CODE-BLOCK:: perl

    polymake> $p = new PermutationAction(GENERATORS => [[1,0,2],[0,2,1]]);
    polymake> $s = stabilizer_of_set($p,new Set<Int>(1,2));
    polymake> print $s->PERMUTATION_ACTION->GENERATORS;
    0 2 1





.. link

.. CODE-BLOCK:: perl

    polymake> print $s->PERMUTATION_ACTION->ORBITS;
    {0}
    {2 1}
    





A note on permutations in polymake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``polymake`` natively uses index notation for permutations: a
permutation g ∈ Sn is an \`Array of length n with entries 0, . . . , n −
1 which corresponds to the second row of the common permutation
notation. For instance, the permutation

::

   0 1 2
   1 0 2

which is equal to (0 1) ∈ S3 in cyclic notation, is represented in
``polymake`` by the integer array [1, 0, 2]. Methods for conversion
between the notation in ``polymake`` and the 1-based cyclic notation as
used, for instance, in GAP are also available:


.. link

.. CODE-BLOCK:: perl

    polymake> $p = new PermutationAction(GENERATORS=>[[1,0,2],[0,2,1]]);
    polymake> print action_to_cyclic_notation($p);
    (1,2),
    (2,3)





.. link

.. CODE-BLOCK:: perl

    polymake> $AGL_1_5 = group_from_cyclic_notation1("(2,3,4,5), (1,2,3,5,4)");
    polymake> print $AGL_1_5->PERMUTATION_ACTION->GENERATORS;
    0 2 3 4 1
    1 2 4 0 3
    





Symmetry groups of polymake objects
-----------------------------------

We switch to the polytope application for the following section:


.. link

.. CODE-BLOCK:: perl

    polymake> application 'polytope';

Polytopes
~~~~~~~~~

There is more than one way to associate a group with any given polytope,
depending on which kind of structural information you want to preserve.
You can find some functions concerning symmetry groups of polytopes
`here <https://polymake.org/release_docs/latest/polytope.html#polytope__Symmetry__36>`__.
It is possibile to attach the group objects described above to polytopes
or cones by using the property ``GROUP``. As there are many possible
groups that operate on a polytope, the property can contain multiple
subobjects; see
`here <https://polymake.org/doku.php/scripting/start#multiple_subobjects>`__
for information on how to handle those.

One interesting group is the group of *combinatorial* automorphisms, the
ones preserving the face lattice. Since the face lattice of a polytope
is atomic and coatomic this group coincides with group of (bipartite)
graph automorphisms of the vertex/facet incidences.


.. link

.. CODE-BLOCK:: perl

    polymake> $c = cube(3);
    polymake> $aut = automorphisms($c->VERTICES_IN_FACETS);
    polymake> print $aut;
    (<0 1 4 5 2 3> <0 1 4 5 2 3 6 7>)
    (<2 3 0 1 4 5> <0 2 1 3 4 6 5 7>)
    (<1 0 2 3 4 5> <1 0 3 2 5 4 7 6>)
    





This says that the combinatorial automorphisms are generated by three
elements, one per line in the output. Each generator is written as a
pair of permutations. The first one gives the action on the FACETS, the
second one gives the action on the VERTICES. Note that ``automorphisms``
does not necessarily output a minimal representation.

Let’s wrap some of this information up in a Group object:


.. link

.. CODE-BLOCK:: perl

    polymake> @g = map { $_->first } @{$aut};
    polymake> $fperm = new group::PermutationAction(GENERATORS=>\@g);
    polymake> $g = new group::Group(FACETS_ACTION=>$fperm);           # note how we use the FACETS_ACTION property this time
    polymake> $g->name = "fullCombinatorialGroupOnFacets";            # is is advisable to give multiple objects a meaningful name
    polymake> $c->add("GROUP",$g);

Now we can, e.g., compute the generators of the action on the vertices
from the action on the facets:


.. link

.. CODE-BLOCK:: perl

    polymake> print $c->GROUP->VERTICES_ACTION->GENERATORS;
    0 1 4 5 2 3 6 7
    0 2 1 3 4 6 5 7
    1 0 3 2 5 4 7 6
    





Many standard constructions of polytopes come with an option to compute
the canonical symmetry group during construction in a more efficient way
than computing the face lattice and then solving the graph automorphism
problem. If you type the name of the function you want to execute and
then hit F1 twice, the available options will be displayed. You will
find a description of the action that will be computed too. For example,
the following creates a cube, but with the action on the facets already
attached:


.. link

.. CODE-BLOCK:: perl

    polymake> $cg = cube(3,group=>1);
    polymake> print $cg->GROUP->FACETS_ACTION->GENERATORS;
    1 0 2 3 4 5
    2 3 0 1 4 5
    0 1 4 5 2 3
    





Orbit polytopes
^^^^^^^^^^^^^^^

Given a group with either a ``COORDINATE_ACTION`` or a
``MATRIX_ACTION``, you can calculate the convex hull of the orbits of a
tuple of points:


.. link

.. CODE-BLOCK:: perl

    polymake> $cg = cube(3,group=>1);
    polymake> print orbit_polytope(new Matrix([[1,1,2,1],[1,5/2,1,0]]), $cg->GROUP->MATRIX_ACTION)->N_VERTICES;
    48
    





See `the
documentation <https://polymake.org/release_docs/latest/polytope.html#polytope__orbit_polytope__319>`__
for more options.

Quotient spaces
~~~~~~~~~~~~~~~

One way of constructing interesting topological spaces is by identifying
points on the boundary of a fundamental region. Polymake can do this in
the case where the fundamental region is a convex polytope. For example,
a cylinder is obtained by identifying opposite sides of a square, and
the
`quarter_turn_manifold() <https://polymake.org/release_docs/latest/polytope.html#polytope__quarter_turn_manifold__238>`__
(see
`here <http://www.math.cornell.edu/~dwh/books/eg99/Ch20/Ch20.html>`__)
is obtained from the boundary of a 3-dimensional cube by identifying
opposite faces by a quarter turn.

For example, to obtain a topological space homeomorphic to a cylinder,
type


.. link

.. CODE-BLOCK:: perl

    polymake> $p = cylinder_2();
    polymake> print $p->QUOTIENT_SPACE->IDENTIFICATION_ACTION->GENERATORS;
    2 3 0 1





.. link

.. CODE-BLOCK:: perl

    polymake> print $p->QUOTIENT_SPACE->IDENTIFICATION_ACTION->ORBITS;
    {0 2}
    {1 3}





.. link

.. CODE-BLOCK:: perl

    polymake> print $p->QUOTIENT_SPACE->FACES;
    {{0} {1}}
    {{0 1} {0 2} {1 3}}
    {{0 1 2 3}}





.. link

.. CODE-BLOCK:: perl

    polymake> print $p->QUOTIENT_SPACE->F_VECTOR;
    2 3 1
    





Thus, vertices 0,2 and vertices 1,3 of a square (a 2-dimensional cube)
are identified, and after identification two vertices, three edges, and
one two-dimensional face remain. In order to get a simplicial complex
without identifications among the vertices, you can calculate the second
barycentric subdivision by asking for the property SIMPLICIAL_COMPLEX:


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->QUOTIENT_SPACE->SIMPLICIAL_COMPLEX->F_VECTOR;
    26 72 48





.. link

.. CODE-BLOCK:: perl

    polymake> print $p->QUOTIENT_SPACE->SIMPLICIAL_COMPLEX->HOMOLOGY;
    ({} 0)
    ({} 0)
    ({} 1)
    





An easy way to make projective spaces is to identify opposite faces in a
centrally symmetric polytope, using the function
`cs_quotient() <https://polymake.org/release_docs/latest/polytope.html#polytope__cs_quotient__239>`__.
For example, to calculate the homology of real 3-dimensional projective
space \**RP3, write


.. link

.. CODE-BLOCK:: perl

    polymake> $m = cs_quotient(cube(3));
    polymake> print $m->QUOTIENT_SPACE->SIMPLICIAL_COMPLEX->HOMOLOGY;
    ({} 0)
    ({(2 1)} 0)
    ({} 0)
    ({} 1)
    





As another example, the `Davis
Manifold <https://people.math.osu.edu/davis.12/old_papers/4-mfld.pdf>`__
is a 4-dimensional hyperbolic manifold obtained by identifying opposite
vertices of a 120-cell:


.. link

.. CODE-BLOCK:: perl

    polymake> $m=davis_manifold();
    polymake> print $m->QUOTIENT_SPACE->F_VECTOR;
    300 600 360 60 1
    





Calculating the homology takes a little bit longer:

polytope > print $m->QUOTIENT_SPACE->SIMPLICIAL_COMPLEX->F_VECTOR;
 94321 1146960 3644640 4320000 1728000
polytope > print $m->QUOTIENT_SPACE->SIMPLICIAL_COMPLEX->HOMOLOGY;
 ({} 0)
 ({(2 1)} 0)
 ({} 0)
 ({(2 1)} 0)
 ({} 0)


Matrix groups
-------------

Let’s switch back to ``group``.


.. link

.. CODE-BLOCK:: perl

    polymake> application 'group';

Polymake can also deal with groups given by matrices that act on the
ambient space. They are stored in the property ``GROUP.MATRIX_ACTION``,
and are paramterized by the number type of the matrices. One way to get
a ``MATRIX_ACTION`` is to convert a permutation action on the vertices
of a polytope:


.. link

.. CODE-BLOCK:: perl

    polymake> $d = polytope::dodecahedron();
    polymake> $d->GROUP->properties();
    type: Group as Polytope<QuadraticExtension<Rational>>::GROUP
        





::

   VERTICES_ACTION
   type: PermutationAction<Int, Rational>


.. link

.. CODE-BLOCK:: perl

    polymake> $d->GROUP->MATRIX_ACTION;




.. link

.. CODE-BLOCK:: perl

    polymake> print $d->GROUP->MATRIX_ACTION->GENERATORS;
    <1 0 0 0
    0 -1 0 0
    0 0 1 0
    0 0 0 1
    >
    <1 0 0 0
    0 1/4-1/4r5 1/2 -1/4-1/4r5
    0 1/2 1/4+1/4r5 -1/4+1/4r5
    0 -1/4-1/4r5 -1/4+1/4r5 1/2
    >
    <1 0 0 0
    0 1 0 0
    0 0 1 0
    0 0 0 -1
    >
    





As we can see, the property ``MATRIX_ACTION`` was calculated on the fly,
specifically by solving matrix equations involving the ``VERTICES`` and
``VERTICES_ACTION->GENERATORS``. Moreover, in this case the matrices are
calculated exactly by adjoining the square root of 5 to the rationals.

Of course, not every combinatorial symmetry group of a concrete point
configuration has a realization as a matrix group, in which case the
above computation will fail. A sure-fire way to get a matrix group is to
calculate the ``REGULAR_REPRESENTATION`` of a permutation group, which
yields the action by permutation matrices on the ambient space of
dimension = number of points.

Orbits
~~~~~~

Once you have a matrix group, you may calculate the orbit of an
arbitrary vector under it:


.. link

.. CODE-BLOCK:: perl

    polymake> $s = symmetric_group(3); 




.. link

.. CODE-BLOCK:: perl

    polymake> $a = $s->REGULAR_REPRESENTATION;




.. link

.. CODE-BLOCK:: perl

    polymake> print orbit($a->GENERATORS, new Vector([1,2,3]));
    {<3 2 1> <1 2 3> <2 1 3> <1 3 2> <3 1 2> <2 3 1>}
    





Invariant polynomials
~~~~~~~~~~~~~~~~~~~~~

Or you can regard the matrices as acting on polynomials, and calculate a
set of invariant polynomials of a given maximum degree. For this, recall
that the action of a matrix on a polynomial is exemplified by

::

   [ 1  1 ]
   [ 1 -1 ]  .  ( x^2 - y^2 )  =  ( x + y )^2 - ( x - y )^2.

You can calculate the polynomials left invariant by the matrices sending
the vertices of a dodecahedron into each other as follows:


.. link

.. CODE-BLOCK:: perl

    polymake> $d = polytope::dodecahedron();




.. link

.. CODE-BLOCK:: perl

    polymake> $d->GROUP->MATRIX_ACTION;




.. link

.. CODE-BLOCK:: perl

    polymake> print join "\n", @{invariant_polynomials($d->GROUP->MATRIX_ACTION, 5)};
    x_0^2 + x_1^2 + x_2^2
    x_0^4 + 2*x_0^2*x_1^2 + 2*x_0^2*x_2^2 + x_1^4 + 2*x_1^2*x_2^2 + x_2^4
    





This is consistent with the Molien series of this action starting out as
1 + x^2 + x^4 + 2x^6 + …, so in particular no invariant of degree
exactly 5 is found. See `this
paper <http://www.ams.org/journals/bull/1979-01-03/S0273-0979-1979-14597-X/S0273-0979-1979-14597-X.pdf>`__
by Stanley for more information.

Decomposition into irreps, and bases of isotypic components
-----------------------------------------------------------

You can calculate

-  the character of a permutation action or matrix action,

-  the decomposition of the action into irreducible representations, and

-  the projection operators to (and vector space bases of) the isotypic
   components.

For ``MATRIX_ACTION``\ s, the character can always be calculated, but
for the rest of these computations the ``CHARACTER_TABLE`` must be
known:


.. link

.. CODE-BLOCK:: perl

    polymake> print $d->GROUP->MATRIX_ACTION->CHARACTER;
    -2 0 2 1 4 1 1/2-1/2r5 3/2-1/2r5 1/2+1/2r5 3/2+1/2r5
        





::

    group > print irreducible_decomposition($d->GROUP->MATRIX_ACTION->CHARACTER, $d->GROUP);
   polymake:  WARNING: available properties insufficient to compute 'CHARACTER_TABLE'

This didn’t work, because the dodecahedron doesn’t (yet) come with a
character table; this might change in future versions, though.

It does work, for instance, for the symmetric group of order 5! (in
fact, up to order 7!):


.. link

.. CODE-BLOCK:: perl

    polymake> $s=symmetric_group(5);




.. link

.. CODE-BLOCK:: perl

    polymake> print $s->CHARACTER_TABLE;
    1 -1 1 1 -1 -1 1
    4 -2 0 1 1 0 -1
    5 -1 1 -1 -1 1 0
    6 0 -2 0 0 0 1
    5 1 1 -1 1 -1 0
    4 2 0 1 -1 0 -1
    1 1 1 1 1 1 1
        





.. link

.. CODE-BLOCK:: perl

    polymake> $s->REGULAR_REPRESENTATION;




.. link

.. CODE-BLOCK:: perl

    polymake> print $s->REGULAR_REPRESENTATION->CHARACTER;
    5 3 1 2 0 1 0
        





.. link

.. CODE-BLOCK:: perl

    polymake> print irreducible_decomposition($s->REGULAR_REPRESENTATION->CHARACTER,$s);
    0 0 0 0 0 1 1
    





So the regular (permutation) representation decomposes into one copy
each of the invariant subspaces associated to the characters in the last
two lines of the character table. The first entries there, 4 and 1, say
that these components should have dimensions 4 and 1, respectively:


.. link

.. CODE-BLOCK:: perl

    polymake> print isotypic_basis($s, $s->REGULAR_REPRESENTATION, 5);
    4/5 -1/5 -1/5 -1/5 -1/5
    -1/5 4/5 -1/5 -1/5 -1/5
    -1/5 -1/5 4/5 -1/5 -1/5
    -1/5 -1/5 -1/5 4/5 -1/5
        





.. link

.. CODE-BLOCK:: perl

    polymake> print isotypic_basis($s, $s->REGULAR_REPRESENTATION, 6);
    1/5 1/5 1/5 1/5 1/5
    






