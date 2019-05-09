.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Getting started with the application ``fan``
============================================

Besides the name giving polyhedral fans this application covers a few
other big objects and related functions. An overview can be found in the
`documentation <https://polymake.org/release_docs/latest/fan.html>`__ or
the `interactive
help <https://polymake.org/doku.php/user_guide/intro_tutorial#getting_help>`__.

However, this tutorial focuses on
`PolyhedralFan <https://polymake.org/release_docs/latest/fan.html#fan__PolyhedralFan__27>`__
objects. Subdivisions have their own tutorial
`here <https://polymake.org/doku.php/user_guide/tutorials/regular_subdivisions>`__
and some notes on polyhedral complexes can be found
`here <https://polymake.org/doku.php/user_guide/tutorials/pcom>`__.

Most of the following code snippets will only work in your polymake
shell after switching to the application ``fan`` with the command


::

    polymake> application 'fan';

Polyhedral fans
---------------

Construction from scratch
~~~~~~~~~~~~~~~~~~~~~~~~~

A primal description containing rays and rays-cones incidence relations
can be passed to the constructor like this:


::

    polymake> $f = new PolyhedralFan(INPUT_RAYS=>[[1,0],[0,1],[-1,0],[0,-1],[2,0]], INPUT_CONES=>[[0,1,4],[1,2],[2,3],[3,0],[0]]);

Former are assigned to
`INPUT_RAYS <https://polymake.org/release_docs/latest/fan.html#fan__INPUT_RAYS__161>`__
as an array of row vectors (which is a matrix). All input rays must not
be zero but redundancies are allowed. The latter are assigned to
`INPUT_CONES <https://polymake.org/release_docs/latest/fan.html#fan__INPUT_CONES__160>`__
and encoded as an array of index sets. Each index set refers to a subset
of ``INPUT_RAYS`` that forms a cone in the fan, indexing starts with
zero. Input rays that do not belong to any of the input cones are
ignored. Input cones do not need to be inclusion-wise maximal. Subcones
of input cones are, however, implicitly included. Indeed, for our fan
``$f`` we obtain:


::

    polymake> print $f->CONES;
    <{0}
    {1}
    {2}
    {3}
    >
    <{0 1}
    {1 2}
    {2 3}
    {0 3}
    >





.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package lrs
      Implementation of the reverse search algorithm of Avis and Fukuda.
      Copyright by David Avis.
      http://cgm.cs.mcgill.ca/~avis/C/lrs.html
    
    polymake: used package cdd
      cddlib
      Implementation of the double description method of Motzkin et al.
      Copyright by Komei Fukuda.
      http://www-oldurls.inf.ethz.ch/personal/fukudak/cdd_home/
    
    </pre>
    </details>




You can specify a fan with lineality by additionally passing
`INPUT_LINEALITY <https://polymake.org/release_docs/latest/fan.html#fan__INPUT_LINEALITY__162>`__.
Nevertheless, a fan given by input rays and input cones can have
lineality as well. Please remind yourself, that all cones in a fan share
the same lineality space.

The properties
`RAYS <https://polymake.org/release_docs/latest/fan.html#fan__RAYS__176>`__,
`MAXIMAL_CONES <https://polymake.org/release_docs/latest/fan.html#fan__MAXIMAL_CONES__150>`__
and
`LINEALITY_SPACE <https://polymake.org/release_docs/latest/fan.html#fan__LINEALITY_SPACE__180>`__
are giving a **non-redundant** primal description:


::

    polymake> print rows_labeled($f->RAYS),"\n";
    ........> print $f->MAXIMAL_CONES,"\n";
    ........> print "lineality dimensions: ", $f->LINEALITY_SPACE->rows() ."x". $f->LINEALITY_SPACE->cols();
    0:1 0
    1:0 1
    2:-1 0
    3:0 -1
    
    {0 1}
    {1 2}
    {2 3}
    {0 3}
    
    lineality dimensions: 0x2




Note that, even though ``LINEALITY_SPACE`` is an empty matrix, its
number of columns is equal to the ambient dimension of ``$f``.

Instead of the input properties, you may right away use ``RAYS``,
``MAXIMAL_CONES`` and ``LINEALITY_SPACE`` for construction purposes but
keep in mind:

Unlike input rays and input cones, only providing rays and maximal cones
may not describe a fan with lineality. In this case polymake assumes an
empty lineality space. All given rays must be non-redundant and in case
of non-pointed fans ``LINEALITY_SPACE`` stores a basis of the lineality
space.


The dual description
^^^^^^^^^^^^^^^^^^^^

The following properties give rise to a dual description:


::

    polymake> print rows_labeled($f->FACET_NORMALS),"\n";
    ........> print rows_labeled($f->MAXIMAL_CONES_FACETS);
    0:1 0
    1:0 1
    
    0:1 1
    1:-1 1
    2:-1 -1
    3:1 -1





Where ``FACET_NORMALS`` is an array of row vectors, the facet normals of
all maximal cones. Incidence relations between them are stored in the
sparse matrix ``MAXIMAL_CONES_FACETS``. Each row corresponds to a
maximal cone and each column to a facet normal. Its entries are 0, 1 or
-1 encoding either no incidence, an inner or and outer facet normal of
the cone, respectively. For example, the second row of
``MAXIMAL_CONES_FACETS`` shows that the first one is an outer and the
second one is an inner facet normal of the second maximal cone.

The dual description requires additional information on the linear span
of each maximal cone. This is stored in ``LINEAR_SPAN_NORMALS`` and
``MAXIMAL_CONES_LINEAR_SPAN_NORMALS``. An empty index set in the latter
corresponds to a full dimensional maximal cone. Check out the
`documentation <https://polymake.org/release_docs/latest/fan.html#fan__MAXIMAL_CONES_LINEAR_SPAN_NORMALS__172>`__
for more informations. All maximal cones in ``$f`` are full dimensional,
hence ``LINEAR_SPAN_NORMALS`` is empty:


::

    polymake> print $f->LINEAR_SPAN_NORMALS->rows."\n\n";
    ........> print $f->MAXIMAL_CONES_LINEAR_SPAN_NORMALS;
    0
    
    {}
    {}
    {}
    {}





Construction from a set of cones
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example one can extract the second and fourth maximal cone of
``$f``:


::

    polymake> $c1 = $f->cone(1);
    ........> $c3 = $f->cone(3);

and pass them to the user method
`check_fan_objects <https://polymake.org/release_docs/latest/fan.html#fan__check_fan_objects__54>`__,
which returns the corresponding ``PolyhedralFan`` object if and only if
the set of provided cones defines a valid polyhedral fan, id est
satisfies the intersection property.


::

    polymake> $checkedfan = check_fan_objects($c1,$c3);
    ........> print $checkedfan->MAXIMAL_CONES;
    {0 1}
    {2 3}





Construction from other objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Polymake provides several clients doing this job.

`normal_fan <https://polymake.org/release_docs/latest/fan.html#fan__normal_fan__45>`__
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The inner normal fan of a polytope can be produced with this client. For
example the normal fan of the 3-dimensional +/-1 cube:


::

    polymake> $nf = normal_fan(cube(3));

Normal fans of bounded feasible polytopes always satisfy the following
properties:


::

    polymake> foreach my $prop (qw(regular pure complete full_dim)) {
    ........>     print ucfirst($prop),": ", $nf->give(uc($prop)),"\n";
    ........> }
    Regular: true
    Pure: true
    Complete: true
    Full_dim: true





If the given polytope is not full-dimensional, its normal fan will have
lineality.

`face_fan <https://polymake.org/release_docs/latest/fan.html#fan__face_fan__44>`__
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Face fans of polytopes are always constructed with respect to a certain
point in the polytopes relative interior. Providing it is optional if
the polytope is centered. Zero will be used as default. If the polytope
is not centered you have to pass such a point as a second argument (in
homogeneous coordinates). For example:


::

    polymake> $v = new Vector([1,0,0,1/2]);
    ........> $ff = face_fan(cross(3), $v);

`k_skeleton <https://polymake.org/release_docs/latest/fan.html#fan__k_skeleton__46>`__
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This client can be used to obtain a subfan consisting of all cones up to
a certain dimension. As an example we construct the skeleton of ``$nf``
with `k=2`:


::

    polymake> $nf2skel = k_skeleton($nf,2);

By taking a look at the f-vectors one can see that the latter has no
cones of dimension 3.


::

    polymake> print "normal fan: ",$nf->F_VECTOR,"\n";
    ........> print "skeleton:   ",$nf2skel->F_VECTOR;
    normal fan: 6 12 8
    skeleton:   6 12




This can also be seen in the Hasse diagram of the skeleton.

Note that the Hasse diagram of a polyhedral fan will always contain an
artifical node at the top which is marked in black and does not
correspond to any cone.


::

    polymake> svg($nf2skel->HASSE_DIAGRAM->VISUAL);
    requires PDFLaTeX and a PDF viewer;
    please specify the output File option or call reconfigure("common::pdfviewer.rules");



