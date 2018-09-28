.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Getting started with the application ``fan``
============================================

Besides the name giving polyhedral fans this application covers a few
other big objects and related functions. An overview can be found in the
`documentation <https://polymake.org/release_docs/latest/fan.html>`__ or
the `interactive
help <https://polymake.org/doku.php/tutorial/intro_tutorial#getting_help>`__.

However, his tutorial focuses on
`PolyhedralFan <https://polymake.org/release_docs/latest/fan.html#fan__PolyhedralFan__27>`__
objects. Subdivisions have their own tutorial
`here <https://polymake.org/doku.php/tutorial/regular_subdivisions>`__
and some notes on polyhedral complexes can be found
`here <https://polymake.org/doku.php/tutorial/pcom>`__.

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
    polymake: used package cdd
      cddlib
      Implementation of the double description method of Motzkin et al.
      Copyright by Komei Fukuda.
      http://www.ifor.math.ethz.ch/~fukuda/cdd_home/cdd.html
    
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
`documentation </release_docs/latest/fan.html#fan__MAXIMAL_CONES_LINEAR_SPAN_NORMALS__172>`__
for more informations. All maximal cones in ``$f`` are full dimensional,
hence ``LINEAR_SPAN_NORMALS`` is empty and treated like a false:


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
    Regular: 1
    Pure: 1
    Complete: 1
    Full_dim: 1




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
    polymake: used package SVG
       Generated using the Perl SVG Module
       by Ronan Oger




.. raw:: html

    <?xml version="1.0" encoding="UTF-8" standalone="yes"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
    <svg height="841pt" id="document" viewBox="0 -824 570 824" width="595pt" xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
    	<title id="document_title">unnamed</title>
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="297.5" x2="470.561224489796" y1="-66.5" y2="-302.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="297.5" x2="124.438775510204" y1="-66.5" y2="-302.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="297.5" x2="265.776596945832" y1="-66.5" y2="-302.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="297.5" x2="351.581632653061" y1="-66.5" y2="-302.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="297.5" x2="319.132653061225" y1="-66.5" y2="-302.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="297.5" x2="243.418367346939" y1="-66.5" y2="-302.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="470.561224489796" x2="452.083333333333" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="470.561224489796" x2="540.416666666667" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="470.561224489796" x2="496.25" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="470.561224489796" x2="407.916666666667" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="124.438775510204" x2="98.75" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="124.438775510204" x2="187.083333333333" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="124.438775510204" x2="142.916666666667" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="124.438775510204" x2="54.5833333333333" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="265.776596945832" x2="452.083333333333" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="265.776596945832" x2="98.75" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="265.776596945832" x2="275.416666666667" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="265.776596945832" x2="231.25" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="351.581632653061" x2="540.416666666667" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="351.581632653061" x2="187.083333333333" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="351.581632653061" x2="363.75" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="351.581632653061" x2="319.583333333333" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="319.132653061225" x2="496.25" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="319.132653061225" x2="142.916666666667" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="319.132653061225" x2="275.416666666667" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="319.132653061225" x2="363.75" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="243.418367346939" x2="407.916666666667" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="243.418367346939" x2="54.5833333333333" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="243.418367346939" x2="231.25" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="243.418367346939" x2="319.583333333333" y1="-302.5" y2="-538.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="452.083333333333" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="540.416666666667" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="496.25" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="407.916666666667" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="98.75" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="187.083333333333" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="142.916666666667" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="54.5833333333333" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="275.416666666667" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="231.25" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="363.75" x2="297.5" y1="-538.5" y2="-774.5" />
    	<line stroke="rgb(0,0,0)" stroke-width="1" x1="319.583333333333" x2="297.5" y1="-538.5" y2="-774.5" />
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="14.52" x="290.24" y="-74.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="297.5" y="-64"> </text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="9.76" x="465.681224489796" y="-310.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="470.561224489796" y="-300">0</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="9.76" x="119.558775510204" y="-310.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="124.438775510204" y="-300">1</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="9.76" x="260.896596945832" y="-310.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="265.776596945832" y="-300">2</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="9.76" x="346.701632653061" y="-310.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="351.581632653061" y="-300">3</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="9.76" x="314.252653061225" y="-310.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="319.132653061225" y="-300">4</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="9.76" x="238.538367346939" y="-310.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="243.418367346939" y="-300">5</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="442.443333333333" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="452.083333333333" y="-536">0 2</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="530.776666666667" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="540.416666666667" y="-536">0 3</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="486.61" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="496.25" y="-536">0 4</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="398.276666666667" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="407.916666666667" y="-536">0 5</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="89.11" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="98.75" y="-536">1 2</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="177.443333333333" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="187.083333333333" y="-536">1 3</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="133.276666666667" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="142.916666666667" y="-536">1 4</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="44.9433333333333" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="54.5833333333333" y="-536">1 5</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="265.776666666667" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="275.416666666667" y="-536">2 4</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="221.61" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="231.25" y="-536">2 5</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="354.11" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="363.75" y="-536">3 4</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(255,255,255); stroke: rgb(0,0,0); stroke-width: 1" width="19.28" x="309.943333333333" y="-546.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="319.583333333333" y="-536">3 5</text>
    	<rect height="13.5" rx="0" ry="0" style="fill: rgb(0,0,0); stroke: rgb(0,0,0); stroke-width: 1" width="14.52" x="290.24" y="-782.5" />
    	<text font-family="Times-Roman" font-size="10" text-anchor="middle" x="297.5" y="-772"> </text>
    	<!-- 
    	Generated using the Perl SVG Module V2.64
    	by Ronan Oger
    	Info: http://www.roitsystems.com/
     -->
    </svg

