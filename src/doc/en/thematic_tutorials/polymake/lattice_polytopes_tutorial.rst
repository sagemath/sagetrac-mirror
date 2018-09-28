.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Tutorial for Lattice Polytopes
==============================

This page gives a small introduction to lattice polytopes in
``polymake``, some useful external software, and usage hints for it. For
a list of methods and properties applicable to lattice polytopes see
`here <tutorial/lattice_polytopes_doc>`__. For an introduction to the
``polymake`` package see `here <tutorial/start>`__.

``polymake`` always assumes that the lattice used to define a lattice
polytope is the standard lattice Zd. Some rules also require that the
polytope is full dimensional. There are user functions that transform a
polytope sitting in some affine subspace of Rd into a full dimensional
polytope, either in the induced lattice or the lattice spanned by the
vertices, see below.

Dependence on other Software
----------------------------

For some computations ``polymake`` has no built-in commands and passes
the computation to external software. Currently, polymake has an
interface to the following packages that compute various properties of
lattice polytopes.

-  `libnormaliz <http://www.math.uos.de/normaliz/>`__ by Winfried Bruns
   and Bogdan Ichim, bundled with polymake

-  `4ti2 <http://www.4ti2.de/>`__ by the 4ti2 team

-  `LattE macchiato <http://www.math.ucdavis.edu/~mkoeppe/latte/>`__ by
   Matthias Köppe, building on ``LattE`` by Jesus de Loera et. al.

-  (`barvinok <http://freshmeat.net/projects/barvinok>`__ by Sven
   Verdoolaege)

Unless you want to deal with Hilbert bases of cones you don’t need them.
If you do, either the bundled extension ``libnormaliz`` or the external
package ``4ti2`` suffices to do most computations with lattice
polytopes. Computation of Gröbner bases currently requires ``4ti2``.
``LattE`` only counts lattice points in a polytope and computes its
Ehrhart polynomial, but may be faster on that than any other methods
implemented. ``barvinok`` can be used to compute the number of lattice
points and the h-polynomial. Access to barvinok is realized via an
extension which has to be downloaded separately.

For some of the commands in this tutorial you will need at least one of
``bundled:libnormaliz`` enabled or ``4ti2`` installed on your machine.
We’ll remind you at the relevant places.

Lattice Points in Rational Polytopes
------------------------------------

We start by creating a rational polytope using one of ``polymake``\ ’s
standard polytope constructions. We choose the 3-dimensional cube with
coordinates +1 and -1. So we start ``polymake`` at the command line and
assign a cube to the variable $p.


::

    polymake> $p=cube(3);

Suppose we want to know how many lattice points this cube contains. The
answer is of course already known, as the cube has one relative interior
integral point per non-empty face. So we expect to get the answer 27.


::

    polymake> print $p->N_LATTICE_POINTS;
    27




.. raw:: html

    <details>
    <summary>
    <pre><small>Click here for additional output</small></pre>
    </summary>
    <pre>
    polymake: used package latte
      LattE (Lattice point Enumeration) is a computer software dedicated to the 
      problems of counting lattice points and integration inside convex polytopes.
      Copyright by Matthias Koeppe, Jesus A. De Loera and others.
      http://www.math.ucdavis.edu/~latte/
    
    </pre>
    </details>




To satisfy this request, ``polymake`` computes all properties necessary
to call an external program that provides the number of lattice points.
In this case, ``polymake`` has passed the request to ``lattE``, which is
shown by the credit message that appears before the answer. By default,
credits for external software are shown when an external package is used
for the first time. You can change this behavior using the variable
``$Verbose::credits``. If you don’t have a version of ``LattE``, or if
you have set different preferences, then ``polymake`` may choose one of
the other programs. So the credit statement depends on your
configuration.

We can of course also ask ``polymake`` to compute the integral points
for us. For our next computations we are only interested in the integral
points in the interior of the cube, so we ask for


::

    polymake> print $p->INTERIOR_LATTICE_POINTS;
    1 0 0 0





.. raw:: html

    <details>
    <summary>
    <pre><small>Click here for additional output</small></pre>
    </summary>
    <pre>
    polymake: used package libnormaliz
      Normaliz is a tool for computations in affine monoids, vector configurations, lattice polytopes, and rational cones.
      Copyright by Winfried Bruns, Bogdan Ichim, Christof Soeger.
      http://www.math.uos.de/normaliz/
    
    </pre>
    </details>




Internally, ``polymake`` computes the intersection of the polytope with
the integer lattice, and then checks which of the points lies on a facet
of $p. By default, ``polymake`` uses a project-and-lift algorithms to
enumerate the lattice points. Note that our call to ``LattE`` above has
only computed the number of integral points (which is done with an
improved version of Barvinok’s algorithm), so ``polymake`` really has to
compute something here. If we had asked for ``INTERIOR_LATTICE_POINTS``
first, then ``N_LATTICE_POINTS`` would just have counted the rows of a
matrix, which would have been much faster. So computation time can
depend on the history.

You can also ask for the HILBERT_BASIS, though in the case of a cube the
result is not so exciting:


::

    polymake> print $p->HILBERT_BASIS;
    1 -1 -1 -1
    1 -1 -1 0
    1 -1 -1 1
    1 -1 0 -1
    1 -1 0 0
    1 -1 0 1
    1 -1 1 -1
    1 -1 1 0
    1 -1 1 1
    1 0 -1 -1
    1 0 -1 0
    1 0 -1 1
    1 0 0 -1
    1 0 0 0
    1 0 0 1
    1 0 1 -1
    1 0 1 0
    1 0 1 1
    1 1 -1 -1
    1 1 -1 0
    1 1 -1 1
    1 1 0 -1
    1 1 0 0
    1 1 0 1
    1 1 1 -1
    1 1 1 0
    1 1 1 1





``polymake`` has no native method to compute a Hilbert basis, so it has
passed the computation to ``4ti2``. The choice may vary, depending on
what is installed on your computer (and configured for ``polymake``).
You can influence the choice with the appropriate ``prefer`` statement.

Note that so far these commands also work for rational polytopes.

Lattice Polytopes
-----------------

Now we want to do some computations that don’t make sense for polytopes
that have non-integral vertex coordinates. We can let ``polymake`` check
that our cube is indeed a polytope with integral vertices.


::

    polymake> print $p->LATTICE;
    1




A particularly interesting class of lattice polytopes is that of
reflexive polytopes. A polytope is *reflexive* if its polar is agein
alattice polytope. This implies in particular that the origin is the
unique interior lattice point in the polytope. So, as we have seen
above, our cube is a candidate. But this is not sufficient, so we have
to do further checks.

Reflexivity is a property that is not defined for polytopes with
non-integral vertices. So if we ask for it in ``polymake``, then
``polymake`` checks that the entered polytope is indeed a lattice
polytope (i.e. it is **bounded** and has **integral vertices**). In that
case the object will automatically get the specialization
``Polytope::Lattice``.


::

    polymake> print $p->REFLEXIVE;
    1




Lattice polytopes can be used to define toric varieties with an ample
line bundle, and many properties of the variety are reflected by the
polytope. here is an example: The toric variety defined by our cube is
*smooth*, i.e. it is one of the *smooth toric Fano varieties*. In
``polymake``, we can just ask for this property in the following way.


::

    polymake> print $p->SMOOTH;
    1




The number of integral points in the k-th dilate of a polytope is given
by a polynomial of degree d in k. This is the famous *Ehrhart Theorem*.
In ``polymake`` you can obtain the coefficients of this polynomial
(starting with the constant coefficient).


::

    polymake> print $p->EHRHART_POLYNOMIAL_COEFF;
    1 6 12 8




``polymake`` has passed this request to ``LattE`` or ``normaliz``, but
as we have used these programs already the credit message is suppressed
(but if you save the cube to a file, then you will find it in there).
Some coefficients of this polynomial have a geometric interpretation.
E.g., the highest coefficient is the Euclidean volume of the polytope.


::

    polymake> print $p->VOLUME;
    8




By a theorem of Stanley, the generating function for the number of
lattice points can be written as the quotient of a polynomial h(t) by
(1-t)d+1, and this polynomial has non-negative integral coefficients.


::

    polymake> print $p->H_STAR_VECTOR;
    1 23 23 1




::

    polymake> print $p->LATTICE_DEGREE;
    3




::

    polymake> print $p->LATTICE_CODEGREE;
    1




In our case the coefficient vector is symmetric, as the polytope is
reflexive. The *co-degree* of the polytope is d+1 minus the degree of
the h-polynomial. It is the smallest factor by which we have to dilate
the polytope to obtain an interior integral point. In our case, this is
1, as the cube already has an integral point.

We can obtain the volume of our polytope also from the
``H_STAR_VECTOR``: Summing up the coefficients give the *lattice volume*
of the polytope, which is d! times its Euclidean volume.


::

    polymake> print $p->LATTICE_VOLUME;
    48




Let us look at a different example:


::

    polymake> $q=new Polytope(INEQUALITIES=>[[5,-4,0,1],[-3,0,-4,1],[-2,1,0,0],[-4,4,4,-1],[0,0,1,0],[8,0,0,-1],[1,0,-1,0],[3,-1,0,0]]);

This actually defines a lattice polytope, which we can see from the list
of vertices:


::

    polymake> print $q->VERTICES;
    1 3 0 8
    1 2 1 8
    1 3 1 8
    1 2 0 3
    1 3 0 7
    1 3 1 7
    1 2 1 7
    1 2 0 4





.. raw:: html

    <details>
    <summary>
    <pre><small>Click here for additional output</small></pre>
    </summary>
    <pre>
    polymake: used package ppl
      The Parma Polyhedra Library (PPL): A C++ library for convex polyhedra
      and other numerical abstractions.
      http://www.cs.unipr.it/ppl/
    
    </pre>
    </details>




``polymake`` provides basically three methods for convex hull
conversion, double description, reverse search, and beneath beyond. The
first two are provided by the packages ``cdd`` and ``lrs``, the last in
internal. By default, ``cdd`` is chosen, and that is what was used above
(they are bundled with ``polymake``, you don’t have to install them). A
polytope Q is *normal* if every lattice point in the k-th dilate of Q is
the sum of k lattice points in Q. You can check this property via


::

    polymake> print $q->NORMAL;

So our polytope is not normal. We can also find a point that violates
the condition. Being normal is equivalent to the fact, that the Hilbert
basis of the cone C(Q) obtained from Q by embedding the polytope at
height one and the coning over it has all its generators in height one.
The property HILBERT_BASIS computes these generators:


::

    polymake> print $q->HILBERT_BASIS;
    1 2 0 3
    1 2 0 4
    1 2 1 7
    1 2 1 8
    1 3 0 7
    1 3 0 8
    1 3 1 7
    1 3 1 8
    2 5 1 13





The last row is the desired vector: [2,5,1,13] is a vector in 2*Q, but
it is not a sum of lattice points in Q. The cone C(Q) corresponds to an
affine toric variety, and the above tells us that this variety is not
normal. Yet, it is very ample, as we can check with


::

    polymake> print $q->VERY_AMPLE;
    1




.. raw:: html

    <details>
    <summary>
    <pre><small>Click here for additional output</small></pre>
    </summary>
    <pre>
    polymake: used package cdd
      cddlib
      Implementation of the double description method of Motzkin et al.
      Copyright by Komei Fukuda.
      http://www-oldurls.inf.ethz.ch/personal/fukudak/cdd_home/
    
    </pre>
    </details>




Now assume we are particularly interested in the third facet of Q. We
can pick this via


::

    polymake> $f=facet($q,2);

Recall that indexes in ``polymake`` start at 0, so the third facet has
index 2. This is again a very ample polytope:


::

    polymake> print $f->VERY_AMPLE;
    1




The result is no surprise, being very ample is inherited by faces. We
could also be interested in the facet width of the polytope $f. This is
the minimum over the maximal distance of a facet to any other vertex.
``polymake`` knows how to compute this:


::

    polymake> print $f->FACET_WIDTH;


.. raw:: html

    <details>
    <summary>
    <pre><small>Click here for additional output</small></pre>
    </summary>
    <pre>
    polymake:  WARNING: could not compute 'FACET_WIDTH' probably because of unsatisfied preconditions:
    precondition : FULL_DIM ( FACET_WIDTHS : VERTICES , FACETS , CONE_AMBIENT_DIM )
    </pre>
    </details>




Almost. It tells you that it can only do this for a full dimensional
polytope, i.e. for a polytope whose dimension coincides with the ambient
dimension. This is not true for our facet: It lives in the same ambient
space as $q, but has one dimension less. We can remedy this by applying
the following:


::

    polymake> $g=ambient_lattice_normalization($f);
    ........> print $g->FACET_WIDTH;
    1




The function ``ambient_lattice_normalization`` returns a full
dimensional version of the polytope ``$f`` in the lattice induced by the
intersection of the affine space of ``$f`` with Z^n. Now ``$g`` is full
dimensional, and we can compute the facet width. Note that there is also
a function which normalizes in the lattice spanned by the vertices of
the polytope: ``vertex_facet_normalization``. This can also be usefull
for full dimensional polytopes. E.g. consider the cube we defined above.
The sum of the entries of each vertex is odd, so the lattice spannd by
the vertices is a sublattice of the integer lattice:


::

    polymake> $cr=vertex_lattice_normalization($p);
    ........> print $cr->VERTICES;
    (4) (0 1)
    1 1 0 0
    1 0 1 0
    1 1 1 0
    1 0 0 1
    1 1 0 1
    1 0 1 1
    1 1 1 1





``$cr`` is the same cube, but we have reduced the lattice. (The first
line is a *sparse representation* of a vector: it has length 4, and the
only non-zero entry is at position 0 and is 1 (note that indexes start
at 0)).

Toric Varieties
---------------

``polymake`` has only few builtin functions to compute properties of the
variety associated to a fan or lattice polytope. There are two
extensions available that add more properties, both currently at an
early stage:

-  `Toric Varieties and Singular
   interface <https///github.com/lkastner>`__ by Lars Kastner/Benjamin
   Lorenz

-  `ToricVarieties-v0.3 <http://www.mathematik.tu-darmstadt.de/~paffenholz/software.html>`__
   by Andreas Paffenholz. Defines a new property for toric varieties
   associated to a fan and divisors on that variety.

Here we will do some computations that do not require one of the
extensions. We start by defining a fan. We’ll make our live easy and
take the normal fan of our cube:


::

    polymake> application "fan";

::

    polymake> $f = normal_fan($p);
    ........> print $f->SMOOTH_FAN;
    1




With the last line we have verified that our fan defines a smooth toric
variety. Note that switching the application is not strictly necessary,
you can also prepend calls to functions and constructors with ``fan::``.
The fan object $f itself knows its type, and chooses available
properties based on this. Any smooth variety is Gorenstein, so we expect
the following:


::

    polymake> print $f->GORENSTEIN;
    1




Similarly, we could check for Q-Gorensteinness with ``Q_GORENSTEIN``. It
is also a complete fan:


::

    polymake> print $f->COMPLETE;
    1




but currently there is little support to detect completeness in
``polymake``. In our case it was already decided during construction,
normal fans are complete. You can also check standard features of fans,
like their rays. Let us do this for the normal fan of our other example:


::

    polymake> $g=normal_fan($q);
    ........> print $g->RAYS;
    -1 0 1/4
    0 -1 1/4
    1 0 0
    1 1 -1/4
    0 1 0
    0 0 -1
    0 -1 0
    -1 0 0





This is not what we wanted. We would like to see the minimal lattice
generators of the rays. We can fix this using


::

    polymake> print primitive($g->RAYS);
    -4 0 1
    0 -4 1
    1 0 0
    4 4 -1
    0 1 0
    0 0 -1
    0 -1 0
    -1 0 0





Note that the function ``primitive`` returns a copy of the argument, the
RAYS as stored in the fan are unchanged. So you have to apply this
function each time you need the primitive generators, or you store them
in a new variable. The fan `g` is not smooth, but still
Gorenstein:


::

    polymake> print $g->SMOOTH_FAN;

::

    polymake> print $g->GORENSTEIN;
    1




You can also access the maximal cones of the fan via


::

    polymake> print $g->MAXIMAL_CONES;
    {3 4 5 7}
    {2 3 5 6}
    {5 6 7}
    {0 1 2 4}
    {0 4 7}
    {0 1 6 7}
    {1 2 6}
    {2 3 4}





The indices in these list refer to the list of rays. Sometimes you might
be interested in the walls, i.e. the codimension 2 faces of the fan.
Here is one way to get them


::

    polymake> print rows_numbered($g->HASSE_DIAGRAM->FACES);
    0:-1
    1:3 4 5 7
    2:2 3 5 6
    3:5 6 7
    4:0 1 2 4
    5:0 4 7
    6:0 1 6 7
    7:1 2 6
    8:2 3 4
    9:3 5
    10:5 7
    11:4 7
    12:3 4
    13:5 6
    14:2 6
    15:2 3
    16:6 7
    17:0 4
    18:0 1
    19:1 2
    20:2 4
    21:0 7
    22:1 6
    23:3
    24:5
    25:7
    26:4
    27:6
    28:2
    29:0
    30:1
    31:





::

    polymake> print $g->HASSE_DIAGRAM->nodes_of_dim($g->DIM-2);
    {23 24 25 26 27 28 29 30}




where the list of numbers given by the latter are the indices of the
codimension 2 faces in the list of all faces given before. There is a
more concise way to list those, using some simple perl programming:


::

    polymake> print map($g->HASSE_DIAGRAM->FACES->[$_], @{$g->HASSE_DIAGRAM->nodes_of_dim($g->DIM-2)});
    {3}{5}{7}{4}{6}{2}{0}{1}




Visualization
-------------

If the lattice polytope lives in R^2 or R^3, then we can visualize the
polytope together with its lattice points. The picture below has been
made with ```javaview`` <http://www.javaview.de>`__, but
``polymake``\ ’s standard visualization method is now ``jreality``,
which is bundled with ``polymake``.


::

    polymake> $p->VISUAL->LATTICE_COLORED;


.. raw:: html

    <?xml version="1.0" encoding="UTF-8" standalone="yes"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
    <svg height="841pt" id="document" viewBox="0 -615 575 615" width="595pt" xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
    	<title id="document_title">p</title>
    	<g>
    		<polygon points="35,-578.5 35,-53.5 560,-53.5 560,-578.5 " style="fill: rgb(); fill-opacity: 1; stroke: rgb(0,0,0); stroke-width: 1" />
    	</g>
    	<g />
    	<g>
    		<polygon points="35,-53.5 35,-53.5 35,-578.5 35,-578.5 " style="fill: rgb(); fill-opacity: 1; stroke: rgb(0,0,0); stroke-width: 1" />
    		<polygon points="560,-578.5 560,-53.5 560,-53.5 560,-578.5 " style="fill: rgb(); fill-opacity: 1; stroke: rgb(0,0,0); stroke-width: 1" />
    		<polygon points="560,-53.5 35,-53.5 35,-53.5 560,-53.5 " style="fill: rgb(); fill-opacity: 1; stroke: rgb(0,0,0); stroke-width: 1" />
    		<polygon points="35,-578.5 35,-578.5 560,-578.5 560,-578.5 " style="fill: rgb(); fill-opacity: 1; stroke: rgb(0,0,0); stroke-width: 1" />
    		<polygon points="35,-53.5 35,-578.5 560,-578.5 560,-53.5 " style="fill: rgb(); fill-opacity: 1; stroke: rgb(0,0,0); stroke-width: 1" />
    	</g>
    	<g />
    	<circle cx="297.5" cy="-316" r="2" style="fill: rgb(30,250,30)" />
    	<circle cx="35" cy="-53.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="35" cy="-53.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="35" cy="-53.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="35" cy="-316" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="35" cy="-316" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="35" cy="-316" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="35" cy="-578.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="35" cy="-578.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="35" cy="-578.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="297.5" cy="-53.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="297.5" cy="-53.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="297.5" cy="-53.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="297.5" cy="-316" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="297.5" cy="-316" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="297.5" cy="-578.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="297.5" cy="-578.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="297.5" cy="-578.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="560" cy="-53.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="560" cy="-53.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="560" cy="-53.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="560" cy="-316" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="560" cy="-316" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="560" cy="-316" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="560" cy="-578.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="560" cy="-578.5" r="2" style="fill: rgb(70,150,70)" />
    	<circle cx="560" cy="-578.5" r="2" style="fill: rgb(70,150,70)" />
    	<!--
    	Generated using the Perl SVG Module V2.82
    	by Ronan Oger
    	Info: http://www.roitsystems.com/
    	-->
    </svg



.. raw:: html

    <details>
    <summary>
    <pre><small>Click here for additional output</small></pre>
    </summary>
    <pre>
    polymake: used package SVG
       Generated using the Perl SVG Module
       by Ronan Oger
    
    </pre>
    </details>




The output may look similar to the following.

.. figure:: attachment:cube-in-lattice.gif
   :alt: {{:tutorial:cube-in-lattice.gif?298}}

   {{:tutorial:cube-in-lattice.gif?298}}

The command ``LATTICE_COLORED`` sorted the lattice points into three
classes before visualization: lattice points in the interior of the
polytope, lattice points on the boundary, and vertices that are not in
the lattice. These classes are then visualized with different colors
(where we only see two in the above picture, as all vertices of the cube
are in the lattice). If you don’t need this distinction,
``VISUAL->LATTICE`` avoids the additional computations.

External Packages
-----------------

``polymake`` can use ``4ti2`` and ``lattE`` via a file based interface
and ``libnormaliz >= 3.1.0`` as library, (the file based interface to
``normaliz`` has been discontinued) for lattice computations and prints
all available packages during startup. To tell ``polymake`` about a
newly installed program run ``%%polymake --reconfigure%%`` or issue the
command ``reconfigure`` during the interactive session. polymake may ask
you to confirm the paths to the binaries.

::

   Application polytope uses following third-party software (for details: help 'credits';)
   4ti2, cddlib, latte, libnormaliz, lrslib, nauty

The output at this position depends on the software available on your
computer. To see each call to an external program you can either set the
variable ``$Verbose::external=1;`` on the command line or include the
line ``$Polymake::User::Verbose::external=1`` in
``~/.polymake/prefer.pl``. If you just want to see the credit message
instead of the program call, set ``$Verbose::credits=2`` instead. If
this is 1, then a credit is shown when a package is used for the first
time, if 0, then all credits are suppressed (but you can find them in
the file afterwards).


::

    polymake> $Verbose::external=1;

::

    polymake> print $p->EHRHART_POLYNOMIAL_COEFF;
    1 6 12 8




You can ask ``polymake`` to prefer one package over another by setting
``prefer "program";`` where program is one of ``_4ti2``, ``latte`` and
``normaliz2``. Of course, the corresponding package needs to be
installed on your computer.

To prefer one program only for some computations you may append one of
.integer_points, .hilbert, .ehrhartpoly for rules computing
N_LATTICE_POINTS, LATTICE_POINTS, HILBERT_BASIS or
EHRHART_POLYNOMIAL_COEFF. (Or ``prefer_now`` just for the next
computation)


::

    polymake> print cube(2)->N_LATTICE_POINTS;
    9




.. raw:: html

    <details>
    <summary>
    <pre><small>Click here for additional output</small></pre>
    </summary>
    <pre>
    polymake: running latte's count: cd /tmp/poly1567Taaaa0002; count  --ehrhart-polynomial input.ine
    </pre>
    </details>




::

    polymake> prefer_now "libnormaliz";
    ........> print cube(2)->N_LATTICE_POINTS;
    9




::

    polymake> print cube(2)->EHRHART_POLYNOMIAL_COEFF;
    1 4 4




.. raw:: html

    <details>
    <summary>
    <pre><small>Click here for additional output</small></pre>
    </summary>
    <pre>
    polymake: running latte's count: cd /tmp/poly1567Taaaa0003; count  --ehrhart-polynomial input.ine
    </pre>
    </details>


