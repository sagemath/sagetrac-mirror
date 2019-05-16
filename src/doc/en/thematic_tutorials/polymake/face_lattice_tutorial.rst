.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Face lattices (of Polytopes)
============================

By definition the face lattice of a polytope contains all the
combinatorial information about a polytope. Here we want to explore how
to work with this in polymake. Let’s start simple.


.. link

.. CODE-BLOCK:: perl

    polymake> $p = n_gon(5);
    polymake> $HD = $p->HASSE_DIAGRAM;       
    polymake> print $HD->FACES;
    {}
    {0}
    {1}
    {2}
    {3}
    {4}
    {1 2}
    {2 3}
    {3 4}
    {0 4}
    {0 1}
    {0 1 2 3 4}





.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package cdd
      cddlib
      Implementation of the double description method of Motzkin et al.
      Copyright by Komei Fukuda.
      http://www-oldurls.inf.ethz.ch/personal/fukudak/cdd_home/
    
    </pre>
    </details>




The obvious question is: How to interpret that output? Well, the Hasse
diagram of a partially ordered set is implemented as a special kind of a
directed graph. Hence all operations on (directed) graphs work. Each
node of the Hasse diagram represents a face. One way to give a name to
such a face is to list all the vertices contained, and this is the
output above. A key feature is that the faces come sorted by dimension.

Very often just a part of the face lattice is interesting. The following
command lists just the 1-dimensional faces.


.. link

.. CODE-BLOCK:: perl

    polymake> print map { $p->HASSE_DIAGRAM->FACES->[$_] } @{$p->HASSE_DIAGRAM->nodes_of_dim(1)};
    {1 2}{2 3}{3 4}{0 4}{0 1}




Face lattices of polytopes can be huge. So it may be an advantage to
compute only a part. The following computes the 2-skeleton of an
8-dimensional cube. We need to specify 3 as a limit since this function
uses the rank which is one more than the dimension.


.. link

.. CODE-BLOCK:: perl

    polymake> $c=cube(8);
    polymake> $HD_partial = lower_hasse_diagram($c->VERTICES_IN_FACETS,3);
    polymake> print $HD_partial->INVERSE_RANK_MAP;
    {(0 (0 0)) (1 (1 256)) (2 (257 1280)) (3 (1281 3072)) (4 (3073 3073))}




Instead of listing all those thousands of faces here we give only the
pairs indicating the start-node and end-node for each rank.


.. link

.. CODE-BLOCK:: perl

    polymake> print map { $HD_partial->nodes_of_rank($_)->size," " } (1..3);
    256 1024 1792 




Dealing with Large Polytopes
----------------------------

In order to get the most out of the above it is important to understand
that this kind of computation is a two-staged process. In the first
step, polymake determines VERTICES_IN_FACETS, that is, the combinatorial
description of the polytope, and then the computation of HASSE_DIAGRAM
only takes this incidence matrix as its input; no coordinates involved
in the second step.

To get at the incidence matrix typically requires a convex hull
computation. This can be triggered like this. Intentionally, there is no
output.


.. link

.. CODE-BLOCK:: perl

    polymake> $p=rand_sphere(6,100); 
    polymake> $p->VERTICES_IN_FACETS;

Notice that this takes a couple of seconds with the default convex hull
code (via cdd), even on a large machine; and the reason is that the
double description method employed is not best possible for this kind of
input. You can speed up *this* computation as follows.


.. link

.. CODE-BLOCK:: perl

    polymake> prefer_now "beneath_beyond"; $p->VERTICES_IN_FACETS;

Temporarily (that is, only in this command) we change the default convex
hull code to polymake’s built-in beneath-and-beyond method. This is
faster on this particular input.

In general, there is no way to tell ahead of time which convex hull
algorithm works best. So, for your own experiments you will have to try.
To get an idea you might want to look up:

-  Avis, David; Bremner, David; Seidel, Raimund: How good are convex
   hull algorithms? 11th ACM Symposium on Computational Geometry
   (Vancouver, BC, 1995). Comput. Geom. 7 (1997), no. 5-6, 265–301.

-  Joswig, Michael: Beneath-and-beyond revisited. Algebra, geometry, and
   software systems, 1–21, Springer, Berlin, 2003.

The subsequent second stage looks as above; but the difference is that
VERTICES_IN_FACETS is known already.


.. link

.. CODE-BLOCK:: perl

    polymake> $HD_partial = lower_hasse_diagram($p->VERTICES_IN_FACETS,3);
    polymake> print map { $HD_partial->nodes_of_rank($_)->size," " } (1..3);
    100 1965 10402 




The executive summary: While polymake is designed to do all kinds of
things automatically, you might have to guide it a little if you are
computing with large or special input.

One more caveat: A *d*-polytope with *n* vertices has at most
*O(n^(d/2))* facets. This is the consequence of the Upper-Bound-Theorem.

-  McMullen, Peter: The maximum numbers of faces of a convex polytope.
   Mathematika 17 (1970) 179-184. This number is actually attained by
   neighborly polytopes; for example, by the cyclic polytopes.
