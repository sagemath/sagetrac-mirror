.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Voronoi Diagrams
================

Voronoi diagrams are constructed from their sites (given in homogeneous
coordinates).


::

    polymake> $VD = new VoronoiPolyhedron(SITES=>[[1,1,1],[1,0,1],[1,-1,1],[1,1,-1],[1,0,-1],[1,-1,-1]]);
    ........> $VD->VISUAL_VORONOI;

Actually, via lifting to the standard paraboloid, Voronoi diagrams are
derived from ``Polytope``. That’s why they have ``VERTICES``,
``FACETS``, and such.


::

    polymake> print $VD->FACETS;
    2 -2 -2 1
    1 0 -2 1
    2 2 -2 1
    2 -2 2 1
    1 0 2 1
    2 2 2 1
    1 0 0 0
        





::

    polymake> print $VD->VERTICES;
    0 0 1 2
    0 1 0 2
    1 1/2 0 -1
    0 -1 0 2
    0 0 -1 2
    1 -1/2 0 -1
    






