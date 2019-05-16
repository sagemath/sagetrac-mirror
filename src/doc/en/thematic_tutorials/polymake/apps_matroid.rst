.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Application ``matroid``: Short Introduction
===========================================

This tutorial is meant to show the main features for handling matroids
available. To make ``matroid`` your current application start
``polymake`` with the option ``-A matroid`` or use the context switch


.. link

.. CODE-BLOCK:: perl

    polymake> application "matroid";

from within the ``polymake`` shell. A permanent setting can be stored
with ``%%set_custom $default_application="matroid";%%``

In Sage the switching is done as follows::

  sage: polymake.application("matroid")

Constructing a Simple Matroid and Playing Around
------------------------------------------------

This is how to produce a matroid from a vector configuration over the
rationals. The matroid is defined by the linear dependence among subsets
of these vectors.


.. link

.. CODE-BLOCK:: perl

    polymake> $M=new Matroid(VECTORS=>[[1,0,0],[1,0,1],[1,1,0],[1,0,2]]);

.. link

In Sage::

    sage: M = polymake.new_object("Matroid", VECTORS=[[1,0,0],[1,0,1],[1,1,0],[1,0,2]])

If ``matroid`` is not your default application you have to qualify
``Matroid`` as in:

.. link

.. CODE-BLOCK:: perl

    polymake> $M=new matroid::Matroid(VECTORS=>[[1,0,0],[1,0,1],[1,1,0],[1,0,2]]);

.. link

In Sage::

    sage: M = polymake.new_object("matroid::Matroid", VECTORS=[[1,0,0],[1,0,1],[1,1,0],[1,0,2]])

Output of basic statistics.


.. link

.. CODE-BLOCK:: perl

    polymake> print $M->N_BASES, " ", $M->N_ELEMENTS, " ", $M->RANK;
    3 4 3

.. link

In Sage::

    sage: M.N_BASES, M.N_ELEMENTS, M.RANK
    (3, 4, 3)

|{{ :tutorial:matroid_lattice_of_flats_example.png?nolink&200|}}| The
``VECTORS`` are numbered consecutively, starting from zero. The bases
are encoded as sets of these ordinal numbers.

.. |{{ :tutorial:matroid_lattice_of_flats_example.png?nolink&200|}}| image:: attachment:matroid_lattice_of_flats_example.png


.. link

.. CODE-BLOCK:: perl

    polymake> print $M->BASES;
    {0 1 2}
    {0 2 3}
    {1 2 3}

.. link

In Sage::

    sage: M.BASES
    {0 1 2}
    {0 2 3}
    {1 2 3}

Similarly you can compute the circuits and cocircuits.


.. link

.. CODE-BLOCK:: perl

    polymake> print $M->CIRCUITS;
    {0 1 3}

.. link

In Sage::

    sage: M.CIRCUITS
    {0 1 3}

.. link

.. CODE-BLOCK:: perl

    polymake> print $M->COCIRCUITS;
    {2}
    {1 3}
    {0 3}
    {0 1}

.. link

In Sage::

    sage: M.COCIRCUITS
    {2}
    {1 3}
    {0 3}
    {0 1}

You can also compute other properties, like


.. link

.. CODE-BLOCK:: perl

    polymake> print $M->PAVING?"1":"0", " ",
    ........> $M->BINARY?"1":"0", " ",
    ........> $M->SERIES_PARALLEL?"1":"0", " ",
    ........> $M->CONNECTED?"1":"0";
    1 1 0 0

.. link

In Sage::

    sage: M.PAVING, M.BINARY, M.SERIES_PARALLEL, M.CONNECTED
    (true, true, false, false)

.. link

.. CODE-BLOCK:: perl

    polymake> print $M->CONNECTED_COMPONENTS;
    {0 1 3}
    {2}

.. link

In Sage::

    sage: M.CONNECTED_COMPONENTS
    {0 1 3}
    {2}

.. link

.. CODE-BLOCK:: perl

    polymake> print $M->TUTTE_POLYNOMIAL;
    x_0^3 + x_0^2 + x_0*x_1

.. link

In Sage::

    sage: M.TUTTE_POLYNOMIAL
    x_0^3 + x_0^2 + x_0*x_1

Even the lattice of flats could be computed and visualised.


.. link

.. CODE-BLOCK:: perl

    polymake> $lattice=$M->LATTICE_OF_FLATS;
    polymake> foreach (@{$lattice->nodes_of_rank(2)}){print $lattice->FACES->[$_]," "};
    {0 2} {0 1 3} {1 2} {2 3} 

.. link

In Sage::

    sage: lattice = M.LATTICE_OF_FLATS
    sage: [ lattice.FACES[i] for i in lattice.nodes_of_rank(2) ]
    [{0 2}, {0 1 3}, {1 2}, {2 3}]

.. link

.. CODE-BLOCK:: perl

    polymake> print $M->MATROID_HYPERPLANES;
    {0 1 3}
    {0 2}
    {1 2}
    {2 3}

.. link

In Sage::

    sage: M.MATROID_HYPERPLANES
    {0 1 3}
    {0 2}
    {1 2}
    {2 3}

.. link

.. CODE-BLOCK:: perl

    polymake> $M->LATTICE_OF_FLATS->VISUAL;

.. link

In Sage::

    sage: M.LATTICE_OF_FLATS.VISUAL()             # not tested


Matroid Polytopes
-----------------

You can construct a polytope from the bases of a matroid as the convex
hull of the characteristic vectors of the bases. This is the *matroid
polytope* of that matroid, sometimes also called the *matroid bases
polytope*. The matroid polytope of the matroid ``$M`` is a subobject
``POLYTOPE`` of type ``polytope::Polytope``.


.. link

.. CODE-BLOCK:: perl

    polymake> print $M->POLYTOPE->VERTICES;
    1 1 1 1 0
    1 1 0 1 1
    1 0 1 1 1

.. link

In Sage::

    sage: M.POLYTOPE.VERTICES
    1 1 1 1 0
    1 1 0 1 1
    1 0 1 1 1

.. link

.. CODE-BLOCK:: perl

    polymake> print $M->POLYTOPE->F_VECTOR;
    3 3

.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package lrs
      Implementation of the reverse search algorithm of Avis and Fukuda.
      Copyright by David Avis.
      http://cgm.cs.mcgill.ca/~avis/C/lrs.html
    
    </pre>
    </details>

.. link

In Sage::

    sage: M.POLYTOPE.F_VECTOR
    3 3


Other Constructions
-------------------

The vertices of a polytope give rise to a matroid. Here is an example
for the vertices of the three-dimensional regular cube. Notice that
point coordinates in the application ‘polytope’ are given by homogeneous
coordinates. Hence this matroid is defined by the relation of affine
dependence.


.. link

.. CODE-BLOCK:: perl

    polymake> $C=new Matroid(VECTORS=>polytope::cube(3)->VERTICES);

.. link

In Sage::

    sage: polymake.application("polytope")
    sage: C = polymake.new_object("matroid::Matroid", VECTORS=polymake.cube(3).VERTICES)

.. link

.. CODE-BLOCK:: perl

    polymake> print $C->N_BASES;
    58

.. link

In Sage::

    sage: C.N_BASES
    58

The system also allows you to construct a matroid from a graph. The
bases correspond to the spanning trees then. Notice that there is more
than one way to encode a graph in ``polymake``. Read the `tutorial on
graphs <apps_graph>`__ for details.


.. link

.. CODE-BLOCK:: perl

    polymake> $G=matroid_from_graph(polytope::cube(3)->GRAPH);

.. link

In Sage::

    sage: polymake.application("polytope")
    sage: c = polymake.cube(3)
    sage: polymake.application("matroid")
    sage: G = polymake.matroid_from_graph(c.GRAPH)

.. link

.. CODE-BLOCK:: perl

    polymake> print $G->N_BASES;
    384

.. link

In Sage::

    sage: G.N_BASES
    384

It is also possible to derive a new matroid from others.


.. link

.. CODE-BLOCK:: perl

    polymake> # The arguments are two matroids and for each matroid a basepoint. The basepoints will be identified. 
    polymake> $se=series_extension(uniform_matroid(2,3),0,uniform_matroid(1,3),0);

.. link

In Sage::

    sage: u23 = polymake.uniform_matroid(2,3)
    sage: u13 = polymake.uniform_matroid(1,3)
    sage: se = polymake.series_extension(u23, 0, u13, 0)

.. link

.. CODE-BLOCK:: perl

    polymake> print deletion($se,4)->VECTORS;
    1 0 0
    0 1 0
    0 0 1
    1 1 1

.. link

In Sage::

    sage: polymake.deletion(se, 4).VECTORS
    1 0 0
    0 1 0
    0 0 1
    1 1 1

.. link

.. CODE-BLOCK:: perl

    polymake> $pe=parallel_extension(uniform_matroid(1,3),0,uniform_matroid(2,3),0);

.. link

In Sage::

    sage: pe = polymake.parallel_extension(u13, 0, u23, 0)

.. link

.. CODE-BLOCK:: perl

    polymake> print dual(contraction($pe,4))->VECTORS;
    1 1 1
    1 0 0
    0 1 0
    0 0 1

.. link

In Sage::

    sage: polymake.dual(polymake.contraction(pe, 4)).VECTORS
    1 1 1
    1 0 0
    0 1 0
    0 0 1

.. link

.. CODE-BLOCK:: perl

    polymake> print projective_plane(3)->N_BASES;
    234

.. link

In Sage::

    sage: polymake.projective_plane(3).N_BASES
    234

.. link

.. CODE-BLOCK:: perl

    polymake> print fano_matroid()->N_BASES;
    28

.. link

In Sage::

    sage: polymake.fano_matroid().N_BASES
    28

.. link

.. CODE-BLOCK:: perl

    polymake> print direct_sum(projective_plane(3),fano_matroid())->N_BASES," = 234*28";
    6552 = 234*28

.. link

In Sage::

   sage: polymake.direct_sum(polymake.projective_plane(3), polymake.fano_matroid()).N_BASES
   6552
   sage: _ == 234*28
   True

.. link

.. CODE-BLOCK:: perl

    polymake> print two_sum(uniform_matroid(2,4),0,uniform_matroid(2,4),0)->CIRCUITS;
    {0 1 2}
    {3 4 5}
    {0 1 3 4}
    {0 1 3 5}
    {0 1 4 5}
    {0 2 3 4}
    {0 2 3 5}
    {0 2 4 5}
    {1 2 3 4}
    {1 2 3 5}
    {1 2 4 5}

.. link

In Sage::

    sage: polymake.two_sum(polymake.uniform_matroid(2,4), 0, polymake.uniform_matroid(2,4), 0).CIRCUITS
    {0 1 2}
    {3 4 5}
    {0 1 3 4}
    {0 1 3 5}
    {0 1 4 5}
    {0 2 3 4}
    {0 2 3 5}
    {0 2 4 5}
    {1 2 3 4}
    {1 2 3 5}
    {1 2 4 5}

Of course you can also construct your matroid from scratch by
specifying, e.g., its set of bases or non-bases and then compute other
properties. The following constructs the Fano matroid, which is the
simplest matroid that cannot be constructed from a vector configuration
(over a field with a characteristic other than two).


.. CODE-BLOCK:: perl

    polymake> $a=new Array<Set<Int>>([0,1,5],[1,2,6],[0,2,3],[1,3,4],[2,4,5],[3,5,6],[0,4,6]);
    polymake> $m=new Matroid(NON_BASES=>$a,N_ELEMENTS=>7);
    polymake> print $m->COCIRCUITS;
    {0 1 2 4}
    {0 1 3 6}
    {0 2 5 6}
    {0 3 4 5}
    {1 2 3 5}
    {1 4 5 6}
    {2 3 4 6}

In Sage::

    sage: a = polymake.new_object("Array<Set<Int>>", [0,1,5],[1,2,6],[0,2,3],[1,3,4],[2,4,5],[3,5,6],[0,4,6])
    sage: m = polymake.new_object("matroid::Matroid", NON_BASES=a, N_ELEMENTS=7)
    sage: m.COCIRCUITS
    {0 1 2 4}
    {0 1 3 6}
    {0 2 5 6}
    {0 3 4 5}
    {1 2 3 5}
    {1 4 5 6}
    {2 3 4 6}

Note that you have to specify N_ELEMENTS when constructing a matroid in
this way because this is not implicit in BASES, etc.
