.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


A Counter-example to an integer analog to Caratheodory’s Theorem
----------------------------------------------------------------

The construction
~~~~~~~~~~~~~~~~

This tutorial describes the construction of a specific rational cone in
six dimensions which is due to:

-  Bruns, Winfried; Gubeladze, Joseph; Henk, Martin; Martin, Alexander;
   Weismantel, Robert: A counterexample to an integer analogue of
   Carathéodory’s theorem. J. Reine Angew. Math. 510 (1999), 179-185.

The rows of this matrix describe a cone *C*:


.. link

.. CODE-BLOCK:: perl

    polymake> $M = new Matrix<Rational>([[0,1,0,0,0,0],
    polymake> [0,0,1,0,0,0],
    polymake> [0,0,0,1,0,0],
    polymake> [0,0,0,0,1,0],
    polymake> [0,0,0,0,0,1],
    polymake> [1,0,2,1,1,2],
    polymake> [1,2,0,2,1,1],
    polymake> [1,1,2,0,2,1],
    polymake> [1,1,1,2,0,2],
    polymake> [1,2,1,1,2,0]]);
    polymake> $C=new Polytope<Rational>(POINTS=>$M);

From


.. link

.. CODE-BLOCK:: perl

    polymake> print $C->HILBERT_BASIS;
    0 0 0 0 0 1
    0 0 0 0 1 0
    0 0 0 1 0 0
    0 0 1 0 0 0
    1 0 2 1 1 2
    0 1 0 0 0 0
    1 1 1 2 0 2
    1 1 2 0 2 1
    1 2 0 2 1 1
    1 2 1 1 2 0
    





one can see that the given generators of *C* form a Hilbert basis. Now
we consider one particular point *x*. The output of the second command
(all coefficients positive) shows that *x* is contained in the interior
of *C*.


.. link

.. CODE-BLOCK:: perl

    polymake> $x=new Vector<Rational>([9,13,13,13,13,13]);
    polymake> print $C->FACETS * $x;
    8 15 19/2 19/2 17 13 17 13 9 13 13 17 8 19/2 13 17 15 19/2 15 15 19/2 17 11 15 8 8 8
    





The following loop iterates over all invertible 6x6 submatrices of *M*
and computes the unique representation of *x* as a linear combination of
the rows of the submatrix. The output (suppressed as it is too long)
shows that each such linear combination requires at least one negative
or one non-integral coefficient.


.. link

.. CODE-BLOCK:: perl

    polymake> foreach (@{all_subsets_of_k(range(0,9),6)}) {
    polymake>   $B = $M->minor($_,All);
    polymake>   if (det($B)) {
    polymake>     print lin_solve(transpose($B),$x), "\n";
    polymake>   }
    polymake> }

This means that *x* cannot be represented as a non-negative linear
combination of any six of the given generators of *C*.

Analyzing the combinatorics
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following is taken from

-  Michael Joswig, Benjamin Müller, and Andreas Paffenholz: ``polymake``
   and lattice polytopes. In Christian Krattenthaler, Volker Strehl and
   Manuel Kauers (eds.), Proceedings of the 21th International
   Conference on Formal Power Series and Algebraic Combinatoric,
   Hagenberg, Austria, 2009, pp. 493-504.


.. link

.. CODE-BLOCK:: perl

    polymake> print $C->N_VERTICES, " ", $C->DIM;
    polymake> print rows_labeled($C->VERTICES_IN_FACETS);

There are two disjoint facets covering all the vertices. Beware the
numbering of facets depends on the convex hull algorithm employed.


.. link

.. CODE-BLOCK:: perl

    polymake> print $C->VERTICES_IN_FACETS->[8];
    polymake> print $C->VERTICES_IN_FACETS->[22];




.. link

.. CODE-BLOCK:: perl

    polymake> print rows_labeled($M);

Here is another polytope which is somewhat similar but not quite the
same.


.. link

.. CODE-BLOCK:: perl

    polymake> $cross5=cross(5);
    polymake> print isomorphic($C,$cross5);
    polymake> print isomorphic($C->GRAPH->ADJACENCY,$cross5->GRAPH->ADJACENCY);




.. link

.. CODE-BLOCK:: perl

    polymake> print $cross5->F_VECTOR - $C->F_VECTOR;

Look at two facets of the five-dimensional cross polytope and their
positions in the dual graph.


.. link

.. CODE-BLOCK:: perl

    polymake> print $cross5->VERTICES_IN_FACETS->[12];
    polymake> print $cross5->VERTICES_IN_FACETS->[13];
    polymake> print rows_labeled($cross5->DUAL_GRAPH->ADJACENCY);

Now we construct a new graph by manipulating the dual graph of the cross
polytope by contracting a perfect matching.


.. link

.. CODE-BLOCK:: perl

    polymake> $g=new props::Graph($cross5->DUAL_GRAPH->ADJACENCY);
    polymake> $g->contract_edge(12,13);
    polymake> $g->contract_edge(24,26);
    polymake> $g->contract_edge(17,21);
    polymake> $g->contract_edge(3,11);
    polymake> $g->contract_edge(6,22);
    polymake> $g->squeeze;

The last command renumbers the nodes sequentially, starting from 0. This
is necessary to render the graph a valid object.


.. link

.. CODE-BLOCK:: perl

    polymake> print isomorphic($C->DUAL_GRAPH->ADJACENCY,$g);

This finally reveals the combinatorial structure: The cone *C* is a cone
over a 5-polytope which can be obtained from the 5-dimensional cross
polytope by \`\ ``straightening`` five pairs of adjacent (simplex)
facets into bipyramids over 3-simplices.
