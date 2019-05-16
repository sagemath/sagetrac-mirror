.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Matching Polytopes
------------------

In this tutorial we will use ``polymake`` to construct and analyse
matching polytopes.

First we construct a graph, the complete graph on four nodes:


.. link

.. CODE-BLOCK:: perl

    polymake> $K4=new props::Graph(4);

.. link

.. CODE-BLOCK:: perl

    polymake> for (my $i=0; $i<4; ++$i) {
    polymake>   for (my $j=$i+1; $j<4; ++$j) {
    polymake>     $K4->edge($i,$j);
    polymake>   }
    polymake> }

(See also the `Tutorial on Graphs <apps_graph>`__ for more on the
construction of graphs.)

Next we like to have the node-edge-incidence matrix of our graph. Since
the latest release of ``polymake`` does not yet support this, we have to
write the function ourselves:


.. link

.. CODE-BLOCK:: perl

    polymake> sub node_edge_incidences {
    polymake>     my $g=shift;
    polymake>    my $A=new Matrix<Int>($g->nodes, $g->edges);
    polymake>     my $k=0;
    polymake>    for (my $i=0; $i<$g->nodes-1; ++$i) {
    polymake>        foreach (@{$g->adjacent_nodes($i)}) {
    polymake>            if ($_>$i) {
    polymake>                $A->[$i]->[$k]=1;
    polymake>                $A->[$_]->[$k]=1;
    polymake>                ++$k;
    polymake>            }
    polymake>        }
    polymake>    }
    polymake>    return $A;
    polymake> }

Now we can construct the node-edge-incidence matrix of our graph ``K4``:


.. link

.. CODE-BLOCK:: perl

    polymake> $A=node_edge_incidences($K4);
    polymake> print $A;
    1 1 1 0 0 0
    1 0 0 1 1 0
    0 1 0 1 0 1
    0 0 1 0 1 1





With this we can now construct the constraint matrix consisting of an
upper part for the nonnegativity constraints xe

.. raw:: html

   <html>

≥

.. raw:: html

   </html>

0 …


.. link

.. CODE-BLOCK:: perl

    polymake> $I=new Matrix<Int>([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]);
    polymake> $Block1=new Matrix<Int>(new Vector<Int>([0,0,0,0,0,0]) | $I);

… and a lower part for the constraints

.. raw:: html

   <html>

Σ

.. raw:: html

   </html>

e xe

.. raw:: html

   <html>

≤

.. raw:: html

   </html>

1 for each vertex v

.. raw:: html

   <html>

∈

.. raw:: html

   </html>

V, where the sum is over all edges e containing v:


.. link

.. CODE-BLOCK:: perl

    polymake> $Block2=new Matrix<Int>(new Vector<Int>([1,1,1,1]) | -$A);

Now we can put both parts together and define the polytope:


.. link

.. CODE-BLOCK:: perl

    polymake> $Ineqs=new Matrix<Rational>($Block1 / $Block2);
    polymake> $P=new Polytope<Rational>(INEQUALITIES=>$Ineqs);

The matching polytope of ``K4`` is the integer hull of ``P``:


.. link

.. CODE-BLOCK:: perl

    polymake> $P_I=new Polytope<Rational>(POINTS=>$P->LATTICE_POINTS);


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




We can analyse some elementary properties of ``P_I`` …


.. link

.. CODE-BLOCK:: perl

    polymake> print $P_I->POINTS;
    1 0 0 0 0 0 0
    1 0 0 0 0 0 1
    1 0 0 0 0 1 0
    1 0 0 0 1 0 0
    1 0 0 1 0 0 0
    1 0 0 1 1 0 0
    1 0 1 0 0 0 0
    1 0 1 0 0 1 0
    1 1 0 0 0 0 0
    1 1 0 0 0 0 1





.. link

.. CODE-BLOCK:: perl

    polymake> print $P_I->FACETS;
    0 0 0 0 0 0 1
    0 1 0 0 0 0 0
    1 0 0 0 -1 -1 -1
    1 -1 0 0 -1 -1 0
    1 0 -1 0 -1 0 -1
    1 -1 -1 0 -1 0 0
    1 0 0 -1 0 -1 -1
    1 -1 0 -1 0 -1 0
    1 0 -1 -1 0 0 -1
    1 -1 -1 -1 0 0 0
    0 0 0 0 0 1 0
    0 0 1 0 0 0 0
    0 0 0 0 1 0 0
    0 0 0 1 0 0 0





.. link

.. CODE-BLOCK:: perl

    polymake> print $P_I->N_FACETS;
    14




… and compare them with the according properties of the defining
polytope ``P``:


.. link

.. CODE-BLOCK:: perl

    polymake> print $P->VERTICES;
    1 0 0 0 1 0 0
    1 0 1 0 0 0 0
    1 1/2 1/2 0 1/2 0 0
    1 0 0 0 0 0 0
    1 1 0 0 0 0 0
    1 1/2 0 1/2 0 1/2 0
    1 0 1/2 1/2 0 0 1/2
    1 0 0 0 1/2 1/2 1/2
    1 0 0 0 0 1 0
    1 0 0 1 0 0 0
    1 0 0 0 0 0 1
    1 1 0 0 0 0 1
    1 0 1 0 0 1 0
    1 0 0 1 1 0 0





.. link

.. CODE-BLOCK:: perl

    polymake> print $P->VOLUME;
    1/72




.. link

.. CODE-BLOCK:: perl

    polymake> print $P_I->VOLUME;
    1/90




Next we analyse the combinatorics of ``P_I``: |{{
:tutorial:ilp:gale.png?300|The Gale diagram of ``facet0``}}|

.. |{{ :tutorial:ilp:gale.png?300|The Gale diagram of ``facet0``}}| image:: attachment:gale.png


.. link

.. CODE-BLOCK:: perl

    polymake> print $P_I->AMBIENT_DIM, " ", $P_I->DIM;
    6 6




.. link

.. CODE-BLOCK:: perl

    polymake> print $P_I->F_VECTOR;
    10 39 78 86 51 14




.. link

.. CODE-BLOCK:: perl

    polymake> print $P_I->FACET_SIZES;
    8 8 6 6 6 6 6 6 6 6 8 8 8 8




.. link

.. CODE-BLOCK:: perl

    polymake> $facet0=facet($P_I,0);




.. link

.. CODE-BLOCK:: perl

    polymake> print $facet0->AMBIENT_DIM, " ", $facet0->DIM;
    6 5




.. link

.. CODE-BLOCK:: perl

    polymake> print rows_labeled($facet0->VERTICES_IN_FACETS);
    0:0 1 2 3 4 5 6
    1:1 2 4 6 7
    2:2 4 5 6 7
    3:1 3 4 6 7
    4:3 4 5 6 7
    5:0 2 3 4 5 7
    6:0 1 2 3 4 7
    7:0 1 3 5 6 7
    8:0 1 2 5 6 7





.. link

.. CODE-BLOCK:: perl

    polymake> $facet0->GALE;

The Gale diagram of ``facet0`` is depicted on the right.
