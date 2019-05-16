.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Dealing with Graphs
===================

Graphs are ubiquitous in geometric combinatorics. Hence they occur a lot
throughout the ``polymake`` system, explicitly and implicitly. It is
important to understand that the user encounters graphs on two distinct
layers in the object hierarchy. It is the purpose of this tutorial to
explore the various features. For the sake of simplicity here we
restrict our attention to undirected graphs.

Graphs of Polytopes
-------------------

Coming from polytopes the first situation in which a graph occurs is the
vertex-edge graph of such a polytope.


.. link

.. CODE-BLOCK:: perl

    polymake> $p=rand_sphere(3,20);
    polymake> print $p->GRAPH->N_NODES;
    20
    





Here ``GRAPH`` is a property of the polytope object ``$p`` which happens
to be of the object type ``Graph``. The following is a fragment of the
file ``apps/polytopes/rules/polytope_properties.rules``. This is where
all the standard properties of polytopes are declared.

::

   property GRAPH : objects::Graph {

      # Difference of the vertices for each edge (only defined up to signs).
      property EDGE_DIRECTIONS : EdgeMap<Undirected, Vector<Scalar>>;
   }

In fact, this ``objects::Graph`` is the main object class of another
application, called ``graph``. This application ``graph`` is defined in
``apps/graph`` and its subdirectories. Although the application
``graph`` has much fewer features than the application ``polytope`` the
overall mechanism of interaction is the same. In particular, there are
properties, rules, clients, and such. Non-trivial features include the
algorithm for computing the diameter and the visualization of graphs
based on a pseudo-physical model (described in `this
paper <http://front.math.ucdavis.edu/0711.2397>`__).


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->DIAMETER;
    4





.. link

.. CODE-BLOCK:: perl

    polymake> $p->VISUAL;

As polytopes and all other objects in ``polymake``\ ’s object hierarchy
the graphs from the application ``graph`` are immutable objects. It is
not possible to add a node or to delete an edge. It is instructive to
look at the beginning of the file ``apps/graph/graph_properties.rules``
which reads like this.

::

   declare object Graph<Dir=Undirected> {

   # combinatorial description of the Graph in the form of adjacency matrix
   property ADJACENCY : props::Graph<Dir>;

   ...

The key property of a graph object is its ``ADJACENCY``. For each node
all neighbors are listed. Here we want to focus on the type
``props::Graph`` of this property. This refers to C++ class from the
Polymake Template Library named ``Graph``, and this is where the data
structure and most algorithms reside. It is possible to directly
manipulate objects of this type, and these are *not* immutable, they can
be changed. The following shows how one can create a 5-cycle. Calling
the method ``edge`` creates an edge if it did not exist before. The
output is the ordered list of neighbors per node.


.. link

.. CODE-BLOCK:: perl

    polymake> $g=new props::Graph(5);                      
    polymake> for (my $i=0; $i<5; ++$i) { $g->edge($i,($i+1)%5) };
    polymake> print $g;
    {1 4}               
    {0 2}               
    {1 3}               
    {2 4}               
    {0 3}               
    





If a graph has many nodes it is convenient to know which line of the
output refers to which node. If an array of labels is given this could
also be used instead of the numbers which are the default.


.. link

.. CODE-BLOCK:: perl

    polymake> print rows_labeled($g);
    0:1 4                             
    1:0 2                             
    2:1 3                             
    3:2 4                             
    4:0 3                             
    





There are other ways to change such a graph. Contracting the edge
*(x,y)* where *x* is smaller than *y* implies that the node *y* is
destroyed.


.. link

.. CODE-BLOCK:: perl

    polymake> $g->delete_edge(0,1);
    polymake> $g->contract_edge(2,3);
    polymake> $g->squeeze();

However, most of our graph algorithms expect a graph with consecutively
numbered nodes. The function ``squeeze`` takes care of a proper
renumbering, but this takes linear time in the number of nodes.


.. link

.. CODE-BLOCK:: perl

    polymake> print rows_labeled($g);
    0:4
    1:2
    2:1 4
    3:0 2
    





How do I iterate over the adjacent nodes to a given node?

::

   foreach (@{$g->adjacent_nodes(0)}) {
      print "node number $_ is adjacent to node number 0\n";
   }

It is also legal to copy all adjacent nodes to an array as in:

::

   @x = @{$g->adjacent_nodes(0)};

Subsequently, the individial neighbors can be accessed, for instance, as
``$x[1]``. However, for technical reasons too difficult to explain here,
it is *not* legal to write ``$g->adjacent_nodes(0)->[1]``! Usually it is
preferred to avoid copying; so use constructions like ``foreach`` and
``map`` if possible.

Defining a Graph from Scratch
-----------------------------

You can also work with graphs independent of their connection to
polytopes. We will switch to ``application "graph"`` for the following
commands, but this is not strictly necessary. We want to define a new
object of type ``Graph`` in ``polymake``.

The key property of a graph is its adjacency matrix, which is stored in
the property ``ADJACENCY``. It lists the neighbors of each node. We use
again the above example of a 5-cycle C5 with consecutively numbered
nodes. Then one can define C5 by


.. link

.. CODE-BLOCK:: perl

    polymake> application "graph";




.. link

.. CODE-BLOCK:: perl

    polymake> $g=new objects::Graph(ADJACENCY=>[[1,4],[0,2],[1,3],[2,4],[0,3]]);

Note the ``objects::`` in front of the key word ``Graph``, which is not
needed when you define any of the other ``polymake`` objects, like
``Polytope<Rational> or``\ Matroid\ ``.  This is necessary here to distinguish the``\ polymake\ ``object``\ Graph\ ``from the``\ C++\ ``class``\ Graph\ ``that we have used above, and that is accessed with the additional qualification``\ props::`.

The list of edges of the graph is induced by the adjacency matrix
(please note that in a undirected graph each edge appears twice). You
can get an explicit list of the edges with the user function ``EDGES``.


.. link

.. CODE-BLOCK:: perl

    polymake> print $g->EDGES;
    {0 1}
    {1 2}
    {2 3}
    {0 4}
    {3 4}
    





Note however, that this list is not stored in the object, as it is just
a different view on the adjacency matrix.

Most often when you define a graph you would not write it down as a list
of adjacencies, but as a list of edges. For convenience, ``polymake``
provides a way to create a graph from a list of edges. The same 5-cycle
as above could also be defined via


.. link

.. CODE-BLOCK:: perl

    polymake>  $g=graph_from_edges([[0,1],[1,2],[2,3],[0,4],[3,4]]);

The order of the edges, and the order of the nodes for each edge in a
undirected case, is not important. We can check the adjacency matrix,


.. link

.. CODE-BLOCK:: perl

    polymake> print $g->ADJACENCY;
    {1 4}
    {0 2}
    {1 3}
    {2 4}
    {0 3}
    





and continue to work with the graph by e.g. checking its ``DIAMETER``,
``BIPARTITE``-ness or other properties:


.. link

.. CODE-BLOCK:: perl

    polymake> print $g->DIAMETER;
    2





.. link

.. CODE-BLOCK:: perl

    polymake> print $g->BIPARTITE;
    0





.. link

.. CODE-BLOCK:: perl

    polymake> print $g->MAX_CLIQUES;
    {{0 1} {0 4} {1 2} {2 3} {3 4}}
    





Directed Graphs
---------------

By specifying the template parameter ``Directed`` a graph is born as a
directed graph. Properties which make sense for directed graphs work as
expected. A directed graph may have two arcs between any two nodes with
opposite orientations.


.. link

.. CODE-BLOCK:: perl

    polymake> $g=new objects::Graph<Directed>(ADJACENCY=>[[1],[2],[3],[2,4],[0]]);
    polymake> print $g->DIAMETER;
    4
    





Some properties of graphs do not make sense for directed graph. Here is
an example of an undirected graph property which does not make sense for
directed graphs.


.. link

.. CODE-BLOCK:: perl

    polymake> #print $g->MAX_CLIQUES;

::

   polymake:  ERROR: Object Graph<Directed> does not have a property or method MAX_CLIQUES

Graphs with multiple edges/arcs are currently not supported.

Visualizing Graphs
------------------

Like other “big” ``polymake`` objects the ``Graph`` class has a member
(function) ``VISUAL`` which returns an abstract visualization object.
Depending on the configuration it typically uses ``JReality`` or
``JavaView``. Particularly interesting for graph drawing is the
visualization via ``Graphviz``.


.. link

.. CODE-BLOCK:: perl

    polymake> graphviz($g->VISUAL);

Note that the latter starts a postscript viewer with the ``Graphviz``
output. Make sure that the custom variable ``$Postscript::viewer`` is
set to something reasonable (like, e.g., ``/usr/bin/evince``).
