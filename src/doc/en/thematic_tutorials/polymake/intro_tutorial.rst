.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Introduction to polymake
========================

Installation
~~~~~~~~~~~~

For installation instructions, see the `HowTo </howto/start>`__ pages.
### Getting started

Starting up ``polymake`` on the command line gives:

::

   Welcome to polymake version 3.1
   Copyright (c) 1997-2018
   Ewgenij Gawrilow, Michael Joswig (TU Berlin)
   http://www.polymake.org

   This is free software licensed under GPL; see the source for copying conditions.
   There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

   Press F1 or enter 'help;' for basic instructions.

   Application polytope currently uses following third-party software packages:
   cdd, jreality, libnormaliz, lrs, nauty, permlib, povray, ppl, sketch, sympol, threejs, tikz, tosimplex
   For more details:  show_credits;

You will face an input prompt which indicates the start up application
(usually ``polytope >``). After the prompt you can enter (an enriched
form of) perl code. See `this introduction <tutorial/perl_intro>`__ and
`the advanced polymake/perl tutorial <tutorial/perl_continued>`__ for
more info.

Applications
~~~~~~~~~~~~

``polymake``\ s namespace is subdivided into so-called “applications”.
They group functions and objects that belong together mathematically
(see `this article <howto/lingo#Application>`__ for more). For example,
the default application ``polytope`` deals with polyhedra, while
``topaz`` provides tools for topology and ``tropical`` for tropical
geometry.

To switch to another application, do this:

::

   application 'topaz';

There are `introductory
tutorials <tutorial/start#application_tutorials>`__ for all
applications. Also check out the `documentation </release_docs/3.0>`__.

Getting help
~~~~~~~~~~~~

Apart from the `tutorials <start>`__ and `HowTos </howto/start>`__, you
can often find out what to do to reach your goal by using the built-in
help features of the ``polymake`` shell.

The polymake help system provides extensive usage information. In the
shell, type


::

    polymake> help;

to get started. Invoking


::

    polymake> help 'help';

will explain thy syntax of the help system.

If you type part of an expresison in the shell, you can hit ``TAB`` to
display possible ways of completing it. For example, to see what methods
you can invoke on the graph of the 3-cube, enter the following and hit
``TAB``:

::

    polytope > cube(3)->GRAPH->                              # hit TAB!
   add                               disable_rules                     MAX_CLIQUES                       set_as_default
   ADJACENCY                         dont_save                         name                              set_as_default_now
   apply_rule                        EDGE_DIRECTIONS                   N_CONNECTED_COMPONENTS            SIGNATURE
   attach                            EDGES                             N_EDGES                           SIGNED_INCIDENCE_MATRIX
   AVERAGE_DEGREE                    EIGENVALUES_LAPLACIAN             N_NODES                           SQUARED_EDGE_LENGTHS
   BICONNECTED_COMPONENTS            get_attachment                    NODE_DEGREES                      STRONG_COMPONENTS
   BIPARTITE                         get_schedule                      NODE_IN_DEGREES                   STRONGLY_CONNECTED
   CHARACTERISTIC_POLYNOMIAL         give                              NODE_LABELS                       take
   CONNECTED                         LATTICE_ACCUMULATED_EDGE_LENGTHS  NODE_OUT_DEGREES                  TRIANGLE_FREE
   CONNECTED_COMPONENTS              LATTICE_EDGE_LENGTHS              properties                        type
   CONNECTIVITY                      list_attachments                  provide                           VISUAL
   DEGREE_SEQUENCE                   list_names                        remove                            WEAKLY_CONNECTED
   description                       list_properties                   remove_attachment                 WEAKLY_CONNECTED_COMPONENTS
   DIAMETER   

You can conveniently access the ``polymake`` documentation in the shell.
Place the cursor above a thing you want to know about and hit ``F1``.
Hitting it once will display brief info on types or function signatures.
Hitting it a second time will show you the complete documentation,
sometimes even small usage examples:

::

    polytope > simplex(                                      # hit F1 twice!
   functions/Producing a polytope from scratch/simplex:
   simplex(d; scale, Options) -> Polytope

    Produce the standard d-simplex.
    Combinatorially equivalent to a regular polytope corresponding to the Coxeter group of type A<sub>d-1</sub>.
    Optionally, the simplex can be scaled by the parameter scale.

   Arguments:
     Int d the dimension
     Scalar scale default value: 1

   Options: 
     group => Bool 

   Returns Polytope 

   Examples:


   *) To print the vertices (in homogeneous coordinates) of the standard
      2-simplex, i.e. a right-angled isoceles triangle, type this:
       > print simplex(2)->VERTICES;
       (3) (0 1)
       1 1 0
       1 0 1
      The first row vector is sparse and encodes the origin.

   *) To create a 3-simplex and also calculate its symmetry group, type this:
       > simplex(3, group=>1);

If you have any questions, feel free to ask in the
`forum <https///forum.polymake.org/>`__.

--------------

Back to `Tutorial Overview <start>`__
