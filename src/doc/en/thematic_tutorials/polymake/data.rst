.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Save and load data in polymake
==============================

In polymake there are different ways to save and load data depending on
the type and the format of the data. We distinguish between polymake
objects (Polytope, Matroid,…), complex data types (Set, Matrix,
Array<Vector >,…), and data from files in arbitrary formats.

Handling polymake objects
-------------------------

Let us take this nice example object:


.. link

.. CODE-BLOCK:: perl

    polymake> $p = cube(3);

To store polymake objects use the command


.. link

.. CODE-BLOCK:: perl

    polymake> save($p,"myPolyObject.poly");

This silently overwrites existing files.

polymake objects that are stored in polymake’s own XML file format can
be loaded via


.. link

.. CODE-BLOCK:: perl

    polymake> $p=load("myPolyObject.poly");

If you did not start ``polymake`` in the directory containing your
object, it is necessary to add the relative or absolute path, e.g.
$p=load(“MyFolder/myPolyObject.poly”); TAB completion like in a usual
UNIX shell supports you in navigating through the file system.

**Note:** If you load a polymake object and compute new properties,
these properties will automatically be added to the original XML-file at
the end of the session. You can suppress this with the command


.. link

.. CODE-BLOCK:: perl

    polymake> $p->dont_save;

called prior to leaving the session (but after the last computation with
$p).

If you want to store a collection of objects into a single file, there
is an `extra tutorial <.tarballs>`__ for you. ## Handling complex data
types

Apart from the full objects, you can also persistently store arbitrary
data structures like matrices or graphs in XML format via ``save_data``,
e.g.


.. link

.. CODE-BLOCK:: perl

    polymake> $s=new Set<Int>(1,2,3,4);
    polymake> save_data($s, "mySet.poly", "My very own set.");

The description text is optional; it can be an arbitrary text, even
stretching over several lines.

To load such files just type


.. link

.. CODE-BLOCK:: perl

    polymake> $s=load_data("mySet.poly");

Saving visualized objects
-------------------------

Furthermore, most visualization methods provide an option to save the
visualized object in a suitable format. Consult the `F1
help <:user_guide:intro_tutorial#getting_help>`__ for information on the
file format and further options.

To save the cube visualized via JReality in a new file called
``mycube.bsh``, do this:

::

   jreality(cube(3)->VISUAL,File=>"mycube");

To save the cube as a TiKz file named ``mycube.tikz`` that you can
e.g. import in a LaTeX document, do this instead:

::

   tikz(cube(3)->VISUAL,File=>"mycube");

Handling arbitrary files
------------------------

Of course, it is also possible to load data from files in other formats.
For this purpose use the standard Perl functions for reading and
writing. Here is an example:

Assume you want to load some points stored in the file points.txt which
looks like this: 1 0 0 0 1 1 0 0 1 0 1 0 1 1 1 0 1 0 0 1 1 1 0 1 1 0 1 1
1 1 1 1 For the sake of the example, let’s create this file:


.. link

.. CODE-BLOCK:: perl

    polymake> open(my $f, '> points.txt'); print $f "1 0 0 0\n1 1 0 0\n1 0 1 0\n1 1 1 0\n1 0 0 1\n1 1 0 1\n1 0 1 1\n1 1 1 1\n"; close $f;

To read this file try the following:


.. link

.. CODE-BLOCK:: perl

    polymake> open(INPUT, "< points.txt");
    polymake> while(<INPUT>){
    polymake>   print $_;
    polymake> }
    polymake> close(INPUT);

``<INPUT>`` is a perl input iterator reading the file line by line.
Variable ``$_`` refers to the current line within this loop; it has a
plain string value.

A reasonable task could be to store the points from the file as a
matrix. This can be done immediately, because the matrix constructor
called with a list of values interprets each value as a matrix line:


.. link

.. CODE-BLOCK:: perl

    polymake> open(INPUT, "< points.txt");
    polymake> $matrix=new Matrix<Rational>(<INPUT>);
    polymake> close(INPUT);


