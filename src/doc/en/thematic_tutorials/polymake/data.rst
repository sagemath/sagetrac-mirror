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


::

    polymake> $p = cube(3);

To store polymake objects use the command


::

    polymake> save($p,"myPolyObject.poly");

polymake objects that are stored in polymake’s own XML file format can
be loaded via


::

    polymake> $p=load("myPolyObject.poly");

If you did not start ``polymake`` in the directory containing your
object, it is necessary to add the relative path, e.g.
$p=load(“MyFolder/myPolyObject.poly”);

**Note:** If you load a polymake object and compute new properties,
these properties will automatically be added to the original XML-file at
the end of the session. You can suppress this with the command


::

    polymake> $p->dont_save;

called prior to leaving the session (but after the last compuation with
$p).

Handling complex data types
---------------------------

It is also possible to store complex data structures in XML format via
``save_data``, e.g.


::

    polymake> $s=new Set<Int>(1,2,3,4);
    ........> save_data($s,"mySet.poly");

To load such files just type


::

    polymake> $s=load_data("mySet.poly");

Handling arbitrary files
------------------------

Of course, it is also possible to load data from files in other formats.
For this purpose use the standard Perl functions for reading and
writing. Here is an example:

Assume you want to load some points stored in the file points.txt which
looks like this: 1 0 0 0 1 1 0 0 1 0 1 0 1 1 1 0 1 0 0 1 1 1 0 1 1 0 1 1
1 1 1 1 For the sake of the example, let’s create this file:


::

    polymake> open(my $f, '> points.txt'); print $f "1 0 0 0\n1 1 0 0\n1 0 1 0\n1 1 1 0\n1 0 0 1\n1 1 0 1\n1 0 1 1\n1 1 1 1\n"; close $f;

To read this file try the following:


::

    polymake> open(INPUT, "< points.txt");
    ........> while(<INPUT>){
    ........>   print $_;
    ........> }
    ........> close(INPUT);

 is a perl input iterator reading the file line by line. Variable ``$_``
refers to the current line within this loop; it has a plain string
value.

A reasonable task could be to store the points from the file as a
matrix. This can be done immediately, because the matrix constructor
called with a list of values interprets each value as a matrix line:


::

    polymake> open(INPUT, "< points.txt");
    ........> $matrix=new Matrix<Rational>(<INPUT>);
    ........> close(INPUT);


