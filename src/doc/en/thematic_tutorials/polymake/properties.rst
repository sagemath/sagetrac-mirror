.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Tutorial on Properties and Rules
================================

Properties
~~~~~~~~~~

Each object has a list of properties of various types. When an object is
‘born’ it comes with an initial list of properties, and all other
properties will be derived from those. Here we discuss an example from
the ``polytope`` application. The following creates a 3-dimensional
cube.


::

    polymake> $c=cube(3);

The object is defined by calling some function, but how does one find
out what the initial set of properties is? Of course, one could look at
the source code, but the following is the direct way from the
interpreter.


::

    polymake> print join(", ", $c->list_properties);
    AMBIENT_DIM, DIM, FACETS, VERTICES_IN_FACETS, BOUNDED
    





The relevant method, which is defined for any ``polymake`` object, is
called ``list_properties``. It returns an array of strings. The extra
code is just there to print this list nicely. The object is changed if
we ask for a property which has not been computed before.


::

    polymake> print $c->VERTICES;
    1 -1 -1 -1
    1 1 -1 -1
    1 -1 1 -1
    1 1 1 -1
    1 -1 -1 1
    1 1 -1 1
    1 -1 1 1
    1 1 1 1
        





::

    polymake> print join(", ", $c->list_properties);
    AMBIENT_DIM, DIM, FACETS, VERTICES_IN_FACETS, BOUNDED, N_VERTICES, FEASIBLE, SIMPLE, SIMPLE_POLYHEDRON, AFFINE_HULL, VERTICES
    





The property ``VERTICES`` was added, but also a few others. These were
computed on the way. Which properties show up after some computation
depends on the rules applied. What is the set of properties that *can*
be computed for a given object? This depends on your set of rule valid
for the object in question. Here is a short sequence of commands which
lets you find out. The properties listed come in alphabetical ordering.


::

    polymake> $t=$c->type;
    ........> print join(", ", sorted_uniq(sort { $a cmp $b } map { keys %{$_->properties} } $t, @{$t->super}));

Instead of showing the (lengthy) enumeration have a look at the
`documentation <release_docs/latest/polytope.html>`__ for a complete
list of properties known for objects of the application ``polytope``.

Schedules
~~~~~~~~~

[beware: output from branch “cones”]

Let us restart with our cube from scratch.


::

    polymake> $c=cube(3);
    ........> print join(", ", $c->list_properties);
    POLYTOPE_AMBIENT_DIM, POLYTOPE_DIM, FACETS, VERTICES_IN_FACETS, BOUNDED
    





Suppose we want to see which sequence of rules leads to the computation
of the F_VECTOR.


::

    polymake> $schedule=$c->get_schedule("F_VECTOR");
    ........> print join("\n", $schedule->list);
    HASSE_DIAGRAM : RAYS_IN_FACETS
    F_VECTOR : HASSE_DIAGRAM
    





Applying the schedule to the object yields the same as asking for the
property right away.


::

    polymake> $schedule->apply($c);
    ........> print join(", ", $c->list_properties);
    POLYTOPE_AMBIENT_DIM, POLYTOPE_DIM, FACETS, VERTICES_IN_FACETS, BOUNDED, HASSE_DIAGRAM, F_VECTOR
    





It is possible to apply the same schedule to several polytopes. This is
useful for a slight speed up in the total time of the computation.
