.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


General chain complexes intopaz======
=====================================

Apart from being capable of computing integer homology of simplicial
complexes (see this `tutorial <topaz_tutorial>`__ for an introduction),
``polymake`` is able to handle general chain complexes and compute
homology for coefficients from different domains. When experimenting in
the interactive shell, switch to the topology application first:


::

    polymake> application 'topaz';

Constructing a ChainComplex
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can construct a chain complex via its differential matrices. For
example purposes, we use the sparse boundary matrices of a triangulation
of the real projective plane. You can then construct a general chain
complex from it like this:


::

    polymake> $bd1 = real_projective_plane()->boundary_matrix(1);
    ........> $bd2 = real_projective_plane()->boundary_matrix(2);
    ........> $a = new Array<SparseMatrix<Integer>>($bd1,$bd2);   # omit the trivial zeroth differential
    ........> $cc = new ChainComplex<SparseMatrix<Integer>>($a,1);

The template parameter of ``ChainComplex`` denotes the type of the
boundary matrices. It defaults to \`SparseMatrix, as this allows
computation of integer homology. The second parameter of the chain
complex constructor defaults to 0, indicating whether to perform a
sanity check on the matrices (i.e. whether matrix dimensions match and
successive maps compose to the zero map).

You can access the data stored in the object like this:


::

    polymake> print $cc->boundary_matrix(2);
    (15) (0 1) (1 -1) (2 1)
    (15) (0 1) (3 -1) (4 1)
    (15) (5 1) (6 -1) (7 1)
    (15) (1 -1) (5 1) (8 1)
    (15) (3 -1) (6 1) (9 1)
    (15) (7 1) (10 1) (11 -1)
    (15) (4 -1) (10 1) (12 1)
    (15) (2 -1) (11 1) (13 1)
    (15) (8 1) (12 -1) (14 1)
    (15) (9 -1) (13 1) (14 1)
    





Computing integer homology
~~~~~~~~~~~~~~~~~~~~~~~~~~

There is a user function to compute integer homology of your complex.
You can access the documentation by typing the name of the function in
the interactive shell and then pressing F1.


::

    polymake> print homology($cc,0);
    ({} 1)
    ({(2 1)} 0)
    ({} 0)
    





The output rows correspond to the dimensions of your homology modules,
containing the torsion coefficients in curly brackets, and the betti
number. Note that this is non-reduced homology, unlike what gets
computed when using the ``HOMOLOGY`` property of a simplicial complex.

There is an extra function for computing the generators of the homology
modules as well.


::

    polymake> print homology_and_cycles($cc,0);
    (({} 1)
    <(6) (0 1)
    >
    )
    (({(2 1)} 0)
    <(15) (10 1) (11 -1) (12 1) (13 -1) (14 -1)
    >
    )
    (({} 0)
    <>
    )
    





The output pairs the homology module representation with a
representation of the cycles generating the respective modules, where
the indices correspond to the indices in your input matrices.

Computing Betti numbers
~~~~~~~~~~~~~~~~~~~~~~~

If your complex’ differentials do not have ``Integer`` coefficients,
computing integer homology is not possible. You can still (and very
efficiently!) compute the Betti numbers by using the corresponding user
function:

::

   print betti_numbers($cc);
   1 0 0
