.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


If you have not used polymake in a long time…
=============================================

…you might want to read up on some things that are important for
backward compatibility.

Numbers
~~~~~~~

``polymake`` always was a hybrid system written half in C++, half in
Perl, but it is only now that the user can directly take advantage of
C++ data types and interfaces in Perl. For instance, via the interface
to `GMP <http://www.swox.com/gmp/>`__ ``polymake`` can also become your
favorite programmable pocket calculator:


.. link

.. CODE-BLOCK:: perl

    polymake> $f = new Integer(1);




.. link

.. CODE-BLOCK:: perl

    polymake> for (my $i = new Integer(100); $i>0; --$i) { $f *= $i; }




.. link

.. CODE-BLOCK:: perl

    polymake> print $f;
    93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000
    





The input of large integer and rational numbers was kind of subtle in
the past, but now (since version 2.11) it has become quite intuitive:


.. link

.. CODE-BLOCK:: perl

    polymake> $bignum=93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000;




.. link

.. CODE-BLOCK:: perl

    polymake> print $bignum;
    93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000
        





.. link

.. CODE-BLOCK:: perl

    polymake> $ratnum=123456/789012;




.. link

.. CODE-BLOCK:: perl

    polymake> print $ratnum;
    10288/65751
        





Each integer constant being too large to fit into a normal perl scalar
value is automatically converted to an ``Integer`` object; each fraction
of two integer constants is automatically converted to a ``Rational``
object (and canonicalized, as can be seen in the example above).

Stored files
~~~~~~~~~~~~

Suppose you still have a file ``cube.poly``, e.g., from trying out the
tutorial of a previous version. You can still do

::

   polymake cube.poly N_FACETS


   N_FACETS
   6

from the command line as you used to. This is the backward compatibility
mode. While this may give the impression that nothing changed and that
you do not have to adapt to the new, this is plain wrong. There are two
things to keep in mind: 1. The old stand-alone clients (such as
``cube``, e.g.) are gone. 2. Once you used the next generation
``polymake`` on your old files they will be transformed into XML
(keeping all your data). In particular, once you called the next
generation ``polymake`` on your files you will not be able to use any
old version on them later.

Equivalent to calling “``cube c3.poly 3``” as before would now be:

::

   polymake 'save(cube(3),"c3.poly")'

A word of warning: It was rarely legal but always popular to edit files
that ``polymake`` worked on with an ASCII text processor. This is still
possible (if you know what you are doing), but in addition to the
caveats previously in place (which are still valid) you have to pay
attention to producing valid XML.
