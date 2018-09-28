.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Saving and Restoring an Array of Polytopes
==========================================

If you want to deal with a whole family of polytopes at the same time,
it will sometimes be convenient to save and restore them to a single
file. polymake has a simple mechanism for this, storing your array of
polytopes into a single tarball.

The necessary functions for this are contained in the script “tarballs”
that you can load into your polymake session by calling


::

    polymake> script("tarballs");

It provides the two functions ``pack_tarball`` and ``unpack_tarball``.

Storing
-------

Here is a simple example, where we create an array @a containing a cube
and a simplex and save this to a file.


::

    polymake> @a = ();
    ........> $a[0] = cube(3);
    ........> $a[1] = simplex(3);
    ........> pack_tarball("simple_polytopes.tgz",@a);

This creates a file ``simple_polytopes.tgz`` in the current directory
that is a tarred (and gzipped) archive containing two polymake files.
polymake detects whether you want the archive gzipped or not from the
supplied file extension (``*.tar.gz`` or ``*.tgz``). You can verify this
by calling ``tar tvfz simple_polytopes.tgz`` from the command line:

::

   [nightingale]:~/temp>tar tvfz simple_polytopes.tgz 
   -rw------- xxx/yyy 1468 2009-07-01 17:20 1.poly
   -rw------- xxx/yyy  854 2009-07-01 17:20 2.poly
   [nightingale]:~/temp>

If you want to get more descriptive names for your polymake files then
you have to set a name for each polytope first.


::

    polymake> $a[0]->name = "my_cube";
    ........> $a[1]->name = "my_simplex";
    ........> pack_tarball("simple_polytopes.tgz",@a);

sets the names of the files in the tarball to ``my_cube.poly`` and
``my_simplex.poly``:

::

   [nightingale]:~/temp>tar tvfz simple_polytopes.tgz 
   -rw------- xxx/yyy  952 2009-07-01 17:21 my_cube.poly
   -rw------- xxx/yyy  650 2009-07-01 17:21 my_simplex.poly
   [nightingale]:~/temp>

Restoring the array
-------------------

You can restore your saved array by using the function
``unpack_tarball``:


::

    polymake> @a=unpack_tarball("simple_polytopes.tgz");
    ........> print $a[0]->name;
    my_cube
    





If you just want a specific polytope from your tarball, then you can
supply its name in the command:


::

    polymake> @a=unpack_tarball("simple_polytopes.tgz","my_simplex.poly");
    ........> print $a[0]->name;
    my_simplex
    





You may supply more than one filename. However, wildcards are not
supported. Note that changes in the files are not automatically stored
in the archive, you have to call ``pack_tarball`` to update the files.

Packing archives outside polymake
---------------------------------

You can of course apply ``tar`` to your favorite family of polymake
files to create a tarball without using polymake. It can be read by
polymake as long as the files are at the root of the archive (i.e. don’t
pack a whole directory tree). The archive also may not contain
non-polymake files.
