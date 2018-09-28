.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


ILP and Hilbert bases
---------------------

A first example
~~~~~~~~~~~~~~~

First we will construct a new rational polytope:


::

    polymake> $p=new Polytope<Rational>;




::

    polymake> $p->POINTS=<<".";
    ........> 1 0 0 0
    ........> 1 1 0 0
    ........> 1 0 1 0
    ........> 1 1 1 0
    ........> 1 0 0 1
    ........> 1 1 0 1
    ........> 1 0 1 1
    ........> 1 1 1 1
    ........> .

Note that points in ``polymake`` are always given in homogenous
coordinates. I.e., the point (a,b,c) in R3 is represented as ``1 a b c``
in ``polymake``.

Now we can examine some properties of ``$p``. For instance we can
determine the number of facets or whether ``$p`` is simple:


::

    polymake> print $p->N_FACETS;
    6
        





::

    polymake> print $p->SIMPLE;
    1
    





As you might already have noticed, our polytope is just a 3-dimensional
cube. So there would have been an easier way to create it using the
client ``cube``:


::

    polymake> $c = cube(3,0);

(You can check out the details of any function in the ```polymake``
documentation <http://wwwopt.mathematik.tu-darmstadt.de/polymake_doku/2.9.8/>`__.)

And we can also verify that the two polytopes are actually equal:


::

    polymake> print equal_polyhedra($p,$c);
    1
    





Another example
~~~~~~~~~~~~~~~

Now let us proceed with a somewhat more interesting example: The convex
hull of 20 randomly chosen points on the 2-dimensional sphere. |{{
:tutorial:ilp:rand_sphere.png?200|}}|

.. |{{ :tutorial:ilp:rand_sphere.png?200|}}| image:: attachment:rand_sphere.png


::

    polymake> $rs = rand_sphere(3,20);

``polymake`` can of course visualise this polytope:


::

    polymake> $rs->VISUAL;

Now we will create yet another new polytope by scaling our random sphere
by a factor lambda. (Otherwise there are rather few integral points
contained in it.)

To this end, we have to multiply every coordinate (except for the
homogenising 1 in the beginning) of every vertex by lamda. Then we can
create a new polytope by specifying its vertices.

.. figure:: attachment:rand_sphere_lattice.png
   :alt: {{ :tutorial:ilp:rand_sphere_lattice.png?200|}}

   {{ :tutorial:ilp:rand_sphere_lattice.png?200|}}


::

    polymake> $lambda=2;




::

    polymake> $s=new Matrix<Rational>([[1,0,0,0],[0,$lambda,0,0],[0,0,$lambda,0],[0,0,0,$lambda]]);




::

    polymake> print $s;
    1 0 0 0
    0 2 0 0
    0 0 2 0
    0 0 0 2
        





::

    polymake> $scaled_rs=new Polytope<Rational>(VERTICES=>($rs->VERTICES * $s), LINEALITY_SPACE=>[]);

``polymake`` can visualise the polytope together with its lattice
points:


::

    polymake> $scaled_rs->VISUAL->LATTICE_COLORED;

Now will construct the integer hull of ``$scaled_rs`` and visualise it:

.. figure:: attachment:ilp_lattice.png
   :alt: {{ :tutorial:ilp:ilp_lattice.png?200|}}

   {{ :tutorial:ilp:ilp_lattice.png?200|}}


::

    polymake> $integer_hull=new Polytope<Rational>(POINTS=>$scaled_rs->LATTICE_POINTS);




::

    polymake> $integer_hull->VISUAL->LATTICE_COLORED;

In order to obtain the integer hull we simply define a new polytope
``$integer_hull`` as the convex hull of all ``LATTICE_POINTS`` contained
in ``$scaled_rs``.

Note that if we give ``POINTS`` (in contrast to ``VERTICES``)
``polymake`` constructs a polytope that is the convex hull of the given
points regardless of whether they are vertices or not. I.e., redundacies
are allowed here.

If you specify ``VERTICES`` you have to make sure yourself that your
points are actually vertices since ``polymake`` does not check this. You
also need to specify the ``LINEALITY_SPACE``, see `Tutorial on
polytopes <tutorial/apps_polytope>`__.

Linear Programming
~~~~~~~~~~~~~~~~~~

Now that we have constructed a nice integral polytope we want to apply
some linear program to it.

First we define a ``LinearProgram`` with our favourite
``LINEAR_OBJECTIVE``. The linear objective is an given as a vector of
length d+1, d being the dimension of the space. The vector [c0,c1, …,
cd] corresponds to the linear objective c0 + c1x1 + … + cdxd.


::

    polymake> $objective=new LinearProgram<Rational>(LINEAR_OBJECTIVE=>[0,1,1,1]);

Then we define a new polytope, which is a copy of our old one
(``$inter_hull``) with the LP as an additional property.


::

    polymake> $ilp=new Polytope<Rational>(VERTICES=>$integer_hull->VERTICES, LP=>$objective);

|{{ :tutorial:ilp:ilp_min_face.png?200|}}| |{{
:tutorial:ilp:ilp_max_face.png?200|}}|

And now we can perform some computations:

.. |{{ :tutorial:ilp:ilp_min_face.png?200|}}| image:: attachment:ilp_min_face.png
.. |{{ :tutorial:ilp:ilp_max_face.png?200|}}| image:: attachment:ilp_max_face.png


::

    polymake> print $ilp->LP->MAXIMAL_VALUE;
    2
        





::

    polymake> print $ilp->LP->MAXIMAL_FACE;
    {6 9 10}
        





::

    polymake> $ilp->VISUAL->MIN_MAX_FACE;

Hence the LP attains its maximal value 2 on the 2-face spanned by the
vertices 6, 9 and 10.

``polymake`` can visualise the polytope and highlight both its maximal
and minimal face in a different (by default admittedly almost painful
;-) ) colour. Here you see the maximal face ``{6 9 10}`` in red and the
minimal face ``{0 3}`` (on the opposite side of the polytope) in yellow.

Note though that since we started out with a random polytope these
results may vary if we perform the same computations another time on a
different random polytope.


::

    polymake> print $ilp->VERTICES;
    1 -1 0 -1
    1 -1 0 1
    1 -1 1 0
    1 0 -1 -1
    1 0 -1 1
    1 0 1 -1
    1 0 1 1
    1 1 -1 0
    1 1 0 -1
    1 1 0 1
    1 1 1 0
    





Hilbert bases
~~~~~~~~~~~~~

Finally, we can have ``polymake`` compute and print a Hilbert basis for
the cone spanned by ``$ilp``. Notice that this requires normaliz or 4ti2
to be installed in order to work.


::

    polymake> print $ilp->HILBERT_BASIS;
    1 0 0 -1
    1 -1 1 0
    1 1 0 0
    1 0 1 0
    1 0 1 -1
    1 1 1 0
    1 0 1 1
    1 1 0 -1
    1 1 0 1
    1 0 0 0
    1 0 0 1
    1 1 -1 0
    1 -1 0 -1
    1 -1 0 0
    1 -1 0 1
    1 0 -1 -1
    1 0 -1 0
    1 0 -1 1
    






