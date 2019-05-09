.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Polytopes and Linear Maps
-------------------------

polymake works with `homogeneous coordinates <coordinates>`__, which is
why *projective* linear transformations are natural to apply to
polytopes. Affine transformations are a special case. By the way, a
*transformation* is always bijective, by definition.

Transformations
~~~~~~~~~~~~~~~

We start out with a regular 3-cube …


::

    polymake> $c=cube(3);




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





… and a homethetic image:


::

    polymake> $T=new Matrix<Rational>([[1,0,0,0],[0,2,0,0],[0,0,3,0],[0,0,0,4]]);




::

    polymake> $ct=transform($c,$T);




::

    polymake> print $ct->VERTICES;
    1 -2 -3 -4
    1 2 -3 -4
    1 -2 3 -4
    1 2 3 -4
    1 -2 -3 4
    1 2 -3 4
    1 -2 3 4
    1 2 3 4





Our points are row vectors, so (projective) linear transformations are
applied by multiplying the corresponding matrix from the right. In the
above example the first column of the matrix T is the vector [1,0,0,0]
which means that T acts as an affine map on *R³*. Also the first row
reads [1,0,0,0], and this says that T fixes the origin. This is to say,
T acts linearly.

The purpose of the function transform used above is not only to work on
the VERTICES but also on the FACETS (if available).


::

    polymake> print $c->FACETS;
    1 1 0 0
    1 -1 0 0
    1 0 1 0
    1 0 -1 0
    1 0 0 1
    1 0 0 -1





::

    polymake> print $ct->FACETS;
    1 1/2 0 0
    1 -1/2 0 0
    1 0 1/3 0
    1 0 -1/3 0
    1 0 0 1/4
    1 0 0 -1/4





If we also read the FACETS as row vectors then the corresponding action
is given by the transpose of the inverse of T.

Non-Bijective Linear Maps
~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes we are interested in images of polytopes under a linear map
which is not bijective. An interesting case are projections, for
instance, onto a coordinate subspace.


::

    polymake> $A=new Matrix<Rational>([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,0]]);

Using transform would not work in this case:

::

   # polytope > transform($c,$A);
   polymake:  ERROR: matrix not invertible

The above error says that transform is not the proper function to deal
with this situation as the linear map given by A is not invertible.

To produce the image the following command works:


::

    polymake> $ca=new Polytope<Rational>(POINTS=>$c->VERTICES*$A);




::

    polymake> print $ca->VERTICES;
    1 1 -1 0
    1 1 1 0
    1 -1 -1 0
    1 -1 1 0





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




Since we are applying a non-bijective map, the images of VERTICES do not
have to be VERTICES. Moreover, even if this is the case, multiple
VERTICES may be mapped to the same (like two onto one as in the example
above). If a polytope already has a double description, that is, both
VERTICES and FACETS are known, then the VERTICES and FACETS of the image
under a transformation (that is, a bijective map) cane be read off right
away. However, in the non-bijective case a convex hull computation is
required to compute the FACETS of the image.

Special Examples of Linear Maps to Apply
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[to be continued]


