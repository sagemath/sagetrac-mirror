.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


A short note on variable naming up front: You can alter the settings for
the names that are used for polynomial variables in parsing them from
strings or for pretty-printing using the ``set_var_names`` or
``local_var_names`` functions. Refer to their
`documentation <https://polymake.org/release_docs/master/common.html#common__set_var_names__239>`__
for defaults and usage information. To restore the default settings to
your ``customize.pl``, type this:


.. link

.. CODE-BLOCK:: perl

    polymake> reset_custom %polynomial_var_names;

Usage of Polynomials in Perl
----------------------------

Constructors
^^^^^^^^^^^^

The easiest way to create a simple
`Polynomial <https://polymake.org/release_docs/master/common.html#common__Polynomial__339>`__
object is through a string:


.. link

.. CODE-BLOCK:: perl

    polymake> $p = new Polynomial("4 + 3x_1 + x_2^5");

Sometimes it’s convenient to use the constructor that takes a vector of
coefficients and a matrix of exponents:


.. link

.. CODE-BLOCK:: perl

    polymake> $coeff = new Vector([9,-5]);
    polymake> $exp = new Matrix<Int>([[0,4],[8,3]]);
    polymake> $p2 = new Polynomial($coeff, $exp);
    polymake> print $p2;
    -5*x_0^8*x_1^3 + 9*x_1^4




There is a seperate type for univariate polynomials, called
`UniPolynomial <https://polymake.org/release_docs/master/common.html#common__UniPolynomial__342>`__.


.. link

.. CODE-BLOCK:: perl

    polymake> $up = new UniPolynomial("3x + 2x^2 + 4");

Polynomials (and UniPolynomials) are templated by their coefficient and
exponent types, defaulting to Rational for coefficients and Int for
exponents. You can even have polynomials of polynomials (of
polynomials…).


.. link

.. CODE-BLOCK:: perl

    polymake> $pp = new UniPolynomial<UniPolynomial<Rational,Int>,Rational>("(4x^2+5)y3/2 - 5/3x4y2/3");
    polymake> print $pp;
    (4*x^2 + 5)*y^3/2 + (-5/3*x^4)*y^2/3




Computations
^^^^^^^^^^^^

The standard arithmetic functions “+”, “-”, "\*“,”^" are defined for
polynomials of matching type.


.. link

.. CODE-BLOCK:: perl

    polymake> print $p + ($p^2);
    9*x_1^2 + 6*x_1*x_2^5 + 27*x_1 + x_2^10 + 9*x_2^5 + 20




However, note that due to the fact that their precedence is given in
perl, it may be necessary to write more parentheses than expected at
first sight. For example, as above, you always have to write “($p^2)”
because of the lower precedence of the “^” operator…


.. link

.. CODE-BLOCK:: perl

    polymake> print $p + $p^2;
    36*x_1^2 + 24*x_1*x_2^5 + 96*x_1 + 4*x_2^10 + 32*x_2^5 + 64




For UniPolynomials, we even have polynome division:


.. link

.. CODE-BLOCK:: perl

    polymake> print (($up^2)/$up);
    (2*x^2 + 3*x + 4)/(1)




Example: Newton Polynomials
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is one way to produce polytopes from polynomials (as the convex
hull of the exponent vectors of all terms).


.. link

.. CODE-BLOCK:: perl

    polymake> $np = newton($p*($p+$p));
    polymake> print $np->VERTICES;
    1 0 0 0
    1 0 2 0
    1 0 0 10





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




.. link

.. CODE-BLOCK:: perl

    polymake> print equal_polyhedra($np,minkowski_sum(newton($p),newton($p+$p)));
    true




The Newton polytope of the product of two polynomials always equals the
Minkowski sum of the Newton polytopes of the factors.


Example: Toric Degeneration
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following describes how to construct the polynomial which describes
the toric deformation with respect to a point configuration and a height
function. This is the input data:


.. link

.. CODE-BLOCK:: perl

    polymake> $points = new Matrix<Int>([1,0],[0,1]);
    polymake> $height = new Vector<Int>([2,3]);
    polymake> $coefficients = new Vector<Rational>([-1/2,1/3]);

The following is generic (assuming that the dimensions of the objects
above match).


.. link

.. CODE-BLOCK:: perl

    polymake> $p = new Polynomial($coefficients,$height|$points);

Notice that the points are given in Euclidean coordinates; that is, if
applied, e.g., to the VERTICES of a polytope do not forget to strip the
homogenizing coordinate. The output in our example looks like this:


.. link

.. CODE-BLOCK:: perl

    polymake> print $p;
    1/3*x_0^3*x_2 -1/2*x_0^2*x_1




Puiseux Fractions
-----------------

Polymake supports the usage of Puiseux fractions - see for example `this
paper <https://arxiv.org/abs/1507.08092>`__ for reference.

The preferred way of creating a new Puiseux fraction is to create an
ordinary monomial, and then use that to define a new ``PuiseuxFraction``
object:


.. link

.. CODE-BLOCK:: perl

    polymake> $x = monomials<Rational,Rational>(1); # create a list of `1` monomial, with `Rational` coefficients and `Rational` exponents
    polymake> $f = new PuiseuxFraction<Min>(2*($x^(1/3)) + ($x^(5/2)));

If you have the common denominator of all exponents at hand you could
also intermediately set ``$x = $x^(1/N)`` to save yourself some work.


We can compute the valuation of a puiseux fraction:


.. link

.. CODE-BLOCK:: perl

    polymake> print $f->val;
    1/3




Evaluate a puiseux fraction at `2^6`:


.. link

.. CODE-BLOCK:: perl

    polymake> print $f->evaluate(2,6);
    32776




Operators like ``+``, ``-``, ``*``, ``/`` are defined as you’d expect.

Besides, puiseux fractions, similar to rational functions over any
ordered field, have a natural ordering induced by the ordering of the
coefficients (see the above mentioned paper for detals) - polymake
correspondingly overloads the operators ``<``, ``>``, ``<=``, ``>=``:


.. link

.. CODE-BLOCK:: perl

    polymake> $g = new PuiseuxFraction<Min>(3*($x^(3/2)));
    polymake> print $f>$g;
    true




Applications
^^^^^^^^^^^^


One usage example is parametrized polyhedra. As an example we compute a
family of 3 dimensional Klee-Minty cubes:


.. link

.. CODE-BLOCK:: perl

    polymake> $k = klee_minty_cube(3, $f);
    polymake> print "facets:\n", $k->FACETS, "\nvolume:\n", $k->VOLUME;
    facets:
    (0) (1) (0) (0)
    (1) (- 1) (0) (0)
    (0) (-2*x^1/3 - x^5/2) (1) (0)
    (1) (-2*x^1/3 - x^5/2) (- 1) (0)
    (0) (0) (-2*x^1/3 - x^5/2) (1)
    (1) (0) (-2*x^1/3 - x^5/2) (- 1)
    
    volume:
    (1 -4*x^1/3 + 4*x^2/3 -2*x^5/2 + 4*x^17/6 + x^5)




You can even check for (combinatorial) isomorphy:


.. link

.. CODE-BLOCK:: perl

    polymake> print isomorphic($k, cube(3));
    true




.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package nauty
      Computation of automorphism groups of graphs.
      Copyright by Brendan McKay and Adolfo Piperno.
      http://pallini.di.uniroma1.it/
    </pre>
    </details>




As another example related to linear optimization we compute a family of
3 dimensional Goldfarb-Sit cubes (again, see the above mentioned paper,
and consult:


.. link

.. CODE-BLOCK:: perl

    polymake> $l = goldfarb_sit(3, $g, 1/2);
    polymake> print $l->LP->MAXIMAL_VALUE;
    (1)




.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package tosimplex
      Dual simplex algorithm implemented by Thomas Opfer
    
    </pre>
    </details>




.. link

.. CODE-BLOCK:: perl

    polymake> print $l->LP->MAXIMAL_VERTEX;
    (1) (0) (0) (1)




.. link

.. CODE-BLOCK:: perl

    polymake> print $l->VOLUME;
    (27/8*x^9/2 -81/4*x^6 + 243/8*x^15/2)


