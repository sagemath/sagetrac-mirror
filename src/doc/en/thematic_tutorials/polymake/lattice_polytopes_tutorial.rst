.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Tutorial for Lattice Polytopes
==============================

This page gives a small introduction to lattice polytopes in
``polymake``, some useful external software, and usage hints for it. For
a list of methods and properties applicable to lattice polytopes see
`here <:user_guide:lattice_polytopes_doc>`__. For an introduction to the
``polymake`` package see `here <:user_guide:start>`__.

``polymake`` always assumes that the lattice used to define a lattice
polytope is the standard lattice Zd. Some rules also require that the
polytope is full dimensional. There are user functions that transform a
polytope sitting in some affine subspace of Rd into a full dimensional
polytope, either in the induced lattice or the lattice spanned by the
vertices, see below.

Dependence on other Software
----------------------------

For some computations ``polymake`` has no built-in commands and passes
the computation to external software. Currently, polymake has an
interface to the following packages that compute various properties of
lattice polytopes.

-  `libnormaliz <http://www.math.uos.de/normaliz/>`__ by Winfried Bruns
   and Bogdan Ichim, bundled with polymake

-  `4ti2 <http://www.4ti2.de/>`__ by the 4ti2 team

-  `LattE macchiato <http://www.math.ucdavis.edu/~mkoeppe/latte/>`__ by
   Matthias Köppe, building on ``LattE`` by Jesus de Loera et. al.

-  (`barvinok <http://freshmeat.net/projects/barvinok>`__ by Sven
   Verdoolaege)

Unless you want to deal with Hilbert bases of cones you don’t need them.
If you do, either the bundled extension ``libnormaliz`` or the external
package ``4ti2`` suffices to do most computations with lattice
polytopes. Computation of Gröbner bases currently requires ``4ti2``.
``LattE`` only counts lattice points in a polytope and computes its
Ehrhart polynomial, but may be faster on that than any other methods
implemented. ``barvinok`` can be used to compute the number of lattice
points and the h-polynomial. Access to barvinok is realized via an
extension which has to be downloaded separately.

For some of the commands in this tutorial you will need at least one of
``bundled:libnormaliz`` enabled or ``4ti2`` installed on your machine.
We’ll remind you at the relevant places.

Lattice Points in Rational Polytopes
------------------------------------

We start by creating a rational polytope using one of ``polymake``\ ’s
standard polytope constructions. We choose the 3-dimensional cube with
coordinates +1 and -1. So we start ``polymake`` at the command line and
assign a cube to the variable $p.


.. link

.. CODE-BLOCK:: perl

    polymake> $p=cube(3);

Suppose we want to know how many lattice points this cube contains. The
answer is of course already known, as the cube has one relative interior
integral point per non-empty face. So we expect to get the answer 27.


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->N_LATTICE_POINTS;
    27




To satisfy this request, ``polymake`` computes all properties necessary
to call an external program that provides the number of lattice points.
In this case, ``polymake`` has passed the request to ``lattE``, which is
shown by the credit message that appears before the answer. By default,
credits for external software are shown when an external package is used
for the first time. You can change this behavior using the variable
``$Verbose::credits``. If you don’t have a version of ``LattE``, or if
you have set different preferences, then ``polymake`` may choose one of
the other programs. So the credit statement depends on your
configuration.

We can of course also ask ``polymake`` to compute the integral points
for us. For our next computations we are only interested in the integral
points in the interior of the cube, so we ask for


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->INTERIOR_LATTICE_POINTS;
    1 0 0 0





Internally, ``polymake`` computes the intersection of the polytope with
the integer lattice, and then checks which of the points lies on a facet
of $p. By default, ``polymake`` uses a project-and-lift algorithms to
enumerate the lattice points. Note that our call to ``LattE`` above has
only computed the number of integral points (which is done with an
improved version of Barvinok’s algorithm), so ``polymake`` really has to
compute something here. If we had asked for ``INTERIOR_LATTICE_POINTS``
first, then ``N_LATTICE_POINTS`` would just have counted the rows of a
matrix, which would have been much faster. So computation time can
depend on the history.

You can also ask for the HILBERT_BASIS, though in the case of a cube the
result is not so exciting:


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->HILBERT_BASIS;
    1 -1 -1 -1
    1 -1 -1 0
    1 -1 -1 1
    1 -1 0 -1
    1 -1 0 0
    1 -1 0 1
    1 -1 1 -1
    1 -1 1 0
    1 -1 1 1
    1 0 -1 -1
    1 0 -1 0
    1 0 -1 1
    1 0 0 -1
    1 0 0 0
    1 0 0 1
    1 0 1 -1
    1 0 1 0
    1 0 1 1
    1 1 -1 -1
    1 1 -1 0
    1 1 -1 1
    1 1 0 -1
    1 1 0 0
    1 1 0 1
    1 1 1 -1
    1 1 1 0
    1 1 1 1





.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package libnormaliz
      Normaliz is a tool for computations in affine monoids, vector configurations, lattice polytopes, and rational cones.
      Copyright by Winfried Bruns, Bogdan Ichim, Christof Soeger.
      http://www.math.uos.de/normaliz/
    
    </pre>
    </details>




``polymake`` has no native method to compute a Hilbert basis, so it has
passed the computation to ``4ti2``. The choice may vary, depending on
what is installed on your computer (and configured for ``polymake``).
You can influence the choice with the appropriate ``prefer`` statement.

Note that so far these commands also work for rational polytopes.

Lattice Polytopes
-----------------

Now we want to do some computations that don’t make sense for polytopes
that have non-integral vertex coordinates. We can let ``polymake`` check
that our cube is indeed a polytope with integral vertices.


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->LATTICE;
    true




A particularly interesting class of lattice polytopes is that of
reflexive polytopes. A polytope is *reflexive* if its polar is agein
alattice polytope. This implies in particular that the origin is the
unique interior lattice point in the polytope. So, as we have seen
above, our cube is a candidate. But this is not sufficient, so we have
to do further checks.

Reflexivity is a property that is not defined for polytopes with
non-integral vertices. So if we ask for it in ``polymake``, then
``polymake`` checks that the entered polytope is indeed a lattice
polytope (i.e. it is **bounded** and has **integral vertices**). In that
case the object will automatically get the specialization
``Polytope::Lattice``.


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->REFLEXIVE;
    true




Lattice polytopes can be used to define toric varieties with an ample
line bundle, and many properties of the variety are reflected by the
polytope. here is an example: The toric variety defined by our cube is
*smooth*, i.e. it is one of the *smooth toric Fano varieties*. In
``polymake``, we can just ask for this property in the following way.


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->SMOOTH;
    true




The number of integral points in the k-th dilate of a polytope is given
by a polynomial of degree d in k. This is the famous *Ehrhart Theorem*.
In ``polymake`` you can obtain the coefficients of this polynomial
(starting with the constant coefficient).


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->EHRHART_POLYNOMIAL;
    8*x^3 + 12*x^2 + 6*x + 1




``polymake`` has passed this request to ``LattE`` or ``normaliz``, but
as we have used these programs already the credit message is suppressed
(but if you save the cube to a file, then you will find it in there).
Some coefficients of this polynomial have a geometric interpretation.
E.g., the highest coefficient is the Euclidean volume of the polytope.


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->VOLUME;
    8




By a theorem of Stanley, the generating function for the number of
lattice points can be written as the quotient of a polynomial h(t) by
(1-t)d+1, and this polynomial has non-negative integral coefficients.


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->H_STAR_VECTOR;
    1 23 23 1




.. link

.. CODE-BLOCK:: perl

    polymake> print $p->LATTICE_DEGREE;
    3




.. link

.. CODE-BLOCK:: perl

    polymake> print $p->LATTICE_CODEGREE;
    1




In our case the coefficient vector is symmetric, as the polytope is
reflexive. The *co-degree* of the polytope is d+1 minus the degree of
the h-polynomial. It is the smallest factor by which we have to dilate
the polytope to obtain an interior integral point. In our case, this is
1, as the cube already has an integral point.

We can obtain the volume of our polytope also from the
``H_STAR_VECTOR``: Summing up the coefficients give the *lattice volume*
of the polytope, which is d! times its Euclidean volume.


.. link

.. CODE-BLOCK:: perl

    polymake> print $p->LATTICE_VOLUME;
    48




Let us look at a different example:


.. link

.. CODE-BLOCK:: perl

    polymake> $q=new Polytope(INEQUALITIES=>[[5,-4,0,1],[-3,0,-4,1],[-2,1,0,0],[-4,4,4,-1],[0,0,1,0],[8,0,0,-1],[1,0,-1,0],[3,-1,0,0]]);

This actually defines a lattice polytope, which we can see from the list
of vertices:


.. link

.. CODE-BLOCK:: perl

    polymake> print $q->VERTICES;
    1 3 1 7
    1 2 0 3
    1 3 0 7
    1 2 1 7
    1 2 0 4
    1 3 1 8
    1 3 0 8
    1 2 1 8





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




``polymake`` provides basically three methods for convex hull
conversion, double description, reverse search, and beneath beyond. The
first two are provided by the packages ``cdd`` and ``lrs``, the last in
internal. By default, ``cdd`` is chosen, and that is what was used above
(they are bundled with ``polymake``, you don’t have to install them). A
polytope Q is *normal* if every lattice point in the k-th dilate of Q is
the sum of k lattice points in Q. You can check this property via


.. link

.. CODE-BLOCK:: perl

    polymake> print $q->NORMAL;
    false




So our polytope is not normal. We can also find a point that violates
the condition. Being normal is equivalent to the fact, that the Hilbert
basis of the cone C(Q) obtained from Q by embedding the polytope at
height one and the coning over it has all its generators in height one.
The property HILBERT_BASIS computes these generators:


.. link

.. CODE-BLOCK:: perl

    polymake> print $q->HILBERT_BASIS;
    1 2 0 3
    1 2 0 4
    1 2 1 7
    1 2 1 8
    1 3 0 7
    1 3 0 8
    1 3 1 7
    1 3 1 8
    2 5 1 13





The last row is the desired vector: [2,5,1,13] is a vector in 2*Q, but
it is not a sum of lattice points in Q. The cone C(Q) corresponds to an
affine toric variety, and the above tells us that this variety is not
normal. Yet, it is very ample, as we can check with


.. link

.. CODE-BLOCK:: perl

    polymake> print $q->VERY_AMPLE;
    true




Now assume we are particularly interested in the third facet of Q. We
can pick this via


.. link

.. CODE-BLOCK:: perl

    polymake> $f=facet($q,2);

Recall that indexes in ``polymake`` start at 0, so the third facet has
index 2. This is again a very ample polytope:


.. link

.. CODE-BLOCK:: perl

    polymake> print $f->VERY_AMPLE;
    true




The result is no surprise, being very ample is inherited by faces. We
could also be interested in the facet width of the polytope $f. This is
the minimum over the maximal distance of a facet to any other vertex.
``polymake`` knows how to compute this:


.. link

.. CODE-BLOCK:: perl

    polymake> #print $f->FACET_WIDTH;

Almost. It tells you that it can only do this for a full dimensional
polytope, i.e. for a polytope whose dimension coincides with the ambient
dimension. This is not true for our facet: It lives in the same ambient
space as $q, but has one dimension less. We can remedy this by applying
the following:


.. link

.. CODE-BLOCK:: perl

    polymake> $g=ambient_lattice_normalization($f);
    polymake> print $g->FACET_WIDTH;
    1




The function ``ambient_lattice_normalization`` returns a full
dimensional version of the polytope ``$f`` in the lattice induced by the
intersection of the affine space of ``$f`` with Z^n. Now ``$g`` is full
dimensional, and we can compute the facet width. Note that there is also
a function which normalizes in the lattice spanned by the vertices of
the polytope: ``vertex_facet_normalization``. This can also be usefull
for full dimensional polytopes. E.g. consider the cube we defined above.
The sum of the entries of each vertex is odd, so the lattice spannd by
the vertices is a sublattice of the integer lattice:


.. link

.. CODE-BLOCK:: perl

    polymake> $cr=vertex_lattice_normalization($p);
    polymake> print $cr->VERTICES;
    (4) (0 1)
    1 1 0 0
    1 0 1 0
    1 1 1 0
    1 0 0 1
    1 1 0 1
    1 0 1 1
    1 1 1 1





``$cr`` is the same cube, but we have reduced the lattice. (The first
line is a *sparse representation* of a vector: it has length 4, and the
only non-zero entry is at position 0 and is 1 (note that indexes start
at 0)).

Toric Varieties
---------------

``polymake`` has only few builtin functions to compute properties of the
variety associated to a fan or lattice polytope. There are two
extensions available that add more properties, both currently at an
early stage:

-  `Toric Varieties and Singular
   interface <https///github.com/lkastner>`__ by Lars Kastner/Benjamin
   Lorenz

-  `ToricVarieties-v0.3 <http://www.mathematik.tu-darmstadt.de/~paffenholz/software.html>`__
   by Andreas Paffenholz. Defines a new property for toric varieties
   associated to a fan and divisors on that variety.

Here we will do some computations that do not require one of the
extensions. We start by defining a fan. We’ll make our live easy and
take the normal fan of our cube:


.. link

.. CODE-BLOCK:: perl

    polymake> application "fan";

.. link

.. CODE-BLOCK:: perl

    polymake> $f = normal_fan($p);
    polymake> print $f->SMOOTH_FAN;
    true




With the last line we have verified that our fan defines a smooth toric
variety. Note that switching the application is not strictly necessary,
you can also prepend calls to functions and constructors with ``fan::``.
The fan object $f itself knows its type, and chooses available
properties based on this. Any smooth variety is Gorenstein, so we expect
the following:


.. link

.. CODE-BLOCK:: perl

    polymake> print $f->GORENSTEIN;
    true




Similarly, we could check for Q-Gorensteinness with ``Q_GORENSTEIN``. It
is also a complete fan:


.. link

.. CODE-BLOCK:: perl

    polymake> print $f->COMPLETE;
    true




but currently there is little support to detect completeness in
``polymake``. In our case it was already decided during construction,
normal fans are complete. You can also check standard features of fans,
like their rays. Let us do this for the normal fan of our other example:


.. link

.. CODE-BLOCK:: perl

    polymake> $g=normal_fan($q);
    polymake> print $g->RAYS;
    -1 0 1/4
    0 -1 1/4
    1 0 0
    1 1 -1/4
    0 1 0
    0 0 -1
    0 -1 0
    -1 0 0





This is not what we wanted. We would like to see the minimal lattice
generators of the rays. We can fix this using


.. link

.. CODE-BLOCK:: perl

    polymake> print primitive($g->RAYS);
    -4 0 1
    0 -4 1
    1 0 0
    4 4 -1
    0 1 0
    0 0 -1
    0 -1 0
    -1 0 0





Note that the function ``primitive`` returns a copy of the argument, the
RAYS as stored in the fan are unchanged. So you have to apply this
function each time you need the primitive generators, or you store them
in a new variable. The fan `g` is not smooth, but still
Gorenstein:


.. link

.. CODE-BLOCK:: perl

    polymake> print $g->SMOOTH_FAN;
    false




.. link

.. CODE-BLOCK:: perl

    polymake> print $g->GORENSTEIN;
    true




You can also access the maximal cones of the fan via


.. link

.. CODE-BLOCK:: perl

    polymake> print $g->MAXIMAL_CONES;
    {0 1 6 7}
    {0 1 2 4}
    {0 4 7}
    {1 2 6}
    {2 3 4}
    {5 6 7}
    {3 4 5 7}
    {2 3 5 6}





The indices in these list refer to the list of rays. Sometimes you might
be interested in the walls, i.e. the codimension 2 faces of the fan.
Here is one way to get them


.. link

.. CODE-BLOCK:: perl

    polymake> print rows_numbered($g->HASSE_DIAGRAM->FACES);
    0:-1
    1:0 1 6 7
    2:0 1 2 4
    3:0 4 7
    4:1 2 6
    5:2 3 4
    6:5 6 7
    7:3 4 5 7
    8:2 3 5 6
    9:0 1
    10:0 7
    11:1 6
    12:6 7
    13:0 4
    14:1 2
    15:2 4
    16:4 7
    17:2 6
    18:3 4
    19:2 3
    20:5 7
    21:5 6
    22:3 5
    23:0
    24:1
    25:7
    26:6
    27:4
    28:2
    29:3
    30:5
    31:





.. link

.. CODE-BLOCK:: perl

    polymake> print $g->HASSE_DIAGRAM->nodes_of_dim($g->DIM-2);
    {23 24 25 26 27 28 29 30}




where the list of numbers given by the latter are the indices of the
codimension 2 faces in the list of all faces given before. There is a
more concise way to list those, using some simple perl programming:


.. link

.. CODE-BLOCK:: perl

    polymake> print map($g->HASSE_DIAGRAM->FACES->[$_], @{$g->HASSE_DIAGRAM->nodes_of_dim($g->DIM-2)});
    {0}{1}{7}{6}{4}{2}{3}{5}




Visualization
-------------

If the lattice polytope lives in R^2 or R^3, then we can visualize the
polytope together with its lattice points. The picture below has been
made with ```javaview`` <http://www.javaview.de>`__, but
``polymake``\ ’s standard visualization method is now ``jreality``,
which is bundled with ``polymake``.


.. link

.. CODE-BLOCK:: perl

    polymake> $p->VISUAL->LATTICE_COLORED;


.. raw:: html

    <!--
    polymake for andrew
    Thu Apr  4 10:50:43 2019
    p
    -->
    
    
    <html>
       <head>
          <title>p</title>
          <style>
    /*
    // COMMON_CODE_BLOCK_BEGIN
    */
             html{overflow: scroll;}
             body { font-family: Arial, Helvetica, sans-serif}
             strong{font-size: 18px;}
             canvas { z-index: 8; }
             input[type='range'] {}
             input[type='radio'] {margin-left:0;}
             input[type='checkbox'] {margin-right:7px; margin-left: 0px; padding-left:0px;}
             .group{padding-bottom: 40px;}
             .settings * {z-index: 11; }
             .settings{z-index: 10; margin-left: 30px; display: none; width: 14em; height: 90%; border: solid 1px silver; padding: 2px; overflow-y: scroll; background-color: white }
             .indented{margin-left: 20px; margin-top: 15px; padding-bottom: 0px;} 
             .shownObjectsList{overflow: auto; max-width: 150px; max-height: 150px;}
             .showSettingsButton{display: block; z-index: 12; position: absolute }
             .hideSettingsButton{display: none; z-index: 12; position: absolute; opacity: 0.5}
             .resetButton{margin-top: 20px;}
             button{margin-left: 0;}
             img{cursor: pointer;}
             .suboption{padding-top: 30px;}
             .transparency{display: none;}
             .labelsCheckbox{margin-top: 10px;}
    
    
             input[type=range] {
               -webkit-appearance: none;
               padding:0; 
               width:90%; 
               margin-left: auto;
               margin-right: auto;
               margin-top: 20px;
               display: block;	
             }
             input[type=range]:focus {
               outline: none;
             }
             input[type=range]::-webkit-slider-runnable-track {
               height: 4px;
               cursor: pointer;
               animate: 0.2s;
               box-shadow: 0px 0px 0px #000000;
               background: #E3E3E3;
               border-radius: 0px;
               border: 0px solid #000000;
             }
             input[type=range]::-webkit-slider-thumb {
               box-shadow: 1px 1px 2px #B8B8B8;
               border: 1px solid #ABABAB;
               height: 13px;
               width: 25px;
               border-radius: 20px;
               background: #E0E0E0;
               cursor: pointer;
               -webkit-appearance: none;
               margin-top: -5px;
             }
             input[type=range]:focus::-webkit-slider-runnable-track {
               background: #E3E3E3;
             }
             input[type=range]::-moz-range-track {
               height: 4px;
               cursor: pointer;
               animate: 0.2s;
               box-shadow: 0px 0px 0px #000000;
               background: #E3E3E3;
               border-radius: 0px;
               border: 0px solid #000000;
             }
             input[type=range]::-moz-range-thumb {
               box-shadow: 1px 1px 2px #B8B8B8;
               border: 1px solid #ABABAB;
               height: 13px;
               width: 25px;
               border-radius: 20px;
               background: #E0E0E0;
               cursor: pointer;
             }
             input[type=range]::-ms-track {
               height: 4px;
               cursor: pointer;
               animate: 0.2s;
               background: transparent;
               border-color: transparent;
               color: transparent;
             }
             input[type=range]::-ms-fill-lower {
               background: #E3E3E3;
               border: 0px solid #000000;
               border-radius: 0px;
               box-shadow: 0px 0px 0px #000000;
             }
             input[type=range]::-ms-fill-upper {
               background: #E3E3E3;
               border: 0px solid #000000;
               border-radius: 0px;
               box-shadow: 0px 0px 0px #000000;
             }
             input[type=range]::-ms-thumb {
               box-shadow: 1px 1px 2px #B8B8B8;
               border: 1px solid #ABABAB;
               height: 13px;
               width: 25px;
               border-radius: 20px;
               background: #E0E0E0;
               cursor: pointer;
             }
             input[type=range]:focus::-ms-fill-lower {
               background: #E3E3E3;
             }
             input[type=range]:focus::-ms-fill-upper {
               background: #E3E3E3;
             }
    /*
    // COMMON_CODE_BLOCK_END
    */
    		</style>
       </head>
    
    <body>
    
    		<div id='settings_0' class='settings'>
    			<div class=group id='explode_0'>
    				<strong>Explode</strong>
    				<input id='explodeRange_0' type='range' min=0 max=6 step=0.01 value=0>
    				<div class=indented><input id='explodeCheckbox_0' type='checkbox'>Automatic explosion</div>
    				<div class=suboption>Exploding speed</div>
    				<input id='explodingSpeedRange_0' type='range' min=0 max=0.5 step=0.001 value=0.05>
    			</div>
    
    			
    			<div class=group id='transparency_0' class='transparency'>
    				<strong>Transparency</strong>
    				<input id='transparencyRange_0' type='range' min=0 max=1 step=0.01 value=0>
    			</div>
    			
    			<div class=group id='rotation_0'>
    				<strong>Rotation</strong>
    				<div class=indented>
    					<div><input type='checkbox' id='changeRotationX_0'> x-axis</div>
    					<div><input type='checkbox' id='changeRotationY_0'> y-axis</div>
    					<div><input type='checkbox' id='changeRotationZ_0'> z-axis</div>
    					<button id='resetButton_0' class='resetButton' >Reset</button>
    				</div>
    
    				<div class=suboption>Rotation speed</div>
    				<input id='rotationSpeedRange_0' type='range' min=0 max=5 step=0.01 value=2>
    
    			</div>
    
    
    			<div class=group id='display_0'>
    				<strong>Display</strong>
    				<div class=indented>
    					<div id='shownObjectsList_0' class='shownObjectsList'></div>
    					<div class='labelsCheckbox'><input type='checkbox' id='labelsCheckboxInput_0' checked>Labels</div>
    				</div>
    			</div>
    
    
    			<div class=group id='svg_0'>
    				<strong>SVG</strong>
    				<div class=indented>
    					<form>
    						<input type="radio" name='screenshotMode' value='download' id='download_0' checked> Download<br>
    						<input type="radio" name='screenshotMode' value='tab' id='tab_0' > New tab<br>
    					</form>
    					<button id='takeScreenshot_0'>Screenshot</button>
    				</div>
    			</div>
    
    		</div>	<!-- end of settings -->
    		<img id='hideSettingsButton_0' style="display: none" class='hideSettingsButton' src='/kernelspecs/polymake/close.svg' width=20px">
    		<img id='showSettingsButton_0' class='showSettingsButton' src='/kernelspecs/polymake/menu.svg' width=20px">
    <div id="model50695986238"></div>
    
    <script>
    requirejs.config({
      paths: {
        three: '/kernelspecs/polymake/three',
        Detector: '/kernelspecs/polymake/Detector',
        SVGRenderer: '/kernelspecs/polymake/SVGRenderer',
        CanvasRenderer: '/kernelspecs/polymake/CanvasRenderer',
        Projector: '/kernelspecs/polymake/Projector',
        TrackballControls: '/kernelspecs/polymake/TrackballControls'
      },
      shim: {
        'three':
        {
          exports: 'THREE'
        },
        'Detector':
        {
          deps: [ 'three' ],
          exports: 'Detector'
        },
        'SVGRenderer':
        {
          deps: [ 'three' ],
          exports: 'THREE.SVGRenderer'
        },
        'CanvasRenderer':
        {
          deps: [ 'three' ],
          exports: 'THREE.CanvasRenderer'
        },
        'Projector':
        {
          deps: [ 'three' ],
          exports: 'THREE.Projector'
        },
        'TrackballControls':
        {
          deps: [ 'three' ],
          exports: 'THREE.TrackballControls'
        }
      }
    });
    require(['three'],function(THREE){
        window.THREE = THREE;
      require(['Detector','SVGRenderer','CanvasRenderer','Projector','TrackballControls'],function(Detector,SVGRenderer,CanvasRenderer,Projector,TrackballControls){
          THREE.SVGRenderer = SVGRenderer;
          THREE.CanvasRenderer = CanvasRenderer;
          THREE.Projector = Projector;
          THREE.TrackballControls = TrackballControls;
    
    // COMMON_CODE_BLOCK_BEGIN
    	var foldable = false;
       var container = document.getElementById( 'model50695986238' );
       var renderer = Detector.webgl? new THREE.WebGLRenderer({antialias: true}): new THREE.CanvasRenderer({antialias: true});
    	var svgRenderer = new THREE.SVGRenderer({antialias: true});
                var box = document.getElementsByClassName( 'output_subarea' )[0];
             var notebook = document.getElementById( 'notebook_panel' );
    
       var width = box.clientWidth - 25;
       var height = notebook.clientHeight * 0.8;
       renderer.setSize(width, height);
       svgRenderer.setSize(width, height);
       renderer.setClearColor(0xFFFFFF, 1);
       svgRenderer.setClearColor(0xFFFFFF, 1);
    
       container.appendChild(renderer.domElement);
    
       var scene = new THREE.Scene();
       var camera = new THREE.PerspectiveCamera(75, width/height, 0.1, 1000);
    
       var renderid;
    
       camera.position.set(0, 0, 5);
       camera.lookAt(0, 0, 0);
       camera.up.set(0, 1, 0);
    
       // class to allow move points together with labels and spheres
       var PMPoint = function (x,y,z) {
          this.vector = new THREE.Vector3(x,y,z);
          this.sprite = null;
          this.sphere = null;
       }
       PMPoint.prototype.makelabel = function(label) {
          this.sprite = textSprite( label );
          this.sprite.position.copy(this.vector);
       }
       PMPoint.prototype.makesphere = function(radius,material) {
          this.sphere = new THREE.Mesh(new THREE.SphereGeometry(radius), material);
          this.sphere.position.copy(this.vector);
       }
    
       PMPoint.prototype.setX = function(x) {
          this.vector.setX(x);
          if (this.sprite) {
             this.sprite.position.setX(x);
          }
          if (this.sphere) {
             this.sphere.position.setX(x);
          }
       };
       PMPoint.prototype.setY = function(y) {
          this.vector.setY(y);
          if (this.sprite) {
             this.sprite.position.setY(y);
          }
          if (this.sphere) {
             this.sphere.position.setY(y);
          }
       };
       PMPoint.prototype.setZ = function(z) {
          this.vector.setZ(z);
          if (this.sprite) {
             this.sprite.position.setZ(z);
          }
          if (this.sphere) {
             this.sphere.position.setZ(z);
          }
       };
       PMPoint.prototype.set = function(x,y,z) {
          this.vector.set(x,y,z);
          if (this.sprite) {
             this.sprite.position.set(x,y,z);
          }
          if (this.sphere) {
             this.sphere.position.set(x,y,z);
          }
       };
       PMPoint.prototype.add = function(o) {
          if (this.sprite) {
             o.add(this.sprite);
          }
          if (this.sphere) {
             o.add(this.sphere);
          }
       };
    
    
       var controls = new THREE.TrackballControls(camera, container);
    	controls.zoomSpeed = 0.2;
    	controls.rotateSpeed = 4;
    
       var all_objects = [];
       var centroids = [];
       // select the target node
       var target = document.querySelector('#model50695986238');
    
       // create an observer instance
       var observer = new MutationObserver(function(mutations) {
          mutations.forEach(function(mutation) {
             if (mutation.removedNodes && mutation.removedNodes.length > 0) {
                cancelAnimationFrame(renderId);
                observer.disconnect();
                console.log("cancelled frame "+renderId);
             }
          });
       });
    
       // configuration of the observer:
       var config = { childList: true, characterData: true }
    
       // pass in the target node, as well as the observer options
       while (target) {
          if (target.className=="output") {
             observer.observe(target, config);
             break;
          }
          target = target.parentNode;
       }
    
    // COMMON_CODE_BLOCK_END
    
       var objectnames = ["p","Lattice points and vertices of p"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(-1, -1, -1));
       allpoints.push(new PMPoint(1, -1, -1));
       allpoints.push(new PMPoint(-1, 1, -1));
       allpoints.push(new PMPoint(1, 1, -1));
       allpoints.push(new PMPoint(-1, -1, 1));
       allpoints.push(new PMPoint(1, -1, 1));
       allpoints.push(new PMPoint(-1, 1, 1));
       allpoints.push(new PMPoint(1, 1, 1));
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       <!-- Edge style -->
       var line_material = new THREE.LineBasicMaterial ( {color: 0x000000, linewidth: 1.5, } );
    
       line_material.side = THREE.DoubleSide;
       line_material.transparent = true;
    
       <!-- EDGES --> 
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[7].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[6].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(-1, -1, -1));
       allpoints.push(new PMPoint(-1, -1, 0));
       allpoints.push(new PMPoint(-1, -1, 1));
       allpoints.push(new PMPoint(-1, 0, -1));
       allpoints.push(new PMPoint(-1, 0, 0));
       allpoints.push(new PMPoint(-1, 0, 1));
       allpoints.push(new PMPoint(-1, 1, -1));
       allpoints.push(new PMPoint(-1, 1, 0));
       allpoints.push(new PMPoint(-1, 1, 1));
       allpoints.push(new PMPoint(0, -1, -1));
       allpoints.push(new PMPoint(0, -1, 0));
       allpoints.push(new PMPoint(0, -1, 1));
       allpoints.push(new PMPoint(0, 0, -1));
       allpoints.push(new PMPoint(0, 0, 1));
       allpoints.push(new PMPoint(0, 1, -1));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(0, 1, 1));
       allpoints.push(new PMPoint(1, -1, -1));
       allpoints.push(new PMPoint(1, -1, 0));
       allpoints.push(new PMPoint(1, -1, 1));
       allpoints.push(new PMPoint(1, 0, -1));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(1, 0, 1));
       allpoints.push(new PMPoint(1, 1, -1));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(1, 1, 1));
    
       <!-- Vertex style -->
       var materials = [
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
       ];
       for (index = 0; index < materials.length; ++index) {
          materials[index].side = THREE.DoubleSide;
       }
       var points_material = new THREE.MeshFaceMaterial ( materials );
    
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,materials[0]);
       allpoints[1].makesphere(0.02,materials[1]);
       allpoints[2].makesphere(0.02,materials[2]);
       allpoints[3].makesphere(0.02,materials[3]);
       allpoints[4].makesphere(0.02,materials[4]);
       allpoints[5].makesphere(0.02,materials[5]);
       allpoints[6].makesphere(0.02,materials[6]);
       allpoints[7].makesphere(0.02,materials[7]);
       allpoints[8].makesphere(0.02,materials[8]);
       allpoints[9].makesphere(0.02,materials[9]);
       allpoints[10].makesphere(0.02,materials[10]);
       allpoints[11].makesphere(0.02,materials[11]);
       allpoints[12].makesphere(0.02,materials[12]);
       allpoints[13].makesphere(0.02,materials[13]);
       allpoints[14].makesphere(0.02,materials[14]);
       allpoints[15].makesphere(0.02,materials[15]);
       allpoints[16].makesphere(0.02,materials[16]);
       allpoints[17].makesphere(0.02,materials[17]);
       allpoints[18].makesphere(0.02,materials[18]);
       allpoints[19].makesphere(0.02,materials[19]);
       allpoints[20].makesphere(0.02,materials[20]);
       allpoints[21].makesphere(0.02,materials[21]);
       allpoints[22].makesphere(0.02,materials[22]);
       allpoints[23].makesphere(0.02,materials[23]);
       allpoints[24].makesphere(0.02,materials[24]);
       allpoints[25].makesphere(0.02,materials[25]);
       allpoints[26].makesphere(0.02,materials[26]);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       scene.add(obj);
       all_objects.push(obj);
    
    // COMMON_CODE_BLOCK_BEGIN
    var xRotationEnabled = false;
    var yRotationEnabled = false;
    var zRotationEnabled = false;
    var rotationSpeedFactor = 1;
    var settingsShown = false;
    var labelsShown = true;
    var intervals = [];
    var timeouts = [];
    var explodingSpeed = 0.05;
    var explodeScale = 0;
    var XMLS = new XMLSerializer();
    var svgElement;
    var renderId;
    
    	var render = function () {
    
    		renderId = requestAnimationFrame(render);
    
    //		comment in for automatic explosion
    //		explode(updateFactor());
    
    		var phi = 0.02 * rotationSpeedFactor;
    
    		if (xRotationEnabled){
    			scene.rotation.x += phi;
    		}
    		if(yRotationEnabled){
    			scene.rotation.y += phi;
    		}
    		if(zRotationEnabled){
    			scene.rotation.z += phi;
    		}
    
    		controls.update();
    		renderer.render(scene, camera);
    	};
    
    	render();
    
    	function computeCentroid(geom) {
    		centroid = new THREE.Vector3();
    		geom.vertices.forEach(function(v) {
    			centroid.add(v);			
    		});
    		centroid.divideScalar(geom.vertices.length);
    		return centroid;
    	}
    
    	function changeTransparency(event){
    		var opacity = 1-Number(event.currentTarget.value);
    		for (var i=0; i<all_objects.length; i++){
    			for (var j=0; j<all_objects[i].children.length; j++){
    				if (all_objects[i].children[j].material.type == "MultiMaterial") {
    					for (var k=0; k<all_objects[i].children[j].material.materials.length; k++){
    						all_objects[i].children[j].material.materials[k].opacity = opacity;
    						all_objects[i].children[j].material.materials[k].depthWrite = opacity < 0.5 ? false : true;
    						all_objects[i].children[j].material.materials[k].depthTest = opacity < 0.5 ? false : true;
    					}
    				} else if (all_objects[i].children[j].material.transparent && 
    							  all_objects[i].children[j].material.type == "MeshBasicMaterial" &&
    							  all_objects[i].children[j].geometry.type == "Geometry"){
    					all_objects[i].children[j].material.opacity = opacity;
    					all_objects[i].children[j].material.depthWrite = opacity < 0.5 ? false : true;
    					all_objects[i].children[j].material.depthTest = opacity < 0.5 ? false : true;
    				}
    			}
    		}
    	}
    
    	function changeRotationX(event){
    		xRotationEnabled = event.currentTarget.checked;
    	}	
    
    	function changeRotationY(event){
    		yRotationEnabled = event.currentTarget.checked;
    	}	
    
    	function changeRotationZ(event){
    		zRotationEnabled = event.currentTarget.checked;
    	}	
    
    
    	function changeRotationSpeedFactor(event){
    		rotationSpeedFactor = Number(event.currentTarget.value);
    	}
    
    	function resetScene(){
    		scene.rotation.set(0,0,0);
    		camera.position.set(0,0,5);
    		camera.up.set(0,1,0);
    	}
    
    	function showSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_0').style.position = 'absolute';
    		document.getElementById('settings_0').style.display = 'block';
    		document.getElementById('showSettingsButton_0').style.display = 'none';
    		document.getElementById('hideSettingsButton_0').style.display = 'block';
    		settingsShown = true;
    	}
    
    	function hideSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_0').style.display = 'none';
    		document.getElementById('hideSettingsButton_0').style.display = 'none';
    		document.getElementById('showSettingsButton_0').style.display = 'block';
    		settingsShown = false;
    	}
    
    
    
    	var pos = 150* Math.PI;
    
    	function updateFactor() {
    		pos++;
    		return Math.sin(.01*pos)+1;
    	}
    
    	function makelabel(message, x, y, z, params) {
    		var spritey = textSprite( message, params );
    		spritey.position.set(x, y, z);
    		obj.add(spritey);
    	}
    
    	function textSprite(message, parameters)
    	{
    		if ( parameters === undefined ) parameters = {};
    
    		var fontface = "Helvetica";
    
    		var fontsize = parameters.hasOwnProperty("fontsize") ? 
    			parameters["fontsize"] : 18;
    		fontsize = fontsize*10;
    
    		var canvas = document.createElement('canvas');
    		var size = 1024;
    		canvas.width = size;
    		canvas.height = size;
    		var context = canvas.getContext('2d');
    		context.font = fontsize + "px " + fontface;
    
    		// text color
    		context.fillStyle = "rgba(0, 0, 0, 1.0)";
    
    		context.fillText(message, size/2, size/2);
    
    		// canvas contents will be used for a texture
    		var texture = new THREE.Texture(canvas);
    		texture.needsUpdate = true;
    
    		var spriteMaterial = new THREE.SpriteMaterial(
    			{map: texture, useScreenCoordinates: false});
    		var sprite = new THREE.Sprite(spriteMaterial);
    		return sprite;
    	}
    
    	function takeSvgScreenshot(){
    		if (labelsShown){
    			hideLabels();
    		}
    		svgRenderer.render(scene,camera);
    		svgElement = XMLS.serializeToString(svgRenderer.domElement);
    		
    		if (labelsShown){
    			displayLabels();
    		}
    
    		if (document.getElementById('tab_0').checked){
    			//show in new tab
    			var myWindow = window.open("","");
    			myWindow.document.body.innerHTML = svgElement;
    		} else{
    			// download svg file 
    			download("screenshot.svg", svgElement);
    		}
    	}
    		
    
    	function showOrHideObject(event){
    		var nr = Number(event.currentTarget.name);
    		all_objects[nr].visible = event.currentTarget.checked;
    	}
    
    	function displayOrHideOptionsRecursive( obj ) {
    		for (var j=0; j<obj.children.length; j++) {
    			var child = obj.children[j];
    			if (child.material===undefined && child) {
    				displayOrHideOptionsRecursive( child );
    			} else {
    				if (child.material.type == "MultiMaterial") {
    					for (var k=0; k<child.material.materials.length; k++) {
    						if (child.material.materials[k].transparent) {
    							document.getElementById('transparency_0').style.display = 'block';
    							document.getElementById('transparencyRange_0').value = 1 - 
    								child.material.materials[k].opacity;
    							return;
    						}
    					}
    				} else if (	child.material.transparent && 
    								child.material.type == "MeshBasicMaterial" &&
    								child.geometry.type == "Geometry"){
    					document.getElementById('transparency_0').style.display = 'block';
    					return;
    				}
    			}
    		}
    	}
    
    	function displayOrHideOptions() {
    		for (var i=0; i<all_objects.length; i++) {
    			var obj = all_objects[i];
    			displayOrHideOptionsRecursive( obj );
    		}
    	}
    
    	displayOrHideOptions()
    
    
    
    
    // ---------------------- EXPLOSION ------------------------------------------------
    // ---------------------------------------------------------------------------------
    
    	function explode(factor) {
    		var obj, c;
    		var c0 = centroids[0];
    		for (var i = 0; i<centroids.length; ++i) {
    			c = centroids[i];
    			obj = all_objects[all_objects.length - centroids.length + i];
    			obj.position.set(c.x*factor, c.y*factor, c.z*factor);
    		}	
    	}
    
    	function triggerExplode(event){
    		explodeScale = Number(event.currentTarget.value);
    		explode(explodeScale);
    	}
    
    	function setExplodingSpeed(event){
    		explodingSpeed = Number(event.currentTarget.value);
    	}
    
    	function triggerAutomaticExplode(event){
    		if (event.currentTarget.checked){
    			startExploding();
    		} else {
    			clearIntervals();
    		}	
    	}
    
    	function startExploding(){
    		intervals.push(setInterval(explodingInterval, 25));
    	}
    
    
    	function explodingInterval(){
    		explodeScale += explodingSpeed;
    		if (explodeScale <= 6){ 
    			explode(explodeScale);
    		}
    		else{
    			explode(6);
    			explodeScale = 6;
    			clearIntervals();
    			timeouts.push(setTimeout(startUnexploding, 3000));
    		}
    		document.getElementById('explodeRange_0').value = explodeScale;
    	}
    
    
    	function startUnexploding(){
    		intervals.push(setInterval(unexplodingInterval, 25));
    	}
    
    	function unexplodingInterval(){
    		explodeScale -= explodingSpeed;
    		if (explodeScale >= 0){	
    			explode(explodeScale);
    		}
    		else {
    			explode(0);
    			explodeScale = 0;
    			clearIntervals();
    			timeouts.push(setTimeout(startExploding, 3000));
    		}
    		document.getElementById('explodeRange_0').value = explodeScale;
    	}
    
    	function clearIntervals(){
    		intervals.forEach(function(interval){
    			clearInterval(interval);
    		});
    		intervals = [];
    		timeouts.forEach(function(timeout){
    			clearTimeout(timeout);
    		});
    		timeouts = [];
    	}
    
    			
    
    	// append checkboxes for displaying or hiding objects
    	var shownObjectsList = document.getElementById('shownObjectsList_0');
    	for (var i=0; i<all_objects.length; i++){
    		var objNode = document.createElement('span');
    		objNode.innerHTML = objectnames[i] + '<br>';
    		var checkbox = document.createElement('input');
    		checkbox.type = 'checkbox';
    		checkbox.checked = true;
    		checkbox.name = String(i);
    		checkbox.onchange = showOrHideObject;
    		shownObjectsList.appendChild(checkbox);
    		shownObjectsList.appendChild(objNode);
    	}
    
    	function displayLabels(){
    		for (var i=0; i<all_objects.length; i++){
    			for (var j=0; j<all_objects[i].children.length; j++){
    				var child = all_objects[i].children[j];
    				if (child.type == 'Sprite'){
    					child.visible = true;
    				}
    			}
    		}
    	}
    
    	function hideLabels(){
    		for (var i=0; i<all_objects.length; i++){
    			for (var j=0; j<all_objects[i].children.length; j++){
    				var child = all_objects[i].children[j];
    				if (child.type == 'Sprite'){
    					child.visible = false;
    				}
    			}
    		}
    	}
    
    	function displayOrHideLabels(event){
    		if (event.currentTarget.checked){
    			displayLabels();
    			labelsShown = true;
    		} else {
    			hideLabels();
    			labelsShown = false;
    		}
    	}
    
    	function download(filename, text) {
    	  var element = document.createElement('a');
    	  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    	  element.setAttribute('download', filename);
    
    	  element.style.display = 'none';
    	  document.body.appendChild(element);
    
    	  element.click();
    
    	  document.body.removeChild(element);
    	}
    
    var tempobj;
    tempobj = document.getElementById('explodeRange_0');
    if (tempobj) {
       tempobj.oninput = triggerExplode;
       document.getElementById('explodeCheckbox_0').onchange = triggerAutomaticExplode;
       document.getElementById('explodingSpeedRange_0').oninput = setExplodingSpeed;
    }
    tempobj = document.getElementById('foldRange_0');
    if (tempobj) {
       tempobj.oninput = fold;
    }
    document.getElementById('transparencyRange_0').oninput = changeTransparency;
    document.getElementById('changeRotationX_0').onchange = changeRotationX;
    document.getElementById('changeRotationY_0').onchange = changeRotationY;
    document.getElementById('changeRotationZ_0').onchange = changeRotationZ;
    document.getElementById('resetButton_0').onclick = resetScene;
    document.getElementById('rotationSpeedRange_0').oninput = changeRotationSpeedFactor;
    document.getElementById('labelsCheckboxInput_0').onchange = displayOrHideLabels;
    document.getElementById('takeScreenshot_0').onclick = takeSvgScreenshot;
    document.getElementById('showSettingsButton_0').onclick = showSettings;
    document.getElementById('hideSettingsButton_0').onclick = hideSettings;
    
    	
    
    // ------------------ SHORTCUTS --------------------------------------------
    // -------------------------------------------------------------------------
    
    /**
     * http://www.openjs.com/scripts/events/keyboard_shortcuts/
     * Version : 2.01.B
     * By Binny V A
     * License : BSD
     */
    shortcut = {
    	'all_shortcuts':{},//All the shortcuts are stored in this array
    	'add': function(shortcut_combination,callback,opt) {
    		//Provide a set of default options
    		var default_options = {
    			'type':'keydown',
    			'propagate':false,
    			'disable_in_input':false,
    			'target':document,
    			'keycode':false
    		}
    		if(!opt) opt = default_options;
    		else {
    			for(var dfo in default_options) {
    				if(typeof opt[dfo] == 'undefined') opt[dfo] = default_options[dfo];
    			}
    		}
    
    		var ele = opt.target;
    		if(typeof opt.target == 'string') ele = document.getElementById(opt.target);
    		var ths = this;
    		shortcut_combination = shortcut_combination.toLowerCase();
    
    		//The function to be called at keypress
    		var func = function(e) {
    			e = e || window.event;
    			
    			if(opt['disable_in_input']) { //Don't enable shortcut keys in Input, Textarea fields
    				var element;
    				if(e.target) element=e.target;
    				else if(e.srcElement) element=e.srcElement;
    				if(element.nodeType==3) element=element.parentNode;
    
    				if(element.tagName == 'INPUT' || element.tagName == 'TEXTAREA') return;
    			}
    	
    			//Find Which key is pressed
    			if (e.keyCode) code = e.keyCode;
    			else if (e.which) code = e.which;
    			var character = String.fromCharCode(code).toLowerCase();
    			
    			if(code == 188) character=","; //If the user presses , when the type is onkeydown
    			if(code == 190) character="."; //If the user presses , when the type is onkeydown
    
    			var keys = shortcut_combination.split("+");
    			//Key Pressed - counts the number of valid keypresses - if it is same as the number of keys, the shortcut function is invoked
    			var kp = 0;
    			
    			//Work around for stupid Shift key bug created by using lowercase - as a result the shift+num combination was broken
    			var shift_nums = {
    				"`":"~",
    				"1":"!",
    				"2":"@",
    				"3":"#",
    				"4":"$",
    				"5":"%",
    				"6":"^",
    				"7":"&",
    				"8":"*",
    				"9":"(",
    				"0":")",
    				"-":"_",
    				"=":"+",
    				";":":",
    				"'":"\"",
    				",":"<",
    				".":">",
    				"/":"?",
    				"\\":"|"
    			}
    			//Special Keys - and their codes
    			var special_keys = {
    				'esc':27,
    				'escape':27,
    				'tab':9,
    				'space':32,
    				'return':13,
    				'enter':13,
    				'backspace':8,
    	
    				'scrolllock':145,
    				'scroll_lock':145,
    				'scroll':145,
    				'capslock':20,
    				'caps_lock':20,
    				'caps':20,
    				'numlock':144,
    				'num_lock':144,
    				'num':144,
    				
    				'pause':19,
    				'break':19,
    				
    				'insert':45,
    				'home':36,
    				'delete':46,
    				'end':35,
    				
    				'pageup':33,
    				'page_up':33,
    				'pu':33,
    	
    				'pagedown':34,
    				'page_down':34,
    				'pd':34,
    	
    				'left':37,
    				'up':38,
    				'right':39,
    				'down':40,
    	
    				'f1':112,
    				'f2':113,
    				'f3':114,
    				'f4':115,
    				'f5':116,
    				'f6':117,
    				'f7':118,
    				'f8':119,
    				'f9':120,
    				'f10':121,
    				'f11':122,
    				'f12':123
    			}
    	
    			var modifiers = { 
    				shift: { wanted:false, pressed:false},
    				ctrl : { wanted:false, pressed:false},
    				alt  : { wanted:false, pressed:false},
    				meta : { wanted:false, pressed:false}	//Meta is Mac specific
    			};
                            
    			if(e.ctrlKey)	modifiers.ctrl.pressed = true;
    			if(e.shiftKey)	modifiers.shift.pressed = true;
    			if(e.altKey)	modifiers.alt.pressed = true;
    			if(e.metaKey)   modifiers.meta.pressed = true;
                            
    			for(var i=0; k=keys[i],i<keys.length; i++) {
    				//Modifiers
    				if(k == 'ctrl' || k == 'control') {
    					kp++;
    					modifiers.ctrl.wanted = true;
    
    				} else if(k == 'shift') {
    					kp++;
    					modifiers.shift.wanted = true;
    
    				} else if(k == 'alt') {
    					kp++;
    					modifiers.alt.wanted = true;
    				} else if(k == 'meta') {
    					kp++;
    					modifiers.meta.wanted = true;
    				} else if(k.length > 1) { //If it is a special key
    					if(special_keys[k] == code) kp++;
    					
    				} else if(opt['keycode']) {
    					if(opt['keycode'] == code) kp++;
    
    				} else { //The special keys did not match
    					if(character == k) kp++;
    					else {
    						if(shift_nums[character] && e.shiftKey) { //Stupid Shift key bug created by using lowercase
    							character = shift_nums[character]; 
    							if(character == k) kp++;
    						}
    					}
    				}
    			}
    			
    			if(kp == keys.length && 
    						modifiers.ctrl.pressed == modifiers.ctrl.wanted &&
    						modifiers.shift.pressed == modifiers.shift.wanted &&
    						modifiers.alt.pressed == modifiers.alt.wanted &&
    						modifiers.meta.pressed == modifiers.meta.wanted) {
    				callback(e);
    	
    				if(!opt['propagate']) { //Stop the event
    					//e.cancelBubble is supported by IE - this will kill the bubbling process.
    					e.cancelBubble = true;
    					e.returnValue = false;
    	
    					//e.stopPropagation works in Firefox.
    					if (e.stopPropagation) {
    						e.stopPropagation();
    						e.preventDefault();
    					}
    					return false;
    				}
    			}
    		}
    		this.all_shortcuts[shortcut_combination] = {
    			'callback':func, 
    			'target':ele, 
    			'event': opt['type']
    		};
    		//Attach the function with the event
    		if(ele.addEventListener) ele.addEventListener(opt['type'], func, false);
    		else if(ele.attachEvent) ele.attachEvent('on'+opt['type'], func);
    		else ele['on'+opt['type']] = func;
    	},
    
    	//Remove the shortcut - just specify the shortcut and I will remove the binding
    	'remove':function(shortcut_combination) {
    		shortcut_combination = shortcut_combination.toLowerCase();
    		var binding = this.all_shortcuts[shortcut_combination];
    		delete(this.all_shortcuts[shortcut_combination])
    		if(!binding) return;
    		var type = binding['event'];
    		var ele = binding['target'];
    		var callback = binding['callback'];
    
    		if(ele.detachEvent) ele.detachEvent('on'+type, callback);
    		else if(ele.removeEventListener) ele.removeEventListener(type, callback, false);
    		else ele['on'+type] = false;
    	}
    }
    
    shortcut.add("Alt+Left",function() {
    	var event = new Event('click');
    	if (settingsShown){
    		document.getElementById('hideSettingsButton_0').dispatchEvent(event);
    	} else{
    		document.getElementById('showSettingsButton_0').dispatchEvent(event);
    	}
    });
    
    if (foldable) moveToBaryCenter();
    
    
    });});
    // COMMON_CODE_BLOCK_END
    </script>
    
    </body>
    </html>



.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package threejs
       Three.js is a lightweight cross-browser JavaScript library/API used to create and display animated 3D computer graphics on a Web browser.
       See http://github.com/mrdoob for the source code.
    
    </pre>
    </details>




The output may look similar to the following.

.. figure:: attachment:cube-in-lattice.gif
   :alt: {{:tutorial:cube-in-lattice.gif?298}}

   {{:tutorial:cube-in-lattice.gif?298}}

The command ``LATTICE_COLORED`` sorted the lattice points into three
classes before visualization: lattice points in the interior of the
polytope, lattice points on the boundary, and vertices that are not in
the lattice. These classes are then visualized with different colors
(where we only see two in the above picture, as all vertices of the cube
are in the lattice). If you don’t need this distinction,
``VISUAL->LATTICE`` avoids the additional computations.

External Packages
-----------------

``polymake`` can use ``4ti2`` and ``lattE`` via a file based interface
and ``libnormaliz >= 3.1.0`` as library, (the file based interface to
``normaliz`` has been discontinued) for lattice computations and prints
all available packages during startup. To tell ``polymake`` about a
newly installed program run ``%%polymake --reconfigure%%`` or issue the
command ``reconfigure`` during the interactive session. polymake may ask
you to confirm the paths to the binaries.

::

   Application polytope uses following third-party software (for details: help 'credits';)
   4ti2, cddlib, latte, libnormaliz, lrslib, nauty

The output at this position depends on the software available on your
computer. To see each call to an external program you can either set the
variable ``$Verbose::external=1;`` on the command line or include the
line ``$Polymake::User::Verbose::external=1`` in
``~/.polymake/prefer.pl``. If you just want to see the credit message
instead of the program call, set ``$Verbose::credits=2`` instead. If
this is 1, then a credit is shown when a package is used for the first
time, if 0, then all credits are suppressed (but you can find them in
the file afterwards).


.. link

.. CODE-BLOCK:: perl

    polymake> $Verbose::external=1;

.. link

.. CODE-BLOCK:: perl

    polymake> print $p->EHRHART_POLYNOMIAL;
    8*x^3 + 12*x^2 + 6*x + 1




You can ask ``polymake`` to prefer one package over another by setting
``prefer "program";`` where program is one of ``_4ti2``, ``latte`` and
``normaliz2``. Of course, the corresponding package needs to be
installed on your computer.

To prefer one program only for some computations you may append one of
.integer_points, .hilbert, .ehrhartpoly for rules computing
N_LATTICE_POINTS, LATTICE_POINTS, HILBERT_BASIS or EHRHART_POLYNOMIAL.
(Or ``prefer_now`` just for the next computation)


.. link

.. CODE-BLOCK:: perl

    polymake> print cube(2)->N_LATTICE_POINTS;
    9




.. link

.. CODE-BLOCK:: perl

    polymake> prefer_now "libnormaliz";
    polymake> print cube(2)->N_LATTICE_POINTS;
    9




.. link

.. CODE-BLOCK:: perl

    polymake> print cube(2)->EHRHART_POLYNOMIAL;
    4*x^2 + 4*x + 1


