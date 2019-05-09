.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Random Constructions
====================

Random points on the unit sphere
--------------------------------

The easiest way to randomly construct a polytope is by sampling points
on the unit sphere. The following chooses 100 points on the units sphere
in 3-space.


::

    polymake> $p1=rand_sphere(3,100);
    ........> print $p1->SIMPLICIAL;
    true




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




With probability one such polytopes are simplicial.

Random polytopes with are neither simplicial nor simple
-------------------------------------------------------


::

    polymake> ($d,$m,$n) = (4,50,30);
    ........> $p1=rand_sphere($d,$m);
    ........> $p2=polarize($p1);
    ........> $p3=new Polytope(POINTS=>rand_vert($p2->VERTICES,$n));
    ........> print $p3->SIMPLICIAL, " ", $p3->SIMPLE, "\n", $p3->F_VECTOR;
    false false
    30 163 256 123


