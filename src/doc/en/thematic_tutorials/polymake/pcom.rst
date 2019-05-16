.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Polyhedral complexes in polymake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Polyhedral complexes are contained in the application ``fan``, so you
hanve to switch application to access the full functionality.


.. link

.. CODE-BLOCK:: perl

    polymake> application "fan";

To define polyhedral complexes in ``polymake``, you need to provide an
array of input points and a list of polytopes represented as an array of
arrays of point indices.


.. link

.. CODE-BLOCK:: perl

    polymake> $pc1 = new PolyhedralComplex(POINTS=>[[1,0,0],[1,0,1],[1,1,0],[1,1,1]],INPUT_POLYTOPES=>[[0,1,2],[2,3],[1]]);

Since some of the input polytopes may be redundant, you should ask for
the ``MAXIMAL_POLYTOPES``.


.. link

.. CODE-BLOCK:: perl

    polymake> print $pc1->MAXIMAL_POLYTOPES;
    {0 1 2}
    {2 3}





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




{{ :tutorial:pcom2.png }}

Triangulations
^^^^^^^^^^^^^^

Triangulations of polytopes form an important special class of polytopal
complexes. In polymake they are objects of type ``SimplicialComplex``
(and thus belong to the application ``topaz``). However, it is easy to
convert them as follows:


.. link

.. CODE-BLOCK:: perl

    polymake> $c=cube(3);
    polymake> $triangulation=new PolyhedralComplex(VERTICES=>$c->VERTICES,MAXIMAL_POLYTOPES=>$c->TRIANGULATION->FACETS);

Voronoi Diagrams and regular subdivisions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are seperate tutorials for `Voronoi diagrams <voronoi>`__ and
`regluar subdivisions <regular_subdivisions>`__ of point sets.
