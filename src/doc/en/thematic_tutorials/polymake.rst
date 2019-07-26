.. -*- coding: utf-8 -*-

.. _polymake:

=============================
 Thematic Tutorial: Polymake
=============================

polymake is a mature open source software package for research in
polyhedral geometry, by Ewgenij Gawrilow and Michael Joswig and
contributing authors. It deals with polytopes, polyhedra and fans as
well as simplicial complexes, matroids, graphs, tropical
hypersurfaces, and other objects.

This thematic tutorial was obtained by (1) converting from the polymake
tutorials available in form of Jupyter notebooks in polymake release
V3.4, (2) manually translating the Perl code to Python.

.. Converted from the Polymake tutorials in polymake/demo/*.ipynb from
   polymake release V3.4 using (in `sage -sh`)

   for a in `ls *.ipynb | grep -v perl_intro | grep -v apps_polytope`; do
     sage-ipynb2rst-polymake $a $SAGE_ROOT/src/doc/en/thematic_tutorials/polymake/`basename $a .ipynb`.rst;
   done

   The following toctree must be updated manually.

.. toctree::
   :maxdepth: 1

   polymake/apps_fan
   polymake/apps_fulton
   polymake/apps_graph
   polymake/apps_group
   polymake/apps_matroid
   polymake/apps_polytope
   polymake/apps_topaz
   polymake/apps_tropical
   polymake/aut_of_graphs
   polymake/caratheodory
   polymake/chain_complex_homology
   polymake/coordinates
   polymake/data
   polymake/face_lattice_tutorial
   polymake/hyperbolic_surface_tutorial
   polymake/ilp_and_hilbertbases
   polymake/intro_tutorial
   polymake/lattice_polytopes_tutorial
   polymake/legacy
   polymake/matching_polytopes
   polymake/matrix_classes
   polymake/optimization
   polymake/pcom
   polymake/persistent_homology
   polymake/polynomials_tutorial
   polymake/polytope_semantics
   polymake/properties
   polymake/random
   polymake/regular_subdivisions
   polymake/tarballs
   polymake/transformations
   polymake/visual_tutorial
   polymake/voronoi


.. Removed:
.. polymake/perl_intro


