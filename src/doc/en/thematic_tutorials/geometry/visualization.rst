.. -*- coding: utf-8 -*-

.. linkall

.. _polyhedron_visualization:

==================================================
Visualization of polyhedron objects in Sage
==================================================

.. MODULEAUTHOR:: sarah-marie belcastro <smbelcas@toroidalsnark.net>, Jean-Philippe Labbé <labbe@math.fu-berlin.de>


There are different ways to visualize polyhedron object of dimension at most 4.

:code:`render_solid`
==================================================

This plots the polyhedron as a solid. You can also adjust the :code:`opacity`
parameter.

::

    sage: Cube = polytopes.cube()
    sage: Cube.render_solid(opacity=0.7)
    Graphics3d Object

.. end of output

:code:`render_wireframe`
==================================================

This plots the graph (with unbounded edges) of the polyhedron

::

    sage: Cube.render_wireframe()
    Graphics3d Object

.. end of output

:code:`plot`
==================================================

The :code:`plot` method draws the graph, the polygons and vertices of the
polyhedron all together.

::

    sage: Cube.plot()
    Graphics3d Object

.. end of output

:code:`show`
==================================================

This is similar to :code:`plot` but does not return an object that you can
manipulate.


:code:`schlegel_projection`
==================================================

It is possible to visualize 4-dimensional polytopes using a schlegel diagram.

::

    sage: HC = polytopes.hypercube(4)
    sage: HC.schlegel_projection()
    The projection of a polyhedron into 3 dimensions
    sage: HC.schlegel_projection().plot()
    Graphics3d Object

.. end of output

We can see it from a different perspective:

::

    sage: HC.schlegel_projection([2,5,11,17]).plot()
    Graphics3d Object

.. end of output
