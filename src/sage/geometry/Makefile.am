# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = triangulation
include $(top_srcdir)/common.mk
geometrydir = $(pythondir)/sage/geometry
geometry_PYTHON = \
    cone.py polyhedra.py __init__.py fan_isomorphism.py toric_plotter.py \
    fan_morphism.py fan.py all.py hasse_diagram.py lattice_polytope.py polytope.py \
    toric_lattice.py newton_polygon.py pseudolines.py


AM_DEFAULT_SOURCE_EXT = .pyx
geometry_LTLIBRARIES = integral_points.la point_collection.la toric_lattice_element.la
toric_lattice_element_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage

geometry_polyhedrondir = $(geometrydir)/polyhedron
geometry_polyhedron_PYTHON = \
    polyhedron/all.py \
    polyhedron/backend_cdd.py \
    polyhedron/backend_ppl.py \
    polyhedron/base.py \
    polyhedron/base_QQ.py \
    polyhedron/base_RDF.py \
    polyhedron/base_ZZ.py \
    polyhedron/cdd_file_format.py \
    polyhedron/constructor.py \
    polyhedron/face.py \
    polyhedron/__init__.py \
    polyhedron/library.py \
    polyhedron/misc.py \
    polyhedron/palp_database.py \
    polyhedron/parent.py \
    polyhedron/plot.py \
    polyhedron/ppl_lattice_polytope.py \
    polyhedron/representation.py \
    polyhedron/lattice_euclidean_group_element.py \
    polyhedron/ppl_lattice_polygon.py
