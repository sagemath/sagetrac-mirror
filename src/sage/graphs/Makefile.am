# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = base graph_decompositions modular_decomposition
include $(top_srcdir)/common.mk
inst_PYTHON = \
    pq_trees.py graph_database.py dot2tex_utils.py generic_graph.py __init__.py \
    graph_coloring.py graph_list.py graph_bundle.py hypergraph_generators.py \
    linearextensions.py schnyder.py graph.py digraph_generators.py \
    graph_generators.py graph_latex.py all.py digraph.py print_graphs.py \
    bipartite_graph.py graph_editor.py isgci.py graph_plot.py line_graph.py \
    hypergraph.py

AM_DEFAULT_SOURCE_EXT = .pyx
inst_LTLIBRARIES = \
    chrompoly.la cliquer.la comparability.la convexity_properties.la \
    distances_all_pairs.la generic_graph_pyx.la genus.la graph_generators_pyx.la \
    hyperbolicity.la matchpoly.la planarity.la spanning_tree.la trees.la \
    weakly_chordal.la

chrompoly_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
cliquer_la_SOURCES = cliquer.pyx cl.c
cliquer_la_LIBADD = -lcliquer -lcsage
comparability_la_LIBADD = -lcsage
distances_all_pairs_la_LIBADD = -lcsage
generic_graph_pyx_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
matchpoly_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lflint -lcsage
matchpoly_la_CFLAGS = -std=c99 $(AM_CFLAGS)
planarity_la_SOURCES = \
    planarity.pyx planarity_c/graphColorVertices.c \
    planarity_c/graphColorVertices_Extensions.c planarity_c/graphDrawPlanar.c \
    planarity_c/graphDrawPlanar_Extensions.c planarity_c/graphEmbed.c \
    planarity_c/graphExtensions.c planarity_c/graphIO.c planarity_c/graphIsolator.c \
    planarity_c/graphK23Search.c planarity_c/graphK23Search_Extensions.c \
    planarity_c/graphK33Search.c planarity_c/graphK33Search_Extensions.c \
    planarity_c/graphK4Search.c planarity_c/graphK4Search_Extensions.c \
    planarity_c/graphNonplanar.c planarity_c/graphOuterplanarObstruction.c \
    planarity_c/graphPreprocess.c planarity_c/graphTests.c planarity_c/graphUtils.c \
    planarity_c/listcoll.c planarity_c/planarity_c.c \
    planarity_c/planarityCommandLine.c planarity_c/planarityRandomGraphs.c \
    planarity_c/planaritySpecificGraph.c planarity_c/planarityUtils.c \
    planarity_c/stack.c

graphs_generatorsdir = $(instdir)/generators
graphs_generators_PYTHON = \
    generators/basic.py \
    generators/chessboard.py \
    generators/degree_sequence.py \
    generators/families.py \
    generators/__init__.py \
    generators/platonic_solids.py \
    generators/random.py \
    generators/smallgraphs.py \
    generators/world_map.py

trees_la_LIBADD = -lcsage
weakly_chordal_la_LIBADD = -lcsage

noinst_HEADERS = cl.h
EXTRA_DIST = cl.c
