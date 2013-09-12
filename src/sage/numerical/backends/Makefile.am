# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
inst_PYTHON = __init__.py

if USE_LIBCOIN
COND_COIN = coin_backend.la
endif

if USE_LIBCPLEX
COND_CPLEX = cplex_backend.la
endif

if USE_LIBGUROBI
COND_GUROBI = gurobi_backend.la
endif

AM_DEFAULT_SOURCE_EXT = .pyx
inst_LTLIBRARIES = \
	 $(COND_COIN) $(COND_CPLEX) $(COND_GUROBI) \
    generic_backend.la glpk_backend.la \
    glpk_graph_backend.la ppl_backend.la

coin_backend_la_LIBADD = -lcsage -lstdc++ -lCbc -lCbcSolver -lCgl -lClp -lCoinUtils -lOsiCbc -lOsiClp -lOsi
coin_backend_la_SOURCES = coin_backend.pyxx
coin_backend_la_CPPFLAGS = -I@LIBCSAGE_INCLUDES@ $(AM_CPPFLAGS)

cplex_backend_la_LIBADD = -lcsage -lstdc++ -lcplex
cplex_backend_la_CPPFLAGS = -I@LIBCSAGE_INCLUDES@ $(AM_CPPFLAGS)

generic_backend_la_LIBADD = -lcsage -lstdc++
generic_backend_la_CPPFLAGS = -I@LIBCSAGE_INCLUDES@ $(AM_CPPFLAGS)

glpk_backend_la_LIBADD = -lcsage -lstdc++ -lglpk @GMP_OR_MPIR_LIBS@ -lz
glpk_backend_la_SOURCES = glpk_backend.pyxx
glpk_backend_la_CPPFLAGS = -I@LIBCSAGE_INCLUDES@ $(AM_CPPFLAGS)

glpk_graph_backend_la_LIBADD = -lcsage -lstdc++ -lglpk @GMP_OR_MPIR_LIBS@ -lz
glpk_graph_backend_la_SOURCES = glpk_graph_backend.pyxx
glpk_graph_backend_la_CPPFLAGS = -I@LIBCSAGE_INCLUDES@ $(AM_CPPFLAGS)

gurobi_backend_la_LIBADD = -lcsage -lstdc++ -lcplex
gurobi_backend_la_CPPFLAGS = -I@LIBCSAGE_INCLUDES@ $(AM_CPPFLAGS)

ppl_backend_la_LIBADD = -lcsage -lstdc++
ppl_backend_la_CPPFLAGS = -I@LIBCSAGE_INCLUDES@ $(AM_CPPFLAGS)
