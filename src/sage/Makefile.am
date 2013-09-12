# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
ACLOCAL_AMFLAGS = -I m4

include common.mk

# FIXME: fallback to setup.py if desired
SUBDIRS = \
    libs algebras calculus categories coding combinat crypto databases \
    ext finance functions games geometry \
    graphs groups gsl homology \
    interfaces interacts lfunctions logic \
    matrix media misc modular modules monoids numerical \
    parallel plot probability quadratic_forms rings sandpiles \
    schemes sets stats structure symbolic tensor tests server doctest sat \
    dynamics matroids

inst_PYTHON = \
    all_cmdline.py \
    all_notebook.py \
    __init__.py \
    config.py \
    version.py \
    env.py \
    all.py

# runtime environment for bare sagelib installs
# (will be sourced by env.d)
00bare.sh:
	echo "export PYTHONPATH=\"$(pythondir)\$${PYTHONPATH/*/:}\$$PYTHONPATH\"" > $@
	echo "export SAGE_SRC=\"$(pythondir)\"" >> $@
	echo "export PATH=\"${prefix}/bin:\$$PATH\"" >> $@

sageenvdir=$(sysconfdir)/env.d
sageenv_DATA = 00bare.sh



# this implements the currently
# unavaible AM_EXTRA_RECURSIVE_TARGETS, see common.mk
py: @VPATH_TRUE@pys
	for i in $(SUBDIRS); do $(MAKE) -C $$i $@; done

pyall: .
	( cd $(VPATH); find . -name '*.h' -or -name '*.hh' -or -name '*.hpp' -or \
	     -name '*.py' -or -name '*.pxd' -or -name '*.pxi' -or \
	     -name '*.pyx' -or -name '*.pyxx' ) | \
	    sort > $@

pycheck: pyall
	echo $(PYS) $(noinst_HEADERS) | tr ' ' '\n' | \
	      sed 's#^#@abs_builddir@/#' | \
	      sed 's#^@abs_top_builddir@#.#' > pyreg
	for i in $(SUBDIRS); do $(MAKE) -C $$i $@ PYLIST=../pyreg; done
	sort pyreg > pyregs; mv pyregs pyreg
	diff -d pyreg pyall | tee pydiff; exit $${PIPESTATUS[0]}
	rm pyreg pyall pydiff

# in VPATH mode, we need __init__.py files, so cython -I works right
# (otherwise, $(INIT_PY_HERE) is empty)
pys: | $(INIT_PY_HERE)
check-local: | $(INIT_PY_HERE)

# move this up one level when merging core modules
sageconfdir=$(sysconfdir)
sageconf_DATA = default.qepcadrc
default.qepcadrc:
	echo SINGULAR $(dir $(SINGULAR)) > $@

# manually implementing AM_EXTRA_RECURSIVE_TARGETS([py])
# this will be in automake1.12 (?)
.PHONY: py $(PHONYPYS)

# bug or feature: know which files belong to sagelib
noinst_HEADERS = \
    modular/arithgroup/farey.pxd modular/modsym/apply.pxd \
    modular/modsym/p1list.pxd structure/mutability.pxd \
    structure/wrapper_parent.pxd structure/parent_gens.pxd \
    structure/coerce_maps.pxd structure/list_clone.pxd \
    structure/coerce_dict.pxd structure/generators.pxd \
    structure/sage_object.pxd structure/coerce_actions.pxd structure/misc.pxd \
    structure/coerce.pxd structure/parent_old.pxd structure/element.pxd \
    structure/category_object.pxd structure/parent_base.pxd \
    structure/parent.pxd numerical/backends/cplex_backend.pxd \
    numerical/backends/generic_backend.pxd numerical/backends/ppl_backend.pxd \
    numerical/backends/coin_backend.pxd numerical/backends/gurobi_backend.pxd \
    numerical/backends/glpk_backend.pxd \
    numerical/backends/glpk_graph_backend.pxd numerical/mip.pxd \
    numerical/linear_functions.pxd \
    graphs/modular_decomposition/modular_decomposition.pxd graphs/trees.pxd \
    graphs/generic_graph_pyx.pxd graphs/convexity_properties.pxd \
    graphs/cliquer.pxd graphs/graph_decompositions/rankwidth.pxd \
    graphs/distances_all_pairs.pxd graphs/base/dense_graph.pxd \
    graphs/base/c_graph.pxd graphs/base/static_sparse_graph.pxd \
    graphs/base/sparse_graph.pxd libs/ginac.pxd libs/gap/gap_includes.pxd \
    libs/gap/util.pxd libs/gap/element.pxd libs/lcalc/lcalc_Lfunction.pxd \
    libs/fplll/fplll.pxd libs/coxeter3/coxeter.pxd libs/gmp/random.pxd \
    libs/gmp/types.pxd libs/gmp/mpf.pxd libs/gmp/all.pxd libs/gmp/mpq.pxd \
    libs/gmp/mpn.pxd libs/gmp/mpz.pxd libs/pari/gen.pxd libs/m4rie.pxd \
    libs/mpmath/utils.pxd libs/mpmath/ext_impl.pxd libs/mpmath/ext_main.pxd \
    libs/m4ri.pxd libs/ntl/ntl_ZZ_p.pxd libs/ntl/ntl_ZZ_p_decl.pxd \
    libs/ntl/ntl_GF2EX.pxd libs/ntl/ntl_ZZ_pContext_decl.pxd \
    libs/ntl/ntl_GF2X.pxd libs/ntl/ntl_ZZ_decl.pxd \
    libs/ntl/ntl_lzz_pX_decl.pxd libs/ntl/ntl_ZZX.pxd libs/ntl/ntl_ZZ_pE.pxd \
    libs/ntl/ntl_GF2EContext.pxd libs/ntl/ntl_ZZ_pEContext_decl.pxd \
    libs/ntl/ntl_ZZ_pEContext.pxd libs/ntl/ntl_GF2E.pxd \
    libs/ntl/ntl_mat_GF2.pxd libs/ntl/ntl_GF2_decl.pxd \
    libs/ntl/ntl_ZZX_decl.pxd libs/ntl/ntl_GF2X_decl.pxd \
    libs/ntl/ntl_ZZ_pContext.pxd libs/ntl/ntl_ZZ_pE_decl.pxd \
    libs/ntl/ntl_vec_ZZ_p_decl.pxd libs/ntl/ntl_lzz_pContext.pxd \
    libs/ntl/ntl_lzz_pContext_decl.pxd libs/ntl/ntl_mat_ZZ.pxd \
    libs/ntl/ntl_vec_ZZ_pE_decl.pxd libs/ntl/ntl_ZZ_pX_decl.pxd \
    libs/ntl/ntl_vec_GF2_decl.pxd libs/ntl/ntl_ZZ_pEX_decl.pxd \
    libs/ntl/ntl_lzz_p.pxd libs/ntl/ntl_lzz_p_decl.pxd \
    libs/ntl/ntl_mat_GF2E.pxd libs/ntl/ntl_ZZ.pxd libs/ntl/ntl_GF2.pxd \
    libs/ntl/ntl_ZZ_pX.pxd libs/ntl/ntl_ZZ_pEX.pxd libs/ntl/ntl_lzz_pX.pxd \
    libs/linbox/linbox.pxd libs/linbox/echelonform.pxd \
    libs/linbox/modular.pxd libs/linbox/fflas.pxd libs/lrcalc/lrcalc.pxd \
    libs/mpfr.pxd libs/ecl.pxd libs/flint/flint.pxd \
    libs/flint/ntl_interface.pxd libs/flint/nmod_poly.pxd \
    libs/flint/fmpz_poly.pxd libs/flint/fmpq_poly.pxd \
    libs/flint/ulong_extras.pxd libs/ratpoints.pxd \
    libs/cremona/newforms.pxd libs/cremona/homspace.pxd libs/cremona/mat.pxd \
    libs/singular/singular.pxd libs/singular/ring.pxd \
    libs/singular/function.pxd libs/singular/decl.pxd \
    libs/singular/groebner_strategy.pxd libs/singular/polynomial.pxd \
    sets/finite_set_map_cy.pxd sets/disjoint_set.pxd symbolic/expression.pxd \
    symbolic/function.pxd symbolic/getitem.pxd symbolic/constants_c.pxd \
    matrix/matrix_rational_dense.pxd matrix/matrix_complex_double_dense.pxd \
    matrix/matrix2.pxd matrix/matrix_generic_dense.pxd matrix/matrix1.pxd \
    matrix/matrix_mod2e_dense.pxd matrix/matrix_modn_dense_double.pxd \
    matrix/matrix_integer_sparse.pxd matrix/matrix_domain_sparse.pxd \
    matrix/template.pxd matrix/matrix_symbolic_dense.pxd \
    matrix/matrix_integer_2x2.pxd matrix/matrix_sparse.pxd \
    matrix/matrix_generic_sparse.pxd matrix/matrix_mpolynomial_dense.pxd \
    matrix/matrix_real_double_dense.pxd matrix/matrix_modn_dense.pxd \
    matrix/matrix0.pxd matrix/matrix_dense.pxd matrix/matrix_cyclo_dense.pxd \
    matrix/matrix_window_modn_dense.pxd matrix/matrix_modn_dense_float.pxd \
    matrix/matrix_mod2_dense.pxd matrix/matrix_domain_dense.pxd \
    matrix/matrix_modn_sparse.pxd matrix/matrix_window.pxd \
    matrix/matrix_integer_dense.pxd matrix/matrix_rational_sparse.pxd \
    matrix/matrix.pxd matrix/action.pxd matrix/matrix_double_dense.pxd \
    modules/vector_real_double_dense.pxd modules/vector_modn_dense.pxd \
    modules/finite_submodule_iter.pxd modules/vector_integer_dense.pxd \
    modules/free_module_element.pxd modules/vector_complex_double_dense.pxd \
    modules/module.pxd modules/vector_rational_dense.pxd \
    modules/vector_double_dense.pxd modules/vector_mod2_dense.pxd \
    combinat/permutation_cython.pxd combinat/enumeration_mod_permgroup.pxd \
    combinat/combinat_cython.pxd combinat/debruijn_sequence.pxd \
    geometry/triangulation/functions.pxd geometry/triangulation/data.pxd \
    geometry/triangulation/triangulations.pxd groups/group.pxd groups/old.pxd \
    groups/perm_gps/partn_ref/refinement_python.pxd \
    groups/perm_gps/partn_ref/refinement_graphs.pxd \
    groups/perm_gps/partn_ref/double_coset.pxd \
    groups/perm_gps/partn_ref/automorphism_group_canonical_label.pxd \
    groups/perm_gps/partn_ref/canonical_augmentation.pxd \
    groups/perm_gps/partn_ref/refinement_binary.pxd \
    groups/perm_gps/partn_ref/refinement_sets.pxd \
    groups/perm_gps/partn_ref/refinement_matrices.pxd \
    groups/perm_gps/partn_ref/refinement_lists.pxd \
    groups/perm_gps/permgroup_element.pxd groups/libgap_wrapper.pxd \
    algebras/letterplace/free_algebra_element_letterplace.pxd \
    algebras/letterplace/free_algebra_letterplace.pxd \
    algebras/quatalg/quaternion_algebra_element.pxd plot/plot3d/shapes.pxd \
    plot/plot3d/transform.pxd plot/plot3d/parametric_surface.pxd \
    plot/plot3d/index_face_set.pxd plot/plot3d/base.pxd \
    finance/time_series.pxd coding/binary_code.pxd rings/real_mpfr.pxd \
    rings/finite_rings/element_base.pxd rings/finite_rings/integer_mod.pxd \
    rings/finite_rings/element_givaro.pxd \
    rings/finite_rings/element_ntl_gf2e.pxd \
    rings/finite_rings/finite_field_base.pxd rings/residue_field.pxd \
    rings/ring.pxd rings/real_lazy.pxd rings/power_series_poly.pxd \
    rings/morphism.pxd rings/fraction_field_FpT.pxd rings/fast_arith.pxd \
    rings/power_series_mpoly.pxd rings/real_double.pxd \
    rings/complex_double.pxd rings/rational.pxd rings/complex_interval.pxd \
    rings/integer_ring.pxd \
    rings/number_field/number_field_element_quadratic.pxd \
    rings/number_field/totallyreal_data.pxd \
    rings/number_field/number_field_element.pxd \
    rings/number_field/number_field_base.pxd \
    rings/power_series_ring_element.pxd rings/real_mpfi.pxd \
    rings/polynomial/multi_polynomial.pxd \
    rings/polynomial/polynomial_integer_dense_flint.pxd \
    rings/polynomial/polynomial_rational_flint.pxd \
    rings/polynomial/symmetric_reduction.pxd \
    rings/polynomial/polynomial_compiled.pxd rings/polynomial/pbori.pxd \
    rings/polynomial/plural.pxd rings/polynomial/polydict.pxd \
    rings/polynomial/polynomial_zmod_flint.pxd \
    rings/polynomial/polynomial_integer_dense_ntl.pxd \
    rings/polynomial/polynomial_gf2x.pxd \
    rings/polynomial/multi_polynomial_libsingular.pxd \
    rings/polynomial/polynomial_zz_pex.pxd \
    rings/polynomial/laurent_polynomial.pxd rings/polynomial/real_roots.pxd \
    rings/polynomial/polynomial_modn_dense_ntl.pxd \
    rings/polynomial/multi_polynomial_ring_generic.pxd \
    rings/polynomial/multi_polynomial_ideal_libsingular.pxd \
    rings/polynomial/polynomial_element.pxd \
    rings/laurent_series_ring_element.pxd rings/complex_number.pxd \
    rings/integer.pxd rings/complex_mpc.pxd \
    rings/padics/local_generic_element.pxd rings/padics/padic_printing.pxd \
    rings/padics/padic_capped_absolute_element.pxd \
    rings/padics/padic_fixed_mod_element.pxd rings/padics/padic_ext_element.pxd \
    rings/padics/padic_generic_element.pxd \
    rings/padics/padic_base_generic_element.pxd \
    rings/padics/padic_capped_relative_element.pxd \
    rings/padics/pow_computer.pxd rings/padics/padic_base_coercion.pxd \
    rings/padics/padic_ZZ_pX_element.pxd \
    rings/padics/padic_ZZ_pX_CR_element.pxd \
    rings/padics/padic_ZZ_pX_FM_element.pxd \
    rings/padics/padic_ZZ_pX_CA_element.pxd rings/padics/pow_computer_ext.pxd \
    gsl/dwt.pxd gsl/gsl_array.pxd gsl/fft.pxd gsl/ode.pxd \
    gsl/interpolation.pxd stats/hmm/distributions.pxd stats/hmm/util.pxd \
    stats/hmm/hmm.pxd stats/intlist.pxd misc/allocator.pxd \
    misc/nested_class.pxd misc/randstate.pxd misc/fast_methods.pxd \
    misc/lazy_list.pxd misc/misc_c.pxd misc/bitset.pxd misc/cachefunc.pxd \
    misc/classcall_metaclass.pxd misc/binary_tree.pxd ext/multi_modular.pxd \
    ext/fast_callable.pxd ext/fast_eval.pxd sat/solvers/satsolver.pxd \
    sat/solvers/cryptominisat/cryptominisat.pxd \
    sat/solvers/cryptominisat/solverconf.pxd sat/solvers/cryptominisat/decl.pxd \
    categories/morphism.pxd categories/map.pxd categories/functor.pxd \
    categories/category_singleton.pxd categories/action.pxd

noinst_HEADERS+= \
    tests/c_lib.h \
    graphs/modular_decomposition/src/dm_english.h \
    graphs/planarity_c/graphDrawPlanar.private.h \
    graphs/planarity_c/graphK23Search.h \
    graphs/planarity_c/graphFunctionTable.h \
    graphs/planarity_c/graphK33Search.private.h \
    graphs/planarity_c/listcoll.h \
    graphs/planarity_c/graphDrawPlanar.h \
    graphs/planarity_c/graphK4Search.h \
    graphs/planarity_c/graph.h \
    graphs/planarity_c/graphStructures.h \
    graphs/planarity_c/graphColorVertices.h \
    graphs/planarity_c/graphExtensions.h \
    graphs/planarity_c/planarity.h \
    graphs/planarity_c/platformTime.h \
    graphs/planarity_c/graphK33Search.h \
    graphs/planarity_c/graphK23Search.private.h \
    graphs/planarity_c/graphExtensions.private.h \
    graphs/planarity_c/appconst.h \
    graphs/planarity_c/graphK4Search.private.h \
    graphs/planarity_c/graphColorVertices.private.h \
    graphs/planarity_c/stack.h \
    graphs/graph_decompositions/rankwidth_c/rw.h \
    libs/lcalc/lcalc_sage.h \
    libs/pari/misc.h \
    libs/pari/pari_err.h \
    libs/mwrank/wrap.h \
    libs/eclsig.h \
    symbolic/pynac_cc.h \
    schemes/hyperelliptic_curves/hypellfrob/hypellfrob.h \
    schemes/hyperelliptic_curves/hypellfrob/recurrences_zn_poly.h \
    schemes/hyperelliptic_curves/hypellfrob/recurrences_ntl.h \
    combinat/matrices/dancing_links_c.h \
    combinat/partitions_c.h \
    geometry/triangulation/functions.h \
    geometry/triangulation/data.h \
    geometry/triangulation/triangulations.h \
    rings/bernmm/bern_modp_util.h \
    rings/bernmm/bern_rat.h \
    rings/bernmm/bern_modp.h \
    misc/darwin_memory_usage.h \
    ext/python_debug.h \
    ext/solaris_fixes.h \
    ext/mod_int.h \
    ext/mod_int.pxd \
    sat/solvers/cryptominisat/solverconf_helper.h \
    sat/solvers/cryptominisat/cryptominisat_helper.h

# provide a python module, so we can build docs
if BUILD_DOC
DOCPY=py
endif

CLEANFILES+=$(INIT_PY_HERE)

all: $(INIT_PY_HERE)
all-local: $(DOCPY)

# in case (parts of) build system needs to be regenerated
# take autotools from sage (the distribution) if necessary.
# (this is only required on toplevel)
ACLOCAL = export PATH="$(PATH)"; @ACLOCAL@
AUTOCONF = export PATH="$(PATH)"; @AUTOCONF@
AUTOHEADER = export PATH="$(PATH)"; @AUTOHEADER@
