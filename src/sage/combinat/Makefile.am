# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
SUBDIRS = matrices words crystals
include $(top_srcdir)/common.mk
inst_PYTHON = \
    tutorial.py misc.py multichoose_nk.py family.py permutation_nk.py \
    partition_tuple.py sidon_sets.py yamanouchi.py lyndon_word.py __init__.py \
    rsk.py \
    binary_tree.py tableau.py q_analogues.py restricted_growth.py backtrack.py \
    set_partition.py kazhdan_lusztig.py combinatorial_algebra.py skew_tableau.py \
    perfect_matching.py integer_list.py symmetric_group_representations.py \
    combination.py schubert_polynomial.py output.py choose_nk.py tools.py \
    integer_vectors_mod_permgroup.py partition.py cyclic_sieving_phenomenon.py \
    parking_functions.py e_one_star.py ranker.py dlx.py set_partition_ordered.py \
    integer_matrices.py ribbon.py dyck_word.py tableau_tuple.py sloane_functions.py \
    tamari_lattices.py abstract_tree.py partition_algebra.py subword.py \
    combinatorial_map.py necklace.py core.py all.py knutson_tao_puzzles.py \
    alternating_sign_matrix.py subsets_pairwise.py integer_vector_weighted.py \
    tiling.py permutation.py graph_path.py derangements.py \
    gelfand_tsetlin_patterns.py symmetric_group_algebra.py integer_vector.py \
    yang_baxter_graph.py finite_class.py combinat.py tuple.py ribbon_tableau.py \
    cartesian_product.py subset.py free_module.py ordered_tree.py skew_partition.py \
    composition_signed.py split_nk.py composition.py \
    non_decreasing_parking_function.py affine_permutation.py vector_partition.py \
    diagram_algebras.py similarity_class_type.py

AM_DEFAULT_SOURCE_EXT = .pyx
inst_LTLIBRARIES = \
    combinat_cython.la debruijn_sequence.la degree_sequences.la dict_addition.la \
    enumeration_mod_permgroup.la expnums.la partitions.la permutation_cython.la \
    q_bernoulli.la

degree_sequences_la_LIBADD = -lcsage
combinat_cython_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
debruijn_sequence_la_LIBADD = -lcsage
expnums_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
partitions_la_SOURCES = partitions.pyxx partitions_c.cc
partitions_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lmpfr -lcsage

combinat_cluster_algebra_quiverdir = $(instdir)/cluster_algebra_quiver
combinat_cluster_algebra_quiver_PYTHON = \
    cluster_algebra_quiver/all.py \
    cluster_algebra_quiver/cluster_seed.py \
    cluster_algebra_quiver/__init__.py \
    cluster_algebra_quiver/mutation_class.py \
    cluster_algebra_quiver/mutation_type.py \
    cluster_algebra_quiver/quiver_mutation_type.py \
    cluster_algebra_quiver/quiver.py

combinat_posetsdir = $(instdir)/posets
combinat_posets_PYTHON = \
    posets/all.py \
    posets/elements.py \
    posets/hasse_diagram.py \
    posets/__init__.py \
    posets/lattices.py \
    posets/linear_extensions.py \
    posets/poset_examples.py \
    posets/posets.py

combinat_rigged_configurationsdir = $(instdir)/rigged_configurations
combinat_rigged_configurations_PYTHON = \
    rigged_configurations/all.py \
    rigged_configurations/bij_abstract_class.py \
    rigged_configurations/bijection.py \
    rigged_configurations/bij_type_A.py \
    rigged_configurations/bij_type_D.py \
    rigged_configurations/__init__.py \
    rigged_configurations/kleber_tree.py \
    rigged_configurations/kr_tableaux.py \
    rigged_configurations/rigged_configuration_element.py \
    rigged_configurations/rigged_configurations.py \
    rigged_configurations/rigged_partition.py \
    rigged_configurations/tensor_product_kr_tableaux_element.py \
    rigged_configurations/tensor_product_kr_tableaux.py

combinat_root_systemdir = $(instdir)/root_system
combinat_root_system_PYTHON = \
    root_system/all.py \
    root_system/ambient_space.py \
    root_system/associahedron.py \
    root_system/cartan_matrix.py \
    root_system/cartan_type.py \
    root_system/coxeter_group.py \
    root_system/coxeter_matrix.py \
    root_system/dynkin_diagram.py \
    root_system/__init__.py \
    root_system/pieri_factors.py \
    root_system/plot.py \
    root_system/root_lattice_realizations.py \
    root_system/root_space.py \
    root_system/root_system.py \
    root_system/type_A_affine.py \
    root_system/type_affine.py \
    root_system/type_A.py \
    root_system/type_B_affine.py \
    root_system/type_BC_affine.py \
    root_system/type_B.py \
    root_system/type_C_affine.py \
    root_system/type_C.py \
    root_system/type_D_affine.py \
    root_system/type_D.py \
    root_system/type_dual.py \
    root_system/type_E_affine.py \
    root_system/type_E.py \
    root_system/type_F_affine.py \
    root_system/type_F.py \
    root_system/type_G_affine.py \
    root_system/type_G.py \
    root_system/type_H.py \
    root_system/type_I.py \
    root_system/type_reducible.py \
    root_system/type_relabel.py \
    root_system/weight_lattice_realizations.py \
    root_system/weight_space.py \
    root_system/weyl_characters.py \
    root_system/weyl_group.py

combinat_sfdir = $(instdir)/sf
combinat_sf_PYTHON = \
    sf/all.py \
    sf/classical.py \
    sf/dual.py \
    sf/elementary.py \
    sf/hall_littlewood.py \
    sf/homogeneous.py \
    sf/__init__.py \
    sf/jack.py \
    sf/k_dual.py \
    sf/kfpoly.py \
    sf/kschur.py \
    sf/llt.py \
    sf/macdonald.py \
    sf/monomial.py \
    sf/multiplicative.py \
    sf/new_kschur.py \
    sf/ns_macdonald.py \
    sf/orthotriang.py \
    sf/powersum.py \
    sf/schur.py \
    sf/sfa.py \
    sf/sf.py

combinat_designsdir = $(instdir)/designs
combinat_designs_PYTHON = \
    designs/all.py \
    designs/block_design.py \
    designs/covering_design.py \
    designs/design_catalog.py \
    designs/ext_rep.py \
    designs/incidence_structures.py \
    designs/__init__.py \
    designs/steiner_quadruple_systems.py

combinat_ncsf_qsymdir = $(instdir)/ncsf_qsym
combinat_ncsf_qsym_PYTHON = \
    ncsf_qsym/all.py \
    ncsf_qsym/combinatorics.py \
    ncsf_qsym/generic_basis_code.py \
    ncsf_qsym/__init__.py \
    ncsf_qsym/ncsf.py \
    ncsf_qsym/qsym.py \
    ncsf_qsym/tutorial.py

combinat_speciesdir = $(instdir)/species
combinat_species_PYTHON = \
    species/all.py \
    species/characteristic_species.py \
    species/combinatorial_logarithm.py \
    species/composition_species.py \
    species/cycle_species.py \
    species/empty_species.py \
    species/functorial_composition_species.py \
    species/generating_series.py \
    species/__init__.py \
    species/library.py \
    species/linear_order_species.py \
    species/misc.py \
    species/partition_species.py \
    species/permutation_species.py \
    species/product_species.py \
    species/recursive_species.py \
    species/series_order.py \
    species/series.py \
    species/set_species.py \
    species/species.py \
    species/stream.py \
    species/structure.py \
    species/subset_species.py \
    species/sum_species.py
