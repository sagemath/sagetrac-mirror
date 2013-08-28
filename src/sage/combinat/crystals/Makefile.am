# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
combinat_crystalsdir = $(pythondir)/sage/combinat/crystals

combinat_crystals_PYTHON = \
    monomial_crystals.py \
    elementary_crystals.py \
    affine.py \
    alcove_path.py \
    all.py \
    crystals.py \
    direct_sum.py \
    fast_crystals.py \
    generalized_young_walls.py \
    highest_weight_crystals.py \
    infinity_crystals.py \
    __init__.py \
    kirillov_reshetikhin.py \
    littelmann_path.py \
    spins.py \
    tensor_product.py \
    kyoto_path_model.py

AM_DEFAULT_SOURCE_EXT = .pyx
combinat_crystals_LTLIBRARIES = letters.la
letters_la_LIBADD = -lcsage
