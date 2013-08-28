# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
rings_number_fielddir = $(pythondir)/sage/rings/number_field
rings_number_field_PYTHON = \
    all.py class_group.py galois_group.py __init__.py maps.py morphism.py \
    number_field_ideal.py number_field_ideal_rel.py number_field.py \
    number_field_rel.py order.py small_primes_of_degree_one.py totallyreal_phc.py \
    totallyreal_rel.py unit_group.py

AM_DEFAULT_SOURCE_EXT = .pyx
rings_number_field_LTLIBRARIES = \
    number_field_base.la number_field_element.la \
	 number_field_element_quadratic.la \
	 number_field_morphisms.la totallyreal.la totallyreal_data.la

number_field_element_la_SOURCES = number_field_element.pyxx
number_field_element_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage

number_field_element_quadratic_la_SOURCES = number_field_element_quadratic.pyxx
number_field_element_quadratic_la_LIBADD = -lntl @GMP_OR_MPIR_LIBS@ -lcsage

totallyreal_la_LIBADD = -lpari @GMP_OR_MPIR_LIBS@ -lcsage
totallyreal_data_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
