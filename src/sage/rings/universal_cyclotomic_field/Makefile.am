# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
rings_universal_cyclotomic_fielddir = $(pythondir)/sage/rings/universal_cyclotomic_field
rings_universal_cyclotomic_field_PYTHON = __init__.py all.py universal_cyclotomic_field.py

AM_DEFAULT_SOURCE_EXT = .pyx
rings_universal_cyclotomic_field_LTLIBRARIES = universal_cyclotomic_field_c.la
universal_cyclotomic_field_c_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
universal_cyclotomic_field_c_la_CPPFLAGS = @NUMPY_INCLUDES@ $(AM_CPPFLAGS)
