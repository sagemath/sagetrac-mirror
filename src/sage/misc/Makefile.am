# this file is part of Sage
# (c) 2013 Felix Salfelder
# license: gplv3+
#
include $(top_srcdir)/common.mk
miscdir = $(pythondir)/sage/misc
misc_PYTHON = \
    misc.py table.py map_threaded.py sage_unittest.py sage_extension.py \
    memory_info.py nested_class_test.py benchmark.py trace.py bindable_class.py \
    proof.py __init__.py python.py mathml.py latex_macros.py defaults.py \
    classgraph.py remote_file.py superseded.py lazy_string.py persist.py prun.py \
    sage_itertools.py unknown.py prandom.py flatten.py typecheck.py cache.py \
    functional.py edit_module.py exceptions.py package.py sage_input.py \
    method_decorator.py lazy_format.py decorators.py preparser.py \
    interpreter.py bug.py dist.py temporary_file.py \
    object_multiplexer.py inline_fortran.py multireplace.py messaging.py \
    file_to_worksheet.py sagedoc.py all.py sh.py explain_pickle.py \
    abstract_method.py random_testing.py cython.py sage_timeit.py log.py html.py \
    pager.py sage_eval.py displayhook.py banner.py getusage.py sageinspect.py \
    profiler.py lazy_import_cache.py mrange.py viewer.py db.py copying.py latex.py \
    func_persist.py dev_tools.py sage_ostools.py \
    test_class_pickling.py ascii_art.py attached_files.py

AM_DEFAULT_SOURCE_EXT = .pyx
misc_LTLIBRARIES = \
    allocator.la binary_tree.la bitset.la c3.la cachefunc.la citation.la \
    classcall_metaclass.la constant_function.la cython_c.la derivative.la \
    fast_methods.la fpickle.la function_mangling.la lazy_import.la lazy_list.la \
    misc_c.la nested_class.la parser.la pickle_old.la randstate.la \
    readline_extra_commands.la refcount.la reset.la sage_timeit_class.la search.la \
    session.la stopgap.la lazy_attribute.la c3_controlled.la inputhook.la

allocator_la_LIBADD = -lcsage
randstate_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage
#randstate_la_LDFLAGS = $(AM_LDFLAGS)

readline_extra_commands_la_LIBADD = -lreadline

bitset_la_LIBADD = -lcsage
binary_tree_la_LIBADD = -lcsage
misc_c_la_LIBADD = @GMP_OR_MPIR_LIBS@ -lcsage

noinst_HEADERS = bitset_pxd.pxi bitset.pxi
