"""
Test relative lazy imports

TESTS::

    sage: from sage.misc import lazy_import_test
    sage: type(lazy_import_test.temporary_file)
    <type 'sage.misc.lazy_import.LazyImport'>
    sage: lazy_import_test.temporary_file.atomic_write
    <class sage.misc.temporary_file.atomic_write at ...>
    sage: type(lazy_import_test.version)
    <type 'sage.misc.lazy_import.LazyImport'>
    sage: lazy_import_test.version()
    'SageMath version ...'
"""

from __future__ import absolute_import

from .lazy_import import lazyimport
with lazyimport:
    from . import temporary_file
    from .all import version
