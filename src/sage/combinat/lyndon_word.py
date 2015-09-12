"""
Lyndon Word

This module has been moved to sage.combinat.words. It should be called from there.

TESTS::

    sage: from sage.combinat.lyndon_word import LyndonWord
"""

from sage.misc.lazy_import import lazy_import
lazy_import('sage.combinat.words.lyndon_word', '*', deprecation=19150)
