"""
Crystal features that are imported by default in the interpreter namespace
"""
from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.crystals', 'catalog', 'crystals')
