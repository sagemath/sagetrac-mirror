from __future__ import absolute_import
# code exports

from sage.misc.lazy_import import lazy_import

from .fano_variety import CPRFanoToricVariety
lazy_import('sage.schemes.ideal','ToricIdeal')
lazy_import('sage.schemes.library', 'toric_varieties')
lazy_import('sage.schemes.toric.variety',
            ['AffineToricVariety', 'ToricVariety'])
lazy_import('sage.schemes.toric.weierstrass', 'WeierstrassForm')
