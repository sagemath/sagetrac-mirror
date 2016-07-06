from __future__ import absolute_import

import sage.libs.ntl.all as ntl

from sage.libs.pari.all import pari, pari_gen, PariError

import sage.libs.symmetrica.all as symmetrica

from sage.misc.lazy_import import lazy_import
lazy_import('sage.libs.gap.libgap', 'libgap')

lazy_import('sage.libs.eclib.all', ('mwrank_EllipticCurve',
        'mwrank_MordellWeil', 'mwrank_initprimes', 'CremonaModularSymbols'))
lazy_import('sage.libs.eclib.all', 'get_precision', 'mwrank_get_precision')
lazy_import('sage.libs.eclib.all', 'set_precision', 'mwrank_set_precision')
