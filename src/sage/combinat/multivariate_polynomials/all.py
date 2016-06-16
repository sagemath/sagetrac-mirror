from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import
lazy_import('sage.combinat.multivariate_polynomials.bases',
            ['SchubertPolynomialsOnVectors', 'DemazurePolynomials'
             'DemazureHatPolynomials', 'GrothendieckPolynomials'])

lazy_import('sage.combinat.multivariate_polynomials.multivariate_polynomials',
            ['AbstractPolynomialRing'])  # FIXME: Remove this as an import

