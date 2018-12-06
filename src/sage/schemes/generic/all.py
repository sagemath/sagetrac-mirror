from __future__ import absolute_import
# code exports
from sage.misc.lazy_import import lazy_import

lazy_import('sage.schemes.generic.spec', 'Spec')
lazy_import('sage.schemes.generic.hypersurface',
            ['ProjectiveHypersurface', 'AffineHypersurface'])
