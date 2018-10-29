from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import as _lazy_import

_lazy_import("sage.coding.code_constructions", ["permutation_action",
            "walsh_matrix"])

_lazy_import("sage.coding.linear_code", [
            "LinearCode",
            "LinearCodeFromVectorSpace",
            "self_orthogonal_binary_codes"])

# Functions removed from the global namespace

_lazy_import('sage.coding', 'codes_catalog', 'codes')
_lazy_import('sage.coding', 'channels_catalog', 'channels')
