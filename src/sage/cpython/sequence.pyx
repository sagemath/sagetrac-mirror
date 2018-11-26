"""
Utilities for "sequence-like" Python objects.

Here "sequence-like" is defined more closely to the Sequence ABC (abstract base
class) than the C-level ``PySequence`` API in that the object's type must
define ``__getitem__`` and ``__len__``, and should be possible to iterate over
monotonically from ``0`` to ``len(obj)`` like ``for idx in range(len(obj)):
obj[idx]``.
"""

# See src/sage/cpython/sequence.pxd for implementation details
