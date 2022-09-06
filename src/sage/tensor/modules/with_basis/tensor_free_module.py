r"""
Parent of Components
"""

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class TensorFreeModuleWithBasis(UniqueRepresentation, Parent):

    @staticmethod
    def __classcall__(cls, base_ring, frames, start_indices, output_formatter=None, multiplicities=None):
        def list_to_tuple(x):
            if isinstance(x, list):
                return tuple(x)
            return x
        frames = tuple(list_to_tuple(frame) for frame in frames)
        start_indices = tuple(start_indices)
        return super().__classcall__(cls, base_ring, frames, start_indices, output_formatter)

    def __init__(self, base_ring, frames, start_indices, output_formatter):
        self._output_formatter = output_formatter
