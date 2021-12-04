from __future__ import annotations
from typing import Generic, TypeVar
from collections.abc import Callable
from typing_extensions import Literal
from sage.manifolds.differentiable.chart import DiffChart, DiffCoordChange
from sage.manifolds.differentiable.manifold import DifferentiableManifold
from sage.manifolds.differentiable.vectorframe import VectorFrame
from sage.structure.element import Matrix

V = TypeVar("V")


class NaturalLocalData(Generic[V]):
    """
    Local data specified on a differentiable manifold that transforms naturally with respect to diffeomorphisms.
    """

    data: dict[VectorFrame, V]
    manifold: DifferentiableManifold
    action: Callable[[DiffCoordChange, V], V]

    def __init__(
        self,
        manifold: DifferentiableManifold,
        action: Callable[[DiffCoordChange, V], V],
    ):
        """
        Construct the local data.

        INPUT:

        - ``manifold`` -- the manifold on which the local data is defined
        - ``action`` -- a function that transforms the local data with respect to a diffeomorphism
        """
        self.manifold = manifold
        self.data = dict()
        self.action = action

        # TODO: We probably also need an action on frame changes

    def __getitem__(self, frame: VectorFrame) -> V:
        """
        Get the local data on a given frame.

        INPUT:

        - ``frame`` -- the frame on which the local data is defined

        OUTPUT:

        - the local data
        """
        return self.data[frame]

    def __setitem__(self, frame: VectorFrame, data: V):
        """
        Set the local data on a given frame.

        INPUT:

        - ``frame`` -- the frame on which the local data is defined
        - ``data`` -- the local data
        """
        if not frame.domain().is_subset(self.manifold):
            raise ValueError(f"{frame} must be defined on a subset of {self.manifold}")

        self.data[frame] = data

    def set(self, data: dict[VectorFrame, V]):
        """
        Set the local data from a dictionary.

        INPUT:

        - ``data`` -- a dictionary of local data
        """
        self.data = data

    def check_consistency(self) -> bool:
        """
        Check the consistency of the local data.
        """
        for coordChange in self.manifold.coord_changes().values():
            source_frame = coordChange.source().frame()
            target_frame = coordChange.target().frame()

            source_data = self.data[source_frame]
            target_data = self.data[target_frame]

            if (
                source_data
                and target_data
                and self.action(coordChange, source_data) != target_data
            ):
                return False

        # TODO: do the same for different frames

        return True

    def is_globally_defined(self) -> bool:
        """
        Check if the local data defines a global field.
        """
        domains_union = None
        for frame in self.data.keys():
            domain = frame.domain()
            if domains_union is not None:
                domains_union = domain.union(domains_union)
            else:
                domains_union = domain
        return domains_union == self.manifold


class NaturalLocalDataFirstOrder(NaturalLocalData[V]):
    def __init__(
        self, manifold: DifferentiableManifold, action: Callable[[Matrix, V], V]
    ):
        super().__init__(
            manifold, lambda coordChange, value: action(coordChange.jacobian(), value)
        )


OrientationValue = Literal[-1, +1]


class Orientation:
    """
    Orientation of a differentiable manifold.
    """

    data: NaturalLocalDataFirstOrder[OrientationValue]

    def __init__(self, manifold: DifferentiableManifold):
        self.data = NaturalLocalDataFirstOrder(
            manifold, lambda matrix, oriention: matrix.det().sign() * oriention
        )

    def __getitem__(self, frame: VectorFrame) -> OrientationValue:
        """
        Get the orientation of a given frame.

        INPUT:

        - ``frame`` -- the frame on which the orientation is defined

        OUTPUT:

        - the orientation
        """
        return self.data[frame]

    def __setitem__(self, target: VectorFrame | DiffChart, oriention: OrientationValue):
        """
        Set the orientation of a given frame.

        INPUT:

        - ``frame`` -- the frame on which the orientation is defined
        - ``oriention`` -- the orientation
        """
        if isinstance(target, VectorFrame):
            self.data[target] = oriention
        elif isinstance(target, DiffChart):
            self.data[target.frame()] = oriention
        else:
            raise TypeError(f"{target} is neither a vector frame nor chart")

    def set(self, data: dict[VectorFrame, OrientationValue]):
        self.data.set(data)

    def is_globally_defined(self) -> bool:
        return self.data.is_globally_defined()
