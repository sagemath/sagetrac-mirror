r"""
Snake graph

This module implements snake graph and band graph matchings:
When a snake/band graph G comes from a generalized arc/loop sigma and
an initial ideal triangulation T of a bordered surface
(possibly with punctures and possibly with empty boundary),
the sum of the weights of all perfect matchings of G is equal to
the cluster algebra element corresponding to sigma
with respect to the cluster corresponding to T.
(see [MSW_Positivity]_ and [MSW_Bases]_).

REFERENCES:

.. [SchifflerAndTODO] Schiffler and TODO,
    *Snake Graph Algebra*

.. SEEALSO::

    Snake graphs closely interact with
    :class:`~sage.combinat.cluster_algebra_quiver.cluster_triangulation.ClusterTriangulation`

"""

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
#from sage.structure.list_clone import ClonableArray
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.combinatorial_map import combinatorial_map
from sage.structure.element import Element
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.combinat.cluster_algebra_quiver.cluster_snakegraph import ClusterSnakeGraph,\
ClusterSnakeGraphs

class ClusterWeightedSnakeGraph(ClusterSnakeGraph):
    """
    A snake graph.

    The weight of a perfect matching of a snake graph correspond to
    a numerator term
shape, first_tile_orientation=1,edge_weights={}, diagonal_weights={}
    INPUT:

    - ``shape`` -- a tuple/list listing the sizes of the rows of the snake graph

    - ``first_tile_orientation`` -- (default: 1) whether the first tile
    has orientation 1 or -1

    - ``edge_weights`` -- (default: {}) a dictionary of the weights of all edges

    - ``diagonal_weights`` -- (default: {}) a dictionary of the weights of all
    diagonals (of the tiles)

    EXAMPLES::

        sage: SG=ClusterWeightedSnakeGraph((2,1,3),first_tile_orientation=1)
        sage: SG
            -- -- --
           |  |  |  |
            -- -- --
           |  |
         -- --
        |  |  |
         -- --

        sage: #from sage.combinat.cluster_algebra_quiver.cluster_snakegraph_matching import ClusterSnakeGraphMatching
        sage: #ClusterSnakeGraphMatching(snake_graph=SG,tile_toggles=[])
        sage: #The minimal matching of a snake graph of shape (2,1,3)
        sage: #ClusterSnakeGraphMatching(snake_graph_shape=(2,1,3),tile_toggles=[])
        sage: #The minimal matching of a snake graph of shape (2,1,3)

        sage: #ClusterSnakeGraphMatching(snake_graph=SG,tile_toggles=(0,2,4,2))
        sage: #The matching of a snake graph of shape (2,1,3) with a flip sequence (0,4) from the minimal matching
        sage: #ClusterSnakeGraphMatching(snake_graph_shape=(2,1,3),tile_toggles=(0,2,4,2))
        sage: #The matching of a snake graph of shape (2,1,3) with a flip sequence (0,4) from the minimal matching

        sage: #ClusterSnakeGraphMatching([],[])
        sage: #Traceback (most recent call last):
        sage: #...
        sage: #ValueError: The first input must be a SnakeGraph class or a list (or tuple) of shape.
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    def __init__(self, shape, first_tile_orientation=1,edge_weights={}, diagonal_weights={}):
        """
        Initialize ``self``.

        TESTS::

            sage: G = ClusterWeightedSnakeGraph((2,1,1))
            sage: TestSuite(G).run()
        """
        self._shape = shape
        self._first_tile_orientation = 1
        self._edge_weights = edge_weights
        self._diagonal_weights = diagonal_weights
        ClusterSnakeGraph.__init__(self,ClusterSnakeGraphs(sum(shape)),shape)

    def plot(self, color='sign'):
        """
        Return a plot of ``self``.

        INPUT:

        - ``color`` -- can be any of the following:

          * ``4`` - use 4 colors: black, red, blue, and green with each
            corresponding to up, right, down, and left respectively
          * ``2`` - use 2 colors: red for horizontal, blue for vertical arrows
          * ``'sign'`` - use red for right and down arrows, blue for left
            and up arrows
          * a list of 4 colors for each direction
          * a function which takes a direction and a boolean corresponding
            to the sign

        EXAMPLES::

            sage: M = ClusterSnakeGraphs(3)
        """
        return 1

    def edge_weights(self):
        """
        Return the weights of ``self``, i.e. a dictionary with keys the position
        of the tiles of ``self`` and values (bottom,right,top,left), a tuple
        of initial cluster variables and boundary edge variables.

        EXAMPLES::

            sage: G = ClusterWeightedSnakeGraph([2,1,3])
            sage: G.edge_weights()
            {}
        """
        return self._edge_weights

    def diagonal_weights(self):
        """
        Return the weights of ``self``, i.e. a dictionary with keys the position
        of the tiles of ``self`` and values (bottom,right,top,left), a tuple
        of initial cluster variables and boundary edge variables.

        EXAMPLES::

            sage: G = ClusterWeightedSnakeGraph([2,1,3])
            sage: G.edge_weights()
            {}
        """
        return self._diagonal_weights

    def first_tile_orientation(self):
        """
        Return orientation of the first tile of the snake graph. The orientation
        is either 1 or -1.

        EXAMPLES::

            sage: G = ClusterWeightedSnakeGraph([2,1,3])
            sage: G.first_tile_orientation()
            1
        """
        return self._first_tile_orientation