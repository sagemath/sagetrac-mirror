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


class ClusterSnakeGraphMatching(Element):
    """
    A perfect matching of a SnakeGraph.

    The weight of a perfect matching of a snake graph correspond to
    a numerator term

    INPUT:

    - ``snake_graph`` -- (default:``None``) a SnakeGraph class

    - ``snake_graph_shape`` -- (default: []) a list (or tuple) which
    describes the shape of the snake graph

    - ``tile_toggles`` -- (default: []) a list (or tuple) of tiles which is flipped
    to get to self from the minimal matching

    .. WARNING::

        either ``snake_graph`` or ``snake_graph_shape`` have to TODO

    EXAMPLES::

        sage: SG = ClusterSnakeGraph((2,1,3),first_tile_orientation=1)
        sage: SnakeGraphMatching(snake_graph=SG,tile_toggles='Min PM')
        The minimal matching of a snake graph of shape (2,1,3)
        sage: SnakeGraphMatching(snake_graph=SG,tile_toggles=[])
        The minimal matching of a snake graph of shape (2,1,3)
        sage: SnakeGraphMatching(snake_graph_shape=(2,1,3),tile_toggles='Min PM')
        The minimal matching of a snake graph of shape (2,1,3)

        sage: SnakeGraphMatching(snake_graph_shape=(2,1,3),tile_toggles='Min PM')
            -- -- --
           |  |  |  |
            -- -- --
           |  |
         -- --
        |  |  |
         -- --

               --
           |

              |
         --
              |
         --

        sage: SnakeGraphMatching(snake_graph=SG,tile_toggles=(0,2,4,2))
        The matching of a snake graph of shape (2,1,3) with a flip sequence (0,4) from the minimal matching
        sage: SnakeGraphMatching(snake_graph_shape=(2,1,3),tile_toggles=(0,2,4,2))
        The matching of a snake graph of shape (2,1,3) with a flip sequence (0,4) from the minimal matching

        sage: SnakeGraphMatching([],[])
        Traceback (most recent call last):
        ...
        ValueError: The first input must be a SnakeGraph class or a list (or tuple) of shape.
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, snakegraph):
        """
        Create an snake graph.

        EXAMPLES::

            sage: ClusterSnakeGraph([2,1,3])
             ***
             *
            **
        """
        from sage.combinat.composition import Compositions
        from sympy import Sum
        shape = snakegraph._shape
        if not shape in Compositions():
            raise ValueError("The input must be a composition of positive integers")
        SGs = ClusterSnakeGraphMatchings(sum(shape))
        return SGs(shape)

    def __init__(self, parent, matching={}, tile_toggles=[]):
        """
        Initialize ``self``.

        TESTS::

            sage: A = AlternatingSignMatrices(3)
            sage: elt = A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
            sage: TestSuite(elt).run()
        """
        if matching == [] and tile_toggles ==[]:
            raise ValueError("The first input must be a SnakeGraph class or a list (or tuple) of shape.")
        self._snakegraph = snakegraph
        self._matching = matching
        self._tile_toggles=tile_toggles
        #self._minimal_matching = [] # todo _get_minimal_matching()
        self._coefficients = {} # todo
        Element.__init__(self, parent)

    def check(self):
        """
        Check if ``self`` is a valid snake graph.

        EXAMPLES::

            sage: M = ClusterSnakeGraphs(3)
            sage: M[0].check()
        """
        if self not in self.parent():
            raise ValueError("invalid snake graph")

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ClusterSnakeGraphs(5).list()
            [
             --
            |  |
             --    -- --       --                   --
            |  |  |  |  |     |  |                 |  |
             --    -- --    -- --    -- -- --       --       -- --          --
            |  |  |  |     |  |  |  |  |  |  |     |  |     |  |  |        |  |
             --    --       -- --    -- -- --    -- --    -- -- --    -- -- --
            |  |  |  |     |  |     |  |        |  |  |  |  |  |     |  |  |  |
             --    --       --       --          -- --    -- --       -- -- --
            |  |  |  |     |  |     |  |        |  |     |  |        |  |
             -- ,  --    ,  --    ,  --       ,  --    ,  --       ,  --       ,


                               --
                              |  |
                               --       -- --          --
                              |  |     |  |  |        |  |
             -- -- -- --       --       -- --       -- --       -- -- --
            |  |  |  |  |     |  |     |  |        |  |  |     |  |  |  |
             -- -- -- --    -- --    -- --       -- -- --    -- -- -- --
            |  |           |  |  |  |  |  |     |  |  |     |  |  |
             --          ,  -- -- ,  -- --    ,  -- --    ,  -- --       ,


                   --
                  |  |
                   --          -- --             --
                  |  |        |  |  |           |  |
             -- -- --    -- -- -- --    -- -- -- --    -- -- -- -- --
            |  |  |  |  |  |  |  |     |  |  |  |  |  |  |  |  |  |  |
             -- -- -- ,  -- -- --    ,  -- -- -- -- ,  -- -- -- -- --
            ]
        """
        sh = self._shape
        ret = ''
        top_row = sh[-1]
        skips = sum(sh[:-1])-(len(sh)-1)
        white_sp = '   '

        ret += white_sp* skips
        ret +=' -- '
        for i in range(1,top_row):
            ret +='-- '
        ret +='\n' + white_sp * skips + '|  |'
        for i in range(1,top_row):
            ret +='  |'

        for i in range(len(sh)-2,-1,-1):
            r = sh[i]
            skips += -(r-1)

            ret +='\n' + white_sp * skips
            for i in range(0,r+sh[i+1]-1):
                ret +=' --'

            ret +='\n' + white_sp * skips + '|  |'
            for i in range(1,r):
                ret +='  |'

        bottom_row = sh[0]
        ret +='\n' + ' --'
        for i in range(1,bottom_row):
            ret +=' --'

        return ret

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

    def shape(self):
        """
        Return shape of the snake graph.
        """
        return self._shape

class ClusterSnakeGraphMatchings(Parent, UniqueRepresentation):
    """
    INPUT:

    - ``d`` -- the number of tiles

      There are also the following predefined boundary conditions:

    EXAMPLES:

    Here are the two types of snake graphs that can be created::

        sage: M = ClusterSnakeGraphs(2)
        sage: list(M)
        [
                  --
                 |  |
          -- --   --
         |  |  | |  |
          -- -- , --
        ]

    REFERENCES:

    """
    def __init__(self, d):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: M = ClusterSnakeGraphs(2)
            sage: TestSuite(M).run()
        """
        self._d = d
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ClusterSnakeGraphs(2, boundary_conditions='ice')
            The snake graphs with 2 tiles
        """
        return "The snake graphs with {} tiles".format(self._d)

    def _repr_option(self, key):
        """
        Metadata about the ``_repr_()`` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: M = ClusterSnakeGraph(2)
            sage: M._repr_option('element_ascii_art')
            True
        """
        if key == 'element_ascii_art':
            return True
        return Parent._repr_option(self, key)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: M = ClusterSnakeGraphs(6)
            sage: M((2,1,3))
                -- -- --
               |  |  |  |
                -- -- --
               |  |
             -- --
            |  |  |
             -- --
        """
        if isinstance(x, ClusterSnakeGraphMatching):
            if x in self.parent():
                return x
            else:
                print 'print an error message here TODO'
        else:
            return self.element_class(self, x)

    Element = ClusterSnakeGraphMatching

    def __iter__(self):
        """
        Iterate through ``self``.

        EXAMPLES::

            sage: M = ClusterSnakeGraphs(2)
            sage: list(M)
            [
                              --
                             |  |
                      -- --   --
                     |  |  | |  |
                      -- -- , --
            ]
        """
        from sage.combinat.composition import Compositions
        for c in Compositions(self._d):
            yield self.element_class(self,c)

    def number_of_tiles(self):
        """
        Return the boundary conditions of ``self``.

        EXAMPLES::

            sage: M = ClusterSnakeGraphs(4)
            sage: M.number_of_tiles()
            4
        """
        return self._d

