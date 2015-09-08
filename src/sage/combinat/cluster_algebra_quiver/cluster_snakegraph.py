r"""
Snake Graphs

REFERENCES:

.. [CanakciSchiffler] Canakci and Schiffler,
    *Snake graph calculus and cluster algebras from surfaces*
    :arxiv:`abs/1209.4617`

.. [MSW_Positivity] Musiker - Schiffler - Williams,
   *Positivity for Cluster Algebras from Surfaces*,
   :arxiv:`0906.0748`

.. [MSW_Bases] Musiker - Schiffler - Williams,
   *Bases for Cluster Algebras from Surfaces*,
   :arxiv:`1110.4364`

.. [MW_MatrixFormulae] Musiker and Williams,
   *Matrix Formulae and Skein Relations for Cluster Algebras
   from Surfaces*,
   :arxiv:`1108.3382`

.. [FominShapiroThurston] Fomin - Shapiro - Thurston,
   *Cluster algebras and triangulated surfaces. part I: Cluster
   complexes*,
   :arxiv:`math/0608367`

.. SEEALSO::

    Cluster triangulations closely interact with
    :class:`~sage.combinat.cluster_algebra_quiver.cluster_seed.ClusterSeed`,
    :class:`~sage.combinat.cluster_algebra_quiver.quiver.ClusterQuiver`

"""

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.element import Element
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.list_clone import ClonableArray

class SnakeGraph(ClonableArray):
    """
    A snake graph is a connected sequence of square tiles.
    To build a snake graph, start with one tile, then glue a new tile so that
    the new tile is glued to the north or the east of the previous tile.
    See [MSW_Positivity] or [CanakciSchiffler]_.

    Note that the edges of the graph are not labeled. Hence a snake graph is uniquely
    determined by a list of positive integers (``shape``) such that their sum is
    equal to the number of the snake graph's tiles, i.e. a snake graph is uniquely
    determined by a composition of ``d``, where ``d`` is the number of tiles
    of the snake graph.

    INPUT:

    - ``shape`` -- a tuple/list listing the sizes of the rows of the snake graph

    EXAMPLES::

        sage: SnakeGraph((2,1,3))
            -- -- --
           |  |  |  |
            -- -- --
           |  |
         -- --
        |  |  |
         -- --
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, shape):
        """
        Create an snake graph.

        EXAMPLES::

            sage: SnakeGraph([2,1,3])
                -- -- --
               |  |  |  |
                -- -- --
               |  |
             -- --
            |  |  |
             -- --
        """
        from sage.combinat.composition import Compositions
        from sympy import Sum
        if not list(shape) in Compositions():
            raise ValueError("The input must be a composition of positive integers")
        SGs = SnakeGraphs(sum(shape))
        return SGs(shape)#, first_tile_orientation, edge_weights, diagonal_weights)

    def __init__(self, parent, shape):
        """
        Initialize ``self``.

        TESTS::

            sage: G = SnakeGraph((2,1,1))
            sage: TestSuite(G).run()
        """
        self._shape = list(shape)
        ClonableArray.__init__(self, parent, shape)

    def check(self):
        """
        Check if ``self`` is a valid snake graph.

        EXAMPLES::

            sage: M = SnakeGraphs(3)
            sage: M[0].check()
        """
        if self not in self.parent():
            raise ValueError("invalid snake graph")

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SnakeGraphs(5).list()[0:6]
            [
             --
            |  |
             --    -- --       --                   --
            |  |  |  |  |     |  |                 |  |
             --    -- --    -- --    -- -- --       --       -- --
            |  |  |  |     |  |  |  |  |  |  |     |  |     |  |  |
             --    --       -- --    -- -- --    -- --    -- -- --
            |  |  |  |     |  |     |  |        |  |  |  |  |  |
             --    --       --       --          -- --    -- --
            |  |  |  |     |  |     |  |        |  |     |  |
             -- ,  --    ,  --    ,  --       ,  --    ,  --
             ]

            sage: SnakeGraphs(5).list()[6:11]
            [
                                           --
                                          |  |
                   --                      --       -- --          --
                  |  |                    |  |     |  |  |        |  |
             -- -- --    -- -- -- --       --       -- --       -- --
            |  |  |  |  |  |  |  |  |     |  |     |  |        |  |  |
             -- -- --    -- -- -- --    -- --    -- --       -- -- --
            |  |        |  |           |  |  |  |  |  |     |  |  |
             --       ,  --          ,  -- -- ,  -- --    ,  -- --
            ]

            sage: SnakeGraphs(5).list()[12:16]
            [
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

    def shape(self):
        """
        Return the shape of ``self``. The shape is a list of positive integers,
        which corresponds to the number of tiles on each row of the snake graph.
        Recall that the shape of a snake graph uniquely determines a snake graph.

        EXAMPLES::

            sage: Gs = SnakeGraphs(9)
            sage: Gs((3,3,3))
                         -- -- --
                        |  |  |  |
                   -- -- -- -- --
                  |  |  |  |
             -- -- -- -- --
            |  |  |  |
             -- -- --
            sage: Gs((3,3,3)).shape()
            [3, 3, 3]
        """
        return self._shape

    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: Gs = SnakeGraphs(9)
            sage: G = Gs((3,2,4))
            sage: G == Gs([3,2,4])
            True
            sage: G == Gs([3,3,3])
            False
            sage: G == 'I am a string'
            False

        """
        if isinstance(other, SnakeGraph):
            return self._shape == other._shape
        return False

    def __ne__(self, other):
        """
        Check not equals. This is needed because otherwise != gives a wrong result.

        EXAMPLES::

            sage: SnakeGraph([1,1,1])==SnakeGraph([3])
            False
            sage: SnakeGraph([1,1,1])!=SnakeGraph([3])
            True
        """
        return not self.__eq__(other)

    def directions(self):
        """
        Return the list DIRs of directions (either 'up' or 'right').
        This list is of length `len(self)-1` and corresponds to all
        the tiles of ``self`` except for the last tile.

        Recall that we build a snake graph by starting with one tile, then glue
        a new tile so that the new tile is glued to the north or the east of
        the previous tile (see [MSW_Positivity] or [CanakciSchiffler]_).
        The entry in position `k` the list DIRs is `up`
        if the tile `k+1` is glued above tile `k`,
        and `right` if the tile `k+1` is glued to the right of tile `k`.

        EXAMPLES::

            sage: G = SnakeGraph([1,3,3,1,2,4,2])
            sage: G.directions()
            ['up',
            'right',
            'right',
            'up',
            'right',
            'right',
            'up',
            'up',
            'right',
            'up',
            'right',
            'right',
            'right',
            'up',
            'right']

            sage: SnakeGraph([1]).directions()
            []
        """
        temp_shape = self._shape[:]
        temp_shape.reverse()
        DIRs = []
        for i in range(len(self._shape)):
            r = temp_shape.pop()
            DIRs.extend((r-1)*['right'])
            if i == len(self._shape)-1:
                break
            DIRs.append('up')
        return DIRs

    def plot(self, rgb_color=(0,0,0), xy=(0, 0)):
        """
        Return a plot of ``self``.

        INPUT:

        - ``rgb_color`` -- (default:(0,0,0), black) The color as an RGB tuple
        - ``xy`` -- (default:(0,0)) snake graph will be plotted at xy=(a,b)

        EXAMPLES::

            sage: g=SnakeGraph((2,3,1))
            sage: print g.plot().description()
            Line defined by 5 points:   [(1.0, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]
            Line defined by 5 points:   [(2.0, 0.0), (1.0, 0.0), (1.0, 1.0), (2.0, 1.0), (2.0, 0.0)]
            Line defined by 5 points:   [(2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)]
            Line defined by 5 points:   [(3.0, 1.0), (2.0, 1.0), (2.0, 2.0), (3.0, 2.0), (3.0, 1.0)]
            Line defined by 5 points:   [(4.0, 1.0), (3.0, 1.0), (3.0, 2.0), (4.0, 2.0), (4.0, 1.0)]
        """
        from sage.plot.graphics import Graphics
        from sage.plot.line import line

        DIRs = self.directions()[:]

        drawing = Graphics()
        x, y = 0,0
        (x,y)=xy

        for pos in range(0,len(DIRs)):

            tile_drawing = line([(x+1,y+0),(x+0,y+0),(x+0,y+1),(x+1,y+1),(x+1,y+0)],rgbcolor=rgb_color)

            DIR = DIRs[pos]
            if DIR == 'up':
                y=y+1
            else:
                x=x+1

            drawing = drawing + tile_drawing
            drawing.axes(False)
            drawing.set_aspect_ratio(1)

        return drawing

class SnakeGraphs(Parent, UniqueRepresentation):
    """
    Class of all snake graphs with `d` square tiles.

    A snake graph is a connected sequence of square tiles.
    To build a snake graph, start with one tile, then glue a new tile so that
    the new tile is glued to the north or the east of the previous tile.
    See [MSW_Positivity] or [CanakciSchiffler]_.

    Note that the edges of the graph are not labeled. Hence snake graphs with `d`
    tiles are in bijection with :class:`Compositions` (of positive integers)
    with total sum `d`

    .. SEEALSO::

        :class:SnakeGraph

    INPUT:

    - ``d`` -- the number of tiles

    EXAMPLES:

    Here are the two types of snake graphs that can be created::

        sage: M = SnakeGraphs(2)
        sage: list(M)
        [
         --
        |  |
         --    -- --
        |  |  |  |  |
         -- ,  -- --
        ]

    REFERENCES:

    """
    def __init__(self, d):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: M = SnakeGraphs(2)
            sage: TestSuite(M).run()
        """
        self._d = d
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SnakeGraphs(2)
            Snake graphs with 2 tiles
        """
        return "Snake graphs with {} tiles".format(self._d)

    def _repr_option(self, key):
        """
        Metadata about the ``_repr_()`` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: Gs = SnakeGraphs(2)
            sage: Gs._repr_option('element_ascii_art')
            True
        """
        if key == 'element_ascii_art':
            return True
        return Parent._repr_option(self, key)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: Gs = SnakeGraphs(6)
            sage: g = Gs((2,1,3))
            sage: g
                -- -- --
               |  |  |  |
                -- -- --
               |  |
             -- --
            |  |  |
             -- --

            sage: Gs(g)
                -- -- --
               |  |  |  |
                -- -- --
               |  |
             -- --
            |  |  |
             -- --

            sage: Gs((2,2))
            Traceback (most recent call last):
            ...
            ValueError: Input a composition of 6

            sage: Gs((2,2))
            Traceback (most recent call last):
            ...
            ValueError: Input a composition of 6

            sage: Gs(matrix([1,1]))
            Traceback (most recent call last):
            ...
            ValueError: [1 1] is not a SnakeGraph nor a list of positive integers
        """
        if isinstance(x, SnakeGraph):
            if x in self.parent():
                return x
            else:
                raise ValueError("Cannot convert between Snake Graphs of different number of tiles")
        elif isinstance(x, list) or isinstance(x, tuple) or isinstance(x, set):
            if sum(x) == self._d:
                return self.element_class(self, x)
            else:
                raise ValueError("Input a composition of {}".format(self._d))
        else:
            raise ValueError("{} is not a SnakeGraph nor a list of positive integers".format(x))

    Element = SnakeGraph

    def __iter__(self):
        """
        Iterate through ``self``.

        EXAMPLES::

            sage: M = SnakeGraphs(2)
            sage: list(M)
            [
             --
            |  |
             --    -- --
            |  |  |  |  |
             -- ,  -- --
            ]
        """
        from sage.combinat.composition import Compositions
        for c in Compositions(self._d):
            yield self.element_class(self,c)

    def number_of_tiles(self):
        """
        Return the number of tiles of the snake graph.

        EXAMPLES::

            sage: M = SnakeGraphs(4)
            sage: M.number_of_tiles()
            4
        """
        return self._d