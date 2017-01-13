r"""
Snake Graphs

REFERENCES:

.. [CanakciSchiffler] Canakci and Schiffler,
   *Snake graph calculus and cluster algebras from surfaces*
   :arxiv:`1209.4617`

.. [MSW_Positivity] Musiker - Schiffler - Williams,
   *Positivity for Cluster Algebras from Surfaces*,
   :arxiv:`0906.0748`

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
from sage.categories.sets_cat import Sets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.element import Element
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.list_clone import ClonableArray
from sage.rings.all import ZZ
from sage.combinat.composition import composition_iterator_fast

class SnakeGraph(ClonableArray):
    """
    A snake graph is a connected sequence of square tiles.

    A snake graph is a connected sequence of square tiles. A square tile
    is considered as a graph with four vertices and four edges in the
    obvious way, and we use the following ascii art to visualize each tile::

             north
               --
        west  |  | east
               --
             south

    To build a snake graph, start with one tile, then glue a new tile so that
    the new tile is glued to the north or the east of the previous tile.
    See [MSW_Positivity]_ or [CanakciSchiffler]_.

    Note that the edges of the graph are not labeled. Hence a snake graph
    is uniquely determined by a list of positive integers (``shape``) such
    that their sum is equal to the number of the snake graph's tiles, i.e.
    a snake graph is uniquely determined by a composition of ``d``, where
    ``d`` is the number of tiles of the snake graph.

    INPUT:

    - ``shape`` -- a tuple/list listing the sizes of the rows of
      the snake graph

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

        TESTS::

            sage: Gs = SnakeGraphs(9)
            sage: G = Gs((3,2,4))
            sage: G == Gs([3,2,4])
            True
            sage: G == Gs([3,3,3])
            False
            sage: G == 'I am a string'
            False
            sage: SnakeGraph([1,1,1]) == SnakeGraph([3])
            False
            sage: SnakeGraph([1,1,1]) != SnakeGraph([3])
            True

            sage: G = SnakeGraph((2,1,1))
            sage: TestSuite(G).run()
        """
        SGs = SnakeGraphs(sum(shape))
        return SGs(shape)#, first_tile_orientation, edge_weights, diagonal_weights)

    def check(self):
        """
        Check if ``self`` is a valid snake graph.

        EXAMPLES::

            sage: M = SnakeGraphs(3)
            sage: M[0].check()
        """
        if any(x not in ZZ or x <= 0 for x in self):
            raise ValueError("the snake graph must consist of positive integers")

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
        ret = ''
        top_row = self[-1]
        skips = sum(self[:-1])-(len(self)-1)
        white_sp = '   '

        ret += white_sp* skips
        ret +=' -- '
        for i in range(1, top_row):
            ret +='-- '
        ret +='\n' + white_sp * skips + '|  |'
        for i in range(1, top_row):
            ret +='  |'

        for i in range(len(self)-2,-1,-1):
            r = self[i]
            skips += -(r-1)

            ret +='\n' + white_sp * skips
            for i in range(r + self[i+1] - 1):
                ret +=' --'

            ret +='\n' + white_sp * skips + '|  |'
            for i in range(1, r):
                ret +='  |'

        bottom_row = self[0]
        ret +='\n' + ' --'
        for i in range(1, bottom_row):
            ret +=' --'

        return ret

    def shape(self):
        """
        Return the shape of ``self``.

        The shape is a list of positive integers, which corresponds to
        the number of tiles on each row of the snake graph. Recall that
        the shape of a snake graph uniquely determines a snake graph.

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
        return list(self)

    def directions(self):
        """
        Return the list of directions (either 'up' or 'right').

        This list is of length `len(self)-1` and corresponds to all
        the tiles of ``self`` except for the last tile.

        Recall that we build a snake graph by starting with one tile,
        then glue a new tile so that the new tile is glued to the north
        or the east of the previous tile (see [MSW_Positivity]_ or
        [CanakciSchiffler]_). The entry in position `k` the list DIRs
        is ``'up'`` if the tile `k+1` is glued above tile `k`, and
        ``'right'`` if the tile `k+1` is glued to the right of tile `k`.

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
        temp_shape = self[:]
        temp_shape.reverse()
        DIRs = []
        for i in range(len(self)):
            r = temp_shape.pop()
            DIRs.extend((r-1)*['right'])
            if i == len(self)-1:
                break
            DIRs.append('up')
        return DIRs

    def plot(self, rgb_color=(0,0,0), xy=(0,0)):
        """
        Return a plot of ``self``.

        INPUT:

        - ``rgb_color`` -- (default:(0,0,0), black) The color as an RGB tuple
        - ``xy`` -- (default:(0,0)) snake graph will be plotted at xy=(a,b)

        EXAMPLES::

            sage: g=SnakeGraph((2,3,1))
            sage: print g.plot(rgb_color=(1,0,1)).description()
            Line defined by 5 points:       [(1.0, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]
            Line defined by 5 points:       [(2.0, 0.0), (1.0, 0.0), (1.0, 1.0), (2.0, 1.0), (2.0, 0.0)]
            Line defined by 5 points:       [(2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)]
            Line defined by 5 points:       [(3.0, 1.0), (2.0, 1.0), (2.0, 2.0), (3.0, 2.0), (3.0, 1.0)]
            Line defined by 5 points:       [(4.0, 1.0), (3.0, 1.0), (3.0, 2.0), (4.0, 2.0), (4.0, 1.0)]
            Line defined by 5 points:       [(4.0, 2.0), (3.0, 2.0), (3.0, 3.0), (4.0, 3.0), (4.0, 2.0)]
        """
        from sage.plot.graphics import Graphics
        from sage.plot.line import line

        DIRs = self.directions()[:]

        drawing = Graphics()
        x,y = 0,0
        (x,y) = xy

        for pos in range(0,len(DIRs)+1):

            tile_drawing = line([(x+1,y+0),(x+0,y+0),(x+0,y+1),
                                 (x+1,y+1),(x+1,y+0)],rgbcolor=rgb_color)

            if pos < len(DIRs):
                DIR = DIRs[pos]

            if DIR == 'up':
                y = y + 1
            else:
                x = x + 1

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
    See [MSW_Positivity]_ or [CanakciSchiffler]_.

    Note that the edges of the graph are not labeled. Hence snake graphs
    with `d` tiles are in bijection with :class:`Compositions` (of
    positive integers) with total sum `d`.

    .. SEEALSO::

        :class:`SnakeGraph`

    INPUT:

    - ``d`` -- the number of tiles

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
            ValueError: input a composition of 6

            sage: Gs(matrix([1,1]))
            Traceback (most recent call last):
            ...
            ValueError: [1 1] is not a snake graph nor a list of positive integers
        """
        if isinstance(x, SnakeGraph):
            if x in self.parent():
                return x
            else:
                raise ValueError("cannot convert between snake graphs of different number of tiles")
        elif isinstance(x, list) or isinstance(x, tuple) or isinstance(x, set):
            if sum(x) == self._d:
                return self.element_class(self, x)
            else:
                raise ValueError("input a composition of {}".format(self._d))
        else:
            raise ValueError("{} is not a snake graph nor a list of positive integers".format(x))

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
        for c in composition_iterator_fast(self._d):
            yield self.element_class(self, c)

    def number_of_tiles(self):
        """
        Return the number of tiles of the snake graphs of ``self``.

        EXAMPLES::

            sage: M = SnakeGraphs(4)
            sage: M.number_of_tiles()
            4
        """
        return self._d

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: M = SnakeGraphs(4)
            sage: M.cardinality()
            16
        """
        return ZZ(2)**self._d

class LabeledSnakeGraph(SnakeGraph):
    r"""
    A labeled snake graph.

    A *labeled snake graph* is a snake graph in which each edge and each
    tile carries a label or weight [CanakciSchiffler]_. For example, for
    snake graphs arising from cluster algebras from surfaces
    [FominShapiroThurston]_, these labels are cluster variables.
    See [MSW_Positivity]_.

    In some situation we would like to consider the weights of the diagonals of the
    snake graph tiles. The diagonal of a tile is defined to be an edge between the
    northwest (NW) and southeast (SE) corners of a tile::

        NW
          --
         |\ |
         | \|
          --
            SE

    .. SEEALSO::

        :class:`SnakeGraph`

    Note that :class:`LabeledSnakeGraph` differs from :class:`SnakeGraph`
    in that user may specify two optional attributes ``diagonal_weights``
    and ``weights``.

    INPUT:

    - ``shape`` -- a tuple/list listing the sizes of the rows of
      the snake graph
    - ``weights`` -- (default: ``None``) a list/tuple/dictionary giving
      the weight of each edge of the snake graph
    - ``diagonal_weights`` -- (default: ``None``) a list/tuple/dictionary
      giving the weight for the diagonal of each tile
    - ``first_tile_orientation`` -- (default: 1) whether the orientation
      of the first tile is `1` or `-1`

    EXAMPLES::

        sage: LabeledSnakeGraph((2,1,3))
            -- -- --
           |  |  |  |
            -- -- --
           |  |
         -- --
        |  |  |
         -- --
         with edge labels

        sage: LabeledSnakeGraph((2,1),weights={0:('B1','B2','b','d'),\
        ....: 1:('a','c','B3','a'),2:('B3','B4','d','b')})
            --
           |  |
         -- --
        |  |  |
         -- --
        with edge labels

        sage: L = LabeledSnakeGraph([2,1],'my weight')
        Traceback (most recent call last):
        ...
        ValueError: weights must be a dictionary of length 3
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, shape, weights={}, diagonal_weights={},
                 first_tile_orientation=1, from_surface=False):
        return LabeledSnakeGraphs()(shape, weights, diagonal_weights,
                 first_tile_orientation, from_surface)

    def __init__(self, parent, shape, weights, diagonal_weights,
                 first_tile_orientation, from_surface):
        """
        Initialize ``self``.

        .. TODO::

            For the input weights, we should check that if edge `e` belongs
            to two tiles, the weight of `e` is consistent on both tiles.

        TESTS::

            sage: G = LabeledSnakeGraph((2,1,1))
            sage: TestSuite(G).run()
        """
        SnakeGraph.__init__(self, parent, shape)
        self._weights = weights
        self._diagonal_weights = diagonal_weights
        self._first_tile_orientation = 1
        self._from_surface = from_surface # TODO
        if weights:
            num_of_tiles = sum(shape)
            if not isinstance(weights,dict) or len(weights) != num_of_tiles:
                raise ValueError("weights must be a dictionary of length {}".format(num_of_tiles))

            if weights.keys() != range(num_of_tiles) or\
            not all(isinstance(v,(list,tuple,dict)) and len(v)==4 for v in weights.values()):
                raise ValueError("weights must be in the form ",\
                "{0:(S0,E0,N0,W0), 1:(S1,E1,N1,W1),..., d:(Sn,En,Nn,Wn)}, ",\
                " where d is {}".format(num_of_tiles-1))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LabeledSnakeGraph([1,1,2,1], weights={0:(1,2,3,4),1:(5,6,7,8),2:(9,8,7,6),
            ....:                                       3:(1,1,1,1),4:(2,2,2,2)})
                --
               |  |
             -- --
            |  |  |
             -- --
            |  |
             --
            |  |
             --
            with edge labels
        """
        return SnakeGraph._repr_(self) + '\nwith edge labels'

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
        if isinstance(other, LabeledSnakeGraph):
            return (list(self) == list(other)
                    and self._weights == other._weights
                    and self._diagonal_weights == other._diagonal_weights)
        return False

    def __ne__(self, other):
        """
        Check not equals. This is needed because otherwise != gives a wrong result.
        [TODO] maybe this is not needed here because this inherits __ne__from SnakeGraph

        EXAMPLES::

            sage: LabeledSnakeGraph([1,1,1]) == LabeledSnakeGraph([3])
            False
            sage: LabeledSnakeGraph([1,1,1]) != LabeledSnakeGraph([3])
            True

            sage: LabeledSnakeGraph([2,1],{0:('B1','B2','b','d'),1:('a','c','B3','a'),
            ....: 2:('B3','B4','d','b')}) == LabeledSnakeGraph([2,1])
            False
            sage: LabeledSnakeGraph([2,1],{0:('B1','B2','b','d'),1:('a','c','B3','a'),
            ....: 2:('B3','B4','d','b')}) == SnakeGraph([2,1])
            False
            sage: LabeledSnakeGraph([2,1],{0:('B1','B2','b','d'),1:('a','c','B3','a'),
            ....: 2:('B3','B4','d','b')}) != LabeledSnakeGraph([2,1])
            True
            sage: LabeledSnakeGraph([2,1],{0:('B1','B2','b','d'),1:('a','c','B3','a'),
            ....: 2:('B3','B4','d','b')}) != SnakeGraph([2,1])
            True

            sage: LabeledSnakeGraph([2,1]) == SnakeGraph([2,1])
            False
            sage: LabeledSnakeGraph([2,1]) != SnakeGraph([2,1])
            True
        """
        return not self.__eq__(other)

    def weights(self):
        """
        Return the weights for the edges of the snake graph.

        EXAMPLES::

           sage: L = LabeledSnakeGraph([1],{0:('south','east','north','west')})
           sage: L.weights()
           {0: ('south', 'east', 'north', 'west')}
        """
        return self._weights

    def plot(self, rgb_color=(0,0,0), xy=(0,0), draw_weights=True,
             draw_diagonal_weights=True, text_color=(1,0,0)):
        """
        Return a plot of ``self``.

        INPUT:

        - ``rgb_color`` -- (default: ``(0,0,0)``, black) the color as
          an RGB tuple
        - ``xy`` -- (default: ``(0,0)``) the coordinates to start plotting
          the snake graph
        - ``text_color`` -- (default: ``(1,0,0)``, red) the color of the
          edge labels

        EXAMPLES::

            sage: L = LabeledSnakeGraph([2])
            sage: print L.plot().description()
            Line defined by 5 points:       [(1.0, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]
            Line defined by 5 points:       [(2.0, 0.0), (1.0, 0.0), (1.0, 1.0), (2.0, 1.0), (2.0, 0.0)]
            Text '$+$' at the point (0.8,0.8)
            Text '$-$' at the point (1.8,0.8)
            Text '' at the point (0.0,0.5)
            Text '' at the point (0.5,0.0)
            Text '' at the point (0.5,0.5)
            Text '' at the point (0.5,1.0)
            Text '' at the point (1.0,0.5)
            Text '' at the point (1.5,0.0)
            Text '' at the point (1.5,0.5)
            Text '' at the point (1.5,1.0)
            Text '' at the point (2.0,0.5)

            sage: L=LabeledSnakeGraph([2,1],{0:('B1','B2','b','d'),1:('a','c','B3','a'),
            ....: 2:('B3','B4','d','b')})
            sage: print L.plot().description()
            Line defined by 5 points:       [(1.0, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]
            Line defined by 5 points:       [(2.0, 0.0), (1.0, 0.0), (1.0, 1.0), (2.0, 1.0), (2.0, 0.0)]
            Line defined by 5 points:       [(2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)]
            Text '$+$' at the point (0.8,0.8)
            Text '$+$' at the point (1.8,1.8)
            Text '$-$' at the point (1.8,0.8)
            Text '' at the point (0.5,0.5)
            Text '' at the point (1.5,0.5)
            Text '' at the point (1.5,1.5)
            Text 'B1' at the point (0.5,0.0)
            Text 'B2' at the point (1.0,0.5)
            Text 'B3' at the point (1.5,1.0)
            Text 'B4' at the point (2.0,1.5)
            Text 'a' at the point (1.5,0.0)
            Text 'b' at the point (0.5,1.0)
            Text 'b' at the point (1.0,1.5)
            Text 'c' at the point (2.0,0.5)
            Text 'd' at the point (0.0,0.5)
            Text 'd' at the point (1.5,2.0)

            sage: L=LabeledSnakeGraph([2,1],{0:('B1','B2','b','d'),1:('a','c','B3','a'),
            ....: 2:('B3','B4','d','b')},{0:'a',1:'b',2:'c'})
            sage: print L.plot().description()
            Line defined by 5 points:       [(1.0, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]
            Line defined by 5 points:       [(2.0, 0.0), (1.0, 0.0), (1.0, 1.0), (2.0, 1.0), (2.0, 0.0)]
            Line defined by 5 points:       [(2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)]
            Text '$+$' at the point (0.8,0.8)
            Text '$+$' at the point (1.8,1.8)
            Text '$-$' at the point (1.8,0.8)
            Text 'B1' at the point (0.5,0.0)
            Text 'B2' at the point (1.0,0.5)
            Text 'B3' at the point (1.5,1.0)
            Text 'B4' at the point (2.0,1.5)
            Text 'a' at the point (0.5,0.5)
            Text 'a' at the point (1.5,0.0)
            Text 'b' at the point (0.5,1.0)
            Text 'b' at the point (1.0,1.5)
            Text 'b' at the point (1.5,0.5)
            Text 'c' at the point (1.5,1.5)
            Text 'c' at the point (2.0,0.5)
            Text 'd' at the point (0.0,0.5)
            Text 'd' at the point (1.5,2.0)
        """
        if not draw_weights:
            return SnakeGraph(self).plot()

        from sage.plot.graphics import Graphics
        from sage.plot.line import line
        from sage.plot.text import text

        DIRs = self.directions()[:]

        drawing = Graphics()
        x, y = 0, 0
        (x,y) = xy

        for pos in range(0,len(DIRs)+1):
            tile_drawing = line([(x+1,y+0),(x+0,y+0),(x+0,y+1),
                                 (x+1,y+1),(x+1,y+0)], rgbcolor=rgb_color)

            if self._weights:
                floor = str(self._weights[pos][0])
                right_side = str(self._weights[pos][1])
                ceiling = str(self._weights[pos][2])
                left_side = str(self._weights[pos][3])
            else:
                floor, right_side, ceiling, left_side = '', '', '', ''

            if self._diagonal_weights:
                diagonal = str(self._diagonal_weights[pos])
            else:
                diagonal = ''

            if self._first_tile_orientation == 1:
                orientation = (-1)**pos
            elif self._first_tile_orientation == -1:
                orientation = -1*(-1)**pos
            if orientation == 1: orientation='$+$'
            else: orientation='$-$'

            labels = ( text(diagonal, (x+0.5,y+0.5), vertical_alignment='bottom', rgbcolor=text_color)
                      + text(right_side, (x+1,y+0.5), horizontal_alignment='left', rgbcolor=text_color)
                      + text(ceiling, (x+0.5,y+1), vertical_alignment='bottom', rgbcolor=text_color)
                      + text(orientation, (x+0.8, y+0.8)) )

            if pos > 0:
                PREVIOUS_DIR = DIRs[pos-1]
                if PREVIOUS_DIR == 'right': # Then draw the label of the bottom edge
                    labels = labels + text(floor,(x+0.5,y+0),vertical_alignment='bottom', rgbcolor=text_color)
                elif PREVIOUS_DIR == 'up': # Then draw the label of the left edge
                    labels = labels + text(left_side,(x+0,y+0.5),horizontal_alignment='left', rgbcolor=text_color)
            else: # For the first tile, draw labels for both bottom and left edges
                labels = labels + \
                text(floor,(x+0.5,y+0),vertical_alignment='bottom', rgbcolor=text_color)\
                + text(left_side,(x+0,y+0.5),horizontal_alignment='left', rgbcolor=text_color)

            if pos < len(DIRs):
                DIR = DIRs[pos]
            if DIR == 'up':
                y = y + 1
            else:
                x = x + 1

            drawing = drawing + tile_drawing + labels
            drawing.axes(False)
            drawing.set_aspect_ratio(1)

        return drawing

class LabeledSnakeGraphs(Parent, UniqueRepresentation):
    Element = LabeledSnakeGraph

    def __init__(self):
        Parent.__init__(self, category=Sets().Infinite())

    def _element_constructor_(self, shape, weights={}, diagonal_weights={},
                 first_tile_orientation=1, from_surface=False):
        return self.element_class(self, shape, weights, diagonal_weights,
                                  first_tile_orientation, from_surface)

