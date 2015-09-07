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

RIGHT='right'
UP='up'

class SnakeGraph(Element):
    """
    A snake graph is a connected sequence of square tiles which goes north and east,
    see [CanakciSchiffler]_.

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
        Element.__init__(self, parent)

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
        Return shape of the snake graph.

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
            sage: G == AlternatingSignMatrix([[1, 0, 0],[0, 0, 1],[0, 1, 0]])
            False
        """
        if isinstance(other, SnakeGraph):
            return self._shape == other._shape
        return False

    def directions(self):
        """
        Return the list of directions for each tile (except the final tile),
        thinking of the graph as starting from
        the southwest corner to the northeast corner.

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
        """
        temp_shape = self._shape[:]
        temp_shape.reverse()
        DIRs = []
        for i in range(len(self._shape)):
            r = temp_shape.pop()
            DIRs.extend((r-1)*[RIGHT])
            if i == len(self._shape)-1:
                break
            DIRs.append(UP)
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
            if DIR == UP:
                y=y+1
            else:
                x=x+1

            drawing = drawing + tile_drawing
            drawing.axes(False)
            drawing.set_aspect_ratio(1)

        return drawing

class SnakeGraphs(Parent, UniqueRepresentation):
    """
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
            The snake graphs with 2 tiles
        """
        return "The snake graphs with {} tiles".format(self._d)

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

            sage: M = SnakeGraphs(6)
            sage: M((2,1,3))
                -- -- --
               |  |  |  |
                -- -- --
               |  |
             -- --
            |  |  |
             -- --
        """
        if isinstance(x, SnakeGraph):
            if x in self.parent():
                return x
            else:
                print 'print an error message here TODO'
        else:
            return self.element_class(self, x)

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