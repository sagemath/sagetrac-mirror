r"""
Tikz Painter
=========================

The goal of this module is to give some tools to produce drawing in the Tikz
Format.
"""
# *****************************************************************************
#  Copyright (C) 2018,2019 Adrien Boussicault (boussica@labri.fr),
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from __future__ import (
    division, absolute_import, print_function, unicode_literals
)
# from sage.functions.trig import cos, sin
#from sage.functions.other import sqrt
from math import cos, sin, sqrt


class TikzPainter:
    r"""
    Technical class to produce TIKZ drawing.

    This class contains some 2D geometric tools to produce some TIKZ
    drawings.

    With that classes you can use options to set up drawing informations.
    Then the class will produce a drawing by using those informations.

    EXAMPLES::

        sage: dt = TikzPainter()
        sage: dt.draw_line([1, 1], [-1, -1])
        u'\n  \\draw[color=black, line width=1] (1.000000, 1.000000) --
        (-1.000000, -1.000000);'

        sage: fct = lambda vec: [2*vec[0], vec[1]]
        sage: import copy
        sage: opt = copy.deepcopy(TikzPainter.default_options)
        sage: opt['coordinate_transformation'] = fct
        sage: dt = TikzPainter(options=opt)
        sage: dt.draw_line([1, 1], [-1, -1])
        u'\n  \\draw[color=black, line width=1] (2.000000, 1.000000) --
        (-2.000000, -1.000000);'

        sage: opt = copy.deepcopy(TikzPainter.default_options)
        sage: opt['mirror'] = [0,1]
        sage: dt = TikzPainter(opt)
        sage: dt.draw_line([1, 1], [-1, -1])
        u'\n  \\draw[color=black, line width=1] (-1.000000, 1.000000) --
        (1.000000, -1.000000);'
    """
    class identity_transformation:
        def __call__(self,v):
            return v
        def __repr__(self):
            return "lambda v : v"
    default_options = dict(
        scale=1, line_size=1, point_size=3.5, color_line='black',
        color_point='black', translation=[0, 0], rotation=0, mirror=None,
        coordinate_transformation=identity_transformation()
    )
    r"""
    This is the default TIKZ options.

    This option is used to configurate element of a drawing to allow
    TIKZ code generation.
    """

    def __init__(self, options=None):
        r"""
        Construct a drawing tools to produce some TIKZ drawing.

        INPUTS:

        - ``options`` -- drawing options

        EXAMPLES::

            sage: dt = TikzPainter()
            sage: dt.draw_line([1, 1], [-1, -1])
            u'\n  \\draw[color=black, line width=1] (1.000000, 1.000000) --
            (-1.000000, -1.000000);'

            sage: import copy
            sage: opt = copy.deepcopy(TikzPainter.default_options)
            sage: opt['rotation'] = 2*pi/4
            sage: dt = TikzPainter(opt)
            sage: dt.draw_line([1, 1], [-1, -1])
            u'\n  \\draw[color=black, line width=1] (-1.000000, 1.000000) --
            (1.000000, -1.000000);'

            sage: opt = dict(scale=3, line_size=2)
            sage: dt = TikzPainter(opt)
            sage: dt.draw_line([1, 1], [-1, -1])
            u'\n  \\draw[color=black, line width=2] (3.000000, 3.000000) --
            (-3.000000, -3.000000);'

            sage: TikzPainter.default_options
            {'color_line': u'black',
             'color_point': u'black',
             'coordinate_transformation': lambda v : v,
             'line_size': 1,
             'mirror': None,
             'point_size': 3.5,
             'rotation': 0,
             'scale': 1,
             'translation': [0, 0]}
        """
        def get_option( key, options ):
            if options is None or not key in options :
                return TikzPainter.default_options[key]
            else:
                return options[key]
        ct = get_option('coordinate_transformation', options)
        def XY(v):
            v = ct(v)
            return [float(v[0]), float(v[1])]
        self._XY = XY
        self._translation = get_option('translation', options)
        self._mirror = get_option('mirror', options)
        self._rotation = get_option('rotation', options)
        self._color_line = get_option('color_line', options)
        self._line_size = get_option('line_size', options)
        self._point_size = get_option('point_size', options)
        self._color_point = get_option('color_point', options)
        self._scale = get_option('scale', options)

    def XY(self, v):
        r"""
        This function give the image of v by some transformation given by the
        drawing option of ``TikzPainter``.

        The transformation is the composition of rotation, mirror, translation
        and XY user function.

        First we apply XY function, then the translation, then the mirror and
        finaly the rotation.

        INPUT:

        - ``v`` -- The vector to transform.

        OUTPUT:

        A list of 2 floats encoding a vector.

        EXAMPLES::

            sage: dt = TikzPainter()
            sage: dt.XY([1, 1])
            [1.0, 1.0]

            sage: fct = lambda vec: [2*vec[0], vec[1]]
            sage: import copy
            sage: opt = copy.deepcopy(TikzPainter.default_options)
            sage: opt['coordinate_transformation'] = fct
            sage: dt = TikzPainter(options=opt)
            sage: dt.XY([1, 1])
            [2.0, 1.0]

            sage: import copy
            sage: opt = copy.deepcopy(TikzPainter.default_options)
            sage: opt['mirror'] = [0, 1]
            sage: dt = TikzPainter(opt)
            sage: dt.XY([1, 1])
            [-1.0, 1.0]
        """
        def translate(pos, v):
            r"""
            Translate a position with a vector.

            INPUTS:

            - ``pos`` -- The position to translate.

            - ``v`` -- The translation vector.

            OUTPUT:

            The translated position.
            """
            return [pos[0]+v[0], pos[1]+v[1]]

        def rotate(pos, angle):
            r"""
            Rotate by `angle` a position around the origin.

            INPUT:

            - ``pos`` -- The position to rotate.

            - ``angle`` -- The angle of rotation (in radian).

            OUTPUT:

            The rotated position.
            """
            [x, y] = pos
            return [x*cos(angle) - y*sin(angle), x*sin(angle) + y*cos(angle)]

        def scale(pos, scale):
            r"""
            Scale by `scale` a position around the origin.

            INPUT:

            - ``pos`` -- The position to scale.

            - ``scale`` -- The scale of rotation (in radian).

            OUTPUT:

            The scaled position.
            """
            [x, y] = pos
            return [scale*x, scale*y]

        def mirror(pos, axe):
            r"""
            Return the mirror of a position according to a given axe.

            INPUTS:

            - ``pos`` -- The position to mirror.

            - ``axe`` -- The axe vector.

            OUTPUT:

            The mirrored position.
            """
            if axe is None:
                return pos
            if not isinstance(axe, (list, tuple)):
                raise ValueError(
                    "mirror option should be None or a list of two real" +
                    " encoding a 2D vector."
                )
            n = float(sqrt(axe[0]**2 + axe[1]**2))
            axe[0] = float(axe[0]/n)
            axe[1] = float(axe[1]/n)
            sp = (pos[0]*axe[0] + pos[1]*axe[1])
            sn = (- pos[0]*axe[1] + pos[1]*axe[0])
            return [
                sp*axe[0] + sn*axe[1],
                sp*axe[1] - sn*axe[0]
            ]
        # TODO : MAKE tests with roatations AND translation ! 
        return translate(
            mirror(
                rotate(
                    scale(
                        self._XY(v),
                        self._scale
                    ),
                    self._rotation
                ),
                self._mirror
            ),
            self._translation
        )

    def draw_line(self, v1, v2, color=None, size=None):
        r"""
        Return the TIKZ code for a line.

        INPUT:

        - ``v1`` -- point, The first point of the line.

        - ``v2`` -- point, The second point of the line.

        - ``color`` -- string (default:``None``), The color of the line.
          If set to ``None``, the color is choosen according the
          drawing option given by ``TikzPainter``.

        - ``size`` -- integer (default:``None``), The size of the line.
          If set to ``None``, the size is choosen according the
          drawing option given by ``TikzPainter``.

        OUTPUT:

        The code of a line in TIKZ.

        EXAMPLES::

            sage: dt = TikzPainter()
            sage: dt.draw_line([1, 1], [-1, -1])
            u'\n  \\draw[color=black, line width=1] (1.000000, 1.000000) --
            (-1.000000, -1.000000);'

        """
        if color is None:
            color = self._color_line
        if size is None:
            size = self._line_size
        [x1, y1] = self.XY(v1)
        [x2, y2] = self.XY(v2)
        return "\n  \\draw[color=%s, line width=%s] (%f, %f) -- (%f, %f);" % (
            color, size, float(x1), float(y1), float(x2), float(y2)
        )

    def draw_polyline(self, list_of_vertices, color=None, size=None):
        r"""
        Return the TIKZ code for a polyline.

        INPUT:

        - ``list_of_vertices`` -- A list of points

        - ``color`` -- string (default:``None``), The color of the line.
          If set to ``None``, the color is choosen according the
          drawing option given by ``TikzPainter``.

        - ``size`` -- integer (default:``None``), The size of the line.
          If set to ``None``, the size is choosen according the
          drawing option given by ``TikzPainter``.

        OUTPUT:

        The code of a polyline in TIKZ.

        EXAMPLES::

            sage: dt = TikzPainter()
            sage: dt.draw_polyline([[1, 1], [-1, -1], [0,0]])
            u'\n  \\draw[color=black, line width=1] (1.000000, 1.000000) --
            (-1.000000, -1.000000);\n  \\draw[color=black, line width=1]
            (-1.000000, -1.000000) -- (0.000000, 0.000000);'
        """
        res = ""
        for i in range(len(list_of_vertices)-1):
            res += self.draw_line(
                list_of_vertices[i], list_of_vertices[i+1], color, size
            )
        return res

    def draw_point(self, p1, color=None, size=None):
        r"""
        Return the TIKZ code for a point.


        INPUT:

        - ``p1`` -- A point

        - ``color`` -- string (default:``None``), The color of the line.
          If set to ``None``, the color is choosen according the
          drawing option given by ``TikzPainter``.

        - ``size`` -- integer (default:``None``), The size of the line.
          If set to ``None``, the size is choosen according the
          drawing option given by ``TikzPainter``.

        OUTPUT:

        The code of a point in TIKZ.

        EXAMPLES::

            sage: dt = TikzPainter()
            sage: dt.draw_point([1, 1])
            u'\n  \\filldraw[color=black] (1.000000, 1.000000) circle (3.5pt);'

        """
        if color is None:
            color = self._color_point
        if size is None:
            size = self._point_size
        [x1, y1] = self.XY(p1)
        return "\n  \\filldraw[color=%s] (%f, %f) circle (%spt);" % (
            color, float(x1), float(y1), size
        )
