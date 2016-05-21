"""
Text in plots
"""
#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options, rename_keyword
from sage.plot.colors import to_mpl_color

class Text(GraphicPrimitive):
    """
    Base class for Text graphics primitive.

    TESTS:

    We test creating some text::

        sage: text("I like Fibonacci",(3,5))
        Graphics object consisting of 1 graphics primitive


    """
    def __init__(self, string, point, options):
        """
        Initializes base class Text.

        EXAMPLES::

            sage: T = text("I like Fibonacci", (3,5))
            sage: t = T[0]
            sage: t.string
            'I like Fibonacci'
            sage: t.x
            3.0
            sage: t.options()['fontsize']
            10
        """
        self.string = string
        self.x = float(point[0])
        self.y = float(point[1])
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data. Notice
        that, for text, the box is just the location itself.

        EXAMPLES::

            sage: T = text("Where am I?",(1,1))
            sage: t=T[0]
            sage: t.get_minmax_data()['ymin']
            1.0
            sage: t.get_minmax_data()['ymax']
            1.0
        """
        from sage.plot.plot import minmax_data
        return minmax_data([self.x], [self.y], dict=True)

    def _repr_(self):
        """
        String representation of Text primitive.

        EXAMPLES::

            sage: T = text("I like cool constants", (pi,e))
            sage: t=T[0];t
            Text 'I like cool constants' at the point (3.14159265359,2.71828182846)
        """
        return "Text '%s' at the point (%s,%s)"%(self.string, self.x, self.y)

    def _allowed_options(self):
        """
        Return the allowed options for the Text class.

        EXAMPLES::

            sage: T = text("ABC",(1,1),zorder=3)
            sage: T[0]._allowed_options()['fontsize']
            "How big the text is. Either the size in points or a relative size, e.g. 'smaller', 'x-large', etc"
            sage: T[0]._allowed_options()['zorder']
            'The layer level in which to draw'
            sage: T[0]._allowed_options()['rotation']
            'how to rotate the text: angle in degrees, vertical, horizontal'
        """
        return {'fontsize': 'How big the text is. Either the size in points or a relative size, e.g. \'smaller\', \'x-large\', etc',
                'fontstyle': 'either \'normal\', \'italic\' or \'oblic\'',
                'fontweight': 'a numeric value in the range 0-1000 or a string'
                              '\'ultralight\', \'light\', \'normal\', \'regular\', \'book\','
                              '\'medium\', \'roman\', \'semibold\', \'demibold\', \'demi\','
                              '\'bold,\', \'heavy\', \'extra bold\', \'black\'',
                'rgbcolor':'The color as an RGB tuple.',
                'background_color': 'The background color.',
                'bounding_box': 'A dictionary specifying a bounding box',
                'hue':'The color given as a hue.',
                'alpha': 'a float (0.0 transparent through 1.0 opaque)',
                'axis_coords':'Uses axis coordinates -- (0,0) lower left and (1,1) upper right',
                'rotation': 'how to rotate the text: angle in degrees, vertical, horizontal',
                'vertical_alignment': 'how to align vertically: top, center, bottom',
                'horizontal_alignment':'how to align horizontally: left, center, right',
                'zorder':'The layer level in which to draw',
                'clip': 'Whether to clip or not.'}

    def _plot3d_options(self, options=None):
        """
        Translate 2D plot options into 3D plot options.

        EXAMPLES::

            sage: T = text("ABC",(1,1))
            sage: t = T[0]
            sage: t.options()['rgbcolor']
            (0.0, 0.0, 1.0)
            sage: s=t.plot3d()
            sage: s.jmol_repr(s.testing_render_params())[0][1]
            'color atom  [0,0,255]'
        """
        if options is None:
            options = dict(self.options())
        options_3d = {}
        # TODO: figure out how to implement rather than ignore
        for s in ['axis_coords', 'clip', 'fontsize', 'horizontal_alignment',
                'rotation', 'vertical_alignment' ]:
            if s in options:
                del options[s]
        options_3d.update(GraphicPrimitive._plot3d_options(self, options))
        return options_3d

    def plot3d(self, **kwds):
        """
        Plots 2D text in 3D.

        EXAMPLES::

            sage: T = text("ABC",(1,1))
            sage: t = T[0]
            sage: s=t.plot3d()
            sage: s.jmol_repr(s.testing_render_params())[0][2]
            'label "ABC"'
            sage: s._trans
            (1.0, 1.0, 0)
        """
        from sage.plot.plot3d.shapes2 import text3d
        options = self._plot3d_options()
        options.update(kwds)
        return text3d(self.string, (self.x, self.y, 0), **options)

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: t1 = text("Hello",(1,1), vertical_alignment="top", fontsize=30, rgbcolor='black')
            sage: t2 = text("World", (1,1), horizontal_alignment="left",fontsize=20, zorder=-1)
            sage: t1 + t2   # render the sum
            Graphics object consisting of 2 graphics primitives
        """
        options = self.options()
        opts = {}
        opts['color'] = options['rgbcolor']
        opts['verticalalignment'] = options['vertical_alignment']
        opts['horizontalalignment'] = options['horizontal_alignment']
        if 'background_color' in options:
            opts['backgroundcolor'] = options['background_color']
        if 'fontweight' in options:
            opts['fontweight'] = options['fontweight']
        if 'alpha' in options:
            opts['alpha'] = options['alpha']
        if 'fontstyle' in options:
            opts['fontstyle'] = options['fontstyle']
        if 'bounding_box' in options:
            opts['bbox'] = options['bounding_box']
        if 'zorder' in options:
            opts['zorder'] = options['zorder']
        if options['axis_coords']:
            opts['transform'] = subplot.transAxes
        if 'fontsize' in options:
            val = options['fontsize']
            if isinstance(val, str):
                opts['fontsize'] = val
            else:
                opts['fontsize'] = int(val)
        if 'rotation' in options:
            val = options['rotation']
            if isinstance(val, str):
                opts['rotation'] = options['rotation']
            else:
                opts['rotation'] = float(options['rotation'])

        p = subplot.text(self.x, self.y, self.string, clip_on=options['clip'], **opts)
        if not options['clip']:
            self._bbox_extra_artists=[p]

@rename_keyword(color='rgbcolor')
@options(fontsize=10, rgbcolor=(0,0,1), horizontal_alignment='center',
         vertical_alignment='center', axis_coords=False, clip=False)
def text(string, xy, **options):
    r"""
    Returns a 2D text graphics object at the point `(x,y)`.

    Type ``text.options`` for a dictionary of options for 2D text.

    2D OPTIONS:

    - ``fontsize`` - How big the text is. It is either an integer that specifies
      the size in points or a string which specifies a size (one of
      'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large')

    - ``fontstyle`` - A string either 'normal', 'italique' or 'oblique'.

    - ``fontweight`` - A numeric value in the range 0-1000 or a string (one of
      'ultralight', 'light', 'normal', 'regular', 'book',' 'medium', 'roman',
      'semibold', 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black')

    - ``rgbcolor`` - The color as an RGB tuple

    - ``hue`` - The color given as a hue

    - ``rotation`` - How to rotate the text: angle in degrees, vertical, horizontal

    - ``vertical_alignment`` - How to align vertically: top, center, bottom

    - ``horizontal_alignment`` - How to align horizontally: left, center, right

    - ``axis_coords`` - (default: False) if True, use axis coordinates, so that
      (0,0) is the lower left and (1,1) upper right, regardless of the x and y
      range of plotted values.

    - ``bounding_box`` - a dictionary specifying a bounding box. See the
      examples (or the matplotlib documentation).

    EXAMPLES::

        sage: text("Sage graphics are really neat because they use matplotlib!", (2,12))
        Graphics object consisting of 1 graphics primitive

    Larger font, bold, colored red and transparent text::

        sage: text("I had a dream!", (2,12), alpha=0.3, fontsize='large', fontweight='bold', color='red')
        Graphics object consisting of 1 graphics primitive

    By setting ``horizontal_alignment`` to 'left' the text is guaranteed to be in
    the lower left no matter what::

        sage: text("I got a horse and he lives in a tree", (0,0), axis_coords=True, horizontal_alignment='left')
        Graphics object consisting of 1 graphics primitive

    Various rotations::

        sage: text("noitator", (0,0), rotation=45.0, horizontal_alignment='left', vertical_alignment='bottom')
        Graphics object consisting of 1 graphics primitive

        sage: text("Sage is really neat!!",(0,0), rotation="vertical")
        Graphics object consisting of 1 graphics primitive

    You can also align text differently::

        sage: t1 = text("Hello",(1,1), vertical_alignment="top")
        sage: t2 = text("World", (1,0.5), horizontal_alignment="left")
        sage: t1 + t2   # render the sum
        Graphics object consisting of 2 graphics primitives

    You can save text as part of PDF output::

        sage: text("sage", (0,0), rgbcolor=(0,0,0)).save(os.path.join(SAGE_TMP, 'a.pdf'))

    Some examples of bounding box::

        sage: bbox = {'boxstyle':"rarrow,pad=0.3", 'fc':"cyan", 'ec':"b", 'lw':2}
        sage: text("I feel good", (1,2), bounding_box=bbox)
        Graphics object consisting of 1 graphics primitive

        sage: text("So good", (0,0), bounding_box={'boxstyle':'round', 'fc':'w'})
        Graphics object consisting of 1 graphics primitive

    The possible options of the bounding box are 'boxstyle' (one of 'larrow',
    'rarrow', 'round', 'round4', 'roundtooth', 'sawtooth', 'square'), 'fc' or
    'facecolor', 'ec' or 'edgecolor', 'ha' or horizontalalignment', 'va' or
    'verticalalignment', 'lw' or 'linewidth'.

    A text with a background color::

        sage: text("So good", (-2,2), background_color='red')
        Graphics object consisting of 1 graphics primitive

    Text must be 2D (use the text3d command for 3D text)::

        sage: t = text("hi",(1,2,3))
        Traceback (most recent call last):
        ...
        ValueError: use text3d instead for text in 3d
        sage: t = text3d("hi",(1,2,3))

    Extra options will get passed on to show(), as long as they are valid::

        sage: text("MATH IS AWESOME", (0, 0), fontsize=40, axes=False)
        Graphics object consisting of 1 graphics primitive
        sage: text("MATH IS AWESOME", (0, 0), fontsize=40).show(axes=False) # These are equivalent
    """
    try:
        x, y = xy
    except ValueError:
        if isinstance(xy, (list, tuple)) and len(xy) == 3:
            raise ValueError("use text3d instead for text in 3d")
        raise
    from sage.plot.all import Graphics
    options['rgbcolor'] = to_mpl_color(options['rgbcolor'])
    point = (float(x), float(y))
    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options, ignore='fontsize'))
    g.add_primitive(Text(string, point, options))
    return g
