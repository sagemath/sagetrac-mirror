"""
Arcs of circles and ellipses
"""
#*****************************************************************************
#       Copyright (C) 2010 Vincent Delecroix <20100.delecroix@gmail.com>,
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
from sage.plot.colors import to_mpl_color

from sage.plot.misc import options, rename_keyword

from math import fmod, sin, cos, pi, atan

class Arc(GraphicPrimitive):
    """
    Primitive class for the Arc graphics type.  See ``arc?`` for information
    about actually plotting an arc of a circle or an ellipse.

    INPUT:

    - ``x,y`` - coordinates of the center of the arc

    - ``r1``, ``r2`` - lengths of the two radii

    - ``angle`` - angle of the horizontal with width

    - ``sector`` - sector of angle

    - ``options`` - dict of valid plot options to pass to constructor

    EXAMPLES:

    Note that the construction should be done using ``arc``::

        sage: from sage.plot.arc import Arc
        sage: print Arc(0,0,1,1,pi/4,pi/4,pi/2,{})
        Arc with center (0.0,0.0) radii (1.0,1.0) angle 0.785398163397 inside the sector (0.785398163397,1.57079632679)
    """
    def __init__(self, x, y, r1, r2, angle, s1, s2, options):
        """
        Initializes base class ``Arc``.

        EXAMPLES:

            sage: A = arc((2,3),1,1,pi/4,(0,pi))
            sage: A[0].x == 2
            True
            sage: A[0].y == 3
            True
            sage: A[0].r1 == 1
            True
            sage: A[0].r2 == 1
            True
            sage: bool(A[0].angle == pi/4)
            True
            sage: bool(A[0].s1 == 0)
            True
            sage: bool(A[0].s2 == pi)
            True

        TESTS::

            sage: from sage.plot.arc import Arc
            sage: a = Arc(0,0,1,1,0,0,1,{})
            sage: print loads(dumps(a))
            Arc with center (0.0,0.0) radii (1.0,1.0) angle 0.0 inside the sector (0.0,1.0)
        """
        self.x = float(x)
        self.y = float(y)
        self.r1 = float(r1)
        self.r2 = float(r2)
        if self.r1 <= 0 or self.r2 <= 0:
            raise ValueError("the radii must be positive real numbers.")

        self.angle = float(angle)
        self.s1 = float(s1)
        self.s2 = float(s2)
        if self.s2 < self.s1:
            self.s1,self.s2=self.s2,self.s1
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        The bounding box is computed as minimal as possible.

        EXAMPLES:

        An example without angle::

            sage: p = arc((-2, 3), 1, 2)
            sage: d = p.get_minmax_data()
            sage: d['xmin']
            -3.0
            sage: d['xmax']
            -1.0
            sage: d['ymin']
            1.0
            sage: d['ymax']
            5.0

        The same example with a rotation of angle `\pi/2`::

            sage: p = arc((-2, 3), 1, 2, pi/2)
            sage: d = p.get_minmax_data()
            sage: d['xmin']
            -4.0
            sage: d['xmax']
            0.0
            sage: d['ymin']
            2.0
            sage: d['ymax']
            4.0
        """
        from sage.plot.plot import minmax_data

        twopi = 2*pi

        s1 = self.s1
        s2 = self.s2
        s = s2-s1
        s1 = fmod(s1,twopi)
        if s1 < 0: s1 += twopi
        s2 = fmod(s1 + s,twopi)
        if s2 < 0: s2 += twopi

        r1 = self.r1
        r2 = self.r2

        angle = fmod(self.angle,twopi)
        if angle < 0: angle += twopi

        epsilon = float(0.0000001)

        cos_angle = cos(angle)
        sin_angle = sin(angle)

        if cos_angle > 1-epsilon:
            xmin=-r1; ymin=-r2
            xmax=r1; ymax=r2
            axmin = pi; axmax = 0
            aymin = 3*pi/2; aymax = pi/2

        elif cos_angle < -1+epsilon:
            xmin=-r1; ymin=-r2
            xmax=r1; ymax=r2
            axmin=0; axmax=pi
            aymin=pi/2; aymax=3*pi/2

        elif sin_angle > 1-epsilon:
            xmin=-r2; ymin=-r1
            xmax=r2; ymax=r1
            axmin = pi/2; axmax = 3*pi/2
            aymin = pi; aymax = 0

        elif sin_angle < -1+epsilon:
            xmin=-r2; ymin=-r1
            xmax=r2; ymax=r1
            axmin = 3*pi/2; axmax = pi/2
            aymin = 0; aymax = pi

        else:
            tan_angle = sin_angle / cos_angle
            axmax = atan(-r2/r1*tan_angle)
            if axmax < 0: axmax += twopi
            xmax = (
              r1 * cos_angle * cos(axmax) -
              r2 * sin_angle * sin(axmax))
            if xmax < 0:
                xmax = -xmax
                axmax = fmod(axmax+pi,twopi)
            xmin = -xmax
            axmin = fmod(axmax + pi,twopi)

            aymax = atan(r2/(r1*tan_angle))
            if aymax < 0: aymax += twopi
            ymax = (
              r1 * sin_angle * cos(aymax) +
              r2 * cos_angle * sin(aymax))
            if ymax < 0:
                ymax = -ymax
                aymax = fmod(aymax+pi,twopi)
            ymin = -ymax
            aymin = fmod(aymax + pi, twopi)

        if s < twopi-epsilon: # bb determined by the sector
            def is_cyclic_ordered(x1,x2,x3):
                return (
                  (x1 < x2 and x2 < x3) or
                  (x2 < x3 and x3 < x1) or
                  (x3 < x1 and x1 < x2))

            x1 = cos_angle*r1*cos(s1) - sin_angle*r2*sin(s1)
            x2 = cos_angle*r1*cos(s2) - sin_angle*r2*sin(s2)
            y1 = sin_angle*r1*cos(s1) + cos_angle*r2*sin(s1)
            y2 = sin_angle*r1*cos(s2) + cos_angle*r2*sin(s2)

            if is_cyclic_ordered(s1,s2,axmin): xmin = min(x1,x2)
            if is_cyclic_ordered(s1,s2,aymin): ymin = min(y1,y2)
            if is_cyclic_ordered(s1,s2,axmax): xmax = max(x1,x2)
            if is_cyclic_ordered(s1,s2,aymax): ymax = max(y1,y2)

        return minmax_data([self.x + xmin, self.x + xmax],
                           [self.y + ymin, self.y + ymax],
                           dict=True)

    def _allowed_options(self):
        """
        Return the allowed options for the ``Arc`` class.

        EXAMPLES::

            sage: p = arc((3, 3), 1, 1)
            sage: p[0]._allowed_options()['alpha']
            'How transparent the figure is.'
        """
        return {'alpha':'How transparent the figure is.',
                'thickness':'How thick the border of the arc is.',
                'hue':'The color given as a hue.',
                'rgbcolor':'The color',
                'zorder':'2D only: The layer level in which to draw',
                'linestyle':"2D only: The style of the line, which is one of "
                "'dashed', 'dotted', 'solid', 'dashdot', or '--', ':', '-', '-.', "
                "respectively."}

    def _repr_(self):
        """
        String representation of ``Arc`` primitive.

        EXAMPLES::

            sage: from sage.plot.arc import Arc
            sage: print Arc(2,3,2.2,2.2,0,2,3,{})
            Arc with center (2.0,3.0) radii (2.2,2.2) angle 0.0 inside the sector (2.0,3.0)
        """
        return "Arc with center (%s,%s) radii (%s,%s) angle %s inside the sector (%s,%s)" %(self.x,self.y,self.r1,self.r2,self.angle,self.s1,self.s2)

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: A = arc((1,1),3,4,pi/4,(pi,4*pi/3)); A
        """
        import matplotlib.patches as patches
        from sage.plot.misc import get_matplotlib_linestyle


        options = self.options()

        p = patches.Arc(
                (self.x,self.y),
                2.*self.r1,
                2.*self.r2,
                fmod(self.angle,2*pi)*(180./pi),
                self.s1*(180./pi),
                self.s2*(180./pi))
        p.set_linewidth(float(options['thickness']))
        a = float(options['alpha'])
        p.set_alpha(a)
        z = int(options.pop('zorder',1))
        p.set_zorder(z)
        c = to_mpl_color(options['rgbcolor'])
        p.set_linestyle(get_matplotlib_linestyle(options['linestyle'],return_type='long'))
        p.set_edgecolor(c)
        subplot.add_patch(p)

    def plot3d(self):
        r"""
        TESTS::

            sage: from sage.plot.arc import Arc
            sage: Arc(0,0,1,1,0,0,1,{}).plot3d()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

@rename_keyword(color='rgbcolor')
@options(alpha=1, thickness=1, linestyle='solid', zorder=5,rgbcolor='blue',
         aspect_ratio=1.0)
def arc(center, r1, r2=None, angle=0.0, sector=(0.0,2*pi), **options):
    r"""
    An arc (that is a portion of a circle or an ellipse)

    Type ``arc.options`` to see all options.

    INPUT:

    - ``center`` - 2-tuple of real numbers - position of the center.

    - ``r1``, ``r2`` - positive real numbers - radii of the ellipse. If only ``r1``
      is set, then the two radii are supposed to be equal and this function returns
      an arc of of circle.

    - ``angle`` - real number - angle between the horizontal and the axis that
      corresponds to ``r1``.

    - ``sector`` - 2-tuple (default: (0,2*pi))- angles sector in which the arc will
      be drawn.

    OPTIONS:

    - ``alpha`` - float (default: 1) - transparency

    - ``thickness`` - float (default: 1) - thickness of the arc

    - ``color``, ``rgbcolor`` - string or 2-tuple (default: 'blue') - the color
      of the arc

    - ``linestyle`` - string (default: ``'solid'``) - The style of the line,
      which is one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``,
      or ``'--'``, ``':'``, ``'-'``, ``'-.'``, respectively.

    EXAMPLES:

    Plot an arc of circle centered at (0,0) with radius 1 in the sector
    `(\pi/4,3*\pi/4)`::

        sage: arc((0,0), 1, sector=(pi/4,3*pi/4))

    Plot an arc of an ellipse between the angles 0 and `\pi/2`::

        sage: arc((2,3), 2, 1, sector=(0,pi/2))

    Plot an arc of a rotated ellipse between the angles 0 and `\pi/2`::

        sage: arc((2,3), 2, 1, angle=pi/5, sector=(0,pi/2))

    Plot an arc of an ellipse in red with a dashed linestyle::

        sage: arc((0,0), 2, 1, 0, (0,pi/2), linestyle="dashed", color="red")
        sage: arc((0,0), 2, 1, 0, (0,pi/2), linestyle="--", color="red")

    The default aspect ratio for arcs is 1.0::

        sage: arc((0,0), 1, sector=(pi/4,3*pi/4)).aspect_ratio()
        1.0

    It is not possible to draw arcs in 3D::

        sage: A = arc((0,0,0), 1)
        Traceback (most recent call last):
        ...
        NotImplementedError
    """
    from sage.plot.all import Graphics

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = options.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale == 'semilogy' or scale == 'semilogx':
        options['aspect_ratio'] = 'automatic'

    if len(center)==2:
        if r2 is None: r2 = r1
        g = Graphics()
        g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
        if len(sector) != 2:
            raise ValueError("the sector must consist of two angles")
        g.add_primitive(Arc(
            center[0],center[1],
            r1,r2,
            angle,
            sector[0],sector[1],
            options))
        return g
    elif len(center)==3:
        raise NotImplementedError

