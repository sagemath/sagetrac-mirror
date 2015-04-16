r"""
Animated plots

Animations are generated from a list (or other iterable) of graphics
objects.
Images are produced by calling the ``save_image`` method on each input
object, creating a sequence of PNG files.
These are then assembled to various target formats using different
tools.
In particular, the ``convert`` program from ImageMagick_ can be used to
generate an animated GIF file.
FFmpeg_ (with the command line program ``ffmpeg``) provides support for
various video formats, but also an alternative method of generating
animated GIFs.
For `browsers which support it`_, APNG_ can be used as another
alternative which works without any extra dependencies.

.. Warning::

    Note that ImageMagick and FFmpeg are not included with Sage, and
    must be installed by the user.  On unix systems, type ``which
    convert`` at a command prompt to see if ``convert`` (part of the
    ImageMagick suite) is installed.  If it is, you will be given its
    location.  Similarly, you can check for ``ffmpeg`` with ``which
    ffmpeg``.  See the websites of ImageMagick_ or FFmpeg_ for
    installation instructions.

EXAMPLES:

The sine function::

    sage: sines = [plot(c*sin(x), (-2*pi,2*pi), color=Color(c,0,0), ymin=-1, ymax=1) for c in sxrange(0,1,.2)]
    sage: a = animate(sines)
    sage: a         # optional -- ImageMagick
    Animation with 5 frames
    sage: a.show()  # optional -- ImageMagick

Animate using FFmpeg_ instead of ImageMagick::

    sage: f = tmp_filename(ext='.gif')
    sage: a.save(filename=f, use_ffmpeg=True) # optional -- ffmpeg

Animate as an APNG_::

    sage: a.apng()  # long time
    doctest:...: DeprecationWarning: use tmp_filename instead
    See http://trac.sagemath.org/17234 for details.

An animated :class:`sage.plot.graphics.GraphicsArray` of rotating ellipses::

    sage: E = animate((graphics_array([[ellipse((0,0),a,b,angle=t,xmin=-3,xmax=3)+circle((0,0),3,color='blue') for a in range(1,3)] for b in range(2,4)]) for t in sxrange(0,pi/4,.15)))
    sage: str(E)    # animations produced from a generator do not have a known length
    'Animation with unknown number of frames'
    sage: E.show()  # optional -- ImageMagick

A simple animation of a circle shooting up to the right::

    sage: c = animate([circle((i,i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
    ....:               xmin=0,ymin=0,xmax=2,ymax=2,figsize=[2,2])
    sage: c.show() # optional -- ImageMagick


Animations of 3d objects::

    sage: var('s,t')
    (s, t)
    sage: def sphere_and_plane(x):
    ....:     return sphere((0,0,0),1,color='red',opacity=.5)+parametric_plot3d([t,x,s],(s,-1,1),(t,-1,1),color='green',opacity=.7)
    sage: sp = animate([sphere_and_plane(x) for x in sxrange(-1,1,.3)])
    sage: sp[0]      # first frame
    Graphics3d Object
    sage: sp[-1]     # last frame
    Graphics3d Object
    sage: sp.show()  # optional -- ImageMagick

    sage: (x,y,z) = var('x,y,z')
    sage: def frame(t):
    ....:     return implicit_plot3d((x^2 + y^2 + z^2), (x, -2, 2), (y, -2, 2), (z, -2, 2), plot_points=60, contour=[1,3,5], region=lambda x,y,z: x<=t or y>=t or z<=t)
    sage: a = animate([frame(t) for t in srange(.01,1.5,.2)])
    sage: a[0]       # long time
    Graphics3d Object
    sage: a.show()   # optional -- ImageMagick

If the input objects do not have a ``save_image`` method, then the
animation object attempts to make an image by calling its internal
method :meth:`sage.plot.animate.Animation.make_image`.  This is
illustrated by the following example::

    sage: t = var('t')
    sage: a = animate((sin(c*pi*t) for c in sxrange(1,2,.2)))
    sage: a.show()  # optional -- ImageMagick


AUTHORS:

- William Stein
- John Palmieri
- Niles Johnson (2013-12): Expand to animate more graphics objects
- Martin von Gagern (2014-12): Added APNG support

.. REFERENCES (not rendered as a section, but linked inline):

.. _ImageMagick: http://www.imagemagick.org
.. _FFmpeg: http://www.ffmpeg.org
.. _APNG: https://wiki.mozilla.org/APNG_Specification
.. _`browsers which support it`: http://caniuse.com/#feat=apng

"""

############################################################################
#  Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
############################################################################

import os
import struct
import zlib

from sage.structure.sage_object import SageObject
from sage.misc.temporary_file import tmp_dir, tmp_filename, graphics_filename
import plot
import sage.misc.misc
import sage.misc.viewer


def animate(frames, **kwds):
    r"""
    Animate a list of frames by creating a
    :class:`sage.plot.animate.Animation` object.

    EXAMPLES::

        sage: t = var('t')
        sage: a = animate((cos(c*pi*t) for c in sxrange(1,2,.2)))
        sage: a.show()  # optional -- ImageMagick

    See also :mod:`sage.plot.animate` for more examples.
    """
    return Animation(frames, **kwds)

class Animation(SageObject):
    r"""
    Return an animation of a sequence of plots of objects.

    INPUT:


    - ``v`` - iterable of Sage objects. These should preferably be
      graphics objects, but if they aren't then :meth:`make_image` is
      called on them.

    - ``xmin, xmax, ymin, ymax`` - the ranges of the x and y axes.

    - ``**kwds`` - all additional inputs are passed onto the rendering
      command. E.g., use figsize to adjust the resolution and aspect
      ratio.


    EXAMPLES::

        sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.3)],
        ....:                xmin=0, xmax=2*pi, figsize=[2,1])
        sage: a                 # optional -- ImageMagick
        Animation with 21 frames
        sage: a[:5]             # optional -- ImageMagick
        Animation with 5 frames
        sage: a.show()          # optional -- ImageMagick
        sage: a[:5].show()      # optional -- ImageMagick

    The :meth:`show` method takes arguments to specify the
    delay between frames (measured in hundredths of a second, default
    value 20) and the number of iterations (default value 0, which
    means to iterate forever). To iterate 4 times with half a second
    between each frame::

        sage: a.show(delay=50, iterations=4) # optional -- ImageMagick

    An animation of drawing a parabola::

        sage: step = 0.1
        sage: L = Graphics()
        sage: v = []
        sage: for i in srange(0,1,step):
        ....:       L += line([(i,i^2),(i+step,(i+step)^2)], rgbcolor=(1,0,0), thickness=2)
        ....:       v.append(L)
        sage: a = animate(v, xmin=0, ymin=0)
        sage: a.show() # optional -- ImageMagick
        sage: show(L)

    TESTS:

    This illustrates that :trac:`2066` is fixed (setting axes
    ranges when an endpoint is 0)::

        sage: animate([plot(sin, -1,1)], xmin=0, ymin=0)._kwds['xmin']
        0

    We check that :trac:`7981` is fixed::

        sage: a = animate([plot(sin(x + float(k)), (0, 2*pi), ymin=-5, ymax=5)
        ....:              for k in srange(0,2*pi,0.3)])
        sage: a.show() # optional -- ImageMagick

    Do not convert input iterator to a list, but ensure that
    the frame count is known after rendering the frames::

        sage: a = animate((plot(x^p, (x,0,2)) for p in sxrange(1,2,.1)))
        sage: str(a)
        'Animation with unknown number of frames'
        sage: a.png()    # long time
        '.../'
        sage: len(a)     # long time
        10
        sage: a._frames
        <generator object ...

    """
    def __init__(self, v=None, **kwds):
        r"""
        Return an animation of a sequence of plots of objects.  See
        documentation of :func:`animate` for more details and
        examples.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.3)],
            ....:                xmin=0, xmax=2*pi, figsize=[2,1]) # indirect doctest
            sage: a           # optional -- ImageMagick
            Animation with 21 frames
        """
        self._frames = v
        self._kwds = kwds

    def _combine_kwds(self, *kwds_tuple):
        """
        Returns a dictionary which is a combination of the all the
        dictionaries in kwds_tuple. This also does the appropriate thing
        for taking the mins and maxes of all of the x/y mins/maxes.

        EXAMPLES::

            sage: a = animate([plot(sin, -1,1)], xmin=0, ymin=0)
            sage: kwds1 = {'a':1, 'b':2, 'xmin':2, 'xmax':5}
            sage: kwds2 = {'b':3, 'xmin':0, 'xmax':4}
            sage: kwds = a._combine_kwds(kwds1, kwds2)
            sage: list(sorted(kwds.items()))
            [('a', 1), ('b', 3), ('xmax', 5), ('xmin', 0)]

        Test that the bug reported in :trac:`12107` has been fixed::

            sage: kwds3 = {}
            sage: kwds4 = {'b':3, 'xmin':0, 'xmax':4}
            sage: a._combine_kwds(kwds3, kwds4)['xmin']
            0
        """
        new_kwds = {}

        for kwds in kwds_tuple:
            new_kwds.update(kwds)

        import __builtin__
        for name in ['xmin', 'xmax', 'ymin', 'ymax']:
            values = [v for v in [kwds.get(name, None) for kwds in kwds_tuple] if v is not None]
            if values:
                new_kwds[name] = getattr(__builtin__, name[1:])(values)
        return new_kwds


    def __getitem__(self, i):
        """
        Get a frame from an animation or
        slice this animation returning a subanimation.

        EXAMPLES::

            sage: a = animate([circle((i,-i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
            ....:               xmin=0,ymin=-2,xmax=2,ymax=0,figsize=[2,2])
            sage: a           # optional -- ImageMagick
            Animation with 10 frames
            sage: frame2 = a[2]  # indirect doctest
            sage: frame2.show()
            sage: a.show() # optional -- ImageMagick
            sage: a[3:7]   # optional -- ImageMagick   # indirect doctest
            Animation with 4 frames
            sage: a[3:7].show() # optional -- ImageMagick
        """
        if isinstance(i, slice):
            return Animation(self._frames[i], **self._kwds)
        else:
            return self._frames[i]

    def _repr_(self):
        """
        Print representation for an animation.

        EXAMPLES::

            sage: a = animate([circle((i,-i), 1-1/(i+1), hue=i/10) for i in srange(0,2,0.2)],
            ....:               xmin=0,ymin=-2,xmax=2,ymax=0,figsize=[2,2])
            sage: a           # optional -- ImageMagick
            Animation with 10 frames
            sage: a._repr_()
            'Animation with 10 frames'
        """
        try:
            num = len(self)
        except TypeError:
            num = "unknown number of"
        return "Animation with %s frames"%num

    def __add__(self, other):
        """
        Add two animations. This has the effect of superimposing the two
        animations frame-by-frame.

        EXAMPLES::

            sage: a = animate([circle((i,0),1) for i in srange(0,2,0.4)],
            ....:             xmin=0, ymin=-1, xmax=3, ymax=1,
            ....:             figsize=[3,1], ticks=[1,1])
            sage: a.show()        # optional -- ImageMagick
            sage: b = animate([circle((0,i),1,hue=0) for i in srange(0,2,0.4)],
            ....:             xmin=0, ymin=-1, xmax=2, ymax=3,
            ....:             figsize=[1,2], ticks=[1,1])
            sage: b.show()        # optional -- ImageMagick
            sage: s = a+b         # indirect doctest
            sage: len(a), len(b)
            (5, 5)
            sage: len(s)
            5
            sage: s.show()        # optional -- ImageMagick
        """
        if not isinstance(other, Animation):
            other = Animation(other)

        kwds = self._combine_kwds(self._kwds, other._kwds)

        #Combine the frames
        m = max(len(self), len(other))
        frames = [a+b for a,b in zip(self._frames, other._frames)]
        frames += self._frames[m:] + other._frames[m:]

        return Animation(frames, **kwds)

    def __mul__(self, other):
        """
        Multiply two animations. This has the effect of appending the two
        animations (the second comes after the first).

        EXAMPLES::

            sage: a = animate([circle((i,0),1,thickness=20*i) for i in srange(0,2,0.4)],
            ....:                xmin=0, ymin=-1, xmax=3, ymax=1, figsize=[2,1], axes=False)
            sage: a.show()             # optional -- ImageMagick
            sage: b = animate([circle((0,i),1,hue=0,thickness=20*i) for i in srange(0,2,0.4)],
            ....:                xmin=0, ymin=-1, xmax=1, ymax=3, figsize=[1,2], axes=False)
            sage: b.show()             # optional -- ImageMagick
            sage: p = a*b              # indirect doctest
            sage: len(a), len(b)
            (5, 5)
            sage: len(p)
            10
            sage: (a*b).show()         # optional -- ImageMagick
        """
        if not isinstance(other, Animation):
            other = Animation(other)

        kwds = self._combine_kwds(self._kwds, other._kwds)

        return Animation(self._frames + other._frames, **kwds)

    def __len__(self):
        """
        Length of self

        EXAMPLES::

            sage: a = animate([circle((i,0),1,thickness=20*i) for i in srange(0,2,0.4)],
            ....:                xmin=0, ymin=-1, xmax=3, ymax=1, figsize=[2,1], axes=False)
            sage: len(a)
            5
        """
        try:
            return self._num_frames
        except AttributeError:
            return len(self._frames)

    def make_image(self, frame, filename, **kwds):
        r"""
        Given a frame which has no ``save_image()`` method, make a graphics
        object and save it as an image with the given filename.  By default, this is
        :meth:`sage.plot.plot.plot`.  To make animations of other objects,
        override this method in a subclass.

        EXAMPLES::

            sage: from sage.plot.animate import Animation
            sage: class MyAnimation(Animation):
            ....:    def make_image(self, frame, filename, **kwds):
            ....:        P = parametric_plot(frame[0], frame[1], **frame[2])
            ....:        P.save_image(filename,**kwds)

            sage: t = var('t')
            sage: x = lambda t: cos(t)
            sage: y = lambda n,t: sin(t)/n
            sage: B = MyAnimation([([x(t), y(i+1,t)],(t,0,1), {'color':Color((1,0,i/4)), 'aspect_ratio':1, 'ymax':1}) for i in range(4)])

            sage: d = B.png(); v = os.listdir(d); v.sort(); v  # long time
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
            sage: B.show()  # not tested

            sage: class MyAnimation(Animation):
            ....:    def make_image(self, frame, filename, **kwds):
            ....:        G = frame.plot()
            ....:        G.set_axes_range(floor(G.xmin()),ceil(G.xmax()),floor(G.ymin()),ceil(G.ymax()))
            ....:        G.save_image(filename, **kwds)

            sage: B = MyAnimation([graphs.CompleteGraph(n) for n in range(7,11)], figsize=5)
            sage: d = B.png()
            sage: v = os.listdir(d); v.sort(); v
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
            sage: B.show()  # not tested

        """
        p = plot.plot(frame, **kwds)
        p.save_image(filename)

    def png(self, dir=None):
        r"""
        Render PNG images of the frames in this animation, saving them
        in ``dir``.  Return the absolute path to that directory.  If
        the frames have been previously rendered and ``dir`` is
        ``None``, just return the directory in which they are stored.

        When ``dir`` is other than ``None``, force re-rendering of
        frames.

        INPUT:

        - ``dir`` -- Directory in which to store frames.  Default
          ``None``; in this case, a temporary directory will be
          created for storing the frames.

        EXAMPLES::

            sage: a = animate([plot(x^2 + n) for n in range(4)], ymin=0, ymax=4)
            sage: d = a.png(); v = os.listdir(d); v.sort(); v  # long time
            ['00000000.png', '00000001.png', '00000002.png', '00000003.png']
        """
        if dir is None:
            try:
                return self._png_dir
            except AttributeError:
                pass
            dir = tmp_dir()
        i = 0
        for frame in self._frames:
            filename = '%s/%08d.png'%(dir,i)
            try:
                save_image = frame.save_image
            except AttributeError:
                self.make_image(frame, filename, **self._kwds)
            else:
                save_image(filename, **self._kwds)
            i += 1
        self._num_frames = i
        self._png_dir = dir
        return dir

    def graphics_array(self, ncols=3):
        r"""
        Return a :class:`sage.plot.graphics.GraphicsArray` with plots of the
        frames of this animation, using the given number of columns.
        The frames must be acceptable inputs for
        :class:`sage.plot.graphics.GraphicsArray`.


        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: v = [E.change_ring(GF(p)).plot(pointsize=30) for p in [97, 101, 103, 107]]
            sage: a = animate(v, xmin=0, ymin=0, axes=False)
            sage: a        # optional -- ImageMagick
            Animation with 4 frames
            sage: a.show() # optional -- ImageMagick

        Modify the default arrangement of array::

            sage: g = a.graphics_array(); print g
            Graphics Array of size 2 x 3
            sage: g.show(figsize=[6,3]) # optional

        Specify different arrangement of array and save it with a given file name::

            sage: g = a.graphics_array(ncols=2); print g
            Graphics Array of size 2 x 2
            sage: f = tmp_filename(ext='.png')
            sage: g.save(f) # optional

        Frames can be specified as a generator too; it is internally converted to a list::

            sage: t = var('t')
            sage: b = animate((plot(sin(c*pi*t)) for c in sxrange(1,2,.2)))
            sage: g = b.graphics_array(); print g
            Graphics Array of size 2 x 3
            sage: g.show() # optional
        """
        ncols = int(ncols)
        frame_list = list(self._frames)
        n = len(frame_list)
        nrows, rem = divmod(n,ncols)
        if rem > 0:
            nrows += 1
        return plot.graphics_array(frame_list, nrows,  ncols)

    def gif(self, delay=20, savefile=None, iterations=0, show_path=False,
            use_ffmpeg=False):
        r"""
        Returns an animated gif composed from rendering the graphics
        objects in self.

        This method will only work if either (a) the ImageMagick
        software suite is installed, i.e., you have the ``convert``
        command or (b) ``ffmpeg`` is installed.  See the web sites of
        ImageMagick_ and FFmpeg_ for more details.  By default, this
        produces the gif using ``convert`` if it is present.  If this
        can't find ``convert`` or if ``use_ffmpeg`` is True, then it
        uses ``ffmpeg`` instead.

        INPUT:

        -  ``delay`` - (default: 20) delay in hundredths of a
           second between frames

        -  ``savefile`` - file that the animated gif gets saved
           to

        -  ``iterations`` - integer (default: 0); number of
           iterations of animation. If 0, loop forever.

        -  ``show_path`` - boolean (default: False); if True,
           print the path to the saved file

        - ``use_ffmpeg`` - boolean (default: False); if True, use
          'ffmpeg' by default instead of 'convert'.

        If ``savefile`` is not specified: in notebook mode, display the
        animation; otherwise, save it to a default file name.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ....:             xmin=0, xmax=2*pi, ymin=-1, ymax=1, figsize=[2,1])
            sage: td = tmp_dir()
            sage: a.gif()              # not tested
            sage: a.gif(savefile=td + 'my_animation.gif', delay=35, iterations=3)  # optional -- ImageMagick
            sage: with open(td + 'my_animation.gif', 'rb') as f: print '\x21\xf9\x04\x08\x23\x00' in f.read()  # optional -- ImageMagick
            True
            sage: a.gif(savefile=td + 'my_animation.gif', show_path=True) # optional -- ImageMagick
            Animation saved to .../my_animation.gif.
            sage: a.gif(savefile=td + 'my_animation_2.gif', show_path=True, use_ffmpeg=True) # optional -- ffmpeg
            Animation saved to .../my_animation_2.gif.

        .. note::

           If neither ffmpeg nor ImageMagick is installed, you will
           get an error message like this::

              Error: Neither ImageMagick nor ffmpeg appears to be installed. Saving an
              animation to a GIF file or displaying an animation requires one of these
              packages, so please install one of them and try again.

              See www.imagemagick.org and www.ffmpeg.org for more information.

        .. REFERENCES (not rendered as a section, but linked inline):

        .. _ImageMagick: http://www.imagemagick.org
        .. _FFmpeg: http://www.ffmpeg.org
        """
        from sage.misc.sage_ostools import have_program
        have_convert = have_program('convert')
        have_ffmpeg = self._have_ffmpeg()
        if use_ffmpeg or not have_convert:
            if have_ffmpeg:
                self.ffmpeg(savefile=savefile, show_path=show_path,
                            output_format='.gif', delay=delay,
                            iterations=iterations)
            else:
                if not have_convert:
                    msg = """
Error: Neither ImageMagick nor ffmpeg appears to be installed. Saving an
animation to a GIF file or displaying an animation requires one of these
packages, so please install one of them and try again.

See www.imagemagick.org and www.ffmpeg.org for more information."""
                else:
                    msg = """
Error: ffmpeg does not appear to be installed.  Download it from
www.ffmpeg.org, or use 'convert' to produce gifs instead."""
                raise OSError(msg)
        else:
            if not savefile:
                savefile = graphics_filename(ext='.gif')
            if not savefile.endswith('.gif'):
                savefile += '.gif'
            savefile = os.path.abspath(savefile)
            d = self.png()
            cmd = ( 'cd "%s"; sage-native-execute convert -dispose Background '
                    '-delay %s -loop %s *.png "%s"' ) % ( d, int(delay),
                        int(iterations), savefile )
            from subprocess import check_call, CalledProcessError
            try:
                check_call(cmd, shell=True)
                if show_path:
                    print "Animation saved to file %s." % savefile
            except (CalledProcessError, OSError):
                msg = """
Error: Cannot generate GIF animation.  Verify that convert
(ImageMagick) or ffmpeg is installed, and that the objects passed to
the animate command can be saved in PNG image format.

See www.imagemagick.org and www.ffmpeg.org for more information."""
                raise OSError(msg)

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: a = animate([plot(x^2 + n) for n in range(4)], ymin=0, ymax=4)
            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: a._rich_repr_(dm)       # optional -- ImageMagick
            OutputImageGif container
        """
        OutputImageGif = display_manager.types.OutputImageGif
        if OutputImageGif not in display_manager.supported_output():
            return
        return display_manager.graphics_from_save(
            self.save, kwds, '.gif', OutputImageGif)

    def show(self, delay=20, iterations=0):
        r"""
        Show this animation immediately.

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        INPUT:

        -  ``delay`` -- (default: 20) delay in hundredths of a
           second between frames

        -  ``iterations`` -- integer (default: 0); number of
           iterations of animation. If 0, loop forever.

        OUTPUT:

        This method does not return anything. Use :meth:`save` if you
        want to save the figure as an image.

        .. note::

           Currently this is done using an animated gif, though this
           could change in the future. This requires that either
           ffmpeg or the ImageMagick suite (in particular, the
           ``convert`` command) is installed.

        See also the :meth:`ffmpeg` method.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ....:                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: a.show()       # optional -- ImageMagick

        The preceding will loop the animation forever. If you want to show
        only three iterations instead::

            sage: a.show(iterations=3)    # optional -- ImageMagick

        To put a half-second delay between frames::

            sage: a.show(delay=50)        # optional -- ImageMagick

        .. note::

           If you don't have ffmpeg or ImageMagick installed, you will
           get an error message like this::

              Error: Neither ImageMagick nor ffmpeg appears to be installed. Saving an
              animation to a GIF file or displaying an animation requires one of these
              packages, so please install one of them and try again.

              See www.imagemagick.org and www.ffmpeg.org for more information.
        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        dm.display_immediately(self, delay=delay, iterations=iterations)

    def _have_ffmpeg(self):
        """
        Return True if the program 'ffmpeg' is installed.  See
        www.ffmpeg.org to download ffmpeg.

        EXAMPLES::

            sage: a = animate([plot(sin, -1,1)], xmin=0, ymin=0)
            sage: a._have_ffmpeg() # random: depends on whether ffmpeg is installed
            False
        """
        from sage.misc.sage_ostools import have_program
        return have_program('ffmpeg')

    def ffmpeg(self, savefile=None, show_path=False, output_format=None,
               ffmpeg_options='', delay=None, iterations=0, pix_fmt='rgb24'):
        r"""
        Returns a movie showing an animation composed from rendering
        the frames in self.

        This method will only work if ffmpeg is installed.  See
        http://www.ffmpeg.org for information about ffmpeg.

        INPUT:

        -  ``savefile`` - file that the mpeg gets saved to.

        .. warning:

            This will overwrite ``savefile`` if it already exists.

        - ``show_path`` - boolean (default: False); if True, print the
          path to the saved file

        - ``output_format`` - string (default: None); format and
          suffix to use for the video.  This may be 'mpg', 'mpeg',
          'avi', 'gif', or any other format that ffmpeg can handle.
          If this is None and the user specifies ``savefile`` with a
          suffix, say ``savefile='animation.avi'``, try to determine the
          format ('avi' in this case) from that file name.  If no file
          is specified or if the suffix cannot be determined, 'mpg' is
          used.

        - ``ffmpeg_options`` - string (default: ''); this string is
          passed directly to ffmpeg.

        - ``delay`` - integer (default: None); delay in hundredths of a
          second between frames.  The framerate is 100/delay.
          This is not supported for mpeg files: for mpegs, the frame
          rate is always 25 fps.

        - ``iterations`` - integer (default: 0); number of iterations
          of animation. If 0, loop forever.  This is only supported
          for animated gif output and requires ffmpeg version 0.9 or
          later.  For older versions, set ``iterations=None``.

        - ``pix_fmt`` - string (default: 'rgb24'); used only for gif
          output.  Different values such as 'rgb8' or 'pal8' may be
          necessary depending on how ffmpeg was installed.  Set
          ``pix_fmt=None`` to disable this option.

        If ``savefile`` is not specified: in notebook mode, display
        the animation; otherwise, save it to a default file name.  Use
        :func:`sage.misc.misc.set_verbose` with ``level=1`` to see
        additional output.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ....:             xmin=0, xmax=2*pi, ymin=-1, ymax=1, figsize=[2,1])
            sage: td = tmp_dir()
            sage: a.ffmpeg(savefile=td + 'new.mpg')       # optional -- ffmpeg
            sage: a.ffmpeg(savefile=td + 'new.avi')       # optional -- ffmpeg
            sage: a.ffmpeg(savefile=td + 'new.gif')       # optional -- ffmpeg
            sage: a.ffmpeg(savefile=td + 'new.mpg', show_path=True) # optional -- ffmpeg
            Animation saved to .../new.mpg.

        .. note::

           If ffmpeg is not installed, you will get an error message
           like this::

              Error: ffmpeg does not appear to be installed. Saving an animation to
              a movie file in any format other than GIF requires this software, so
              please install it and try again.

              See www.ffmpeg.org for more information.


        TESTS::

            sage: a.ffmpeg(output_format='gif',delay=30,iterations=5)     # optional -- ffmpeg
            doctest:...: DeprecationWarning: use tmp_filename instead
            See http://trac.sagemath.org/17234 for details.
        """
        if not self._have_ffmpeg():
            msg = """Error: ffmpeg does not appear to be installed. Saving an animation to
a movie file in any format other than GIF requires this software, so
please install it and try again."""
            raise OSError(msg)
        else:
            if savefile is None:
                if output_format is None:
                    output_format = '.mpg'
                else:
                    if output_format[0] != '.':
                        output_format = '.'+output_format
                savefile = graphics_filename(ext=output_format)
            else:
                if output_format is None:
                    suffix = os.path.splitext(savefile)[1]
                    if len(suffix) > 0:
                        output_format = suffix
                    else:
                        output_format = '.mpg'
            if not savefile.endswith(output_format):
                savefile += output_format
            early_options = ''
            if output_format == '.gif':
                # We try to set reasonable options for gif output.
                #
                # Older versions of ffmpeg (before 0.9, summer 2011)
                # use the option -loop_output instead of -loop.
                # Setting iterations=None is a way of preventing sage
                # from adding the -loop option.  A separate
                # -loop_output option can be added with the
                # ffmpeg_options argument.
                if iterations is not None:
                    loop_cmd = '-loop {0} '.format(iterations)
                else:
                    loop_cmd = ''
                # A pix_fmt value is required for some but not all
                # ffmpeg installations.  Setting pix_fmt=None will
                # prevent sage from adding this option, and it may be
                # controlled separately through ffmpeg_options.
                if pix_fmt is not None:
                    pix_fmt_cmd = '-pix_fmt {0} '.format(pix_fmt)
                else:
                    pix_fmt_cmd = ''
                ffmpeg_options += ' {0}{1}'.format(pix_fmt_cmd,loop_cmd)
            if delay is not None and output_format != '.mpeg' and output_format != '.mpg':
                early_options += ' -r %s ' % int(100/delay)
            savefile = os.path.abspath(savefile)
            pngdir = self.png()
            pngs = os.path.join(pngdir, "%08d.png")
            # For ffmpeg, it seems that some options, like '-g ... -r
            # ...', need to come before the input file names, while
            # some options, like '-pix_fmt rgb24', need to come
            # afterwards.  Hence 'early_options' and 'ffmpeg_options'
            cmd = 'cd "%s"; sage-native-execute ffmpeg -y -f image2 %s -i %s %s %s' % (pngdir, early_options, pngs, ffmpeg_options, savefile)
            from subprocess import check_call, CalledProcessError, PIPE
            try:
                if sage.misc.misc.get_verbose() > 0:
                    set_stderr = None
                else:
                    set_stderr = PIPE
                sage.misc.misc.verbose("Executing '%s'" % cmd,level=1)
                sage.misc.misc.verbose("\n---- ffmpeg output below ----\n")
                check_call(cmd, shell=True, stderr=set_stderr)
                if show_path:
                    print "Animation saved to file %s." % savefile
            except (CalledProcessError, OSError):
                print "Error running ffmpeg."
                raise

    def apng(self, savefile=None, show_path=False, delay=20, iterations=0):
        r"""
        Creates an animated PNG composed from rendering the graphics
        objects in self. Return the absolute path to that file.

        Notice that not all web browsers are capable of displaying APNG
        files, though they should still present the first frame of the
        animation as a fallback.

        The generated file is not optimized, so it may be quite large.

        Input:

        -  ``delay`` - (default: 20) delay in hundredths of a
           second between frames

        -  ``savefile`` - file that the animated gif gets saved
           to

        -  ``iterations`` - integer (default: 0); number of
           iterations of animation. If 0, loop forever.

        -  ``show_path`` - boolean (default: False); if True,
           print the path to the saved file

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ....:                xmin=0, xmax=2*pi, figsize=[2,1])
            sage: dir = tmp_dir()
            sage: a.apng()  # long time
            sage: a.apng(savefile=dir + 'my_animation.png', delay=35, iterations=3)  # long time
            sage: a.apng(savefile=dir + 'my_animation.png', show_path=True)  # long time
            Animation saved to .../my_animation.png.

        If the individual frames have different sizes, an error will be raised::

            sage: a = animate([plot(sin(x), (x, 0, k)) for k in range(1,4)],
            ....:             ymin=-1, ymax=1, aspect_ratio=1, figsize=[2,1])
            sage: a.apng()  # long time
            Traceback (most recent call last):
            ...
            ValueError: Chunk IHDR mismatch

        """
        pngdir = self.png()
        if savefile is None:
            savefile = graphics_filename('.png')
        with open(savefile, "wb") as out:
            apng = APngAssembler(
                out, len(self),
                delay=delay, num_plays=iterations)
            for i in range(len(self)):
                png = os.path.join(pngdir, "%08d.png" % i)
                apng.add_frame(png)
        if show_path:
            print "Animation saved to file %s." % savefile

    def save(self, filename=None, show_path=False, use_ffmpeg=False, **kwds):
        r"""
        Save this animation.

        INPUT:

        -  ``filename`` - (default: None) name of save file

        -  ``show_path`` - boolean (default: False); if True,
           print the path to the saved file

        - ``use_ffmpeg`` - boolean (default: False); if True, use
          'ffmpeg' by default instead of 'convert' when creating GIF
          files.

        If filename is None, then in notebook mode, display the
        animation; otherwise, save the animation to a GIF file.  If
        filename ends in '.sobj', save to an sobj file.  Otherwise,
        try to determine the format from the filename extension
        ('.mpg', '.gif', '.avi', etc.).  If the format cannot be
        determined, default to GIF.

        For GIF files, either ffmpeg or the ImageMagick suite must be
        installed.  For other movie formats, ffmpeg must be installed.
        An sobj file can be saved with no extra software installed.

        EXAMPLES::

            sage: a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
            ....:             xmin=0, xmax=2*pi, ymin=-1, ymax=1, figsize=[2,1])
            sage: td = tmp_dir()
            sage: a.save()         # not tested
            sage: a.save(td + 'wave.gif')   # optional -- ImageMagick
            sage: a.save(td + 'wave.gif', show_path=True)   # optional -- ImageMagick
            Animation saved to file .../wave.gif.
            sage: a.save(td + 'wave.avi', show_path=True)   # optional -- ffmpeg
            Animation saved to file .../wave.avi.
            sage: a.save(td + 'wave0.sobj')
            sage: a.save(td + 'wave1.sobj', show_path=True)
            Animation saved to file .../wave1.sobj.

        TESTS:

        Ensure that we can pass delay and iteration count to the saved
        GIF image (see :trac:`18176`)::

            sage: a.save(td + 'wave.gif')   # optional -- ImageMagick
            sage: with open(td + 'wave.gif', 'rb') as f: print '!\xf9\x04\x08\x14\x00' in f.read()  # optional -- ImageMagick
            True
            sage: with open(td + 'wave.gif', 'rb') as f: print '!\xff\x0bNETSCAPE2.0\x03\x01\x00\x00\x00' in f.read()  # optional -- ImageMagick
            True
            sage: a.save(td + 'wave.gif', delay=35)   # optional -- ImageMagick
            sage: with open(td + 'wave.gif', 'rb') as f: print '!\xf9\x04\x08\x14\x00' in f.read()  # optional -- ImageMagick
            False
            sage: with open(td + 'wave.gif', 'rb') as f: print '!\xf9\x04\x08\x23\x00' in f.read()  # optional -- ImageMagick
            True
            sage: a.save(td + 'wave.gif', iterations=3)   # optional -- ImageMagick
            sage: with open(td + 'wave.gif', 'rb') as f: print '!\xff\x0bNETSCAPE2.0\x03\x01\x00\x00\x00' in f.read()  # optional -- ImageMagick
            False
            sage: with open(td + 'wave.gif', 'rb') as f: print '!\xff\x0bNETSCAPE2.0\x03\x01\x03\x00\x00' in f.read()  # optional -- ImageMagick
            True
        """
        if filename is None:
            suffix = '.gif'
        else:
            suffix = os.path.splitext(filename)[1]
            if len(suffix) == 0:
                suffix = '.gif'

        if filename is None or suffix == '.gif':
            self.gif(savefile=filename, show_path=show_path,
                     use_ffmpeg=use_ffmpeg, **kwds)
        elif suffix == '.sobj':
            SageObject.save(self, filename)
            if show_path:
                print "Animation saved to file %s." % filename
        else:
            self.ffmpeg(savefile=filename, show_path=show_path, **kwds)


class APngAssembler(object):
    r"""
    Builds an APNG_ (Animated PNG) from a sequence of PNG files.
    This is used by the :meth:`sage.plot.animate.Animation.apng` method.

    This code is quite simple; it does little more than copying chunks
    from input PNG files to the output file. There is no optimization
    involved. This does not depend on external programs or libraries.

    INPUT:

        - ``out`` -- a file opened for binary writing to which the data
          will be written

        - ``num_frames`` -- the number of frames in the animation

        - ``num_plays`` -- how often to iterate, 0 means infinitely

        - ``delay`` -- numerator of the delay fraction in seconds

        - ``delay_denominator`` -- denominator of the delay in seconds

    EXAMPLES::

        sage: from sage.plot.animate import APngAssembler
        sage: def assembleAPNG():
        ....:     a = animate([sin(x + float(k)) for k in srange(0,2*pi,0.7)],
        ....:                 xmin=0, xmax=2*pi, figsize=[2,1])
        ....:     pngdir = a.png()
        ....:     outfile = sage.misc.temporary_file.tmp_filename(ext='.png')
        ....:     with open(outfile, "wb") as f:
        ....:         apng = APngAssembler(f, len(a))
        ....:         for i in range(len(a)):
        ....:             png = os.path.join(pngdir, "{:08d}.png".format(i))
        ....:             apng.add_frame(png, delay=10*i + 10)
        ....:     return outfile
        ....:
        sage: assembleAPNG()  # long time
        '...png'

    .. REFERENCES:

    .. _APNG: https://wiki.mozilla.org/APNG_Specification
    """

    magic = b"\x89PNG\x0d\x0a\x1a\x0a"
    mustmatch = frozenset([b"IHDR", b"PLTE", b"bKGD", b"cHRM", b"gAMA",
                           b"pHYs", b"sBIT", b"tRNS"])

    def __init__(self, out, num_frames,
                 num_plays=0, delay=200, delay_denominator=100):
        r"""
        Initialize for creation of an APNG file.
        """
        self._last_seqno = -1
        self._idx = 0
        self._first = True
        self.out = out
        self.num_frames = num_frames
        self.num_plays = num_plays
        self.default_delay_numerator = delay
        self.default_delay_denominator = delay_denominator
        self._matchref = dict()
        self.out.write(self.magic)

    def add_frame(self, pngfile, delay=None, delay_denominator=None):
        r"""
        Adds a single frame to the APNG file.

        INPUT:

        - ``pngfile`` -- file name of the PNG file with data for this frame

        - ``delay`` -- numerator of the delay fraction in seconds

        - ``delay_denominator`` -- denominator of the delay in seconds

        If the delay is not specified, the default from the constructor
        applies.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: from StringIO import StringIO
            sage: buf = StringIO()
            sage: apng = APngAssembler(buf, 2)
            sage: fn = APngAssembler._testData("input1", True)
            sage: apng.add_frame(fn, delay=0x567, delay_denominator=0x1234)
            sage: fn = APngAssembler._testData("input2", True)
            sage: apng.add_frame(fn)
            sage: len(buf.getvalue())
            217
            sage: buf.getvalue() == APngAssembler._testData("anim12", False)
            True
            sage: apng.add_frame(fn)
            Traceback (most recent call last):
            ...
            RuntimeError: Already reached the declared number of frames

        """
        if self._idx == self.num_frames:
            raise RuntimeError("Already reached the declared number of frames")
        self.delay_numerator = self.default_delay_numerator
        self.delay_denominator = self.default_delay_denominator
        self._actl_written = False
        self._fctl_written = False
        if delay is not None:
            self.delay_numerator = delay
        if delay_denominator is not None:
            self.delay_denominator = delay_denominator
        self._add_png(pngfile)
        self._idx += 1
        if self._idx == self.num_frames:
            self._chunk(b"IEND", b"")

    def set_default(self, pngfile):
        r"""
        Adds a default image for the APNG file.

        This image is used as a fallback in case some application does
        not understand the APNG format.  This method must be called
        prior to any calls to the ``add_frame`` method, if it is called
        at all.  If it is not called, then the first frame of the
        animation will be the default.

        INPUT:

        - ``pngfile`` -- file name of the PNG file with data
          for the default image

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: from StringIO import StringIO
            sage: buf = StringIO()
            sage: apng = APngAssembler(buf, 1)
            sage: fn = APngAssembler._testData("input1", True)
            sage: apng.set_default(fn)
            sage: fn = APngAssembler._testData("input2", True)
            sage: apng.add_frame(fn, delay=0x567, delay_denominator=0x1234)
            sage: len(buf.getvalue())
            179
            sage: buf.getvalue() == APngAssembler._testData("still1anim2", False)
            True
            sage: apng.add_frame(fn)
            Traceback (most recent call last):
            ...
            RuntimeError: Already reached the declared number of frames

        """
        if self._idx != 0:
            raise RuntimeError("Default image must precede all animation frames")
        self._actl_written = False
        self._fctl_written = True
        self._add_png(pngfile)

    def _add_png(self, pngfile):
        r"""
        Add data from one PNG still image.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_add_png", reads=False)
            enter _add_png('...png')
              write _current_chunk = ('\x00\x00\x00\r', 'IHDR', '\x00\x00\x00\x03\x00\x00\x00\x02\x08\x00\x00\x00\x00', '\xb8\x1f9\xc6')
              call _copy() -> None
              call _first_IHDR('\x00\x00\x00\x03\x00\x00\x00\x02\x08\x00\x00\x00\x00') -> None
              write _current_chunk = ('\x00\x00\x00\x04', 'gAMA', '\x00\x01\x86\xa0', '1\xe8\x96_')
              call _copy() -> None
              write _current_chunk = ('\x00\x00\x00\x07', 'tIME', '\x07\xde\x06\x1b\x0b&$', '\x1f0z\xd5')
              write _current_chunk = ('\x00\x00\x00\x08', 'IDAT', 'img1data', '\xce\x8aI\x99')
              call _first_IDAT('img1data') -> None
              write _current_chunk = ('\x00\x00\x00\x00', 'IEND', '', '\xaeB`\x82')
              write _first = False
            exit _add_png -> None
            enter _add_png('...png')
              write _current_chunk = ('\x00\x00\x00\r', 'IHDR', '\x00\x00\x00\x03\x00\x00\x00\x02\x08\x00\x00\x00\x00', '\xb8\x1f9\xc6')
              write _current_chunk = ('\x00\x00\x00\x04', 'gAMA', '\x00\x01\x86\xa0', '1\xe8\x96_')
              write _current_chunk = ('\x00\x00\x00\x04', 'IDAT', 'img2', '\x0ei\xab\x1d')
              call _next_IDAT('img2') -> None
              write _current_chunk = ('\x00\x00\x00\x04', 'IDAT', 'data', 'f\x94\xcbx')
              call _next_IDAT('data') -> None
              write _current_chunk = ('\x00\x00\x00\x00', 'IEND', '', '\xaeB`\x82')
              write _first = False
            exit _add_png -> None
        """
        with open(pngfile, 'rb') as png:
            if png.read(8) != self.magic:
                raise ValueError("{} is not a PNG file".format(pngfile))
            while True:
                chead = png.read(8)
                if len(chead) == 0:
                    break
                clen, ctype = struct.unpack(">L4s", chead)
                cdata = png.read(clen)
                ccrc = png.read(4)
                utype = ctype.decode("ascii")
                self._current_chunk = (chead[:4], ctype, cdata, ccrc)
                if ctype in self.mustmatch:
                    ref = self._matchref.get(ctype)
                    if ref is None:
                        self._matchref[ctype] = cdata
                        self._copy()
                    else:
                        if cdata != ref:
                            raise ValueError("Chunk {} mismatch".format(utype))
                met = ("_first_" if self._first else "_next_") + utype
                try:
                    met = getattr(self, met)
                except AttributeError:
                    pass
                else:
                    met(cdata)
        self._first = False

    def _seqno(self):
        r"""
        Generate next sequence number.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: from StringIO import StringIO
            sage: buf = StringIO()
            sage: apng = APngAssembler(buf, 1)
            sage: apng._seqno()
            '\x00\x00\x00\x00'
            sage: apng._seqno()
            '\x00\x00\x00\x01'
            sage: apng._seqno()
            '\x00\x00\x00\x02'
        """
        self._last_seqno += 1
        return struct.pack(">L", self._last_seqno)

    def _first_IHDR(self, data):
        r"""
        Remember image size.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_first_IHDR")
            enter _first_IHDR('\x00\x00\x00\x03\x00\x00\x00\x02\x08\x00\x00\x00\x00')
              write width = 3
              write height = 2
            exit _first_IHDR -> None
        """
        w, h, d, ctype, comp, filt, ilace = struct.unpack(">2L5B", data)
        self.width = w
        self.height = h

    def _first_IDAT(self, data):
        r"""
        Write acTL and fcTL, then copy as IDAT.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_first_IDAT")
            enter _first_IDAT('img1data')
              call _actl() -> None
              call _fctl() -> None
              call _copy() -> None
            exit _first_IDAT -> None
        """
        self._actl()
        self._fctl()
        self._copy()

    def _next_IDAT(self, data):
        r"""
        Write fcTL, then convert to fdAT.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_next_IDAT")
            enter _next_IDAT('img2')
              call _fctl() -> None
              call _seqno() -> '\x00\x00\x00\x02'
              call _chunk('fdAT', '\x00\x00\x00\x02img2') -> None
            exit _next_IDAT -> None
            enter _next_IDAT('data')
              call _fctl() -> None
              call _seqno() -> '\x00\x00\x00\x03'
              call _chunk('fdAT', '\x00\x00\x00\x03data') -> None
            exit _next_IDAT -> None
        """
        self._fctl()
        maxlen = 0x7ffffffb
        while len(data) > maxlen:
            self._chunk(b"fdAT", self._seqno() + data[:maxlen])
            data = data[maxlen:]
        self._chunk(b"fdAT", self._seqno() + data)

    def _copy(self):
        r"""
        Copy an existing chunk without modification.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_copy")
            enter _copy()
              read _current_chunk = ('\x00\x00\x00\r', 'IHDR', '\x00\x00\x00\x03\x00\x00\x00\x02\x08\x00\x00\x00\x00', '\xb8\x1f9\xc6')
              read out = <StringIO.StringIO instance at ...
              read out = <StringIO.StringIO instance at ...
              read out = <StringIO.StringIO instance at ...
              read out = <StringIO.StringIO instance at ...
            exit _copy -> None
            enter _copy()
              read _current_chunk = ('\x00\x00\x00\x04', 'gAMA', '\x00\x01\x86\xa0', '1\xe8\x96_')
            ...
              read _current_chunk = ('\x00\x00\x00\x08', 'IDAT', 'img1data', '\xce\x8aI\x99')
            ...
            exit _copy -> None
        """
        for d in self._current_chunk:
            self.out.write(d)

    def _actl(self):
        r"""
        Write animation control data (acTL).

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_actl")
            enter _actl()
              read _actl_written = False
              read num_frames = 2
              read num_plays = 0
              call _chunk('acTL', '\x00\x00\x00\x02\x00\x00\x00\x00') -> None
              write _actl_written = True
            exit _actl -> None
        """
        if self._actl_written:
            return
        data = struct.pack(">2L", self.num_frames, self.num_plays)
        self._chunk(b"acTL", data)
        self._actl_written = True

    def _fctl(self):
        r"""
        Write frame control data (fcTL).

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1("_fctl")
            enter _fctl()
              read _fctl_written = False
              read width = 3
              read height = 2
              read delay_numerator = 1383
              read delay_denominator = 4660
              call _seqno() -> '\x00\x00\x00\x00'
              call _chunk('fcTL', '\x00\x00\x00\x00\x00\x00\x00\x03\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x00\x05g\x124\x01\x00') -> None
              write _fctl_written = True
            exit _fctl -> None
            enter _fctl()
              read _fctl_written = False
              read width = 3
              read height = 2
              read delay_numerator = 200
              read delay_denominator = 100
              call _seqno() -> '\x00\x00\x00\x01'
              call _chunk('fcTL', '\x00\x00\x00\x01\x00\x00\x00\x03\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x00\x00\xc8\x00d\x01\x00') -> None
              write _fctl_written = True
            exit _fctl -> None
            enter _fctl()
              read _fctl_written = True
            exit _fctl -> None
        """
        if self._fctl_written:
            return
        data = struct.pack(
            ">4L2H2B",
            self.width, self.height, 0, 0,
            self.delay_numerator, self.delay_denominator,
            1, 0)
        self._chunk(b"fcTL", self._seqno() + data)
        self._fctl_written = True

    def _chunk(self, ctype, cdata):
        r"""
        Write a new (or modified) chunk of data

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: from StringIO import StringIO
            sage: buf = StringIO()
            sage: apng = APngAssembler(buf, 1)
            sage: buf.getvalue()
            '\x89PNG\r\n\x1a\n'
            sage: apng._chunk("abcd", "efgh")
            sage: buf.getvalue()
            '\x89PNG\r\n\x1a\n\x00\x00\x00\x04abcdefgh\xae\xef*P'
        """
        ccrc = struct.pack(">L", zlib.crc32(ctype + cdata) & 0xffffffff)
        clen = struct.pack(">L", len(cdata))
        for d in [clen, ctype, cdata, ccrc]:
            self.out.write(d)

    @classmethod
    def _hex2bin(cls, h):
        r"""
        Convert hex data to binary.

        This is a helper method used for testing.
        Most data is given as lower-case hex digits,
        possibly intermixed with whitespace.
        A dot causes the next four bytes to be copied verbatim
        even if they look like hex digits. This is used for chunk types.
        Other characters which are not hex digits are passed verbatim.

        EXAMPLE::

            sage: from sage.plot.animate import APngAssembler
            sage: h2b = APngAssembler._hex2bin
            sage: h2b("0123") == b"\x01\x23"
            True
            sage: h2b(" 01 \n 23 ") == b"\x01\x23"
            True
            sage: h2b(".abcdef") == b"abcd\xef"
            True
            sage: h2b("PNG") == b"PNG"
            True
        """
        b = []
        while h:
            if h[0] in ' \n': # ignore whitespace
                h = h[1:]
            elif h[0] in '0123456789abcdef': # hex byte
                b.append(int(h[:2], 16))
                h = h[2:]
            elif h[0] == '.': # for chunk type
                b.extend(ord(h[i]) for i in range(1,5))
                h = h[5:]
            else: # for PNG magic
                b.append(ord(h[0]))
                h = h[1:]
        return ''.join(map(chr,b))

    @classmethod
    def _testData(cls, name, asFile):
        r"""
        Retrieve data for test cases.

        INPUT:

        - ``name``: The name of the file content.

        - ``asFile``: Whether to return a binary string of the named data
                      or the path of a file containing that data.

        EXAMPLE::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testData("input1", False)
            '\x89PNG\r\n\x1a\n\x00...'
            sage: APngAssembler._testData("input2", True)
            '...png'
        """
        data = {

            # Input 1: one PNG image, except the data makes no real sense
            "input1": """89 PNG 0d0a1a0a
            0000000d.IHDR 00000003000000020800000000 b81f39c6
            00000004.gAMA 000186a0 31e8965f
            00000007.tIME 07de061b0b2624 1f307ad5
            00000008.IDAT 696d673164617461 ce8a4999
            00000000.IEND ae426082""",

            # Input 2: slightly different, data in two chunks
            "input2": """89 PNG 0d0a1a0a
            0000000d.IHDR 00000003000000020800000000 b81f39c6
            00000004.gAMA 000186a0 31e8965f
            00000004.IDAT 696d6732 0e69ab1d
            00000004.IDAT 64617461 6694cb78
            00000000.IEND ae426082""",

            # Expected output 1: both images as frames of an animation
            "anim12": """89 PNG 0d0a1a0a
            0000000d.IHDR 00000003000000020800000000 b81f39c6
            00000004.gAMA 000186a0 31e8965f
            00000008.acTL 0000000200000000 f38d9370
            0000001a.fcTL 000000000000000300000002
                          0000000000000000056712340100 b4f729c9
            00000008.IDAT 696d673164617461 ce8a4999
            0000001a.fcTL 000000010000000300000002
                          000000000000000000c800640100 1b92eb4d
            00000008.fdAT 00000002696d6732 9cfb89a3
            00000008.fdAT 0000000364617461 c966c076
            00000000.IEND ae426082""",

            # Expected output 2: first image as fallback, second as animation
            "still1anim2": """89 PNG 0d0a1a0a
            0000000d.IHDR 00000003000000020800000000 b81f39c6
            00000004.gAMA 000186a0 31e8965f
            00000008.acTL 0000000100000000 b42de9a0
            00000008.IDAT 696d673164617461 ce8a4999
            0000001a.fcTL 000000000000000300000002
                          0000000000000000056712340100 b4f729c9
            00000008.fdAT 00000001696d6732 db5bf373
            00000008.fdAT 0000000264617461 f406e9c6
            00000000.IEND ae426082""",

        }
        d = cls._hex2bin(data[name])
        if asFile:
            from sage.misc.temporary_file import tmp_filename
            fn = tmp_filename(ext=".png")
            with open(fn, 'wb') as f:
                f.write(d)
            return fn
        return d

    @classmethod
    def _testCase1(cls, methodToTrace=None, **kwds):
        r"""
        Run common test case.

        This test case is one animation of two frames.
        The named method (if not None) will be traced during execution.
        This will demonstrate the role of each method in the doctests.

        TESTS::

            sage: from sage.plot.animate import APngAssembler
            sage: APngAssembler._testCase1()
        """
        from sage.doctest.fixtures import trace_method
        from StringIO import StringIO
        buf = StringIO()
        apng = cls(buf, 2)
        if methodToTrace is not None:
            trace_method(apng, methodToTrace, **kwds)
        apng.add_frame(cls._testData("input1", True),
                       delay=0x567, delay_denominator=0x1234)
        apng.add_frame(cls._testData("input2", True))
        out = buf.getvalue()
        assert len(out) == 217
        assert out == cls._testData("anim12", False)
