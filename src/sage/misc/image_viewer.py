r"""
Determination of programs for viewing images and PDF files.

The function :func:`default_image_viewer` defines reasonable defaults for
these programs.  To use something else, use ``image_viewer``.  First
import it::

    sage: from sage.misc.image_viewer import image_viewer

On OS X, PDFs are opened by default using the 'open' command, which
runs whatever has been designated as the PDF viewer in the OS.  To
change this to use 'Adobe Reader'::

    sage: image_viewer.pdf_viewer('open -a /Applications/Adobe\ Reader.app') # not tested

Similarly, you can set ``image_viewer.pdf_viewer(...)``, ``image_viewer.dvi_viewer(...)``,
and ``image_viewer.png_viewer(...)``.  You can make this change permanent by adding
lines like these to your :file:`SAGE_STARTUP_FILE` (which is
:file:`$HOME/.sage/init.sage` by default)::

    from sage.misc.image_viewer import image_viewer
    viewer.pdf_viewer('open -a /Applications/Adobe\ Reader.app')

Functions and classes
---------------------
"""
import os
import webbrowser

from sage.structure.sage_object import SageObject
from sage.features.xdgopen import Open, XdgOpen, Cygstart

VIEWERS = ['dvi', 'pdf', 'png']

def default_image_viewer(viewer=None):
    """
    Set up default programs for opening web pages, PDFs, PNGs, and DVI files.

    INPUT:

    - ``viewer``: ``None`` or a string. Currently ignored.

    EXAMPLES::

        sage: from sage.misc.image_viewer import default_image_viewer
        sage: default_image_viewer() # random -- depends on OS, etc.
        'open'
    """
    if 'SAGE_IMAGE_VIEWER' in os.environ:
        VIEWER = os.environ['SAGE_IMAGE_VIEWER']

    elif XdgOpen().is_present(): # linux
        # On other OS'es try xdg-open if present.
        # See http://portland.freedesktop.org/xdg-utils-1.0.
        VIEWER = 'xdg-open'

    elif Cygstart().is_present():  # Cygwin
        # Windows is easy, since it has a system for
        # determining what opens things.  However, on Cygwin we
        # should access this through the 'cygstart' program rather
        # than trying to run rundll32 directly, which on newer Windows versions
        # has security implications
        VIEWER = 'cygstart'

    elif Open().is_present(): # Darwin
        # Check for "open" last, in case it is present on other
        # systems but doesn't do what we need.
        #
        # Simple on OS X, since there is an open command that opens
        # anything, using the user's preferences.
        VIEWER = 'open'

    else:
        raise NotImplementedError("unable to determine image viewing program for this system; try running sage.misc.image_viewer.image_viewer or setting the environment variable SAGE_IMAGE_VIEWER")

    return VIEWER


# _viewer_prefs: a dictionary holding global preferences for viewers.
_viewer_prefs = {}

class ImageViewer(SageObject):
    """
    Set defaults for various viewing applications: a web browser, a
    dvi viewer, a pdf viewer, and a png viewer.

    EXAMPLES::

        sage: from sage.misc.image_viewer import image_viewer
        sage: old_png = image_viewer.png_viewer()  # indirect doctest
        sage: image_viewer.png_viewer('open -a /Applications/Firefox.app')
        sage: image_viewer.png_viewer()
        'open -a /Applications/Firefox.app'
        sage: image_viewer.png_viewer(old_png) # restore old value
    """
    def _set(self, app=None, TYPE='png'):
        r"""
        Change the default viewer. Return the current setting if the
        argument ``app`` is ``None``.

        INPUT:

        - ``app`` -- ``None`` or a string, the program to use
        - ``TYPE`` -- a string, must be in the list ``VIEWERS`` defined in
          :mod:`sage.misc.image_viewer`.  Default 'png'.

        EXAMPLES::

            sage: from sage.misc.image_viewer import image_viewer
            sage: old_pdf = image_viewer.pdf_viewer()
            sage: image_viewer.pdf_viewer('open -a /Applications/Firefox.app') # indirect doctest
            sage: image_viewer.pdf_viewer()
            'open -a /Applications/Firefox.app'
            sage: image_viewer.pdf_viewer(old_pdf) # restore old value
        """
        TYPE = TYPE.lower()
        if TYPE not in VIEWERS:
            raise ValueError('Unrecognized type of viewer: {}'.format(TYPE))
        if app is None:
            try:
                return _viewer_prefs[TYPE]
            except KeyError:
                return default_image_viewer(TYPE)
        elif app == 'default':
            try:
                del _viewer_prefs[TYPE]
            except KeyError:
                pass
        else:
            _viewer_prefs[TYPE] = app

    def dvi_viewer(self, app=None):
        r"""
        Change the default dvi viewer. Return the current setting if arg
        is ``None``, which is the default.

        INPUT:

        - ``app`` -- ``None`` or a string, the program to use

        EXAMPLES::

            sage: from sage.misc.image_viewer import image_viewer
            sage: old_dvi_app = image_viewer.dvi_viewer()
            sage: image_viewer.dvi_viewer('/usr/bin/xdvi') # indirect doctest
            sage: image_viewer.dvi_viewer()
            '/usr/bin/xdvi'
            sage: image_viewer.dvi_viewer(old_dvi_app) # restore old value
        """
        return self._set(app, TYPE='dvi')

    def pdf_viewer(self, app=None):
        r"""
        Change the default pdf viewer. Return the current setting if arg
        is ``None``, which is the default.

        INPUT:

        - ``app`` -- ``None`` or a string, the program to use

        EXAMPLES::

            sage: from sage.misc.image_viewer import image_viewer
            sage: old_pdf_app = image_viewer.pdf_viewer()
            sage: image_viewer.pdf_viewer('/usr/bin/pdfopen') # indirect doctest
            sage: image_viewer.pdf_viewer()
            '/usr/bin/pdfopen'
            sage: image_viewer.pdf_viewer(old_pdf_app) # restore old value
        """
        return self._set(app, TYPE='pdf')

    def png_viewer(self, app=None):
        r"""
        Change the default png viewer. Return the current setting if arg
        is ``None``, which is the default.

        INPUT:

        - ``app`` -- ``None`` or a string, the program to use

        EXAMPLES::

            sage: from sage.misc.image_viewer import image_viewer
            sage: old_png_app = image_viewer.png_viewer()
            sage: image_viewer.png_viewer('display') # indirect doctest
            sage: image_viewer.png_viewer()
            'display'
            sage: image_viewer.png_viewer(old_png_app) # restore old value
        """
        return self._set(app, TYPE='png')

    def __call__(self, x=None):
        """
        Return the current setting for a viewer program.

        INPUT:

        - ``x`` -- string

        If ``x`` is ``None``, then return the png viewer app. If ``x``
        starts with 'dvi', return the current dvi viewer, and
        similarly if ``x`` starts with 'pdf'.

        EXAMPLES::

            sage: from sage.misc.image_viewer import image_viewer
            sage: image_viewer('pdf') # random -- depends on OS, etc.
            'xdg-open'
            sage: image_viewer('png') == image_viewer()
            True
        """
        if isinstance(x, str):
            x = x.lower()

        if x is None or x.startswith('png'):
            return self.png_viewer()
        elif x.startswith('dvi'):
            return self.dvi_viewer()
        elif x.startswith('pdf'):
            return self.pdf_viewer()

image_viewer = ImageViewer()

def dvi_viewer():
    """
    Return the program used to display a dvi file.  By default, the
    program used depends on the platform and other factors, like
    settings of certain environment variables.  To use a different
    program, call ``image_viewer.dvi_viewer('PROG')``, where 'PROG' is
    the desired program.

    EXAMPLES::

        sage: from sage.misc.image_viewer import dvi_viewer
        sage: dvi_viewer() # random -- depends on OS, etc.
        'open'
    """
    image_viewer()
    return image_viewer.dvi_viewer()

def pdf_viewer():
    """
    Return the program used to display a pdf file.  By default, the
    program used depends on the platform and other factors, like
    settings of certain environment variables.  To use a different
    program, call ``viewer.pdf_viewer('PROG')``, where 'PROG' is the
    desired program.

    EXAMPLES::

        sage: from sage.misc.image_viewer import pdf_viewer, image_viewer
        sage: old_pdf_app = image_viewer.pdf_viewer()
        sage: image_viewer.pdf_viewer('acroread')
        sage: pdf_viewer()
        'acroread'
        sage: image_viewer.pdf_viewer('old_pdf_app')
    """
    image_viewer()
    return image_viewer.pdf_viewer()

def png_viewer():
    """
    Return the program used to display a png file. By default, the
    program used depends on the platform and other factors, like
    settings of certain environment variables.  To use a different
    program, call ``image_viewer.png_viewer('PROG')``, where 'PROG' is
    the desired program.

    EXAMPLES::

        sage: from sage.misc.image_viewer import png_viewer
        sage: png_viewer() # random -- depends on OS, etc.
        'xdg-open'
    """
    image_viewer()
    return image_viewer.png_viewer()
