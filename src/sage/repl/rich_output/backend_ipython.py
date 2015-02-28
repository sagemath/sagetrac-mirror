# -*- encoding: utf-8 -*-
"""
IPython Backend for the Sage Rich Output System

This module defines the IPython backends for
:mod:`sage.repl.rich_output`.
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
from sage.repl.rich_output.backend_base import BackendBase
from sage.repl.rich_output.output_catalog import *


class BackendIPython(BackendBase):
    """
    Common base for the IPython UIs

    EXAMPLES::

        sage: from sage.repl.rich_output.backend_ipython import BackendIPython
        sage: BackendIPython()._repr_()
        Traceback (most recent call last):
        NotImplementedError: derived classes must implement this method
    """

    def install(self, **kwds):
        """
        Switch the Sage rich output to the IPython backend

        INPUT:

        - ``shell`` -- keyword argument. The IPython shell.

        No tests since switching away from the doctest rich output
        backend will break the doctests.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: from sage.repl.rich_output.backend_ipython import BackendIPython
            sage: backend = BackendIPython()
            sage: shell = get_test_shell();
            sage: backend.install(shell=shell)
            sage: shell.run_cell('1+1')
            2
        """
        shell = kwds['shell']
        from sage.repl.display.formatter import SageDisplayFormatter
        shell.display_formatter = SageDisplayFormatter(parent=shell)
        shell.configurables.append(shell.display_formatter)
    
    def set_underscore_variable(self, obj):
        """
        Set the ``_`` builtin variable.
        
        Since IPython handles the history itself, this does nothing.

        INPUT:

        - ``obj`` -- anything.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: from sage.repl.rich_output.backend_ipython import BackendIPython
            sage: backend = BackendIPython()
            sage: backend.set_underscore_variable(123)
            sage: _
            0
        """
        pass

    
class BackendIPythonCommandline(BackendIPython):
    """
    Backend for the IPython Command Line

    EXAMPLES::

        sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
        sage: BackendIPythonCommandline()
        IPython command line
    """

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend._repr_()
            'IPython command line'
        """
        return 'IPython command line'
    
    def supported_output(self):
        """
        Return the outputs that are supported by the IPython commandline backend.

        OUTPUT:

        Iterable of output container classes, that is, subclass of
        :class:`~sage.repl.rich_output.output_basic.OutputBase`).
        The order is ignored.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: supp = backend.supported_output();  supp     # random output
            set([<class 'sage.repl.rich_output.output_graphics.OutputImageGif'>, 
                 ...,
                 <class 'sage.repl.rich_output.output_graphics.OutputImagePng'>])
            sage: from sage.repl.rich_output.output_basic import OutputLatex
            sage: OutputLatex in supp
            True
        """
        return set([
            OutputPlainText, OutputAsciiArt, OutputLatex,
            OutputImagePng, OutputImageGif,
            OutputImagePdf, OutputImageDvi,
            OutputSceneJmol, OutputSceneWavefront,
        ])

    def displayhook(self, plain_text, rich_output):
        """
        Backend implementation of the displayhook
        
        INPUT:

        - ``plain_text`` -- instance of
          :class:`~sage.repl.rich_output.output_basic.OutputPlainText`. The
          plain text version of the output.

        - ``rich_output`` -- instance of an output container class
          (subclass of
          :class:`~sage.repl.rich_output.output_basic.OutputBase`). Guaranteed
          to be one of the output containers returned from
          :meth:`supported_output`, possibly the same as
          ``plain_text``.

        OUTPUT:

        The IPython commandline display hook returns the IPython
        display data, a pair of dictionaries. The first dictionary
        contains mime types as keys and the respective output as
        value. The second dictionary is metadata.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend.displayhook(plain_text, plain_text)
            ({u'text/plain': 'Example plain text output'}, {})
        """
        if isinstance(rich_output, OutputPlainText):
            return ({u'text/plain': rich_output.text.get()}, {})
        elif isinstance(rich_output, OutputAsciiArt):
            return ({u'text/plain': rich_output.ascii_art.get()}, {})
        elif isinstance(rich_output, OutputLatex):
            return ({u'text/plain': rich_output.latex.get()}, {})
        elif isinstance(rich_output, OutputImagePng):
            msg = self.launch_viewer(
                rich_output.png.filename(ext='png'), plain_text.text.get())
            return ({u'text/plain': msg}, {})
        elif isinstance(rich_output, OutputImageGif):
            msg = self.launch_viewer(
                rich_output.gif.filename(ext='gif'), plain_text.text.get())
            return ({u'text/plain': msg}, {})
        elif isinstance(rich_output, OutputImagePdf):
            msg = self.launch_viewer(
                rich_output.pdf.filename(ext='pdf'), plain_text.text.get())
            return ({u'text/plain': msg}, {})
        elif isinstance(rich_output, OutputImageDvi):
            msg = self.launch_viewer(
                rich_output.dvi.filename(ext='dvi'), plain_text.text.get())
            return ({u'text/plain': msg}, {})
        elif isinstance(rich_output, OutputSceneJmol):
            msg = self.launch_jmol(rich_output, plain_text.text.get())
            return ({u'text/plain': msg}, {})
        elif isinstance(rich_output, OutputSceneWavefront):
            msg = self.launch_sage3d(rich_output, plain_text.text.get())
            return ({u'text/plain': msg}, {})
        else:
            raise TypeError('rich_output type not supported')

    def display_immediately(self, plain_text, rich_output):
        """
        Show output without going back to the command line prompt.

        This method is similar to the rich output :meth:`displayhook`,
        except that it can be invoked at any time. On the Sage command
        line it launches viewers just like :meth:`displayhook`.
        
        INPUT:

        Same as :meth:`displayhook`.

        OUTPUT:

        This method does not return anything.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend.display_immediately(plain_text, plain_text)
            Example plain text output
        """
        formatdata, metadata = self.displayhook(plain_text, rich_output)
        print(formatdata[u'text/plain'])

    def launch_viewer(self, image_file, plain_text):
        """
        Launch external viewer for the graphics file.

        INPUT:

        - ``image_file`` -- string. File name of the image file.

        - ``plain_text`` -- string. The plain text representation of
          the image file.

        OUTPUT:

        String. Human-readable message indicating whether the viewer
        was launched successfully.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: backend.launch_viewer('/path/to/foo.bar', 'Graphics object')
            'Launched bar viewer for Graphics object'
        """
        base, dot_ext = os.path.splitext(image_file)
        ext = dot_ext.lstrip(os.path.extsep)
        from sage.misc.viewer import viewer
        command = viewer(ext)
        if not command:
            command = viewer.browser()
        from sage.doctest import DOCTEST_MODE
        if not DOCTEST_MODE:
            os.system('{0} {1} 2>/dev/null 1>/dev/null &'
                      .format(command, image_file))
        return 'Launched {0} viewer for {1}'.format(ext, plain_text)

    def launch_jmol(self, output_jmol, plain_text):
        """
        Launch the stand-alone jmol viewer

        INPUT:

        - ``output_jmol`` --
          :class:`~sage.repl.rich_output.output_graphics3d.OutputSceneJmol`. The
          scene to launch Jmol with.

        - ``plain_text`` -- string. The plain text representation.

        OUTPUT:

        String. Human-readable message indicating that the viewer was launched.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: from sage.repl.rich_output.output_graphics3d import OutputSceneJmol
            sage: backend.launch_jmol(OutputSceneJmol.example(), 'Graphics3d object')
            'Launched jmol viewer for Graphics3d object'
        """
        from sage.doctest import DOCTEST_MODE
        from sage.interfaces.jmoldata import JmolData
        jdata = JmolData()
        if not jdata.is_jvm_available() and not DOCTEST_MODE:
            raise RuntimeError('jmol cannot run, no suitable java version found')
        launch_script = output_jmol.launch_script_filename()
        from sage.env import SAGE_LOCAL
        jmol_cmd = os.path.join(SAGE_LOCAL, 'bin', 'jmol')
        if not DOCTEST_MODE:
            os.system('{0} {1} 2>/dev/null 1>/dev/null &'
                      .format(jmol_cmd, launch_script))
        return 'Launched jmol viewer for {0}'.format(plain_text)

    def launch_sage3d(self, output_wavefront, plain_text):
        """
        Launch the stand-alone java3d viewer

        INPUT:

        - ``output_wavefront`` --
          :class:`~sage.repl.rich_output.output_graphics3d.OutputSceneWavefront`. The
          scene to launch Java3d with.

        - ``plain_text`` -- string. The plain text representation.

        OUTPUT:

        String. Human-readable message indicating that the viewer was launched.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
            sage: backend = BackendIPythonCommandline()
            sage: from sage.repl.rich_output.output_graphics3d import OutputSceneWavefront
            sage: backend.launch_sage3d(OutputSceneWavefront.example(), 'Graphics3d object')
            'Launched Java 3D viewer for Graphics3d object'
        """
        from sage.env import SAGE_LOCAL
        sage3d = os.path.join(SAGE_LOCAL, 'bin', 'sage3d')
        obj = output_wavefront.obj_filename()
        from sage.doctest import DOCTEST_MODE
        if not DOCTEST_MODE:
            os.system('{0} {1} 2>/dev/null 1>/dev/null &'
                      .format(sage3d, obj))
        return 'Launched Java 3D viewer for {0}'.format(plain_text)

    
class BackendIPythonNotebook(BackendIPython):
    """
    Backend for the IPython Notebook

    EXAMPLES::

        sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
        sage: BackendIPythonNotebook()
        IPython notebook
    """

    def _repr_(self):
        """
        Return string representation of the backend

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
            sage: backend = BackendIPythonNotebook()
            sage: backend._repr_()
            'IPython notebook'
        """
        return 'IPython notebook'
    
    def supported_output(self):
        """
        Return the outputs that are supported by the IPython notebook backend.

        OUTPUT:

        Iterable of output container classes, that is, subclass of
        :class:`~sage.repl.rich_output.output_basic.OutputBase`).
        The order is ignored.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
            sage: backend = BackendIPythonNotebook()
            sage: supp = backend.supported_output();  supp     # random output
            set([<class 'sage.repl.rich_output.output_graphics.OutputPlainText'>, 
                 ...,
                 <class 'sage.repl.rich_output.output_graphics.OutputImagePdf'>])
            sage: from sage.repl.rich_output.output_basic import OutputLatex
            sage: OutputLatex in supp
            True

        The IPython notebook cannot display gif images, see
        https://github.com/ipython/ipython/issues/2115 ::

            sage: from sage.repl.rich_output.output_graphics import OutputImageGif
            sage: OutputImageGif in supp
            False
        """
        return set([
            OutputPlainText, OutputAsciiArt, OutputLatex,
            OutputImagePng, OutputImageJpg,
            OutputImageSvg, OutputImagePdf, 
        ])

    def displayhook(self, plain_text, rich_output):
        """
        Backend implementation of the displayhook
        
        INPUT:

        - ``plain_text`` -- instance of
          :class:`~sage.repl.rich_output.output_basic.OutputPlainText`. The
          plain text version of the output.

        - ``rich_output`` -- instance of an output container class
          (subclass of
          :class:`~sage.repl.rich_output.output_basic.OutputBase`). Guaranteed
          to be one of the output containers returned from
          :meth:`supported_output`, possibly the same as
          ``plain_text``.

        OUTPUT:

        The IPython notebook display hook returns the IPython
        display data, a pair of dictionaries. The first dictionary
        contains mime types as keys and the respective output as
        value. The second dictionary is metadata.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
            sage: backend = BackendIPythonNotebook()
            sage: backend.displayhook(plain_text, plain_text)
            ({u'text/plain': 'Example plain text output'}, {})
        """
        if isinstance(rich_output, OutputPlainText):
            return ({u'text/plain': rich_output.text.get()}, {})
        elif isinstance(rich_output, OutputAsciiArt):
            return ({u'text/plain': rich_output.ascii_art.get()}, {})
        elif isinstance(rich_output, OutputLatex):
            return ({u'text/html':  rich_output.mathjax(),
                     u'text/plain': plain_text.text.get(),
            }, {})
        elif isinstance(rich_output, OutputImagePng):
            return ({u'image/png':  rich_output.png.get(),
                     u'text/plain': plain_text.text.get(),
            }, {})
        elif isinstance(rich_output, OutputImageJpg):
            return ({u'image/jpeg':  rich_output.jpg.get(),
                     u'text/plain':  plain_text.text.get(),
            }, {})

        elif isinstance(rich_output, OutputImageSvg):
            return ({u'image/svg+xml': rich_output.svg.get(),
                     u'text/plain':    plain_text.text.get(),
            }, {})
        elif isinstance(rich_output, OutputImagePdf):
            return ({u'image/png':  rich_output.png.get(),
                     u'text/plain': plain_text.text.get(),
            }, {})
        else:
            raise TypeError('rich_output type not supported')

        
    def display_immediately(self, plain_text, rich_output):
        """
        Show output immediately.

        This method is similar to the rich output :meth:`displayhook`,
        except that it can be invoked at any time.

        .. TODO::

            This does not work currently.
        
        INPUT/OUTPUT:

        Same as :meth:`displayhook`.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
            sage: backend = BackendIPythonNotebook()
            sage: backend.display_immediately(plain_text, plain_text)
            ({u'text/plain': 'Example plain text output'}, {})
        """
        return self.displayhook(plain_text, rich_output)
