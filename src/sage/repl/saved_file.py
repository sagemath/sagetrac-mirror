# -*- encoding: utf-8 -*-
r"""
Saved File Proxy

The purpose of :class:`SavedFile` is somewhat cryptic, it is the
return value of `save()` methods. If the save command is the last
statement in a repl command or notebook, its return value is seen by
the display hook. We use this to automatically download the saved file
from a remote server.

EXAMPLES::

    sage: filename = tmp_filename(ext='filename.png')
    sage: saved_file = plot(sin(x), x, 0, 10).save(filename)
    sage: type(saved_file)
    <class 'sage.repl.saved_file.SavedFileWithShow'>
    sage: saved_file
    Saved file at /...filename.png
    sage: saved_file.load()
    '\x89PNG...'
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
import imghdr
from sage.structure.sage_object import SageObject


class RichReprHolder(SageObject):

    def __init__(self, plain_text, rich_output):
        """
        Proxy object for Rich Output

        This object holds a rich output container instance, its only
        purpose is to return that as rich output when it is pretty
        printed.

        INPUT

        Instance of :class:`~sage.repl.rich_output.output_base.OutputBase`

        EXAMPLES::
        
            sage: from sage.repl.saved_file import RichReprHolder
            sage: from sage.repl.rich_output.output_catalog import OutputAsciiArt
            sage: ascii_art = OutputAsciiArt.example()
            sage: holder = RichReprHolder('ASCII Art', ascii_art)
            sage: holder     # displays the stored ascii art
            [                        *   *   *    * ]
            [      **   **   *    *  *   *  *    *  ]
            [ ***, * , *  , **, ** , *, * , * , *   ]
        """
        if self._can_display(rich_output):
            self._plain_text = plain_text
            self._rich_output = rich_output
        else:
            self._plain_text = 'Your user interfaces cannot display {0}' \
                .format(self._plain_text)
            self._rich_output = None

    def _can_display(self, rich_output):
        """
        Test whether the user interface can display ``rich_output``

        INPUT
        
        - ``rich_output`` --
          :class:`~sage.repl.rich_output.output_base.OutputBase`. The
          rich output to (potentially) display.

        OUTPUT:

        Boolean.

        EXAMPLES::
        
            sage: from sage.repl.saved_file import RichReprHolder
            sage: from sage.repl.rich_output.output_catalog import OutputAsciiArt
            sage: ascii_art = OutputAsciiArt.example()
            sage: holder = RichReprHolder('ASCII Art', ascii_art)
            sage: holder._can_display(ascii_art)
            True
        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        return any(isinstance(rich_output, output_container)
                   for output_container in dm.supported_output())
        
    def _repr_(self):
        """
        Return the plain text representation
        
        OUTPUT:

        String.

            sage: from sage.repl.saved_file import RichReprHolder
            sage: from sage.repl.rich_output.output_catalog import OutputAsciiArt
            sage: ascii_art = OutputAsciiArt.example()
            sage: holder = RichReprHolder('ASCII Art', ascii_art)
            sage: holder._repr_()
            'ASCII Art'
            sage: str(holder)
            'ASCII Art'
            sage: repr(holder)
            'ASCII Art'
        """
        return self._plain_text
        
    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::
        
            sage: from sage.repl.saved_file import RichReprHolder
            sage: from sage.repl.rich_output.output_catalog import OutputAsciiArt
            sage: ascii_art = OutputAsciiArt.example()
            sage: holder = RichReprHolder('ASCII Art', ascii_art)
            sage: dm = sage.repl.rich_output.get_display_manager()
            sage: holder._rich_repr_(dm) is ascii_art
            True
        """
        if any(isinstance(self._rich_output, output_container)
               for output_container in display_manager.supported_output()):
            return self._rich_output
        from sage.repl.rich_output.output_basic import OutputPlainText
        return OutputPlainText(self._plain_text)


class SavedFile(SageObject):

    def __init__(self, filename, rich_output=None):
        """
        Saved File Proxy

        This object holds data about where we just saved a file.

        INPUT:

        - ``filename`` -- string. The just-saved filename.

        - ``rich_output`` -- a rich output container type, see
          :mod:`~sage.repl.rich_output.output_catalog` (optional,
          default: ``None``). If set, the :meth:`show` method might be
          able to display the just-saved file, depending on the
          capabilities of the user interface. If not set, some effort
          is made to identify the file type but no promises.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImagePng
            sage: png_file = OutputImagePng.example().png.filename(ext='filename.png')
            sage: from sage.repl.saved_file import SavedFile
            sage: SavedFile(png_file)
            Saved file at /...filename.png
        """
        self._filename = filename
        if rich_output is not None:
            self._rich_output = rich_output
        else:
            from sage.repl.rich_output.guess_output import guess_output
            self._rich_output = guess_output(filename)
        if self._rich_output is not None:
            self.__class__ = SavedFileWithShow    
                
    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.saved_file import SavedFile
            sage: from sage.repl.rich_output.output_catalog import OutputImagePng
            sage: png_file = OutputImagePng.example().png.filename(ext='.png')
            sage: SavedFile(png_file)._repr_()
            'Saved file at /....png'
        """
        return 'Saved file at {0}'.format(self.abspath())

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::
        
            sage: from sage.repl.saved_file import SavedFile
            sage: from sage.repl.rich_output.output_catalog import OutputImagePng
            sage: png_file = OutputImagePng.example().png.filename(ext='.png')
            sage: display_manager = sage.repl.rich_output.get_display_manager()
            sage: SavedFile(png_file)._rich_repr_(display_manager)
            OutputSavedFile container
            sage: SavedFile(png_file)   # indirect doctest        
            Saved file at /....png
        """
        from sage.repl.rich_output.output_special import OutputSavedFile
        if OutputSavedFile in display_manager.supported_output():
            return OutputSavedFile(self.abspath())

    def basename(self):
        """
        Return the file name excluding path

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputImagePng
            sage: png_file = OutputImagePng.example().png.filename(ext='.png')
            sage: from sage.repl.saved_file import SavedFile
            sage: '/' in SavedFile(png_file).basename()
            False
        """
        return os.path.basename(self._filename)

    def abspath(self):
        """
        Return the normalized absolute filename

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.saved_file import SavedFile
            sage: from sage.repl.rich_output.output_catalog import OutputImagePng
            sage: png_file = OutputImagePng.example().png.filename(ext='.png')
            sage: SavedFile(png_file).abspath()
            '/....png'
        """
        return os.path.abspath(self._filename)

    def load(self):
        r"""
        Return the content of the saved file

        OUTPUT:

        Bytes (string in Python 2)

        EXAMPLES::

            sage: from sage.repl.saved_file import SavedFile
            sage: from sage.repl.rich_output.output_catalog import OutputImagePng
            sage: png_file = OutputImagePng.example().png.filename(ext='.png')
            sage: SavedFile(png_file).load()
            '\x89PNG\r\n\x1a\n\...\x00IEND\xaeB`\x82'
        """
        with open(self.abspath()) as f:
            return f.read()

    
class SavedFileWithShow(SavedFile):
    
    def show(self):
        """
        Display the saved file

        This method lets you view the just-saved file if we can either
        automatically figure out the file type or if the rich output
        type was specified to :class:`SavedFile`.

        EXAMPLES::

            sage: f = tmp_filename(ext='.png')
            sage: line2d([(-1,1), (1,-1)]).save(f).show()
        """
        from sage.repl.rich_output.buffer import OutputBuffer
        buf = OutputBuffer.from_file(self.abspath())
        rich_output = self._rich_output(buf)
        from sage.repl.rich_output import get_display_manager
        get_display_manager().display_immediately(
            RichReprHolder(self.basename(), rich_output))
