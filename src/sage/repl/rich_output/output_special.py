# -*- encoding: utf-8 -*-
r"""
Special Output Types

This module defines rich output types that do not correspond to a
particular type of media file, but instead are hooks for the
displayhook to perform certain actions.
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

from sage.repl.rich_output.output_basic import OutputBase



class OutputSavedFile(OutputBase):

    def __init__(self, filename, auto_download=True):
        """
        A Saved File

        This is the rich representation type of
        :class:`~sage.repl.saved_file.SavedFile`. A user interface
        (usually, web browser) for Sage running on a remote computer
        shall react to this rich output container by downloading the
        file. Ideally the download is initiated without further user
        actions, but showing a download link is also acceptable.

        INPUT:

        - ``filename`` -- string. The filename. Must point to an
          existing file.

        - ``auto_download`` -- boolean (default: ``True``). Whether
          remote files should be downloaded automatically.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSavedFile
            sage: OutputSavedFile.example()
            OutputSavedFile container
        """
        self.filename = os.path.abspath(filename)

    @classmethod
    def example(cls):
        r"""
        Construct a sample saved file container

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT:

        An instance of :class:`OutputSavedFile`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputSavedFile
            sage: OutputSavedFile.example()
            OutputSavedFile container
            sage: OutputSavedFile.example().filename
            '/.../local/share/sage/ext/doctest/rich_output/example.png'
        """
        from sage.env import SAGE_EXTCODE
        filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.png')
        return OutputSavedFile(filename)

