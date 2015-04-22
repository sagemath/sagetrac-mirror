# -*- encoding: utf-8 -*-
"""
Portable Document Format (PDF) Wrapper

This object is a wrapper around PDF documents, usually generated from
(pdf)latex. The wrapper object can then be manipulated on the command
line, for example it allows you to preview pages as svg or bitmap
graphics.
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
from sage.structure.sage_object import SageObject
from sage.repl.rich_output.buffer import OutputBuffer
from sage.misc.cachefunc import cached_method

from sage.interfaces.cmdline import Tool, ToolMissingException
from sage.interfaces.cmdline.pdf import (
    pdftocairo_svg, 
    pdf2svg, pdfcrop,
    convert_svg, convert_png,
    ghostscript_png,
)


class PortableDocumentFormat(SageObject):

    def __init__(self, pdf):
        """
        Wrapper for a PDF Document

        This wraps a in-memory pdf document. See also
        :meth:`from_file` if you want to wrap a pdf document that is
        already saved in a file.

        INPUT:

        - ``pdf`` -- bytes (string). A PDF document.

        EXAMPLES::

            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: pdf = PortableDocumentFormat._example();  pdf
            PDF document (4285 bytes)
            sage: type(pdf)
            <class 'sage.repl.portable_document_format.PortableDocumentFormat'>
        """
        self._pdf = OutputBuffer(pdf)
        if not self._pdf.get().startswith('%PDF-'):
            raise ValueError('not a PDF document')

    @classmethod
    def from_file(cls, filename):
        """
        Construct Wrapper for a PDF Document

        INPUT:

        - ``filename`` -- string. The name of a PDF file to load.

        OUTPUT:

        :class:`PortableDocumentFormat` instance.

        EXAMPLES::

            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: pdf = PortableDocumentFormat._example()
            sage: PortableDocumentFormat.from_file(pdf.filename())
            PDF document (4285 bytes)
        """
        return cls(OutputBuffer.from_file(filename))

    @classmethod
    def _example(cls):
        """
        Construct an Example PDF
        
        OUTPUT:
        
        :class:`PortableDocumentFormat`.
         
        EXAMPLES::
        
            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: PortableDocumentFormat._example()
            PDF document (4285 bytes)
        """
        from sage.env import SAGE_EXTCODE
        filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.pdf')
        return cls.from_file(filename)

    @cached_method
    def filename(self):
        """
        Return the filename of the PDF document

        If necessary, the document is saved first.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: pdf = PortableDocumentFormat._example()
            sage: pdf.filename()
            '/.../local/share/sage/ext/doctest/rich_output/example.pdf'
        """
        return self._pdf.filename(ext='pdf')

    def save(self, filename):
        """
        Save the PDF document

        INPUT:

        - ``filename`` -- string. The filename to save under.

        OUTPUT:

        A :class:`sage.repl.saved_file.SavedFile` instance.

        EXAMPLES::

            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: pdf = PortableDocumentFormat._example()
            sage: pdf.save(tmp_filename())
            Saved file at /...

        The returned wrapper object lets you access the file
        metadata::

            sage: filename = tmp_filename(ext='filename.pdf')
            sage: saved_file = pdf.save(filename)
            sage: saved_file.abspath()
            '/...filename.pdf'
        """
        self._pdf.save_as(filename)
        from sage.repl.saved_file import SavedFile
        return SavedFile(filename)
        
    def _repr_(self):
        """
        Return a string representation
        
        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: pdf = PortableDocumentFormat._example()
            sage: pdf._repr_()
            'PDF document (4285 bytes)'
        """
        return 'PDF document ({0} bytes)'.format(len(self._pdf.get()))

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: pdf = PortableDocumentFormat._example()
            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: pdf._rich_repr_(dm)
            OutputImagePdf container
        """
        dpi = kwds.get('dpi', 150)
        if display_manager.preferences.graphics == 'disable':
            return
        OutputImagePdf = display_manager.types.OutputImagePdf
        if OutputImagePdf in display_manager.supported_output():
            return OutputImagePdf(self._pdf)
        ### This does not work properly yet because of id collisions
        #OutputImageSvg = display_manager.types.OutputImageSvg
        #if OutputImageSvg in display_manager.supported_output():
        #    try:
        #        svg = self.to_svg()
        #    except ToolMissingException as exc:
        #        exc.convert_to_warning()
        #    return OutputImageSvg(svg)
        OutputImagePng = display_manager.types.OutputImagePng
        if OutputImagePng in display_manager.supported_output():
            try:
                png = self.to_png(dpi=dpi)
            except ToolMissingException as exc:
                exc.convert_to_warning()
            return OutputImagePng(png)

    def show(self, dpi=150):
        """
        Show this document immediately.

        This method attempts to display the document immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        INPUT:

        - ``dpi`` -- integer. The resolution at which to raster the
          pdf file if necessary.

        EXAMPLES::

            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: pdf = PortableDocumentFormat._example()
            sage: pdf.show(dpi=75)
        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        dm.display_immediately(self, dpi=dpi)

    def trim(self):
        """
        Cut off white margins

        OUTPUT:

        A new :class:`PortableDocumentFormat` with the margins
        trimmed.

        EXAMPLES::
        
            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: pdf = PortableDocumentFormat._example()
            sage: pdf.trim()   # optional - pdfcrop
            Launched pdf viewer for PDF document (4863 bytes)
        """
        tool = Tool.require(pdfcrop)
        trimmed = tool(self._pdf.get())
        return self.__class__(trimmed)
        
    def to_svg(self, trim=False):
        r"""
        Convert to SVG

        INPUT:

        - ``trim`` -- boolean (optional). Whether to cut off
          white margins. Ignored if the ``pdfcrop`` command line tool
          is not installed.

        OUTPUT:

        String. The first page of the PDF document converted to SVG.

        EXAMPLES::

            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: pdf = PortableDocumentFormat._example()
            sage: pdf.to_svg()
            '<?xml...</svg>\n'
        """
        if trim not in [True, False]:
            raise ValueError('trim must be a boolean, got {0}'.format(trim))
        pdf = self._pdf.get()
        if trim and pdfcrop.is_installed():
            pdf = pdfcrop(pdf)
        tool = Tool.require(pdftocairo_svg, pdf2svg, convert_svg)
        return tool(pdf)

    def to_png(self, trim=False, dpi=75):
        r"""
        Convert to PNG

        INPUT:

        - ``trim`` -- boolean (default: ``True``). Whether to cut off
          white margins. Ignored if the ``pdfcrop`` command line tool
          is not installed.

        - ``dpi`` -- integer (default: 75). Resolution of the raster
          graphics output in pixels per inch.

        OUTPUT:

        String. The first page of the PDF document converted to PNG.

        EXAMPLES::

            sage: from sage.repl.portable_document_format import PortableDocumentFormat
            sage: pdf = PortableDocumentFormat._example()
            sage: pdf.to_png()
            '\x89PNG...'
        """
        if trim not in [True, False]:
            raise ValueError('trim must be a boolean, got {0}'.format(trim))
        pdf = self._pdf.get()
        if trim and pdfcrop.is_installed():
            pdf = pdfcrop(pdf)
        tool = Tool.require(convert_png, ghostscript_png)
        return tool(pdf, dpi=dpi)
    
    

