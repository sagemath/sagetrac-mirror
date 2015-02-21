# -*- encoding: utf-8 -*-
r"""
Base class for Backends

The display backends are the commandline, the SageNB notebook, the
ipython notebook, the Emacs sage mode, the Sage doctester, .... All of
these have different capabilities for what they can display.

To implement a new display backend, you need to subclass
:class:`BackendBase`. All backend-specific handlig of rich output
should be in :meth:`~BackendBase.displayhook` and
:meth:`~BackendBase.display_immediately`. See :class:`BackendSimple`
for an absolutely minimal example of a functioning backend.

You declare the types of rich output that your backend can handle in
:meth:`~BackendBase.supported_output`. There are two ways to then
display specific output types in your own backend.

* Directly use one of the existing output containers listed in
  :mod:`sage.repl.rich_output.output_catalog`. That is, test for the
  rich output type in :meth:`~BackendBase.display_immediately` and
  handle it.

* Subclass the rich output container to attach your backend-specific
  functionality. Then :meth:`~BackendBase.display_immediately` will
  receive instances of your subclass. See
  :class:`~sage.repl.rich_output.backend_test.BackendTest` for an
  example of how this is done.

You can also mix both ways of implementing different rich output types.

EXAMPLES::

    sage: from sage.repl.rich_output.backend_base import BackendSimple
    sage: backend = BackendSimple()
    sage: plain_text = backend.plain_text_formatter(range(10));  plain_text
    OutputPlainText container
    sage: backend.displayhook(plain_text, plain_text)
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject


class BackendBase(SageObject):

    def _repr_(self):
        """
        Return string representation of the backend

        Every backend must implement this method.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend._repr_()
            Traceback (most recent call last):
            ...
            NotImplementedError: derived classes must implement this method
        """
        raise NotImplementedError('derived classes must implement this method')

    def get_display_manager(self):
        """
        Return the display manager singleton

        This is a convenience method to access the display manager
        singleton.

        OUTPUT:

        The unique
        :class:`~sage.repl.rich_output.display_manager.DisplayManager`
        instance.

        EXAMPLES::
        
            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.get_display_manager()
            The Sage display manager using the doctest backend
        """
        from sage.repl.rich_output import get_display_manager
        return get_display_manager()

    def install(self, **kwds):
        """
        Hook that will be called once before the backend is used for the
        first time.

        The default implementation does nothing.

        INPUT:

        - ``kwds`` -- optional keyword arguments that are passed
          through by the
          :meth:`~sage.repl.rich_output.display_manager.DisplayManager.switch_backend`
          method.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.install()
        """
        pass

    def uninstall(self):
        """
        Hook that will be called once right before the backend is removed.

        The default implementation does nothing.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.uninstall()
        """
        pass

    def default_preferences(self):
        """
        Return the backend's display preferences

        Override this method to change the default preferences when
        using your backend.

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.preferences.DisplayPreferences`.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.default_preferences()
            Display preferences:
            * graphics is not specified
            * text is not specified
        """
        from sage.repl.rich_output.preferences import DisplayPreferences
        return DisplayPreferences()

    def supported_output(self):
        """
        Return the outputs that are supported by the backend.

        Subclasses must implement this method.

        OUTPUT:

        Iterable of output container classes, that is, subclass of
        :class:`~sage.repl.rich_output.output_basic.OutputBase`).  May
        be a list/tuple/set/frozenset. The order is ignored. Only used
        internally by the display manager.

        You may return backend-specific subclasses of existing output
        containers. This allows you to attach backend-specifc
        functionality to the output container.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.supported_output()
            Traceback (most recent call last):
            ...
            NotImplementedError: derived classes must implement this method
        """
        raise NotImplementedError('derived classes must implement this method')

    def max_width(self):
        """
        Return the number of characters that fit into one output line

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.max_width()
            79
        """
        return 79

    def newline(self):
        r"""
        Return the newline string.

        OUTPUT:

        String for starting a new line of output.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.newline()
            '\n'
        """
        return '\n'

    def _apply_pretty_printer(self, pretty_printer_class, obj):
        """
        Helper method to format ``obj`` as text

        INPUT:

        - ``pretty_printer_class`` -- subclass of
          :class:`sage.repl.display.pretty_print.SagePrettyPrinter`.

        - ``obj`` -- anything.

        OUTPUT:

        String.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: from sage.repl.display.pretty_print import SagePrettyPrinter
            sage: backend._apply_pretty_printer(SagePrettyPrinter, 1/2)
            '1/2'
        """
        import StringIO
        stream = StringIO.StringIO()
        printer = pretty_printer_class(
            stream, self.max_width(), self.newline())
        printer.pretty(obj)
        printer.flush()
        return stream.getvalue()

    def plain_text_formatter(self, obj):
        r"""
        Hook to override how plain text is being formatted.

        If the object does not have a ``_rich_repr_`` method, or if it
        does not return a rich output object
        (:class:`~sage.repl.rich_output.output_basic.OutputBase`),
        then this method is used to generate plain text output.

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.output_basic.OutputPlainText`
        containing the string representation of the object.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: out = backend.plain_text_formatter(range(30))
            sage: out
            OutputPlainText container
            sage: out.text
            buffer containing 139 bytes
            sage: out.text.get()
            '[0,\n 1,\n 2,\n 3,\n 4,\n 5,\n 6,\n 7,\n 8,\n 9,\n 
            10,\n 11,\n 12,\n 13,\n 14,\n 15,\n 16,\n 17,\n 18,\n 
            19,\n 20,\n 21,\n 22,\n 23,\n 24,\n 25,\n 26,\n 27,\n 
            28,\n 29]'
        """
        from sage.repl.display.pretty_print import SagePrettyPrinter
        plain_text = self._apply_pretty_printer(SagePrettyPrinter, obj)
        from sage.repl.rich_output.output_basic import OutputPlainText
        return OutputPlainText(plain_text)

    def ascii_art_formatter(self, obj):
        """
        Hook to override how ascii art is being formatted.

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.output_basic.OutputAsciiArt`
        containing the ascii art string representation of the object.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: out = backend.ascii_art_formatter(range(30))
            sage: out
            OutputAsciiArt container
            sage: out.ascii_art
            buffer containing 228 bytes
            sage: print(out.ascii_art.get())
            [                                                                              
            [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
            <BLANKLINE>
                                            ]
             22, 23, 24, 25, 26, 27, 28, 29 ]
        """
        from sage.repl.display.pretty_print import AsciiArtPrettyPrinter
        ascii_art = self._apply_pretty_printer(AsciiArtPrettyPrinter, obj)
        from sage.repl.rich_output.output_basic import OutputAsciiArt
        return OutputAsciiArt(ascii_art)

    def latex_formatter(self, obj):
        r"""
        Hook to override how Latex is being formatted.

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        Instance of
        :class:`~sage.repl.rich_output.output_basic.OutputLatex`
        containing the latex string representation of the object.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: out = backend.latex_formatter(1/2)
            sage: out
            OutputLatex container
            sage: out.latex
            buffer containing 45 bytes
            sage: out.latex.get()
            '\\newcommand{\\Bold}[1]{\\mathbf{#1}}\\frac{1}{2}'
            sage: out.mathjax()
            '<html><script type="math/tex; mode=display">\\newcommand{\\Bold}[1]{\\mathbf{#1}}\\frac{1}{2}</script></html>'
        """
        from sage.misc.latex import MathJax
        mathjax = MathJax().eval(obj, mode='plain', combine_all=True)
        from sage.repl.rich_output.output_basic import OutputLatex
        return OutputLatex(str(mathjax))

    def set_underscore_variable(self, obj):
        """
        Set the ``_`` builtin variable.

        By default, this sets the special ``_`` variable. Backends
        that organize the history differently (e.g. IPython) can
        override this method.

        INPUT:

        - ``obj`` -- result of the most recent evaluation.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.set_underscore_variable(123)
            sage: _
            123

            sage: 'foo'
            'foo'
            sage: _     # indirect doctest
            'foo'
        """
        import __builtin__
        __builtin__._ = obj
    
    def displayhook(self, plain_text, rich_output):
        """
        Backend implementation of the displayhook
        
        The value of the last statement on a REPL input line or
        notebook cell are usually handed to the Python displayhook and
        shown on screen.  By overriding this method you define how
        your backend handles output. The difference to the usual
        displayhook is that Sage already converted the value to the
        most suitable rich output container.

        Derived classes must implement this method.
        
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

        This method may return something, which is then returned from
        the display manager's
        :meth:`~sage.repl.rich_output.display_manager.DisplayManager.displayhook`
        method.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.displayhook(plain_text, plain_text)
            Traceback (most recent call last):
            ...
            NotImplementedError: derived classes must implement this method
        """
        return self.display_immediately(plain_text, rich_output)

    def display_immediately(self, plain_text, rich_output):
        """
        Show output without going back to the command line prompt.

        This method is similar to the rich output :meth:`displayhook`,
        except that it can be invoked at any time. Typically, it ends
        up being called by :meth:`sage.plot.graphics.Graphics.show`.

        Derived classes must implement this method.
        
        INPUT:

        Same as :meth:`displayhook`.

        OUTPUT:

        This method may return something so you can implement
        :meth:`displayhook` by calling this method. However, when
        called by the display manager any potential return value is
        discarded: There is no way to return anything without
        returning to the command prompt.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_base import BackendBase
            sage: backend = BackendBase()
            sage: backend.display_immediately(plain_text, plain_text)
            Traceback (most recent call last):
            ...
            NotImplementedError: derived classes must implement this method
        """
        raise NotImplementedError('derived classes must implement this method')


class BackendSimple(BackendBase):
    """
    Simple Backend

    This backend only supports plain text.
    
    EXAMPLES::

        sage: from sage.repl.rich_output.backend_base import BackendSimple
        sage: BackendSimple()
        simple
    """
    
    def _repr_(self):
        r"""
        Return string representation of the backend

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendSimple
            sage: backend = BackendSimple()
            sage: backend._repr_()
            'simple'
        """
        return 'simple'

    def supported_output(self):
        """
        Return the outputs that are supported by the backend.

        OUTPUT:

        Iterable of output container classes, that is, subclass of
        :class:`~sage.repl.rich_output.output_basic.OutputBase`).
        This backend only supports the plain text output container.

        EXAMPLES::

            sage: from sage.repl.rich_output.backend_base import BackendSimple
            sage: backend = BackendSimple()
            sage: backend.supported_output()
            {<class 'sage.repl.rich_output.output_basic.OutputPlainText'>}
        """
        from sage.repl.rich_output.output_basic import OutputPlainText
        return set([OutputPlainText])

    def display_immediately(self, plain_text, rich_output):
        """
        Show output without going back to the command line prompt.

        INPUT:

        Same as :meth:`~BackendBase.displayhook`.

        OUTPUT:

        This backend returns nothing, it just prints to stdout.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: from sage.repl.rich_output.backend_base import BackendSimple
            sage: backend = BackendSimple()
            sage: backend.display_immediately(plain_text, plain_text)
            Example plain text output
        """
        print(rich_output.text.get())

