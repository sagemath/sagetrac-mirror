# -*- encoding: utf-8 -*-
"""
The ``pretty_print`` command.

Works similar to the ``print`` function, except that it always tries
to use a rich output for an object. Only if that is not available it
is falling back on the plain text.

EXAMPLES::

    sage: pretty_print(1, 2, 3)
    <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}1 2 3</script></html>

Printing a graphics object just prints a string, whereas
:func:`pretty_print` does not print anything and just shows the
graphics instead::

    sage: print(plot(sin))
    Graphics object consisting of 1 graphics primitive
    sage: pretty_print(plot(sin))
"""


#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import types
import collections

from sage.structure.sage_object import SageObject
from sage.repl.rich_output import get_display_manager


class SequencePrettyPrinter(SageObject):

    def __init__(self, *args, **kwds):
        """
        Pretty Printer for Muliple Arguments.

        INPUT/OUTPUT:

        Same as :func:`pretty_print`, except that the number of
        arguments must be >= 2. Otherwise its not a sequence of things
        to print.

        EXAMPLES::

            sage: pretty_print(1, 2, 3)   # indirect doctest
            <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}1 2 3</script></html>
            sage: from sage.repl.rich_output.pretty_print import SequencePrettyPrinter
            sage: SequencePrettyPrinter(1, 2, 3).pretty_print()
            1 2 3
        """
        self.args = args
        assert len(self.args) >= 2
        self.kwds = kwds

    def is_homogeneous(self, common_type):
        """
        Return whether the pretty print items are homogeneous

        INPUT:
        
        - ``common_type`` -- a type.

        OUTPUT:

        Boolean. Whether all items to be pretty printed are of said
        type.

        EXAMPLES::

            sage: from sage.repl.rich_output.pretty_print import SequencePrettyPrinter
            sage: seq = SequencePrettyPrinter(1, 2, 3)
            sage: seq.is_homogeneous(Integer)
            True
            sage: seq.is_homogeneous(str)
            False
        """
        return all(isinstance(arg, common_type) for arg in self.args)
        
    def _concatenate_graphs(self):
        """
        Plot multiple graphs into a single plot

        OUTPUT:

        A graphics object.

        EXAMPLES::

            sage: from sage.repl.rich_output.pretty_print import SequencePrettyPrinter
            sage: plt = SequencePrettyPrinter(*list(graphs(3)))._concatenate_graphs()
            sage: type(plt)
            <class 'sage.plot.graphics.GraphicsArray'>
            sage: plt
            Graphics Array of size 2 x 4
        """
        import sage.graphs.graph_list as graphs_list
        return graphs_list.to_graphics_array(self.args, **self.kwds)

    def _concatenate_graphics(self):
        """
        Combine multiple graphics objects into one graphics array

        OUTPUT:

        A graphics array.

        EXAMPLES::

            sage: from sage.repl.rich_output.pretty_print import SequencePrettyPrinter
            sage: ga = SequencePrettyPrinter(*[Graphics()]*5)._concatenate_graphics()
            sage: type(ga)
            <class 'sage.plot.graphics.GraphicsArray'>
            sage: ga.nrows(), ga.ncols()
            (2, 4)
        """
        from sage.plot.plot import graphics_array
        return graphics_array(self.args, ncols=4, **self.kwds)
    
    def pretty_print(self):
        """
        Actually do the pretty print.

        EXAMPLES::

            sage: from sage.repl.rich_output.pretty_print import SequencePrettyPrinter
            sage: SequencePrettyPrinter(1, 2, 3).pretty_print()
            1 2 3

        The keyword arguments are only used the first time graphics
        output is generated::

            sage: seq = SequencePrettyPrinter(Graph(), Graph(), edge_labels=True)
            sage: seq.pretty_print()   # does not pass edge_labels to graphics object
            sage: seq._concatenate_graphs().show(edge_labels=True)
            Traceback (most recent call last):
            ...
            TypeError: matplotlib() got an unexpected keyword argument 'edge_labels'
        """
        from sage.plot.plot import Graphics
        from sage.graphs.graph import GenericGraph
        if self.is_homogeneous(GenericGraph):
            args = self._concatenate_graphs()
            kwds = dict()
        elif self.is_homogeneous(Graphics):
            args = self._concatenate_graphics()
            kwds = dict()
        else:
            args = self.args
            kwds = dict(self.kwds)
            kwds['concatenate'] = True
        get_display_manager().display_immediately(args, **kwds)


def pretty_print(*args, **kwds):
    r"""
    Pretty print the arguments in an intelligent way.
    
    For a single positional argument, this function chooses the
    highest-quality output supported by the user interface.
    
    For certain homogeneous multiple positional arguments a suitable
    combined graphical output is generated. In particular, graphs and
    plots are treated special.

    Otherwise this function will concatenate the textual
    representations. Latex output is preferred if none is specified
    via
    :meth:`~sage.repl.rich_output.display_manager.DisplayManager.preferences`.

    INPUT:

    - ``*args`` -- any number of positional arguments. The objects to
      pretty print. If the single argument is an iterator/generator
      then it is expanded.

    - ``**kwds`` -- optional keyword arguments that are passed to the
      rich representation. Examples include:

        - ``dpi`` - dots per inch

        - ``figsize``- [width, height] (same for square aspect)

        - ``axes`` - (default: True)

        - ``fontsize`` - positive integer

        - ``frame`` - (default: False) draw a MATLAB-like frame around
          the image

    EXAMPLES::

        sage: pretty_print(ZZ)
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}</script></html>

        sage: pretty_print("Integers = ", ZZ) # trac 11775
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}\verb|Integers|\phantom{\verb!x!}\verb|=| \Bold{Z}</script></html>

    To typeset LaTeX code as-is, use :class:`LatexExpr`::

        sage: pretty_print(LatexExpr(r"\frac{x^2 + 1}{x - 2}"))
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}\frac{x^2 + 1}{x - 2}</script></html>

    Iterators and generators are unwrapped::

        sage: iterator = iter(range(3));  iterator
        <listiterator object at 0x...>
        sage: pretty_print(iterator)
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}0 1 2</script></html>

    TESTS::

        sage: plt = plot(sin)
        sage: pretty_print(plt)             # graphics output
        sage: pretty_print(ZZ, 123, plt)    # latex output
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z} 123 \verb|Graphics|\phantom{\verb!x!}\verb|object|\phantom{\verb!x!}\verb|consisting|\phantom{\verb!x!}\verb|of|\phantom{\verb!x!}\verb|1|\phantom{\verb!x!}\verb|graphics|\phantom{\verb!x!}\verb|primitive|</script></html>
        sage: pretty_print(plt, plt)        # graphics output
    """
    if len(args) == 1 and isinstance(args[0], (types.GeneratorType, collections.Iterator)):
        args = tuple(args[0])

    # Support deprecation trac #18292
    if len(args) == 1:
        import sage.misc.html
        if sage.misc.html.WarnIfNotPrinted.skip_pretty_print(args[0]):
            return

    dm = get_display_manager()
    old_preferences_text = dm.preferences.text
    try:
        if dm.preferences.text is None:
            dm.preferences.text = 'latex'
        if len(args) == 0:
            pass
        elif len(args) == 1:
            dm.display_immediately(*args, **kwds)
        else:
            SequencePrettyPrinter(*args, **kwds).pretty_print()
    finally:
        dm.preferences.text = old_preferences_text
    

def show(*args, **kwds):
    """
    Alias for ``pretty_print``

    This function is an alias for :meth:`pretty_print`.

    INPUT/OUTPUT:

    See :meth:`pretty_print`. Except if the argument is a graph, in
    which case it is plotted instead.

    EXAMPLES::

        sage: show(1)
        <html><script type="math/tex">\newcommand{\Bold}[1]{\mathbf{#1}}1</script></html>
    """
    from sage.graphs.generic_graph import GenericGraph
    if len(args) == 1 and isinstance(args[0], GenericGraph):
        # Graphs are special, they ride the short bus...
        # Please, somebody help me get rid of this! #18289
        args[0].show()
        return
    pretty_print(*args, **kwds)
