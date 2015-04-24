"""
HTML Fragments

This module defines a HTML fragment class, which holds a piece of
HTML. This is primarily used in browser-based notebooks, though it
might be useful for creating static pages as well.
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import warnings
from sage.misc.latex import latex
from sage.misc.sage_eval import sage_eval
from sage.structure.sage_object import SageObject
from sage.misc.superseded import deprecation


class HtmlFragment(str, SageObject):
    r"""
    A HTML fragment.

    This is a piece of HTML, usually not a complete document.  For
    example, just a ``<div>...</div>`` piece and not the entire
    ``<html>...</html>``.

    .. automethod:: _rich_repr_
    """

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: h = html(1/2)
            sage: h._rich_repr_(dm)    # the doctest backend does not suppot html
        """
        OutputHtml = display_manager.types.OutputHtml
        if OutputHtml in display_manager.supported_output():
            return OutputHtml(self)


def math_parse(s):
    r"""
    Replace TeX-``$`` with Mathjax equations.
    
    Turn the HTML-ish string s that can have \$\$ and \$'s in it into
    pure HTML.  See below for a precise definition of what this means.

    INPUT:

    - ``s`` -- a string

    OUTPUT:

    A :class:`HtmlFragment` instance.

    Do the following:

    * Replace all ``\$ text \$``\'s by
      ``<script type="math/tex"> text </script>``
    * Replace all ``\$\$ text \$\$``\'s by
      ``<script type="math/tex; mode=display"> text </script>``
    * Replace all ``\ \$``\'s by ``\$``\'s.  Note that in
      the above two cases nothing is done if the ``\$``
      is preceeded by a backslash.
    * Replace all ``\[ text \]``\'s by
      ``<script type="math/tex; mode=display"> text </script>``

    EXAMPLES::

        sage: sage.misc.html.math_parse('This is $2+2$.')
        'This is <script type="math/tex">2+2</script>.'
        sage: sage.misc.html.math_parse('This is $$2+2$$.')
        'This is <script type="math/tex; mode=display">2+2</script>.'
        sage: sage.misc.html.math_parse('This is \\[2+2\\].')
        'This is <script type="math/tex; mode=display">2+2</script>.'
        sage: sage.misc.html.math_parse(r'This is \[2+2\].')
        'This is <script type="math/tex; mode=display">2+2</script>.'

    TESTS::

        sage: sage.misc.html.math_parse(r'This \$\$is $2+2$.')
        'This $$is <script type="math/tex">2+2</script>.'
    """
    # first replace \\[ and \\] by <script type="math/tex; mode=display">
    # and </script>, respectively.
    while True:
        i = s.find('\\[')
        if i == -1:
            break
        else:
            s = s[:i] + '<script type="math/tex; mode=display">' + s[i+2:]
            j = s.find('\\]')
            if j == -1:  # missing right-hand delimiter, so add one
                s = s + '</script>'
            else:
                s = s[:j] + '</script>' + s[j+2:]

    # Below t always has the "parsed so far" version of s, and s is
    # just the part of the original input s that hasn't been parsed.
    t = ''
    while True:
        i = s.find('$')
        if i == -1:
            # No dollar signs -- definitely done.
            return HtmlFragment(t + s)
        elif i > 0 and s[i-1] == '\\':
            # A dollar sign with a backslash right before it, so
            # we ignore it by sticking it in the parsed string t
            # and skip to the next iteration.
            t += s[:i-1] + '$'
            s = s[i+1:]
            continue
        elif i+1 < len(s) and s[i+1] == '$':
            # Found a math environment. Double dollar sign so display mode.
            disp = '; mode=display'
        else:
            # Found math environment. Single dollar sign so default mode.
            disp = ''

        # Now find the matching $ sign and form the html string.

        if len(disp) > 0:
            j = s[i+2:].find('$$')
            if j == -1:
                j = len(s)
                s += '$$'
            else:
                j += i + 2
            txt = s[i+2:j]
        else:
            j = s[i+2:].find('$')
            if j == -1:
                j = len(s)
                s += '$'
            else:
                j += i + 2
            txt = s[i+1:j]

        t += s[:i] + '<script type="math/tex%s">%s</script>'%(disp,
                      ' '.join(txt.splitlines()))
        s = s[j+1:]
        if len(disp) > 0:
            s = s[1:]
    return HtmlFragment(t)


class HTMLFragmentFactory(SageObject):

    def _repr_(self):
        """
        Return string representation

        OUTPUT:

        String.

        EXAMPLES::
        
            sage: html
            Create HTML output (see html? for details)
        """
        return 'Create HTML output (see html? for details)'
    
    def __call__(self, obj):
        r"""
        Construct a HTML fragment
     
        INPUT:
     
        - ``obj`` -- anything. An object for which you want a HTML
          representation.
     
        OUTPUT:
     
        A :class:`HtmlFragment` instance.
     
        EXAMPLES::
     
            sage: html('<hr>')
            '<hr>'
            sage: type(html('<hr>'))
            <class 'sage.misc.html.HtmlFragment'>

            sage: html(1/2)
            '<script type="math/tex">\\frac{1}{2}</script>'

            sage: html('<a href="http://sagemath.org">sagemath</a>')
            '<a href="http://sagemath.org">sagemath</a>'
        """
        # Prefer dedicated _html_() method
        try:
            result = obj._html_()
        except AttributeError:
            pass
        else:
            if not isinstance(result, HtmlFragment):
                warnings.warn('_html_() did not return a HtmlFragment')
                return HtmlFragment(result)
            else:
                return result
        # Otherwise: convert latex to html
        try:
            result = obj._latex_()
        except AttributeError:
            pass
        else:
            return math_parse('${0}$'.format(obj._latex_()))
        # If all else fails
        return math_parse(HtmlFragment(obj))
         
    def eval(self, s, locals=None):
        r"""
        Evaluate embedded <sage> tags

        INPUT:

        - ``s`` -- string.

        - ``globals`` -- dictionary. The global variables when
          evaluating ``s``. Default: the current global variables.

        OUTPUT:
     
        A :class:`HtmlFragment` instance.
     
        EXAMPLES::

            sage: a = 123
            sage: html.eval('<sage>a</sage>')
            '<script type="math/tex">123</script>'
            sage: html.eval('<sage>a</sage>', locals={'a': 456})
            '<script type="math/tex">456</script>'
        """
        if hasattr(s, '_html_'):
            deprecation(18292, 'html.eval() is for strings, use html() for sage objects')
            return s._html_()
        if locals is None:
            from sage.repl.user_globals import get_globals
            locals = get_globals()
        s = str(s)
        s = math_parse(s)
        t = ''
        while len(s) > 0:
            i = s.find('<sage>')
            if i == -1:
                 t += s
                 break
            j = s.find('</sage>')
            if j == -1:
                 t += s
                 break
            t += s[:i] + '<script type="math/tex">%s</script>'%\
                     latex(sage_eval(s[6+i:j], locals=locals))
            s = s[j+7:]
        return HtmlFragment(t)

    def table(self, x, header=False):
        r"""
        Generate a HTML table.  

        See :class:`~sage.misc.table.table`.

        INPUT:

        - ``x`` -- a list of lists (i.e., a list of table rows)

        - ``header`` -- a row of headers.  If ``True``, then the first
          row of the table is taken to be the header.

        OUTPUT:
     
        A :class:`HtmlFragment` instance.

        EXAMPLES::

            sage: output = html.table([(i, j, i == j) for i in [0..1] for j in [0..1]])
            sage: output
            '<div.../div>'
            sage: print(output)
            <div class="notruncate">
            <table class="table_form">
            <tbody>
            <tr class ="row-a">
            <td><script type="math/tex">0</script></td>
            <td><script type="math/tex">0</script></td>
            <td><script type="math/tex">\mathrm{True}</script></td>
            </tr>
            <tr class ="row-b">
            <td><script type="math/tex">0</script></td>
            <td><script type="math/tex">1</script></td>
            <td><script type="math/tex">\mathrm{False}</script></td>
            </tr>
            <tr class ="row-a">
            <td><script type="math/tex">1</script></td>
            <td><script type="math/tex">0</script></td>
            <td><script type="math/tex">\mathrm{False}</script></td>
            </tr>
            <tr class ="row-b">
            <td><script type="math/tex">1</script></td>
            <td><script type="math/tex">1</script></td>
            <td><script type="math/tex">\mathrm{True}</script></td>
            </tr>
            </tbody>
            </table>
            </div>

            sage: output = html.table([(x,n(sin(x), digits=2)) for x in [0..3]], 
            ....:                     header=["$x$", "$\sin(x)$"])
            sage: print(output)
            <div class="notruncate">
            <table class="table_form">
            <tbody>
            <tr>
            <th><script type="math/tex">x</script></th>
            <th><script type="math/tex">\sin(x)</script></th>
            </tr>
            <tr class ="row-a">
            <td><script type="math/tex">0</script></td>
            <td><script type="math/tex">0.00</script></td>
            </tr>
            <tr class ="row-b">
            <td><script type="math/tex">1</script></td>
            <td><script type="math/tex">0.84</script></td>
            </tr>
            <tr class ="row-a">
            <td><script type="math/tex">2</script></td>
            <td><script type="math/tex">0.91</script></td>
            </tr>
            <tr class ="row-b">
            <td><script type="math/tex">3</script></td>
            <td><script type="math/tex">0.14</script></td>
            </tr>
            </tbody>
            </table>
            </div>
        """
        from table import table
        return table(x, header_row=header)._html_()

    def iframe(self, url, height=400, width=800):
        r"""
        Generate an iframe HTML fragment

        INPUT:

        - ``url`` -- string. A url, either with or without URI scheme
          (defaults to "http"), or an absolute file path.

        - ``height`` -- the number of pixels for the page height.
          Defaults to 400.

        - ``width`` -- the number of pixels for the page width.
          Defaults to 800.

        OUTPUT:
     
        A :class:`HtmlFragment` instance.

        EXAMPLES::

            sage: html.iframe("sagemath.org")
            '<iframe height="400" width="800"
            src="http://sagemath.org"></iframe>'
            sage: html.iframe("http://sagemath.org",30,40)
            '<iframe height="30" width="40"
            src="http://sagemath.org"></iframe>'
            sage: html.iframe("https://sagemath.org",30)
            '<iframe height="30" width="800"
            src="https://sagemath.org"></iframe>'
            sage: html.iframe("/home/admin/0/data/filename")
            '<iframe height="400" width="800"
            src="file:///home/admin/0/data/filename"></iframe>'
            sage: html.iframe('data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA'
            ... 'AUAAAAFCAYAAACNbyblAAAAHElEQVQI12P4//8/w38GIAXDIBKE0DHxgljNBA'
            ... 'AO9TXL0Y4OHwAAAABJRU5ErkJggg=="')
            '<iframe height="400" width="800" 
            src="http://data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAUAAAAFCAYAAACNbyblAAAAHElEQVQI12P4//8/w38GIAXDIBKE0DHxgljNBAAO9TXL0Y4OHwAAAABJRU5ErkJggg==""></iframe>'
        """
        if url.startswith('/'):
            url = 'file://{0}'.format(url)
        elif '://' not in url:
            url = 'http://{0}'.format(url)
        return HtmlFragment('<iframe height="{0}" width="{1}" src="{2}"></iframe>'
                            .format(height, width, url))


html = HTMLFragmentFactory()
html.__doc__ = HTMLFragmentFactory.__call__.__doc__
