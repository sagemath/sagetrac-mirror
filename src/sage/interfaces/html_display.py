r"""
HTML display
------------

The :class:`HTMLDisplay` can be used to display ``sage`` objects in a web
browser.  In other words, this class provides pretty printing on wheels.  This
can be useful for ``sage`` objects which are hard to display in the terminal,
either because of their size of because they require special characters or
extensive use of colour. The class is primarily intended to be used by
developers, rather than being called from the command line.

The classes defined below are intended to be use by other classes to construct
HTML rather than being calle directly. The following three classes are
implemented:

- :class:`HTMLElement`   a base class which implements the basic functionality for the
                         classes which can be added to sage.
- :class:`HTMLTable`     a class for turning an matrix or array into a
                         (scrollable) html table with prescribed row and column labels
- :class:`HTMLDisplay`   a container for constructing an HTML page, section by
                         section, in the ``sage`` documentation format.

EXAMPLES::

    sage: from sage.interfaces.html_display import HTMLDisplay
    sage: page=HTMLDisplay('Prime numbers')
    sage: page
    <html_display: page>Prime numbers
    sage: page.add_section('<B>Theorem</B> There are infinitely many prime numbers<br>','The number of primes')
    sage: page.add_section('<ul>'+'\n'.join('<li>%s</li>'%p for p in primes(100))+'</ul>','Prime numbers less than 100')
    sage: page.html_display()      # not tested

The derived classes of :class:`HTMLElement` can all be added to an HTML page.
Now we add a table to the web page::

    sage: from sage.interfaces.html_display import HTMLTable
    sage: table=HTMLTable(table=matrix(15), title="A matrix", row_label=range(14), col_label=range(15))
    sage: page.add_section(table)
    sage: page
    <html_display: page>Prime numbers
       0. The number of primes
       1. Prime numbers less than 100
       2. A matrix
    sage: page.html_display()  # not tested

Here is a more realistic example which wraps a character table for a finite
group::

    sage: from sage.interfaces.html_display import HTMLTable
    sage: def CharacterTable(G):
    ...       alphabet='abcdefghijklnmopqrstuvwxyz'
    ...       cc_orders=[c.order() for c in G.conjugacy_classes_representatives()]
    ...       cc_labels=['%d%s' %(cc_orders[c], alphabet[cc_orders[:c].count(cc_orders[c])]) for c in xrange(len(cc_orders))]
    ...       return HTMLTable(title='Character table of %s' % G,
    ...                        table=G.character_table(),
    ...                        col_label=cc_labels,
    ...                        row_label=['X%d'%d for d in xrange(len(cc_labels))])
    sage: ct=CharacterTable(SymmetricGroup(6))
    sage: ct.html_display()   # not tested



..NOTE::

Unlike other ``sage`` interfaces``, :class:`HTMLDisplay` is a one way interface:
the web pages opened by :class:`HTMLDisplay` cannot interact with ``sage``. It
is intended to be used from the command line rather than from inside the
:class:`notebook`.

AUTHORS:

- Andrew Mathas (2013-12-17): initial version

"""

#*****************************************************************************
#  Copyright (C) 2012 Andrew Mathas andrw dot mathas at sydney dot edu dot au
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from __builtin__ import str
from os import environ
from time import localtime,strftime

from sage.misc.flatten import flatten
from sage.misc.latex import latex
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.element import is_Matrix
from sage.structure.sage_object import SageObject
from sage.version import version as sage_version

class HTMLElement(SageObject):
    def __init__(self, html_tag=None, **options):
        r"""
        Creates the :class:`HTMLElement` ``self`` and sets all of its 
        options as specifies in ``options``. The ``html_tag`` describes the type
        of html block that ``self`` represents. The html_tag is mandatory!

        EXAMPLES::

            sage: from sage.interfaces.html_display import HTMLTable, HTMLElement
            sage: table=HTMLTable(title='A table', scrollable_rows=True, table=[[1,2,3],[4,5,6]])  # indirect doctest
            sage: table.html_tag
            'table'
            sage: table
            <html_display: table>A table
            sage: HTMLElement()
            Traceback (most recent call last):
            ...
            ValueError: An html_tag must be specified!
        """
        if html_tag is None:
            raise ValueError('An html_tag must be specified!')
        self.html_tag=html_tag # identifies this html tag/object

        # Every HTMLElement has the following dictionaries which keep track of
        # the header information that they need. When one HTMLElement is added
        # to another any new header information is also added and these
        # dictionaries ensure that we do not duplicate any of this data.
        self._css={}        # dictionary of css style definitions
        self._javascript={} # dictionary of javascript code
        self._head={}       # dictionary of header lines other than css and javascript
        self._document_ready={} # dictionary for javascript to be executed at "document ready"

        # To process options we first set all of the default options and
        # then overwrite them with the specified options
        if not hasattr(self, '_default_options'):
            self._default_options={'title':self.html_tag}
        self.options(**self._default_options)

        if options!={}: self.options(**options)


    def __getitem__(self,key):
        r'''
        A hack to allow easy formatting of the form:
            print '%s(tag|indent)' % self
        which will print the tag with the specified indentation, which is
        optional.

        EXAMPLES::

            sage: from sage.interfaces.html_display import HTMLDisplay
            sage: page=HTMLDisplay(title="A nice page")
            sage: print '<title>%(title)</title>' % page   # indirect doctest
            <title>A nice page</title>
        '''
        if '|' in key:
            tag, indent = key.split('|')
        else:
            tag, indent = key, '0'

        # make the indent into an indent string of the righht length
        try:
            indent='  '*int(indent)
        except ValueError:
            raise ValueError('%s is not a valid key for %s\n' % (key, self) )

        # change the tag to _'_+tag and pick the corresponding attribute/method
        tag='_'+tag
        try:
            tag_value=getattr(self,tag)
        except AttributeError:
            raise AttributeEror('%s is not an attribute of %s\n' % (tag, self))

        if callable(tag_value):
            return tag_value(indent)
        elif isinstance(tag_value, list):
            return indent.join('%s\n'% line for line in tag_value) if tag_value!=[] else ''
        else:
            return indent+tag_value


    def html_display(self):
        r"""
        Any :class:`HTMLElement` can be displayed on a web page by applying its
        method :meth:`html_display`. It does this by wrapping itself inside a
        :class:`HTMLDisplay` element as a single section.

        EXAMPLES:

            sage: from sage.interfaces.html_display import HTMLTable
            sage: table=HTMLTable(table=[[1,2],[3,4]])
            sage: table.html_display()    # not tested
        """
        page=HTMLDisplay(title=self._title)
        page.add_section(self)
        page.html_display()

    def options(self, *get_option, **set_option):
        r"""
        Sets and displays the options for the array. When the options are set
        we also do some post processing to ensure that they do what we want.

        When called with no arguments :meth:`options` prints the list of the
        current options of ``self``. The options corresponding to the strings in
        ``get_option`` are returned - as a single value if ``get_option`` has
        length one and otherwise as a dictionary of the requested options. Any
        options in ``set_option`` are set to the corresponding value in this
        dictionary.

        EXAMPLES::

            sage: from sage.interfaces.html_display import HTMLTable
            sage: table=HTMLTable(table=[[1,2],[3,4]])
            sage: table.options()
            Current options for <html_display: table>Table:
              - alternate_colours = True
              - cell_heights = 1em
              - cell_widths = None
              - col_label = None
              - col_label_math_mode = True
              - column_footer = False
              - column_header = True
              - entries_math_mode = True
              - row_column_header = 
              - row_header = True
              - row_label = None
              - row_label_math_mode = True
              - scrollable_columns = True
              - scrollable_rows = True
              - table_height = 90%
              - table_width = 90%
              - title = Table
            sage: table.options('title', 'scrollable_rows')
            {'scrollable_rows': True, 'title': 'Table'}
            sage: table.options(title='A nice web page', scrollable_rows=True)
        """
        # print the current settings for the options in ``option``
        if get_option==() and set_option=={}:
            print 'Current options for %s:' % self
            for opt in sorted(self._default_options.keys()):
                print '  - %s = %s' % (opt, getattr(self, '_'+opt))

        # print the current settings for the options in ``option``
        if len(get_option)>0:
            for opt in get_option:
                if not opt in self._default_options:
                    raise ValueError('%s is not a valid option for %s\n'%(opt,self))
            if len(get_option)==1:
                return getattr(self,'_'+get_option[0])
            else:
                return {option: getattr(self,'_'+option) for option in get_option}

        # set the specified options, but process them in alphabetical order as
        # some options need to be set before others and conveniently
        # alphabetical order is compatible with this.
        for opt in sorted(set_option.keys()):
            if opt in self._default_options:
                setattr(self, '_'+opt, set_option[opt])
            else:
                raise ValueError('%s is not a valid option for %s\n'%(opt,self))


    def _repr_(self):
        r"""
        A representation of the :class:`HTMLElement` ``self``.

        EXAMPLES::

            sage: from sage.interfaces.html_display import HTMLTable
            sage: table=HTMLTable(table=[[1,2],[3,4]])
            sage: table   #indirect doc test
            <html_display: table>Table
        """
        return '<html_display: %s>%s' % (self.html_tag, self._title)


    def body(self):
        r"""
        This is a place holder as very :class:`HTMLElement` must have a method
        :meth:`body` which returns a string for displaying the element on
        the page.

        EXAMPLE::

            sage: from sage.interfaces.html_display import HTMLElement
            sage: HTMLElement('ul').body()
            Traceback (most recent call last):
            ...
            NotImplementedError: body has not been implemented for <html_display: ul>ul
        """
        raise NotImplementedError('body has not been implemented for %s' %self)


    def add_meta_data(self, element):
        r"""
        Adds all of the meta data from ``element`` to ``self``.

        EXAMPLES::

            sage: from sage.interfaces.html_display import HTMLDisplay, HTMLTable
            sage: page=HTMLDisplay(); table=HTMLTable(table=[[1,2,3]])
            sage: table.add_to_meta_data('css','table_class','table_class.td {text-align: center;}')
            sage: page.add_meta_data(table)
            sage: page._css
            {'table_class': ['table_class.td {text-align: center;}']}
        """
        for name in element._head:
            self.add_to_meta_data('head', name, element._head[name])
        for name in  element._css:
            self.add_to_meta_data('css', name, element._css[name])
        for name in  element._javascript:
            self.add_to_meta_data('javascript', name, element._javascript[name])
        for name in  element._document_ready:
            self.add_to_meta_data('document_ready', name, element._document_ready[name])


    def add_to_meta_data(self, meta_type, name, data):
        r"""
        Add meta data to the :class:`HTMLElement` ``self``. Here::
        - ``meta_type`` is one of ``css``, ``head``, ``javascript`` or ``document_ready``
        - ``name`` is a string which uniquely identifies the data (so that it
              will not be added twice)
        - ``data`` is the actual meta data.

        The ``css`` and ``javascript`` ``meta_types`` are clear. Data of type
        ``head`` are typically links and meta data put into the html header and
        ``document_ready`` is javascript to be executed inside the html documents
        call to "$(document).ready(...)".

        EXAMPLES:

            sage: from sage.interfaces.html_display import HTMLDisplay
            sage: page=HTMLDisplay()
            sage: page.add_to_meta_data('css', 'stuff', 'td { font-weight: bold; }')
            sage: page._css
            {'stuff': ['td { font-weight: bold; }']}
            sage: page.add_to_meta_data('javascript','stuff',"alert('page ready!');")
            sage: page._javascript
            {'stuff': ["alert('page ready!');"]}

        As the example indicates, the ``name`` identifies only need to be unique
        within each ``meta_type``.
        """
        meta=getattr(self,'_'+meta_type)
        if not name in meta:
            if isinstance(data, str):
                meta[name]=data.split('\n')
            else:
                meta[name]=data


class HTMLTable(HTMLElement):
    r"""
    This class constructs a HTML table from a matrix or two dimensional array,
    which is suitable for inclusion in a :class:`HTMLDisplay`.

    As an example, we construct a :class:`HTMLTable` for the character table of
    a group.

    EXAMPLES::

        sage: from sage.interfaces.html_display import HTMLTable
        sage: def CharacterTable(G):
        ...       alphabet='abcdefghijklnmopqrstuvwxyz'
        ...       cc_orders=[c.order() for c in G.conjugacy_classes_representatives()]
        ...       cc_labels=['%d%s' %(cc_orders[c], alphabet[cc_orders[:c].count(cc_orders[c])]) for c in xrange(len(cc_orders))]
        ...       return HTMLTable(title='Character table of %s' % G,
        ...                        table=G.character_table(),
        ...                        col_label=cc_labels,
        ...                        row_label=['X%d'%d for d in xrange(len(cc_labels))])
        sage: ct=CharacterTable(CyclicPermutationGroup(12))
        sage: ct.html_display()   # not tested
        sage: ct.options()
        Current options for <html_display: table>Character table of Cyclic group of order 12 as a permutation group:
          - alternate_colours = True
          - cell_heights = 1em
          - cell_widths = None
          - col_label = ['1a', '12a', '6a', '4a', '3a', '12b', '2a', '12c', '3b', '4b', '6b', '12d']
          - col_label_math_mode = True
          - column_footer = False
          - column_header = True
          - entries_math_mode = True
          - row_column_header = 
          - row_header = True
          - row_label = ['X0', 'X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'X11']
          - row_label_math_mode = True
          - scrollable_columns = True
          - scrollable_rows = True
          - table_height = 90%
          - table_width = 90%
          - title = Character table of Cyclic group of order 12 as a permutation group
    """
    def __init__(self, **options):
        # default table options
        self._default_options=dict(
                alternate_colours=True,
                cell_heights='1em',
                cell_widths=None,
                col_label=None,
                col_label_math_mode=True,
                column_footer=False,
                column_header=True,
                entries_math_mode=True,
                row_header=True,
                row_label=None,
                row_label_math_mode=True,
                row_column_header='',
                scrollable_columns=True,
                scrollable_rows=True,
                table_height='90%',
                table_width='90%',
                title='Table'
        )
        try:
            self._table=options.pop('table')
        except KeyError:
            raise ValueError('a table must be specified!\n')

        super(HTMLTable, self).__init__('table', **options)


    def body(self, url_include):
        r""" 
        Return a list which contains the lines in the table body and, at the
        same time ensure that the meta data is set correctly for the table.

        The format for the scrollable table is:
            <table id="Open_Text_Arrays" class="dataTable">
                <thead>
                   <tr><th>Col1</th><th>Col2</th>...</tr>
                </thead>
                <tbody>
                   <tr><th>Col1</th><td>Col2</td>...</tr>
                </tbody>
                <tfoot>
                   <tr><th>Col1</th>...</tr>
                </tfoot>
            </table>
        """

        table_id='123456'
        body=[]
        if self._scrollable_rows or self._scrollable_columns:
            # add javascript for scrolling table
            self.add_to_meta_data('css','scrollable_table',css_for_scrollable_table)
            self.add_to_meta_data('document_ready','scrollable_table',javascript_for_scrollable_table)
            # load the necessary javascript
            head ='<link rel="stylesheet" type="text/css" href="{url_include}_static/jquery.dataTables.css" />\n'
            head+='<script type="text/javascript" src="{url_include}_static/jquery.dataTables.min.js"></script>\n'
            head+='<script type="text/javascript" src="{url_include}_static/FixedColumns.min.js"></script>\n'
            self.add_to_meta_data('head','scrolling table', head.format(url_include=url_include))
        else:
            # ensure that javascript for scrolling tables is not being added
            if 'scrollable_table' in self._javascript:
                self._javascript.pop('scrollable_table')

        body.append('<table id="html_table_display_%s" class="dataTable">' % table_id)

        # column headings
        if self._column_header:
            body.append('    <thead>')
            body.append('      <tr>')
            body.append('        <th>%s</th>' % self._row_column_header)
            if self._col_label_math_mode:
                col_label=lambda col: '        <th> \( %s \) </th>' % latex(col)
            else:
                col_label=lambda col: '        <th> %s </th>' % col
            body+=[col_label(col) for col in self._col_label]
            body.append('    </tr>')
            body.append('  </thead>')

        # we are ready to dump the actual table, but first some helpers
        body.append('  <tbody>')

        # row headings
        if self._row_header:
            if self._row_label_math_mode:
                row_label=lambda label: '      <td>\( %s \) </td>' % latex(self._row_label[label])
            else:
                row_label=lambda label: '      <td>%s </td>' % self._row_label[label]
        else: row_label=lambda label:''

        # now define a short-cut for printing the table entries. This
        # depends upon whether self._table is a matrix or an array and
        # whether the entries should be typeset using latex
        if is_Matrix(self._table):
            if self._entries_math_mode:
                entry=lambda table, row, col: '      <td>\( %s \) </td>' % latex(table[row,col])
            else:
                entry=lambda table, row, col: '      <td>%s </td>' % table[row,col]
        else:
            if self._entries_math_mode:
                entry=lambda table, row, col: '      <td>\( %s \) </td>' % latex(table[row][col])
            else:
                entry=lambda table, row, col: '      <td>%s </td>' % table[row][col]

        # the table contents
        for row in xrange(len(self._row_label)):
            body.append('    <tr>')
            body.append(row_label(row))
            for col in xrange(len(self._col_label)):
                body.append(entry(self._table,row,col))
            body.append('    </tr>')
        body.append('  </tbody>')

        # column footer
        if self._column_footer:
            body.append('  <tfoot>')
            body.append('    <tr>')
            if self._col_label_math_mode:
                col_label=lambda col: '      <th> \( %s \) </th>' % latex(col)
            else:
                col_label=lambda col: '      <th> %s </th>' % col
            body+=[col_label(col) for col in self._col_label]
            body.append('    </tr>')
            body.append('  </tfoot>')

        # finally, close off the table and the div
        body.append('</table>     <!-- close  html_table_display_%s -->' % table_id)

        return body


class HTMLDisplay(HTMLElement):
    r"""
    The wrapper class for creating an HTML page with the ``sage`` logo.

    Options
    -------
    - ``title``       the title of the web page
    - ``url_include`` the default base URL
    - ``export``      if ``True`` a web exportable file is generated (default: ``False``)
    - ``html_file``   name of the html file (defaults to temporary file)

    All of these options default to
    """

    def __init__(self, title='', **options):
        self._default_options=dict(
                export=False,
                html_file=None,
                url_include=None
        )
        super(HTMLDisplay, self).__init__('page', **options)

        if self._url_include is None:
            if self._export:
                self._url_include='http://www.sagemath.org/doc/output/html/en'
            else:
                self._url_include='//file:%s/devel/sage/doc/output/html/en/' % environ['SAGE_ROOT']

        self._title=title   # page title
        self._navigation=[] # a list of navigation items to include on the pages
        # Finally, an ordered list of the sections on the page. Each entry is
        # a dictionary consisting of the section title and the text it contains.
        # These can later be modified by add_section when redisplaying the page.
        self._sections=[]

        # the default header and footer data
        self.add_to_meta_data('head','generic',generic_meta_data.format(
                           url_include=self._url_include, version=sage_version))
        self.add_to_meta_data('document_ready','sphinx_sidebar', sphinx_sidebar_toggle)


    def _repr_(self):
        r"""
        HTMLDisplay pages print showing their list of sections.

        EXAMPLES::

            sage: from sage.interfaces.html_display import HTMLDisplay
            sage: page=HTMLDisplay('Prime numbers')
            sage: page.add_section('<B>Theorem</B> There are infinitely many prime numbers<br>','The number of primes')
            sage: page.add_section('<ul>'+'\n'.join('<li>%s</li>'%p for p in primes(100))+'</ul>','Prime numbers less than 100')
            sage: page  # indirect doctest
            <html_display: page>Prime numbers
               0. The number of primes
               1. Prime numbers less than 100
        """
        return '<html_display: page>%s\n%s' %(self._title,
                  '\n'.join('  %2d. %s'%(d, self._sections[d]['title']) for d in xrange(len(self._sections))))


    def _html_navigation(self, indent):
        nav =indent+'<div class="related">\n'
        nav+=indent+'  <h3>Navigation</h3>\n'
        nav+=indent+'  <ul>\n'
        nav+=indent.join('    %s\n' for line in self._navigation)
        nav+=indent+'    <li><a href="%(url_include)sindex.html">\n' % self
        nav+=indent+'      <img src="%(url_include)s_static/sagelogo.png" style="vertical-align:middle" alt="Sage Logo" title="Sage Logo">\n' % self
        nav+=indent+'    </a></li>\n'
        nav+=indent+'    <li> <a href="#">HTML display</a>&raquo; %(title)s</li>\n' % self
        nav+=indent+'  </ul>\n'
        nav+=indent+'</div>   <!-- close related -->\n'
        return nav


    def _html_meta(self, indent):
        r"""
        Returns a string which defines al of the header data for inclusion on
        the HTML page. The integer ``indent``defines the number of leading
        spaces on each line.
        """
        head=[]
        if len(self._head)>0:
            head+=[ indent+'%s\n'% line for values in self._head.values() for line in values ]

        if len(self._css)>0:
            head.append(indent+'<style type="text/css">\n')
            head+=[ indent+'  %s\n'%line for values in self._css.values() for line in values ]
            head.append(indent+'</style>\n')

        if len(self._javascript)>0:
            head.append(indent+'<script type="text/javascript">\n')
            head+=[ indent+'  %s\n'%line for values in self._javascript.values() for line in values ]
            head.append(indent+'</script>\n')

        return ''.join(head)


    def _html_document_ready(self, indent):
        doc_ready=indent+'<script type="text/javascript">\n'
        doc_ready+=indent+'  $(document).ready(function() {\n'
        doc_ready+=''.join(indent+'    %s\n'% line for values in self._document_ready.values() for line in values )
        doc_ready+=indent+'  })\n'
        doc_ready+=indent+'</script>\n'
        return doc_ready


    def _html_toc(self, indent):
        r"""
        Return a string for the table of contents sidebar.
        """

        # only return a toc if tere are sections
        if self._sections==[]: return ''

        toc =indent+'<div class="sphinxsidebar">\n'
        toc+=indent+'  <div class="sphinxsidebarwrapper">\n'
        toc+=indent+'    <h3><a href="#">Table Of Contents</a></h3>\n'
        toc+=indent+'    <ul>\n'
        toc+=indent+''.join('      <li><a class="reference internal" href="#%s">%s</a></li>\n'%(
          section['href'],section['title']) for section in self._sections if section['title']!=''
        )
        toc+=indent+'    </ul>\n'
        toc+=indent+'  </div>   <!-- close sphinxsidebarwrapper -->\n'
        toc+=indent+'</div>   <!-- close sphinxsidebar -->\n'
        return toc


    def _html_sections(self, indent):
        r"""
        Return a string for displaying all of the sections on the page as stored
        in ``self._section.
        """
        sections=''
        for section in self._sections:
            sections+='\n'+indent+'<div class="section" id="%s">\n' % section['href']
            sections+=indent+'  <h1>%s<a class="headerlink" href="#%s" title="Permalink to this headline">&para;</a></h1>\n' % (
                    section['title'], section['href'] )
            sections+=''.join(indent+'  %s\n'%line for line in section['text'])
            sections+=indent+'</div>    <!-- close section %s -->\n' % section['href']
        return sections


    def _html_(self):
        """
        Return an html string for the body of the web page ``self``.

        TODO:  Need to rewrite so that the substitions are only done via a call
        to :meth:`_html_` as currently changing the``export`` option after the
        table has been initialised will not have any effect.
        """
        return html_page_template.format(page=self,year=strftime('%Y',localtime()))


    def __html_file(self):
        r"""
        Return a (temporary) file name (and throw away the associated handler).
        This will be used to write the html file for. We make this a lazy
        attribute so that it is only necessary to reload this file when the
        underlying html document changes.
        """
        if self._html_file is None:
            # generate a filename for a temporary system file
            import tempfile
            tmpfile=tempfile.NamedTemporaryFile(suffix='.html',delete=False)
            tmpfile.close()
            self._html_file=tmpfile.name
        return self._html_file


    def __delitem__(self):
        r"""
        (Attempt to) remove the temporary file ``self._html_file`` when ``self``
        is deleted. May causse havoc if the web page is still being used?
        """
        import os.remove
        try:
            remove(self._html_file)
        except OSError:
            pass


    def html_display(self):
        r"""
        Display the page in a web browser.
        """
        import webbrowser
        html=open(self.__html_file(),'w')
        html.write(self._html_())
        html.close()
        webbrowser.open("file:///" + self._html_file, new=0)


    def add_section(self, section, heading=None, section_number=None):
        r"""
        Adds a new section to the page with title ``heading`` and text ``text``.
        All of the section titles on the page must be unique, so if ``heading``
        is equal to a pre-existing title then this section will overwrite the
        previous one. If ``section_number`` is not ``None`` then this section
        will be inserted into this position in the page - and any previous
        section of the same heading being deleted.
        """
        if heading is None and isinstance(section, HTMLElement):
            heading='%(title)s' % section

        if heading is None or heading=='':
            raise ValueError('every section needs a heading\n')

        headings=[sec['title'] for sec in self._sections]
        if heading in headings:
            sec=headings.index(heading)
            if section_number is not None and sec!=section_number:
                self._sections.pop(sec)  # remove old section
                sec=section_number
                self._sections.insert(sec,None) # make space for the new section
        else:
            sec=section_number if section_number is not None else len(self._sections)
            self._sections.insert(sec, {}) # make space for the new section

        # this section will become self._sections[sec]
        anchor=heading.lower().replace(' ','_')  # use for href's
        if isinstance(section, HTMLElement):
            self._sections[sec]=dict(title=heading, href=anchor, text=section.body(self._url_include))
            self.add_meta_data(section)
        else:
            self._sections[sec]=dict(title=heading, href=anchor, text=section.split('\n'))


    def table(self):
        r"""
        Returns a :class:`HTMLTable` which can be used to construct an HTML
        table which can then be added back into the page. This is really a
        convenient shorthand for adding tables to the HTML page rather than
        an actual method of the class:`HTMLDisplay` ``self``.
        """
        return HTMLTable()

#-----------------------------------------------------------------------
#Various constant strings/templates which specify parts of the web pages
#-----------------------------------------------------------------------

css_for_scrollable_table='''div.dataTables_scroll { clear: both; }
  table.dataTable {
    table-layout: fixed;
    width:100%;
  }
  table.dataTable th, table.dataTable td {
    text-align: center;
    padding: 0;
  }
  table.dataTable tr {
    height: 1.5em;
    padding: 1px;
  }
'''

javascript_for_scrollable_table='''$(".dataTable").each(function() {
    var Id = $(this).get(0).id;

    var oTable=$("#" + Id).dataTable({
        "sScrollY": ( 0.9 * $(window).height() ),
        "sScrollX": "100%", //columns will be mis-aligned if using sSCrollX
        "sScrollXInner": "150%",
        "bScrollCollapse": true,
        "aaSorting": [], // Disable initial sort
        "bDeferRender": true,
        "bPaginate": false,
        "bRetrieve": false,
        "bFilter": false,
        "bInfo": false,
        "bSort": false,
        "sPaginationType": "full_numbers",
        "sHeightMatch": "none",
        "aoColumnDefs" :[
          { aTargets: ["_all"], bSortable: false },
          { aTargets: [1], bSortable: true },
          { "sContentPadding": "mmm", "aTargets": ['_all'] }
        ]
    });

    new FixedColumns( oTable, {
        "iLeftColumns": 1 // one fixed column on the left
    } );
    oTable.fnAdjustColumnSizing();
});
'''


generic_meta_data='''<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="Generator" content="Sage HTML Display version {version}" />
<meta name="Copyright"
      content="Distributed under the terms of the GNU General Public License (GPLh) http://www.gnu.org/licenses/sage_html" />
<link rel="stylesheet" href="{url_include}_static/sage.css" type="text/css" />
<link rel="stylesheet" href="{url_include}_static/pygments.css" type="text/css" />
<link rel="shortcut icon" href="{url_include}_static/favicon.ico"/>
<link rel="icon" href="{url_include}_static/sageicon.png" type="image/x-icon" />
<script type="text/javascript" src="{url_include}_static/jquery.js"></script>
<script type="text/javascript" src="{url_include}_static/underscore.js"></script>
<script type="text/javascript" src="{url_include}_static/doctools.js"></script>
<script type="text/javascript" src="{url_include}_static/MathJax.js?config=TeX-AMS_HTML-full,../mathjax_sage.js"></script>
<script type="text/javascript">
  var DOCUMENTATION_OPTIONS = {{
      URL_ROOT: '{url_include}',
      VERSION: '{version}',
      COLLAPSE_INDEX: true,
      FILE_SUFFIX: '.html',
      HAS_SOURCE: false
  }};
  // start with the contents pane hidden as it will usually be redundant
  bar = jQuery('div.sphinxsidebar');
  bar.hide();
</script>'''


# This template determines the html string for an :class:`HTMLDisplay` web page.
# Via the :meth:`__getitem__` method defined in :class:`HTMLElement` this
# template triggers calls to various methods in :class:`HTMLDisplay`.
# Specifically, below ``page[method|indent]`` triggers page._method(width)
html_page_template='''<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <title>{page[title]}</title>
{page[html_meta|2]}
  </head>

  <body>

{page[html_navigation|2]}

    <div class="document">
      <div class="documentwrapper">
        <div class="bodywrapper">
          <div class="body">
{page[html_sections|6]}

          </div>   <!-- close body -->
        </div>   <!-- close bodywrapper -->
      </div>   <!-- close documentwrapper -->

{page[html_toc|3]}

      <div class="clearer"></div>

    </div>   <!-- close document -->

{page[html_navigation|2]}
    <div class="footer">
        &copy; Copyright 2005--{year}, The Sage Development Team.
    </div>
{page[html_document_ready|2]}
  </body>
</html>
'''


sphinx_sidebar_toggle='''/*global jQuery, window */
/* Sphinx sidebar toggle.  Putting this code at the end of the body
* enables the toggle for the live, static, and offline docs.  Note:
* sage.misc.html.math_parse() eats jQuery's dollar-sign shortcut. */
var bar, bod, bg, fg, key, tog, wid_old, wid_new, resize, get_state, set_state;
bod =jQuery('div.bodywrapper');
bar =jQuery('div.sphinxsidebar');
tog =jQuery('<div class="sphinxsidebartoggle"></div>');

/* Delayed resize helper.  Not perfect but good enough. */
resize = function () {
    setTimeout(function () {
        tog.height(bod.height());
    }, 100);
};
jQuery(window).resize(function () {
    resize();
});

/* Setup and add the toggle. See Sphinx v0.5.1 default.css. */
fg =jQuery('div.sphinxsidebar p a').css('color') || 'rgb(152, 219, 204)';
bg =jQuery('div.document').css('background-color') || 'rgb(28, 78, 99)';
wid_old = '230px';
wid_new = '5px';
tog.css('background-color', bg)
    .css('border-width', '0px')
    .css('border-right', wid_new + ' ridge ' + bg)
    .css('cursor', 'pointer')
    .css('position', 'absolute')
    .css('left', '-' + wid_new)
    .css('top', '0px')
    .css('width', wid_new);
bod.css('position', 'relative');
bod.prepend(tog);
resize()

/* Cookie helpers. */
key = 'sphinxsidebar=';
set_state = function (s) {
    var date = new Date();
    /* Expiry in 7 days. */
    date.setTime(date.getTime() + (7 * 24 * 3600 * 1000));
    document.cookie = key + encodeURIComponent(s) + '; expires=' +
        date.toUTCString() + '; path=/';
};
get_state = function () {
    var i, c, crumbs = document.cookie.split(';');
    for (i = 0; i < crumbs.length; i += 1) {
        c = crumbs[i].replace(/^\s+/, '');
        if (c.indexOf(key) === 0) {
            return decodeURIComponent(c.substring(key.length, c.length));
        }
    }
    return null;
};

/* Event handlers. */
tog.mouseover(function (ev) {
    tog.css('border-right-color', fg);
}).mouseout(function (ev) {
    tog.css('border-right-color', bg);
}).click(function (ev) {
    if (bod.hasClass('wide')) {
        bod.removeClass('wide');
        bod.css('margin-left', wid_old);
        bar.css('width', wid_old);
        bar.show();
        set_state('visible');
    } else {
        set_state('hidden');
        bar.hide();
        bar.css('width', '0px');
        bod.css('margin-left', wid_new);
        bod.addClass('wide');
    }
    resize();
});

/* Hide the normally visible sidebar? */
if (get_state() === 'hidden') {
    tog.trigger('click');
} else {
    set_state('visible');
}
'''
