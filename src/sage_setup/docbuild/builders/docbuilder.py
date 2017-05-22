class DocBuilder(object):
    def __init__(self, name, lang='en'):
        """
        INPUT:

        - ``name`` - the name of a subdirectory in SAGE_DOC_SRC, such as
          'tutorial' or 'bordeaux_2008'

        - ``lang`` - (default "en") the language of the document.
        """
        doc = name.split(os.path.sep)

        if doc[0] in LANGUAGES:
            lang = doc[0]
            doc.pop(0)

        self.name = os.path.join(*doc)
        self.lang = lang
        self.dir = os.path.join(SAGE_DOC_SRC, self.lang, self.name)

    def _output_dir(self, type):
        """
        Returns the directory where the output of type ``type`` is stored.
        If the directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: from sage_setup.docbuild import DocBuilder
            sage: b = DocBuilder('tutorial')
            sage: b._output_dir('html')
            '.../html/en/tutorial'
        """
        d = os.path.join(SAGE_DOC, type, self.lang, self.name)
        sage_makedirs(d)
        return d

    def _doctrees_dir(self):
        """
        Returns the directory where the doctrees are stored.  If the
        directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: from sage_setup.docbuild import DocBuilder
            sage: b = DocBuilder('tutorial')
            sage: b._doctrees_dir()
            '.../doctrees/en/tutorial'
        """
        d = os.path.join(SAGE_DOC, 'doctrees', self.lang, self.name)
        sage_makedirs(d)
        return d

    def _output_formats(self):
        """
        Returns a list of the possible output formats.

        EXAMPLES::

            sage: from sage_setup.docbuild import DocBuilder
            sage: b = DocBuilder('tutorial')
            sage: b._output_formats()
            ['changes', 'html', 'htmlhelp', 'inventory', 'json', 'latex', 'linkcheck', 'pickle', 'web']

        """
        #Go through all the attributes of self and check to
        #see which ones have an 'is_output_format' attribute.  These
        #are the ones created with builder_helper.
        output_formats = []
        for attr in dir(self):
            if hasattr(getattr(self, attr), 'is_output_format'):
                output_formats.append(attr)
        output_formats.sort()
        return output_formats

    def pdf(self):
        """
        Builds the PDF files for this document.  This is done by first
        (re)-building the LaTeX output, going into that LaTeX
        directory, and running 'make all-pdf' (or for the special case of
        the ja docs, 'all-pdf-ja(ex,to run platex)' there.

        EXAMPLES::

            sage: from sage_setup.docbuild import DocBuilder
            sage: b = DocBuilder('tutorial')
            sage: b.pdf() #not tested
        """
        self.latex()
        tex_dir = self._output_dir('latex')
        pdf_dir = self._output_dir('pdf')
        make_target = "cd '%s' && $MAKE %s && mv -f *.pdf '%s'"
        error_message = "failed to run $MAKE %s in %s"
        MB_LANG = {'ja': 'all-pdf-ja'} # language name : the modified target

        # Replace the command for languages that require special processing
        if self.lang in MB_LANG:
            command = MB_LANG[self.lang]
        else:
            command = 'all-pdf'

        if subprocess.call(make_target%(tex_dir, command, pdf_dir), shell=True):
            raise RuntimeError(error_message%(command, tex_dir))
        logger.warning("Build finished.  The built documents can be found in %s", pdf_dir)

    def clean(self, *args):
        shutil.rmtree(self._doctrees_dir())
        output_formats = list(args) if args else self._output_formats()
        for format in output_formats:
            shutil.rmtree(self._output_dir(format), ignore_errors=True)

    html = builder_helper('html')
    pickle = builder_helper('pickle')
    web = pickle
    json = builder_helper('json')
    htmlhelp = builder_helper('htmlhelp')
    latex = builder_helper('latex')
    changes = builder_helper('changes')
    linkcheck = builder_helper('linkcheck')
    # import the customized builder for object.inv files
    inventory = builder_helper('inventory')
