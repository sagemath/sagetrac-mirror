r"""
This is Sage's version of the sphinx-build script

We redirect stdout to our own logger, and remove some unwanted chatter.
"""

import os, sys, re, sphinx

# override the fancy multi-line formatting
def term_width_line(text):
    return text + '\n'

sphinx.util.console.term_width_line = term_width_line


class SageSphinxLogger(object):
    r"""
    This implements the file object interface to serve as sys.stdout
    replacement.
    """
    ansi_color = re.compile(r'\x1b\[[0-9;]*m')
    ansi_reset = re.compile(r'\x1b\[39;49;00m')
    prefix_len = 9

    def __init__(self, stream, prefix):
        self._init_chatter()
        self._stream = stream
        self._color = stream.isatty()
        prefix = prefix[0:self.prefix_len]
        prefix = ('[{0:'+str(self.prefix_len)+'}]').format(prefix)
        self._is_stdout = (stream.fileno() == 1)
        self._is_stderr = (stream.fileno() == 2)
        if self._is_stdout:
            color = 'darkgreen'
        elif self._is_stderr:
            color = 'red'
        else:
            color = 'lightgray'
        self._prefix = sphinx.util.console.colorize(color, prefix)
        # When we see an error in the log, we store it here and raise it at the
        # end of the file (sometimes the lines following the error still
        # contain valuable information.)
        self._error = None

    def _init_chatter(self):
        # useless_chatter: regular expressions to be filtered from
        # Sphinx output.
        self.useless_chatter = (
            re.compile('^$'),
            re.compile('^Running Sphinx v'),
            re.compile('^loading intersphinx inventory from '),
            re.compile('^Compiling a sub-document'),
            re.compile('^updating environment: 0 added, 0 changed, 0 removed'),
            re.compile('^looking for now-outdated files... none found'),
            re.compile('^building \[.*\]: targets for 0 source files that are out of date'),
            re.compile('^loading pickled environment... done'),
            re.compile('^loading cross citations... done \([0-9]* citations\).'),
            re.compile('WARNING: favicon file \'favicon.ico\' does not exist'),
            re.compile('WARNING: html_static_path entry .* does not exist'),
            re.compile('WARNING: while setting up extension'),
            re.compile('WARNING: Any IDs not assiend for figure node'),
            re.compile('WARNING: .* is not referenced'),
            re.compile('WARNING: Build finished'),
            re.compile('language "hu" not supported'),
            )

        # replacements: pairs of regular expressions and their replacements,
        # to be applied to Sphinx output.
        self.replacements = [(re.compile('build succeeded, [0-9]+ warning[s]?.'),
                              'build succeeded.')]

        if 'inventory' in sys.argv:
            # When building the inventory, ignore warnings about missing
            # citations and the search index.
            self.useless_chatter += (
                re.compile('WARNING: citation not found:'),
                re.compile('WARNING: search index couldn\'t be loaded, but not all documents will be built: the index will be incomplete.')
                )

        # Regular expressions indicating a problem with docbuilding. Raise an
        # exception if any of these occur.
        self._error_patterns = (re.compile('Segmentation fault'),
                    re.compile('SEVERE'),
                    re.compile('ERROR'),
                    re.compile('^make.*Error'),
                    re.compile('Exception occurred'),
                    re.compile('Sphinx error'))

        # We want all warnings to actually be errors.
        # Exceptions:
        # - warnings upon building the LaTeX documentation
        # - undefined labels upon the first pass of the compilation: some
        #   cross links may legitimately not yet be resolvable at this point.
        if 'latex' not in sys.argv:
            if 'multidoc_first_pass=1' in sys.argv:
                # Catch all warnings except 'WARNING: undefined label'
                self._error_patterns += (re.compile('WARNING: (?!undefined label)'),)
            else:
                self._error_patterns += (re.compile('WARNING:'),)

    def _filter_out(self, line):
        if self._error and self._is_stdout:
            # swallow non-errors after an error occurred
            return True
        line = re.sub(self.ansi_color, '', line)
        for regex in self.useless_chatter:
            if regex.search(line) is not None:
                return True
        return False

    def _check_errors(self, line):
        r"""
        Search for errors in line.

        EXAMPLES::

            sage: from sys import stdout
            sage: from sage_setup.docbuild.sphinxbuild import SageSphinxLogger
            sage: logger = SageSphinxLogger(stdout, "doctesting")
            sage: logger._log_line("Segmentation fault!\n") # indirect doctest
            [doctestin] Segmentation fault!
            sage: logger.raise_errors()
            Traceback (most recent call last):
            ...
            OSError: [doctestin] Segmentation fault!

        """
        if self._error:
            return  # we already have found an error
        for regex in self._error_patterns:
            if regex.search(line) is not None:
                self._error = line
                return

    def _log_line(self, line):
        r"""
        Write ``line`` to the output stream with some mangling.

        EXAMPLES::

            sage: from sys import stdout
            sage: from sage_setup.docbuild.sphinxbuild import SageSphinxLogger
            sage: logger = SageSphinxLogger(stdout, "doctesting")
            sage: logger._log_line("building documentation…\n")
            [doctestin] building documentation…

        """
        if self._filter_out(line):
            return
        for (old, new) in self.replacements:
            line = old.sub(new, line)
        line = self._prefix + ' ' + line.rstrip() + '\n'
        if not self._color:
            line = self.ansi_color.sub('', line)
        self._stream.write(line)
        self._stream.flush()
        self._check_errors(line)

    def raise_errors(self):
        r"""
        Raise an exceptions if any errors have been found while parsing the
        Sphinx output.

        EXAMPLES::

            sage: from sys import stdout
            sage: from sage_setup.docbuild.sphinxbuild import SageSphinxLogger
            sage: logger = SageSphinxLogger(stdout, "doctesting")
            sage: logger._log_line("This is a SEVERE error\n")
            [doctestin] This is a SEVERE error
            sage: logger.raise_errors()
            Traceback (most recent call last):
            ...
            OSError: [doctestin] This is a SEVERE error

        """
        if self._error:
            raise OSError(self._error)

    _line_buffer = ''

    def _write(self, string):
        self._line_buffer += string
        lines = self._line_buffer.splitlines()
        for i, line in enumerate(lines):
            last = (i == len(lines)-1)
            if last and not self._line_buffer.endswith('\n'):
                self._line_buffer = line
                return
            self._log_line(line)
        self._line_buffer = ''


    # file object interface follows

    closed = False
    encoding = None
    mode = 'w'
    name = '<log>'
    newlines = None
    softspace = 0

    def isatty(self):
        return True

    def close(self):
        if self._line_buffer != '':
            self._log_line(self._line_buffer)
            self._line_buffer = ''

    def flush(self):
        self._stream.flush()

    def write(self, str):
        try:
            self._write(str)
        except OSError:
            raise
        except Exception:
            import traceback
            traceback.print_exc(file=self._stream)

    def writelines(self, sequence):
        for line in sequence:
            self.write(line)


def runsphinx():
    output_dir = sys.argv[-1]

    saved_stdout = sys.stdout
    saved_stderr = sys.stderr

    try:
        sys.stdout = SageSphinxLogger(sys.stdout, os.path.basename(output_dir))
        sys.stderr = SageSphinxLogger(sys.stderr, os.path.basename(output_dir))
        sphinx.cmdline.main(sys.argv)
        sys.stderr.raise_errors()
        sys.stdout.raise_errors()
    finally:
        sys.stdout = saved_stdout
        sys.stderr = saved_stderr
        sys.stdout.flush()
        sys.stderr.flush()

    if exception is not None:
        raise exception
