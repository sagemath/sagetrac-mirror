# -*- encoding: utf-8 -*-
"""
Wrapper around Commandline Tools

The instances of this class are in separate files inside the
``sage.interfaces.cmdline`` package.
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
import subprocess
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.sage_ostools import have_program


class ToolMissingException(Exception):

    def __init__(self, *tools):
        """
        Exception to be thrown if command-line tools are not installed

        INPUT:

        - ``*tools`` -- any number of :class:`Tool` instances. The
          tools that are not installed.

        EXAMPLES::

            sage: from sage.interfaces.cmdline import ToolMissingException
            sage: from sage.interfaces.cmdline.posix import cat
            sage: raise ToolMissingException(cat)
            Traceback (most recent call last):
            ...
            ToolMissingException: One or more of the following required tools is not installed:
            * cat: Part of GNU coreutils
            <BLANKLINE>
            You might be able to install it with:
            * sudo yum install coreutils
        """
        self._tools = tools
        msg = ['One or more of the following required tools is not installed:']
        for tool in self._tools:
            msg.append('* {0}: {1}'.format(tool._name, tool._help))
        msg.append('')
        is_not_None = lambda x: x is not None
        rpm = [tool._package.get('rpm', None) for tool in self._tools]
        rpm = ' '.join(filter(is_not_None, rpm))
        deb = [tool._package.get('deb', None) for tool in self._tools]
        deb = ' '.join(filter(is_not_None, deb))
        if rpm or deb:
            msg.append('You might be able to install it with:')
        if rpm:
            msg.append('* sudo yum install {0}'.format(rpm))
        if deb:
            msg.append('* sudo apt-get install {0}'.format(deb))
        super(ToolMissingException, self).__init__('\n'.join(msg))
        
    def convert_to_warning(self):
        """
        Display warning with the same exception message

        EXAMPLES::

            sage: from sage.interfaces.cmdline.tool import ToolMissingException
            sage: from sage.interfaces.cmdline.posix import cat
            sage: exc = ToolMissingException(cat)
            sage: exc.convert_to_warning()
            doctest:...: UserWarning: One or more of the following required tools is not installed:
            * cat: Part of GNU coreutils
            <BLANKLINE>
            You might be able to install it with:
            * sudo ...
        """
        import warnings
        warnings.warn(self.message)
    
        
class CalledToolError(subprocess.CalledProcessError):

    def __str__(self):
        """
        Return a string representation

        OUTPUT:

        String.
        
        EXAMPLES::

            sage: from sage.interfaces.cmdline.tool import CalledToolError
            sage: exc = CalledToolError(123, 'foo', output='stdout is here')
            sage: print str(exc)
            Command 'foo' returned non-zero exit status 123
            ----[ output follows ]---------------------------------------------------------
            stdout is here
            -------------------------------------------------------------------------------
        """
        s = [super(CalledToolError, self).__str__()]
        s.append('----[ output follows ]' + '-' * 57)
        s.append(self.output)
        s.append('-' * 79)
        return '\n'.join(s)
        
        
class Tool(SageObject):

    _instances = set()

    def __init__(self, name, args, help, test, rpm=None, deb=None):
        """
        Command line tool

        All command line tools are supposed to work as filters where
        you pipe in the pdf file (or other input), and the tool writes
        the output to stdout.

        INPUT:

        - ``name`` -- string. The name of the tool. Need not be
          unique.

        - ``args`` -- list of strings. The name of the executable and
          the necessary command line arguments.

        - ``help`` -- string. Human-readable hint instructions for
          where the tool can be found. Should be short and not refer
          to distribution specifics.

        - ``test`` -- a function that can test the proper function of
          the tool. Must take the tool as the only argument. Will only
          be called if the tool is installed.

        - ``rpm`` -- optional string. The rpm package name in
          Redhat-derived distributions.

        - ``deb`` -- optional string. The deb package name in
          Debian-derived distributions.

        EXAMPLES::

            sage: from sage.interfaces.cmdline.posix import cat
            sage: type(cat)
            <class 'sage.interfaces.cmdline.tool.Tool'>
            sage: cat.is_installed()
            True
        """
        self._name = name
        self._args = tuple(args)
        self._help = help
        self._test = test
        self._package = dict(rpm=rpm, deb=deb)
        assert not isinstance(Tool._instances, frozenset), \
            'tools can only be defined within Sage'
        Tool._instances.add(self)

    def __hash__(self):
        """
        Make tool instances hashable

        OUTPUT:

        Integer. The hash of the tool instance.

        EXAMPLES::

            sage: from sage.interfaces.cmdline.posix import cat
            sage: cat.__hash__()   # random output
            -4675830580846295185
        """
        return hash(self._name) ^ hash(self._args)

    def __eq__(self, other):
        """
        Define Equality of tool instances

        INPUT:

        - ``other`` -- another tool instance.
        
        OUTPUT:

        Boolean. Whether ``other`` equals this tool.

        EXAMPLES::

            sage: from sage.interfaces.cmdline.posix import cat, cp
            sage: cat.__eq__(cat)
            True
            sage: cat.__eq__(cp)
            False
            sage: cat == cp
            False
            sage: cat == None
            False
        """
        if type(self) != type(other):
            return False
        return (self._name == other._name) and (self._args == other._args)

    @classmethod
    def _disallow_further_instances(cls):
        """
        Called from ``__init__.py`` to disallow further instances

        EXAMPLES::

            sage: from sage.interfaces.cmdline import Tool
            sage: Tool._disallow_further_instances()
            sage: Tool('foo', ['foo'], 'help', 'test')
            Traceback (most recent call last):
            ...
            AssertionError: tools can only be defined within Sage
        """
        cls._instances = frozenset(cls._instances)
    
    @classmethod
    def all(cls):
        """
        Return a list of all PDF tools

        OUTPUT:

        List of :class:`Tool` instances.

        EXAMPLES::

            sage: from sage.interfaces.cmdline import Tool
            sage: Tool.all()
            [cat,
             convert dvi->pdf,
             convert pdf->png,
             convert pdf->svg,
             cp,
             dvipdf,
             gcc,
             ghostscript pdf->png,
             latex,
             pdf2svg,
             pdfcrop,
             pdfjam,
             pdflatex,
             pdftocairo pdf->svg,
             xelatex]
        """
        return sorted(cls._instances)
    
    @classmethod
    def _test_all(cls):
        """
        Test the different tools

        EXAMPLES::

            sage: from sage.interfaces.cmdline import Tool
            sage: Tool._test_all()
            testing cat
            testing convert dvi->pdf
            testing convert pdf->png
            testing convert pdf->svg
            testing cp
            testing dvipdf
            testing gcc
            testing ghostscript pdf->png
            testing latex
            testing pdf2svg
            testing pdfcrop
            testing pdfjam
            testing pdflatex
            testing pdftocairo pdf->svg
            testing xelatex
        """
        for tool in sorted(cls.all()):
            print('testing {0}'.format(tool))
            if tool.is_installed():
                tool._test(tool)
            
    def is_installed(self):
        """
        Return whether the tool is installed
    
        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.interfaces.cmdline.posix import cat
            sage: type(cat)
            <class 'sage.interfaces.cmdline.tool.Tool'>
            sage: cat.is_installed()
            True
        """
        return have_program(self._args[0])

    @classmethod
    def require(cls, *tools):
        """
        Check that at least one of the tools is installed

        INPUT:
        
        - ``*tools`` -- :class:`Tool`. The command-line tools that are
          required.

        OUTPUT:

        :class:`Tool`. The first tool listed in the arguments that is
        installed. Raises :class:`ToolMissingException` if none of them is
        installed.

        EXAMPLES::

            sage: from sage.interfaces.cmdline import Tool
            sage: from sage.interfaces.cmdline.posix import cat, cp
            sage: Tool.require(cat, cp)
            cat
        """
        for tool in tools:
            if tool.is_installed():
                return tool
        raise ToolMissingException(*tools)
                       
    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.interfaces.cmdline.posix import cat
            sage: cat._repr_()
            'cat'
        """
        return self._name
    
    def __lt__(self, other):
        """
        Make tools sortable

        INPUT:

        - ``other`` -- another tool.

        OUTPUT:

        Boolean. 

        EXAMPLES::

            sage: from sage.interfaces.cmdline.posix import cp, cat
            sage: cat.__lt__(cp)
            True
            sage: cat < cp
            True
        """
        return self._name < other._name

    def __call__(self, stdin, **kwds):
        """
        Pipe data through the tool
        
        INPUT:

        - ``stdin`` -- bytes (string). Data to pipe into the command
          line tool.

        OUTPUT:

        Bytes (string). The output of the command line tool.

        EXAMPLES::

            sage: from sage.interfaces.cmdline.posix import cat
            sage: cat.__call__('test')
            'test'
            sage: cat('test')
            'test'
        """
        args = tuple(arg.format(**kwds) for arg in self._args)
        process = subprocess.Popen(
            ('sage-native-execute',) + args,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate(stdin)
        retcode = process.poll()
        if retcode:
            raise CalledToolError(
                retcode, ' '.join(args), output=stdout+stderr)
        return stdout


class ToolAdapter(Tool):

    def __init__(self, *args, **kwds):
        """
        Adapter for command line tools that do not accept stdin

        The stdout / stderr of such tools is generally useless chatter
        and ignored.

        INPUT:

        - ``file_ext`` -- pair of strings (optional keyword
          argument). The file extension for input and output file.

        All other argument are the same as in :meth:`Tool.__init__`.

        EXAMPLES::

            sage: from sage.interfaces.cmdline.posix import cp
            sage: type(cp)
            <class 'sage.interfaces.cmdline.tool.ToolAdapter'>
        """
        self._input_ext, self._output_ext = kwds.pop('file_ext', ['', ''])
        super(ToolAdapter, self).__init__(*args, **kwds)
    
    def __call__(self, stdin, **kwds):
        """
        Simulate piping data through the tool
        
        INPUT:

        - ``stdin`` -- bytes (string). Input data.

        OUTPUT:

        Bytes (string). The output of the command line tool.

        EXAMPLES::

            sage: from sage.interfaces.cmdline.posix import cp
            sage: cp.__call__('test')
            'test'
            sage: cp('test')
            'test'
        """
        from sage.misc.temporary_file import tmp_filename
        input_file = tmp_filename(ext=self._input_ext)
        with open(input_file, 'wb') as f:
            f.write(stdin)
        output_file = tmp_filename(ext=self._output_ext)
        kwds = dict(kwds)
        kwds['input'] = input_file
        kwds['output'] = output_file
        args = tuple(arg.format(**kwds) for arg in self._args)
        process = subprocess.Popen(
            ('sage-native-execute',) + args,
            stdin=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
        output, _ = process.communicate()
        retcode = process.poll()
        if retcode:
            raise CalledToolError(retcode, ' '.join(args), output=output)
        with open(output_file, 'rb') as f:
            return f.read()



        

class CompilerTool(Tool):

    def __init__(self, *args, **kwds):
        """
        Adapter for command line compilers

        The assumption is that the tool produces multiple output
        files, however only one that we are intersted in. It is run in
        a temporary directory.

        INPUT:

        - ``filename`` -- pair of strings (required keyword
          argument). The file name for input and output file.

        All other argument are the same as in :meth:`Tool.__init__`.

        EXAMPLES::

            sage: from sage.interfaces.cmdline.gcc import gcc
            sage: type(gcc)
            <class 'sage.interfaces.cmdline.tool.CompilerTool'>
        """
        self._input_file, self._output_file = kwds.pop('filename')
        super(CompilerTool, self).__init__(*args, **kwds)

    def __call__(self, stdin, **kwds):
        r"""
        Simulate piping data through a compiler
        
        INPUT:

        - ``stdin`` -- bytes (string). Input data.

        OUTPUT:

        Bytes (string). The output of the command line tool.

        EXAMPLES::

            sage: src = 'int main(int argc, char **argv) { return 0; }'
            sage: from sage.interfaces.cmdline.gcc import gcc
            sage: gcc(src)    # random output
            '\x7fELF...'
        """
        from sage.misc.temporary_file import tmp_dir
        directory = tmp_dir()
        input_file = os.path.join(directory, self._input_file)
        output_file = os.path.join(directory, self._output_file)
        with open(input_file, 'wb') as f:
            f.write(stdin)
        args = tuple(arg.format(**kwds) for arg in self._args)
        process = subprocess.Popen(
            ('sage-native-execute',) + args,
            stdin=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            cwd=directory,
        )
        output, _ = process.communicate()
        retcode = process.poll()
        import shutil
        try:
            if retcode:
                raise CalledToolError(retcode, ' '.join(args), output=output)
            with open(output_file, 'rb') as f:
                return f.read()
        finally:
            shutil.rmtree(directory)
    
