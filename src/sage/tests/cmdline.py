r"""
This file contains some tests that Sage command line options actually
do something.

We test the following command line options:

test.py
/path/to/test.py
test.sage
/path/to/test.sage
test.spyx
/path/to/test.spyx
--advanced
-c
--cython
--dev
--ecl
--experimental
--fixdoctests
--gap
--gdb
--gp
-h
--help
--info
--ipython
--kash
--lisp
--maxima
--min
--mwrank
--optional
--preparse
--python
-q
--R
--root
--rst2txt
--rst2sws
--sh
--singular
--sqlite3
--standard
--startuptime
-t
-v
--zzfoobar (illegal option)


AUTHORS:

- Jeroen Demeyer (2010-11-20): initial version (#10300)

"""

from subprocess import *
import os, select


def test_executable(args, input="", timeout=100.0, **kwds):
    r"""
    Run the program defined by ``args`` using the string ``input`` on
    the standard input.

    INPUT:

    - ``args`` -- a list of program arguments, the first being the
      executable.

    - ``input`` -- a string serving as standard input.  Usually, this
      should end with a newline.

    - ``timeout`` -- if the program produces no output for ``timeout``
      seconds, a RuntimeError is raised.

    - ``**kwds`` -- Additional keyword arguments passed to the
      :class:`Popen` constructor.

    OUTPUT: a tuple ``(out, err, ret)`` with the standard output,
    standard error and exitcode of the program run.


    EXAMPLES::

        sage: from sage.tests.cmdline import test_executable
        sage: (out, err, ret) = test_executable(["cat"], "Hello World!")
        sage: out
        'Hello World!'
        sage: err
        ''
        sage: ret
        0

    We test the timeout option::

        sage: (out, err, ret) = test_executable(["sleep", "1"], timeout=0.1)
        Traceback (most recent call last):
        ...
        RuntimeError: timeout in test_executable()

    TESTS:

    Run Sage itself with various options::

        sage: (out, err, ret) = test_executable(["sage"])
        sage: out.find(version()) >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage"], "3^33\n")
        sage: out.find(version()) >= 0
        True
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "-q"], "3^33\n")
        sage: out.find(version()) >= 0
        False
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "-c", "print 3^33"])
        sage: print out
        5559060566555523
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--min", "-c", "print 3^33"])
        sage: print out
        5559060566555523
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--startuptime"])
        sage: out.find("Slowest module import") >= 0
        True
        sage: err
        ''
        sage: ret
        0

    Test help::

        sage: (out, err, ret) = test_executable(["sage", "-h"])
        sage: out.find("Optional arguments:") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--help"])
        sage: out.find("Optional arguments:") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--advanced"])
        sage: out.find("search through the Sage documentation") >= 0
        True
        sage: err
        ''
        sage: ret
        0

    Basic information about the Sage installation::

        sage: (out, err, ret) = test_executable(["sage", "-v"])
        sage: out.find(version()) >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--root"])
        sage: len(out) >= 2   # at least one character + newline
        True
        sage: err
        ''
        sage: ret
        0

    Test ``sage --info [packages]`` and the equivalent
    ``sage -p --info --info [packages]`` (the doubling of ``--info``
    is intentional, that option should be idempotent)::

        sage: out, err, ret = test_executable(["sage", "--info", "sqlite"])
        sage: print out
        Found local metadata for sqlite-...
        = SQLite =
        ...
        SQLite is a software library that implements a self-contained,
        serverless, zero-configuration, transactional SQL database engine.
        ...
        sage: err
        ''
        sage: ret
        0

        sage: out, err, ret = test_executable(["sage", "-p", "--info", "--info", "sqlite"])
        sage: print out
        Found local metadata for sqlite-...
        = SQLite =
        ...
        SQLite is a software library that implements a self-contained,
        serverless, zero-configuration, transactional SQL database engine.
        ...
        sage: err
        ''
        sage: ret
        0

    Test ``sage-run`` on a Python file, both with an absolute and with a relative path::

        sage: dir = tmp_dir(); name = 'python_test_file.py'
        sage: fullname = os.path.join(dir, name)
        sage: F = open(fullname, 'w')
        sage: F.write("print 3^33\n")
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", fullname])
        sage: print out
        34
        sage: err
        ''
        sage: ret
        0
        sage: (out, err, ret) = test_executable(["sage", name], cwd=dir)
        sage: print out
        34
        sage: err
        ''
        sage: ret
        0

    The same as above, but now with a ``.sage`` file.  This indirectly
    also tests the preparser::

        sage: dir = tmp_dir(); name = 'sage_test_file.sage'
        sage: fullname = os.path.join(dir, name)
        sage: F = open(fullname, 'w')
        sage: F.write("k.<a> = GF(5^3); print a^124\n")
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", fullname])
        sage: print out
        1
        sage: err
        ''
        sage: ret
        0
        sage: (out, err, ret) = test_executable(["sage", name], cwd=dir)
        sage: print out
        1
        sage: err
        ''
        sage: ret
        0

    Test running a ``.spyx`` file::

        sage: dir = tmp_dir(); name = 'sage_test_file.spyx'
        sage: fullname = os.path.join(dir, name)
        sage: F = open(fullname, 'w')
        sage: F.write("from sage.rings.integer cimport Integer\ncdef long i, s = 0\nsig_on()\nfor i in range(1000): s += i\nsig_off()\nprint Integer(s)")
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", fullname])
        sage: print out
        499500
        sage: err
        'Compiling ...spyx...'
        sage: ret
        0
        sage: (out, err, ret) = test_executable(["sage", name], cwd=dir)
        sage: print out
        499500
        sage: err
        'Compiling ...spyx...'
        sage: ret
        0

    Testing ``sage --preparse FILE`` and ``sage -t FILE``.  First create
    a file and preparse it::

        sage: s = "'''\nThis is a test file.\n'''\ndef my_add(a,b):\n    '''\n    Add a to b.\n\n        EXAMPLES::\n\n            sage: my_add(2,2)\n            4\n        '''\n    return a + b\n"
        sage: script = os.path.join(tmp_dir(), 'my_script.sage')
        sage: script_py = script + '.py'
        sage: F = open(script, 'w')
        sage: F.write(s)
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", "--preparse", script])
        sage: ret
        0
        sage: os.path.isfile(script_py)
        True

    Now test my_script.sage and the preparsed version my_script.sage.py::

        sage: (out, err, ret) = test_executable(["sage", "-t", script])
        sage: ret
        0
        sage: out.find("All tests passed!") >= 0
        True
        sage: (out, err, ret) = test_executable(["sage", "-t", script_py])
        sage: ret
        0
        sage: out.find("All tests passed!") >= 0
        True

    Now for a file which should fail tests::

        sage: s = s.replace('4', '5') # (2+2 != 5)
        sage: F = open(script, 'w')
        sage: F.write(s)
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", "-t", script])
        sage: ret
        1
        sage: out.find("1 item had failures:") >= 0
        True

    Test ``sage -t --debug -p 2`` on a ReST file, the ``-p 2`` should
    be ignored. In Pdb, we run the ``help`` command::

        sage: s = "::\n\n    sage: assert True == False\n    sage: 2 + 2\n    5"
        sage: script = tmp_filename(ext='.rst')
        sage: F = open(script, 'w')
        sage: F.write(s)
        sage: F.close()
        sage: (out, err, ret) = test_executable([
        ....:     "sage", "-t", "--debug", "-p", "2", "--warn-long", "0", script], "help")
        sage: print out
        Debugging requires single-threaded operation, setting number of threads to 1.
        Running doctests with ID...
        Doctesting 1 file.
        sage -t ...
        **********************************************************************
        File "...", line 3, in ...
        Failed example:
            assert True == False
        Exception raised:
            Traceback (most recent call last):
            ...
            AssertionError
        > <doctest ...>(1)<module>()
        -> assert True == False
        (Pdb)
        Documented commands (type help <topic>):
        ========================================
        ...
        **********************************************************************
        File "...", line 4, in ...
        Failed example:
            2 + 2
        Expected:
            5
        Got:
            4
        **********************************************************************
        Previously executed commands:
            s...: assert True == False
        debug:
        <BLANKLINE>
        Returning to doctests...
        **********************************************************************
        1 item had failures:
           2 of   3 in ...
            [2 tests, 2 failures, ...]
        ...
        sage: ret
        1

    Check that Sage refuses to run doctests from a directory whose
    permissions are too loose.  We create a world-writable directory
    inside a safe temporary directory to test this::

        sage: d = os.path.join(tmp_dir(), "test")
        sage: os.mkdir(d)
        sage: os.chmod(d, 0o777)
        sage: (out, err, ret) = test_executable(["sage", "-t", "nonexisting.py"], cwd=d)
        sage: print err
        Traceback (most recent call last):
        ...
        RuntimeError: refusing to run doctests...
        sage: (out, err, ret) = test_executable(["sage", "-tp", "1", "nonexisting.py"], cwd=d)
        sage: print err
        Traceback (most recent call last):
        ...
        RuntimeError: refusing to run doctests...

    Now run a test for the fixdoctests script and, in particular, check that the
    issues raised in :trac:`10589` are fixed. We have to go to slightly silly
    lengths to doctest the output.::

        sage: test='r\"\"\"Add a doc-test for the fixdoctest command line option and, in particular, check that\n:trac:`10589` is fixed.\n\nEXAMPLES::\n\n    sage: 1+1              # incorrect output\n    3\n    sage: m=matrix(ZZ,3)   # output when none is expected\n    [0 0 0]\n    [0 0 0]\n    [1 0 0]\n    sage: (2/3)*m          # no output when it is expected\n    sage: mu=PartitionTuple([[4,4],[3,3,2,1],[1,1]])   # output when none is expected\n    [4, 4, 3, 3, 2, 1, 1]\n    sage: mu.pp()          # uneven indentation\n    ****\n    ****\n    sage: PartitionTuples.global_options(convention="French")\n    sage: mu.pp()         # fix doctest with uneven indentation\n    sage: PartitionTuples.global_options.reset()\n\"\"\"\n'
        sage: test_file = os.path.join(tmp_dir(), 'test_file.py')
        sage: F = open(test_file, 'w')
        sage: F.write(test)
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", "--fixdoctests", test_file])
        sage: with open(test_file, 'r') as f:
        ....:     fixed_test = f.read()
        sage: import difflib
        sage: list(difflib.unified_diff(test.splitlines(), fixed_test.splitlines()))[2:-1]
        ['@@ -4,18 +4,23 @@\n',
         ' EXAMPLES::',
         ' ',
         '     sage: 1+1              # incorrect output',
         '-    3',
         '+    2',
         '     sage: m=matrix(ZZ,3)   # output when none is expected',
         '+    sage: (2/3)*m          # no output when it is expected',
         '     [0 0 0]',
         '     [0 0 0]',
         '-    [1 0 0]',
         '-    sage: (2/3)*m          # no output when it is expected',
         '+    [0 0 0]',
         '     sage: mu=PartitionTuple([[4,4],[3,3,2,1],[1,1]])   # output when none is expected',
         '-    [4, 4, 3, 3, 2, 1, 1]',
         '     sage: mu.pp()          # uneven indentation',
         '-    ****',
         '-    ****',
         '+       ****   ***   *',
         '+       ****   ***   *',
         '+              **',
         '+              *',
         '     sage: PartitionTuples.global_options(convention="French")',
         '     sage: mu.pp()         # fix doctest with uneven indentation',
         '+    *',
         '+    **',
         '+    ****   ***   *',
         '+    ****   ***   *',
         '     sage: PartitionTuples.global_options.reset()']

    Test external programs being called by Sage::

        sage: (out, err, ret) = test_executable(["sage", "--sh"], "echo Hello World\nexit 42\n")
        sage: out.find("Hello World\n") >= 0
        True
        sage: ret
        42

        sage: (out, err, ret) = test_executable(["sage", "--sh", "-c", "echo Hello World; exit 42"])
        sage: out.find("Hello World\n") >= 0
        True
        sage: ret
        42

        sage: (out, err, ret) = test_executable(["sage", "--ipython"], "\n3**33\n")
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--python"], "print 3^33\n")
        sage: out
        '34\n'
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--cython"])
        sage: print err
        Cython (http://cython.org) is a compiler for code written in the
        Cython language.  Cython is based on Pyrex by Greg Ewing.
        ...

        sage: def has_tty():
        ....:     try:
        ....:         os.open(os.ctermid(), os.O_RDONLY)
        ....:         return True
        ....:     except OSError:
        ....:         return False 
        sage: (out, err, ret) = test_executable(["sage", "--ecl"], "(* 12345 54321)\n")
        sage: out.find("Embeddable Common-Lisp") >= 0
        True
        sage: out.find("670592745") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--lisp"], "(* 12345 54321)\n")
        sage: out.find("Embeddable Common-Lisp") >= 0
        True
        sage: out.find("670592745") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--gap", "-q"], "Size(SymmetricGroup(5));\n")
        sage: out
        '120\n'
        sage: err.replace('gap: halving pool size.', '').strip()
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--gdb"], 'quit\n')  # optional - gdb
        sage: out.find('(gdb) ') >= 0  # optional - gdb
        True
        sage: ret  # optional - gdb
        0

        sage: (out, err, ret) = test_executable(["sage", "--kash", "-b", "3^33;\n"])  # optional - kash
        sage: out.find("5559060566555523") >= 0  # optional - kash
        True
        sage: err  # optional - kash
        ''
        sage: ret  # optional - kash
        0

        sage: (out, err, ret) = test_executable(["sage", "--mwrank", "-v0", "-q", "-o"], "0 0 1 -7 6 0 0 0 0 0\n")
        sage: out
        'Curve [0,0,1,-7,6] :\tRank = 3\n[[3],[[1,-1],[-2,3],[-7/4,25/8]]]\n\n\n'
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--singular"], "12345*54321;\n")
        sage: out.find("A Computer Algebra System for Polynomial Computations") >= 0
        True
        sage: out.find("670592745") >= 0
        True
        sage: err
        ''
        sage: ret
        0

    Check that ``sage-location`` did its job in making Python scripts
    relative.  We test it on the ``ipython`` script::

        sage: open(os.path.join(SAGE_LOCAL, "bin", "ipython")).readline()
        '#!/usr/bin/env python\n'

    Test GP using the ``-f`` option which prevents the reading of a ``.gprc``
    configuration file::

        sage: (out, err, ret) = test_executable(["sage", "--gp", "-f"], "3^33\nquit(42)\n")
        sage: out.find("PARI/GP") >= 0
        True
        sage: out.find("5559060566555523") >= 0
        True
        sage: err
        ''
        sage: ret
        42

    Some programs of which we check functionality using only ``--version``::

        sage: (out, err, ret) = test_executable(["sage", "--maxima", "--version"])
        sage: out.find("Maxima ") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--R", "--version"])
        sage: out.find("R version ") >= 0
        True
        sage: err
        ''
        sage: ret
        0

        sage: (out, err, ret) = test_executable(["sage", "--sqlite3", "--version"])
        sage: out.startswith("3.")
        True
        sage: err
        ''
        sage: ret
        0

    Check some things requiring an internet connection::

        sage: (out, err, ret) = test_executable(["sage", "--standard"])  # optional - internet
        sage: out.find("atlas") >= 0  # optional - internet
        True
        sage: err  # optional - internet
        ''
        sage: ret  # optional - internet
        0

        sage: (out, err, ret) = test_executable(["sage", "--optional"])  # optional - internet
        sage: out.find("database_cremona_ellcurve") >= 0  # optional - internet
        True
        sage: err  # optional - internet
        ''
        sage: ret  # optional - internet
        0

        sage: (out, err, ret) = test_executable(["sage", "--experimental"])  # optional - internet
        sage: out.find("macaulay2") >= 0  # optional - internet
        True
        sage: err  # optional - internet
        ''
        sage: ret  # optional - internet
        0

    Check an illegal command line option.  This outputs an error to stdout,
    but we allow stderr in case this changes in the future::

        sage: (out, err, ret) = test_executable(["sage", "--zzfoobar"])
        sage: (out+err).find("unknown option: --zzfoobar") >= 0
        True
        sage: ret > 0
        True

    Test ``sage --rst2txt file.rst`` on a ReST file::

        sage: s = "::\n\n    sage: 2^10\n    1024\n    sage: 2 + 2\n    4"
        sage: input = tmp_filename(ext='.rst')
        sage: F = open(input, 'w')
        sage: F.write(s)
        sage: F.close()
        sage: (out, err, ret) = test_executable(["sage", "--rst2txt", input])
        sage: print out
        {{{id=0|
        2^10
        ///
        1024
        }}}
        <BLANKLINE>
        {{{id=1|
        2 + 2
        ///
        4
        }}}
        sage: err
        ''
        sage: ret
        0

    Test ``sage --rst2txt file.rst file.txt`` on a ReST file::

        sage: s = "::\n\n    sage: 2^10\n    1024\n    sage: 2 + 2\n    4"
        sage: input = tmp_filename(ext='.rst')
        sage: output = tmp_filename(ext='.txt')
        sage: F = open(input, 'w')
        sage: F.write(s)
        sage: F.close()
        sage: test_executable(["sage", "--rst2txt", input, output])
        ('', '', 0)
        sage: print open(output, 'r').read()
        {{{id=0|
        2^10
        ///
        1024
        }}}
        <BLANKLINE>
        {{{id=1|
        2 + 2
        ///
        4
        }}}

    Test ``sage --rst2sws file.rst file.sws`` on a ReST file::

        sage: s = "Thetitle\n--------\n\n::\n\n    sage: 2^10\n    1024\n    sage: 2 + 2\n    4"
        sage: input = tmp_filename(ext='.rst')
        sage: output = tmp_filename(ext='.sws')
        sage: F = open(input, 'w')
        sage: F.write(s)
        sage: F.close()
        sage: test_executable(["sage", "--rst2sws", input, output])
        ('', '', 0)
        sage: import tarfile
        sage: f = tarfile.open(output, 'r')
        sage: print f.extractfile('sage_worksheet/worksheet.html').read()
        <h1 class="title">Thetitle</h1>
        <BLANKLINE>
        {{{id=0|
        2^10
        ///
        1024
        }}}
        <BLANKLINE>
        {{{id=1|
        2 + 2
        ///
        4
        }}}
        sage: print f.extractfile('sage_worksheet/worksheet.txt').read()
        Thetitle
        system:sage
        <BLANKLINE>
        <BLANKLINE>
        <h1 class="title">Thetitle</h1>
        <BLANKLINE>
        {{{id=0|
        2^10
        ///
        1024
        }}}
        <BLANKLINE>
        {{{id=1|
        2 + 2
        ///
        4
        }}}
    """
    pexpect_env = dict(os.environ)
    try:
        del pexpect_env["TERM"]
    except KeyError:
        pass
    p = Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, env=pexpect_env, **kwds)
    if input:
        p.stdin.write(input)
    p.stdin.close()
    fdout = p.stdout.fileno()
    fderr = p.stderr.fileno()
    out = ""; err = ""

    while True:
        # Try reading from fdout and fderr
        rfd = []
        if fdout: rfd.append(fdout)
        if fderr: rfd.append(fderr)
        if len(rfd) == 0: break
        rlist = select.select(rfd, [], [], timeout)[0]

        if len(rlist) == 0:
            # Timeout!
            p.terminate()
            raise RuntimeError("timeout in test_executable()")
        if fdout in rlist:
            s = os.read(fdout, 1024)
            if s == "": fdout = None   # EOF
            out += s
        if fderr in rlist:
            s = os.read(fderr, 1024)
            if s == "": fderr = None   # EOF
            err += s

    return (out, err, p.wait())
