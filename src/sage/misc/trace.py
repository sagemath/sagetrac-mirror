"""
Interactively tracing execution of a command
"""
from __future__ import print_function


def trace(code, preparse=True):
    r"""
    Evaluate Sage code using the interactive tracer and return the
    result. The string ``code`` must be a valid expression
    enclosed in quotes (no assignments - the result of the expression
    is returned). In the Sage notebook this just raises a
    NotImplementedException.

    INPUT:


    -  ``code`` - str

    -  ``preparse`` - bool (default: True); if True, run
       expression through the Sage preparser.


    REMARKS: This function is extremely powerful! For example, if you
    want to step through each line of execution of, e.g.,
    ``factor(100)``, type

    ::

        sage: trace("factor(100)")             # not tested

    then at the (Pdb) prompt type ``s`` (or ``step``), then press return
    over and over to step through every line of Python that is called
    in the course of the above computation. Type ``?`` at any time for
    help on how to use the debugger (e.g., ``l`` lists 11 lines around
    the current line; ``bt`` gives a back trace, etc.).

    Setting a break point: If you have some code in a file and would
    like to drop into the debugger at a given point, put the following
    code at that point in the file:

    ``import pdb; pdb.set_trace()``

    For an article on how to use the Python debugger, see
    http://www.onlamp.com/pub/a/python/2005/09/01/debugger.html

    TESTS:

    The only real way to test this is via pexpect spawning a
    sage subprocess that uses IPython.

    ::

        sage: import pexpect
        sage: s = pexpect.spawn('sage')
        sage: _ = s.sendline("trace('print(factor(10))'); print(3+97)")
        sage: _ = s.sendline("s"); _ = s.sendline("c");
        sage: _ = s.expect('100', timeout=90)

    Seeing the ipdb prompt and the 2 \* 5 in the output below is a
    strong indication that the trace command worked correctly.

    ::

        sage: print(s.before[s.before.find('--'):])
        --...
        ipdb> c
        2 * 5

    We test what happens in notebook embedded mode::

        sage: sage.plot.plot.EMBEDDED_MODE = True
        sage: trace('print(factor(10))')
        Traceback (most recent call last):
        ...
        NotImplementedError: the trace command is not implemented in the Sage notebook; you must use the command line.
    """
    from sage.plot.plot import EMBEDDED_MODE
    if EMBEDDED_MODE:
        raise NotImplementedError("the trace command is not implemented in the Sage notebook; you must use the command line.")

    from IPython.core.debugger import Pdb
    pdb = Pdb()

    try:
        ipython = get_ipython()
    except NameError:
        raise NotImplementedError("the trace command can only be run from the Sage command-line")

    from sage.repl.preparse import preparse
    code = preparse(code)
    return pdb.run(code, ipython.user_ns)

    # this could also be useful; it drops
    # us into a debugger in an except block:
    #     import pdb; pdb.post_mortem(sys.exc_info()[2])
