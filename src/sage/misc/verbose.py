r"""
Verbosity System and Logging in SageMath

Howto: Logging with Python's Logging Module
===========================================

Basic
^^^^^

Import it (is done automatically on start-up)::

    sage: import logging
    sage: logging.basicConfig()

Setting the level::

    sage: logging.getLogger().setLevel(logging.INFO)

Log something::

    sage: logger = logging.getLogger(__name__)
    sage: logger.info('Hello. I am talking to you.')
    INFO:__main__:Hello. I am talking to you.

If we haven't set the logging level to ``logging.INFO``, then the previous
wouldn't have been shown.
::

    sage: logger.debug('Hello. I am really talking a lot.')

The latter is not shown as the current logging level is only
``logging.INFO`` and not ``logging.DEBUG``.

Reset the level::

    sage: logging.getLogger().setLevel(logging.WARNING)

Warnings are still shown at this default level (``logging.WARNING``)::

    sage: logger.warn('Hello. I am warning you.')
    WARNING:__main__:Hello. I am warning you.

And that's all.

There are a lot more features, see
:python:`Logging facility for Python<library/logging.html>`.


More Advanced
^^^^^^^^^^^^^

Say, one wants to suppress the warnings in the following example::

    sage: R.<x,y> = CC[]
    sage: I = ideal([x^2+y^2-1,x*y-1])
    sage: I.variety()
    WARNING:sage.rings.polynomial.multi_polynomial_ideal:Warning: computations
    in the complex field are inexact;
    variety may be computed partially or incorrectly.
    WARNING:sage.rings.polynomial.multi_polynomial_ideal:Warning: falling back
    to very slow toy implementation.
    [{y: -0.86602540378443... - 0.500000000000000*I}, ...]

Simply get the corresponding logger and disable it::

    sage: logger = logging.getLogger('sage.rings.polynomial.multi_polynomial_ideal')
    sage: logger.disabled = True

Then::

    sage: R.<x,y> = CC[]
    sage: I = ideal([x^2+y^2-1,x*y-1])
    sage: I.variety()
    [{y: -0.866025403784439 - 0.500000000000000*I}, ...]

Reset the suppression::

    sage: logger.disabled = False


Various
=======

AUTHORS:

- Daniel Krenn (2016)


Functions
=========
"""
#*****************************************************************************
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import logging


def _SageMath_to_Python_level_(level):
    r"""
    Convert SageMath's verbose level to Python's logging level.

    INPUT:

    - ``level`` -- an integer.

    OUTPUT:

    An integer.

    TESTS::

        sage: from sage.misc.verbose import _SageMath_to_Python_level_
        sage: _SageMath_to_Python_level_(-3)
        50
        sage: _SageMath_to_Python_level_(-2)
        50
        sage: _SageMath_to_Python_level_(-1)
        40
        sage: _SageMath_to_Python_level_(0)
        30
        sage: _SageMath_to_Python_level_(1)
        20
        sage: _SageMath_to_Python_level_(2)
        10
        sage: _SageMath_to_Python_level_(3)
        9
    """
    from sage.misc.functional import round
    L = round(level)

    if L <= -2:
        return logging.CRITICAL
    elif L == -1:
        return logging.ERROR
    elif L == 0:
        return logging.WARNING
    elif L == 1:
        return logging.INFO
    elif L == 2:
        return logging.DEBUG
    else:
        return logging.DEBUG - 1


def _Python_to_SageMath_level_(level):
    r"""
    Convert Python's logging level to SageMath's verbose level.

    INPUT:

    - ``level`` -- an integer.

    OUTPUT:

    An integer.

    TESTS::

        sage: from sage.misc.verbose import _Python_to_SageMath_level_
        sage: import logging
        sage: _Python_to_SageMath_level_(logging.CRITICAL)
        -2
        sage: _Python_to_SageMath_level_(logging.ERROR)
        -1
        sage: _Python_to_SageMath_level_(logging.WARNING)
        0
        sage: _Python_to_SageMath_level_(logging.INFO)
        1
        sage: _Python_to_SageMath_level_(logging.DEBUG)
        2
    """
    from sage.misc.functional import round
    from sage.rings.rational_field import QQ
    L = round(QQ(level) / QQ(10)) * 10
    
    if L >= logging.CRITICAL:
        return -2
    elif L == logging.ERROR:
        return -1
    elif L == logging.WARNING:
        return 0
    elif L == logging.INFO:
        return 1
    elif L == logging.DEBUG:
        return 2
    else:
        return 3


def verbose(mesg="", t=None, level=1, caller_name=None, **kwds):
    """
    Print a message if the current verbosity is at least ``level``.

    This is part of SageMath's Verbosity System. We encourage you
    to use to Python's logging
    (:python:`Logging facility for Python<library/logging.html>`,
    see also :mod:`~sage.misc.verbose`).

    INPUT:

    -  ``mesg`` -- a string, a message to print.

    -  ``t`` -- (default: ``None``) an integer. If included, will also print
       :func:`cputime(t) <sage.misc.misc.cputime>`,
       which is the time since time `t`. Thus `t` should have
       been obtained with ``t=cputime()``.

    -  ``level`` -- (default: ``1``) an integer specifying
       the verbosity level of what we are printing.
       This level is mapped to Python's logging levels
       (see note below for details).

    -  ``caller_name`` -- (default: None) a string, the name
       of the calling function; in most cases Python can deduce this, so
       it need not be provided.

    - Keyword arguments can be used in the format string ``mesg``.

    OUTPUT:

    The current :func:`cputime(t) <sage.misc.misc.cputime>`.
    Possibly prints a message to Python's logging module.

    .. NOTE::

        Logging Levels of SageMath and Python

        .. csv-table::
            :class: contentstable
            :widths: 20, 20
            :delim: |

            SageMath | Python
            `-2` | ``logging.CRITICAL``
            `-1` | ``logging.ERROR``
            `0` | ``logging.WARNING``
            `1` | ``logging.INFO``
            `2` | ``logging.DEBUG``

    EXAMPLE::

        sage: from sage.misc.verbose import verbose, set_verbose, get_verbose

        sage: set_verbose(1)
        sage: t = verbose("This is SageMath.", level=0)
        WARNING:root:This is SageMath.
        sage: t = verbose("This is SageMath.", level=1)
        INFO:root:This is SageMath.
        sage: t = verbose("This is SageMath.", level=2)

    ::

        sage: set_verbose(1)
        sage: t = cputime()
        sage: t = verbose("This is SageMath.", t, level=1, caller_name="me")
        INFO:me:This is SageMath. (time = ...)
        sage: set_verbose(0)
    """
    import inspect
    from sage.misc.misc import cputime

    if caller_name is None:
        frame = inspect.stack()[1]
        module = inspect.getmodule(frame[0])
        if module is None:
            caller_name = None
        else:
            caller_name = module.__name__

    logger = logging.getLogger(caller_name)

    if t is not None:
        if mesg:
            mesg += ' (time = {})'.format(t)
        else:
            mesg = "Finished."

    logger.log(_SageMath_to_Python_level_(level), mesg, **kwds)
    return cputime()


def set_verbose(level):
    """
    Set the global SageMath verbosity level.

    This is part of SageMath's Verbosity System. We encourage you
    to use to Python's logging
    (:python:`Logging facility for Python<library/logging.html>`,
    see also :mod:`~sage.misc.verbose`).

    INPUT:

    -  ``level`` -- an integer specifying the verbosity level.
       This level is mapped to Python's logging levels
       (see :func:`~sage.misc.verbose.verbose` for details).

    EXAMPLES::

        sage: from sage.misc.verbose import verbose, set_verbose, get_verbose
        sage: set_verbose(2)
        sage: t = verbose("This is SageMath.", level=0)
        WARNING:root:This is SageMath.
        sage: t = verbose("This is SageMath.", level=1)
        INFO:root:This is SageMath.
        sage: t = verbose("This is SageMath.", level=2)
        DEBUG:root:This is SageMath.
        sage: t = verbose("This is Sage.", level=3)
        sage: set_verbose(0)
    """
    return logging.getLogger().setLevel(_SageMath_to_Python_level_(level))


def get_verbose():
    """
    Return the global SageMath verbosity level.

    This is part of SageMath's Verbosity System. We encourage you
    to use to Python's logging
    (:python:`Logging facility for Python<library/logging.html>`,
    see also :mod:`~sage.misc.verbose`).

    OUTPUT:

    An integer specifying the current verbosity level.
    This level is mapped to Python's logging levels
    (see :func:`~sage.misc.verbose.verbose` for details).

    EXAMPLES::

        sage: from sage.misc.verbose import verbose, set_verbose, get_verbose
        sage: get_verbose()
        0
        sage: set_verbose(2)
        sage: get_verbose()
        2
        sage: set_verbose(0)
    """
    return _Python_to_SageMath_level_(logging.getLogger().getEffectiveLevel())


# this is a version of Python's StreamHandler which prints log
# messages to the stream *currently* pointed to by sys.stderr (not the
# one when StreamHandler is set up).  This is useful in a Sage notebook, 
# where every cell has its own set of streams.

class DynamicStdErrStreamHandler(logging.StreamHandler):
    """
    A handler class which writes logging records, appropriately formatted,
    to a stream. Note that this class does not close the stream, as
    sys.stdout or sys.stderr may be used.
    """

    def __init__(self):
        import sys
        logging.StreamHandler.__init__(self, sys.stderr)
        self.parent_class = logging.StreamHandler           # save in object because name logging.StreamHandler is not available at exit

    def flush(self):
        import sys
        try:
            self.stream = sys.stderr
        except NameError:                                   # happens at exit in terminal
            pass
        self.parent_class.flush(self)

    def emit(self, record):
        import sys
        try:
            self.stream = sys.stderr
        except NameError:                                   # happens at exit in terminal
            pass
        self.parent_class.emit(self, record)


def basic_sage_logging_config(**kwargs):
    # Adapted from Python's basicConfig.
    """
    Do basic configuration for the logging system suitable for Sage.

    This function does nothing if the root logger already has handlers
    configured.

    It creates a `DynamicStdErrStreamHandler` which writes to the stream
    that is the dynamic value of `sys.stderr`, sets a formatter using the 
    BASIC_FORMAT format string, and add the handler to the root logger.
    """
    root = logging.getLogger()
    if not root.handlers:
        fs = kwargs.get("format", logging.BASIC_FORMAT)
        dfs = kwargs.get("datefmt", None)
        fmt = logging.Formatter(fs, dfs)
        hdlr = DynamicStdErrStreamHandler()
        hdlr.setFormatter(fmt)
        root = logging.getLogger()
        root.addHandler(hdlr)
        level = kwargs.get("level")
        if level is not None:
            root.setLevel(level)
