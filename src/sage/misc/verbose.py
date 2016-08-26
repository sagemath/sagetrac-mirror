r"""
Verbosity System via Python's Logging Module

Howto: Using Python's Logging Module in SageMath
================================================

Import it::

    sage: import logging
    sage: logging.basicConfig()  # only needed once

Setting the level::

    sage: logging.getLogger().setLevel(logging.INFO)

Log something::

    sage: logger = logging.getLogger(__name__)
    sage: logger.info('Hello. I am talking to you.')
    INFO:__main__:Hello. I am talking to you.

Reset the level::

    sage: logging.getLogger().setLevel(logging.WARNING)

And that's all.

There is also ``logger.warning`` and ``logger.debug``, and a lot
more features, see
:python:`Logging facility for Python<library/logging.html>`.

Alternatively, this module provides
:func:`verbose`, :func:`set_verbose`, :func:`get_verbose`.


Logging Levels of SageMath and Python
=====================================

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


def verbose(msg="", t=None, level=1, caller_name=None, **kwds):
    """
    Print a message if the current verbosity is at least ``level``.

    INPUT:

    -  ``msg`` -- a string, a message to print.

    -  ``t`` -- (default: ``None``) an integer. If included, will also print
       :func:`cputime(t) <sage.misc.misc.cputime>`,
       which is the time since time `t`. Thus `t` should have
       been obtained with ``t=cputime()``.

    -  ``level`` -- (default: ``1``) an integer specifying
       the verbosity level of what we are printing.
       This level is mapped to Python's logging levels
       (see :mod:`~sage.misc.verbose` for details).

    -  ``caller_name`` -- (default: None) a string, the name
       of the calling function; in most cases Python can deduce this, so
       it need not be provided.

    - Keyword arguments can be used in the format string ``msg``.

    OUTPUT:

    The current :func:`cputime(t) <sage.misc.misc.cputime>`.
    Possibly prints a message to Python's logging module.

    EXAMPLE::

        sage: from sage.misc.verbose import verbose, set_verbose, get_verbose
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
            caller_name = '?'
        else:
            caller_name = module.__name__

    logger = logging.getLogger(caller_name)

    if t is not None:
        if msg:
            msg += ' (time = {})'.format(t)
        else:
            msg = "Finished."

    logger.log(_SageMath_to_Python_level_(level), msg, **kwds)
    return cputime()


def set_verbose(level):
    """
    Set the global SageMath verbosity level.

    INPUT:

    -  ``level`` -- an integer specifying the verbosity level.
       This level is mapped to Python's logging levels
       (see :mod:`~sage.misc.verbose` for details).

    EXAMPLES::

        sage: from sage.misc.verbose import verbose, set_verbose, get_verbose
        sage: set_verbose(2)
        sage: t = verbose("This is SageMath.", level=0)
        WARNING:?:This is SageMath.
        sage: t = verbose("This is SageMath.", level=1)
        INFO:?:This is SageMath.
        sage: t = verbose("This is SageMath.", level=2)
        DEBUG:?:This is SageMath.
        sage: t = verbose("This is Sage.", level=3)  # random
        sage: set_verbose(0)
    """
    return logging.getLogger().setLevel(_SageMath_to_Python_level_(level))


def get_verbose():
    """
    Return the global SageMath verbosity level.

    OUTPUT:

    An integer specifying the current verbosity level.
    This level is mapped to Python's logging levels
    (see :mod:`~sage.misc.verbose` for details).

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
