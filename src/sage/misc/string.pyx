# -*- encoding: utf-8 -*-
"""
String conversion and encoding/decoding utilities, in particular for Python 2/3
compatibility.
"""

#*****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.bray@lri.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import

import locale
import six
import sys

from libc.string cimport strlen

from cpython.bytes cimport PyBytes_AsString as PyString_AsString
if six.PY2:
    from cpython.string cimport PyString_FromString
else:
    from cpython.bytes cimport PyBytes_FromString as PyString_FromString
from cpython.unicode cimport (PyUnicode_Decode, PyUnicode_AsEncodedString)
from cpython.version cimport PY_MAJOR_VERSION

cdef extern from "Python.h":
    # Missing from cpython.unicode
    char* PyUnicode_AsUTF8(object unicode)


DEFAULT_ENCODING = locale.getpreferredencoding()
FS_ENCODING = sys.getfilesystemencoding()


cdef inline str char_to_str(char* c, encoding=DEFAULT_ENCODING):
    """
    Converts a C ``char`` array to a Python `str` object.

    On Python 3 this requires an encoding to be specified with which to
    decode the bytes in the ``char`` array.
    """

    if PY_MAJOR_VERSION <= 2:
        # <str> is needed here to "trick" Cython into thinking we expect
        # this to return a str (on Python 2 it has no problem with this,
        # but on Python 3 it balks)
        # It doesn't matter that this cast doesn't make sense since it
        # will never happen on Python 3
        return <str>PyString_FromString(c)
    else:
        return PyUnicode_Decode(c, strlen(c), PyUnicode_AsUTF8(encoding),
                                "surrogateescape")


cpdef inline str bytes_to_str(bytes b, encoding=DEFAULT_ENCODING):
    """
    Convertes `bytes` to `str`.

    On Python 2 this is a no-op since ``bytes is str``.  On Python 3
    this decodes the given `bytes` to a Python 3 unicode `str` using the
    specified encoding.

    EXAMPLES::

        sage: from six import PY2; from sage.misc import six
        sage: six.DEFAULT_ENCODING = 'utf-8'
        sage: s = six.bytes_to_str(b'\xe2\x98\x83')
        sage: if PY2:
        ....:     s == b'\xe2\x98\x83'
        ....: else:
        ....:     s == u'☃'
        True
    """

    return char_to_str(PyString_AsString(b), encoding=encoding)


cpdef inline bytes str_to_bytes(str s, encoding=DEFAULT_ENCODING):
    """
    Convertes `str` to `bytes`.

    On Python 2 this is a no-op since ``str is bytes``.  On Python 3
    this encodes the given `str` to a Python 3 `bytes` using the
    specified encoding.

    EXAMPLES::

        sage: from six import PY2; from sage.misc import six
        sage: six.DEFAULT_ENCODING = 'utf-8'
        sage: if PY2:
        ....:     b = six.str_to_bytes('\xe2\x98\x83')
        ....: else:
        ....:     b = six.str_to_bytes(u'☃')
        sage: b == b'\xe2\x98\x83'
        True
    """
    if PY_MAJOR_VERSION <= 2:
        return <bytes>s
    else:
        return PyUnicode_AsEncodedString(s, PyUnicode_AsUTF8(encoding),
                                         "surrogateescape")
