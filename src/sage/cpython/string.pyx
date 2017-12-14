from __future__ import absolute_import

from six import *



cpdef string_to_bytes(x):
    r"""
    Convert ``x`` to bytes using the utf-8 encoding.

    INPUT:

    - ``x`` -- a string (bytes or unicode)

    OUTPUT:

    an object of type ``bytes``

    Python2 behaviour:

    If input is str, returns the input.

    If input is unicode, convert to bytes using utf8-encoding.

    Python3 behaviour:

    If input is str, convert to bytes using utf8-encoding.

    If input is bytes, returns the input.

    EXAMPLES::

        sage: from sage.cpython.string import string_to_bytes
        sage: string_to_bytes("500 €")
        '500 \xe2\x82\xac'
        sage: string_to_bytes(u"500 \u20ac")
        '500 \xe2\x82\xac'
    """
    if isinstance(x, text_type):  # py2 unicode and py3 str
        return x.encode("utf-8")
    if isinstance(x, bytes):
        return x
    raise TypeError('input has no conversion to unicode')
