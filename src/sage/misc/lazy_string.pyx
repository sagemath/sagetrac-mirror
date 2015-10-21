# -*- coding: utf-8 -*-
"""
Lazy strings

Based on speaklater: https://github.com/mitsuhiko/speaklater.

A lazy string is an object that behaves almost exactly like a string
but where the value is not computed until needed.  To define a lazy
string you specify a function that produces a string together with the
appropriate arguments for that function.  Sage uses lazy strings in
:mod:`sage.misc.misc` so that the filenames for SAGE_TMP (which
depends on the pid of the process running Sage) are not computed when
importing the Sage library.  This means that when the doctesting code
imports the Sage library and then forks, the variable SAGE_TMP depends
on the new pid rather than the old one.

EXAMPLES::

    sage: from sage.misc.lazy_string import lazy_string
    sage: L = []
    sage: s = lazy_string(lambda x: str(len(x)), L)
    sage: L.append(5)
    sage: s
    l'1'

Note that the function is recomputed each time::

    sage: L.append(6)
    sage: s
    l'2'
"""

#Copyright (c) 2009 by Armin Ronacher.
#
#Some rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are
#met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#
#    * The names of the contributors may not be used to endorse or
#      promote products derived from this software without specific
#      prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from cpython.object cimport PyObject_Call, PyObject_RichCompare

import types

def is_lazy_string(obj):
    """
    Checks if the given object is a lazy string.

    EXAMPLES::

        sage: from sage.misc.lazy_string import lazy_string, is_lazy_string
        sage: f = lambda: "laziness"
        sage: s = lazy_string(f)
        sage: is_lazy_string(s)
        True
    """
    return isinstance(obj, _LazyString)

def lazy_string(f, *args, **kwargs):
    """
    Creates a lazy string.

    INPUT:

    - ``f``, either a callable or a (format) string
    - positional arguments that are given to ``f``, either by calling or by
      applying it as a format string
    - named arguments, that are forwarded to ``f`` if it is not a string


    EXAMPLES::

        sage: from sage.misc.lazy_string import lazy_string
        sage: f = lambda x: "laziness in "+str(x)
        sage: s = lazy_string(f, ZZ); s
        l'laziness in Integer Ring'

    Here, we demonstrate that the evaluation is postponed until the value is
    needed, and that the result is not cached::

        sage: class C:
        ....:     def __repr__(self):
        ....:         print "determining string representation"
        ....:         return "a test"
        sage: c = C()
        sage: s = lazy_string("this is %s", c)
        sage: s
        determining string representation
        l'this is a test'
        sage: s == 'this is a test'
        determining string representation
        True
        sage: unicode(s)
        determining string representation
        u'this is a test'

    """
    return _LazyString(f, args, kwargs)

def _make_lazy_string(ftype, fpickle, args, kwargs):
    """
    Used for pickling.

    EXAMPLES::

        sage: from sage.misc.lazy_string import _make_lazy_string
        sage: s = _make_lazy_string(None, lambda: "laziness", (), {})
        sage: s
        l'laziness'
    """
    if ftype == 'func':
        from sage.misc.fpickle import unpickle_function
        f = unpickle_function(fpickle)
    else:
        f = fpickle
    return _LazyString(f, args, kwargs)

cdef class _LazyString(object):
    """
    Lazy class for strings created by a function call or a format string.

    INPUT:

    - ``f``, either a callable or a (format) string
    - ``args``, a tuple of arguments that are given to ``f``, either by calling
      or by applying it as a format string
    - ``kwargs``, a dictionary of optional arguments, that are forwarded to ``f``
      if it is a callable.

    .. NOTE::

        Evaluation of ``f`` is postponed until it becomes necessary, e.g., for
        comparison. The result of evaluation is not cached. The proxy
        implementation attempts to be as complete as possible, so that the
        lazy objects should mostly work as expected, for example for sorting.

        The function :func:`lazy_string` creates lazy strings in a slightly more
        convenient way, because it is then not needed to provide the arguments as
        tuple and dictionary.

    EXAMPLES::

        sage: from sage.misc.lazy_string import lazy_string, _LazyString
        sage: f = lambda x: "laziness in the " + repr(x)
        sage: s = lazy_string(f, ZZ); s
        l'laziness in the Integer Ring'
        sage: lazy_string("laziness in the %s", ZZ)
        l'laziness in the Integer Ring'

    Here, we demonstrate that the evaluation is postponed until the value is
    needed, and that the result is not cached. Also, we create a lazy string directly,
    without calling :func:`lazy_string`::

        sage: class C:
        ....:     def __repr__(self):
        ....:         print "determining string representation"
        ....:         return "a test"
        sage: c = C()
        sage: s = _LazyString("this is %s", (c,), {})
        sage: s
        determining string representation
        l'this is a test'
        sage: s == 'this is a test'
        determining string representation
        True
        sage: unicode(s)
        determining string representation
        u'this is a test'

    """

    def __init__(self, f, args, kwargs):
        """
        INPUT:

        - ``f``, either a callable or a (format) string
        - ``args``, a tuple of arguments that are given to ``f``, either by calling
          or by applying it as a format string
        - ``kwargs``, a dictionary of optional arguments, that are forwarded to ``f``
          if it is a callable.

        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda x: "laziness" + repr(x)
            sage: s = lazy_string(f, 5); s
            l'laziness5'
            sage: lazy_string("This is %s", ZZ)
            l'This is Integer Ring'
            sage: lazy_string(u"This is %s", ZZ)
            lu'This is Integer Ring'
        """
        self.func = f
        self.args = <tuple?>args
        self.kwargs = <dict?>kwargs

    cdef value(self):
        cdef f = self.func
        if isinstance(f, basestring):
            return f % self.args
        return PyObject_Call(f, self.args, self.kwargs)

    property value:
        def __get__(self):
            """
            Return the value of this lazy string, as an ordinary string.

            EXAMPLES::

                sage: from sage.misc.lazy_string import lazy_string
                sage: f = lambda: "laziness"
                sage: lazy_string(f).value
                'laziness'

            ::

                sage: from sage.misc.lazy_string import lazy_string
                sage: lazy_string("%s", "laziness").value
                'laziness'
            """
            return self.value()

    def __contains__(self, key):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: 'zi' in s
            True
            sage: 'ni' in s
            False
        """
        return key in self.value()

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: bool(lazy_string(f))
            True
            sage: f = lambda: ""
            sage: bool(lazy_string(f))
            False
        """
        return bool(self.value())

    def __dir__(self):
        """
        We assume that the underlying value provides the methods of a
        unicode string.

        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: "split" in dir(s) # indirect doctest
            True
        """
        return dir(unicode)

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: "".join(list(s)) # indirect doctest
            'laziness'
        """
        return iter(self.value())

    def __len__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: len(s)
            8
        """
        return len(self.value())

    def __str__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: str(s) # indirect doctest
            'laziness'
        """
        return str(self.value())

    def __unicode__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: unicode(s) # indirect doctest
            u'laziness'
        """
        return unicode(self.value())

    def __add__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s + " supreme"
            'laziness supreme'
        """
        if isinstance(self, _LazyString):
            return (<_LazyString>self).value() + other
        else:
            return self + (<_LazyString>other).value()

    def __mod__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laz%sss"
            sage: s = lazy_string(f)
            sage: s % "ine"
            'laziness'
            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "ine"
            sage: s = lazy_string(f)
            sage: "laz%sss" % s
            'laziness'
        """
        if isinstance(self, _LazyString):
            return (<_LazyString>self).value() % other
        else:
            return self % (<_LazyString>other).value()

    def __mul__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s * 2
            'lazinesslaziness'
            sage: 2 * s
            'lazinesslaziness'
        """
        if isinstance(self, _LazyString):
            return (<_LazyString>self).value() * other
        else:
            return self * (<_LazyString>other).value()

    def __richcmp__(self, other, int op):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s < 'laziness'
            False
            sage: s < 'azi'
            False
            sage: s < s
            False
            sage: s <= 'laziness'
            True
            sage: s <= 'azi'
            False
            sage: s <= s
            True
            sage: s == 'laziness'
            True
            sage: s == 'azi'
            False
            sage: s == s
            True
            sage: s != 'laziness'
            False
            sage: s != 'azi'
            True
            sage: s != s
            False
            sage: s > 'laziness'
            False
            sage: s > 'azi'
            True
            sage: s > s
            False
            sage: s >= 'laziness'
            True
            sage: s >= 'azi'
            True
            sage: s >= s
            True
        """
        self = (<_LazyString?>self).value()
        return PyObject_RichCompare(self, other, op)

    def __getattr__(self, name):
        """
        We pass attribute lookup through to the underlying value.

        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s.endswith('ess')
            True
            sage: s.find('i')
            3
        """
        if name == '__members__':
            return self.__dir__()
        return getattr(self.value(), name)

    def __reduce__(self):
        """
        Pickling.

        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: TestSuite(s).run() # indirect doctest
        """
        if isinstance(self.func, types.FunctionType):
            from sage.misc.fpickle import pickle_function
            f = pickle_function(self.func)
            ftype = 'func'
        else:
            f = self.func
            ftype = None
        return _make_lazy_string, (ftype, f, self.args, self.kwargs)

    def __getitem__(self, key):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s[4]
            'n'
        """
        return self.value()[key]

    def __copy__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: copy(s) is s
            True
        """
        return self

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s # indirect doctest
            l'laziness'
        """
        try:
            return 'l' + repr(self.value())
        except Exception:
            return '<%s broken>' % self.__class__.__name__

    cpdef update_lazy_string(self, args, kwds):
        """
        Change this lazy string in-place.

        INPUT:

        - ``args``, a tuple
        - ``kwds``, a dict

        .. NOTE::

            Lazy strings are not hashable, and thus an in-place change is
            allowed.

        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda op,A,B:"unsupported operand parent(s) for '%s': '%s' and '%s'"%(op,A,B)
            sage: R = GF(5)
            sage: S = GF(3)
            sage: D = lazy_string(f, '+', R, S)
            sage: D
            l"unsupported operand parent(s) for '+': 'Finite Field of size 5' and 'Finite Field of size 3'"
            sage: D.update_lazy_string(('+', S, R), {})

        Apparently, the lazy string got changed in-place::

            sage: D
            l"unsupported operand parent(s) for '+': 'Finite Field of size 3' and 'Finite Field of size 5'"

        TESTS::

            sage: D.update_lazy_string(None, None)
            Traceback (most recent call last):
            ...
            TypeError: Expected tuple, got NoneType
        """
        self.args = <tuple?>args
        self.kwargs = <dict?>kwds
