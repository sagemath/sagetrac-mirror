"""
Sage customizations to Cython
"""
#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage_setup.find import find_library_files


class distutils_kwds_hook(object):
    """
    A hook function to process extension module metadata with
    Sage-specific customizations:

    - sort libraries according to ``library_order``

    - add libraries (not just header files) as dependencies

    INPUT:

    - ``library_order`` -- A :class:`dict` with ``library:order``
      pairs, where the libraries with the smallest order (which can
      be negative) are put first. Libraries which do not appear in
      the dict are assumed to have order 0.
    """
    def __init__(self, library_order):
        """
        TESTS::

            sage: from sage_setup.cython import distutils_kwds_hook
            sage: hook = distutils_kwds_hook({})
            sage: hook.library_order
            {}
        """
        self.library_order = library_order

    def sort_key(self, lib):
        """
        A key to use for sorting libraries.

        EXAMPLES::

            sage: from sage_setup.cython import distutils_kwds_hook
            sage: order = {'first': -1, 'last': 1}
            sage: hook = distutils_kwds_hook(order)
            sage: hook.sort_key("first")
            -1
            sage: hook.sort_key("unknown")
            0
            sage: sorted(["last", "first", "other"], key=hook.sort_key)
            ['first', 'other', 'last']
        """
        return self.library_order.get(lib, 0)

    def __call__(self, kwds):
        """
        Hook for Cython's ``create_extension_modules()``.

        EXAMPLES::

            sage: from sage_setup.cython import distutils_kwds_hook
            sage: order = {'first': -1, 'last': 1}
            sage: hook = distutils_kwds_hook(order)
            sage: hook({})
            {}
            sage: libs = ["last", "first", "other"]
            sage: hook(dict(libraries=libs, sources=["ext.c"]))
            {'depends': [], 'libraries': ['first', 'other', 'last'], 'sources': ['ext.c']}
        """
        try:
            libs = kwds['libraries']
        except KeyError:
            pass
        else:
            libs = sorted(libs, key=self.sort_key)

            try:
                # Make a copy because the "depends" for different
                # modules might be identical Python objects
                depends = kwds['depends'][:]
            except KeyError:
                depends = []
            for lib in libs:
                depends += find_library_files(lib)

            kwds['libraries'] = libs
            kwds['depends'] = depends
        return kwds
