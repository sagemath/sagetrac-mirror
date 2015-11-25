"""
Optional extensions

An "optional extension" is a Cython extension which is always
cythonized (i.e. converted to a .c or .cpp file), but which is only
compiled depending on some condition. Typically, this condition is a
package which must be installed.
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


import os
from distutils.extension import Extension
from sage.misc.package import is_package_installed


class CythonizeExtension(Extension):
    """
    A class for extensions which are only cythonized, but not built.

    The file ``src/setup.py`` contains some logic to check the
    ``skip_build`` attribute of extensions.

    EXAMPLES::

        sage: from sage_setup.optional_extension import CythonizeExtension
        sage: ext = CythonizeExtension("foo", ["foo.c"])
        sage: ext.skip_build
        True
    """
    skip_build = True


c_code_if_package_not_installed = """
#include <Python.h>

PyMODINIT_FUNC init{name}(void)
<%
    PyObject* modpackage = PyImport_ImportModule("sage.misc.package");
    if (!modpackage) return;

    PyObject* exc = PyObject_GetAttrString(modpackage, "PackageNotFoundError");
    if (!exc) return;

    PyErr_SetString(exc, "{package}");
%>
"""

class PackageNotInstalledExtension(Extension):
    """
    Abstract base class for extensions which depend on a package.
    Due to the way that Cython works, this must be subclassed, where
    the subclass must add a ``package`` attribute.

    EXAMPLES::

        sage: from sage_setup.optional_extension import CythonizeExtension
        sage: ext = CythonizeExtension("foo", ["foo.c"])
        sage: ext.skip_build
        True
    """
    @classmethod
    def package_class(cls, package):
        """
        Return a class for extensions which depend on package
        ``package``.

        EXAMPLES::

            sage: from sage_setup.optional_extension import PackageNotInstalledExtension
            sage: PackageNotInstalledExtension.package_class("mypkg")
            <class 'sage_setup.optional_extension.mypkgNotInstalledExtension'>
        """
        name = package + "NotInstalledExtension"
        return type(cls)(name, (cls, object), dict(package=package))

    def prepare(self, distutils_build_ext):
        """
        Convert the extension to an extension to compile a C file
        raising ``PackageNotFoundError``.

        EXAMPLES::

            sage: from sage_setup.optional_extension import PackageNotInstalledExtension
        """
        name = self.name.split('.')[-1]
        build_py = distutils_build_ext.get_finalized_command('build')
        temp_dir = build_py.build_temp

        csource = os.path.join(temp_dir, name + ".c")
        with open(csource, "w") as f:
            f.write(c_code_if_package_not_installed.format(
                    name=name, package=self.package))

        self.sources = [csource]
        self.depends = [__file__]
        self.libraries = []

    def __reduce__(self):
        """
        Pickle ``self``.

        We need a custom ``__reduce__`` method because the class is
        created at run-time.

        See :func:`make_PackageNotInstalledExtension`

        EXAMPLES::

            sage: from sage_setup.optional_extension import PackageNotInstalledExtension
            sage: cls = PackageNotInstalledExtension.package_class("mypkg")
            sage: ext = cls("foo", ["foo.c"])
            sage: ext.__reduce__()
            (<function make_PackageNotInstalledExtension at 0x...>,
             ('mypkg',),
             {'define_macros': [],
              'depends': [],
              'export_symbols': [],
              'extra_compile_args': [],
              'extra_link_args': [],
              'extra_objects': [],
              'include_dirs': [],
              'language': None,
              'libraries': [],
              'library_dirs': [],
              'name': 'foo',
              'runtime_library_dirs': [],
              'sources': ['foo.c'],
              'swig_opts': [],
              'undef_macros': []})
        """
        return make_PackageNotInstalledExtension, (self.package,), self.__dict__


def make_PackageNotInstalledExtension(package):
    """
    Helper function for unpickling.

    EXAMPLES::

        sage: from sage_setup.optional_extension import PackageNotInstalledExtension
        sage: cls = PackageNotInstalledExtension.package_class("mypkg")
        sage: ext = cls("foo", ["foo.c"])
        sage: ext.name, ext.package
        ('foo', 'mypkg')
        sage: ext = loads(dumps(ext))
        sage: ext.name, ext.package
        ('foo', 'mypkg')
    """
    t = PackageNotInstalledExtension.package_class(package)
    return t.__new__(t)


def OptionalExtension(*args, **kwds):
    """
    If some condition (see INPUT) is satisfied, return an ``Extension``.
    Otherwise, return a ``PackageNotInstalledExtension`` or a
    ``CythonizeExtension``.

    Typically, the condition is some optional package or something
    depending on the operating system.

    INPUT:

    - ``condition`` -- (boolean) the actual condition

    - ``package`` -- (string) the condition is that this package is
      installed (only used if ``condition`` is not given)

    EXAMPLES::

        sage: from sage_setup.optional_extension import OptionalExtension
        sage: ext = OptionalExtension("foo", ["foo.c"], condition=False)
        sage: print ext.__class__
        sage_setup.optional_extension.CythonizeExtension
        sage: ext = OptionalExtension("foo", ["foo.c"], condition=True)
        sage: print ext.__class__
        distutils.extension.Extension
        sage: ext = OptionalExtension("foo", ["foo.c"], package="no_such_package")
        sage: print ext.__class__
        <class 'sage_setup.optional_extension.no_such_packageNotInstalledExtension'>
        sage: ext = OptionalExtension("foo", ["foo.c"], package="pari")
        sage: print ext.__class__
        distutils.extension.Extension
    """
    cls = Extension
    try:
        if not kwds.pop("condition"):
            cls = CythonizeExtension
    except KeyError:
        pkg = kwds.pop("package")
        if not is_package_installed(pkg):
            cls = PackageNotInstalledExtension.package_class(pkg)

    return cls(*args, **kwds)
