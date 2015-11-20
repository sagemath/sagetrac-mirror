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


def wrap_create_extension_list(cython_create_extension_list, library_sort_key):
    """
    Decorator for Cython's ``create_extension_list()`` to add
    Sage-specific customizations:

    - sort libraries according to ``library_sort_key``

    - add libraries (not just header files) as dependencies
    """
    def create_extension_list(*args, **kwds):
        module_list, module_metadata = cython_create_extension_list(*args, **kwds)

        for ext in module_list:
            libs = sorted(ext.libraries, key=library_sort_key)
            # Make a copy because the "depends" for different modules
            # might be identical Python objects
            depends = ext.depends[:]
            for lib in libs:
                depends += find_library_files(lib)

            ext.libraries = libs
            ext.depends = depends
            meta = module_metadata[ext.name]
            meta['libraries'] = libs
            meta['depends'] = depends

        return module_list, module_metadata
    return create_extension_list


def init_cython(library_order):
    """
    Monkeypatch ``create_extension_list``
    """
    library_sort_key = lambda x: library_order.get(x, 0)
    import Cython.Build.Dependencies as Deps
    Deps.create_extension_list = wrap_create_extension_list(
            Deps.create_extension_list, library_sort_key)
