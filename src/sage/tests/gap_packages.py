"""
Test the optional GAP packages

TESTS::

    sage: from sage.tests.gap_packages import all_installed_packages
    sage: for name in all_installed_packages():       # optional: gap_packages
    ....:     stat = libgap.eval('LoadPackage("{0}")'.format(name))
    ....:     prefix = ' ' if str(stat) == 'true' else 'x'
    ....:     print("{2} : {0: <10} got '{1}' when loading.".format(name,stat,prefix))
      : Alnuth     got 'true' when loading.
      : GAPDoc     got 'true' when loading.
    x : HAPcryst   got 'fail' when loading.
      : Hap        got 'true' when loading.
      : autpgrp    got 'true' when loading.
      : braid      got 'true' when loading.
      : crime      got 'true' when loading.
      : ctbllib    got 'true' when loading.
      : design     got 'true' when loading.
      : factint    got 'true' when loading.
      : grape      got 'true' when loading.
      : guava      got 'true' when loading.
      : laguna     got 'true' when loading.
      : polycyclic got 'true' when loading.
      : polymaking got 'true' when loading.
      : sonata     got 'true' when loading.
      : toric      got 'true' when loading.
"""

import os
import os.path

from sage.libs.gap.libgap import libgap



def all_installed_packages():
    """
    Return list of all installed packages.

    OUTPUT:
    
    Tuple of strings in alphabetic order.
    
    EXAMPLES::

        sage: from sage.tests.gap_packages import all_installed_packages
        sage: all_installed_packages()
        (...'GAPDoc',...)
    """
    packages = []
    for path in libgap.eval('GAP_ROOT_PATHS').sage():
        if path.endswith('/.gap/'):
            continue
        pkg_dir = os.path.join(path, 'pkg')
        if not os.path.exists(pkg_dir):
            continue
        for subdir in os.listdir(pkg_dir):
            if not os.path.isdir(os.path.join(pkg_dir, subdir)):
                continue
            packages.append(subdir.rstrip('-.0123456789'))
    packages.sort()
    return tuple(packages)
