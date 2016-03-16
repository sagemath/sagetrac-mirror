"""
Detecting external softwares

This module makes up a list of external softwares that Sage interfaces. Availability
of each software is tested only when necessary.

AUTHORS:

- Kwankyu Lee (2016-03-09) -- initial version, based on the codes by Robert Bradshaw
                              and Nathann Cohen
"""

#*****************************************************************************
#       Copyright (C) 2016 KWANKYU LEE <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from multiprocessing import Array
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt

# Functions in this module whose name is of the form 'has_xxx' tests if the
# software xxx is available to Sage.
prefix = 'has_'

def has_internet():
    """
    Return ``True`` if Internet is available.

    It attempts to connect to the site "http://www.sagemath.org". Failure
    of doing this within a second is regarded as internet being not available.

    EXAMPLES::

        sage: from sage.doctest.external import has_internet
        sage: has_internet() # random
        True
    """
    from six.moves import urllib
    try:
        alarm(1)
        urllib.request.urlopen("http://www.sagemath.org")
        cancel_alarm()
        return True
    except (Exception, AlarmInterrupt):
        return False

def has_latex():
    """
    Return ``True`` if Latex is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_latex
        sage: has_latex() # random
        True
    """
    from sage.misc.latex import _run_latex_, _latex_file_
    from sage.misc.temporary_file import tmp_filename
    try:
        f = tmp_filename(ext='.tex')
        O = open(f, 'w') 
        O.write(_latex_file_('2+3'))
        O.close() 
        _run_latex_(f) 
        return True
    except Exception:
        return False

def has_magma():
    """
    Return ``True`` if Magma is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_magma
        sage: has_magma() # random
        True
    """
    from sage.interfaces.magma import magma
    try:
        magma('2+3')
        return True
    except Exception:
        return False

def has_matlab():
    """
    Return ``True`` if Matlab is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_matlab
        sage: has_matlab() # random
        True
    """
    from sage.interfaces.matlab import matlab
    try:
        matlab('2+3')
        return True
    except Exception:
        return False

def has_mathematica():
    """
    Return ``True`` if Mathematica is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_mathematica
        sage: has_mathematica() # random
        True
    """
    from sage.interfaces.mathematica import mathematica
    try:
        mathematica('2+3')
        return True
    except Exception:
        return False

def has_maple():
    """
    Return ``True`` if Maple is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_maple
        sage: has_maple() # random
        True
    """
    from sage.interfaces.maple import maple
    try:
        maple('2+3')
        return True
    except Exception:
        return False

def has_macaulay2():
    """
    Return ``True`` if Macaulay2 is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_macaulay2
        sage: has_macaulay2() # random
        True
    """
    from sage.interfaces.macaulay2 import macaulay2
    try:
        macaulay2('2+3')
        return True
    except Exception:
        return False

def has_octave():
    """
    Return ``True`` if Octave is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_octave
        sage: has_octave() # random
        True
    """
    from sage.interfaces.octave import octave
    try:
        octave('2+3')
        return True
    except Exception:
        return False

def has_scilab():
    """
    Return ``True`` if Scilab is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_scilab
        sage: has_scilab() # random
        True
    """
    from sage.interfaces.scilab import scilab
    try:
        scilab('2+3')
        return True
    except Exception:
        return False

def has_cplex():
    """
    Return ``True`` if CPLEX is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_cplex
        sage: has_cplex() # random
        True
    """
    from sage.numerical.mip import MixedIntegerLinearProgram
    try:
        MixedIntegerLinearProgram(solver='cplex')
        return True
    except Exception:
        return False

def has_gurobi():
    """
    Return ``True`` if Gurobi is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_gurobi
        sage: has_gurobi() # random
        True
    """
    from sage.numerical.mip import MixedIntegerLinearProgram
    try:
        MixedIntegerLinearProgram(solver='gurobi')
        return True
    except Exception:
        return False

def external_softwares():
    """
    Return the alphabetical list of external softwares supported by this module.
    """
    supported = list()
    for func in globals():
        if func.startswith(prefix):
            supported.append(func[len(prefix):])
    return sorted(supported)

external_softwares = external_softwares()

def _lookup(software):
    """
    Return true if the software is available on the system.
    """
    if software in external_softwares:
        func = globals().get(prefix + software)
        return func()
    else:
        return False

class AvailableSoftwares(object):
    """
    Class of the set of external softwares whose availability is detected lazily.

    EXAMPLES::

        sage: from sage.doctest.external import external_softwares,available_softwares
        sage: external_softwares
        ['cplex',
         'gurobi',
         'internet',
         'latex',
         'macaulay2',
         'magma',
         'maple',
         'mathematica',
         'matlab',
         'octave',
         'scilab']
        sage: 'internet' in available_softwares # random
        True
        sage: available_softwares.issuperset(set(['internet','latex','magma'])) # random
        True
    """
    def __init__(self):
        """
        For multiprocessing of doctests, the data ``self.seen`` should be shared among
        subprocesses, so we use :class:`Array` class from the ``multiprocessing`` module.
        """ 
        self._seen = Array('i', len(external_softwares)) # initialized to zeroes

    def __contains__(self, item):
        try:
            idx = external_softwares.index(item)
        except Exception:
            return False
        if not self._seen[idx]:
            if _lookup(item):
                self._seen[idx] = 1 # available
            else:
                self._seen[idx] = -1 # not available
        if self._seen[idx] == 1:
            return True
        elif self._seen[idx] == -1:
            return False
        else:
            raise AssertionError("Invalid value for self.seen")

    def issuperset(self, other):
        """
        Return ``True`` if ``other`` is a subset of ``self``. 
        """
        for item in other:
            if item not in self:
                return False
        return True

    def seen(self):
        """
        Return the list of detected external softwares by now.
        """
        return [external_softwares[i] for i in range(len(external_softwares)) if self._seen[i] > 0]

available_softwares = AvailableSoftwares()
