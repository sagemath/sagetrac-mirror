# Copyright (C) 2014 Martin Raum

# Author: Martin Raum (martin@raum-brothers.eu)

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from contextlib import contextmanager
import re
from sage.env import SAGE_ROOT

@contextmanager
def citation_record(record):
    import cProfile, pstats
    from sage.misc.citation_items import systems

    profile = cProfile.Profile()

    profile.enable()
    yield
    profile.disable()

    stats = pstats.Stats(profile)
    strings = [a[0].replace(SAGE_ROOT, "") + " " + a[2] for a in stats.stats.keys()]

    #Remove trivial functions
    bad_res = [re.compile(r'is_.*Element')]
    for bad_re in bad_res:
        i = 0
        while i < len(strings):
            if bad_re.findall(strings[i]):
                strings.pop(i)
            else:
                i += 1

    #Check to see which systems appear in the profiled run
    for system in systems:
        if any([(r in s) or (r.replace('.','/') in s) for r in systems[system] for s in strings]):
            record.append(system)

def get_systems(cmd):
    """
    Returns a list of the systems used in running the command
    cmd.  Note that the results can sometimes include systems
    that did not actually contribute to the computation. Due
    to caching and the inability to follow all C calls, it
    could miss some dependencies as well.

    INPUT:

    - ``cmd`` - a string to run

    EXAMPLES::

        sage: from sage.misc.citation import get_systems
        sage: s = get_systems('integrate(x^2, x)'); #priming coercion model
        sage: get_systems('integrate(x^2, x)')
        ['ginac', 'Maxima']
        sage: R.<x,y,z> = QQ[]
        sage: I = R.ideal(x^2+y^2, z^2+y)
        sage: get_systems('I.primary_decomposition()')
        ['Singular']

        sage: a = var('a')
        sage: get_systems('((a+1)^2).expand()')
        ['ginac', 'GMP']
    """
    from sage.misc.superseded import deprecation
    ## TODO: insert ticket number
    deprecation(0, 'get_sytems is replaced by citation_record')


    import cProfile, pstats, re

    if not isinstance(cmd, basestring):
        raise TypeError("command must be a string")

    cmd = preparse(cmd)

    #Run the command and get the stats
    filename = tmp_filename()
    cProfile.runctx(cmd, globals(), {}, filename)
    stats = pstats.Stats(filename)

    #Strings is a list of method names and modules which get run
    strings = [a[0].replace(SAGE_ROOT, "") + " " + a[2] for a in stats.stats.keys()]

    #Remove trivial functions
    bad_res = [re.compile(r'is_.*Element')]
    for bad_re in bad_res:
        i = 0
        while i < len(strings):
            if bad_re.findall(strings[i]):
                strings.pop(i)
            else:
                i += 1

    #Check to see which systems appear in the profiled run
    systems_used = []
    for system in systems:
        if any([(r in s) or (r.replace('.','/') in s) for r in systems[system] for s in strings]):
            systems_used.append(system)
    return systems_used
