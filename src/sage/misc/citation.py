###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from contextlib import contextmanager
import os, re

last_citations = []

@contextmanager
def citation_record(record = None):
    import cProfile, pstats
    import sage.misc.gperftools as gProfiler

    from sage.misc.citation_items.all import citation_items


    try:
        old_cpuprofile_frequency = os.environ['CPUPROFILE_FREQUENCY']
    except KeyError:
        old_cpuprofile_frequency = None
    os.environ["CPUPROFILE_FREQUENCY"] = "1000"

    cprofiler = cProfile.Profile()
    gprofiler = gProfiler.Profiler()

    cprofiler.enable()
    try:
        gprofiler.start()
    except ImportError:
        gprofiler = None

    yield

    cprofiler.disable()
    if gprofiler:
        gprofiler.stop()

    if old_cpuprofile_frequency:
        os.environ["CPUPROFILE_FREQUENCY"] = old_cpuprofile_frequency
    else:
        del os.environ["CPUPROFILE_FREQUENCY"]

    cprofiler_calls = map(cprofile_stat_to_function_string,
                          pstats.Stats(cprofiler).stats.keys())
    if gprofiler:
        gprofiler_calls = gperftools_top_to_functions(gprofiler.top(print_top=False))
    else:
        gprofiler_calls = []


    #Remove trivial functions
    bad_res = map(re.compile,
                  [r'is_.*Element'])
    calls = [c for c in cprofiler_calls + gprofiler_calls
             if all(r.match(c) is None for r in bad_res)]


    #Check to see which citations appear in the profiled run
    called_items = [item for item in citation_items
                    if any(r.match(c) is not None for r in item.re() for c in calls)]

    if record is None:
        import sage
        sage.misc.citation.last_citations = called_items

        print "The computation used the following components."
        print "They are stored as a list in sage.misc.citation.last_citations."
        print "  ", ", ".join(map(repr, called_items))
    else:
        record.extend(called_items)


def cprofile_stat_to_function_string(stat_key):
    if stat_key[0] == '~':
        if stat_key[2].startswith("<method "):
            object_start = stat_key[2].find("of '") + len("of '")
            object_end = stat_key[2].find("'", object_start)
            module_part = stat_key[2][object_start:object_end]

            function_start = stat_key[2].find("method '") + len("method '")
            function_end = stat_key[2].find("'", function_start)
            function_part = stat_key[2][function_start:function_end]
        else:
            module_part = None
            function_part = stat_key[2][1:-1]
    else:
        module_part = (stat_key[0].split("site-packages/")[-1]
                       .split('.')[0]
                       .replace("/", "."))
        function_part = stat_key[2]

    if module_part:
        return module_part + "." + function_part
    else:
        return function_part
    
def gperftools_top_to_functions(top):
    split = re.compile("_[0123456789]+").split
    lines = [l.rstrip().split()[-1] for l in top.splitlines()]

    return [".".join(split(l)[1:])
            for l in lines
            if l.startswith("__pyx")]
            
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
    deprecation(1, 'get_sytems is replaced by citation_record')


    import cProfile, inspect, pstats, re
    from sage.misc.all import preparse, tmp_filename
    from sage.env import SAGE_ROOT
    from sage.misc.citation_items.all import citation_items as systems

    if not isinstance(cmd, basestring):
        raise TypeError("command must be a string")

    cmd = preparse(cmd)

    #Run the command and get the stats
    filename = tmp_filename()
    cProfile.runctx(cmd, inspect.stack()[1][0].f_globals, {}, filename)
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
        res = [r.replace("^", "").replace(".*", "") for r in system._re]
        if any([(r in s) or (r.replace('.','/') in s) for r in res for s in strings]):
            systems_used.append(repr(system))
    return systems_used
