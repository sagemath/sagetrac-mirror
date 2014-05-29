"""
Combinatorial maps
"""
#*****************************************************************************
#       Copyright (C) 2011 Christian Stump <christian.stump at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.homset import Hom

_CMD_is_active = False
_constructions = []

def _unwrap(x, call=True):
    r"""
    if ``x`` is a tuple: then return _unwrap(x[0])(*x[1], **x[2])

    if ``x`` is a string: then we assume that it is a path to a function (or a
    class) and the value of the function is returned

    otherwise return ``x`` itself.

    EXAMPLES::

        sage: from sage.combinat.combinatorial_map import _unwrap
        sage: _unwrap('sage.combinat.partition.Partitions')
        Partitions
        sage: _unwrap(('sage.combinat.partition.Partitions', (5,)))
        Partitions of the integer 5
        sage: _unwrap(('sage.combinat.composition.Compositions', (19,), {'length':4}))
        Compositions of the integer 19 satisfying constraints length=4
    """
    if isinstance(x,tuple):
        f = _unwrap(x[0], call=False)
        args = ()
        if len(x) >= 2:
            args = x[1]
        kwds = {}
        if len(x) >= 3:
            kwds = x[2]

        if call:
            return f(*args, **kwds)
        else:
            return f,args,kwds

    elif isinstance(x,str):
        i = x.rfind('.')
        modname = x[:i]
        objname = x[i+1:]
        f = getattr(__import__(modname, {}, {}, [objname]), objname)
        if call:
            return f()
        else:
            return f

    return x

def build_map(f, domain, codomain, **kwds):
    domain = _unwrap(domain)
    codomain = _unwrap(codomain)

    M = Hom(domain, codomain)(f)

    # TODO: most of the time it is just impossible to rename the map
    #if 'name' in kwds:
    #    M.rename(kwds['name'])
    return M

class CombinatorialMapDecorator(object):
    r"""
    The decorator for combinatorial maps.
    """
    def __init__(self, *args, **kwds):
        self._args = args
        self._kwds = kwds

    def __call__(self, f):
        global _constructions, _CMD_is_active

        if 'domain' in self._kwds and 'codomain' in self._kwds:
            _constructions.append((build_map, (f,), self._kwds))
            if _CMD_is_active:
                M = build_map(f, **self._kwds)
                CombinatorialMapDatabase().add_map(M)
        else:
            print "Warning: %s in %s still uses the old combinatorial_map framework"%(f.__name__, f.__module__)
        return f

# Do we change the name?
combinatorial_map = CombinatorialMapDecorator

class CombinatorialMapDatabase(UniqueRepresentation):
    r"""
    The database of combinatorial maps.

    TODO::

        We need to set ``_CMD_is_active`` to ``False`` when the object is
        destroyed ! Is there a destructor in Python?

    TESTS:

    The following sanity check must pass (the only trouble right now seems to be
    the pickling)::

        sage: CMD = CombinatorialMapDatabase()
        sage: for m in CMD.map_iterator():
        ....:     TestSuite(m).run()
        ...
    """
    def __init__(self):
        self._maps_from = {}
        self._maps_to = {}

        global _CMD_is_active
        _CMD_is_active = True

        self._load_constructions()

    def maps(self):
        r"""
        Return the set of maps.
        """
        return list(self.map_iterator())

    def map_iterator(self):
        r"""
        Return an iterator over the set of maps.
        """
        for l in self._maps_from.itervalues():
            for m in l:
                yield m

    def num_maps(self):
        r"""
        Return the number of maps in the database.
        """
        from sage.rings.integer import Integer
        return Integer(sum(len(x) for x in self._maps_from.itervalues()))

    def __repr__(self):
        r"""
        String representations.
        """
        return "Combinatorial map database with %d maps"%self.num_maps()

    def _load_constructions(self):
        r"""
        Load the maps stored in the global variable _constructions
        """
        global _constructions
        for f,args,kwds in _constructions:
            try:
                M = f(*args, **kwds)
            except (ImportError,TypeError,ValueError,AttributeError),msg:
                raise ValueError("the following construction failed:\n f=%s\n args=%s\n kwds=%s\nerror message: %s"%(f,args,kwds,msg))

            self.add_map(M)

    def add_map(self, M):
        r"""
        Add the map ``M`` to the database.
        """
        if M.domain() not in self._maps_from:
            self._maps_from[M.domain()] = []
        self._maps_from[M.domain()].append(M)

        if M.codomain() not in self._maps_to:
            self._maps_to[M.codomain()] = []
        self._maps_to[M.codomain()].append(M)


