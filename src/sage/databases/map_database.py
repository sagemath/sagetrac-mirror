"""
Map database
"""
#*****************************************************************************
#       Copyright (C) 2014 Vincent Delecroix <20100.delecroix at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.sage_object import SageObject
from sage.categories.morphism import Morphism

#TODO: to be moved in sage/categories/morphism.pyx
class MorphismFromElementMethod(Morphism):
    r"""
    A morphism defined from a method to be applied to an element.
    """
    def __init__(self, parent, method_name):
        Morphism.__init__(self, parent)
        self._method_name = method_name

    def _call_(self, elt):
        return getattr(elt, self._method_name)()

# the following list is built from the decorator. It will contain a list of
# tuples of the form (domain, codomain, category, method, dict_of_extra_args)
_methods_as_maps = []

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

class RegisterMethodAsMapDecorator(object):
    r"""
    A decorator for methods that actually correspond to a map.

    Each method tagged with this decorator fills a list that can be subsequently
    used by :class:`MapDatabase`.

    EXAMPLES::

        sage: pass
    """
    def __init__(self, *args, **kwds):
        self._args = args
        self._kwds = kwds

    def __call__(self, f):
        global _methods_as_maps

        kwds = self._kwds.copy()

        category = None
        if 'category' in kwds:
            category = kwds['category']
            del kwds['category']

        domain = None
        if 'domain' in kwds:
            domain = kwds['domain']
            del kwds['domain']

        codomain = None
        if 'codomain' in kwds:
            codomain = kwds['codomain']
            del kwds['codomain']

        _methods_as_maps.append((domain, codomain, category, f, kwds))

        return f

register_method_as_map = RegisterMethodAsMapDecorator

# Do we get rid of combinatorial_map... not yet
combinatorial_map = register_method_as_map

class MapDatabase(SageObject):
    r"""
    Database of maps.

    TESTS:

    The following sanity check must pass (the only trouble right now seems to be
    the pickling)::

        sage: CMD = MapDatabase()
        sage: for m in CMD.map_iterator():
        ....:     TestSuite(m).run()
    """
    def __init__(self, load_constructions=True):
        self._maps_from = {}
        self._maps_to = {}

        if load_constructions:
            self._load_constructions()

    def __contains__(self, elt):
        r"""
        Check wether the database contains the morphism ``elt``.

        EXAMPLES::

            sage: m = QQ.convert_map_from(ZZ); m
            sage: m in MapDatabase()
            False
        """
        if not isinstance(elt, Morphism):
            return False

        return elt in self._maps_from[elt.domain()]

    def maps(self):
        r"""
        Return the list of maps of this database.
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

    def _load_constructions(self, modules=None):
        r"""
        Load the maps stored in the global variable ``_methods_as_maps`` that
        were registered from the decorator ``register_method_as_map``.

        If ``modules`` is set to a list of modules (i.e. string of the form
        ['sage.combinat.partition', 'sage.combinat.permutation', ...]) then only
        the maps that belong to these modules are loaded.
        """
        global _methods_as_maps
        for domain,codomain,category,method,extra_kwds in _methods_as_maps:
            if modules and method.__module__ not in modules:
                continue

            if domain is None or codomain is None:
                print "Warning: %s in %s do not specify domain or codomain and will not be present in the database"%(
                        method.__name__, method.__module__)
                continue

            from sage.categories.homset import Hom
            try:
                domain = _unwrap(domain)
            except (ImportError,TypeError,ValueError,AttributeError),msg:
                raise ValueError("Error while unwrapping the domain %s for the method %s in %s... original error: %s"%(
                             domain, method.__name__, method.__module__, msg))

            try:
                codomain = _unwrap(codomain)
            except (ImportError,TypeError,ValueError,AttributeError),msg:
                raise ValueError("Error while unwrapping the codomain %s for the method %s in %s... original error: %s"%(
                             codomain, method.__name__, method.__module__, msg))

            try:
                H = Hom(domain,codomain,category)
                M = MorphismFromElementMethod(H, method.__name__)
            except (ImportError,TypeError,ValueError,AttributeError),msg:
                raise ValueError("the following construction failed:\n domain=%s\n codomain=%s\n method_name=%s\nerror message: %s"%(domain,codomain,method.__name__))

            if 'name' in extra_kwds:
                if domain == codomain:
                    map_name = 'Endomorphism %s of %s'%(extra_kwds['name'], domain)
                else:
                    map_name = 'Morphism %s\n  From: %s\n  To:   %s'%(extra_kwds['name'],
                            domain, codomain)
                M.rename(map_name)

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


