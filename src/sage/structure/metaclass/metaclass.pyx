"""
SageMath metaclass framework

Any metaclass used in Sage is supposed to be instance (not sub-class)
of SageMetaclass. This is achieved by putting `__metaclass__ = sage_metaclass`
into the metaclass definition.
"""

try:
    import copyreg
except:
    import copy_reg as copyreg

from sage.misc.cachefunc import cached_function

class SageMetaclass(type):
    def __call__(mcls, name, bases, namespace):
        basemetacls = [type(X) for X in bases]
        basemetacls = tuple(sorted(set([mcls]+[X for X in basemetacls if
                                             not(issubclass(mcls,X))]),
                                   cmp=lambda x,y: cmp(x.__name__, y.__name__)))
        if len(basemetacls)>1:
            mcls = sage_metaclass(None, basemetacls, {})
        out = mcls.__new__(mcls,name,bases,namespace)
        if mcls.__init__ is not type.__init__:
            mcls.__init__(out, name, bases, namespace)
        return out

class DynamicSageMetaclass(SageMetaclass):
    # This is automatically created. Pickling onls works
    # under the assumption that all the methods are
    # inherited from the bases.
    def __reduce__(mcls):
        return sage_metaclass, (mcls.__name__, mcls.__bases__, {})

copyreg.pickle(DynamicSageMetaclass, DynamicSageMetaclass.__reduce__)


_cache = {}
def sage_metaclass(name, bases, namespace):
    #print "new", name, bases, namespace
    # There are two cases: If the namespace is empty,
    # then we consider a dynamically created metaclass.
    # It is uniquely determined by name and bases.
    # If it isn't empty, then we assume that the metaclass
    # is defined in a module.
    if not namespace:
        if name is None:
            name = ''.join([X.__name__.split('Metaclass')[0] for X in bases])+'Metaclass'
        try:
            return _cache[name,bases]
        except KeyError:
            pass
        mcls = DynamicSageMetaclass(name, bases, namespace)
        for X in reversed(mcls.__mro__):
            try:
                metainit = X.__dict__['__metainit__']
            except (AttributeError, KeyError):
                continue
            metainit(mcls, name, bases, namespace)
        _cache[name, bases] = mcls
        return mcls
    else:
        if not bases:
            bases = (type,)
        return SageMetaclass(name, bases, namespace)

###################################################################
#  Examples
###################################################################

#
# Unique representation

class ClasscallMetaclass:
    __metaclass__ = sage_metaclass
    def __call__(cls, *args, **opts):
        try:
            classcall = cls.__classcall__
        except AttributeError:
            return type.__call__(cls, *args, **opts)
        return classcall(cls, *args, **opts)

class MyUniqueRepresentation:
    __metaclass__ = ClasscallMetaclass
    @cached_function
    def __classcall__(cls, *args, **opts):
        out = super(cls,cls).__new__(cls, *args, **opts)
        cls.__init__(out,*args,**opts)
        out._reduction = (type(out), args, opts)
        return out
    def __hash__(self):
        return id(self)
    def __reduce__(self):
        return unreduce, self._reduction

def unreduce(cls, args, opts):
    #print cls, args, opts
    return cls(*args,**opts)

#
# Nested pickling
from sage.misc.nested_class import nested_pickle
class NestedClassMetaclass:
    __metaclass__ = sage_metaclass
    def __init__(cls, name, bases, namespace):
        nested_pickle(cls)

class MyNestedClass:
    __metaclass__ = NestedClassMetaclass
    def __init__(self, a,b):
        self.a = a
        self.b = b
    class Test(MyUniqueRepresentation):
        def bla(self):
            return 5

#
# Dynamic class
class DynamicClassMetaclass:
    __metaclass__ = sage_metaclass
    def __init__(cls, name, bases, namespace):
        if '_reduction_data' not in namespace:
            cls._reduction_data = (generic_unpickle, (cls.__module__,
                                                      cls.__name__))
    @staticmethod
    def __metainit__(mcls, name, bases, namespace):
        copyreg.pickle(mcls, DynamicClassMetaclass.__reduce__)
    def __reduce__(cls):
        red_data = getattr(cls,'_reduction_data', None)
        if red_data is None:
            clsmod = __import__(cls.__module__)
            return getattr, (clsmod, cls.__name__)
        return red_data

def generic_unpickle(mod,name):
    clsmod = __import__(mod)
    return getattr(clsmod, name)

copyreg.pickle(DynamicClassMetaclass, DynamicClassMetaclass.__reduce__)

dynamic_cache = {}
def my_dynamic_class(name, bases, extra_class=None, reduction=None):
    try:
        return dynamic_cache[name, bases,extra_class,reduction]
    except KeyError:
        pass
    if reduction is None:
        reduction = (my_dynamic_class, (name, bases, extra_class, reduction))
    if extra_class is not None:
        methods = dict(extra_class.__dict__)
        if methods.has_key("__dict__"):
            methods.__delitem__("__dict__")
        bases = extra_class.__bases__ + bases
    else:
        assert bases != ()
        extra_class = bases[0]
        methods = {}
    methods['_reduction_data'] = reduction
    methods['__module__'] = extra_class.__module__
    dynamic_cache[name, bases,extra_class,reduction]=DynamicClassMetaclass(name, bases, methods)
    return dynamic_cache[name, bases,extra_class,reduction]
