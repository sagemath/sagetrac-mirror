# -*- coding: utf-8 -*-
"""
Monkey patching 

AUTHOR:

    - Nicolas M. Thiéry
"""
class Dummy: pass
classobj = type(Dummy)

def add_entries_to_cls(cls, entries):
    for (key, value) in entries.iteritems():
        if isinstance(value, (type, classobj)):
            add_entries_to_cls(cls.__dict__[key], value.__dict__)
        else:
            setattr(cls, key, value)

class MonkeyPatchMetaclass:

    def __init__(cls, a, bases, entries):
        if bases == ():
            # Special case: initialization of the MonkeyPatch class
            assert a == "MonkeyPatch"
            return
        assert len(bases)==2
        assert bases[0] is MonkeyPatch
        add_entries_to_cls(bases[1], entries)

class MonkeyPatch:
    __metaclass__ = MonkeyPatchMetaclass

"""
    sage: load("monkey_patch.py")

    sage: class A(object):
    ...       def f(self):
    ...           return "calling A.f"
    ...       def g(self):
    ...           return "calling A.g"
    ...       class Nested:
    ...           pass

    sage: a = A()
    sage: a.f()
    'calling A.f'
    sage: a.g()
    'calling A.g'

    sage: class _(MonkeyPatch, A):
    ...       def f(self):
    ...           return "calling AMonkeyPatch.f"
    ...       class Nested:
    ...           def f(self):
    ...               return "calling AMonkeyPatch.Nested.f"
    ...

    sage: a.f()
    'calling AMonkeyPatch.f'
    sage: a.g()
    'calling A.g'
    sage: a_nested = A.Nested()
    sage: a_nested.f()
    'calling AMonkeyPatch.Nested.f'

    sage: a = A()
    sage: a.f()
    'calling AMonkeyPatch.f'
    sage: a.g()
    'calling A.g'
    sage: a_nested = A.Nested()
    sage: a_nested.f()
    'calling AMonkeyPatch.Nested.f'
"""
