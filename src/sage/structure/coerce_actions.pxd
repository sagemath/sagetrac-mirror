
from element cimport Element, RingElement, ModuleElement
from parent cimport Parent

from sage.categories.action cimport Action
from sage.categories.map cimport Map


cdef class ModuleAction(Action):
    cdef Map connecting
    cdef extended_base

cdef class LeftModuleAction(ModuleAction):
    pass

cdef class RightModuleAction(ModuleAction):
    pass

cdef class PyScalarAction(Action):
    cdef Action _action

cdef class IntegerMulAction(Action):
    pass
