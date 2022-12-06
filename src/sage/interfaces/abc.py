r"""
Abstract base classes for interface elements
"""
from functools import update_wrapper
def _add_doc(c):
        c.__doc__ = "\n Abstract base class for :class:`~sage.interfaces.gap." + c.__name__
        c.__doc__ += """
    
  This class is defined for the purpose of ``isinstance`` tests.  It should not be
  instantiated.

  EXAMPLES:

  By design, there is a unique direct subclass::

  """
        c.__doc__ += "    sage: len(sage.interfaces.abc."+c.__name__+".__subclasses__()) <= 1\n      True"
        return c


@_add_doc
class GapElement: pass


@_add_doc
class GpElement: pass


class Macaulay2Element:
    r"""
    Abstract base class for :class:`~sage.interfaces.macaulay2.Macaulay2Element`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES:

    By design, there is a unique direct subclass::

        sage: len(sage.interfaces.abc.Macaulay2Element.__subclasses__()) <= 1
        True
    """
    pass


class SingularElement:
    r"""
    Abstract base class for :class:`~sage.interfaces.singular.SingularElement`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES:

    By design, there is a unique direct subclass::

        sage: len(sage.interfaces.abc.SingularElement.__subclasses__()) <= 1
        True
    """
    pass
