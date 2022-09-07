# -*- coding: utf-8 -*-
r"""
Feature for testing the presence of ``open``, ``xdg-open``, ``cygstart``

These are programs which are used to open image files according to the
system preferences.
"""
from . import Executable, FeatureTestResult

class Open(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``open``

    This executable is available on OS X.

    EXAMPLES::

        sage: from sage.features.xdgopen import Open
        sage: Open().is_present()  # optional - open
        FeatureTestResult('open', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.xdgopen import Open
            sage: isinstance(Open(), Open)
            True
        """
        Executable.__init__(self, "open", executable="open",)

class XdgOpen(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``xdg-open``

    This executable is available on typical linux systems.

    EXAMPLES::

        sage: from sage.features.xdgopen import XdgOpen
        sage: XdgOpen().is_present()  # optional - xdg-open
        FeatureTestResult('xdg-open', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.xdgopen import Open
            sage: isinstance(Open(), Open)
            True
        """
        Executable.__init__(self, "xdg-open", executable="xdg-open",
                            url="https://freedesktop.org/wiki/Software/xdg-utils/")

class Cygstart(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``cygstart``

    This executable is available on Cygwin.

    EXAMPLES::

        sage: from sage.features.xdgopen import Cygstart
        sage: Cygstart().is_present()  # optional - cygstart
        FeatureTestResult('cygstart', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.xdgopen import Cygstart
            sage: isinstance(Cygstart(), Cygstart)
            True
        """
        Executable.__init__(self, "cygstart", executable="cygstart")

def all_features():
    return [Open(), XdgOpen(), Cygstart()]
