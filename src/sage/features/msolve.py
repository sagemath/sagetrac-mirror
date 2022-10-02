# -*- coding: utf-8 -*-
r"""
Feature for testing the presence of msolve

`msolve <https://msolve.lip6.fr/>`_ is a multivariate polynomial system solver.

.. SEEALSO::

    - :mod:`sage.rings.polynomial.msolve`
"""

import subprocess
from . import Executable
from . import FeatureTestResult

class msolve(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of msolve

    EXAMPLES::

        sage: from sage.features.msolve import msolve
        sage: msolve().is_present() # optional - msolve
        FeatureTestResult('msolve', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.msolve import msolve
            sage: isinstance(msolve(), msolve)
            True
        """
        Executable.__init__(self, "msolve", executable="msolve",
                            url="https://msolve.lip6.fr/")

    def is_functional(self):
        r"""
        Test if our installation of msolve is working

        EXAMPLES::

            sage: from sage.features.msolve import msolve
            sage: msolve().is_functional() # optional - msolve
            FeatureTestResult('msolve', True)
        """
        msolve_out = subprocess.run(["msolve", "-h"], capture_output=True)

        if msolve_out.returncode != 0:
            return FeatureTestResult(self, False, reason="msolve -h returned "
                                f"non-zero exit status {msolve_out.returncode}")
        elif (msolve_out.stdout[:46] !=
              b'\nmsolve library for polynomial system solving\n'):
            return FeatureTestResult(self, False,
                                     reason="output of msolve -h not recognized")
        return FeatureTestResult(self, True)

def all_features():
    return [msolve()]
