# -*- coding: utf-8 -*-
r"""
Check various graph generator programs
"""

import os
import subprocess

from . import Executable, FeatureTestResult


class Plantri(Executable):
    r"""
    A class:`sage.features.graph_generators.Feature` which checks for the ``plantri``
    binary.

    EXAMPLES::

        sage: from sage.features.graph_generators import Plantri
        sage: Plantri().is_present()  # optional: plantri
        FeatureTestResult('plantri', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graph_generators import Plantri
            sage: isinstance(Plantri(), Plantri)
            True
        """
        Executable.__init__(self, name="plantri", spkg="plantri", executable="plantri", url="http://users.cecs.anu.edu.au/~bdm/plantri/")

    def is_functional(self):
        r"""
        Check whether ``plantri`` works on trivial input.

        EXAMPLES::

            sage: from sage.features.graph_generators import Plantri
            sage: Plantri().is_functional()  # optional: plantri
            FeatureTestResult('plantri', True)
        """
        command = ["plantri", "4"]
        try:
            lines = subprocess.check_output(command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                    reason="Call `{command}` failed with exit code {e.returncode}".format(command=" ".join(command), e=e))

        expected = "1 triangulations written"
        if lines.find(expected) == -1:
            return FeatureTestResult(self, False,
                    reason = "Call `{command}` did not produce output which contains `{expected}`".format(command=" ".join(command), expected=expected))

        return FeatureTestResult(self, True)


class Buckygen(Executable):
    r"""
    A class:`sage.features.graph_generators.Feature` which checks for the ``buckygen``
    binary.

    EXAMPLES::

        sage: from sage.features.graph_generators import Buckygen
        sage: Buckygen().is_present()  # optional: buckygen
        FeatureTestResult('Buckygen', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graph_generators import Buckygen
            sage: isinstance(Buckygen(), Buckygen)
            True
        """
        Executable.__init__(self, name="Buckygen", spkg="buckygen", executable="buckygen", url="http://caagt.ugent.be/buckygen/")

    def is_functional(self):
        r"""
        Check whether ``buckygen`` works on trivial input.

        EXAMPLES::

            sage: from sage.features.graph_generators import Buckygen
            sage: Buckygen().is_functional()  # optional: buckygen
            FeatureTestResult('Buckygen', True)
        """
        command = ["buckygen", "-d", "22d"]
        try:
            lines = subprocess.check_output(command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                    reason = "Call `{command}` failed with exit code {e.returncode}".format(command=" ".join(command), e=e))

        expected = "Number of fullerenes generated with 13 vertices: 0"
        if lines.find(expected) == -1:
            return FeatureTestResult(self, False,
                    reason = "Call `{command}` did not produce output which contains `{expected}`".format(command=" ".join(command), expected=expected))

        return FeatureTestResult(self, True)


class Benzene(Executable):
    r"""
    A class:`sage.features.graph_generators.Feature` which checks for the ``benzene``
    binary.

    EXAMPLES::

        sage: from sage.features.graph_generators import Benzene
        sage: Benzene().is_present()  # optional: benzene
        FeatureTestResult('Benzene', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graph_generators import Benzene
            sage: isinstance(Benzene(), Benzene)
            True
        """
        Executable.__init__(self, name="Benzene", spkg="benzene", executable="benzene", url="http://www.grinvin.org/")

    def is_functional(self):
        r"""
        Check whether ``benzene`` works on trivial input.

        EXAMPLES::

            sage: from sage.features.graph_generators import Benzene
            sage: Benzene().is_functional()  # optional: benzene
            FeatureTestResult('Benzene', True)
        """
        devnull = open(os.devnull, 'wb')
        command = ["benzene", "2", "p"]
        try:
            lines = subprocess.check_output(command, stderr=devnull)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                    reason="Call `{command}` failed with exit code {e.returncode}".format(command=" ".join(command), e=e))

        expected = ">>planar_code<<"
        if not lines.startswith(expected):
            return FeatureTestResult(self, False,
                    reason="Call `{command}` did not produce output that started with `{expected}`.".format(command=" ".join(command), expected=expected))

        return FeatureTestResult(self, True)
