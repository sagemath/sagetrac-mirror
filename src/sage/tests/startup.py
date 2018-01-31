r"""
Ensure that certain modules are not loaded on startup.

EXAMPLES::

    sage: 'sage.libs.gap.libgap' in sys.modules
    False

Check that IPython is not imported at startup (:trac:`18726`). It is
imported by the doctest framework, so the simple test like above would
not work. Instead, we test this by starting a new Python process::

    sage: from sage.tests.cmdline import test_executable
    sage: cmd = "from sage.all import *\nprint('IPython' in sys.modules)\n"
    sage: print(test_executable(["sage", "--python"], cmd)[0])  # long time
    False

Check that numpy (:trac:`11714`) and pyparsing are not imported on startup
as they increase the startup time. Since :trac:`23696` those are imported
by the doctest framework via a matplotlib import. Again the simple test would
not work, but this time we don't have to avoid to load ipython.

    sage: from sage.tests.cmdline import test_executable
    sage: cmd = "print('numpy' in sys.modules)\n"
    sage: print(test_executable(["sage", "-c"], cmd)[0])  # long time
    False
    sage: cmd = "print('pyparsing' in sys.modules)\n"
    sage: print(test_executable(["sage", "-c"], cmd)[0])  # long time
    False
"""
