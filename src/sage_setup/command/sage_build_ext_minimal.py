from setuptools.command.build_ext import build_ext
from sage.env import THREAD_COUNT


class sage_build_ext_minimal(build_ext):
    """
    In contrast to :func:`~sage_setup.sage_build_ext.sage_build_ext`, this build extension is designed
    to be used in combination with Cython's cythonize function.
    Thus, we only take care of some options and letting Cython do the main work.
    """

    def initialize_options(self):
        build_ext.initialize_options(self)
        self.parallel = THREAD_COUNT
