import os
import unittest
import importlib
import doctest
import traceback
from glob import glob

from develop.config import config




class TesterABC(object):

    def root_path(self):
        return os.path.abspath(config.HOME_DIR)

    def file_to_module(self, filename):
        """
        Convert filename to module name
        """
        cwd = self.root_path()
        fqn = os.path.abspath(filename)
        assert fqn.startswith(cwd)
        fqn = fqn[len(cwd) : -len('.py')].lstrip(os.path.sep)
        module = fqn.replace(os.path.sep, '.')
        if module.endswith('.__init__'):
            module = module[:-len('.__init__')]
        return module

    def relative_path(self, filename):
        cwd = self.root_path()
        fqn = os.path.abspath(filename)
        assert fqn.startswith(cwd)
        fqn = fqn[len(cwd):].lstrip(os.path.sep)
        return fqn

    def add_file(self, filename):
        self.status_message(self.relative_path(filename))
        #print 'Doctesting ', self.relative_path(filename)
        #doctest.testfile(self.relative_path(filename))
        module_name = self.file_to_module(filename)
        try:
            module = importlib.import_module(module_name)
        except Exception:
            traceback.print_exc()
            return
        self.add_module(module)

    def add_dir(self, *test_paths):
        for test_path in test_paths:
            for filename in glob(os.path.join(test_path, '*.py')):
                self.add_file(filename)

    def add_all(self):
        self.add_dir(self.test_path())

    def status_message(self, filename):
        pass

    def test_path(self):
        raise NotImplementedError

    def add_module(self, module):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError


class DocTester(TesterABC):
    
    def __init__(self):
        self._modules = []

    def test_path(self):
        return os.path.join(config.APP_DIR)

    def status_message(self, filename):
        print('Doc testing {}'.format(filename))

    def add_module(self, module):
        self._modules.append(module)

    def run(self):
        for m in self._modules:
            doctest.testmod(m, globs=self.global_variables())

    def global_variables(self):
        from sage_pkg.config import config
        from sage_pkg.app import Application
        from sage_pkg.tiny_git import GitRepository
        from sage_pkg.package_list import loader
        app = Application()
        return dict(
            config=config,
            app=app,
            git=GitRepository(config.path.dot_git),
            loader=loader,
        )


class UnitTester(TesterABC):

    def __init__(self):
        self.suite = unittest.TestSuite()

    def test_path(self):
        return config.TEST_DIR

    def status_message(self, filename):
        print('Unit testing {}'.format(filename))

    def add_module(self, module):
        module_suite = unittest.defaultTestLoader.loadTestsFromModule(module)
        self.suite.addTests(module_suite)

    def run(self):
        verbosity = 0
        unittest.TextTestRunner(verbosity=verbosity).run(self.suite)
