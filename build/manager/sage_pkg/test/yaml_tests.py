"""
Run the PyYaml tests
"""

import sys
import os
import unittest

from develop.config import config




class YamlTestCase(unittest.TestCase):

    def setUp(self):
        self.PY_YAML = os.environ.get('PY_YAML', None)
        if self.PY_YAML is None:
            return
        PY2 = bool(sys.version_info[0] == 2)
        if PY2:
            compat = os.path.join(config.APP_DIR, 'compat', 'python-2')
            tests = os.path.join(self.PY_YAML, 'tests', 'lib')
        else:
            compat = os.path.join(config.APP_DIR, 'compat', 'python-3')
            tests = os.path.join(self.PY_YAML, 'tests', 'lib3')
        self._old_sys_path = sys.path
        sys.path = [compat, tests] + list(sys.path)
        __import__('test_all')
        self._cwd = os.getcwd()
        os.chdir(self.PY_YAML)

    def tearDown(self):
        if self.PY_YAML is None:
            return        
        sys.path = self._old_sys_path
        os.chdir(self._cwd)

    def test_yaml(self):
        if self.PY_YAML is None:
            print('skipping the PyYAML test suite (set the PY_YAML=... environment variable)')
            return
        sys.modules['test_all'].main([])



if __name__ == '__main__':
    unittest.main()
