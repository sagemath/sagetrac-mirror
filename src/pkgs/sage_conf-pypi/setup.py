import os
import sys
import shutil
import sysconfig

from setuptools import setup
from distutils.command.build_scripts import build_scripts as distutils_build_scripts
from distutils.command.build_py import build_py as distutils_build_py
from distutils.errors import (DistutilsSetupError, DistutilsModuleError,
                              DistutilsOptionError)

class build_py(distutils_build_py):

    def run(self):
        DOT_SAGE = os.environ.get('DOT_SAGE', os.path.join(os.environ.get('HOME'), '.sage'))
        # config.status and other configure output has to be writable.
        SAGE_ROOT = os.path.join(DOT_SAGE, 'sage-{}-{}'.format(self.distribution.version,
                                                               sysconfig.get_config_var('SOABI')))
        SAGE_LOCAL = os.path.join(SAGE_ROOT, 'local')
        if os.path.exists(os.path.join(SAGE_ROOT, 'config.status')):
            print('Reusing {}'.format(SAGE_ROOT))
        else:
            try:
                shutil.copytree('sage_root', SAGE_ROOT)  # will fail if already exists
            except Exception:
                raise DistutilsSetupError("the directory SAGE_ROOT={} already exists but it is not configured.  Please remove it and try again.".format(SAGE_ROOT))
            cmd = "cd {} && ./configure --prefix={} --with-python={} --with-system-python3=force".format(SAGE_ROOT, SAGE_LOCAL, sys.executable)
            print("Running {}".format(cmd))
            if os.system(cmd) != 0:
                raise DistutilsSetupError("configure failed")

        cmd = "cd {} && make V=0 build-local".format(SAGE_ROOT)
        if os.system(cmd) != 0:
            raise DistutilsSetupError("make build-local failed")

        # Install configuration
        shutil.copyfile(os.path.join(SAGE_ROOT, 'build', 'pkgs', 'sage_conf', 'src', 'sage_conf.py'),
                        'sage_conf.py')
        if not self.distribution.py_modules:
            self.py_modules = self.distribution.py_modules = []
        self.distribution.py_modules.append('sage_conf')
        shutil.copyfile(os.path.join(SAGE_ROOT, 'src', 'bin', 'sage-env-config'),
                        os.path.join('bin', 'sage-env-config'))
        distutils_build_py.run(self)

class build_scripts(distutils_build_scripts):

    def run(self):
        self.distribution.scripts.append(os.path.join('bin', 'sage-env-config'))
        if not self.distribution.entry_points:
            self.entry_points = self.distribution.entry_points = dict()
        # if 'console_scripts' not in self.distribution.entry_points:
        #     self.distribution.entry_points['console_scripts'] = []
        # self.distribution.entry_points['console_scripts'].append('sage-config=sage_conf:_main')
        distutils_build_scripts.run(self)

setup(
    cmdclass=dict(build_py=build_py, build_scripts=build_scripts)
)
