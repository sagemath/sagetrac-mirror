import os
import sys
import glob
import shutil
import sysconfig

from setuptools import setup
from distutils.command.build_scripts import build_scripts as distutils_build_scripts
from setuptools.command.build_py import build_py as setuptools_build_py
from distutils.errors import (DistutilsSetupError, DistutilsModuleError,
                              DistutilsOptionError)

class build_py(setuptools_build_py):

    def run(self):
        HERE = os.path.dirname(__file__)
        with open(os.path.join(HERE, 'VERSION.txt')) as f:
            sage_version = f.read().strip()
        SETENV = ':'
        # Until pynac is repackaged as a pip-installable package (#30534), SAGE_LOCAL still has to be specific to
        # the Python version.  Note that as of pynac-0.7.26.sage-2020-04-03, on Cygwin, pynac is linked through
        # to libpython; whereas on all other platforms, it is not linked through, so we only key it to the SOABI.
        soabi = sysconfig.get_config_var('SOABI')
        if sys.platform == 'cygwin':
            libdir_tag = sysconfig.get_config_var('LIBDIR').replace(' ', '-').replace('\\', '-').replace('/', '-')
            ldversion = sysconfig.get_config_var('LDVERSION')
            python_tag = f'{libdir_tag}-{ldversion}'
        else:
            python_tag = soabi

        SAGE_ROOT = os.path.join('/var/tmp', f'sage-{sage_version}-{python_tag}')
        SAGE_LOCAL = os.path.join(SAGE_ROOT, 'local')
        if os.path.exists(os.path.join(SAGE_ROOT, 'config.status')):
            print('Reusing SAGE_ROOT={SAGE_ROOT}')
        else:
            # config.status and other configure output has to be writable.
            # So (until the Sage distribution supports VPATH builds - #21469), we have to make a copy of sage_root.
            try:
                shutil.copytree(os.path.join(HERE, 'sage_root_source'), SAGE_ROOT)  # will fail if already exists
            except Exception:
                raise DistutilsSetupError(f"the directory SAGE_ROOT={SAGE_ROOT} already exists but it is not configured.  Please remove it and try again.")
            # Use our copy of the sage_conf template, which contains the relocation logic
            shutil.copyfile(os.path.join(HERE, 'sage_conf.py.in'),
                            os.path.join(SAGE_ROOT, 'build', 'pkgs', 'sage_conf', 'src', 'sage_conf.py.in'))

            cmd = f"cd {SAGE_ROOT} && {SETENV} && ./configure --prefix={SAGE_LOCAL} --with-python={sys.executable} --with-system-python3=force --without-system-gmp --without-system-mpfr"
            print(f"Running {cmd}")
            if os.system(cmd) != 0:
                raise DistutilsSetupError("configure failed")

        # build-local only builds the non-Python packages of the Sage distribution.
        # It still makes an (empty) venv in SAGE_LOCAL, which is unused by default;
        # but a user could use "make build-venv" to build compatible wheels for all Python packages.
        # TODO: A target to only build wheels of tricky packages
        # (that use native libraries shared with other packages).
        SETMAKE = 'if [ -z "$MAKE" ]; then export MAKE="make -j$(PATH=build/bin:$PATH build/bin/sage-build-num-threads | cut -d" " -f 2)"; fi'
        TARGETS = 'build-local'
        #TARGETS = 'base-toolchain giac'
        cmd = f'cd {SAGE_ROOT} && {SETENV} && {SETMAKE} && $MAKE V=0 {TARGETS}'
        if os.system(cmd) != 0:
            raise DistutilsSetupError("make build-local failed")

        # Install configuration
        shutil.copyfile(os.path.join(SAGE_ROOT, 'build', 'pkgs', 'sage_conf', 'src', 'sage_conf.py'),
                        os.path.join(HERE, 'sage_conf.py'))
        if not self.distribution.py_modules:
            self.py_modules = self.distribution.py_modules = []
        self.distribution.py_modules.append('sage_conf')
        shutil.copyfile(os.path.join(SAGE_ROOT, 'src', 'bin', 'sage-env-config'),
                        os.path.join(HERE, 'bin', 'sage-env-config'))
        # Install built SAGE_ROOT as package data
        if not self.packages:
            self.packages = self.distribution.packages = ['']
        if not self.distribution.package_data:
            self.package_data = self.distribution.package_data = {}
        os.symlink(SAGE_ROOT, os.path.join(HERE, 'sage_root'))
        # We do not include lib64 (a symlink) because all symlinks are followed,
        # causing another copy to be installed.
        self.distribution.package_data[''] = (
            glob.glob('sage_root/*')
            + glob.glob('sage_root/build/**', recursive=True)
            + glob.glob('sage_root/local/*')
            + glob.glob('sage_root/local/bin/**', recursive=True)
            + glob.glob('sage_root/local/include/**', recursive=True)
            + glob.glob('sage_root/local/lib/**', recursive=True)
            + glob.glob('sage_root/local/share/**', recursive=True)
            + glob.glob('sage_root/local/var/**', recursive=True)
            )
        #
        setuptools_build_py.run(self)

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
    cmdclass=dict(build_py=build_py, build_scripts=build_scripts),
    # Do not mark the wheel as pure
    has_ext_modules=lambda: True
)
