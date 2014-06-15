"""
Builder for spkg-install scripts
"""

import os
import shutil
import subprocess

from sage_pkg.config import config
from sage_pkg.chdir_context import chdir
from sage_pkg.utils import cached_property
from .base import PackageBase
from .sage_mirror_mixin import SageMirrorMixin
from .sage_environment_mixin import SageEnvironmentMixin



class SpkgInstallScript(SageEnvironmentMixin, SageMirrorMixin, PackageBase):

    @cached_property
    def install_script(self):
        scriptname = self._config.builder.install_script
        return os.path.join(config.path.packages, self.name, scriptname)

    @cached_property
    def check_script(self):
        try:
            scriptname = self._config.builder.check_script
        except AttributeError:
            return None
        return os.path.join(config.path.packages, self.name, scriptname)

    def get_environment(self):
        """
        Set up the build-time environment for Sage packages
        """
        env = super(SpkgInstallScript, self).get_environment()
        env['PKG_NAME'] = self.name
        return env

    def prepare(self):
        """
        Copy spkg-install to the build dir
        """
        self._copy_spkg_install()
        self._copy_patches()
        super(SpkgInstallScript, self).prepare()

    def _copy_spkg_install(self):
        spkg_install = os.path.join(self.build_dir, 'spkg-install')
        shutil.copy(self.install_script, spkg_install)
        os.chmod(spkg_install, 0o755)

    def _copy_patches(self):
        patches = os.path.join(config.path.packages, self.name, 'patches')
        if os.path.isdir(patches):
            dst = os.path.join(self.build_dir, 'patches')
            shutil.copytree(patches, dst)

    def install(self):        
        env = self.get_environment()
        with chdir(self.build_dir):
            subprocess.check_call('sage-spkg-sage_pkg', env=env, shell=True)
        super(SpkgInstallScript, self).install()

    def check(self):
        if not (self.want_check() and self.check_script):
            return
        env = self.get_environment()
        with chdir(self.build_dir):
            subprocess.check_call(self.check_script)
        super(SpkgInstallScript, self).check()

