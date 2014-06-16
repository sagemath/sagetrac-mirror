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

    def _validate(self):
        super(SpkgInstallScript, self)._validate()
        self._require(
            'source.version',
            'builder.install_script',
        )

    @property
    def version(self):
        """
        Human-readable version.

        May only be used for the UI, but never for dependency
        calculations. See :meth:`version_stamp`.

        EXAMPLES::
  
            >>> loader.get('foo').version
            1.3
        """
        return self._config.source.version

    @cached_property
    def install_script(self):
        return self._config.builder.install_script

    @cached_property
    def check_script(self):
        return self._config.builder('check_script', default=None)

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
        self._copy_pkg_directory_contents()
        self._make_scripts_executable()
        super(SpkgInstallScript, self).prepare()

    def _copy_pkg_directory_contents(self):
        src = os.path.join(config.path.packages, self.name)
        dst = os.path.join(self.build_dir)
        for dirent in os.listdir(src):
            fqn = os.path.join(src, dirent)
            if os.path.isdir(fqn):
                shutil.copytree(fqn, os.path.join(dst, dirent))
            else:
                shutil.copyfile(fqn, os.path.join(dst, dirent))
                
    def _make_scripts_executable(self):
        def make_executable(script):
            if script is not None:
                fname = os.path.join(self.build_dir, script)
                os.chmod(fname, 0o700)
        make_executable(self.install_script)
        make_executable(self.check_script)

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

