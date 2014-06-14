"""
Builder for spkg-install scripts
"""

import os
import subprocess

from sage_pkg.config import config
from sage_pkg.chdir_context import chdir
from .base import PackageBase
from .sage_mirror_mixin import SageMirrorMixin


class SpkgInstallScript(PackageBase, SageMirrorMixin):

    @property
    def install_script(self):
        scriptname = self._config.builder.install_script
        return os.path.join(config.path.packages, self.name, scriptname)

    @property
    def check_script(self):
        try:
            scriptname = self._config.builder.check_script
        except AttributeError:
            return None
        return os.path.join(config.path.packages, self.name, scriptname)

    def install(self):
        with chdir(self.build_dir):
            subprocess.check_call(self.install_script)

    def check(self):
        if not self.check_script:
            return
        with chdir(self.build_dir):
            subprocess.check_call(self.check_script)

