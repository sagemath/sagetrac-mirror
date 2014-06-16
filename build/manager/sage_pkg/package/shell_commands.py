"""
Builder for Inline Shell Commands
"""

import os
import subprocess

from sage_pkg.config import config
from sage_pkg.chdir_context import chdir
from sage_pkg.utils import cached_property
from sage_pkg.git.repo import GitRepository
from sage_pkg.logger import logger

from .base import PackageBase
from .source_version_stamp_mixin import SourceVersionStampMixin
from .sage_environment_mixin import SageEnvironmentMixin


class ShellCommands(SourceVersionStampMixin, SageEnvironmentMixin, PackageBase):

    def _validate(self):
        super(ShellCommands, self)._validate()
        self._config._require(
            'builder.workdir',
            'builder.commands',
        )

    @property
    def version(self):
        return '== Sage Version'

    def install(self):
        env = self.get_environment()
        builder = self._config.builder
        work_dir = os.path.join(config.path.root, builder.workdir)
        for cmd in builder.commands:
            logger.debug('%s: executing %s', self.name, cmd)
            subprocess.check_call(cmd, env=env, shell=True, cwd=work_dir)
        super(ShellCommands, self).install()
