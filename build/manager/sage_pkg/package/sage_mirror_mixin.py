"""
Mixin to download from the Sage mirror network
"""

import os
import shutil
import subprocess

from sage_pkg.config import config
from sage_pkg.chdir_context import chdir


class SageMirrorMixin(object):

    @property
    def tarball(self):
        """
        The tarball name (how upstream names their stuff)
        """
        return self._config.source.tarball

    @property
    def tarball_cache(self):
        return os.path.join(config.path.download_cache, self._config.source.tarball)

    def tarball_url_iter(self):
        try:
            mirrors = config.mirrors
        except AttributeError:
            raise ValueError('you need to define at least one mirror')
        for url in mirrors:
            url = url.rstrip('/')
            yield '/'.join(url, self.name, self.tarball)

    def download(self):
        pass

    def unpack(self):
        build_dir = self.build_dir
        if os.path.exists(build_dir):
            shutil.rmtree(build_dir)
        os.makedirs(build_dir)
        cmd = ['tar', '-xf', self.tarball_cache, '-C', build_dir]
        subprocess.check_call(cmd)

        
