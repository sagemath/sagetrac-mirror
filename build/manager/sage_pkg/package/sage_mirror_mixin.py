"""
Mixin to download from the Sage mirror network
"""

import os
import shutil
import subprocess

from sage_pkg.config import config
from sage_pkg.logger import logger
from sage_pkg.chdir_context import chdir



class SageMirrorMixin(object):

    @classmethod
    def _doctest_example(cls):
        from sage_pkg.package.package_config import PackageConfig
        package_yaml = os.path.join(config.path.test_data, 'spkg', 'package.yaml')
        pc = PackageConfig(package_yaml)
        from sage_pkg.package.all import SpkgInstallScript
        return SpkgInstallScript(pc, 'version', False)

    @property
    def tarball(self):
        """
        The tarball name (how upstream names their stuff)

        EXAMPLES::

            >>> from sage_pkg.package.sage_mirror_mixin import SageMirrorMixin
            >>> pkg = SageMirrorMixin._doctest_example()
            >>> pkg
            spkg_example
            >>> pkg.tarball
            'spkg_example-4.7.4.tar.bz2'
        """
        return self._config.source.tarball.name

    @property
    def tarball_cache(self):
        """
        The tarball cache (file where "download" stores the tarball)

        EXAMPLES::

            >>> from sage_pkg.package.sage_mirror_mixin import SageMirrorMixin
            >>> pkg = SageMirrorMixin._doctest_example()
            >>> pkg.tarball_cache    # doctest: +ELLIPSIS
            '/.../upstream/spkg_example-4.7.4.tar.bz2'
        """
        return os.path.join(config.path.download_cache, self.tarball)

    def tarball_url_iter(self):
        """
        Iterator for potential download URLs

        EXAMPLES::

            >>> from sage_pkg.package.sage_mirror_mixin import SageMirrorMixin
            >>> pkg = SageMirrorMixin._doctest_example()
            >>> list(pkg.tarball_url_iter())
            ['http://download.example.com/packages/spkg_example/spkg_example-4.7.4.tar.bz2']
        """
        try:
            mirrors = config.mirrors
        except AttributeError:
            raise ValueError('you need to define at least one mirror')
        for url in mirrors:
            url = url.rstrip('/')
            yield '/'.join([url, self.name, self.tarball])

    def download(self):
        super(SageMirrorMixin, self).download()

    def unpack(self):
        super(SageMirrorMixin, self).unpack()
        build_dir = self.build_dir
        logger.info('unpacking %s to %s', self.tarball, build_dir)
        if os.path.exists(build_dir):
            shutil.rmtree(build_dir)
        os.makedirs(build_dir)
        os.chmod(build_dir, 0o700)
        cmd = ['tar', '-xf', self.tarball_cache, '-C', build_dir, '--no-same-owner']
        subprocess.check_call(cmd)
        self._make_src_symlink()

    def _make_src_symlink(self):
        """
        Symlink ``src`` to the top-level directory from within the tarball
        """
        build_dir = self.build_dir
        src = os.path.join(build_dir, 'src')
        assert not os.path.exists(src)    # we just deleted the old build_dir
        tarball_dir = None
        for dirent in os.listdir(build_dir):
            fullname = os.path.join(build_dir, dirent)
            if not os.path.isdir(fullname):
                continue
            if dirent.lower().startswith(self.name):
                tarball_dir = dirent
        if tarball_dir is None:
            raise ValueError('no suitable directory as target for the "src" symlink')
        logger.info('symlinking %s to %s', tarball_dir, src)
        os.symlink(tarball_dir, src)
        

