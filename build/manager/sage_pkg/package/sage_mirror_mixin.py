"""
Mixin to download from the Sage mirror network
"""

import os
import shutil
import subprocess

from sage_pkg.config import config
from sage_pkg.logger import logger
from sage_pkg.chdir_context import chdir
from sage_pkg.utils import cached_property
from sage_pkg.mirror_network import MirrorList, Tarball



class SageMirrorMixin(object):

    def _validate(self):
        super(SageMirrorMixin, self)._validate()
        self._require(
            'source.tarball.name',
            'source.tarball.sha1',
        )

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

    @cached_property
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

    def download(self):
        sha1 = self._config.source.tarball('sha1', default=None)
        tarball = Tarball(self.tarball, sha1)
        if not tarball.is_cached():
            logger.info('downloading tarball %s', self.tarball)
            mirror_list = MirrorList()
            tarball.download(mirror_list)
        super(SageMirrorMixin, self).download()

    def unpack(self):
        build_dir = self.build_dir
        logger.info('unpacking %s to %s', self.tarball, build_dir)
        if os.path.exists(build_dir):
            shutil.rmtree(build_dir)
        os.makedirs(build_dir)
        os.chmod(build_dir, 0o700)
        cmd = ['tar', '-xf', self.tarball_cache, '-C', build_dir, '--no-same-owner']
        subprocess.check_call(cmd)
        self._make_src_symlink()
        super(SageMirrorMixin, self).unpack()

    def _make_src_symlink(self):
        """
        Symlink ``src`` to the top-level directory from within the tarball
        """
        build_dir = self.build_dir
        src = os.path.join(build_dir, 'src')
        if os.path.exists(src):
            logger.info('tarball %s already contains a "src" directory', self.tarball)
            return
        tarball_dir = None
        for dirent in os.listdir(build_dir):
            fullname = os.path.join(build_dir, dirent)
            if not os.path.isdir(fullname):
                continue
            if tarball_dir is not None:
                raise ValueError('ambiguous "src" symlink: multiple directories in %s', self.tarball)
            tarball_dir = dirent
        if tarball_dir is None:
            raise ValueError('no suitable directory in %s as target for the "src" symlink', self.tarball)
        logger.info('symlinking %s to %s', tarball_dir, src)
        os.symlink(tarball_dir, src)
        

