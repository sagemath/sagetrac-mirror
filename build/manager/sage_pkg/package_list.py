"""
Package Data Reader
"""

import os
import yaml

from sage_pkg.config import config
from sage_pkg.package.package_config import PackageConfig
from sage_pkg.tiny_git import GitRepository


def load_config(pkgname):
    """
    Return the configuration for a package

    EXAMPLES::

        >>> from sage_pkg.package_list import load_config
        >>> load_config('foo')
        Configuration:
        - config.builder.install_script = spkg-install
        - config.builder.type = SpkgInstallScript
        - config.category = standard
        - config.name = foo
        - config.source.tarball.name = Foo-1.3.tar.gz
        - config.source.tarball.sha1 = 99533b2233a27132ee65bcd2992f7bb464388b55
        - config.source.tarball.url = http://www.example.com/1.3.tar.gz
        - config.source.version = 1.3
    """
    app_yaml = os.path.join(config.path.packages, pkgname, 'app.yaml')
    return PackageConfig(app_yaml)


def load_packages():
    git = GitRepository(config.path.dot_git)
    for name in os.listdir(package_dir):
        if not os.path.isdir(name):
            continue
        app_yaml = os.path.join(package_dir, name, 'app.yaml')
        if not os.path.isfile(name):
            continue
        yield _load(app_yaml, git)

