"""
Package Data Reader
"""

import os
import yaml

from sage_pkg.config import config
from sage_pkg.package.package_config import PackageConfig
from sage_pkg.tiny_git import GitRepository


def load_config(pkg_name):
    """
    Return the configuration for a package

    INPUT:

    - ``pkg_name`` -- string. Name of a subdirectory of the packages
      directory.

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
    app_yaml = os.path.join(config.path.packages, pkg_name, 'app.yaml')
    return PackageConfig(app_yaml)


def make_random_sha1():
    """
    EXAMPLES::

        >>> from sage_pkg.package_list import make_random_sha1
        >>> sha1 = make_random_sha1()
        >>> sha1  # doctest: +SKIP
        '3835bb5604b33160a94f47ee8d4262b9471c0017'
        >>> len(sha1) == 40 and int(sha1, 16) >= 0
        True
    """
    import datetime
    now = str(datetime.datetime.utcnow())
    import sha
    return sha.sha(now).hexdigest()


def make_package(config, version):
    from sage_pkg.packages import all
    try:
        cls = getattr(all, config.builder.type)
    except AttributeError:
        raise ValueError('unknown builder type: ' + config.builder.type)
    return cls(config, version)
    

def load_packages():
    """
    Return the list of all packages
    """
    git = GitRepository(config.path.dot_git)
    tree = git.get_tree(config.path.packages)
    git_version = dict()
    for mode, name, sha1 in git.get_tree(config.path.packages).ls_dirs():
        pkg_tree = git.get(sha1)
        if pkg_tree.is_clean():
            git_version[name] = sha1
    # Packages that are not clean or not checked in get a random version stamp to force rebuild
    random_version = make_random_sha1()
    result = []
    for name in os.listdir(config.path.packages):
        if not os.path.isdir(name):
            continue
        config = load_config(name)
        if config.name.lower() != name.lower():
            raise ValueError('The name in app.yaml must match the directory name')
        try:
            version = git_version[name]
        except KeyError:
            print('Package {0} has changes, forcing rebuild')
            version = random_version
        pkg = make_package(config, version)
        result.append(pkg)
    return sorted(result)
