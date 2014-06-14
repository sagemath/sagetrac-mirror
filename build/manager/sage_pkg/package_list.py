"""
Package Data Reader
"""

import os
import yaml

from sage_pkg.config import config
from sage_pkg.package.package_config import PackageConfig
from sage_pkg.git.repo import GitRepository
from sage_pkg.task_queue import TaskQueue


class InvalidPackage(Exception):
    pass
    

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
        - config.builder.type = TestPackage
        - config.category = standard
        - config.name = foo
        - config.source.tarball.name = Foo-1.3.tar.gz
        - config.source.tarball.sha1 = 99533b2233a27132ee65bcd2992f7bb464388b55
        - config.source.tarball.url = http://www.example.com/1.3.tar.gz
        - config.source.version = 1.3
    """
    package_yaml = os.path.join(config.path.packages, pkg_name, 'package.yaml')
    if not os.path.isfile(package_yaml):
        raise InvalidPackage('missing package.yaml')
    return PackageConfig(package_yaml)


class PackageLoader(object):

    def _make_random_sha1(self):
        """
        EXAMPLES::

            >>> sha1 = loader._make_random_sha1()
            >>> sha1  # doctest: +SKIP
            '3835bb5604b33160a94f47ee8d4262b9471c0017'
            >>> len(sha1) == 40 and int(sha1, 16) >= 0
            True
        """
        import datetime
        now = str(datetime.datetime.utcnow())
        import sha
        return sha.sha(now).hexdigest()
    
    def _make_package(self, pkg_config, version):
        from sage_pkg.package import all
        try:
            cls = getattr(all, pkg_config.builder.type)
        except AttributeError:
            raise ValueError('unknown builder type: ' + pkg_config.builder.type)
        return cls(pkg_config, version)
        
    
    def get_all(self):
        """
        Return the list of all packages
    
        EXAMPLES::

            >>> loader.get_all()
            (bar, baz, foo)
        """
        if hasattr(self, '_packages'):
            return self._packages
        git = GitRepository(config.path.dot_git)
        tree = git.get_tree(config.path.packages)
        git_version = dict()
        for mode, name, sha1 in git.get_tree(config.path.packages).ls_dirs():
            pkg_dir_name = os.path.join(config.path.packages, name)
            if git.is_clean_dir(pkg_dir_name):
                git_version[name] = sha1
        # Packages that are not clean or not checked in get a random version stamp to force rebuild
        random_version = self._make_random_sha1()
        result = []
        for name in os.listdir(config.path.packages):
            fullname = os.path.join(config.path.packages, name)
            if not os.path.isdir(fullname):
                continue
            try:
                pkg_config = load_config(name)
            except InvalidPackage as err:
                print('Package {0} invalid: {1}'.format(name, err))
                continue
            if pkg_config.name.lower() != name.lower():
                raise ValueError('The name in package.yaml must match the directory name')
            try:
                version = git_version[name]
            except KeyError:
                print('Package {0} has changes, forcing rebuild'.format(name))
                version = random_version
            pkg = self._make_package(pkg_config, version)
            result.append(pkg)
        self._packages = tuple(sorted(result))
        return self._packages
    
    def get(self, name):
        """
        Load a single package
        
        This is just for doctesting.
        """
        for pkg in self.get_all():
            if pkg.name == name:
                return pkg
        raise ValueError('unknown package: ' + name)

    def queue(self):
        """
        Return packages as queue

        The packages don't have any work associated (i.e. they are not
        callable). So you can only use it the output to iterate
        through the packages (in a way that satisfies dependencies).

        EXAMPLES::

            >>> queue = loader.queue()
            >>> while not queue.is_finished():
            ...     task = queue.next(reproducible=True)
            ...     print task
            ...     queue.finish(task)
            Task foo
            Task bar
            Task baz
        """
        pkg_dict = dict((pkg.name, pkg) for pkg in self.get_all())
        queue = TaskQueue()
        for pkg in self.get_all():
            queue.add_work(pkg, [pkg_dict[dep] for dep in pkg.get_all_dependencies()])
        return queue

    def build_queue(self):
        """
        Return the queue of build tasks

        The build task queue is a refinement of the :meth:`queue` of
        packages. Each package is broken up into different build
        steps.

        EXAMPLES::

            >>> build = loader.build_queue()
            >>> build.run_serial(reproducible=True)    # doctest: +ELLIPSIS
            bar: downloading
            bar: unpacking
            bar: preparing
            baz: downloading
            ...
            baz: installing
        """
        dependencies = dict()
        packages = self.queue()
        queue = TaskQueue()
        while not packages.is_finished():
            pkg = packages.next()
            task_list = pkg.work.build_tasks(dependencies)
            queue.add_task(*task_list)
            packages.finish(pkg)
        return queue


loader = PackageLoader()
