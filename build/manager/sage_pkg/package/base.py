"""
Base class for Packages

Packages define 7 steps:

* Download of the sources (including verification of checksums)

* Unpacking of the sources

* Preparation of the source (typically apply patches)

* Configuring (typically run ``./configure``)

* Compile (typically run ``make``)

* Check (typically run ``make check``)

* Installation (typically run ``make install``). Note: this is the
  only step that is allowed to put files in the destination directory.

Each of these step is run after (and only if) the previous step finished
successfully, that is, did not raise an exception.
"""


import os
import copy
import shutil
import tempfile
import yaml
import datetime

from sage_pkg.config import config
from sage_pkg.config_yaml import ConfigYAML
from sage_pkg.logger import logger
from sage_pkg.task_queue import TaskQueue
from sage_pkg.utils import cached_property, pretty_age
from .dependency_accumulator import CommutingSha1Accumulator


class PackageBase(object):
    
    def __init__(self, app_config, version_stamp=None, longterm=False):
        """
        INPUT:

        - ``app_config`` -- instance of
          :class:`~sage_pkg.package.package_config.PackageConfig`. The
          package configuration.

        - ``version_stamp`` -- string. Uniquely characterizes
          ``app_config``.

        - ``longterm`` -- boolean. Whether the same ``version_stamp``
          can possibly reappear in a future compilation. If true, the
          compiled package might be cached.
        """
        self._config = app_config
        self._version_stamp = version_stamp
        self._longterm = longterm
        self._metadata = None

    def _init_dependencies(self, pkg_dict):
        """
        Called after loading all packages

        In a second initialization step, we replace dependency names
        with the actual package objects. As a side effect, this checks
        for typos in the dependencies.

        INPUT:

        - ``pkg_dict`` -- a dictionary whose keys are all package
          namse and values the associated :class:`PackageBase`
          instance.
        """
        def deps(dep_type):
            try:
                dependencies = self._config('depends', dep_type)
            except KeyError:
                return tuple()
            logger.debug('initializating %s dependencies for %s: %s', dep_type, self.name, dependencies)
            if dependencies:
                try:
                    return tuple(pkg_dict[dep] for dep in dependencies)
                except KeyError as err:
                    raise ValueError('invalid {0} dependency in {1}: {2}'
                                     .format(dep_type, self.name, err))
            else:
                return tuple()
        self._hard_dependencies = deps('hard')
        self._soft_dependencies = deps('soft')
        self._test_dependencies = deps('test')
        self._build_dependencies = deps('build')

    def __eq__(self, other):
        return self._config.name == other._config.name

    def __lt__(self, other):
        return self._config.name < other._config.name

    def __repr__(self):
        return self.name

    @property
    def config(self):
        """
        The package configuration (read from ``package.yaml``)

        EXAMPLES::
  
            >>> loader.get('foo').config
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
        return self._config

    @property
    def name(self):
        """
        Package name

        EXAMPLES::
  
            >>> loader.get('foo').name
            'foo'
        """
        return self._config.name

    @property
    def category(self):
        """
        Return the package category
        """
        return self._config.category

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
        return 'unknown'

    @property
    def version_stamp(self):
        """
        Unique machine-readable version.

        This is the version relevant for dependency
        calculations. Should not be used in the UI.

        OUTPUT:

        String (ASCII SHA1).

        EXAMPLES::
  
            >>> loader.get('foo').version_stamp    # doctest: +SKIP
            'a71d5accb4ed818985dfcd796090b75afb83885c'
        """
        return self._version_stamp

    @cached_property
    def dependency_stamp(self):
        """
        Dependency version.

        This versions the dependency, that is: If and only if
        this property changes then the package must be recompiled.

        It combines:

        * Our own :meth:`version_stamp`

        * The :meth:`dependency_stamp` of the hard and soft
          dependencies.

        Hence, it recursively includes the version stamps of all
        dependencies.

        Note that the test and build dependencies are not included; If
        you managed to successfully build/test the package with older
        versions then there is no need to recompile.

        OUTPUT:

        String (ASCII SHA1).

        EXAMPLES::

            >>> loader.get('foo').dependency_stamp    # doctest: +SKIP
            'a71d5accb4ed818985dfcd796090b75afb83885c'
        """
        acc = CommutingSha1Accumulator()
        acc += self.version_stamp
        for dep in self.get_hard_dependencies() + self.get_soft_dependencies():
            acc += dep.dependency_stamp
        return str(acc)

    @property
    def is_up_to_date(self):
        """
        Whether the package needs to be built (or rebuilt).

        OUTPUT:

        Boolean

            >>> loader.get('foo').is_up_to_date
            False
        """
        return self.dependency_stamp == self.metadata.dependency_stamp

    @cached_property
    def status(self):
        """
        Human-readable package status.

        OUTPUT:

        String. Indicates whether the package is up to date or not.
        """
        if not self.metadata.previous:
            return 'not installed'
        elif self.dependency_stamp == self.metadata.dependency_stamp:
            return 'up to date'
        elif self.version_stamp == self.metadata.version_stamp:
            return 'dependencies changed'
        else:
            return 'old version'

    @cached_property
    def age(self):
        if not self.metadata.previous:
            return 'NA'
        return pretty_age(self.metadata.time)
        
    @cached_property
    def build_dir(self):
        """
        The build directory. 
        
        The unpack step must create/overwrit this directory with the
        source tree.

        EXAMPLES::

            >>> loader.get('foo').build_dir   # doctest: +ELLIPSIS
            '/.../test_data/local/var/tmp/sage-build/foo'
        """
        return os.path.join(config.path.build, self.name)

    def get_all_dependencies(self):
        """
        Return all direct dependencies.

        EXAMPLES::

            >>> loader.get('baz').get_all_dependencies()    
            (bar,)
        """
        return \
            self.get_hard_dependencies() + \
            self.get_test_dependencies() + \
            self.get_build_dependencies()

    def get_hard_dependencies(self):
        return self._hard_dependencies
        
    def get_soft_dependencies(self):
        return self._soft_dependencies
        
    def get_build_dependencies(self):
        return self._build_dependencies

    def get_test_dependencies(self):
        """
        Return the dependencies for running self-tests
        """
        return self._test_dependencies

    def want_check(self):
        try:
            return config('check', self.name)
        except KeyError:
            return config('check', 'default', default=False)
        
    def get_environment(self):
        """
        Get a dictionary of environment variables

        EXAMPLES::

            >>> foo = loader.get('foo')
            >>> foo.get_environment()   # doctest: +ELLIPSIS
            {'...'}
        """
        return copy.copy(os.environ)

    def build_tasks(self, dependencies, stop_at='install'):
        """
        Each build step is its own task to be finer grained.

        Note that we try to push dependencies as far back as
        possible. Download works without anything, unpack and prepare
        only needs build tools (like bzip2 and GNU patch).

        INPUT:

        - ``dependencies`` -- dict. Which build task (= value) satisfies which
          package dependency (= key).

        OUTPUT:

        A list of build tasks.

        EXAMPLES::

            >>> foo = loader.get('foo')
            >>> deps = dict()
            >>> tasks = foo.build_tasks(deps);  tasks
            [foo-download, foo-unpack, foo-prepare, foo-configure, foo-compile, foo-install]
            >>> for task in tasks:  print task.work
            <bound method TestPackage.download of foo>
            <bound method TestPackage.unpack of foo>
            <bound method TestPackage.prepare of foo>
            <bound method TestPackage.configure of foo>
            <bound method TestPackage.compile of foo>
            <bound method TestPackage.install of foo>

        The ``deps`` is modified to record what task satisifes the
        dependency on the package name::

            >>> deps
            {foo: foo-install}

        You can also stop before installation, in which case ``deps``
        is never modified::

            >>> foo.build_tasks(deps, stop_at='download')
            [foo-download]
            >>> foo.build_tasks(deps, stop_at='unpack')
            [foo-download, foo-unpack]
            >>> foo.build_tasks(deps, stop_at='prepare')
            [foo-download, foo-unpack, foo-prepare]

            >>> foo.build_tasks(deps, stop_at='unknownvalue')
            Traceback (most recent call last):
            ...
            ValueError: unknown value for the stop_at parameter
        """
        from sage_pkg.task_queue import Task
        dependencies[self] = None   # default if you skip installation
        tasklist = []
        task_download  = prev = Task(self.download,  [],     name='{0}-download' .format(self.name))
        tasklist.append(prev)
        if stop_at == 'download':
            return tasklist
        build = [prev] + [dependencies[dep] for dep in self.get_build_dependencies()]
        task_unpack    = prev = Task(self.unpack,    build,  name='{0}-unpack'   .format(self.name))
        tasklist.append(prev)
        if stop_at == 'unpack':
            return tasklist
        task_prepare   = prev = Task(self.prepare,   [prev], name='{0}-prepare'  .format(self.name))
        tasklist.append(prev)
        if stop_at == 'prepare':
            return tasklist
        hard = [prev] + [dependencies[dep] for dep in self.get_hard_dependencies()]
        task_configure = prev = Task(self.configure, hard,   name='{0}-configure'.format(self.name))
        tasklist.append(prev)
        if stop_at == 'configure':
            return tasklist
        task_compile   = prev = Task(self.compile,   [prev], name='{0}-compile'  .format(self.name))
        tasklist.append(prev)
        if stop_at == 'compile':
            return tasklist
        if self.want_check():
            test = [prev] + [dependencies[dep] for dep in self.get_test_dependencies()]
            task_check     = prev = Task(self.check,     test,   name='{0}-check'    .format(self.name))
            tasklist.append(prev)
            if stop_at == 'check':
                return tasklist
        task_install   = prev = Task(self.install,   [prev], name='{0}-install'  .format(self.name))
        tasklist.append(prev)
        if stop_at != 'install':
            raise ValueError('unknown value for the stop_at parameter')
        dependencies[self] = task_install
        return tasklist

    def build_queue(self, stop_at='install'):
        """
        Return the queue of build tasks for this package

        This method is analogous to
        :meth:`~sage_pkg.package_list.PackageLoader.build_queue`
        except that it only contains build tasks for this one
        package. Inter-package dependencies are ignored.

        EXAMPLES::

            >>> pkg = loader.get('baz')
            >>> pkg.get_all_dependencies()
            (bar,)
            >>> queue = pkg.build_queue()
            >>> queue.run_serial()
            baz: downloading
            baz: unpacking
            baz: preparing
            baz: configuring
            baz: installing
        """
        # Ignore dependencies
        dependencies = dict()
        for dep in self.get_all_dependencies():
            dependencies[dep] = None
        queue = TaskQueue()
        tasks = self.build_tasks(dependencies, stop_at=stop_at)
        queue.add_task(*tasks)
        return queue

    def _success(self):
        """
        Installation finished successfully.

        This method is called when the installation finished
        successfully as part of the install step.
        """
        result = dict(
            version=self.version,
            version_stamp=self.version_stamp,
            dependency_stamp=self.dependency_stamp,
            time=datetime.datetime.utcnow(),
            install_dir=config.path.install,
        )
        self._save_metadata(result)
        
    @cached_property
    def _metadata_filename(self):
        """
        Where the package installation metadata is saved
        """
        return os.path.join(config.path.install_metadata, self.name + '.yaml')

    @cached_property
    def metadata(self):
        """
        Return the metadata of previous installation attempts.
        """
        defaults = dict(
            filename=self._metadata_filename,
            version_stamp='0'*40,
            dependency_stamp='0'*40,
        )
        if os.path.exists(self._metadata_filename):
            return ConfigYAML(defaults, dict(previous=True), self._metadata_filename, name='metadata')
        else:
            return ConfigYAML(defaults, dict(previous=False), name='metadata')
            
    def _save_metadata(self, value_dict):
        """
        Save the metadata.

        This should only be called after successful installation. The
        new metadata will then be found the next time we run the
        package manager.

        EXAMPLES::

            >>> pkg = loader.get('baz')
            >>> pkg.metadata    # doctest: +ELLIPSIS
            Configuration:
            - metadata.dependency_stamp = 0000000000000000000000000000000000000000
            - metadata.filename = /.../test_data/local/var/lib/sage/pkgs/baz.yaml
            - metadata.previous = False
            - metadata.version_stamp = 0000000000000000000000000000000000000000
            >>> pkg._metadata_filename    # doctest: +ELLIPSIS
            '/.../test_data/local/var/lib/sage/pkgs/baz.yaml'
        """
        metadata = yaml.dump(value_dict, default_flow_style=False)
        # atomically save result
        tmpfile = tempfile.NamedTemporaryFile(dir=config.path.install_metadata, mode='w', delete=False)
        try:
            tmpfile.write(metadata)
            tmpfile.close()
            os.rename(tmpfile.name, self._metadata_filename)
        except BaseException as error:
            try:
                tmpfile.close()
                os.remove(tmpfile.name)
            except OSError:
                pass
            raise error
        logger.info('saved %s metadata at %s', self.name, self._metadata_filename)

    def download(self):
        """
        Download of the sources.

        This includes verification of checksums, if any.

        To be implemented in a derived class.
        """
        pass

    def unpack(self):
        """
        Unpacking of the sources.

        To be implemented in a derived class
        """
        pass

    def prepare(self):
        """
        Preparation of the source.

        To be implemented in a derived class.

        Typically, this step applies patches.
        """
        pass

    def configure(self):
        """
        Configuring.

        To be implemented in a derived class.

        Typically, this runs ``./configure``.
        """
        pass

    def compile(self):
        """
        Compile.
        
        To be implemented in a derived class

        Typically, this runs ``make``.
        """
        pass

    def check(self):
        """
        Run self-tests
        
        To be implemented in a derived class

        Typically, this runs ``make check``.
        """
        pass

    def install(self):
        """
        Installation

        To be implemented in a derived class. Make sure to call the
        superclass :meth:`install` at the end.

        Typically, this runs ``make install``. 

        .. warning:: 

            This is the only step that is allowed to put files in the
            installation directory.
        """
        self._success()
