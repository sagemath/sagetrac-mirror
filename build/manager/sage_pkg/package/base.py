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

from sage_pkg.config import config
from sage_pkg.task_queue import TaskQueue


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

    def __eq__(self, other):
        return self._config.name == other._config.name

    def __lt__(self, other):
        return self._config.name < other._config.name

    def __repr__(self):
        return self.name

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
    def version(self):
        """
        Human-readable version.

        May only be used for the UI, but never for dependency
        calculations.

        EXAMPLES::
  
            >>> loader.get('foo').version
            1.3
        """
        return self._config.source.version

    @property
    def version_stamp(self):
        """
        Unique machine-readable version.

        This is the version relevant for dependency
        calculations. Should not be used in the UI.

        EXAMPLES::
  
            >>> loader.get('foo').version_stamp    # doctest: +SKIP
            'a71d5accb4ed818985dfcd796090b75afb83885c'
        """
        return self._version_stamp

    @property
    def build_dir(self):
        """
        The build directory. 
        
        The unpack step must create/overwrit this directory with the
        source tree.

        EXAMPLES::

            >>> loader.get('foo').build_dir
            '/var/tmp/sage-build/foo'
        """
        return os.path.join(config.path.build, self.name)

    def get_all_dependencies(self):
        """
        Return all direct dependencies.

        EXAMPLES::

            >>> loader.get('baz').get_all_dependencies()    
            ['bar']
        """
        return \
            self.get_hard_dependencies() + \
            self.get_test_dependencies() + \
            self.get_build_dependencies()

    def get_hard_dependencies(self):
        try:
            return self._config.depends.hard
        except AttributeError:
            return []
        
    def get_build_dependencies(self):
        try:
            return self._config.depends.build
        except AttributeError:
            return []

    def get_test_dependencies(self):
        """
        Return the dependencies for running self-tests
        """
        try:
            return self._config.depends.test
        except AttributeError:
            return []
        
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

        - ``dependencies`` -- dict. Which build task satisfies which
          package dependency.

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
            {'foo': foo-install}

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
        dependencies[self.name] = None   # default if you skip installation
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
        dependencies[self.name] = task_install
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
            ['bar']
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

        To be implemented in a derived class.

        Typically, this runs ``make install``. 

        .. warning:: 

            This is the only
            step that is allowed to put files in the installation directory.
        """
        pass
