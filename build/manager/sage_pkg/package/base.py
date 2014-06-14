"""
Base class for Packages

Packages define 7 steps:

* Download of the sources

* Unpacking of the sources

* Preparation of the source (typically apply patches)

* Configuring (typically run ``./configure``)

* Build (typically run ``make``)

* Check (typically run ``make check``)

* Installation (typically run ``make install``). Note: this is the
  only step that is allowed to put files in the destination directory.

Each of these step is run after (and only if) the previous step finished
successfully, that is, did not raise an exception.
"""


import os

from sage_pkg.config import config


class PackageBase(object):
    
    def __init__(self, app_config, version_stamp=None):
        """
        INPUT:

        - ``app_config`` -- instance of
          :class:`~sage_pkg.package.package_config.PackageConfig`. The
          package configuration.

        - ``version_stamp`` -- string that uniquely characterizes
          ``app_config``.
        """
        self._config = app_config
        self._version_stamp = version_stamp

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
        return self.get_hard_dependencies() + self.get_build_dependencies()

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
        
    def build_tasks(self, dependencies):
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
            >>> foo.build_tasks(deps)
            [foo-download, foo-unpack, foo-prepare, foo-configure, foo-build, foo-check, foo-install]

        The ``deps`` is modified to record what task satisifes the
        dependency on the package name::

            >>> deps
            {'foo': foo-install}
        """
        from sage_pkg.task_queue import Task
        task_download  = prev = Task(self.download,  [],     name='{0}-download' .format(self.name))
        build = [prev] + [dependencies[dep] for dep in self.get_build_dependencies()]
        task_unpack    = prev = Task(self.unpack,    build,  name='{0}-unpack'   .format(self.name))
        task_prepare   = prev = Task(self.prepare,   [prev], name='{0}-prepare'  .format(self.name))
        hard = [prev] + [dependencies[dep] for dep in self.get_hard_dependencies()]
        task_configure = prev = Task(self.configure, hard,   name='{0}-configure'.format(self.name))
        task_build     = prev = Task(self.build,     [prev], name='{0}-build'    .format(self.name))
        task_check     = prev = Task(self.check,     [prev], name='{0}-check'    .format(self.name))
        task_install   = prev = Task(self.install,   [prev], name='{0}-install'  .format(self.name))
        dependencies[self.name] = task_install
        return [task_download, task_unpack, task_prepare, 
                task_configure, task_build, task_check, task_install]

    def download(self):
        pass

    def unpack(self):
        pass

    def prepare(self):
        pass

    def configure(self):
        pass

    def build(self):
        pass

    def check(self):
        pass

    def install(self):
        pass
