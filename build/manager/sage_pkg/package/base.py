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

    def all_dependencies(self):
        """
        EXAMPLES::

            >>> loader.get('baz').all_dependencies()    
            ['foo', 'bar']
        """
        result = []
        try:
            depends = self._config.depends
        except AttributeError:
            return result
        for deps in depends._c.values():
            result.extend(deps)
        return result

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
