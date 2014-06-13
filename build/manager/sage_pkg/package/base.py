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
    
    def __init__(self, app_config, version=None):
        """
        INPUT:

        - ``app_config`` -- instance of
          :class:`~sage_pkg.package.package_config.PackageConfig`. The
          package configuration.

        - ``version`` -- string that uniquely characterizes
          ``app_config``.
        """
        self._verify(app_config)
        self._config = app_config
        self._version = version

    def __lt__(self, other):
        return self._config.name < other._config.name

    def __repr__(self):
        return repr(self._config)

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
