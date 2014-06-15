"""
Type of the configuration data for the package manager
"""


import os

from sage_pkg.config_yaml import ConfigYAML

class ManagerConfig(ConfigYAML):
    """
    The Configuration for the Sage Package Manager itself

    EXAMPLES::

        >>> config    # doctest: +ELLIPSIS
        Configuration:
        - config.mirrors.extra = ['http://download.example.com/packages']
        - config.mirrors.first = ['http://try.this.first/packages']
        - config.mirrors.list = None
        - config.path.build = /.../test_data/local/var/tmp/sage-build
        - config.path.configuration = []
        - config.path.dot_git = /.../.git
        - config.path.download_cache = /.../upstream
        - config.path.install = /.../test_data/local
        - config.path.install_metadata = /.../test_data/local/var/lib/sage/pkgs
        - config.path.packages = /.../build/manager/test_data/test_pkgs
        - config.path.root = /...
        - config.path.sage_pkg = /.../build/manager
        - config.path.test_data = /.../build/manager/test_data
    """
    def __init__(self, *args, **kwds):
        self._base_dir = kwds.pop('base_dir', '.')
        super(ManagerConfig, self).__init__(*args, **kwds)

    def _normalize(self, config):
        super(ManagerConfig, self)._normalize(config)
        path = config.get('path', dict())
        self._normalize_paths(path)

    def _validate(self):
        super(ManagerConfig, self)._validate()
        self._require(
            'path.configuration',
            'path.dot_git',
            'path.packages',
            'path.sage_pkg',
            'path.download_cache',
            'mirrors.list',
        )

    def _normalize_paths(self, path, root=None):
        if root is None:
            root = path['root']
            root = os.path.abspath(os.path.join(self._base_dir, root))
            path['root'] = root
    
        def normalize(value):
            value = os.path.expanduser(value)
            if not os.path.isabs(value):
                value = os.path.join(root, value)
            return os.path.abspath(value)

        for key in path.keys():
            value = path[key]
            if isinstance(value, list):
                path[key] = map(normalize, value)
            elif isinstance(value, dict):
                self._normalize_paths(value, root)
            else:
                path[key] = normalize(value)
       
 
