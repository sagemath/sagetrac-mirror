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
        - config.mirrors = ['http://download.example.com/packages']
        - config.path.build = /var/tmp/sage-build
        - config.path.configuration = []
        - config.path.dot_git = /.../.git
        - config.path.download_cache = /.../upstream
        - config.path.install = /usr/local
        - config.path.packages = /.../build/manager/test_pkgs
        - config.path.root = /...
        - config.path.sage_pkg = /.../build/manager
    """
    
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
        )

    def _normalize_paths(self, path, root=None):
        if root is None:
            root = path['root']
            root = os.path.abspath(root)
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
       
 
