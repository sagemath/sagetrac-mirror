"""
Type of the configuration data for individual packges
"""


import os
import copy

from sage_pkg.config_yaml import ConfigYAML

class PackageConfig(ConfigYAML):
    """
    The Configuration for individual packages.

    EXAMPLES::

        >>> app.info('baz')    # doctest: +ELLIPSIS
        Configuration:
        - config.builder.check_script = spkg-check
        - config.builder.install_script = spkg-install
        - config.builder.type = TestPackage
        - config.category = standard
        - config.depends.build = None
        - config.depends.hard = ['bar']
        - config.depends.soft = None
        - config.depends.test = None
        - config.name = baz
        - config.source.git.commit = refs/tags/2.0.1
        - config.source.git.repo = git://github.com/example/baz.git
        - config.source.tarball.name = baz_2.0.1.tar.gz
        - config.source.tarball.sha1 = 99533b2233a27132ee65bcd2992f7bb464388b55
        - config.source.tarball.url = http://www.example.com/2.0.1.tar.gz
        - config.source.version = 2.0.1
        Status: not installed
        Installed when: NA
    """
    
    def _normalize(self, config):
        super(PackageConfig, self)._normalize(config)
        formatting = dict((key, ConfigYAML(value) if isinstance(value, dict) else value) 
                          for key, value in config.items())
        self._normalize_format(config, formatting)

    def _validate(self):
        super(PackageConfig, self)._validate()
        self._require(
            'name',
            'category',
            'source',
            'builder.type',
        )
        if not self.name.islower():
            raise ValueError('package name must be all lower case, got ' + self.name)
        if self.category not in ['standard', 'optional', 'experimental', 'huge']:
            raise ValueError('package category unknown, got ' + self.category)

    def _normalize_format(self, config, formatting):
        def recurse_into(config):
            if isinstance(config, basestring):
                return config.format(**formatting)
            elif isinstance(config, list):
                return map(recurse_into, config)
            elif isinstance(config, dict):
                for key in config.keys():
                    config[key] = recurse_into(config[key])
                return config
            else:
                return config
        recurse_into(config)
       
 
