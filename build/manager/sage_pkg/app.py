"""
The main routines
"""

import yaml

from .tiny_git import GitRepository
from .package import load_package, all_packages


class Application(object):

    def __init__(self, configfile):
        self._load_config(configfile)
        self.git = GitRepository(self.dot_git)
        
    def _load_config(self, configfile):
        with open(configfile, 'r') as f:
            self.config = yaml.load_all(f).next()
        dirname = os.path.abspath(os.dirname(configfile))
        def path(config_key):
            p = self.config[config_key]
            if os.path.isabs(p):
                p = os.path.abspath(p)
            else:
                p = os.path.abspath(os.path.join(dirname, p))
            if not os.path.exists(p):
                raise ValueError('unknown path: {0}'.format(p))
            return p
        self.dot_git = path('dot_git')
        self.package_dir = path('package_dir')

    def table(self):
        """
        Print a table of packages 
        """
        for pkg in all_packages(self.package_dir, self.git):
            print(pkg)

    def info(self, pkg=None):
        """
        Print a table of packages 
        """
        if pkg is None:
            return self.table()
        pkg = load_package(pkg, self.git)
        
        

    
