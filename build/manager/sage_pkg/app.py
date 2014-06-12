"""
The main routines
"""

from sage_pkg.config import config

from sage_pkg.tiny_git import GitRepository
from sage_pkg.package_list import load_packages


class Application(object):

    def __init__(self):
        """
        The Python constructor
        """
        self.git = GitRepository(config.path.dot_git)
        
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
        
        

    
