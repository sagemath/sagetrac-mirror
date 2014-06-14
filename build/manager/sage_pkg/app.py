"""
The main routines
"""

from sage_pkg.config import config

from sage_pkg.package_list import load_config, loader


class Application(object):

    def __init__(self):
        """
        The Python constructor
        """
        
    def table(self):
        """
        Print a table of packages 
        """
        for pkg in loader.get_all():
            print(pkg)

    def info(self, pkg):
        """
        Print information about a package
        """
        pkg_config = load_config(pkg)
        print(pkg_config)

    def list_all(self, category=None):
        """
        Print a table of packages 
        """
        for pkg in loader.get_all():
            print(pkg)
        
    def upgrade_v1(self):
        """
        Generate package.yaml from old-style directories
        """
        from sage_pkg.upgrade_v1 import upgrade_all
        upgrade_all()
    
