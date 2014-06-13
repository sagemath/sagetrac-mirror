"""
The main routines
"""

from sage_pkg.config import config

from sage_pkg.package_list import load_packages, load_config


class Application(object):

    def __init__(self):
        """
        The Python constructor
        """
        
    def table(self):
        """
        Print a table of packages 
        """
        for pkg in load_packages():
            print(pkg)

    def info(self, pkg):
        """
        Print information about a package
        """
        pkg_config = load_config(pkg)
        print(pkg_config)

    def list(self, category=None):
        """
        Print a table of packages 
        """
        for pkg in load_packages():
            print(pkg)
        
        

    
