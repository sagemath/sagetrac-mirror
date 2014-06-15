"""
The main routines
"""

from sage_pkg.config import config
from sage_pkg.package_list import load_config, loader
from sage_pkg.directories import make_required_directories


def build_queue(pkg=None, stop_at='install'):
    if pkg is None:
        packages = loader.get_all()
    else:
        packages = loader.get(pkg)
    return packages.build_queue(stop_at=stop_at)


class Application(object):

    def __init__(self):
        """
        The Python constructor
        """
        make_required_directories()
        
    def list_all(self, category=None):
        """
        Print a table of packages 

        EXAMPLES::

            >>> app.list_all()
        """
        def print_line(name, category, version, status, age):
            print('{0:<27} {1:>7}  {2:>16}  {3:>15}  {4:>11}'.format(
                name, category, version, status, age))
        print_line(
            'Package name',
            'Category',
            'Version',
            'Status',
            'Age'
        )
        print('-' * 79)
        for pkg in loader.get_all():
            print_line(pkg.name, pkg.category, pkg.version, pkg.status, pkg.age)

    def info(self, pkg):
        """
        Print information about a package
        """
        pkg = loader.get(pkg)
        print(pkg._config)
        print('Status: ' + pkg.status)
        print('Installed when: ' + pkg.age)
        
    def upgrade_v1(self):
        """
        Generate package.yaml from old-style directories
        """
        from sage_pkg.upgrade_v1 import upgrade_all
        upgrade_all()
    

    def _build_until_step(self, pkg, step_name):
        queue = build_queue(pkg, stop_at=step_name)
        queue.run_serial()

    def download(self, pkg=None):
        """
        Stop after downloading the sources
        """
        self._build_until_step(pkg, 'download')

    def unpack(self, pkg=None):
        """
        Stop after finishing the unpack step.
        """
        self._build_until_step(pkg, 'unpack')

    def prepare(self, pkg=None):
        """
        Stop after finishing the prepare step.
        """
        self._build_until_step(pkg, 'prepare')

    def configure(self, pkg=None):
        """
        Stop after finishing the configure step.
        """
        self._build_until_step(pkg, 'configure')
        pass

    def compile(self, pkg=None):
        """
        Stop after finishing the compile step.
        """
        self._build_until_step(pkg, 'compile')
        pass

    def check(self, pkg=None):
        """
        Run self-tests but don't install.
        """
        self._build_until_step(pkg, 'check')

    def install(self, pkg=None):
        """
        Install package.
        """
        self._build_until_step(pkg, 'install')
        
