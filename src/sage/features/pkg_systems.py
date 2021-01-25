r"""
Package systems as features
"""

from . import Feature

_cache_package_systems = None

def package_systems():
    """
    Return a list of ``PackageSystem`` objects representing the available package systems.

    The list is ordered by decreasing preference.

    EXAMPLES::

        sage: from sage.features import package_systems
        sage: package_systems()    # random
        [Feature('homebrew'), Feature('sage_spkg'), Feature('pip')]
    """
    # The current implementation never returns more than one system.
    from subprocess import run, CalledProcessError, PIPE
    global _cache_package_systems
    if _cache_package_systems is None:
        _cache_package_systems = []
        # Try to use scripts from SAGE_ROOT (or an installation of sage_bootstrap)
        # to obtain system package advice.
        try:
            proc = run('sage-guess-package-system', shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True, check=True)
            system_name = proc.stdout.strip()
            if system_name != 'unknown':
                _cache_package_systems = [PackageSystem(system_name)]
        except CalledProcessError:
            pass
        more_package_systems = [SagePackageSystem(), PipPackageSystem()]
        _cache_package_systems += [ps for ps in more_package_systems if ps.is_present()]

    return _cache_package_systems

class PackageSystem(Feature):
    r"""
    A feature describing a system package manager.

    EXAMPLES::

        sage: from sage.features import PackageSystem
        sage: PackageSystem('conda')
        Feature('conda')
    """
    def _is_present(self):
        r"""
        Test whether ``self`` appears in the list of available package systems.

        EXAMPLES::

            sage: from sage.features import PackageSystem
            sage: debian = PackageSystem('debian')
            sage: debian.is_present()  # indirect doctest, random
            True
        """
        return self in package_systems()

    def spkg_installation_hint(self, spkgs, *, prompt="  !", feature=None):
        r"""
        Return a string that explains how to install ``feature``.

        EXAMPLES::

            sage: from sage.features import PackageSystem
            sage: homebrew = PackageSystem('homebrew')
            sage: homebrew.spkg_installation_hint('openblas')  # optional - SAGE_ROOT
            'To install openblas using the homebrew package manager, you can try to run:\n!brew install openblas'
        """
        if isinstance(spkgs, (tuple, list)):
            spkgs = ' '.join(spkgs)
        if feature is None:
            feature = spkgs
        return self._spkg_installation_hint(spkgs, prompt, feature)

    def _spkg_installation_hint(self, spkgs, prompt, feature):
        r"""
        Return a string that explains how to install ``feature``.

        Override this method in derived classes.

        EXAMPLES::

            sage: from sage.features import PackageSystem
            sage: fedora = PackageSystem('fedora')
            sage: fedora.spkg_installation_hint('openblas')  # optional - SAGE_ROOT
            'To install openblas using the fedora package manager, you can try to run:\n!sudo yum install openblas-devel'
        """
        from subprocess import run, CalledProcessError, PIPE
        lines = []
        system = self.name
        try:
            proc = run(f'sage-get-system-packages {system} {spkgs}',
                       shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True, check=True)
            system_packages = proc.stdout.strip()
            print_sys = f'sage-print-system-package-command {system} --verbose --sudo --prompt="{prompt}"'
            command = f'{print_sys} update && {print_sys} install {system_packages}'
            proc = run(command, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True, check=True)
            command = proc.stdout.strip()
            if command:
                lines.append(f'To install {feature} using the {system} package manager, you can try to run:')
                lines.append(command)
                return '\n'.join(lines)
        except CalledProcessError:
            pass
        return f'No equivalent system packages for {system} are known to Sage.'

class SagePackageSystem(PackageSystem):
    r"""
    The feature describing the Sage package manager.

    EXAMPLES::

        sage: from sage.features import SagePackageSystem
        sage: SagePackageSystem()
        Feature('sage_spkg')
    """
    @staticmethod
    def __classcall__(cls):
        r"""
        Normalize initargs.

        TESTS::

            sage: from sage.features import SagePackageSystem
            sage: SagePackageSystem() is SagePackageSystem()  # indirect doctest
            True
        """
        return PackageSystem.__classcall__(cls, "sage_spkg")

    def _is_present(self):
        r"""
        Test whether ``sage-spkg`` is available.

        EXAMPLES::

            sage: from sage.features import SagePackageSystem
            sage: bool(SagePackageSystem().is_present())  # indirect doctest, optional - sage_spkg
            True
        """
        from subprocess import run, DEVNULL, CalledProcessError
        try:
            # "sage -p" is a fast way of checking whether sage-spkg is available.
            run('sage -p', shell=True, stdout=DEVNULL, stderr=DEVNULL, check=True)
            return True
        except CalledProcessError:
            return False

    def _spkg_installation_hint(self, spkgs, prompt, feature):
        r"""
        Return a string that explains how to install ``feature``.

        EXAMPLES::

            sage: from sage.features import SagePackageSystem
            sage: print(SagePackageSystem().spkg_installation_hint(['foo', 'bar'], prompt="### ", feature='foobarability'))  # indirect doctest
            To install foobarability using the Sage package manager, you can try to run:
            ### sage -i foo bar
        """
        lines = []
        lines.append(f'To install {feature} using the Sage package manager, you can try to run:')
        lines.append(f'{prompt}sage -i {spkgs}')
        return '\n'.join(lines)

class PipPackageSystem(PackageSystem):
    r"""
    The feature describing the Pip package manager.

    EXAMPLES::

        sage: from sage.features import PipPackageSystem
        sage: PipPackageSystem()
        Feature('pip')
    """
    @staticmethod
    def __classcall__(cls):
        r"""
        Normalize initargs.

        TESTS::

            sage: from sage.features import PipPackageSystem
            sage: PipPackageSystem() is PipPackageSystem()  # indirect doctest
            True
        """
        return PackageSystem.__classcall__(cls, "pip")

    def _is_present(self):
        r"""
        Test whether ``pip`` is available.

        EXAMPLES::

            sage: from sage.features import PipPackageSystem
            sage: bool(PipPackageSystem().is_present())    # indirect doctest
            True
        """
        from subprocess import run, DEVNULL, CalledProcessError
        try:
            run('sage -pip --version', shell=True, stdout=DEVNULL, stderr=DEVNULL, check=True)
            return True
        except CalledProcessError:
            return False
