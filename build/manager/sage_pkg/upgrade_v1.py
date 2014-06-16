"""
Tools for upgrading Sage to the first version of the package manager
"""




PACKAGE_YAML_TEMPLATE = """
name:
    {name}

category:
    {category}

source:
    version:  
        {version}
    tarball:
        name:     {tarball}
        sha1:     {sha1}

builder:
    type:            SpkgInstallScript
    install_script:  spkg-install
""".strip('\n')

CHECK_ADDON = """
    check_script:    spkg-check
""".strip('\n')

DEPS_PREFIX = """
depends:
    hard:
""".strip('\n')

DEPS_TEMPLATE = """
        - {dep}
""".strip('\n')


DEPENDENCIES = {
    'sqlite': ['readline', ],
    'bzip2': ['pkgconf'],
    'patch': ['bzip2'],
    'atlas': ['python'],
    'boehm_gc': ['pkgconf'],
    'boost_cropped': [],
    'cliquer': [],
    'ncurses': ['pkgconf'],
    'readline': ['ncurses'],
    'iconv': [],
    'docutils': ['python'],
    'elliptic_curves': ['python', 'sqlite'],
    'conway': [],
    'graphs': [],
    'glpk': ['mpir', 'zlib'],
    'python': ['zlib', 'bzip2', 'pkgconf', 'readline', 'sqlite', 'libpng'],
    'mpir': ['iconv'],
    'gsl': ['atlas'],
    'gf2x': [],
    'ntl': ['mpir', 'gf2x'],
    'libfplll': ['mpir', 'mpfr'],
    'pari': ['readline', 'mpir', 'pari_galdata', 'pari_seadata_small'],
    'pari_galdata': [],
    'pari_seadata_small': [],
    'polybori': ['python', 'ipython', 'scons', 'boost_cropped', 'libm4ri', 'gd'],
    'polytopes_db': [],
    'ppl': ['mpir', 'glpk'],
    'mpc': ['mpir', 'mpfr'],
    'mpfr': ['mpir'],
    'mpfi': ['mpir', 'mpfr'],
    'givaro': ['mpir'],
    'git': ['zlib', 'python'],
    'fflas_ffpack': ['mpir', 'givaro', 'gsl', 'atlas'],
    'linbox': [
        'mpir', 'ntl', 'givaro',
        'mpfr', 'libfplll', 'iml',
        'libm4ri', 'libm4rie', 'fflas_ffpack'],
    'iml': ['mpir', 'gsl', 'atlas'],
    'genus2reduction': ['pari'],
    'palp': [],
    'lcalc': ['pari', 'mpfr'],
    'lrcalc': [],
    'pynac': ['python'],
    'sympow': [],
    'symmetrica': [],
    'gap': ['ncurses', 'readline', 'mpir'],
    'libgap': ['gap'],
    'ipython': ['python'],
    'pexpect': ['python'],
    'gd': ['libpng', 'freetype', 'iconv'],
    'gdmodule': ['python', 'gd', 'iconv'],
    'scons': ['python'],
    'rubiks': [],
    'sqlite': ['readline'],
    'sagetex': ['python',
                'maxima', 'scipy',
                'matplotlib', 'pillow', 'tachyon'],
    'setuptools': ['python'],
    'singular': ['mpir', 'ntl', 'readline', 'mpfr'],
    'pycrypto': ['python'],
    'networkx': ['python'],
    'mpmath': ['python'],
    'zlib': [],
    'jmol': ['sagenb'],
    'freetype': ['libpng'],
    'libpng': ['zlib'],
    'six': ['python'],
    'dateutil': ['python', 'six', 'setuptools'],
    'pyparsing': ['python'],
    'tornado': ['python'],
    'matplotlib': [
        'python', 'numpy',
        'freetype', 'libpng',
        'gdmodule', 'dateutil',
        'pkgconf', 'pyparsing',
        'setuptools'],
    'cddlib': ['mpir'],
    'gfan': ['mpir', 'cddlib'],
    'tachyon': ['libpng'],
    'ecm': ['mpir'],
    'ratpoints': ['mpir'],
    'ecl': ['mpir', 'readline', 'boehm_gc'],
    'maxima': ['ecl'],
    'r': ['atlas', 'iconv', 'readline', 'pkgconf'],
    'rpy': ['python', 'r'],
    'sympy': ['python'],
    'cython': ['python'],
    'flintqs': ['mpir'],
    'flint': ['mpir', 'mpfr', 'ntl'],
    'eclib': ['pari', 'ntl', 'flint'],
    'libm4ri': ['libpng', 'pkgconf'],
    'libm4rie': ['libm4ri', 'givaro', 'ntl'],
    'znpoly': ['mpir', 'python'],
    'sagenb': [
        'python', 'setuptools', 'pexpect',
        'jinja2', 'sphinx', 'docutils'],
    'sqlalchemy': ['python', 'setuptools'],
    'sphinx': [
        'python', 'setuptools', 'docutils',
        'jinja2', 'pygments'],
    'jinja2': ['python', 'setuptools', 'docutils'],
    'pygments': ['python', 'setuptools'],
    'gcc': ['mpir', 'mpfr', 'mpc', 'zlib'],
    'pillow': ['python', 'setuptools'],
    'pkgconf': [],
    'pkgconfig': ['python', 'setuptools'],
    'numpy': ['python', 'atlas', 'pkgconf'],
    'scipy': ['atlas', 'numpy'],
    'cvxopt': [
        'numpy',
        'atlas', 'cephes',
        'gsl', 'glpk', 'matplotlib'],
    'cephes': [],
}


import os
import re
from sage_pkg.config import config


def upgrade_all():
    """
    Upgrade all packages
    """
    for name in os.listdir(config.path.packages):
        fullname = os.path.join(config.path.packages, name)
        if not os.path.isdir(fullname):
            continue
        if not os.path.exists('package-version.txt'):
            continue
        upgrade(name, fullname)


def old_version(pkg_dir):
    with open(os.path.join(pkg_dir, 'package-version.txt'), 'r') as f:
        version_patchlevel = f.read().strip()
    return re.sub(r'\.p[0-9]*', '', version_patchlevel)


def old_sha1(pkg_dir):
    with open(os.path.join(pkg_dir, 'checksums.ini'), 'r') as f:
        for line in f.readlines():
            if line.startswith('tarball='):
                tarball = line[8:].strip()
            if line.startswith('sha1='):
                sha1 = line[5:].strip()
    return tarball, sha1


def make_package_yaml(name, pkg_dir):
    """
    Package.yaml from old-style directory
    
    EXAMPLES::

        >>> pkg_dir = os.path.join(config.path.test_data, 'oldpkg')
        >>> from sage_pkg.upgrade_v1 import make_package_yaml
        >>> print(make_package_yaml('sqlite', pkg_dir))
        name:
            sqlite
        <BLANKLINE>
        category:
            standard
        <BLANKLINE>
        source:
            version:  
                4.7.4
            tarball:
                name:     gap-{source.version}.tar.bz2
                sha1:     e61699cedf56b3b1668dec7387437eeb08fba3e7
        <BLANKLINE>
        builder:
            type:            SpkgInstallScript
            install_script:  spkg-install
            check_script:    spkg-check
        <BLANKLINE>
        depends:
            hard:
                - readline
        <BLANKLINE>
    """
    if name in DEPENDENCIES:
        category = 'standard'
        dependencies = DEPENDENCIES[name]
    else:
        category = 'optional'
    tarball, sha1 = old_sha1(pkg_dir)
    tarball = tarball.replace('VERSION', '{source.version}')
    result = PACKAGE_YAML_TEMPLATE.format(
        name=name,
        version=old_version(pkg_dir),
        sha1=sha1,
        tarball=tarball,
        category=category
    )
    result += '\n'
    if os.path.exists(os.path.join(pkg_dir, 'spkg-check')):
        result += CHECK_ADDON
    result += '\n'
    if category == 'standard':
        result += '\n'
        result += DEPS_PREFIX + '\n'
        for dep in dependencies:
            result += DEPS_TEMPLATE.format(dep=dep) + '\n'
    return result


def upgrade(name, pkg_dir):
    """
    Upgrade a single package
    """
    with open(os.path.join(pkg_dir, 'package.yaml'), 'w') as yaml:
        package_yaml = make_package_yaml(name, pkg_dir)
        yaml.write(package_yaml)
