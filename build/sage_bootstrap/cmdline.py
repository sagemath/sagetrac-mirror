# -*- coding: utf-8 -*-
"""
Commandline handling

Note that argparse is not part of Python 2.6, so we cannot rely on it here.
"""


#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import sys
import re
from textwrap import dedent
import logging
log = logging.getLogger()

from sage_bootstrap.env import SAGE_DISTFILES
from sage_bootstrap.package import Package
from sage_bootstrap.download import Download
from sage_bootstrap.mirror_list import MirrorList
from sage_bootstrap.tarball import Tarball


class CmdlineSubcommands(object):

    KEYWORD_RE = re.compile('--(?P<option>[A-Za-z][A-Za-z0-9]*)=(?P<value>.*)')
    
    def __init__(self, argv=None):
        if argv is None:
            argv = sys.argv
        if len(argv) == 1:
            self.print_help()
            sys.exit(0)
        self.subcommand = argv[1]
        self._init_argparse(argv[2:])

    def _init_argparse(self, cmdline_args):
        args = []
        kwds = {}
        for cmd in cmdline_args:
            match = self.KEYWORD_RE.match(cmd)
            if match is None:
                args.append(cmd)
            else:
                kwds[match.group('option')] = match.group('value')
        self.extra_args = args
        self.extra_kwds = kwds
            
    def print_help(self):
        print(dedent(self.__doc__).lstrip())
        print('Usage:')
        for method in sorted(dir(self)):
            if not method.startswith('run_'):
                continue
            doc = dedent(getattr(self, method).__doc__).lstrip().splitlines()
            print('')
            print('* ' + doc[0])
            for line in doc[1:]:
                print('  ' + line)
        
    def run(self):
        try:
            method = getattr(self, 'run_{0}'.format(self.subcommand))
        except AttributeError:
            log.error('unknown subcommand: {0}'.format(self.subcommand))
            sys.exit(1)
        try:
            method(*self.extra_args, **self.extra_kwds)
        except TypeError as err:
            log.error(str(err))
            log.error('invalid arguments to the {0} subcommand: {1} {2}'
                      .format(self.subcommand, self.extra_args, self.extra_kwds))
            sys.exit(1)

            
class SagePkgApplication(CmdlineSubcommands):
    """
    sage-package
    --------
    
    The package script is used to manage third-party tarballs.
    """
    
    def run_config(self):
        """
        config: Print the configuration

        $ sage --package config
        Configuration:
          * log = info
          * interactive = True
        """
        from sage_bootstrap.config import Configuration
        print(Configuration())

    def run_list(self):
        """
        list: Print a list of all available packages

        $ sage --package list | sort
        4ti2
        arb
        atlas
        autotools
        [...]
        zn_poly
        """
        for pkg in Package.all():
            print(pkg.name)

    def run_name(self, tarball_filename):
        """
        name: Find the package name given a tarball filename
    
        $ sage --package name pari-2.8-1564-gdeac36e.tar.gz
        pari
        """
        tarball = Tarball(os.path.basename(tarball_filename))
        print(tarball.package.name)

    def run_tarball(self, package_name):
        """
        tarball: Find the tarball filename given a package name
    
        $ sage --package tarball pari
        pari-2.8-1564-gdeac36e.tar.gz
        """
        package = Package(package_name)
        print(package.tarball.filename)

    def run_apropos(self, incorrect_name):
        """
        apropos: Find up to 5 package names that are close to the given name

        $ sage --package apropos python
        Did you mean: cython, ipython, python2, python3, patch?
        """
        from sage_bootstrap.levenshtein import Levenshtein, DistanceExceeded
        levenshtein = Levenshtein(5)
        names = []
        for pkg in Package.all():
            try:
                names.append([levenshtein(pkg.name, incorrect_name), pkg.name])
            except DistanceExceeded:
                pass
        if names:
            names = sorted(names)[:5]
            print('Did you mean: {0}?'.format(', '.join(name[1] for name in names)))
        else:
            print('There is no package similar to {0}'.format(incorrect_name))
            print('You can find further packages at http://files.sagemath.org/spkg/')

    def run_update(self, package_name, new_version, url=None):
        """
        update: Update a package. This modifies the Sage sources. 
    
        $ sage --package update pari 2015 --url=http://localhost/pari/tarball.tgz
        """
        from sage_bootstrap.updater import PackageUpdater
        update = PackageUpdater(package_name, new_version)
        if url is not None:
            update.download_upstream(url)
        update.fix_checksum()


class SageDownloadFileApplication(object):    
    """
    USAGE:
    
        sage-download-file --print-fastest-mirror
    
    Print out the fastest mirror. All further arguments are ignored in
    that case.
    
        sage-download-file [--quiet] url-or-tarball [destination]
    
    The single mandatory argument can be a http:// url or a tarball
    filename. In the latter case, the tarball is downloaded from the
    mirror network and its checksum is verified.
    
    If the destination is not specified:
    * a url will be downloaded and the content written to stdout
    * a tarball will be saved under {SAGE_DISTFILES}
    """
    
    def run(self):
        progress = True
        url = None
        print_fastest_mirror = None
        destination = None
        for arg in sys.argv[1:]:
            if arg.startswith('--print-fastest-mirror'):
                url = ""
                print_fastest_mirror = True
                continue
            if arg.startswith('--quiet'):
                progress = False
                continue
            if url is None:
                url = arg
                continue
            if destination is None:
                destination = arg
                continue
            raise ValueError('too many arguments')
        if url is None:
            print(dedent(self.__doc__.format(SAGE_DISTFILES=SAGE_DISTFILES)))
            sys.exit(2)

        try:
            if url.startswith('http://') or url.startswith('https://') or url.startswith('ftp://'):
                Download(url, destination, progress=progress, ignore_errors=True).run()
            elif print_fastest_mirror:
                url = "fastest mirror"  # For error message
                print(MirrorList().fastest)
            else:
                # url is a tarball name
                tarball = Tarball(url)
                tarball.download()
                if destination is not None:
                    tarball.save_as(destination)
        except BaseException:
            try:
                stars = '*' * 72 + '\n'
                sys.stderr.write(stars)
                try:
                    import traceback
                    traceback.print_exc(file=sys.stderr)
                    sys.stderr.write(stars)
                except:
                    pass
                sys.stderr.write("Error downloading %s\n"%(url,))
                sys.stderr.write(stars)
            finally:
                sys.exit(1)
