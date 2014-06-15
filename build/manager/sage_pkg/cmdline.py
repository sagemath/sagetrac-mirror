## -*- encoding: utf-8 -*-
"""
Handle Command Line Options
"""

##############################################################################
#  Sage Package Manager
#  Copyright (C) 2014  Volker Braun <vbraun.name@gmail.com>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################


import sys
import os
import warnings

from .logger import logger
from .config import set_configuration



def debug_shell(app, package, parser):
    from IPython.frontend.terminal.ipapp import TerminalIPythonApp
    ip = TerminalIPythonApp.instance()
    ip.initialize(argv=[])
    from develop.config import config as dev_config
    ip.shell.user_global_ns['dev_config'] = dev_config
    from .config import config
    ip.shell.user_global_ns['app'] = app
    ip.shell.user_global_ns['log'] = logger
    ip.shell.user_global_ns['config'] = config
    from sage_pkg.git.repo import GitRepository
    ip.shell.user_global_ns['git'] = GitRepository(config.path.dot_git)
    from sage_pkg.package_list import loader
    ip.shell.user_global_ns['loader'] = loader
    ip.shell.user_global_ns[package] = loader.get(package)
    def ipy_import(module_name, identifier):
        import importlib
        module = importlib.import_module(module_name)
        ip.shell.user_global_ns[identifier] = getattr(module, identifier) 
    # ipy_import('sage_pkg.git_interface', 'GitInterface')
    print('Type "{0}" to access the {0} package'.format(package))
    ip.start()



description = \
"""
The Sage Package Manager
"""



def launch(DEFAULT_CONFIG):
    from argparse import ArgumentParser
    parser = ArgumentParser(description=description)
    parser.add_argument('--config', dest='config', type=str, default=None, 
                        help='builder configuration file')
    parser.add_argument('--log', dest='log', default=None,
                        help='one of [DEBUG, INFO, ERROR, WARNING, CRITICAL]')
    subparsers = parser.add_subparsers(dest='subcommand')

    # sage-pkg info <package>
    parser_info = subparsers.add_parser('info', help='Print information about package')
    parser_info.add_argument('package', type=str, help='Package name')

    # sage-pkg debug <package>
    parser_debug = subparsers.add_parser('debug', help='Debug package')
    parser_debug.add_argument('package', type=str, help='Package name')

    # sage-pkg info <package>
    parser_list = subparsers.add_parser('list', help='List all packages')

    # sage-pkg pkg-upgrade-v1
    parser_upgrade_v1 = subparsers.add_parser('pkg-upgrade-v1', help='Upgrade package descriptions')

    # sage-pkg help
    parser_help = subparsers.add_parser('help', help='Get help')

    # The different build steps
    parser_download = subparsers.add_parser('download', help='Build up to the "download" step')
    parser_download.add_argument('package', type=str, help='Package name')

    parser_unpack = subparsers.add_parser('unpack', help='Build up to the "unpack" step')
    parser_unpack.add_argument('package', type=str, help='Package name')

    parser_prepare = subparsers.add_parser('prepare', help='Build up to the "prepare" step')
    parser_prepare.add_argument('package', type=str, help='Package name')

    parser_configure = subparsers.add_parser('configure', help='Build up to the "configure" step')
    parser_configure.add_argument('package', type=str, help='Package name')

    parser_compile = subparsers.add_parser('compile', help='Build up to the "compile" step')
    parser_compile.add_argument('package', type=str, help='Package name')

    parser_check = subparsers.add_parser('check', help='Build up to the "check" step')
    parser_check.add_argument('package', type=str, help='Package name')

    parser_install = subparsers.add_parser('install', help='Build and install')
    parser_install.add_argument('package', type=str, help='Package name')

    # Build the whole thing
    parser_build = subparsers.add_parser('build', help='Build everything')

    args = parser.parse_args()
    if args.log is not None:
        import logging
        level = getattr(logging, args.log.upper())
        logger.setLevel(level=level)

    set_configuration(args.config if args.config else DEFAULT_CONFIG)

    from .app import Application
    app = Application()

    if args.subcommand == 'debug':
        debug_shell(app, args.package, parser)
    elif args.subcommand == 'info':
        app.info(args.package)
    elif args.subcommand == 'list':
        app.list_all()
    elif args.subcommand == 'pkg-upgrade-v1':
        app.upgrade_v1()
    elif args.subcommand == 'download':
        app.download(args.package)
    elif args.subcommand == 'unpack':
        app.unpack(args.package)
    elif args.subcommand == 'prepare':
        app.prepare(args.package)
    elif args.subcommand == 'configure':
        app.configure(args.package)
    elif args.subcommand == 'compile':
        app.compile(args.package)
    elif args.subcommand == 'check':
        app.check(args.package)
    elif args.subcommand == 'install':
        app.install(args.package)
    elif args.subcommand == 'build':
        app.build()
    elif args.subcommand == 'help':
        parser.print_help()
    else:
        print('Unknown subcommand "{0}"'.format(args.subcommand))
        parser.print_help()
