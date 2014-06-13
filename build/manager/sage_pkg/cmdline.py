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



def debug_shell(app, parser):
    from IPython.frontend.terminal.ipapp import TerminalIPythonApp
    ip = TerminalIPythonApp.instance()
    ip.initialize(argv=[])
    from develop.config import config as dev_config
    ip.shell.user_global_ns['dev_config'] = dev_config
    from .config import config
    ip.shell.user_global_ns['app'] = app
    ip.shell.user_global_ns['log'] = logger
    ip.shell.user_global_ns['config'] = config
    from sage_pkg.tiny_git import GitRepository
    ip.shell.user_global_ns['git'] = GitRepository(config.path.dot_git)
    from sage_pkg.package_list import loader
    ip.shell.user_global_ns['loader'] = loader
    def ipy_import(module_name, identifier):
        import importlib
        module = importlib.import_module(module_name)
        ip.shell.user_global_ns[identifier] = getattr(module, identifier) 
    # ipy_import('sage_pkg.git_interface', 'GitInterface')
    ip.start()



description = \
"""
The Sage Package Manager
"""



def launch(DEFAULT_CONFIG):
    from argparse import ArgumentParser
    parser = ArgumentParser(description=description)
    parser.add_argument('--debug', dest='debug', action='store_true',
                        default=False, 
                        help='debug')
    parser.add_argument('--config', dest='config', type=str, default=None, 
                        help='builder configuration file')
    parser.add_argument('--log', dest='log', default=None,
                        help='one of [DEBUG, INFO, ERROR, WARNING, CRITICAL]')
    subparsers = parser.add_subparsers(dest='subcommand')

    # sage-pkg info <package>
    parser_info = subparsers.add_parser('info', help='Print information')
    parser_info.add_argument('package', type=str, help='Package name')

    # sage-pkg help
    parser_help = subparsers.add_parser('help', help='Get help')

    args = parser.parse_args()
    if args.log is not None:
        import logging
        level = getattr(logging, args.log)
        logger.setLevel(level=level)

    set_configuration(args.config if args.config else DEFAULT_CONFIG)

    from .app import Application
    app = Application()

    if args.debug:
        debug_shell(app, parser)
    elif args.subcommand == 'info':
        app.info(args.package)
    elif args.subcommand == 'help':
        parser.print_help()
    else:
        print('Unknown subcommand "{0}"'.format(args.subcommand))
        parser.print_help()
