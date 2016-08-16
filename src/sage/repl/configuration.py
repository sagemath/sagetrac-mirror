r"""
Sage's IPython Configuration
"""

#*****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import

import sys
import copy
from traitlets.config.loader import Config

from sage.repl.prompts import SagePrompts


# Name of the Sage IPython extension
SAGE_EXTENSION = 'sage'


class SageIpythonConfiguration(object):

    def _doctest_mode(self):
        """
        Whether we are in doctest mode

        This returns ``True`` during doctests.

        EXAMPLES::
        
            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config._doctest_mode()
            True
        """
        from sage.doctest import DOCTEST_MODE
        return DOCTEST_MODE
    
    def _allow_ansi(self):
        """
        Whether to allow ANSI escape sequences

        This returns ``False`` during doctests to avoid ANSI escape 
        sequences.

        EXAMPLES::
        
            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config._allow_ansi()
            False
        """
        return (not self._doctest_mode()) and sys.stdout.isatty()

    def colors(self):
        """
        Return the IPython color palette

        This returns ``'NoColor'`` during doctests to avoid ANSI escape 
        sequences.

        EXAMPLES::
        
            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config.simple_prompt()
            True
        """
        return 'LightBG' if self._allow_ansi() else 'NoColor'

    def simple_prompt(self):
        """
        Return whether to use the simple prompt

        This returns ``True`` during doctests to avoid ANSI escape sequences.

        EXAMPLES::
        
            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config.simple_prompt()
            True
        """
        return not self._allow_ansi()

    def term_title(self):
        """
        Return whether to set the terminal title

        This returns false during doctests to avoid ANSI escape sequences.

        EXAMPLES::
        
            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config.term_title()
            False
        """
        return self._allow_ansi()

    def default(self):
        """
        Return a new default configuration object

        EXAMPLES::
        
            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config.default()
            {'InteractiveShell': {'colors': ...
        """
        from sage.repl.interpreter import SageTerminalInteractiveShell
        cfg = Config(
            TerminalIPythonApp=Config(
                display_banner=False,
                verbose_crash=True,
                test_shell=False,
                shell_class=SageTerminalInteractiveShell,
            ),
            InteractiveShell=Config(
                prompts_class=SagePrompts,
                ast_node_interactivity='all',
                colors=self.colors(),
                simple_prompt=self.simple_prompt(),
                term_title=self.term_title(),
                confirm_exit=False,
                separate_in=''
            ),
            InteractiveShellApp=Config(extensions=[SAGE_EXTENSION]),
        )
        if self._doctest_mode():
            # Using the file-backed history causes problems in parallel tests
            cfg.HistoryManager = Config(hist_file=':memory:')
        return cfg

    def copy(self):
        """
        Return a copy of the current configuration

        EXAMPLES::
        
            sage: from sage.repl.configuration import sage_ipython_config
            sage: sage_ipython_config.copy()
            {'InteractiveShell': {'colors': ...
        """
        try:
            return copy.deepcopy(get_ipython().config)
        except NameError:
            return self.default()


sage_ipython_config = SageIpythonConfiguration()
