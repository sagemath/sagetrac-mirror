"""
The Sage Input Hook

This is a hook into the IPython input prompt and will be called
periodically (every 100ms) while Python is sitting idle. We use it to
reload attached files if they have changed. See :class:`SageInputhook`
for more details.

IPython has analogous code to set an input hook, but we are not using
their implementation. For once, it unsets signal handlers which will
disable Ctrl-C.
"""

###########################################################################
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

include 'sage/ext/stdsage.pxi'
include 'sage/ext/interrupt.pxi'

cdef extern from 'pythonrun.h':
    int (*PyOS_InputHook)() nogil except *

import sage.libs.readline as readline
from sage.misc.attached_files import reload_attached_files_if_modified


cdef int c_sage_inputhook() nogil except *:
    """
    This is the C function that is installed as PyOS_InputHook
    """
    with gil:
        try:
            sig_check()
            return sage_inputhook()
        except KeyboardInterrupt:
            # The user pressed Ctrl-C while at the prompt; We match the normal
            # Python behavior for consistency
            print '\nKeyboardInterrupt'
            readline.initialize()
            readline.forced_update_display()
        return 0

def install():
    """
    Install the Sage input hook

    EXAMPLES::

        sage: from sage.misc.inputhook import install
        sage: install()
    """
    global PyOS_InputHook
    PyOS_InputHook = c_sage_inputhook

def uninstall():
    """
    Uninstall the Sage input hook

    EXAMPLES::

        sage: from sage.misc.inputhook import uninstall
        sage: uninstall()
    """
    global PyOS_InputHook
    PyOS_InputHook = NULL



sage_inputhook = None

class SageInputhook(object):

    def __init__(self):
        """
        The Sage Inputhook.

        The inputhook lets you run code while Sage is otherwise
        idle. There is no parallelism / threads involved, so you have
        to return quickly in case the user wants to type something.

        You can add your own callbacks using :meth:`add`. Valid
        callbacks are any callable without arguments. It will be
        called shortly after Sage returns to the input prompt. If the
        callable returns ``True``, it will be kept and called again
        after about 100ms. Once the callable returns ``False``, it
        will be removed from the callback list and not called
        again. There is no guarantee about the order of multiple
        callbacks.

        Callbacks are weakly referenced, so they are also removed if
        the callable is garbage collected.

        EXAMPLES::
         
            sage: from sage.misc.inputhook import sage_inputhook
            sage: sage_inputhook
            <sage.misc.inputhook.SageInputhook object at 0x...>
        """
        assert sage_inputhook is None   # singleton object
        self._callbacks = set()
        
    def add(self, callback):
        """
        Add an idle callback.

        INPUT:

        - ``callback`` -- a callable that can be weak referenced. Can
          be a function or a bound method.

        EXAMPLES::
        
            sage: def callme(): print 'idle'; return False
            sage: from sage.misc.inputhook import sage_inputhook
            sage: sage_inputhook.add(callme) 
            sage: sage_inputhook()    
            idle
            0
            sage: sage_inputhook()    
            0
        """
        self._callbacks.add(callback)

    def discard(self, callback):
        """
        Discard a callback.

        INPUT:

        - ``callback`` -- a callable to be discarded.

        EXAMPLES::
        
            sage: def callme(): print 'idle'; return False
            sage: from sage.misc.inputhook import sage_inputhook
            sage: sage_inputhook.add(callme)
            sage: sage_inputhook.discard(callme)
            sage: sage_inputhook()
            0

        Note that it is allowed to discard callables that are not
        currently in the list of callbacks. The following is an no-op
        and does not raise an exception::

            sage: sage_inputhook.discard(integral)
        """
        self._callbacks.discard(callback)
            
    def __call__(self):
        """
        The input hook.
   
        This object will be called every 100ms when IPython is idle at
        the command prompt.
    
        EXAMPLES::
    
            sage: from sage.misc.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: tmp = tmp_filename(ext='.py')
            sage: f = open(tmp, 'w'); f.write('a = 2\n'); f.close()
            sage: shell.run_cell('%attach ' + tmp)
            sage: shell.run_cell('a')
            2
            sage: sleep(1)  # filesystem timestamp granularity
            sage: f = open(tmp, 'w'); f.write('a = 3\n'); f.close()
    
        Note that the doctests are never really at the command prompt, so
        we call the input hook manually::
    
            sage: shell.run_cell('from sage.misc.inputhook import sage_inputhook')
            sage: shell.run_cell('sage_inputhook()')
            ### reloading attached file tmp_....py modified at ... ###
            0
    
            sage: shell.run_cell('a')
            3
            sage: shell.run_cell('detach({0})'.format(repr(tmp)))
            sage: shell.run_cell('attached_files()')
            []
        """
        reload_attached_files_if_modified()
        finished = []
        for callback in self._callbacks:
            rc = callback()
            if not rc:
                finished.append(callback)
        self._callbacks.difference_update(finished)
        return 0
        
sage_inputhook = SageInputhook()



