"""
Interface to the Sage cleaner

Triva Note: For the name "sage-cleaner", think of the
"The Cleaner" from Pulp Fiction:
http://www.frankjankowski.de/quiz/illus/keitel.jpg
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

from sage.misc.misc import SAGE_TMP

def cleaner(pid, cmd=''):
    """
    Write a line to the ``spawned_processes`` file with the given
    ``pid`` and ``cmd``.
    """
    if cmd != '':
        cmd = cmd.strip().split()[0]
    # This is safe, since only this process writes to this file.
    F = os.path.join(SAGE_TMP, 'spawned_processes')
    try:
        with open(F, 'a') as o:
            o.write('%s %s\n'%(pid, cmd))
    except IOError:
        pass
    else:
        start_cleaner()


def start_cleaner():
    """
    Start ``sage-cleaner`` in a new process group.
    """
    if not os.fork():
        os.setpgid(os.getpid(), os.getpid())

        # Redirect stdin, stdout, stderr to /dev/null
        with open(os.devnull, "r+") as f:
            os.dup2(f.fileno(), 0)
            os.dup2(f.fileno(), 1)
            os.dup2(f.fileno(), 2)

        # Close all other file descriptors
        try:
            maxopenfiles = os.sysconf("SC_OPEN_MAX")
            if maxopenfiles <= 0:
                raise ValueError
        except ValueError:
            maxopenfiles = 1024
        os.closerange(3, maxopenfiles)

        os.execlp("sage-cleaner", "sage-cleaner")
