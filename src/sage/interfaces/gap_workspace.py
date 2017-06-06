r"""
Support for (lib)GAP workspace files
"""

#*****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import time
from sage.env import DOT_SAGE, SAGE_LOCAL


def gap_workspace_file(system="gap", name="workspace", dir=None):
    r"""
    Return the filename for the GAP workspace.

    INPUT:

    - ``system`` -- the name of the system, either ``"gap"`` or
      ``"libgap"``

    - ``name`` -- the kind of workspace, usually ``"workspace"`` but
      libGAP also uses other files

    - ``dir`` -- the directory where the workspaces should be stored.
      By default, this is ``DOT_SAGE/gap``

    If needed, this creates ``dir``. It also removes any obsolete
    workspaces from ``dir``.

    EXAMPLES::

        sage: from sage.interfaces.gap_workspace import gap_workspace_file
        sage: gap_workspace_file("foo", "bar", "/somewhere")
        '/somewhere/foo-bar-...'

    TESTS::

        sage: from sage.env import DOT_SAGE
        sage: D = gap_workspace_file()
        sage: D.startswith(os.path.join(DOT_SAGE, "gap", "gap-workspace-"))
        True
    """
    if dir is None:
        dir = os.path.join(DOT_SAGE, 'gap')

    h = hex(abs(hash(SAGE_LOCAL)))
    return os.path.join(dir, '%s-%s-%s' % (system, name, h))


def prepare_workspace_dir(dir=None):
    r"""
    Create and clean up the directory for GAP workspaces.

    INPUT:

    - ``dir`` -- the directory where the workspaces should be stored.
      By default, this is ``DOT_SAGE/gap``

    OUTPUT: the actual workspace directory

    EXAMPLES::

        sage: from sage.interfaces.gap_workspace import prepare_workspace_dir
        sage: prepare_workspace_dir()
        '.../gap'

    TESTS::

        sage: prepare_workspace_dir(os.path.join(tmp_dir(), "new"))
        '.../new'
    """
    if dir is None:
        dir = os.path.join(DOT_SAGE, 'gap')

    # Make sure that the workspace directory exists
    try:
        os.makedirs(dir)
    except OSError:
        if not os.path.isdir(dir):
            raise
    else:
        # Directory was created, add a README file
        with open(os.path.join(dir, 'README.txt'), 'w') as f:
            f.write("It is OK to delete all these cache files.  They will be recreated as needed.\n")

    # Delete all gap workspaces that haven't been used in the last
    # week, to avoid needless cruft.  I had an install on sage.math
    # with 90 of these, since I run a lot of different versions of
    # Sage, and it totalled 1.3GB of wasted space!  See trac #4936.
    # We only do this after creating a new workspace, since this cruft
    # issue is only a problem if workspaces get created every so
    # often.  We don't want to have to do this on every startup.
    now = time.time()
    for F in os.listdir(dir):
        if "workspace" in F:
            W = os.path.join(dir, F)
            try:
                age = now - os.path.getatime(W)
                if age >= 604800:    # 1 week in seconds
                    os.unlink(W)
            except OSError:
                # It's not a problem if W doesn't exist, everything
                # else is an error.
                if os.path.exists(W):
                    raise

    return dir
