"""
LibGAP Workspace Support

The single purpose of this module is to provide the location of the
libgap saved workspace and a time stamp to invalidate saved
workspaces.
"""

import os
import glob
from sage.interfaces.gap_workspace import gap_workspace_file


def timestamp():
    """
    Return a time stamp for libgap

    OUTPUT:

    Float. Unix timestamp of the most recently changed LibGAP file.

    EXAMPLES::

        sage: from sage.libs.gap.saved_workspace import timestamp
        sage: timestamp()   # random output
        1406642467.25684
        sage: type(timestamp())
        <... 'float'>
    """
    libgap_dir = os.path.dirname(__file__)
    files = glob.glob(os.path.join(libgap_dir, '*'))
    if len(files) == 0:
        print('Unable to find LibGAP files.')
        return float('inf')
    return max(map(os.path.getmtime, files))


def workspace(name='workspace'):
    """
    Return the filename of the gap workspace and whether it is up to date.

    INPUT:

    - ``name`` -- string. A name that will become part of the
      workspace filename.

    OUTPUT:

    Pair consisting of a string and a boolean. The string is the
    filename of the saved libgap workspace (or that it should have if
    it doesn't exist). The boolean is whether the workspace is
    up-to-date. You may use the workspace file only if the boolean is
    ``True``.

    EXAMPLES::

        sage: from sage.libs.gap.saved_workspace import workspace
        sage: ws, up_to_date = workspace()
        sage: ws
        '/.../gap/libgap-workspace-...'
        sage: isinstance(up_to_date, bool)
        True
    """
    workspace = gap_workspace_file("libgap", name)
    try:
        workspace_mtime = os.path.getmtime(workspace)
    except OSError:
        # workspace does not exist
        return (workspace, False)
    return (workspace, workspace_mtime >= timestamp())
