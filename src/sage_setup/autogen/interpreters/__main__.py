import os
from sage.env import SAGE_SRC

from . import rebuild


rebuild(os.path.join(SAGE_SRC, "sage", "ext", "interpreters"))
