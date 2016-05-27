import os
from sage.env import SAGE_SRC

def autogen_all():
    """
    Regenerate the automatically generated files of the Sage library.
    """
    from . import pari
    pari.rebuild()

    from . import interpreters
    interpreters.rebuild(os.path.join(SAGE_SRC, "sage", "ext", "interpreters"))
