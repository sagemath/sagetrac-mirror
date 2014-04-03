from __future__ import absolute_import
# Makes sage.combinat.root_system? equivalent to sage.combinat.root_system.root_system?
from .root_system import __doc__

# currently needed to activate the backward compatibility register_unpickle_override
from . import type_A
from . import type_B
from . import type_C
from . import type_D
from . import type_E
from . import type_F
from . import type_G

from . import all
