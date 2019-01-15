from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

lazy_import("sage.games.sudoku_backtrack", 'backtrack_all', deprecation=27066)

from .sudoku_backtrack import backtrack_all
from .sudoku import Sudoku, sudoku
from .hexad import Minimog

del absolute_import
