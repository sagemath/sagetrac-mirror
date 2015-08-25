r"""
Tableaux features that are imported by default in the interpreter namespace
"""
from abstract_tableaux  import AbstractTableaux
from abstract_tableau   import AbstractTableau
from bad_shape_tableaux import BadShapeTableaux
from bad_shape_tableau  import BadShapeTableauFactory          as BadShapeTableau

from skew_tableaux      import SkewTableaux
from skew_tableaux      import StandardSkewTableauxFactory     as StandardSkewTableaux
from skew_tableaux      import SemistandardSkewTableauxFactory as SemistandardSkewTableaux
from skew_tableau       import SkewTableauFactory              as SkewTableau
from skew_tableau       import SemistandardSkewTableauFactory  as SemistandardSkewTableau
from skew_tableau       import StandardSkewTableauFactory      as StandardSkewTableau

# TODO: organize all tableaux under a "tableaux" namespace like graphs are and
#   add a zillion deprecation warnings
