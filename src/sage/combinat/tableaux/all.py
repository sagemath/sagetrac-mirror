r"""
Tableaux features that are imported by default in the interpreter namespace
"""
from abstract_tableaux  import AbstractTableaux
from abstract_tableau   import AbstractTableau
from bad_shape_tableaux import BadShapeTableaux
from bad_shape_tableau  import BadShapeTableau
from skew_tableaux      import SkewTableaux     as SSkewTableaux # change back
from skew_tableau       import SkewTableau      as SSkewTableau

# TODO: organize all tableaux under a "tableaux" namespace like graphs are and
#   add a zillion deprecation warnings
