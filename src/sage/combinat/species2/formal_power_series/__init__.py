from sage.categories.formal_power_series import FormalPowerSeries
from sage.structure.category_object import CategoryObject
from sage.structure.dynamic_class import dynamic_class
from sage.structure.unique_representation import UniqueRepresentation


class FPS(UniqueRepresentation, CategoryObject):

    def __init__(self, category=FormalPowerSeries()):
        CategoryObject.__init__(self, category=category)

        # ##### CUSTOM ####### used to inherit the parent methods
        base = self.__class__.__base__
        self.__class__ = dynamic_class("%s_with_category" % base.__name__,
                                       (self.__class__, self.category().parent_class, ),
                                       doccls=base)
        #####################
