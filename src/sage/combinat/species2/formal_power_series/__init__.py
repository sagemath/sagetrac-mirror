from sage.categories.formal_power_series import FormalPowerSeries
from sage.misc.abstract_method import abstract_method
from sage.rings.infinity import Infinity
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

class ValuationFPS:

    def __init__(self):
        self._listeners = set()
        self._valuation_init_()

    def _valuation_registration_(self, listeners):
        for f in listeners:
            if isinstance(f, ValuationFPS):
                f._valuation_add_listener_(self)

    def _valuation_add_listener_(self, listener):
        self._listeners.add(listener)

    def _valuation_init_(self):
        self._valuation = Infinity

    def _valuation_update_(self):
        tmp = self._valuation
        self._valuation = self._valuation_compute_()
        if tmp != self._valuation:
            map(lambda listener: listener._valuation_update_(), self._listeners)

    @abstract_method(optional=False)
    def _valuation_compute_(self):
        """
        The way to compute the valuation

        :return: an Integer

        (For example, `val(f+g) = min(val(f), val(g))`.)
        """

    def valuation(self):
        return self._valuation