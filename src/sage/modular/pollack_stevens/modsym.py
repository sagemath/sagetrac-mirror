from sage.structure.element import ModuleElement

class ModularSymbolElement(ModuleElement):
    def __init__(self, map, parent):
        ModuleElement.__init__(self, parent)
        self._map = map

    def _add_(self, right):
        return self.__class__(self._map + right._map, self.parent())

    def _lmul_(self, right):
        return self.__class__(self._map * right, self.parent())

    def _sub_(self, right):
        return self.__class__(self._map - right._map, self.parent())

    
