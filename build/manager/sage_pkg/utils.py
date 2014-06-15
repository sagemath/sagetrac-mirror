"""
Utilities
"""


class cached_property(object):
    """
    Memoized ``@property``
    """

    def __init__(self, func, readonly=True):
        self.func = func
        self.__name__ = func.__name__
        self.__doc__ = func.__doc__
        self.is_readonly = readonly

    def __get__(self, obj, type=None):
        if obj is None:
            return self
        try:
            return obj.__dict__[self.__name__]
        except KeyError:
            value = self.func(obj)
            obj.__dict__[self.__name__] = value
            return value

    def __set__(self, obj, value):
        if self.is_readonly:
            raise TypeError('read only property')
        obj.__dict__[self.__name__] = value
