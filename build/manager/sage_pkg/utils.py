"""
Utilities
"""

import datetime


def random_sha1():
    """
    Return a random SHA1

    OUTPUT:

    String containing 40 hex digits.

    EXAMPLES::
    
        >>> from sage_pkg.utils import random_sha1
        >>> sha1 = random_sha1()
        >>> sha1  # doctest: +SKIP
        '3835bb5604b33160a94f47ee8d4262b9471c0017'
        >>> is_valid_sha1(sha1)
        True
    """
    import datetime
    now = str(datetime.datetime.utcnow())
    import sha
    return sha.sha(now).hexdigest()


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



def pretty_age(time_utc):
    """
    Return string for age
    
    OUTPUT:

    A pretty string like 'an hour ago', 'Yesterday', '3 months ago',
    'just now', etc.
    """
    now = datetime.datetime.utcnow()
    diff = now - time_utc
    second_diff = diff.seconds
    day_diff = diff.days
    if day_diff < 0:
        return 'from future?'
    if day_diff == 0:
        if second_diff < 10:
            return "just now"
        if second_diff < 60:
            return str(second_diff) + " seconds ago"
        if second_diff < 120:
            return  "a minute ago"
        if second_diff < 3600:
            return str( second_diff / 60 ) + " minutes ago"
        if second_diff < 7200:
            return "an hour ago"
        if second_diff < 86400:
            return str( second_diff / 3600 ) + " hours ago"
    if day_diff == 1:
        return "yesterday"
    if day_diff < 7:
        return str(day_diff) + " days ago"
    if day_diff < 31:
        return str(day_diff/7) + " weeks ago"
    if day_diff < 365:
        return str(day_diff/30) + " months ago"
    return str(day_diff/365) + " years ago"
