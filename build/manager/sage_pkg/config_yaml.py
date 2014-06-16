"""
Configuration management with YAML
"""

import os
import collections
import yaml



def recursive_merge(dst, src):
    """
    Recursively merge dictionaries

    INPUT:

    - ``dst`` -- a dictionary containing dictionaries and so on.

    - ``src`` -- a dictionary containing dictionaries and so on.

    OUTPUT:

    The ``src`` dictionary merged into the ``src``
    dictionary. (Sub-)keys of the same name are overwritten with the
    value in ``src``. The ``src`` dictionary is overwritten.

    EXAMPLES::

        >>> from sage_pkg.config_yaml import recursive_merge
        >>> src = {'a': 1, 'b': 1}
        >>> dst = {'b': 2, 'c': 2}
        >>> recursive_merge(dst, src) == {'a': 1, 'b': 1, 'c': 2}
        True
    """
    if isinstance(src, dict) and isinstance(dst, dict):
        for key, value in dst.iteritems():
            if key not in src:
                src[key] = value
            else:
                src[key] = recursive_merge(value, src[key])
    return src


def recursive_list_attrs(prefix, obj):
    for key in sorted(dir(obj)):
        if key.startswith('_'):
            continue
        value = getattr(obj, key)
        if isinstance(value, ConfigYAML):
            for attr in recursive_list_attrs('{0}.{1}'.format(prefix, key), value):
                yield attr
        else:
            yield '{0}.{1} = {2}'.format(prefix, key, value)
            

class FrozenDictProxy(collections.Mapping):

    def __init__(self, dictionary):
        self._dictionary = dictionary

    def __getitem__(self, key):
        return self._dictionary.__getitem__(key)

    def __iter__(self):
        return self._dictionary.__iter__()

    def __len__(self):
        return self._dictionary.__len__()

    def __repr__(self):
        return 'frozen' + repr(self._dictionary)


class ConfigYAML(object):
    
    def __init__(self, *args, **kwds):
        self._merge(*args)
        self._init_attributes()
        self._validate()
        self._name = kwds.pop('name', 'config')
        if kwds:
            raise ValueError('unknown keyword arguments: {0}'.format(kwds))

    def _normalize(self, config):
        """
        Hook for normalizing the values in the configuration dictionary

        This is run before the object is initialized, so you cannot
        access any attributes of ``self``.
        
        INPUT:
        
        - ``config`` -- deserialized yaml. Nested
          dictionary/list/values.

        The return value is discarded, you are supposed to change
        values in place.
        """
        pass

    def _validate(self):
        """
        Hook for validating the configuration
        
        Unlike :meth:`_normalize`, this method is run *after* the
        object initialization. You can still raise a ``ValueError`` if
        the cofiguration is invalid, but you cannot change it any
        more.
        """
        pass

    def _require(self, *args):
        """
        Require the presence of configuration keys

        INPUT:

        - ``*args`` -- each argument must be a string, or a list of
          strings. Hierarchical sub-keys can also be joined by dots.

        EXAMPLES::

            >>> config._require('path.configuration')
            >>> config.path.configuration
            []
            >>> config._require('path.foo.bar')
            Traceback (most recent call last):
            ...
            ValueError: required key is not defined: path.foo.bar
        """
        for arg in args:
            if isinstance(arg, basestring):
                arg = arg.split('.')
            value = self._c
            for key in arg:
                try:
                    value = value[key]
                except TypeError:
                    raise ValueError('parent entity is list, need dict: ' + '.'.join(arg))
                except KeyError:
                    raise ValueError('required key is not defined: ' + '.'.join(arg))

    def _merge(self, *args):
        """
        Merge configurations

        INPUT:

        - ``args`` -- each argument can either be a dictionary/mapping
          collection or a string. The latter must be the name of a
          yaml configuration file.

        - ``normalize``

        OUTPUT:

        A dictionary obtained by merging the dictionaries from the
        left to the right.
        """
        c = dict()
        for arg in args:
            if isinstance(arg, basestring):
                with open(arg, 'r') as f:
                    additional_c = yaml.load_all(f).next()
            else:
                additional_c = arg
            c = recursive_merge(c, additional_c)
            #print c
        self._normalize(c)
        self._c = FrozenDictProxy(c)

    def _init_attributes(self):
        for key, value in self._c.iteritems():
            if isinstance(value, dict):
                value = ConfigYAML(value)
            setattr(self, key, value)
                
    def __repr__(self):
        """
        Return a string representation
        """
        result = ['Configuration:']
        prefix = '- ' + self._name
        for line in recursive_list_attrs(prefix, self):
            result.append(line)
        return '\n'.join(result)
        
    def __getitem__(self, key):
        """
        Get nested values

        EXAMPLES::
        
            >>> config.mirrors.extra
            ['http://download.example.com/packages']
            >>> config.mirrors['extra']
            ['http://download.example.com/packages']
        """
        return self._c[key]

    def __call__(self, *args, **kwds):
        """
        Get nested values

        EXAMPLES::

            >>> config('path', 'install')   # doctest: +ELLIPSIS
            '/.../test_data/local'
            >>> config('path', 'bar')
            Traceback (most recent call last):
            ...
            KeyError: 'bar'
            >>> config('path', 'bar', default='baz')
            'baz'
        """
        try:
            default = kwds.pop('default')
            have_default = True
        except KeyError:
            have_default = False
        value = self._c
        try:
            for key in args:
                if not isinstance(value, (dict, FrozenDictProxy)):
                    raise KeyError('parent is not a dict, cannot get key = ' + key)
                value = value[key]
        except KeyError:
            if have_default:
                return default
            else:
                raise
        return value
