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


class ConfigYAML(object):
    
    def __init__(self, *args):
        self._merge(*args)
        self._init_attributes()
        #print self

    def _normalize(self, config):
        """
        Hook for normalizing the values in the configuration dictionary
        """
        pass

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
        for line in recursive_list_attrs('- config', self):
            result.append(line)
        return '\n'.join(result)
        



class SagePkgConfig(ConfigYAML):
    """
    The Configuration for the Sage Package Manager itself

    EXAMPLES::

        >>> config    # doctest: +ELLIPSIS
        Configuration:
        - config.path.configuration = []
        - config.path.dot_git = /.../.git
        - config.path.packages = /.../build/manager/test_pkg
        - config.path.root = /...
        - config.path.sage_pkg = /.../build/manager
    """
    
    def _normalize(self, config):
        path = config.get('path', dict())
        self._normalize_paths(path)

    def _normalize_paths(self, path, root=None):
        if root is None:
            root = path['root']
            root = os.path.abspath(root)
            path['root'] = root
    
        def normalize(value):
            value = os.path.expanduser(value)
            if not os.path.isabs(value):
                value = os.path.join(root, value)
            return os.path.abspath(value)

        for key in path.keys():
            value = path[key]
            if isinstance(value, list):
                path[key] = map(normalize, value)
            elif isinstance(value, dict):
                self._normalize_paths(value, root)
            else:
                path[key] = normalize(value)
        
