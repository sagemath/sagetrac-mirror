r"""
Prefix Coding

This module implements functionalities relating to prefix coding.

AUTHORS:

- Jan Wabbersen (2015-04-27): Initial version based on concepts and
  ideas of the Huffman module by Nathann Cohen.
"""

#*****************************************************************************
#       Copyright (C) 2015 Jan Wabbersen <jan.wabbersen@googlemail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from abc import ABCMeta, abstractmethod
from itertools import product

from six import string_types, iteritems

from sage.structure.sage_object import SageObject
from sage.graphs.digraph import DiGraph
from .misc import frequency_table

# Key for code length of the longest code in a decoding table.
MAX_KEY_LEN = "key_len"

class PrefixCoding(SageObject):
    r"""
    This class implements basic functionalities for prefix codes.
    
    Note that this class must not be initiated. It is a base class for
    homomorphic coding methods, which produce prefix codes.
    """
    __metaclass__ = ABCMeta
    
    def __init__(self, source, verbose=False, q=2, char_per_symbol=1,
                 decoding_table_key_len=8):
        r"""
        Constructor for prefix coding methods.

        Note that the object initialization depends on _build_code()
        which is not implemented in this class. Thus, this class
        cannot be used on its own.

        INPUT:

        - ``source`` -- can be either

            - A string from which the encoding table should be created.

            - A dictionary that associates to each symbol of an
              alphabet a numeric value. If we consider the frequency
              of each alphabetic symbol, then ``source`` is considered
              as the frequency table of the alphabet with each numeric
              (non-negative integer) value being the number of
              occurrences of a symbol. The numeric values can also
              represent weights of the symbols. In that case, the
              numeric values are not necessarily integers, but can be
              real numbers.
        
        - ``verbose`` -- (default: False) if True, print intermediate
          data of building the code.
          
        - ``q`` -- (default: 2) the number of symbols in the code
          alphabet. Supported is 2 <= q <= 10.
          
        - ``char_per_symbol`` -- (default: 1) the number of characters
          that define one symbol.
          
        - ``decoding_table_key_len`` -- (default: 8) the length of each
          key in the resulting decoding table. Longer keys may result in
          faster decoding. Decoding tables with shorter keys consume
          less memory, but decoding may be slower.
        """
        if type(self) == PrefixCoding:
            raise Exception("PrefixCoding must not be initiated. See the "
                            "documentation for more information.")
        
        # Encoding table
        self._character_to_code = {}
        self._decoding_table = {}
        # Saves how many characters define a symbol.
        self._char_per_symbol = char_per_symbol
        
        if isinstance(source, string_types):
            if (self._char_per_symbol != 1 and
                len(source) % self._char_per_symbol != 0):
                raise ValueError("The passed string does not match with the "
                                 "passed value for char_per_symbol.")
            if q == 2:
                self._build_code(
                    frequency_table(source, self._char_per_symbol), verbose)
            else:
                self._build_code(
                    frequency_table(source, self._char_per_symbol), verbose, q)
        elif isinstance(source, dict):
            keys = source.keys()
            self._char_per_symbol = len(keys[0])
            if any(len(k) != self._char_per_symbol for k in keys):
                raise ValueError("All keys in the passed dictionary must have "
                                 "the same length.")
            if q == 2:
                self._build_code(source, verbose)
            else:
                self._build_code(source, verbose, q)
        else:
            raise ValueError("Input must be either a string or a dictionary.")
        self._decoding_table = self._generate_decoding_table(
            self._character_to_code, q, decoding_table_key_len)
        
    @abstractmethod
    def _build_code(self, dic, verbose=False, q=2):
        r"""
        Construct a code corresponding to an alphabet with the
        given weight table.

        INPUT:

        - ``dic`` -- a dictionary that associates to each symbol of an
          alphabet a numeric value. If we consider the frequency of
          each alphabetic symbol, then ``dic`` is considered as the
          frequency table of the alphabet with each numeric
          (non-negative integer) value being the number of occurrences
          of a symbol. The numeric values can also represent weights of
          the symbols. In that case, the numeric values are not
          necessarily integers, but can be real numbers. In general,
          we refer to ``dic`` as a weight table.
          
        - ``verbose`` -- (default: False) if True, print the results
          of each step of building the code.
          
        - ``q`` -- (default: 2) number of symbols in the code alphabet.
          Supported is 2 <= q <= 10.
        """
        raise NotImplementedError("This method needs to be implemented by a "
                                  "specific coding method.")
    
    def encode(self, string):
        r"""
        Encode the given string based on the current encoding table.
        
        INPUT:

        - ``string`` -- a string of symbols over an alphabet.
          
        OUTPUT:

        - An encoding of ``string``.
        
        EXAMPLES:
        
        See one of the specific prefix coding methods (e.g., Shannon).
        """
        if self._character_to_code:
            if self._char_per_symbol == 1:
                return "".join([self._character_to_code[x] for x in string])
            else:
                if len(string) % self._char_per_symbol != 0:
                    raise ValueError("The passed string contains unknown "
                                     "symbols.")
                return "".join(
                    [self._character_to_code[x] for x in
                    [string[i:i + self._char_per_symbol]
                    for i in range(0, len(string), self._char_per_symbol)]])
    
    def decode(self, string):
        r"""
        Decode the given string based on the current decoding table.
        
        INPUT:

        - ``string`` -- a string of encoded symbols.
          
        OUTPUT:

        - The decoding of ``string``.
        
        EXAMPLES:
        
        See one of the specific prefix coding methods (e.g., Shannon).
            
        ALGORITHM:
        
        Uses an adapted version of [GA]_.
        
        REFERENCES:
        
        .. [GA] J.-L. Gailly and M. Adler. Decompression algorithm
           (inflate). http://www.gzip.org/algorithm.txt (accessed
           on February 10, 2015)
        """
        output = []
        table = self._decoding_table
        window = table[MAX_KEY_LEN]
        cur_window = window
        start = 0
        v = None
        while start < len(string):
            try:
                v = table[string[start:start + window]]
                while isinstance(v, dict):
                    start += cur_window
                    cur_window = v[MAX_KEY_LEN]
                    v = v[string[start:start + cur_window]]
                output.append(v[0])
                start += v[1]
                cur_window = window
            except KeyError:
                if not isinstance(v, dict):
                    cur_window = window
                    v = table
                last_bits = string[start:start + cur_window]
                v = v.get(last_bits + (cur_window - len(last_bits)) * '0',
                          None)
                if v is not None:
                    output.append(v[0])
                    break;
                else:
                    raise ValueError("Input must be a string encoded by this "
                                     "instance.")
        return "".join(output)
    
    def _generate_decoding_table(self, char_code, q=2, max_len=8):
        r"""
        Generate a decoding table with the given encoding table.
        
        INPUT:

        - ``char_code`` -- encoding table.
        
        - ``q`` -- (default: 2) number of symbols in the code alphabet.
          Supported is 2 <= q <= 10.
          
        - ``max_len`` -- (default: 8) the max. length of each
          key in the resulting decoding table. Longer keys may result in
          faster decoding. Decoding tables with shorter keys consume
          less memory, but decoding may be slower.
          
        OUTPUT:

        - A decoding table for ``char_code``.
        
        EXAMPLES::
        
            sage: from sage.coding.source_coding.shannon import Shannon
            sage: str = "decoding table"
            sage: s = Shannon(str)
            sage: encoding_table = s.encoding_table()
            sage: s._generate_decoding_table(encoding_table, max_len=3)
            {'000': ('e', 3),
             '001': ('d', 3),
             '010': {'0': ('a', 1), '1': (' ', 1), 'key_len': 1},
             '011': {'0': ('c', 1), '1': ('b', 1), 'key_len': 1},
             '100': {'1': ('g', 1), 'key_len': 1},
             '101': {'0': ('i', 1), '1': ('l', 1), 'key_len': 1},
             '110': {'0': ('o', 1), '1': ('n', 1), 'key_len': 1},
             '111': {'0': ('t', 1), 'key_len': 1},
             'key_len': 3}
             
        Note that this is a private method. Thus, it is not designed
        for end users.
            
        ALGORITHM:
        
        Uses an adapted version of [GA]_.
        
        REFERENCES:
        
        .. [GA] J.-L. Gailly and M. Adler. Decompression algorithm
           (inflate). http://www.gzip.org/algorithm.txt (accessed
           on February 10, 2015)
        """
        if isinstance(char_code, dict):
            char_code = char_code.items()
        char_code.sort(key=lambda x: len(x[1]))
        split = next(
            (i for i, (a, c) in enumerate(char_code) if len(c) > max_len),
            None)
        # Create a table for this level and table(s) for the next level.
        if split is not None:
            decoding_table = self._fill_up_decoding_table(
                char_code[:split], q, max_len)
            next_level = char_code[split:]
            next_level.sort(key=lambda x: x[1])
            old_s, old_c = next_level[0]
            next_table = [(old_s, old_c[max_len:])]
            for s, c in next_level[1:]:
                if c[:max_len] == old_c[:max_len]:
                    next_table.append((s, c[max_len:]))
                else:
                    decoding_table[old_c[:max_len]] = self \
                        ._generate_decoding_table(next_table, q, max_len)
                    next_table = [(s, c[max_len:])]
                old_c = c
            # Create last table of this level.
            decoding_table[old_c[:max_len]] = self._generate_decoding_table(
                next_table, q, max_len)
            return decoding_table
        # Only create a table for this level.
        else:
            return self._fill_up_decoding_table(char_code, q=q)
    
    def _fill_up_decoding_table(self, table, q=2, key_len=None):
        r"""
        Generate one table of the decoding table.
        
        Switch key and value and fill up the new keys with all possible
        permutations of the symbols 0,...,q-1 so that all keys have
        the same (minimum) length. Note that this is just a helper
        method for _generate_decoding_table().
        
        INPUT:

        - ``table`` -- table to be filled up.
        
        - ``q`` -- (default: 2) number of symbols for the permutations.
          Supported is 2 <= q <= 10.
          
        - ``key_len`` -- (default: None) the key length of the resulting
          table.
          
        OUTPUT:

        - A filled up table for ``table``.
        
        TESTS::
        
            sage: from sage.coding.source_coding.shannon import Shannon
            sage: s = Shannon("example")
            sage: encoding_table = s.encoding_table()
            sage: s._fill_up_decoding_table(encoding_table.items())
            {'000': ('e', 2),
             '001': ('e', 2),
             '010': ('a', 3),
             '011': ('m', 3),
             '100': ('l', 3),
             '101': ('p', 3),
             '110': ('x', 3),
             'key_len': 3}
        """
        digits = [str(i) for i in range(0, q)]
        d = {}
        if key_len is None:
            # Length of longest code.
            key_len = len(max(table, key=lambda x: len(x[1]))[1])
        perm = []
        old_suffix_len = -1  # for first iteration
        for s, c in table:
            suffix_len = key_len - len(c)
            if not suffix_len == old_suffix_len:
                perm = ["".join(x) for x in product(digits, repeat=suffix_len)]
            # Create all entries for the current code.
            for p in perm:
                d[c + p] = (s, len(c))
            old_suffix_len = suffix_len
        # Adding length of the keys for decoding purposes.
        d[MAX_KEY_LEN] = key_len
        return d
    
    def encoding_table(self):
        r"""
        Return the current encoding table.

        INPUT:

        - None.

        OUTPUT:

        - A dictionary associating an alphabetic symbol to an encoding.

        EXAMPLES:
        
        See one of the specific prefix coding methods (e.g., Shannon).
        """
        return self._character_to_code.copy()
    
    def tree(self):
        r"""
        Return the tree corresponding to the current encoding.

        INPUT:

        - None.

        OUTPUT:

        - The tree representing the code.

        EXAMPLES:
        
        See one of the specific prefix coding methods (e.g., Shannon).
        """
        g = DiGraph()
        edges = set()
        for s, c in iteritems(self._character_to_code):
            duplicate = False
            parent = c[:-1]
            if parent:
                edges.add((parent, s + ": " + c))
            else:
                edges.add(('root', s + ": " + c))
                continue
            for _ in range(len(parent) - 1):
                if (parent[:-1], parent) in edges:
                    duplicate = True
                    break
                edges.add((parent[:-1], parent))
                parent = parent[:-1]
            if not duplicate:
                edges.add(('root', parent))
        g.add_edges(list(edges))
        return g
