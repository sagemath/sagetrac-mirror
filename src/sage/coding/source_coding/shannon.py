r"""
Shannon Coding

This module implements functionalities relating to Shannon encoding and
decoding.

AUTHORS:

- Jan Wabbersen (2015-04-27): initial version based on concepts and
  ideas of the Huffman module by Nathann Cohen
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
from __future__ import print_function, absolute_import

from math import ceil, log
from operator import itemgetter

from .prefix_coding import PrefixCoding 
from .misc import SimpleTable, first_binary_dec_places

class Shannon(PrefixCoding):
    r"""
    This class implements basic functionalities of Shannon coding.
    
    It can build a Shannon code from a given string, or from the
    information of a dictionary associating to each key (the elements
    of the alphabet) a weight (most of the time, a probability value
    or a number of occurrences).
    
    INPUT:

        - ``source`` -- can be either

            - A string from which the Shannon encoding should be
              created.
    
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
          
        - ``char_per_symbol`` -- (default: 1) the number of characters
          that define one symbol.
          
        - ``decoding_table_key_len`` -- (default: 8) the length of each
          key in the resulting decoding table. Longer keys may result in
          faster decoding. Decoding tables with shorter keys consume
          less memory, but decoding may be slower.
    
    In order to construct a Shannon encoding for an alphabet, we use
    exactly one of the following methods:

    #. Let ``source`` be a string of symbols over an alphabet and feed
       ``source`` to the constructor of this class. Based on the input
       string and ``char_per_code``, a frequency table is constructed 
       that contains the frequency of each unique symbol (defined by
       ``char_per_symbol``) in ``source``. The alphabet in question is
       then all the unique symbols in ``source``. A significant
       implication of this is that any subsequent string that we want to
       encode must contain only symbols that can be found in ``source``.

    #. Let ``source`` be the frequency table of an alphabet. We can
       feed this table to the constructor of this class. The table
       ``source`` can be a table of frequencies or a table of weights.
       
    EXAMPLES::
    
        sage: from sage.coding.source_coding.shannon import Shannon
        sage: s1 = Shannon("Encode me!")
        sage: for symbol, code in s1.encoding_table().iteritems():
        ....:     print("'" + symbol + "' : " + code)
        ....:     
        '!' : 0011
        ' ' : 0100
        'c' : 0110
        'e' : 000
        'd' : 1000
        'm' : 1001
        'o' : 1011
        'n' : 1100
        'E' : 1110

    Now, we have a Shannon code for encoding strings consisting of
    symbols that are in the encoding table::
    
        sage: encoded = s1.encode("Encode me!"); encoded
        '11101100011010111000000010010010000011'
        
    And we can decode the encoded string in the following way::
    
        sage: s1.decode(encoded)
        'Encode me!'

    If we try to encode a string with ``s1`` that consists of symbols
    which are not part of the encoding table, we of course get an
    error::
    
        sage: s1.encode("This string contains other symbols!")
        Traceback (most recent call last):
        ...
        KeyError: 'T'

    If we want to encode multiple strings over an alphabet, it might be
    a good idea to generate a Shannon code for this alphabet and a
    given probability distribution for the symbols::
    
        sage: d = {'a': 0.5, 'b': 0.25, 'c': 0.2, 'd': 0.05}
        sage: s2 = Shannon(d)
        sage: for symbol, code in s2.encoding_table().iteritems():
        ....:     print("'" + symbol + "' : " + code)
        ....:     
        'a' : 0
        'c' : 110
        'b' : 10
        'd' : 11110

    We can also use frequencies or in general weights for the symbols::
    
        sage: d2 = {'a': 2.5, 'b': 1.25, 'c': 1, 'd': 0.25}
        sage: s3 = Shannon(d2)
        sage: for symbol, code in s3.encoding_table().iteritems():
        ....:     print("'" + symbol + "' : " + code)
        ....:     
        'a' : 0
        'c' : 110
        'b' : 10
        'd' : 11110
    
    If you are interested in Shannon coding, you can output intermediate
    data of the generation of the Shannon code::
    
        sage: s4 = Shannon("code", verbose=True)
        | symbol | p    | cum_sum | code |
        | c      | 0.25 | 0       | 00   |
        | e      | 0.25 | 0.25    | 01   |
        | d      | 0.25 | 0.5     | 10   |
        | o      | 0.25 | 0.75    | 11   |
    
    With ``char_per_symbol`` we are able to define symbols consisting
    of multiple characters. If we choose a dictionary for ``source``,
    we can just use keys/symbols consisting of multiple characters. But
    we have to make sure they all have the same length. If we choose a
    string for ``source``, we have to specify how many characters build
    one symbol::
    
        sage: str = "Split me into symbols consisting of 3 characters"
        sage: s5 = Shannon(str, char_per_symbol=3)
        sage: for symbol, code in s5.encoding_table().iteritems():
        ....:     print("'" + symbol + "' : " + code)
        ....:     
        'nsi' : 0000
        'me ' : 0011
        'int' : 0110
        '3 c' : 0111
        'it ' : 1001
        'ols' : 0100
        'of ' : 0001
        'ng ' : 0010
        ' co' : 1111
        'ers' : 0101
        'sti' : 1000
        'act' : 1010
        'ymb' : 1011
        'Spl' : 1100
        'o s' : 1101
        'har' : 1110
        
    Obviously, we have to make sure that the length of the string is a
    multiple of ``char_per_symbol`` or else we get an error::
    
        sage: s6 = Shannon("four", char_per_symbol=3)
        Traceback (most recent call last):
        ...
        ValueError: The passed string does not match with the passed value for char_per_symbol.
    
    When decoding long strings with a code consisting of long codewords,
    you can try to adapt the parameter ``decoding_table_key_len`` for
    a faster decoding. Or if you want to save memory when dealing with
    long codewords, you can set it to a smaller value::
        
        sage: s7 = Shannon("Let's call it a day!", decoding_table_key_len=3)
    """
    
    def __init__(self, source, verbose=False, char_per_symbol=1,
                 decoding_table_key_len=8):
        r"""
        Constructor for Shannon.

        See the docstring of this class for full documentation.

        EXAMPLES::

            sage: from sage.coding.source_coding.shannon import Shannon
            sage: str = "Give me an example!"
            sage: s = Shannon(str, decoding_table_key_len=6)

        TESTS:

        Feeding anything else than a string or a dictionary::

            sage: Shannon(Graph())
            Traceback (most recent call last):
            ...
            ValueError: Input must be either a string or a dictionary.
        """
        PrefixCoding.__init__(self, source, verbose, 2, char_per_symbol,
                              decoding_table_key_len)
        
    def _build_code(self, dic, verbose=False):
        r"""
        Construct a Shannon code corresponding to an alphabet with
        the given weight table.

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
          
        - ``verbose`` -- (default: False) if True, print intermediate
          data of building the code.
          
        EXAMPLES::
        
            sage: from sage.coding.source_coding.shannon import Shannon
            sage: d = {'c': 4, 'o': 3, 'd': 2, 'e':1}
            sage: s = Shannon(d)
            sage: s._build_code(d)
        """
        total = sum(dic.values())
        dic.update((s, w * 1.0 / total) for (s, w) in dic.items())
        
        # Sort characters by decreasing probability.
        sorted_list = sorted(dic.items(),
                             key=itemgetter(1),
                             reverse=True)
        
        # Calculate cumulative probabilities.
        cumsum = [0]
        for i, (s, p) in enumerate(sorted_list[:-1]):
            cumsum.append(cumsum[i] + p)

        # Calculate bit sequences for characters.
        bin_list = []
        for i, cum_p in enumerate(cumsum):
            bits = int(ceil(-log(sorted_list[i][1], 2)))
            bin_cumsum = first_binary_dec_places(cum_p, bits)
            bin_list.append(bin_cumsum)
            
        symbols = [s for (s, p) in sorted_list]
        self._character_to_code = dict(zip(symbols, bin_list))
                
        if verbose:
            probs = [str(p) for (s, p) in sorted_list]
            t = SimpleTable(['symbol', 'p', 'cum_sum', 'code'])
            for i in range(len(symbols)):
                t.add_row([symbols[i], probs[i], str(cumsum[i]), bin_list[i]])
            t.print_me()
            
    def encode(self, string):
        r"""
        Encode the given string based on the current encoding table.
        
        INPUT:

        - ``string`` -- a string of symbols over an alphabet.
          
        OUTPUT:

        - A Shannon encoding of ``string``.
        
        EXAMPLES:
        
        This example illustrates how to encode a string::
        
            sage: from sage.coding.source_coding.shannon import Shannon
            sage: str = "Encode me!"
            sage: s = Shannon(str)
            sage: encoded = s.encode(str); encoded
            '11101100011010111000000010010010000011'

        """
        return PrefixCoding.encode(self, string)
    
    def decode(self, string):
        r"""
        Decode the given string based on the current decoding table.
        
        INPUT:

        - ``string`` -- a string of encoded symbols.
          
        OUTPUT:

        - The Shannon decoding of ``string``.
        
        EXAMPLES:
        
        This example illustrates how to encode and then decode a
        string::
        
            sage: from sage.coding.source_coding.shannon import Shannon
            sage: str = "Encode me!"
            sage: s = Shannon(str)
            sage: encoded = s.encode(str); encoded
            '11101100011010111000000010010010000011'
            sage: s.decode(encoded)
            'Encode me!'
        
        TESTS:
            
        It is an error to try to decode an encoding from another
        instance. This may raise an exception or result in
        a wrong decoding::
        
            sage: str2 = "I'm another string."
            sage: s2 = Shannon(str2)
            sage: encoded_s2 = s2.encode(str2); encoded_s2
            '101001001010111000001101000111001010110101011110011000011100010100111101000011000011110'
            sage: s.decode(encoded_s2)
            Traceback (most recent call last):
            ...
            ValueError: Input must be a string encoded by this instance.
        """
        return PrefixCoding.decode(self, string)
    
    def encoding_table(self):
        r"""
        Return the current encoding table.

        INPUT:

        - None.

        OUTPUT:

        - A dictionary associating an alphabetic symbol to an encoding.

        EXAMPLES::
        
            sage: from sage.coding.source_coding.shannon import Shannon
            sage: str = "Show my encoding table!"
            sage: s = Shannon(str)
            sage: t = sorted(s.encoding_table().items())
            sage: for symbol, code in t:
            ....:     print(symbol, code)
            ....:     
              000
            ! 11010
            S 11001
            a 01100
            b 01111
            c 01101
            d 10000
            e 0010
            g 10010
            h 10100
            i 10011
            l 10111
            m 10110
            n 0100
            o 0011
            t 11011
            w 11101
            y 11110
        """
        return PrefixCoding.encoding_table(self)
    
    def tree(self):
        r"""
        Return the tree corresponding to the current encoding.

        INPUT:

        - None.

        OUTPUT:

        - The binary tree representing a Shannon code.

        EXAMPLES::
        
            sage: from sage.coding.source_coding.shannon import Shannon
            sage: str = "Binary tree"
            sage: s = Shannon(str)
            sage: t = s.tree(); t
            Digraph on 22 vertices
            sage: t.show(layout='tree', vertex_size=1800, figsize=[8,8])
        """
        return PrefixCoding.tree(self)
