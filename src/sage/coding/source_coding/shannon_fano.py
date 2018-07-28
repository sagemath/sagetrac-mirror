r"""
Shannon-Fano Coding

This module implements functionalities relating to Shannon-Fano encoding
and decoding.

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
from __future__ import print_function

from operator import itemgetter

from prefix_coding import PrefixCoding
from misc import SimpleTable

class ShannonFano(PrefixCoding):
    r"""
    This class implements basic functionalities of Shannon-Fano coding.
    
    It can build a Shannon-Fano code from a given string, or from the
    information of a dictionary associating to each key (the elements
    of the alphabet) a weight (most of the time, a probability value
    or a number of occurrences).
    
    INPUT:

        - ``source`` -- can be either

            - A string from which the Shannon-Fano encoding should be
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
    
    In order to construct a Shannon-Fano encoding for an alphabet, we
    use exactly one of the following methods:

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
    
        sage: from sage.coding.source_coding.shannon_fano import ShannonFano
        sage: sf1 = ShannonFano("Encode me!")
        sage: for symbol, code in sf1.encoding_table().iteritems():
        ....:     print("'" + symbol + "' : " + code)
        ....:     
        '!' : 010
        ' ' : 0110
        'c' : 0111
        'e' : 00
        'd' : 100
        'm' : 101
        'o' : 110
        'n' : 1110
        'E' : 1111

    Now, we have a Shannon-Fano code for encoding strings consisting of
    symbols that are in the encoding table::
    
        sage: encoded = sf1.encode("Encode me!"); encoded
        '11111110011111010000011010100010'
        
    And we can decode the encoded string in the following way::
    
        sage: sf1.decode(encoded)
        'Encode me!'

    If we try to encode a string with ``sf1`` that consists of symbols
    which are not part of the encoding table, we of course get an
    error::
    
        sage: sf1.encode("This string contains other symbols!")
        Traceback (most recent call last):
        ...
        KeyError: 'T'

    If we want to encode multiple strings over an alphabet, it might be
    a good idea to generate a Shannon-Fano code for this alphabet and a
    given probability distribution for the symbols::
    
        sage: d = {'a': 0.5, 'b': 0.25, 'c': 0.2, 'd': 0.05}
        sage: sf2 = ShannonFano(d)
        sage: for symbol, code in sf2.encoding_table().iteritems():
        ....:     print("'" + symbol + "' : " + code)
        ....:     
        'a' : 0
        'c' : 110
        'b' : 10
        'd' : 111

    We can also use frequencies or in general weights for the symbols::
    
        sage: d2 = {'a': 2.5, 'b': 1.25, 'c': 1, 'd': 0.25}
        sage: sf3 = ShannonFano(d2)
        sage: for symbol, code in sf3.encoding_table().iteritems():
        ....:     print("'" + symbol + "' : " + code)
        ....:     
        'a' : 0
        'c' : 110
        'b' : 10
        'd' : 111
    
    If you are interested in Shannon-Fano coding, you can output
    intermediate data of the generation of the Shannon-Fano code::
    
        sage: sf4 = ShannonFano("code", verbose=True)
        | level | partition 1 | partition 2 | code |
        | 1     | c, e        | d, o        |      |
        | 2     | c           | e           | 0    |
        |       | c           |             | 00   |
        |       | e           |             | 01   |
        | 2     | d           | o           | 1    |
        |       | d           |             | 10   |
        |       | o           |             | 11   |
        
    With ``char_per_symbol`` we are able to define symbols consisting
    of multiple characters. If we choose a dictionary for ``source``,
    we can just use keys/symbols consisting of multiple characters. But
    we have to make sure they all have the same length. If we choose a
    string for ``source``, we have to specify how many characters build
    one symbol::
    
        sage: str = "Split me into symbols consisting of 3 characters"
        sage: sf5 = ShannonFano(str, char_per_symbol=3)
        sage: for symbol, code in sf5.encoding_table().iteritems():
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
    
        sage: sf6 = ShannonFano("four", char_per_symbol=3)
        Traceback (most recent call last):
        ...
        ValueError: The passed string does not match with the passed value for char_per_symbol.

    When decoding long strings with a code consisting of long codewords,
    you can try to adapt the parameter ``decoding_table_key_len`` for
    a faster decoding. Or if you want to save memory when dealing with
    long codewords, you can set it to a smaller value::
        
        sage: sf7 = ShannonFano("Let's call it a day!", decoding_table_key_len=3)
    """
    
    def __init__(self, source, verbose=False, char_per_symbol=1,
                 decoding_table_key_len=8):
        r"""
        Constructor for ShannonFano.

        See the docstring of this class for full documentation.

        EXAMPLES::

            sage: from sage.coding.source_coding.shannon_fano import ShannonFano
            sage: str = "Give me an example!"
            sage: sf = ShannonFano(str, decoding_table_key_len=6)

        TESTS:

        Feeding anything else than a string or a dictionary::

            sage: ShannonFano(Graph())
            Traceback (most recent call last):
            ...
            ValueError: Input must be either a string or a dictionary.
        """
        self._v_table = SimpleTable(["level", "partition 1",
                                     "partition 2", "code"])
        PrefixCoding.__init__(self, source, verbose, 2, char_per_symbol,
                              decoding_table_key_len)
        
    def _sf_split(self, s_w, code, verbose):#
        r"""
        Helper method for building the Shannon-Fano code.

        INPUT:

        - ``s_w`` -- a list of tuples containing a symbol and a
          corresponding weight, sorted by weight.
          
        - ``code`` -- the codeword for the partition ``s_w``
          
        - ``verbose`` -- (default: False) if True, print intermediate
          data of building the code.
          
        EXAMPLES::
        
            sage: from sage.coding.source_coding.shannon_fano import ShannonFano
            sage: d = {'a': 5, 'b': 3.5, 'c': 1, 'd':0.4}
            sage: sf = ShannonFano(d)
            sage: l = [('a', 5), ('b', 3.5), ('c', 1), ('d', 0.4)]
            sage: sf._sf_split(l, '', False)
        """
        if len(s_w) == 1:
            self._character_to_code[s_w[0][0]] = code
            if verbose:
                self._v_table.add_row(["", s_w[0][0], "", code])
            return
        
        total_sum = sum(w for s, w in s_w)
        min_diff = total_sum
        min_diff_i = 0
        for i in range(1, len(s_w)):
            first_sum = sum(w for s, w in s_w[:i])
            second_sum = total_sum - first_sum
            diff = abs(first_sum - second_sum)
            if diff < min_diff:
                min_diff = diff
                min_diff_i = i
                
        if verbose:
            left = ', '.join(s for s, w in s_w[:min_diff_i])
            right = ', '.join(s for s, w in s_w[min_diff_i:])
            level = str(len(code) + 1)
            self._v_table.add_row([level, left, right, code])
        
        self._sf_split(s_w[:min_diff_i], code + '0', verbose)
        self._sf_split(s_w[min_diff_i:], code + '1', verbose)
        
    def _build_code(self, dic, verbose=False):
        r"""
        Construct a Shannon-Fano code corresponding to an alphabet with
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
        
            sage: from sage.coding.source_coding.shannon_fano import ShannonFano
            sage: d = {'c': 4, 'o': 3, 'd': 2, 'e':1}
            sage: sf = ShannonFano(d)
            sage: sf._build_code(d)
        """
        # Sort characters by decreasing probability.
        s_w = sorted(dic.items(), key=itemgetter(1), reverse=True)
        self._sf_split(s_w, '', verbose)
        if verbose:
            self._v_table.print_me()
        
    def encode(self, string):
        r"""
        Encode the given string based on the current encoding table.
        
        INPUT:

        - ``string`` -- a string of symbols over an alphabet.
          
        OUTPUT:

        - A Shannon-Fano encoding of ``string``.
        
        EXAMPLES:
        
        This example illustrates how to encode a string::
        
            sage: from sage.coding.source_coding.shannon_fano import ShannonFano
            sage: str = "Encode me!"
            sage: sf = ShannonFano(str)
            sage: encoded = sf.encode(str); encoded
            '11111110011111010000011010100010'
        """
        return PrefixCoding.encode(self, string)
    
    def decode(self, string):
        r"""
        Decode the given string based on the current decoding table.
        
        INPUT:

        - ``string`` -- a string of encoded symbols.
          
        OUTPUT:

        - The Shannon-Fano decoding of ``string``.
        
        EXAMPLES:
        
        This example illustrates how to encode and then decode a
        string::
        
            sage: from sage.coding.source_coding.shannon_fano import ShannonFano
            sage: str = "Encode me!"
            sage: sf = ShannonFano(str)
            sage: encoded = sf.encode(str); encoded
            '11111110011111010000011010100010'
            sage: sf.decode(encoded)
            'Encode me!'
        
        TESTS:
            
        It is an error to try to decode an encoding from another
        instance. This may raise an exception or result in
        a wrong decoding::
        
            sage: str2 = "I'm another string."
            sage: sf2 = ShannonFano(str2)
            sage: encoded_sf2 = sf2.encode(str2); encoded_sf2
            '10110101011000000111001110101101011110000100001111001100101110001100111111'
            sage: sf.decode(encoded_sf2)
            'mm!oeecenmmcdedeEeo!neocn'
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
        
            sage: from sage.coding.source_coding.shannon_fano import ShannonFano
            sage: str = "Show my encoding table!"
            sage: sf = ShannonFano(str)
            sage: t = sorted(sf.encoding_table().items())
            sage: for symbol, code in t:
            ....:     print(symbol, code)
            ....:     
              000
            ! 11011
            S 11010
            a 01110
            b 1000
            c 01111
            d 10010
            e 001
            g 10011
            h 10110
            i 1010
            l 1100
            m 10111
            n 0110
            o 010
            t 1110
            w 11110
            y 11111
        """
        return PrefixCoding.encoding_table(self)
    
    def tree(self):
        r"""
        Return the tree corresponding to the current encoding.

        INPUT:

        - None.

        OUTPUT:

        - The binary tree representing a Shannon-Fano code.

        EXAMPLES::
        
            sage: from sage.coding.source_coding.shannon_fano import ShannonFano
            sage: str = "Binary tree"
            sage: sf = ShannonFano(str)
            sage: t = sf.tree(); t
            Digraph on 17 vertices
            sage: t.show(layout='tree', vertex_size=1800, figsize=[8,8])
        """
        return PrefixCoding.tree(self)
