r"""
Huffman Coding

This module implements functionalities relating to Huffman encoding
and decoding.

AUTHORS:

- Nathann Cohen (2010-05): initial version.

- Jan Wabbersen (2015-04-27): adapted to new super class and generalized
  to q-nary huffman.
"""

#*****************************************************************************
#       Copyright (C) 2015 Jan Wabbersen <jan.wabbersen@googlemail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# http://www.gnu.org/licenses/
###########################################################################
from __future__ import print_function

from heapq import heapify, heappush, heappop

from prefix_coding import PrefixCoding
from misc import SimpleTable

class Huffman(PrefixCoding):
    r"""
    This class implements basic functionalities of Huffman coding.
    
    It can build a q-nary Huffman code from a given string, or from the
    information of a dictionary associating to each key (the elements
    of the alphabet) a weight (most of the time, a probability value
    or a number of occurrences).
    
    INPUT:

        - ``source`` -- can be either

            - A string from which the Huffman encoding should be
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
          
        - ``q`` -- (default: 2) the number of symbols in the code
          alphabet. Supported is 2 <= q <= 10.
          
        - ``char_per_symbol`` -- (default: 1) the number of characters
          that define one symbol.
          
        - ``decoding_table_key_len`` -- (default: 8) the length of each
          key in the resulting decoding table. Longer keys may result in
          faster decoding. Decoding tables with shorter keys consume
          less memory, but decoding may be slower.
    
    In order to construct a Huffman encoding for an alphabet, we
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
    
        sage: from sage.coding.source_coding.huffman import Huffman
        sage: h1 = Huffman("Encode me!")
        sage: for symbol, code in h1.encoding_table().iteritems():
        ....:     print "'" + symbol + "' : " + code
        ....:     
        '!' : 1101
        ' ' : 1100
        'c' : 1111
        'e' : 10
        'd' : 000
        'm' : 001
        'o' : 011
        'n' : 010
        'E' : 1110

    Now, we have a Huffman code for encoding strings consisting of
    symbols that are in the encoding table::
    
        sage: encoded = h1.encode("Encode me!"); encoded
        '11100101111011000101100001101101'
        
    And we can decode the encoded string in the following way::
    
        sage: h1.decode(encoded)
        'Encode me!'

    If we try to encode a string with ``h1`` that consists of symbols
    which are not part of the encoding table, we of course get an
    error::
    
        sage: h1.encode("This string contains other symbols!")
        Traceback (most recent call last):
        ...
        KeyError: 'T'

    If we want to encode multiple strings over an alphabet, it might be
    a good idea to generate a Huffman code for this alphabet and a
    given probability distribution for the symbols::
    
        sage: d = {'a': 0.5, 'b': 0.25, 'c': 0.2, 'd': 0.05}
        sage: h2 = Huffman(d)
        sage: for symbol, code in h2.encoding_table().iteritems():
        ....:     print "'" + symbol + "' : " + code
        ....:     
        'a' : 1
        'c' : 001
        'b' : 01
        'd' : 000

    We can also use frequencies or in general weights for the symbols::
    
        sage: d2 = {'a': 2.5, 'b': 1.25, 'c': 1, 'd': 0.25}
        sage: h3 = Huffman(d2)
        sage: for symbol, code in h3.encoding_table().iteritems():
        ....:     print "'" + symbol + "' : " + code
        ....:     
        'a' : 1
        'c' : 001
        'b' : 01
        'd' : 000
    
    If you are interested in Huffman coding, you can output intermediate
    data of the generation of the Huffman tree::
    
        sage: h4 = Huffman("huffman", verbose=True)
        | merge | tree 0      | tree 1      |
        | 1     | a           | h           |
        | 2     | m           | n           |
        | 3     | u           | [a, h]      |
        | 4     | [m, n]      | f           |
        | 5     | [u, [a, h]] | [[m, n], f] |
        
    If we want to encode with a larger code alphabet, we can use the
    parameter ``q``::
        
        sage: str = "Sometimes a larger code alphabet is needed."
        sage: h5 = Huffman(str, q=3)
        sage: h5.encode(str)
        '001122121222021021212220121112112011200011222002101012202222111120101012110022220221102201211002222022202000'
        sage: h6 = Huffman(str, q=10)
        sage: h6.encode(str)
        '9110849808376799629482793158769997956928479837968858590'
        
    With ``char_per_symbol`` we are able to define symbols consisting
    of multiple characters. If we choose a dictionary for ``source``,
    we can just use keys/symbols consisting of multiple characters. But
    we have to make sure they all have the same length. If we choose a
    string for ``source``, we have to specify how many characters build
    one symbol::
    
        sage: str2 = "Split me into symbols consisting of 3 characters"
        sage: h6 = Huffman(str2, char_per_symbol=3)
        sage: for symbol, code in h6.encoding_table().iteritems():
        ....:     print "'" + symbol + "' : " + code
        ....:     
        'nsi' : 1010
        'me ' : 1000
        'int' : 0110
        '3 c' : 0001
        'it ' : 0111
        'ols' : 1101
        'of ' : 1100
        'ng ' : 1001
        ' co' : 0000
        'ers' : 0100
        'sti' : 1110
        'act' : 0011
        'ymb' : 1111
        'Spl' : 0010
        'o s' : 1011
        'har' : 0101
        
    Obviously, we have to make sure that the length of the string is a
    multiple of ``char_per_symbol`` or else we get an error::
    
        sage: h7 = Huffman("four", char_per_symbol=3)
        Traceback (most recent call last):
        ...
        ValueError: The passed string does not match with the passed value for char_per_symbol.

    When decoding long strings with a code consisting of long codewords,
    you can try to adapt the parameter ``decoding_table_key_len`` for
    a faster decoding. Or if you want to save memory when dealing with
    long codewords, you can set it to a smaller value::
        
        sage: h8 = Huffman("Let's call it a day!", decoding_table_key_len=3)
    """
    
    def __init__(self, source, verbose=False, q=2, char_per_symbol=1,
                 decoding_table_key_len=8):
        r"""
        Constructor for Huffman.

        See the docstring of this class for full documentation.

        EXAMPLES::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Give me an example!"
            sage: h = Huffman(str, q=3, decoding_table_key_len=6)

        TESTS:

        Feeding anything else than a string or a dictionary::

            sage: Huffman(Graph())
            Traceback (most recent call last):
            ...
            ValueError: Input must be either a string or a dictionary.
        """
        PrefixCoding.__init__(self, source, verbose, q, char_per_symbol,
                              decoding_table_key_len)
        
    def _build_code_from_tree(self, tree, q, prefix):
        r"""
        Build the Huffman code corresponding to a given tree, prefix
        and code alphabet size.

        INPUT:

        - ``tree`` -- list of list with maximum q strings

        - ``q`` -- the code alphabet size

        - ``prefix`` (string) -- string which is the prefix
          of any element of the tree

        EXAMPLES::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: h = Huffman("Some string.")
            sage: tree = [['a', 'b'], 'c', 'd']  
            sage: h._build_code_from_tree(tree, 3, "")
        """
        # This is really a recursive construction of a Huffman code. By
        # feeding this class a sufficiently large alphabet, it is possible to
        # exceed the maximum recursion depth and hence result in a
        # RuntimeError.
        try:
            if isinstance(tree, list):
                for i in range(q):
                    self._build_code_from_tree(tree[i], q, prefix + str(i))
            else:
                self._character_to_code[tree] = prefix
        # For q > 2, there can be a branch with fewer than q leaves.
        except IndexError:
            return

    def _build_code(self, dic, verbose=False, q=2):
        r"""
        Constructs a Huffman code corresponding to an alphabet with the
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
          
        - ``verbose`` -- (default: False) if True, print intermediate
          data of building the code.
          
        - ``q`` -- (default: 2) the number of symbols in the code
          alphabet. Supported is 2 <= q <= 10.

        EXAMPLES::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: d = {'c': 4, 'o': 3, 'd': 2, 'e':1}
            sage: h = Huffman(d)
            sage: h._build_code(d)
        """
        if verbose:
            headings = ["merge"] + ["tree " + str(i) for i in range(q)]
            table = SimpleTable(headings)
            merge = 0
            
        # Switch s and w for heapify.
        heap = [(w, s) for s, w in dic.items()]
        heapify(heap)
        # Number of trees to merge.
        trees = q if q == 2 else len(dic) % (q-1)
        if q != 2 and trees <= 1:
            trees = q - 1 + trees
        # Construct a Huffman tree.
        while len(heap) > 1:
            weight = 0
            tree_list = []
            for _ in range(trees):
                w, t = heappop(heap)
                weight += w
                tree_list.append(t)
            heappush(heap, (weight, tree_list))
            if verbose:
                merge += 1
                row = [str(merge)] + \
                      [str(t).replace("'", "") for t in tree_list] + \
                      ["" for _ in range(trees, q)]
                table.add_row(row)   
            trees = q
        tree = heap[0][1]
        self._build_code_from_tree(tree, q, prefix="")
        if verbose:
            table.print_me()
        
    def encode(self, string):
        r"""
        Encode the given string based on the current encoding table.
        
        INPUT:

        - ``string`` -- a string of symbols over an alphabet.
          
        OUTPUT:

        - A Huffman encoding of ``string``.
        
        EXAMPLES:
        
        This example illustrates how to encode a string::
        
            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Encode me!"
            sage: h = Huffman(str) 
            sage: h.encode(str)
            '11100101111011000101100001101101'
        """
        return PrefixCoding.encode(self, string)
    
    def decode(self, string):
        r"""
        Decode the given string based on the current decoding table.
        
        INPUT:

        - ``string`` -- a string of encoded symbols.
          
        OUTPUT:

        - The Huffman decoding of ``string``.
        
        EXAMPLES:
        
        This example illustrates how to encode and then decode a
        string::
        
            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Encode me!"
            sage: h = Huffman(str)
            sage: encoded = h.encode(str); encoded
            '11100101111011000101100001101101'
            sage: h.decode(encoded)
            'Encode me!'
        
        TESTS:
            
        It is an error to try to decode an encoding from another
        instance. This may raise an exception or result in
        a wrong decoding::
        
            sage: str2 = "I'm another string."
            sage: h2 = Huffman(str2)
            sage: encoded_h2 = h2.encode(str2); encoded_h2
            '11100110100010010111010110011101000011110100010110010110000010111111111011'
            sage: h.decode(encoded_h2)
            'Eonmmoee EedcnmomodmocE '
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
        
            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Show my encoding table!"
            sage: h = Huffman(str)
            sage: t = sorted(h.encoding_table().items())
            sage: for symbol, code in t:
            ....:     print symbol, code
            ....:     
              100
            ! 10100
            S 10101
            a 10110
            b 10111
            c 11000
            d 11001
            e 001
            g 11010
            h 11011
            i 11100
            l 11101
            m 11110
            n 010
            o 011
            t 11111
            w 0000
            y 0001
        """
        return PrefixCoding.encoding_table(self)
    
    def tree(self):
        r"""
        Return the Huffman tree corresponding to the current encoding.

        INPUT:

        - None.

        OUTPUT:

        - The tree representing a Huffman code.

        EXAMPLES::
        
            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Binary tree"
            sage: h = Huffman(str, q=3)
            sage: t = h.tree(); t
            Digraph on 13 vertices
            sage: t.show(layout='tree', vertex_size=1800, figsize=[8,8])
        """
        return PrefixCoding.tree(self)
