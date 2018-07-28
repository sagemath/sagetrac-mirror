r"""
Miscellaneous classes and functions

This module implements miscellaneous classes and functions for source
coding.

AUTHORS:

- Jan Wabbersen (2015-04-27): Initial version.
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

def frequency_table(string, char_per_symbol=1):
    r"""
    Return the frequency table corresponding to the given string.

    INPUT:

    - ``string`` -- a string of symbols over some alphabet.
    
    - ``char_per_symbol`` -- the number of characters
      that define one symbol.

    OUTPUT:

    - A table of frequencies of each unique symbol in ``string``. If
      ``string`` is an empty string, return an empty table.

    EXAMPLES:

    The frequency table of a non-empty string::

        sage: from sage.coding.source_coding.misc import frequency_table
        sage: str = "Stop counting my characters!"
        sage: t = sorted(frequency_table(str).items())
        sage: for symbol, code in t:
        ....:     print(symbol, code)
        ....:     
          3
        ! 1
        S 1
        a 2
        c 3
        e 1
        g 1
        h 1
        i 1
        m 1
        n 2
        o 2
        p 1
        r 2
        s 1
        t 3
        u 1
        y 1

    The frequency of an empty string::

        sage: frequency_table("")
        {}
        
    The frequency table of a string with two characters defined as one
    symbol::
    
        sage: t2 = sorted(frequency_table(str, char_per_symbol=2).items())
        sage: for symbol, code in t2:
        ....:     print(symbol, code)
        ....:     
         c 2
        St 1
        ct 1
        er 1
        g  1
        ha 1
        in 1
        my 1
        nt 1
        op 1
        ou 1
        ra 1
        s! 1
    """
    if char_per_symbol != 1 and len(string) % char_per_symbol != 0:
                raise ValueError("The passed string does not match with the "
                                 "passed value for char_per_symbol.")
    d = {} 
    for i in range(0, len(string), char_per_symbol):
        s = string[i:i+char_per_symbol]
        d[s] = d.get(s, 0) + 1
    return d

def cum_sum(seq):
    r"""
    Return the cumulative sums of the given sequence.

    INPUT:

    - ``string`` -- sequence of numbers.

    OUTPUT:

    - A list of the cumulative sums of ``seq``.

    EXAMPLES::

        sage: from sage.coding.source_coding.misc import cum_sum
        sage: l = [3, 5, 5, 1]
        sage: cum_sum(l)
        [3, 8, 13, 14]
    """
    cum_sum = [seq[0]]
    for i, a in enumerate(seq[1:]):
        cum_sum.append(cum_sum[i] + a)
    return cum_sum

def first_binary_dec_places(number, x):
    r"""
    Return the first ``x`` binary fractional digits of ``number``.

    INPUT:

    - ``number`` -- the number between 0 and 1 inclusive.
    
    - ``x`` -- the number of binary places after the comma.

    OUTPUT:

    - The first ``x`` binary fractional digits of ``number``.

    EXAMPLES::

        sage: from sage.coding.source_coding.misc import first_binary_dec_places
        sage: number = 0.25
        sage: first_binary_dec_places(number, 5)
        '01000'
    """
    b = ''
    for _ in range(x):
        number *= 2
        if number < 1:
            b += '0'
        else:
            b += '1'
            number -= 1
    return b


class SimpleTable():
    r"""
    Saves elements organized in a table with columns and rows.
    
    INPUT:

        - ``cols`` -- a list of column names.
       
    EXAMPLES::
    
        sage: from sage.coding.source_coding.misc import SimpleTable
        sage: l = ["column one", "column two", "column three"]
        sage: t = SimpleTable(l)
        sage: row = [33, "alpha", 1]
        sage: t.add_row(row)
        sage: t.print_me()
        | column one | column two | column three |
        | 33         | alpha      | 1            |
    """
    def __init__(self, cols):
        r"""
        Constructor for SimpleTable.

        See the docstring of this class for full documentation.

        EXAMPLES::

            sage: from sage.coding.source_coding.misc import SimpleTable
            sage: l = ["column one", "column two", "column three"]
            sage: t = SimpleTable(l)
        """
        self._number_of_columns = len(cols)
        self._cols = [cols]
        
    def add_row(self, row):
        r"""
        Add a row.
        
        INPUT:

        - ``row`` -- the row to be added.
        
        EXAMPLES::
        
            sage: from sage.coding.source_coding.misc import SimpleTable
            sage: l = ["column one", "column two", "column three"]
            sage: t = SimpleTable(l)
            sage: row = [33, "alpha", 1]
            sage: t.add_row(row)
        """
        if self._number_of_columns != len(row):
            raise ValueError("The number of elements in the passed row does "
                             "not match with the table!")
        self._cols.append([str(s) for s in row])
            
    def print_me(self):
        r"""
        Print the current table.
        
        EXAMPLES::
        
            sage: from sage.coding.source_coding.misc import SimpleTable
            sage: l = ["column one", "column two", "column three"]
            sage: t = SimpleTable(l)
            sage: row = [33, "alpha", 1]
            sage: t.add_row(row)
            sage: t.print_me()
            | column one | column two | column three |
            | 33         | alpha      | 1            |
        """
        col_width = [max(len(e) for e in col) for col in zip(*self._cols)]
        for line in self._cols:
            print("| " + " | ".join("{:{}}".format(e, col_width[i])
                                    for i, e in enumerate(line)) + " |")
