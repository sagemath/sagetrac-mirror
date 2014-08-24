"""
Complete Intersection Calabi-Yau Manifolds.

Complete Intersection Calabi-Yau Manifolds (CICY) in Products of
Projective Spaces are a rather simple class of Calabi-Yau manifolds.
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import re
import os
from sage.env import SAGE_SHARE
from sage.structure.sage_object import SageObject
from sage.rings.all import ZZ, QQ
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method


MMA_TEST_1 = """
{{Num -> 1, Size -> {7, 6}, DimList -> {1, 1, 1, 1, 1, 2, 2}, Chi -> 0, 
  AllInds -> {0, 38, 1494, 32262, -5472}, HodgeNos -> {15, 15}, 
  C2 -> {24, 24, 24, 24, 24, 36, 36}, 
  Conf -> {{1, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 1}, {0, 0, 0, 0, 1, 1}, 
    {1, 0, 0, 1, 0, 0}, {1, 0, 0, 0, 0, 1}, {0, 0, 1, 2, 0, 0}, 
    {0, 1, 0, 0, 2, 0}}}, {Num -> 2, Size -> {7, 6}, 
  DimList -> {1, 1, 1, 1, 1, 2, 2}, Chi -> 0, 
  AllInds -> {0, 38, 1494, 32154, -5472}, HodgeNos -> {15, 15}, 
  C2 -> {24, 24, 24, 24, 24, 36, 36}, 
  Conf -> {{1, 1, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 1}, {0, 0, 1, 0, 1, 0}, 
    {0, 0, 1, 1, 0, 0}, {1, 0, 1, 0, 0, 0}, {0, 0, 0, 2, 0, 1}, 
    {0, 1, 0, 0, 2, 0}}}, {Num -> 3, Size -> {6, 6}, 
  DimList -> {1, 1, 1, 2, 2, 2}, Chi -> 0, 
  AllInds -> {0, 39, 1638, 36348, -5616}, HodgeNos -> {15, 15}, 
  C2 -> {24, 24, 24, 36, 36, 36}, 
  Conf -> {{1, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 1}, {0, 0, 0, 1, 1, 0}, 
    {2, 0, 0, 1, 0, 0}, {0, 1, 0, 0, 0, 2}, {0, 0, 1, 0, 2, 0}}}}
"""

MMA_TEST_2 = """
{{1, 1, 1, {{6}}, 2610, False, 1, 0, 426, 1752, True}, 
 {2, 1, 2, {{2, 5}}, 2190, False, 1, 0, 356, 1472, True}, 
 {3, 1, 2, {{3, 4}}, 1476, False, 1, 0, 237, 996, True}, 
 {4, 1, 3, {{2, 2, 4}}, 1632, False, 1, 0, 263, 1100, True}, 
 {5, 1, 3, {{2, 3, 3}}, 1206, False, 1, 0, 192, 816, True}}
"""


class MMaList(SageObject):

    NUMBER_RE = re.compile('^[+-]?[0-9]+$')      # only deal with integers for now
    IDENTIFIER_RE = re.compile('^[a-zA-Z][a-zA-Z0-9]*$')
    
    def __init__(self, stream):
        r"""
        Read Mathematica list and iterate over entries.

        INPUT:

        - ``stream`` -- input stream.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import (
            ....:     MMaList, MMA_TEST_1, MMA_TEST_2)

            sage: ii = iter(MMaList(MMA_TEST_1))
            sage: next(ii)    # random output
            {'AllInds': [0, 38, 1494, 32262, -5472], 'Chi': 0, 'HodgeNos': [15, 15], 
             'DimList': [1, 1, 1, 1, 1, 2, 2], 'Num': 1, 
             'Conf': [[1, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 1, 1], 
                      [1, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 1], [0, 0, 1, 2, 0, 0], 
                      [0, 1, 0, 0, 2, 0]], 
             'C2': [24, 24, 24, 24, 24, 36, 36], 
             'Size': [7, 6]}

            sage: ii = iter(MMaList(MMA_TEST_2))
            sage: next(ii)
            [1, 1, 1, [[6]], 2610, False, 1, 0, 426, 1752, True]
            sage: next(ii)
            [2, 1, 2, [[2, 5]], 2190, False, 1, 0, 356, 1472, True]
            sage: next(ii)
            [3, 1, 2, [[3, 4]], 1476, False, 1, 0, 237, 996, True]
        """
        from StringIO import StringIO
        if isinstance(stream, basestring):
            self.stream = StringIO(stream)
        else:
            self.stream = stream
        self.curr = None
        self.ch = None
        self.tok = None
        
    def getchar(self):
        """
        Return the next char from the stream and advance.

        OUTPUT:

        A single character. Empty string if end of file.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import MMaList
            sage: lst = MMaList('{test123}')
            sage: lst.getchar()
            '{'
            sage: lst.getchar()
            't'
            sage: lst.getchar()
            'e'
            sage: lst.getchar()
            's'
        """
        if self.ch is not None:
            try:
                return self.ch
            finally:
                self.ch = None
        else:
            return self.stream.read(1)

    def peek(self):
        """
        Return the next char from the stream without advancing.

        OUTPUT:

        A single character. Empty string if end of file.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import MMaList
            sage: lst = MMaList('{test123}')
            sage: lst.peek()
            '{'
            sage: lst.peek()
            '{'
            sage: lst.peek()
            '{'
        """
        if self.ch is None:
            self.ch = self.stream.read(1)
        return self.ch

    def get_token(self):
        """
        Get the next token.

        OUTPUT:

        Returns the next token as string. Also stores the token in the
        ``self.tok`` attribute.
        
        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import (
            ....:     MMaList, MMA_TEST_1)
            sage: lst = MMaList(MMA_TEST_1)
            sage: lst.get_token()
            '{'
            sage: lst.get_token()
            '{'
            sage: lst.get_token()
            'Num'
            sage: lst.get_token()
            '->'
        """
        while True:
            tok = self.getchar()
            if len(tok) == 0:
                raise ValueError('unexpected end of stream')
            if tok not in [' ', '\n']:
                break
        if tok in ['{', '}', ',']:
            self.tok = tok
            return tok
        while self.peek() not in ['}', ' ', ',', '\n']:
            ch = self.getchar()
            if len(ch) == 0:
                raise ValueError('unexpected end of stream')            
            tok += ch
        self.tok = tok
        return tok

    def get_entry(self):
        """
        Return the next container entry.

        OUTPUT:
        
        Pair consisting of a key and a value, corresponding to the MMa
        notation ``key -> value``. If there is no arrow in the entry
        then only the value is returned (key is ``None``).

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import MMaList
            sage: lst = MMaList('{123, True, foo -> 456}')
            sage: lst.get_token()
            '{'
            sage: lst.tok
            '{'
            sage: lst.get_entry()
            (None, 123)
            sage: lst.tok
            ','
            sage: lst.get_entry()
            (None, True)
            sage: lst.tok
            ','
            sage: lst.get_entry()
            ('foo', 456)
            sage: lst.tok
            '}'
        """
        assert self.tok in ['{', ',', '->'], 'not at list entry starting point'
        tok = self.get_token()
        if tok == 'True':
            self.get_token()
            return None, True
        if tok == 'False':
            self.get_token()
            return None, False
        if tok == 'Null':
            self.get_token()
            return None, None
        if tok == '{':
            lst = self.get_container()
            assert self.tok == '}'
            self.get_token()
            return None, lst
        if self.NUMBER_RE.match(tok):
            self.get_token()
            return None, int(tok)
        if self.IDENTIFIER_RE.match(tok):
            next_tok = self.get_token()
            if next_tok == '->':
                key, value = self.get_entry()
                if key is not None:
                    raise ValueError('nested keys not allowed: {0} -> {1} -> ...'
                                     .format(tok, key))
                return tok, value
            else:
                raise ValueError('unknown identifier: ' + tok)

    def iter_container(self):
        """
        Return the next list/dict entry.

        OUTPUT:

        Iterate over pairs ``(key, value)`` corresponding to the
        container entries.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import MMaList
            sage: lst = MMaList('{123, True}')
            sage: lst.get_token()
            '{'
            sage: ii = lst.iter_container()
            sage: ii.next()
            (None, 123)
            sage: ii.next()
            (None, True)
            sage: ii.next()
            Traceback (most recent call last):
            ...
            StopIteration
            sage: lst.tok
            '}'
        """
        assert self.tok == '{'
        while True:
            yield self.get_entry()
            if self.tok == '}':
                raise StopIteration
        
    def get_container(self):
        """
        Return the next container as list or dict.

        OUTPUT:

        A Python list or dict.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import MMaList
            sage: lst = MMaList('{123, True, 456}')
            sage: lst.get_token()
            '{'
            sage: lst.get_container()
            [123, True, 456]
            sage: lst.tok
            '}'

            sage: dct = MMaList('{foo -> 123, bar -> True, baz -> 456}')
            sage: dct.get_token()
            '{'
            sage: result = dct.get_container();   result   # random output
            sage: type(result)
            <type 'dict'>
            sage: sorted(result.items())
            [('bar', True), ('baz', 456), ('foo', 123)]
            sage: dct.tok
            '}'
        """
        keys = []
        values = []
        assert self.tok == '{'
        for key, value in self.iter_container():
            keys.append(key)
            values.append(value)
        if all(k is None for k in keys):
            return values
        if any(k is None for k in keys):
            raise ValueError('cannot mix list and dict')
        return dict(zip(keys, values))

    def __iter__(self):
        """
        Iterate over the top-level list entries.

        OUTPUT:
        
        Yields the list entries converted to Python objects.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import MMaList
            sage: iter(MMaList('{{foo -> {bar ->\n -123}}, 1}')).next()
            {'foo': {'bar': -123}}        
        """
        tok = self.get_token()
        for key, value in self.iter_container():
            if key is not None:
                raise ValueError('top-level is not a list but a dict')
            yield value


class Database(SageObject):
    
    DATABASE_3D = os.path.join(SAGE_SHARE, 'complete_intersection_calabi_yau', 'cicy3list.txt')
    DATABASE_4D = os.path.join(SAGE_SHARE, 'complete_intersection_calabi_yau', 'cicy4list.zip')

    def __init__(self, dimension, output='complete_intersection', base_ring=QQ):
        r"""
        Read a database of CICYs

        INPUT:

        - ``dimension`` -- integer. The Calabi-Yau dimension of
          interest.

        - ``output`` -- string. How to output each CICY. Allowed
          values are

          * ``'complete_intersection'`` -- (default). Output
            :class:`~sage.schemes.product_projective.complete_intersection.GenericCompleteIntersection`
            objects.

          * ``'matrix'`` -- output matrix. Standard notation: rows are
            projective space factors and columns are equations.

          * ``'list'`` -- output list of list of integers. Fastest
            output.

        - ``base_ring`` -- commutative ring (default: `\QQ`). The base
          ring if ``output`` is ``'complete_intersection'``.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import Database
            sage: db = Database(5)        # optional - database_cicy
            Traceback (most recent call last):
            ...
            ValueError: only have databases for 3d and 4d
        """
        if dimension == 3:
            if not os.path.exists(self.DATABASE_3D):
                raise RuntimeError('The database_cicy optional package is not installed')
            self._iter = self._iter_3d 
        elif dimension == 4:
            if not os.path.exists(self.DATABASE_4D):
                raise RuntimeError('The database_cicy optional package is not installed')
            self._iter = self._iter_4d
        else:
            raise ValueError('only have databases for 3d and 4d')
        if output == 'complete_intersection':
            self._output = self._output_complete_intersection
        elif output == 'list':
            self._output = self._output_list
        elif output == 'matrix':
            self._output = self._output_matrix
        self._base_ring = base_ring

    def _iter_3d(self):
        """
        Iterate over the CICY treefolds.
        
        OUTPUT:
        
        Iterator over pairs consisting of a configuration matrix and
        additional data.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import Database
            sage: db = Database(3)        # optional - database_cicy
            sage: db._iter_3d().next()    # optional - database_cicy    # random output
            ([[1, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 1, 1], [1, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 1], [0, 0, 1, 2, 0, 0], [0, 1, 0, 0, 2, 0]], {'AllInds': [0, 38, 1494, 32262, -5472], 'Chi': 0, 'HodgeNos': [15, 15], 'DimList': [1, 1, 1, 1, 1, 2, 2], 'Num': 1, 'Conf': [[1, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 1, 1], [1, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 1], [0, 0, 1, 2, 0, 0], [0, 1, 0, 0, 2, 0]], 'C2': [24, 24, 24, 24, 24, 36, 36], 'Size': [7, 6]})
        """
        import io
        db = io.open(self.DATABASE_3D, 'rb')
        for cicy in MMaList(db):
            yield cicy['Conf'], cicy
        db.close()

    def _iter_4d(self):
        """
        Iterate over the CICY treefolds.
        
        OUTPUT:
        
        Iterator over pairs consisting of a configuration matrix and
        additional data.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import Database
            sage: db = Database(4)         # optional - database_cicy
            sage: db._iter_4d().next()     # optional - database_cicy    # random output
            ([[6]], [1, 1, 1, [[6]], 2610, False, 1, 0, 426, 1752, True])
        """
        from zipfile import ZipFile
        z = ZipFile(self.DATABASE_4D)
        db = z.open('cicy4list.txt')
        for cicy in MMaList(db):
            yield cicy[3], cicy
        db.close()
    
    def _output_list(self, listlist):
        """
        Construct ``'list'``-style output.

        INPUT:

        - ``listlist`` -- list of lists of integers. The CICY
          configuration.

        OUTPUT:

        Same.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import Database
            sage: db = Database(3)                       # optional - database_cicy
            sage: db._output_list([[1,2,3], [3,4,5]])    # optional - database_cicy
            [[1, 2, 3], [3, 4, 5]]
        """
        return listlist

    def _output_matrix(self, listlist):
        """
        Construct ``'matrix'``-style output.

        INPUT:

        - ``listlist`` -- list of lists of integers. The CICY
          configuration.

        OUTPUT:

        Matrix.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import Database
            sage: db = Database(3)                         # optional - database_cicy
            sage: db._output_matrix([[1,2], [2,1]])    # optional - database_cicy
            [1 2]
            [2 1]
        """
        return matrix(ZZ, listlist)

    @cached_method
    def _product_projective_spaces(self, dimensions):
        """
        Construct product of projective spaces.

        INPUT:

        - ``dimensions`` -- tuple of integers. The projective space
          dimensions.

        OUTPUT:

        A
        :class:`~sage.schemes.product_projective.space.ProductProjectiveSpaces_ring`
        instance.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import Database
            sage: db = Database(3)                          # optional - database_cicy
            sage: db._product_projective_spaces((1,2,3))    # optional - database_cicy
            Product of projective spaces P^1 x P^2 x P^3 over Rational Field
        """
        assert isinstance(dimensions, tuple)
        from sage.schemes.product_projective.space import ProductProjectiveSpaces
        return ProductProjectiveSpaces(dimensions, self._base_ring)

    def _output_complete_intersection(self, listlist):
        """
        Construct ``'complete_intersection'``-style output.

        INPUT:

        - ``listlist`` -- list of lists of integers. The CICY
          configuration.

        OUTPUT:

        A
        :class:`~sage.schemes.product_projective.complete_intersection.GenericCompleteIntersection`
        instance.

        EXAMPLES::

            sage: from sage.schemes.product_projective.cicy import Database
            sage: db = Database(3)                                    # optional - database_cicy
            sage: db._output_complete_intersection([[1,2], [2,1]])    # optional - database_cicy
            Complete intersection in Product of projective spaces P^2 x P^2 over Rational Field
              P^2   |   1   2
              P^2   |   2   1
        """
        dims = tuple(sum(degrees) - 1 for degrees in listlist)
        ambient = self._product_projective_spaces(dims)
        return ambient.complete_intersection(*listlist)

    def __iter__(self):
        """
        Iterate over the entries of the database.

            sage: from sage.schemes.product_projective.cicy import Database
            sage: db = Database(3)     # optional - database_cicy
            sage: next(iter(db))       # optional - database_cicy
            [[1, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 1, 1], 
             [1, 0, 0, 1, 0, 0], [1, 0, 0, 0, 0, 1], [0, 0, 1, 2, 0, 0], 
             [0, 1, 0, 0, 2, 0]]

            sage: db = Database(4)     # optional - database_cicy
            sage: next(iter(db))       # optional - database_cicy
            [[6]]
        """
        for conf, data in self._iter():
            yield self._output(conf)
        
