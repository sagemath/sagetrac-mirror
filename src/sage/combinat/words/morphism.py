# coding=utf-8
r"""
Word morphisms/substitutions

This modules implements morphisms over finite and infinite words.

AUTHORS:

- Sebastien Labbe (2007-06-01): initial version
- Sebastien Labbe (2008-07-01): merged into sage-words
- Sebastien Labbe (2008-12-17): merged into sage
- Sebastien Labbe (2009-02-03): words next generation
- Sebastien Labbe (2009-11-20): allowing the choice of the
  datatype of the image. Doc improvements.
- Stepan Starosta (2012-11-09): growing letters

EXAMPLES:

Creation of a morphism from a dictionary or a string::

    sage: n = WordMorphism({0:[0,2,2,1],1:[0,2],2:[2,2,1]})

::

    sage: m = WordMorphism('x->xyxsxss,s->xyss,y->ys')

::

    sage: n
    WordMorphism: 0->0221, 1->02, 2->221
    sage: m
    WordMorphism: s->xyss, x->xyxsxss, y->ys

The codomain may be specified::

    sage: WordMorphism({0:[0,2,2,1],1:[0,2],2:[2,2,1]}, codomain=Words([0,1,2,3,4]))
    WordMorphism: 0->0221, 1->02, 2->221

Power of a morphism::

    sage: n^2
    WordMorphism: 0->022122122102, 1->0221221, 2->22122102

Image under a morphism::

    sage: m('y')
    word: ys
    sage: m('xxxsy')
    word: xyxsxssxyxsxssxyxsxssxyssys

Iterated image under a morphism::

    sage: m('y', 3)
    word: ysxyssxyxsxssysxyssxyss

Infinite fixed point of morphism::

    sage: fix = m.fixed_point('x')
    sage: fix
    word: xyxsxssysxyxsxssxyssxyxsxssxyssxyssysxys...
    sage: fix.length()
    +Infinity

Incidence matrix::

    sage: matrix(m)
    [2 3 1]
    [1 3 0]
    [1 1 1]

Many other functionalities...::

    sage: m.is_identity()
    False
    sage: m.is_endomorphism()
    True
"""
#*****************************************************************************
#       Copyright (C) 2008 Sebastien Labbe <slabqc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import itertools
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.sets.set import Set
from sage.rings.all import QQ
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer import Integer
from sage.modules.free_module_element import vector
from sage.matrix.constructor import Matrix
from sage.combinat.words.word import FiniteWord_class
from sage.combinat.words.words import Words_all, Words

class CallableDict(dict):
    r"""
    Wrapper of dictionary that makes it callable.

    EXAMPLES::

        sage: from sage.combinat.words.morphism import CallableDict
        sage: d = CallableDict({1:'one', 2:'zwei', 3:'trois'})
        sage: d(1), d(2), d(3)
        ('one', 'zwei', 'trois')
    """
    def __call__(self, key):
        r"""
        Returns the value with key ``key``.

        EXAMPLES::

            sage: from sage.combinat.words.morphism import CallableDict
            sage: d = CallableDict({'one': 1, 'zwei': 2, 'trois': 3})
            sage: d('one'), d('zwei'), d('trois')
            (1, 2, 3)
        """
        return self[key]

def get_cycles(f, domain=None):
    r"""
    Return the cycle of the function ``f`` on the finite set domain. It is
    assumed that f is an endomorphism.

    INPUT:

    - ``f`` - function.

    - ``domain`` - set (default: None) - the domain of ``f``. If none, then
      tries to use ``f.domain()``.

    EXAMPLES::

        sage: from sage.combinat.words.morphism import get_cycles
        sage: get_cycles(lambda i: (i+1)%3, domain=[0,1,2])
        [(0, 1, 2)]
        sage: get_cycles(lambda i: [0,0,0][i], domain=[0,1,2])
        [(0,)]
        sage: get_cycles(lambda i: [1,1,1][i], domain=[0,1,2])
        [(1,)]
    """
    if domain is None:
        try:
            domain = f.domain()
        except AttributeError:
            raise ValueError("you should specify the domain of the function f")
    cycles = []
    not_seen = dict((letter,True) for letter in domain)
    for a in not_seen:
        if not_seen[a]:
            not_seen[a] = False
            cycle = [a]
            b = f(a)
            while not_seen[b]:
                not_seen[b] = False
                cycle.append(b)
                b = f(b)
            if b in cycle:
                cycles.append(tuple(cycle[cycle.index(b):]))

    return cycles

class WordMorphism(SageObject):
    r"""
    WordMorphism class

    EXAMPLES::

        sage: n = WordMorphism({0:[0,2,2,1],1:[0,2],2:[2,2,1]})
        sage: m = WordMorphism('x->xyxsxss,s->xyss,y->ys')

    Power of a morphism::

        sage: n^2
        WordMorphism: 0->022122122102, 1->0221221, 2->22122102

    Image under a morphism::

        sage: m('y')
        word: ys
        sage: m('xxxsy')
        word: xyxsxssxyxsxssxyxsxssxyssys

    Iterated image under a morphism::

        sage: m('y', 3)
        word: ysxyssxyxsxssysxyssxyss

    See more examples in the documentation of the call method
    (``m.__call__?``).

    Infinite fixed point of morphism::

        sage: fix = m.fixed_point('x')
        sage: fix
        word: xyxsxssysxyxsxssxyssxyxsxssxyssxyssysxys...
        sage: fix.length()
        +Infinity

    Incidence matrix::

        sage: matrix(m)
        [2 3 1]
        [1 3 0]
        [1 1 1]

    Many other functionalities...::

        sage: m.is_identity()
        False
        sage: m.is_endomorphism()
        True

    TESTS::

        sage: wm = WordMorphism('a->ab,b->ba')
        sage: wm == loads(dumps(wm))
        True
    """
    def __init__(self, data, codomain=None):
        r"""
        Construction of the morphism.

        EXAMPLES:

        1. If data is a str::

            sage: WordMorphism('a->ab,b->ba')
            WordMorphism: a->ab, b->ba
            sage: WordMorphism('a->ab,b->ba')
            WordMorphism: a->ab, b->ba
            sage: WordMorphism('a->abc,b->bca,c->cab')
            WordMorphism: a->abc, b->bca, c->cab
            sage: WordMorphism('a->abdsf,b->hahdad,c->asdhasd')
            WordMorphism: a->abdsf, b->hahdad, c->asdhasd
            sage: WordMorphism('(->(),)->)(')
            WordMorphism: (->(), )->)(
            sage: WordMorphism('a->53k,b->y5?,$->49i')
            WordMorphism: $->49i, a->53k, b->y5?

        An erasing morphism::

            sage: WordMorphism('a->ab,b->')
            WordMorphism: a->ab, b->

        Use the arrows ('->') correctly::

            sage: WordMorphism('a->ab,b-')
            Traceback (most recent call last):
            ...
            ValueError: The second and third characters must be '->' (not '-')
            sage: WordMorphism('a->ab,b')
            Traceback (most recent call last):
            ...
            ValueError: The second and third characters must be '->' (not '')
            sage: WordMorphism('a->ab,a-]asdfa')
            Traceback (most recent call last):
            ...
            ValueError: The second and third characters must be '->' (not '-]')

        Each letter must be defined only once::

            sage: WordMorphism('a->ab,a->ba')
            Traceback (most recent call last):
            ...
            ValueError: The image of 'a' is defined twice.

        2. From a dictionary::

            sage: WordMorphism({"a":"ab","b":"ba"})
            WordMorphism: a->ab, b->ba
            sage: WordMorphism({2:[4,5,6],3:[1,2,3]})
            WordMorphism: 2->456, 3->123
            sage: WordMorphism({'a':['a',6,'a'],6:[6,6,6,'a']})
            WordMorphism: 6->666a, a->a6a

        The image of a letter can be a set, but the order is not
        preserved::

            sage: WordMorphism({2:[4,5,6],3:set([4,1,8])}) #random results
            WordMorphism: 2->456, 3->814

        If the image of a letter is not iterable, it is considered as a
        letter::

            sage: WordMorphism({0:1, 1:0})
            WordMorphism: 0->1, 1->0
            sage: WordMorphism({0:123, 1:789})
            WordMorphism: 0->123, 1->789
            sage: WordMorphism({2:[4,5,6], 3:123})
            WordMorphism: 2->456, 3->123

        3. From a WordMorphism::

            sage: WordMorphism(WordMorphism('a->ab,b->ba'))
            WordMorphism: a->ab, b->ba

        TESTS::

            sage: WordMorphism(',,,a->ab,,,b->ba,,')
            WordMorphism: a->ab, b->ba
        """
        if isinstance(data, WordMorphism):
            self._domain = data._domain
            self._codomain = data._codomain
            self._morph = data._morph
        else:
            if isinstance(data, str):
                data = self._build_dict(data)
            elif not isinstance(data, dict):
                raise NotImplementedError

            if codomain is None:
                codomain = self._build_codomain(data)

            if not isinstance(codomain,Words_all):
                raise TypeError("the codomain must be a Words domain")
            self._codomain = codomain

            self._morph = {}

            dom_alph = list()
            for (key,val) in data.iteritems():
                dom_alph.append(key)
                if val in codomain.alphabet():
                    self._morph[key] = codomain([val])
                else:
                    self._morph[key] = codomain(val)

            dom_alph.sort()
            self._domain = Words(dom_alph)

    def _build_dict(self, s):
        r"""
        Parse the string input to WordMorphism and build the dictionary
        it represents.

        TESTS::

            sage: wm = WordMorphism('a->ab,b->ba')
            sage: wm._build_dict('a->ab,b->ba') == {'a': 'ab', 'b': 'ba'}
            True
            sage: wm._build_dict('a->ab,a->ba')
            Traceback (most recent call last):
            ...
            ValueError: The image of 'a' is defined twice.
            sage: wm._build_dict('a->ab,b>ba')
            Traceback (most recent call last):
            ...
            ValueError: The second and third characters must be '->' (not '>b')
        """
        tmp_dict = {}
        for fleche in s.split(','):
            if len(fleche) == 0:
                continue

            if len(fleche) < 3 or fleche[1:3] != '->':
                raise ValueError("The second and third characters must be '->' (not '%s')"%fleche[1:3])

            lettre = fleche[0]
            image  = fleche[3:]

            if lettre in tmp_dict:
                raise ValueError("The image of %r is defined twice." %lettre)

            tmp_dict[lettre] = image
        return tmp_dict

    def _build_codomain(self, data):
        r"""
        Returns a Words domain containing all the letter in the keys of
        data (which must be a dictionary).

        TESTS:

        If the image of all the letters are iterable::

            sage: wm = WordMorphism('a->ab,b->ba')
            sage: wm._build_codomain({'a': 'ab', 'b': 'ba'})
            Words over {'a', 'b'}
            sage: wm._build_codomain({'a': 'dcb', 'b': 'a'})
            Words over {'a', 'b', 'c', 'd'}
            sage: wm._build_codomain({2:[4,5,6],3:[1,2,3]})
            Words over {1, 2, 3, 4, 5, 6}
            sage: wm._build_codomain({2:[4,5,6],3:set([4,1,8])})
            Words over {1, 4, 5, 6, 8}

        If the image of a letter is not iterable, it is considered as
        a letter::

            sage: wm._build_codomain({2:[4,5,6],3:123})
            Words over {4, 5, 6, 123}
            sage: wm._build_codomain({0:1, 1:0, 2:2})
            Words over {0, 1, 2}
        """
        codom_alphabet = set()
        for key,val in data.iteritems():
            try:
                it = iter(val)
            except Exception:
                it = [val]
            codom_alphabet.update(it)
        return Words(sorted(codom_alphabet))

    def __eq__(self, other):
        r"""
        Returns ``True`` if ``self`` is equal to ``other``.

        EXAMPLES::

            sage: n = WordMorphism('a->a,b->aa,c->aaa')
            sage: n**3 == n**1
            True
            sage: WordMorphism('b->ba,a->ab') == WordMorphism('a->ab,b->ba')
            True
            sage: WordMorphism('b->ba,a->ab') == WordMorphism({"a":"ab","b":"ba"})
            True
            sage: m = WordMorphism({0:[1,2,3],1:[4,5,6]}); m
            WordMorphism: 0->123, 1->456
            sage: o = WordMorphism('0->123,1->456'); o
            WordMorphism: 0->123, 1->456
            sage: m == o
            False

        TESTS:

        Check that equality depends on the codomain::

            sage: m = WordMorphism('a->a,b->aa,c->aaa')
            sage: n = WordMorphism('a->a,b->aa,c->aaa', codomain=Words('abc'))
            sage: m == n
            False
        """
        if not isinstance(other, WordMorphism):
            return False
        return self._morph == other._morph and self._codomain == other._codomain

    def __ne__(self, other):
        r"""
        Returns whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->baba')
            sage: n = WordMorphism('a->ab,b->baba')
            sage: o = WordMorphism('a->ab,b->bab')
            sage: m != n
            False
            sage: n != o
            True

        This solves :trac:`12475`::

            sage: s = WordMorphism('1->121,2->131,3->4,4->1')
            sage: s == s.reversal()
            True
            sage: s != s.reversal()
            False

        """
        return not self.__eq__(other)

    def __repr__(self):
        r"""
        Returns the string representation of the morphism.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->ba')
            WordMorphism: a->ab, b->ba
            sage: WordMorphism({0:[0,1],1:[1,0]})
            WordMorphism: 0->01, 1->10

        TESTS::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: repr(s)
            'WordMorphism: a->ab, b->ba'
        """
        return "WordMorphism: %s" % str(self)

    def __str__(self):
        r"""
        Returns the morphism in str.

        EXAMPLES::

            sage: print WordMorphism('a->ab,b->ba')
            a->ab, b->ba
            sage: print WordMorphism({0:[0,1],1:[1,0]})
            0->01, 1->10

        The output is sorted to make it unique::

            sage: print WordMorphism('b->ba,a->ab')
            a->ab, b->ba

        The str method is used for string formatting::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: "Here is a map : %s" % s
            'Here is a map : a->ab, b->ba'

        ::

            sage: s = WordMorphism({1:[1,2],2:[1]})
            sage: s.dual_map()
            E_1^*(1->12, 2->1)

        TESTS::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: str(s)
            'a->ab, b->ba'
        """
        L = [str(lettre) + '->' + image.string_rep() for lettre,image in self._morph.iteritems()]
        return ', '.join(sorted(L))

    def __call__(self, w, order=1, datatype='iter'):
        r"""
        Returns the image of ``w`` under self to the given order.

        INPUT:

        -  ``w`` - word or sequence in the domain of self
        -  ``order`` - integer or plus ``Infinity`` (default: 1)
        - ``datatype`` - (default: ``'iter'``) ``'list'``, ``'str'``,
          ``'tuple'``, ``'iter'``. The datatype of the output
          (note that only list, str and tuple allows the word to be
          pickled and saved).

        OUTPUT:

        -  ``word`` - order-th iterated image under self of ``w``

        EXAMPLES:

        The image of a word under a morphism:

        1. The image of a finite word under a morphism::

            sage: tm = WordMorphism ('a->ab,b->ba')
            sage: tm('a')
            word: ab
            sage: tm('aabababb')
            word: ababbaabbaabbaba

        2. The iterated image of a word::

            sage: tm('a', 2)
            word: abba
            sage: tm('aba', 3)
            word: abbabaabbaababbaabbabaab

        3. The infinitely iterated image of a letter::

            sage: tm('a', oo)
            word: abbabaabbaababbabaababbaabbabaabbaababba...

        4. The image of an infinite word::

            sage: t = words.ThueMorseWord()
            sage: n = WordMorphism({0:[0, 1], 1:[1, 0]})
            sage: n(t)
            word: 0110100110010110100101100110100110010110...
            sage: n(t, 3)
            word: 0110100110010110100101100110100110010110...
            sage: n(t)[:1000] == t[:1000]
            True

        The Fibonacci word::

            sage: w = words.FibonacciWord()
            sage: m = WordMorphism({0:'a', 1:'b'})
            sage: m(w)
            word: abaababaabaababaababaabaababaabaababaaba...
            sage: f = words.FibonacciWord('ab')
            sage: f[:1000] == m(w)[:1000]
            True

        ::

            sage: w = words.FibonacciWord("ab")
            sage: m = WordMorphism('a->01,b->101')
            sage: m(w)
            word: 0110101011010110101011010101101011010101...

        The default datatype of the output is an iterable which
        can be saved (for finite word only)::

            sage: m = WordMorphism('a->ab,b->ba')
            sage: w = m('aabb')
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_iter_with_caching'>
            sage: w == loads(dumps(w))
            True
            sage: save(w, filename=os.path.join(SAGE_TMP, 'test.sobj'))

        One may impose the datatype of the resulting word::

            sage: w = m('aaab',datatype='list')
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_list'>
            sage: w = m('aaab',datatype='str')
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: w = m('aaab',datatype='tuple')
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_tuple'>

        To use str datatype for the output word, the domain and codomain
        alphabet must consist of str objects::

            sage: m = WordMorphism({0:[0,1],1:[1,0]})
            sage: w = m([0],4); type(w)
            <class 'sage.combinat.words.word.FiniteWord_iter_with_caching'>
            sage: w = m([0],4,datatype='list'); type(w)
            <class 'sage.combinat.words.word.FiniteWord_list'>
            sage: w = m([0],4,datatype='str')
            Traceback (most recent call last):
            ...
            ValueError: 0 not in alphabet!
            sage: w = m([0],4,datatype='tuple'); type(w)
            <class 'sage.combinat.words.word.FiniteWord_tuple'>

        The word must be in the domain of self::

            sage: tm('0021')
            Traceback (most recent call last):
            ...
            KeyError: '0'

        The order must be a positive integer or plus Infinity::

            sage: tm('a', -1)
            Traceback (most recent call last):
            ...
            TypeError: order (-1) must be a positive integer or plus Infinity
            sage: tm('a', 6.7)
            Traceback (most recent call last):
            ...
            TypeError: order (6.70000000000000) must be a positive integer or plus Infinity

        Only the first letter is considered for infinitely iterated image of
        a word under a morphism::

            sage: tm('aba',oo)
            word: abbabaabbaababbabaababbaabbabaabbaababba...

        The morphism self must be prolongable on the given letter for infinitely
        iterated image::

            sage: m = WordMorphism('a->ba,b->ab')
            sage: m('a', oo)
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on a

        The empty word is fixed by any morphism for all natural
        powers::

            sage: phi = WordMorphism('a->ab,b->a')
            sage: phi(Word())
            word:
            sage: phi(Word(), oo)
            word:
            sage: it = iter([])
            sage: phi(it, oo)
            word:

        TESTS::

            sage: for i in range(6):
            ...     tm('a', i)
            ...
            word: a
            word: ab
            word: abba
            word: abbabaab
            word: abbabaabbaababba
            word: abbabaabbaababbabaababbaabbabaab
            sage: m = WordMorphism('a->,b->')
            sage: m('')
            word:
        """
        if order == 1:
            if isinstance(w, (tuple,str,list)):
                length = 'finite'
            elif isinstance(w, FiniteWord_class):
                #Is it really a good thing to precompute the length?
                length = sum(self._morph[a].length() * b for (a,b) in w.evaluation_dict().iteritems())
            elif hasattr(w, '__iter__'):
                length = Infinity
                datatype = 'iter'
            elif w in self._domain.alphabet():
                return self._morph[w]
            else:
                raise TypeError("Don't know how to handle an input (=%s) that is not iterable or not in the domain alphabet."%w)
            return self.codomain()((x for y in w for x in self._morph[y]), length=length, datatype=datatype)
        elif order is Infinity:
            if isinstance(w, (tuple,str,list,FiniteWord_class)):
                if len(w) == 0:
                    return self.codomain()()
                else:
                    letter = w[0]
            elif hasattr(w, '__iter__'):
                try:
                    letter = next(w)
                except StopIteration:
                    return self.codomain()()
            elif w in self._domain.alphabet():
                letter = w
            else:
                raise TypeError("Don't know how to handle an input (=%s) that is not iterable or not in the domain alphabet."%w)
            return self.fixed_point(letter=letter)
        elif isinstance(order, (int,Integer)) and order > 1:
            return self(self(w, order-1),datatype=datatype)
        elif order == 0:
            return self._domain(w)
        else:
            raise TypeError("order (%s) must be a positive integer or plus Infinity" % order)

    def latex_layout(self, layout=None):
        r"""
        Get or set the actual latex layout (oneliner vs array).

        INPUT:

        - ``layout`` - string (default: ``None``), can take one of the
          following values:

          - ``None`` - Returns the actual latex layout. By default, the
            layout is ``'array'``
          - ``'oneliner'`` - Set the layout to ``'oneliner'``
          - ``'array'`` - Set the layout to ``'array'``

        EXAMPLES::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: s.latex_layout()
            'array'
            sage: s.latex_layout('oneliner')
            sage: s.latex_layout()
            'oneliner'
        """
        if layout is None:
            # return the layout
            if not hasattr(self, '_latex_layout'):
                self._latex_layout = 'array'
            return self._latex_layout
        else:
            # change the layout
            self._latex_layout = layout

    def _latex_(self):
        r"""
        Returns the latex representation of the morphism.

        Use :method:`latex_layout` to change latex layout (oneliner vs
        array). The default is an latex array.

        EXAMPLES::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: s._latex_()
            \begin{array}{l}
            a \mapsto ab\\
            b \mapsto ba
            \end{array}

        Change the latex layout to a one liner::

            sage: s.latex_layout('oneliner')
            sage: s._latex_()
            a \mapsto ab,b \mapsto ba

        TESTS:

        Unknown latex style::

            sage: s.latex_layout('tabular')
            sage: s._latex_()
            Traceback (most recent call last):
            ...
            ValueError: unknown latex_layout(=tabular)

        """
        from sage.misc.latex import LatexExpr
        A = self.domain().alphabet()
        latex_layout = self.latex_layout()
        if latex_layout == 'oneliner':
            L = [r"%s \mapsto %s" % (a, self.image(a)) for a in A]
            return LatexExpr(r','.join(L))
        elif latex_layout == 'array':
            s =  r""
            s += r"\begin{array}{l}" + '\n'
            lines = []
            for a in A:
                lines.append(r"%s \mapsto %s"% (a, self.image(a)))
            s += '\\\\\n'.join(lines)
            s += '\n' + "\end{array}"
            return LatexExpr(s)
        else:
            raise ValueError('unknown latex_layout(=%s)' % latex_layout)

    def __mul__(self, other):
        r"""
        Returns the morphism ``self``\*``other``.

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->ba')
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo*m
            WordMorphism: a->aba, b->aab
            sage: fibo*fibo
            WordMorphism: a->aba, b->ab
            sage: m*fibo
            WordMorphism: a->abba, b->ab

        ::

            sage: n = WordMorphism('a->a,b->aa,c->aaa')
            sage: p1 = n*m
            sage: p1
            WordMorphism: a->aaa, b->aaa
            sage: p1.domain()
            Words over {'a', 'b'}
            sage: p1.codomain()
            Words over {'a'}

        ::

            sage: p2 = m*n
            sage: p2
            WordMorphism: a->ab, b->abab, c->ababab
            sage: p2.domain()
            Words over {'a', 'b', 'c'}
            sage: p2.codomain()
            Words over {'a', 'b'}

        ::

            sage: m = WordMorphism('0->a,1->b')
            sage: n = WordMorphism('a->c,b->e',codomain=Words('abcde'))
            sage: p = n * m
            sage: p.codomain()
            Words over {'a', 'b', 'c', 'd', 'e'}

        TESTS::

            sage: m = WordMorphism('a->b,b->c,c->a')
            sage: WordMorphism('')*m
            Traceback (most recent call last):
            ...
            KeyError: 'b'
            sage: m * WordMorphism('')
            WordMorphism:
        """
        return WordMorphism(dict((key, self(w)) for (key, w) in other._morph.iteritems()), codomain=self.codomain())

    def __pow__(self, exp):
        r"""
        Returns the power of ``self`` with exponent = ``exp``.

        INPUT:

        -  ``exp`` - a positive integer

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->ba')
            sage: m^1
            WordMorphism: a->ab, b->ba
            sage: m^2
            WordMorphism: a->abba, b->baab
            sage: m^3
            WordMorphism: a->abbabaab, b->baababba

        The exponent must be a positive integer::

            sage: m^1.5
            Traceback (most recent call last):
            ...
            ValueError: exponent (1.50000000000000) must be an integer
            sage: m^-2
            Traceback (most recent call last):
            ...
            ValueError: exponent (-2) must be strictly positive

        When ``self`` is not an endomorphism::

            sage: n = WordMorphism('a->ba,b->abc')
            sage: n^2
            Traceback (most recent call last):
            ...
            KeyError: 'c'
        """
        #If exp is not an integer
        if not isinstance(exp, (int,Integer)):
            raise ValueError("exponent (%s) must be an integer" %exp)

        #If exp is negative
        elif exp <= 0:
            raise ValueError("exponent (%s) must be strictly positive" %exp)

        #Base of induction
        elif exp == 1:
            return self

        else:
            nexp = int(exp / 2)
            over = exp % 2
            res = (self * self) ** nexp
            if over == 1:
                res *= self
            return res

    def extend_by(self, other):
        r"""
        Returns ``self`` extended by ``other``.

        Let `\varphi_1:A^*\rightarrow B^*` and `\varphi_2:C^*\rightarrow D^*`
        be two morphisms. A morphism `\mu:(A\cup C)^*\rightarrow (B\cup D)^*`
        corresponds to `\varphi_1` *extended by* `\varphi_2` if
        `\mu(a)=\varphi_1(a)` if `a\in A` and `\mu(a)=\varphi_2(a)` otherwise.

        INPUT:

        -  ``other`` - a WordMorphism.

        OUTPUT:

        WordMorphism

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->ba')
            sage: n = WordMorphism({0:1,1:0,'a':5})
            sage: m.extend_by(n)
            WordMorphism: 0->1, 1->0, a->ab, b->ba
            sage: n.extend_by(m)
            WordMorphism: 0->1, 1->0, a->5, b->ba
            sage: m.extend_by(m)
            WordMorphism: a->ab, b->ba

        TESTS::

            sage: m.extend_by(WordMorphism({})) == m
            True
            sage: m.extend_by(WordMorphism('')) == m
            True

        ::

            sage: m.extend_by(4)
            Traceback (most recent call last):
            ...
            TypeError: other (=4) is not a WordMorphism
        """
        if not isinstance(other, WordMorphism):
            raise TypeError("other (=%s) is not a WordMorphism"%other)

        nv = dict(other._morph)
        for k,v in self._morph.iteritems():
            nv[k] = v
        return WordMorphism(nv)

    def restrict_domain(self, alphabet):
        r"""
        Returns a restriction of ``self`` to the given alphabet.

        INPUT:

        - ``alphabet`` - an iterable

        OUTPUT:

        WordMorphism

        EXAMPLES::

            sage: m = WordMorphism('a->b,b->a')
            sage: m.restrict_domain('a')
            WordMorphism: a->b
            sage: m.restrict_domain('')
            WordMorphism:
            sage: m.restrict_domain('A')
            WordMorphism:
            sage: m.restrict_domain('Aa')
            WordMorphism: a->b

        The input alphabet must be iterable::

            sage: m.restrict_domain(66)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
        """
        return WordMorphism(dict((a, self(a)) for a in alphabet if a in self.domain().alphabet()))

    def _matrix_(self, R=None):
        r"""
        Returns the incidence matrix of the morphism over the specified ring.

        EXAMPLES::

            sage: fibo = WordMorphism('a->ab,b->a')
            sage: tm = WordMorphism('a->ab,b->ba')
            sage: Mfibo = matrix(fibo); Mfibo     # indirect doctest
            [1 1]
            [1 0]
            sage: Mtm = matrix(tm); Mtm
            [1 1]
            [1 1]
            sage: Mtm * Mfibo == matrix(tm*fibo)   # indirect doctest
            True
            sage: Mfibo * Mtm == matrix(fibo*tm)   # indirect doctest
            True
            sage: Mfibo.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: p = Mfibo.charpoly(); p
            x^2 - x - 1
            sage: p.roots(ring=RR, multiplicities=False)
            [-0.618033988749895, 1.61803398874989]
        """
        if R is None:
            return self.incidence_matrix()
        else:
            return self.incidence_matrix().change_ring(R)

    def incidence_matrix(self):
        r"""
        Returns the incidence matrix of the morphism. The order of the rows
        and column are given by the order defined on the alphabet of the
        domain and the codomain.

        The matrix returned is over the integers.  If a different ring is
        desired, use either the ``change_ring`` function or the ``matrix``
        function.

        EXAMPLES::

            sage: m = WordMorphism('a->abc,b->a,c->c')
            sage: m.incidence_matrix()
            [1 1 0]
            [1 0 0]
            [1 0 1]
            sage: m = WordMorphism('a->abc,b->a,c->c,d->abbccccabca,e->abc')
            sage: m.incidence_matrix()
            [1 1 0 3 1]
            [1 0 0 3 1]
            [1 0 1 5 1]
        """
        L = []
        domain_alphabet = self.domain().alphabet()
        codomain_alphabet = self.codomain().alphabet()
        for b in domain_alphabet:
            w = self._morph[b]
            ev_dict = w.evaluation_dict()
            L.append([ev_dict.get(a,0) for a in codomain_alphabet])
        M = Matrix(IntegerRing(), L).transpose()
        return M

    def domain(self):
        r"""
        Returns domain of ``self``.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').domain()
            Words over {'a', 'b'}
            sage: WordMorphism('b->ba,a->ab').domain()
            Words over {'a', 'b'}
            sage: WordMorphism('6->ab,y->5,0->asd').domain()
            Words over {'0', '6', 'y'}
        """
        return self._domain

    def codomain(self):
        r"""
        Returns the codomain of ``self``.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').codomain()
            Words over {'a', 'b'}
            sage: WordMorphism('6->ab,y->5,0->asd').codomain()
            Words over {'5', 'a', 'b', 'd', 's'}
        """
        return self._codomain

    def is_endomorphism(self):
        r"""
        Returns ``True`` if the codomain is a subset of the domain.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').is_endomorphism()
            True
            sage: WordMorphism('6->ab,y->5,0->asd').is_endomorphism()
            False
            sage: WordMorphism('a->a,b->aa,c->aaa').is_endomorphism()
            False
            sage: Wabc = Words('abc')
            sage: m = WordMorphism('a->a,b->aa,c->aaa',codomain = Wabc)
            sage: m.is_endomorphism()
            True

        We check that :trac:`8674` is fixed::

            sage: P = WordPaths('abcd')
            sage: m = WordMorphism('a->adab,b->ab,c->cbcd,d->cd', codomain=P)
            sage: m.is_endomorphism()
            True
        """
        return self.codomain() == self.domain()

    def image(self, letter):
        r"""
        Return the image of a letter

        INPUT:

        - ``letter`` - a letter in the domain alphabet

        OUTPUT:

        word

        ..NOTE::

            The letter is assumed to be in the domain alphabet
            (no check done). Hence, this method is faster
            than the ``__call__`` method suitable for words input.

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->ac,c->a')
            sage: m.image('b')
            word: ac

        ::

            sage: s = WordMorphism({('a', 1):[('a', 1), ('a', 2)], ('a', 2):[('a', 1)]})
            sage: s.image(('a',1))
            word: ('a', 1),('a', 2)

        ::

            sage: s = WordMorphism({0:[1,2], 'a':(2,3,4), ():[9,8,7]})
            sage: s.image(0)
            word: 12
            sage: s.image('a')
            word: 234
            sage: s.image(())
            word: 987
        """
        return self._morph[letter]

    def images(self):
        r"""
        Returns the list of all the images of the letters of the alphabet
        under ``self``.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').images()
            [word: ab, word: a]
            sage: WordMorphism('6->ab,y->5,0->asd').images()
            [word: 5, word: asd, word: ab]
        """
        return self._morph.values()

    def reversal(self):
        r"""
        Returns the reversal of ``self``.

        EXAMPLES::

            sage: WordMorphism('6->ab,y->5,0->asd').reversal()
            WordMorphism: 0->dsa, 6->ba, y->5
            sage: WordMorphism('a->ab,b->a').reversal()
            WordMorphism: a->ba, b->a
        """
        return WordMorphism(dict((key, w.reversal()) for (key, w) in self._morph.iteritems()),codomain=self._codomain)

    def is_empty(self):
        r"""
        Returns ``True`` if the cardinality of the domain is zero and
        ``False`` otherwise.

        EXAMPLES::

            sage: WordMorphism('').is_empty()
            True
            sage: WordMorphism('a->a').is_empty()
            False
        """
        return len(self._morph) == 0

    def is_erasing(self):
        r"""
        Returns ``True`` if ``self`` is an erasing morphism, i.e. the image of a
        letter is the empty word.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').is_erasing()
            False
            sage: WordMorphism('6->ab,y->5,0->asd').is_erasing()
            False
            sage: WordMorphism('6->ab,y->5,0->asd,7->').is_erasing()
            True
            sage: WordMorphism('').is_erasing()
            False
        """
        for image in self.images():
            if image.is_empty():
                return True
        return False

    def is_identity(self):
        r"""
        Returns ``True`` if ``self`` is the identity morphism.

        EXAMPLES::

            sage: m = WordMorphism('a->a,b->b,c->c,d->e')
            sage: m.is_identity()
            False
            sage: WordMorphism('a->a,b->b,c->c').is_identity()
            True
            sage: WordMorphism('a->a,b->b,c->cb').is_identity()
            False
            sage: m = WordMorphism('a->b,b->c,c->a')
            sage: (m^2).is_identity()
            False
            sage: (m^3).is_identity()
            True
            sage: (m^4).is_identity()
            False
            sage: WordMorphism('').is_identity()
            True
            sage: WordMorphism({0:[0],1:[1]}).is_identity()
            True

        We check that #8618 is fixed::

            sage: t = WordMorphism({'a1':['a2'], 'a2':['a1']})
            sage: (t*t).is_identity()
            True
        """
        if self.domain() != self.codomain():
            return False

        for letter in self.domain().alphabet():
            img = self.image(letter)
            if img.length() != 1:
                return False
            elif img[0] != letter:
                return False
        return True

    def partition_of_domain_alphabet(self):
        r"""
        Returns a partition of the domain alphabet.

        Let `\varphi:\Sigma^*\rightarrow\Sigma^*` be an involution. There
        exists a triple of sets `(A, B, C)` such that

         -  `A \cup B \cup C =\Sigma`;
         -  `A`, `B` and `C` are mutually disjoint and
         -  `\varphi(A)= B`, `\varphi(B)= A`, `\varphi(C)= C`.

        These sets are not unique.

        INPUT:

        - ``self`` - An involution.

        OUTPUT:

        A tuple of three sets

        EXAMPLES::

            sage: m = WordMorphism('a->b,b->a')
            sage: m.partition_of_domain_alphabet() #random ordering
            ({'a'}, {'b'}, {})
            sage: m = WordMorphism('a->b,b->a,c->c')
            sage: m.partition_of_domain_alphabet() #random ordering
            ({'a'}, {'b'}, {'c'})
            sage: m = WordMorphism('a->a,b->b,c->c')
            sage: m.partition_of_domain_alphabet() #random ordering
            ({}, {}, {'a', 'c', 'b'})
            sage: m = WordMorphism('A->T,T->A,C->G,G->C')
            sage: m.partition_of_domain_alphabet() #random ordering
            ({'A', 'C'}, {'T', 'G'}, {})
            sage: I = WordMorphism({0:oo,oo:0,1:-1,-1:1,2:-2,-2:2,3:-3,-3:3})
            sage: I.partition_of_domain_alphabet() #random ordering
            ({0, -1, -3, -2}, {1, 2, 3, +Infinity}, {})

        TESTS::

            sage: m = WordMorphism('a->b,b->a,c->a')
            sage: m.partition_of_domain_alphabet()
            Traceback (most recent call last):
            ...
            TypeError: self (=a->b, b->a, c->a) is not an endomorphism
        """
        if not self.is_involution():
            raise TypeError("self is not an involution")

        A = set(); B = set(); C = set()
        for a in self.domain().alphabet():
            if a == self(a)[0]:
                C.add(a)
            elif not (a in A or a in B):
                A.add(a)
                B.add(self(a)[0])

        return Set(A), Set(B), Set(C)

    def is_involution(self):
        r"""
        Returns ``True`` if ``self`` is an involution, i.e. its square
        is the identity.

        INPUT:

        - ``self`` - an endomorphism

        EXAMPLES::

            sage: WordMorphism('a->b,b->a').is_involution()
            True
            sage: WordMorphism('a->b,b->ba').is_involution()
            False
            sage: WordMorphism({0:[1],1:[0]}).is_involution()
            True

        TESTS::

            sage: WordMorphism('').is_involution()
            True
            sage: WordMorphism({0:1,1:0,2:3}).is_involution()
            Traceback (most recent call last):
            ...
            TypeError: self (=0->1, 1->0, 2->3) is not an endomorphism
        """
        if not self.is_endomorphism():
            raise TypeError("self (=%s) is not an endomorphism"%self)

        return (self*self).is_identity()

    def pisot_eigenvector_right(self):
        r"""
        Returns the right eigenvector of the incidence matrix associated
        to the largest eigenvalue (in absolute value).

        Unicity of the result is guaranteed when the multiplicity of the
        largest eigenvalue is one, for example when self is a Pisot
        irreductible substitution.

        A substitution is Pisot irreducible if the characteristic
        polynomial of its incidence matrix is irreducible over `\QQ` and
        has all roots, except one, of modulus strictly smaller than 1.

        INPUT:

        - ``self`` - a Pisot irreducible substitution.

        EXAMPLES::

            sage: m = WordMorphism('a->aaaabbc,b->aaabbc,c->aabc')
            sage: matrix(m)
            [4 3 2]
            [2 2 1]
            [1 1 1]
            sage: m.pisot_eigenvector_right()
            (1, 0.5436890126920763?, 0.2955977425220848?)
        """
        eig = self.incidence_matrix().eigenvectors_right()
        return max(eig, key=lambda x:abs(x[0]))[1][0]

    def pisot_eigenvector_left(self):
        r"""
        Returns the left eigenvector of the incidence matrix associated
        to the largest eigenvalue (in absolute value).

        Unicity of the result is guaranteed when the multiplicity of the
        largest eigenvalue is one, for example when self is a Pisot
        irreductible substitution.

        A substitution is Pisot irreducible if the characteristic
        polynomial of its incidence matrix is irreducible over `\QQ` and
        has all roots, except one, of modulus strictly smaller than 1.

        INPUT:

        - ``self`` - a Pisot irreducible substitution.

        EXAMPLES::

            sage: m = WordMorphism('a->aaaabbc,b->aaabbc,c->aabc')
            sage: matrix(m)
            [4 3 2]
            [2 2 1]
            [1 1 1]
            sage: m.pisot_eigenvector_left()
            (1, 0.8392867552141611?, 0.5436890126920763?)
        """
        eig = self.incidence_matrix().eigenvectors_left()
        return max(eig, key=lambda x:abs(x[0]))[1][0]

    def _check_primitive(self):
        r"""
        Returns ``True`` if all the letters of the domain appear in all the
        images of letters of the domain.

        INPUT:

        - ``self`` - the codomain must be an instance of Words

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->ba')
            sage: m._check_primitive()
            True
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo._check_primitive()
            False
            sage: WordMorphism({2:[4,5,6],3:[4,1,8]})
            WordMorphism: 2->456, 3->418
            sage: WordMorphism({2:[4,5,6],3:[4,1,8]})._check_primitive()
            False

        """
        if not isinstance(self.codomain(),Words_all):
            raise TypeError("codomain of self(=%s) must be an instance of Words"%self)

        dom_alphabet = set(self.domain().alphabet())

        for image in self.images():
            if not dom_alphabet <= set(image):
                return False
        else:
            return True

    def is_primitive(self):
        r"""
        Returns ``True`` if ``self`` is primitive.

        A morphism `\varphi` is *primitive* if there exists
        an positive integer `k` such that for all `\alpha\in\Sigma`,
        `\varphi^k(\alpha)` contains all the letters of `\Sigma`.

        INPUT:

        - ``self`` - an endomorphism

        ALGORITHM:

            Exercices 8.7.8, p.281 in [1] :
            (c) Let `y(M)` be the least integer `e` such that `M^e` has all
            positive entries. Prove that, for all primitive matrices `M`,
            we have `y(M) \leq (d-1)^2 + 1`.
            (d) Prove that the bound `y(M)\leq (d-1)^2+1` is best possible.

        EXAMPLES::

            sage: tm = WordMorphism('a->ab,b->ba')
            sage: tm.is_primitive()
            True
            sage: fibo = WordMorphism('a->ab,b->a');
            sage: fibo.is_primitive()
            True
            sage: m = WordMorphism('a->bb,b->aa')
            sage: m.is_primitive()
            False
            sage: f = WordMorphism({0:[1],1:[0]})
            sage: f.is_primitive()
            False

        ::

            sage: s = WordMorphism('a->b,b->c,c->ab')
            sage: s.is_primitive()
            True
            sage: s = WordMorphism('a->b,b->c,c->d,d->e,e->f,f->g,g->h,h->ab')
            sage: s.is_primitive()
            True

        TESTS::

            sage: m = WordMorphism('a->bb,b->aac')
            sage: m.is_primitive()
            Traceback (most recent call last):
            ...
            TypeError: self (=a->bb, b->aac) is not an endomorphism
            sage: m = WordMorphism('a->,b->',codomain=Words('ab'))
            sage: m.is_primitive()
            False
            sage: m = WordMorphism('a->,b->')
            sage: m.is_primitive()
            Traceback (most recent call last):
            ...
            TypeError: self (=a->, b->) is not an endomorphism

        REFERENCES:

        - [1] Jean-Paul Allouche and Jeffrey Shallit, Automatic Sequences:
          Theory, Applications, Generalizations, Cambridge University Press,
          2003.
        """
        if not self.is_endomorphism():
            raise TypeError("self (=%s) is not an endomorphism"%self)
        m = self.incidence_matrix()
        power = m
        order = 1
        dim = self.domain().size_of_alphabet()
        max_order = (dim-1)**2 + 1
        while True:
            l = power.list()
            if len(l) == 0:
                return False
            try:
                l.index(0)
            except ValueError:
                return True
            if order > max_order:
                return False
            power *= power
            order += order

    def is_prolongable(self, letter):
        r"""
        Returns ``True`` if ``self`` is prolongable on ``letter``.

        A morphism `\varphi` is prolongable on a letter `a`
        if `a` is a prefix of `\varphi(a)`.

        INPUT:

        - ``self`` - its codomain must be an instance of Words
        - ``letter`` - a letter in the domain alphabet

        OUTPUT:

        Boolean

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').is_prolongable(letter='a')
            True
            sage: WordMorphism('a->ab,b->a').is_prolongable(letter='b')
            False
            sage: WordMorphism('a->ba,b->ab').is_prolongable(letter='b')
            False
            sage: (WordMorphism('a->ba,b->ab')^2).is_prolongable(letter='b')
            True
            sage: WordMorphism('a->ba,b->').is_prolongable(letter='b')
            False
            sage: WordMorphism('a->bb,b->aac').is_prolongable(letter='a')
            False

        We check that #8595 is fixed::

            sage: s = WordMorphism({('a', 1) : [('a', 1), ('a', 2)], ('a', 2) : [('a', 1)]})
            sage: s.is_prolongable(('a',1))
            True

        TESTS::

            sage: WordMorphism('a->ab,b->b,c->ba').is_prolongable(letter='d')
            Traceback (most recent call last):
            ...
            TypeError: letter (=d) is not in the domain alphabet (={'a', 'b', 'c'})

        ::

            sage: n0, n1 = matrix(2,[1,1,1,0]), matrix(2,[2,1,1,0])
            sage: n = {'a':n0, 'b':n1}
            sage: WordMorphism(n).is_prolongable(letter='a') #todo: not implemented
            Traceback (most recent call last):
            ...
            TypeError: codomain of self must be an instance of Words
        """
        if not isinstance(self.codomain(), Words_all):
            raise TypeError("codomain of self must be an instance of Words")

        if letter not in self.domain().alphabet():
            raise TypeError("letter (=%s) is not in the domain alphabet (=%s)"\
                                %(letter, self.domain().alphabet()))
        image = self.image(letter)
        return not image.is_empty() and letter == image[0]

    def is_uniform(self, k=None):
        r"""
        Returns True if self is a `k`-uniform morphism.

        Let `k` be a positive integer. A morphism `\phi` is called `k`-uniform
        if for every letter `\alpha`, we have `|\phi(\alpha)| = k`. In other
        words, all images have length `k`. A morphism is called uniform if it
        is `k`-uniform for some positive integer `k`.

        INPUT:

        - ``k`` - a positive integer or None. If set to a positive integer,
          then the function return True if self is `k`-uniform. If set to
          None, then the function return True if self is uniform.

        EXAMPLES::

            sage: phi = WordMorphism('a->ab,b->a')
            sage: phi.is_uniform()
            False
            sage: phi.is_uniform(k=1)
            False
            sage: tau = WordMorphism('a->ab,b->ba')
            sage: tau.is_uniform()
            True
            sage: tau.is_uniform(k=1)
            False
            sage: tau.is_uniform(k=2)
            True
        """
        if k is None:
            return len(set(w.length() for w in self.images())) == 1
        else:
            return all(w.length() == k for w in self.images())

    def _fixed_point_iterator(self, letter):
        r"""
        Returns an iterator of the letters of the fixed point of ``self``
        starting with ``letter``.

        If w is the iterated word, then this iterator: outputs the elements
        of morphism[ w[i] ], appends morphism[ w[i+1] ] to w, increments i.

        INPUT:

        - ``self`` - an endomorphism, must be prolongable on
           letter

        - ``letter`` - a letter in the domain of ``self``

        OUTPUT:

        - iterator of the fixed point

        EXAMPLES::

            sage: m = WordMorphism('a->abc,b->,c->')
            sage: list(m._fixed_point_iterator('a'))
            ['a', 'b', 'c']

        The morphism must be prolongable on the letter or the iterator will
        be empty::

            sage: list(m._fixed_point_iterator('b'))
            []

        The morphism must be an endomorphism::

            sage: m = WordMorphism('a->ac,b->aac')
            sage: list(m._fixed_point_iterator('a'))
            Traceback (most recent call last):
            ...
            KeyError: 'c'

        We check that #8595 is fixed::

            sage: s = WordMorphism({('a', 1):[('a', 1), ('a', 2)], ('a', 2):[('a', 1)]})
            sage: it = s._fixed_point_iterator(('a',1))
            sage: next(it)
            ('a', 1)

        This shows that ticket :trac:`13668` has been resolved::

            sage: s = WordMorphism({1:[1,2],2:[2,3],3:[4],4:[5],5:[6],6:[7],7:[8],8:[9],9:[10],10:[1]})
            sage: s7 = s^7
            sage: s7r = s7.reversal()
            sage: s7r10 = s7r^10
            sage: s7r10.fixed_point(1)
            word: 1,10,9,8,7,6,5,4,3,2,10,9,8,7,6,5,4,3,2,9,8,7,6,5,...
        """
        w = iter(self.image(letter))
        while True:
            for a in self.image(next(w)):
                yield a
            else:
                next_w = next(w)
                w = itertools.chain([next_w], w, self.image(next_w))

    def fixed_point(self, letter):
        r"""
        Returns the fixed point of ``self`` beginning by the given ``letter``.

        A fixed point of morphism `\varphi` is a word `w` such that
        `\varphi(w) = w`.

        INPUT:

        -  ``self`` - an endomorphism, must be prolongable on ``letter``

        -  ``letter`` - in the domain of ``self``, the first letter
           of the fixed point.

        OUTPUT:

        - ``word`` - the fixed point of ``self`` beginning with ``letter``.

        EXAMPLES:

        1. Infinite fixed point::

            sage: WordMorphism('a->ab,b->ba').fixed_point(letter='a')
            word: abbabaabbaababbabaababbaabbabaabbaababba...
            sage: WordMorphism('a->ab,b->a').fixed_point(letter='a')
            word: abaababaabaababaababaabaababaabaababaaba...
            sage: W = Words('abc')
            sage: WordMorphism('a->ab,b->b,c->ba', codomain=W).fixed_point(letter='a')
            word: abbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb...

        2. Infinite fixed point of an erasing morphism::

            sage: WordMorphism('a->ab,b->,c->ba', codomain=W).fixed_point(letter='a')
            word: ab

        3. Finite fixed point::

            sage: WordMorphism('a->ab,b->b,c->ba', codomain=W).fixed_point(letter='b')
            word: b

        4. Finite fixed point of an erasing morphism::

            sage: m = WordMorphism('a->abc,b->,c->')
            sage: fp = m.fixed_point('a'); fp
            word: abc
            sage: m = WordMorphism('a->ba,b->')
            sage: m('ba')
            word: ba
            sage: m.fixed_point('a') #todo: not implemented
            word: ba

        5. Fixed point of a power of a morphism::

            sage: m = WordMorphism('a->ba,b->ab')
            sage: (m^2).fixed_point(letter='a')
            word: abbabaabbaababbabaababbaabbabaabbaababba...

        TESTS::

            sage: WordMorphism('a->ab,b->,c->ba', codomain=W).fixed_point(letter='b')
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on b
            sage: WordMorphism('a->ab,b->,c->ba', codomain=W).fixed_point(letter='c')
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on c
            sage: WordMorphism('a->ab,b->,c->ba', codomain=W).fixed_point(letter='d')
            Traceback (most recent call last):
            ...
            TypeError: letter (=d) is not in the domain alphabet (={'a', 'b', 'c'})
            sage: WordMorphism('a->aa,b->aac').fixed_point(letter='a')
            Traceback (most recent call last):
            ...
            TypeError: self (=a->aa, b->aac) is not an endomorphism
        """
        if not self.is_endomorphism():
            raise TypeError("self (=%s) is not an endomorphism"%self)

        if not self.is_prolongable(letter=letter):
            raise TypeError("self must be prolongable on %s"%letter)

        image = self.image(letter)

        if image.length() == 1:
            return image

        # Construct the word.
        w = self.codomain()(self._fixed_point_iterator(letter), datatype='iter')
        return w

    def fixed_points(self):
        r"""
        Returns the list of all fixed points of ``self``.

        EXAMPLES::

            sage: f = WordMorphism('a->ab,b->ba')
            sage: for w in f.fixed_points(): print w
            abbabaabbaababbabaababbaabbabaabbaababba...
            baababbaabbabaababbabaabbaababbaabbabaab...

            sage: f = WordMorphism('a->ab,b->c,c->a')
            sage: for w in f.fixed_points(): print w
            abcaababcabcaabcaababcaababcabcaababcabc...

            sage: f = WordMorphism('a->ab,b->cab,c->bcc')
            sage: for w in f.fixed_points(): print w
            abcabbccabcabcabbccbccabcabbccabcabbccab...

        This shows that ticket :trac:`13668` has been resolved::

            sage: d = {1:[1,2],2:[2,3],3:[4],4:[5],5:[6],6:[7],7:[8],8:[9],9:[10],10:[1]}
            sage: s = WordMorphism(d)
            sage: s7 = s^7
            sage: s7.fixed_points()
            [word: 12232342..., word: 2,3,4,5,6,7,8...]
            sage: s7r = s7.reversal()
            sage: s7r.fixed_points()
            []

        This shows that ticket :trac:`13668` has been resolved::

            sage: s = "1->321331332133133,2->133321331332133133,3->2133133133321331332133133"
            sage: s = WordMorphism(s)
            sage: (s^2).fixed_points()
            []

        """
        L = []
        for letter in self.domain().alphabet():
            if self.is_prolongable(letter=letter):
                L.append(self.fixed_point(letter=letter))
        return L

    def periodic_point(self, letter):
        r"""
        Return the periodic point of self that starts with ``letter``.

        EXAMPLES::

            sage: f = WordMorphism('a->bab,b->aba')
            sage: f.periodic_point('a')
            word: abababababababababababababababababababab...
            sage: f.fixed_point('a')
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on a
        """
        if self.is_erasing():
            raise NotImplementedError("self should be non erasing")

        cycle = [letter]
        a = self(letter)[0]
        while a not in cycle:
            cycle.append(a)
            a = self(a)[0]
        if a != letter:
            raise ValueError("there is no periodic point starting with letter (=%s)"%letter)
        return (self**len(cycle)).fixed_point(letter)

    def periodic_points(self):
        r"""
        Return the periodic points of ``f`` as a list of tuples where each tuple is
        a periodic orbit of ``f``.

        EXAMPLES::

            sage: f = WordMorphism('a->aba,b->baa')
            sage: for p in f.periodic_points():
            ...      print len(p), ',', p[0]
            1 , ababaaababaaabaabaababaaababaaabaabaabab...
            1 , baaabaabaababaaabaababaaabaababaaababaaa...

            sage: f = WordMorphism('a->bab,b->aa')
            sage: for p in f.periodic_points():
            ...       print len(p), ',', p[0]
            2 , aababaaaababaababbabaababaababbabaababaa...
            sage: f.fixed_points()
            []

        This shows that ticket :trac:`13668` has been resolved::

            sage: d = {1:[1,2],2:[2,3],3:[4],4:[5],5:[6],6:[7],7:[8],8:[9],9:[10],10:[1]}
            sage: s = WordMorphism(d)
            sage: s7 = s^7
            sage: s7r = s7.reversal()
            sage: s7r10 = s7r^10
            sage: for p in s7r10.periodic_points(): p
            [word: 1,10,9,8,7,6,5,4,3,2,10,9,8,7,6,5,4,3,2,...]
            [word: 2,1,1,10,9,8,7,6,5,4,3,2,1,10,9,8,7,6,5,...]
            [word: 3,2,2,1,2,1,1,10,9,8,7,6,5,4,3,2,2,1,1,1...]
            [word: 4,3,2,3,2,2,1,3,2,2,1,2,1,1,10,9,8,7,6,5...]
            [word: 5,4,3,2,4,3,2,3,2,2,1,4,3,2,3,2,2,1,3,2,...]
            [word: 6543254324323221543243232214323221322121...]
            [word: 7654326543254324323221654325432432322154...]
            [word: 8765432765432654325432432322176543265432...]
            [word: 9876543287654327654326543254324323221876...]
            [word: 10,9,8,7,6,5,4,3,2,9,8,7,6,5,4,3,2,8,7,6...]
        """
        assert self.is_endomorphism(), "f should be an endomorphism"

        if self.is_erasing():
            raise NotImplementedError("f should be non erasing")

        A = self.domain().alphabet()
        d = dict((letter,self(letter)[0]) for letter in A)

        res = []
        for cycle in get_cycles(CallableDict(d),A):
            g = self**len(cycle)
            res.append([g.fixed_point(letter) for letter in cycle])

        return res

    def conjugate(self, pos):
        r"""
        Returns the morphism where the image of the letter by ``self``
        is conjugated of parameter ``pos``.

        INPUT:

        - ``pos`` - integer

        EXAMPLES::

            sage: m = WordMorphism('a->abcde')
            sage: m.conjugate(0) == m
            True
            sage: m.conjugate(1)
            WordMorphism: a->bcdea
            sage: m.conjugate(3)
            WordMorphism: a->deabc
            sage: WordMorphism('').conjugate(4)
            WordMorphism:
            sage: m = WordMorphism('a->abcde,b->xyz')
            sage: m.conjugate(2)
            WordMorphism: a->cdeab, b->zxy
        """
        return WordMorphism(dict((key, w.conjugate(pos)) for (key, w) in self._morph.iteritems()))

    def has_left_conjugate(self):
        r"""
        Returns ``True`` if all the non empty images of ``self`` begins with
        the same letter.

        EXAMPLES::

            sage: m = WordMorphism('a->abcde,b->xyz')
            sage: m.has_left_conjugate()
            False
            sage: WordMorphism('b->xyz').has_left_conjugate()
            True
            sage: WordMorphism('').has_left_conjugate()
            True
            sage: WordMorphism('a->,b->xyz').has_left_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb').has_left_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb,c->').has_left_conjugate()
            True
        """
        I = itertools.ifilterfalse(FiniteWord_class.is_empty, self.images())

        try:
            letter = next(I)[0]
        except StopIteration:
            return True

        #Compare the first letter of all the non empty images
        for image in I:
            if image[0] != letter:
                return False

        return True

    def has_right_conjugate(self):
        r"""
        Returns ``True`` if all the non empty images of ``self`` ends with the
        same letter.

        EXAMPLES::

            sage: m = WordMorphism('a->abcde,b->xyz')
            sage: m.has_right_conjugate()
            False
            sage: WordMorphism('b->xyz').has_right_conjugate()
            True
            sage: WordMorphism('').has_right_conjugate()
            True
            sage: WordMorphism('a->,b->xyz').has_right_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb').has_right_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb,c->').has_right_conjugate()
            True
        """
        return self.reversal().has_left_conjugate()

    def list_of_conjugates(self):
        r"""
        Returns the list of all the conjugate morphisms of ``self``.

        DEFINITION:

        Recall from Lothaire [1] (Section 2.3.4)
        that `\varphi` is *right conjugate* of `\varphi'`,
        noted `\varphi\triangleleft\varphi'`, if there exists
        `u \in \Sigma^*` such that

        .. math::

            \varphi(\alpha)u = u\varphi'(\alpha),

        for all `\alpha \in \Sigma`, or equivalently that
        `\varphi(x)u = u\varphi'(x)`, for all words `x \in \Sigma^*`.
        Clearly, this relation is not
        symmetric so that we say that two morphisms `\varphi` and
        `\varphi'` are *conjugate*, noted
        `\varphi\bowtie\varphi'`, if
        `\varphi\triangleleft\varphi'` or
        `\varphi'\triangleleft\varphi`. It is easy to see that
        conjugacy of morphisms is an equivalence relation.

        REFERENCES:

        - [1] M. Lothaire, Algebraic Combinatorics on words, Cambridge
          University Press, 2002.

        EXAMPLES::

            sage: m = WordMorphism('a->abbab,b->abb')
            sage: m.list_of_conjugates()
            [WordMorphism: a->babba, b->bab,
            WordMorphism: a->abbab, b->abb,
            WordMorphism: a->bbaba, b->bba,
            WordMorphism: a->babab, b->bab,
            WordMorphism: a->ababb, b->abb,
            WordMorphism: a->babba, b->bba,
            WordMorphism: a->abbab, b->bab]
            sage: m = WordMorphism('a->aaa,b->aa')
            sage: m.list_of_conjugates()
            [WordMorphism: a->aaa, b->aa]
            sage: WordMorphism('').list_of_conjugates()
            [WordMorphism: ]
            sage: m = WordMorphism('a->aba,b->aba')
            sage: m.list_of_conjugates()
            [WordMorphism: a->baa, b->baa,
            WordMorphism: a->aab, b->aab,
            WordMorphism: a->aba, b->aba]
            sage: m = WordMorphism('a->abb,b->abbab,c->')
            sage: m.list_of_conjugates()
            [WordMorphism: a->bab, b->babba, c->,
            WordMorphism: a->abb, b->abbab, c->,
            WordMorphism: a->bba, b->bbaba, c->,
            WordMorphism: a->bab, b->babab, c->,
            WordMorphism: a->abb, b->ababb, c->,
            WordMorphism: a->bba, b->babba, c->,
            WordMorphism: a->bab, b->abbab, c->]
        """
        if self.is_empty():
            return [self]

        #Construire la liste c des morphismes conjugues
        c = []
        m = self
        c.append(m)
        while(m.has_left_conjugate()):
            m = m.conjugate(1)
            if m == self:
                break
            c.append(m)
        m = self
        while(m.has_right_conjugate()):
            m = m.conjugate(-1)
            if m == self:
                break
            c.insert(0, m)

        #Construire la liste d des morphismes distincts
        d = []
        for m in c:
            if m not in d:
                d.append(m)
        return d

    def is_in_classP(self, f=None):
        r"""
        Returns ``True`` if ``self`` is in class `P` (or `f`-`P`).

        DEFINITION : Let `A` be an alphabet. We say that a
        primitive substitution `S` is in the *class P* if there
        exists a palindrome `p` and for each `b\in A` a
        palindrome `q_b` such that `S(b)=pq_b` for all
        `b\in A`. [1]

        Let `f` be an involution on `A`. "We say that a morphism
        `\varphi` is in class `f`-`P` if there exists an
        `f`-palindrome `p` and for each `\alpha \in A`
        there exists an `f`-palindrome `q_\alpha` such
        that `\varphi(\alpha)=pq_\alpha`. [2]

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of ``self``.
           It must be callable on letters as well as words (e.g. WordMorphism).

        REFERENCES:

        - [1] Hof, A., O. Knill et B. Simon, Singular continuous
          spectrum for palindromic Schrödinger operators,
          Commun. Math. Phys.  174 (1995) 149-159.

        - [2] Labbe, Sebastien. Proprietes combinatoires des
          `f`-palindromes, Memoire de maitrise en Mathematiques,
          Montreal, UQAM, 2008, 109 pages.

        EXAMPLES::

            sage: WordMorphism('a->bbaba,b->bba').is_in_classP()
            True
            sage: tm = WordMorphism('a->ab,b->ba')
            sage: tm.is_in_classP()
            False
            sage: f = WordMorphism('a->b,b->a')
            sage: tm.is_in_classP(f=f)
            True
            sage: (tm^2).is_in_classP()
            True
            sage: (tm^2).is_in_classP(f=f)
            False
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo.is_in_classP()
            True
            sage: fibo.is_in_classP(f=f)
            False
            sage: (fibo^2).is_in_classP()
            False
            sage: f = WordMorphism('a->b,b->a,c->c')
            sage: WordMorphism('a->acbcc,b->acbab,c->acbba').is_in_classP(f)
            True
        """
        if self.is_empty():
            return True

        #Compute the longest common prefix of all the images of letters
        images = self.images()
        lcp = images[0]
        for image in images:
            lcp = lcp.longest_common_prefix(image)

        #Find a common palindrome prefix
        for i in range(lcp.length()+1):
            if lcp[:i].is_palindrome(f=f):

                #If all the suffixes are palindromes,
                for image in images:
                    if not image[i:].is_palindrome(f=f):
                        break
                else:
                    return True

        return False

    def has_conjugate_in_classP(self, f=None):
        r"""
        Returns ``True`` if ``self`` has a conjugate in class `f`-`P`.

        DEFINITION : Let `A` be an alphabet. We say that a
        primitive substitution `S` is in the *class P* if there
        exists a palindrome `p` and for each `b\in A` a
        palindrome `q_b` such that `S(b)=pq_b` for all
        `b\in A`. [1]

        Let `f` be an involution on `A`. We say that a morphism
        `\varphi` is in class `f`-`P` if there exists an
        `f`-palindrome `p` and for each `\alpha \in A`
        there exists an `f`-palindrome `q_\alpha` such
        that `\varphi(\alpha)=pq_\alpha`. [2]

        INPUT:

        -  ``f`` - involution (default: None) on the alphabet of ``self``.
           It must be callable on letters as well as words (e.g. WordMorphism).

        REFERENCES:

        - [1] Hof, A., O. Knill et B. Simon, Singular continuous
          spectrum for palindromic Schrödinger operators,
          Commun. Math. Phys.  174 (1995) 149-159.

        - [2] Labbe, Sebastien. Proprietes combinatoires des
          `f`-palindromes, Memoire de maitrise en Mathematiques,
          Montreal, UQAM, 2008, 109 pages.

        EXAMPLES::

            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo.has_conjugate_in_classP()
            True
            sage: (fibo^2).is_in_classP()
            False
            sage: (fibo^2).has_conjugate_in_classP()
            True
        """
        for k in self.list_of_conjugates():
            if k.is_in_classP(f=f):
                return True
        return False

    def dual_map(self, k=1):
        r"""
        Return the dual map `E_k^*` of self (see [1]).

        .. NOTE::

            It is acually implemented only for `k=1`.

        INPUT:

        - ``self`` - unimodular endomorphism defined on integers
          ``1, 2, \ldots, d``
        - ``k`` - integer (optional, default: 1)

        OUTPUT:

            an instance of E1Star - the dual map

        EXAMPLES::

            sage: sigma = WordMorphism({1:[2],2:[3],3:[1,2]})
            sage: sigma.dual_map()
            E_1^*(1->2, 2->3, 3->12)

        ::

            sage: sigma.dual_map(k=2)
            Traceback (most recent call last):
            ...
            NotImplementedError: The dual map E_k^* is implemented only for k = 1 (not 2)

        REFERENCES:

        - [1] Sano, Y., Arnoux, P. and Ito, S., Higher dimensional
          extensions of substitutions and their dual maps, Journal
          d'Analyse Mathematique 83 (2001), 183-206.
        """
        if k == 1:
            from sage.combinat.e_one_star import E1Star
            return E1Star(self)
        else:
            raise NotImplementedError("The dual map E_k^*" +
                 " is implemented only for k = 1 (not %s)" % k)

    @cached_method
    def rauzy_fractal_projection(self, eig=None, prec=53):
        r"""
        Returns a dictionary giving the projection of the canonical basis.

        See the method :meth:`rauzy_fractal_plot` for more details about the projection.

        INPUT:

        - ``eig`` - a real element of ``QQbar`` of degree >= 2 (default: ``None``).
          The eigenvalue used for the projection.
          It must be an eigenvalue of ``self.incidence_matrix()``.
          The one used by default is the maximal eigenvalue of
          ``self.incidence_matrix()`` (usually a Pisot number),
          but for substitutions with more than 3 letters
          other interesting choices are sometimes possible.

        - ``prec`` - integer (default: ``53``).
          The number of bits used in the floating point representations
          of the coordinates.

        OUTPUT:

            dictionary, letter -> vector, giving the projection

        EXAMPLES:

        The projection for the Rauzy fractal of the Tribonacci substitution
        is::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: s.rauzy_fractal_projection()
            {'1': (1.00000000000000, 0.000000000000000),
             '2': (-1.41964337760708, -0.606290729207199),
             '3': (-0.771844506346038, 1.11514250803994)}

        TESTS::

            sage: t = WordMorphism('1->12,2->3,3->45,4->5,5->6,6->7,7->8,8->1')
            sage: E = t.incidence_matrix().eigenvalues()
            sage: x = [x for x in E if -0.8 < x < -0.7][0]
            sage: t.rauzy_fractal_projection(prec=10)
            {'1': (1.0, 0.00),
             '2': (-1.7, -0.56),
             '3': (0.79, 1.3),
             '4': (1.9, -0.74),
             '5': (-1.7, -0.56),
             '6': (0.79, 1.3),
             '7': (0.21, -1.3),
             '8': (-0.88, 0.74)}
            sage: t.rauzy_fractal_projection(eig=x, prec=10)
            {'1': (1.0, 0.00),
             '2': (-0.12, -0.74),
             '3': (-0.66, -0.56),
             '4': (-0.46, -0.18),
             '5': (-0.54, 0.18),
             '6': (-0.34, 0.56),
             '7': (0.12, 0.74),
             '8': (0.66, 0.56)}

        AUTHOR:

            Timo Jolivet (2012-06-16)
        """
        alphabet = self.domain().alphabet()
        size_alphabet = len(alphabet)

        # Eigenvalues
        if eig is None:
            beta = max(self.incidence_matrix().eigenvalues(), key=abs)
        else:
            beta = eig

        # Test is deg(beta) >= 2
        if beta.degree() < 2:
            raise ValueError("The algebraic degree of ``eig`` must be at least two.")

        # Algebraic conjugates of beta
        from sage.rings.qqbar import QQbar
        beta_conjugates = beta.minpoly().roots(QQbar, multiplicities=False)
        if not beta.imag():
            beta_conjugates.remove(beta)
        for x in beta_conjugates:
            if x.imag():
                beta_conjugates.remove(x.conjugate())

        # Left eigenvector vb in the number field Q(beta)
        from sage.rings.number_field.number_field import NumberField
        K = NumberField(beta.minpoly(), 'b')
        vb = (self.incidence_matrix()-K.gen()).kernel().basis()[0]

        # Projections of canonical base vectors from R^size_alphabet to C, using vb
        from sage.modules.free_module import VectorSpace
        canonical_basis = VectorSpace(K,size_alphabet).basis()
        canonical_basis_proj = {}

        from sage.rings.real_mpfr import RealField
        RealField_prec = RealField(prec)
        for a, x in zip(alphabet, canonical_basis):
            v = []
            for y in beta_conjugates:
                # if y has nonzero imaginary part
                if y.imag():
                    z = (vb*x).lift()(y)
                    z1, z2 = z.real(), z.imag()
                    v += [RealField_prec(z1), RealField_prec(z2)]
                # if y is real
                else:
                    z = (vb*x).lift()(y)
                    v += [RealField_prec(z)]
            canonical_basis_proj[a] = vector(v)

        return canonical_basis_proj

    def rauzy_fractal_points(self, n=None, exchange=False, eig=None, translate=None, prec=53):
        r"""
        Returns a dictionary of list of points associated with the pieces
        of the Rauzy fractal of ``self``.

        INPUT:

            See the method :meth:`rauzy_fractal_plot` for a description
            of the options and more examples.

        OUTPUT:

            dictionary of list of points

        EXAMPLES:

        The Rauzy fractal of the Tribonacci substitution and the number of
        points in the piece of the fractal associated with ``'1'``, ``'2'``
        and ``'3'`` are respectively::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: D = s.rauzy_fractal_points(n=100)
            sage: len(D['1'])
            54
            sage: len(D['2'])
            30
            sage: len(D['3'])
            16

        TESTS::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: D = s.rauzy_fractal_points(n=100, exchange=True, translate=[(3,1,-2), (5,-33,8)], prec=40)
            sage: len(D['1'])
            108

        AUTHOR:

            Timo Jolivet (2012-06-16)
        """
        alphabet = self.domain().alphabet()
        canonical_basis_proj = self.rauzy_fractal_projection(eig=eig, prec=prec)

        # if exchange, set the projection to its opposite
        if exchange:
            for a in canonical_basis_proj:
                canonical_basis_proj[a] = - canonical_basis_proj[a]

        # Compute a fixed point u
        if exchange:
            u = iter(self.reversal().periodic_points()[0][0])
        else:
            u = iter(self.periodic_points()[0][0])

        # Manage various options in function of dimension
        if n is None:
            dim_fractal = len(canonical_basis_proj[alphabet[0]])
            if dim_fractal == 1:
                n = 1000
            elif dim_fractal == 2:
                n = 50000
            elif dim_fractal == 3:
                n = 5000
            else:
                n = 50000

        # Compute orbit points to plot
        S = 0
        orbit_points = dict([(a,[]) for a in alphabet])
        for _ in xrange(n):
            a = next(u)
            S += canonical_basis_proj[a]
            orbit_points[a].append(S)

        # Manage translated copies
        from sage.rings.real_mpfr import RealField
        RealField_prec = RealField(prec)
        if translate is not None:

            if isinstance(translate, dict):
                for a in translate.keys():
                    translate[a] = [vector(RealField_prec, v) for v in translate[a]]

            else:
                translate = [vector(RealField_prec, v) for v in translate]

            for a in alphabet:
                translated_copies = dict([(i,[]) for i in alphabet])

                if isinstance(translate, list):
                    to_treat = translate

                elif isinstance(translate, dict):
                    try:
                        to_treat = translate[a]
                    except KeyError:
                        to_treat = []

                for x in to_treat:
                    v = 0
                    for i,z in zip(alphabet,x):
                        v += z*canonical_basis_proj[i]
                    translated_copies[a] += [vector(v) + w for w in orbit_points[a]]

                orbit_points[a] = translated_copies[a]

        return orbit_points

    def rauzy_fractal_plot(self, n=None, exchange=False, eig=None, translate=None, prec=53, \
                           colormap='hsv', opacity=None, plot_origin=None, plot_basis=False, point_size=None):
        r"""
        Returns a plot of the Rauzy fractal associated with a substitution.

        The substitution does not have to be irreducible.
        The usual definition of a Rauzy fractal requires that
        its dominant eigenvalue is a Pisot number but the present method
        doesn't require this, allowing to plot some interesting pictures
        in the non-Pisot case (see the examples below).

        For more details about the definition of the fractal and the
        projection which is used, see Section 3.1 of [1].

        Plots with less than 100,000 points take a few seconds,
        and several millions of points can be plotted in reasonable time.

        Other ways to draw Rauzy fractals (and more generally projections of paths)
        can be found in :meth:`sage.combinat.words.paths.FiniteWordPath_all.plot_projection`
        or in :meth:`sage.combinat.e_one_star`.

        OUTPUT:

        A Graphics object.

        INPUT:

        - ``n`` - integer (default: ``None``)
          The number of points used to plot the fractal.
          Default values: ``1000`` for a 1D fractal,
          ``50000`` for a 2D fractal, ``10000`` for a 3D fractal.

        - ``exchange`` - boolean (default: ``False``).
          Plot the Rauzy fractal with domain exchange.

        - ``eig`` - a real element of ``QQbar`` of degree >= 2 (default: ``None``).
          The eigenvalue used to plot the fractal.
          It must be an eigenvalue of ``self.incidence_matrix()``.
          The one used by default the maximal eigenvalue of
          ``self.incidence_matrix()`` (usually a Pisot number),
          but for substitutions with more than 3 letters
          other interesting choices are sometimes possible.

        - ``translate`` - a list of vectors of ``RR^size_alphabet``,
          or a dictionary from the alphabet to lists of vectors (default: ``None``).
          Plot translated copies of the fractal.
          This option allows to plot tilings easily.
          The projection used for these vectors is the same as
          the projection used for the canonical basis to plot the fractal.
          If the input is a list, all the pieces will be translated and plotted.
          If the input is a dictionary, each piece will be translated and plotted
          accordingly to the vectors associated with each letter in the dictionary.
          Note: by default, the Rauzy fractal placed at the origin
          is not plotted with the ``translate`` option;
          the vector ``(0,0,...,0)`` has to be added manually.

        - ``prec`` - integer (default: ``53``).
          The number of bits used in the floating point representations
          of the points of the fractal.

        - ``colormap`` - color map or dictionary (default: ``'hsv'``).
          It can be one of the following :

           - ``string`` - a coloring map. For available coloring map names type:
             ``sorted(colormaps)``

           - ``dict`` - a dictionary of the alphabet mapped to colors.

        - ``opacity`` - a dictionary from the alphabet to the real interval [0,1] (default: ``None``).
          If none is specified, all letters are plotted with opacity ``1``.

        - ``plot_origin`` - a couple ``(k,c)`` (default: ``None``).
          If specified, mark the origin by a point of size ``k`` and color ``c``.

        - ``plot_basis`` - boolean (default: ``False``).
          Plot the projection of the canonical basis with the fractal.

        - ``point_size`` - float (default: ``None``).
          The size of the points used to plot the fractal.

        EXAMPLES:

        #. The Rauzy fractal of the Tribonacci substitution::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: s.rauzy_fractal_plot()     # long time
            Graphics object consisting of 3 graphics primitives

        #. The "Hokkaido" fractal. We tweak the plot using the plotting options
           to get a nice reusable picture, in which we mark the origin by a black dot::

            sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
            sage: G = s.rauzy_fractal_plot(n=100000, point_size=3, plot_origin=(50,"black"))  # not tested
            sage: G.show(figsize=10, axes=false) # not tested

        #. Another "Hokkaido" fractal and its domain exchange::

            sage: s = WordMorphism({1:[2], 2:[4,3], 3:[4], 4:[5,3], 5:[6], 6:[1]})
            sage: s.rauzy_fractal_plot()                  # not tested takes > 1 second
            sage: s.rauzy_fractal_plot(exchange=True)     # not tested takes > 1 second

        #. A three-dimensional Rauzy fractal::

            sage: s = WordMorphism('1->12,2->13,3->14,4->1')
            sage: s.rauzy_fractal_plot()     # not tested takes > 1 second

        #. A one-dimensional Rauzy fractal (very scattered)::

            sage: s = WordMorphism('1->2122,2->1')
            sage: s.rauzy_fractal_plot().show(figsize=20)     # not tested takes > 1 second

        #. A high resolution plot of a complicated fractal::

            sage: s = WordMorphism('1->23,2->123,3->1122233')
            sage: G = s.rauzy_fractal_plot(n=300000)  # not tested takes > 1 second
            sage: G.show(axes=false, figsize=20)      # not tested takes > 1 second

        #. A nice colorful animation of a domain exchange::

            sage: s = WordMorphism('1->21,2->3,3->4,4->25,5->6,6->7,7->1')
            sage: L = [s.rauzy_fractal_plot(), s.rauzy_fractal_plot(exchange=True)]     # not tested takes > 1 second
            sage: animate(L, axes=false).show(delay=100)     # not tested takes > 1 second

        #. Plotting with only one color::

            sage: s = WordMorphism('1->12,2->31,3->1')
            sage: s.rauzy_fractal_plot(colormap={'1':'black', '2':'black', '3':'black'})     # not tested takes > 1 second

        #. Different fractals can be obtained by choosing another (non-Pisot) eigenvalue::

            sage: s = WordMorphism('1->12,2->3,3->45,4->5,5->6,6->7,7->8,8->1')
            sage: E = s.incidence_matrix().eigenvalues()
            sage: x = [x for x in E if -0.8 < x < -0.7][0]
            sage: s.rauzy_fractal_plot()          # not tested takes > 1 second
            sage: s.rauzy_fractal_plot(eig=x)     # not tested takes > 1 second

        #. A Pisot reducible substitution with seemingly overlapping tiles::

            sage: s = WordMorphism({1:[1,2], 2:[2,3], 3:[4], 4:[5], 5:[6], 6:[7], 7:[8], 8:[9], 9:[10], 10:[1]})
            sage: s.rauzy_fractal_plot()     # not tested takes > 1 second

        #. A non-Pisot reducible substitution with a strange Rauzy fractal::

            sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
            sage: s.rauzy_fractal_plot()     # not tested takes > 1 second

        #. A substitution with overlapping tiles. We use the options
           ``colormap`` and ``opacity`` to study how the tiles overlap::

            sage: s = WordMorphism('1->213,2->4,3->5,4->1,5->21')
            sage: s.rauzy_fractal_plot()                                   # not tested takes > 1 second
            sage: s.rauzy_fractal_plot(colormap={'1':'red', '4':'purple'})     # not tested takes > 1 second
            sage: s.rauzy_fractal_plot(opacity={'1':0.1,'2':1,'3':0.1,'4':0.1,'5':0.1}, n=150000)     # not tested takes > 1 second

        #. Funny experiments by playing with the precision of the float numbers used to plot the fractal::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: s.rauzy_fractal_plot(prec=6)      # not tested
            sage: s.rauzy_fractal_plot(prec=9)      # not tested
            sage: s.rauzy_fractal_plot(prec=15)     # not tested
            sage: s.rauzy_fractal_plot(prec=19)     # not tested
            sage: s.rauzy_fractal_plot(prec=25)     # not tested

        #. Using the ``translate`` option to plot periodic tilings::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: s.rauzy_fractal_plot(n=10000, translate=[(0,0,0),(-1,0,1),(0,-1,1),(1,-1,0),(1,0,-1),(0,1,-1),(-1,1,0)])     # not tested takes > 1 second

           ::

            sage: t = WordMorphism("a->aC,b->d,C->de,d->a,e->ab")   # substitution found by Julien Bernat
            sage: V = [vector((0,0,1,0,-1)), vector((0,0,1,-1,0))]
            sage: S = set(map(tuple, [i*V[0] + j*V[1] for i in [-1,0,1] for j in [-1,0,1]]))
            sage: t.rauzy_fractal_plot(n=10000, translate=S, exchange=true)     # not tested takes > 1 second

        #. Using the ``translate`` option to plot arbitrary tilings with the fractal pieces.
           This can be used for example to plot the self-replicating tiling of the Rauzy fractal::

            sage: s = WordMorphism({1:[1,2], 2:[3], 3:[4,3], 4:[5], 5:[6], 6:[1]})
            sage: s.rauzy_fractal_plot()     # not tested takes > 1 second
            sage: D = {1:[(0,0,0,0,0,0), (0,1,0,0,0,0)], 3:[(0,0,0,0,0,0), (0,1,0,0,0,0)], 6:[(0,1,0,0,0,0)]}
            sage: s.rauzy_fractal_plot(n=30000, translate=D)     # not tested takes > 1 second

        #. Plot the projection of the canonical basis with the fractal::

            sage: s = WordMorphism({1:[2,1], 2:[3], 3:[6,4], 4:[5,1], 5:[6], 6:[7], 7:[8], 8:[9], 9:[1]})
            sage: s.rauzy_fractal_plot(plot_basis=True)     # not tested takes > 1 second

        TESTS::

            sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
            sage: s.rauzy_fractal_plot(n=1000, colormap='Set1', opacity={'a':0.5,'b':1,'c':0.7,'d':0,'e':0.2}, plot_origin=(100,"black"), plot_basis=True, point_size=2.5)
            Graphics object consisting of 10 graphics primitives

        REFERENCES:

        - [1] Valerie Berthe and Anne Siegel,
          Tilings associated with beta-numeration and substitutions,
          Integers 5 (3), 2005.
          http://www.integers-ejcnt.org/vol5-3.html

        AUTHOR:

            Timo Jolivet (2012-06-16)
        """
        alphabet = self.domain().alphabet()
        size_alphabet = len(alphabet)

        orbit_points = self.rauzy_fractal_points(n=n, exchange=exchange, eig=eig, translate=translate, prec=prec)

        dim_fractal = len(orbit_points[alphabet[0]][0])

        # Manage colors and opacity
        if isinstance(colormap, dict):
            col_dict = colormap

        elif isinstance(colormap, str):
            from matplotlib import cm

            if not colormap in cm.datad.keys():
                raise RuntimeError("Color map %s not known (type sorted(colors) for valid names)" % colormap)

            colormap = cm.__dict__[colormap]
            col_dict = {}
            for i, a in enumerate(alphabet):
                col_dict[a] = colormap(float(i)/float(size_alphabet))[:3]

        else:
            raise TypeError("Type of option colormap (=%s) must be dict or str" % colormap)

        if opacity is None:
            opacity = dict([(a,1) for a in alphabet])

        elif not isinstance(opacity, dict):
            raise TypeError("Type of option opacity (=%s) must be dict" % opacity)

        # Plot points size
        if point_size is None:
            if dim_fractal == 1 or dim_fractal == 2:
                point_size = 1
            elif dim_fractal == 3:
                point_size = 8

        # Make graphics
        from sage.plot.plot import Graphics
        G = Graphics()

        from sage.plot.point import points

        # 1D plots
        if dim_fractal == 1:
            from sage.all import plot
            for a in col_dict.keys():
                # We plot only the points with a color in col_dict and with positive opacity
                if (a in col_dict.keys()) and (opacity[a] > 0):
                    G += plot([x[0] for x in orbit_points[a]], color=col_dict[a], alpha=opacity[a], thickness=point_size)
            if plot_basis:
                from matplotlib import cm
                from sage.plot.arrow import arrow
                canonical_basis_proj = self.rauzy_fractal_projection(eig=eig, prec=prec)
                for i,a in enumerate(alphabet):
                    x = canonical_basis_proj[a]
                    G += arrow((-1.1,0), (-1.1,x[0]), color=cm.__dict__["gist_gray"](0.75*float(i)/float(size_alphabet))[:3])

        # 2D or 3D plots
        else:
            if point_size is None and dim_fractal == 2:
                point_size = 1
            elif point_size is None and dim_fractal == 3:
                point_size = 8

            for a in col_dict.keys():
                # We plot only the points with a color in col_dict and with positive opacity
                if (a in col_dict.keys()) and (opacity[a] > 0):
                    G += points(orbit_points[a], color=col_dict[a], alpha=opacity[a], size=point_size)

            if plot_basis:
                from matplotlib import cm
                from sage.plot.arrow import arrow
                canonical_basis_proj = self.rauzy_fractal_projection(eig=eig, prec=prec)
                for i,a in enumerate(alphabet):
                    x = canonical_basis_proj[a]
                    G += arrow([0]*dim_fractal, x, color=cm.__dict__["gist_gray"](0.75*float(i)/float(size_alphabet))[:3])

        if plot_origin:
            G += points([(0,0)], size=plot_origin[0], color=plot_origin[1])

        if dim_fractal == 1 or dim_fractal == 2:
            G.set_aspect_ratio(1)

        return G

    def is_growing(self, letter=None):
        r"""
        Return ``True`` if ``letter`` is a growing letter.

        A letter `a` is *growing* for the morphism `s` if the length of the
        iterates of `| s^n(a) |` tend to infinity as `n` goes to infinity.

        INPUT:

        - ``letter`` -- ``None`` or a letter in the domain of ``self``

        .. NOTE::

            If letter is ``None``, this returns ``True`` if ``self`` is
            everywhere growing, i.e., all letters are growing letters (see
            [CassNic10]_), and that ``self`` **must** be an endomorphism.

        EXAMPLES::

            sage: WordMorphism('0->01,1->1').is_growing('0')
            True
            sage: WordMorphism('0->01,1->1').is_growing('1')
            False
            sage: WordMorphism('0->01,1->10').is_growing()
            True
            sage: WordMorphism('0->1,1->2,2->01').is_growing()
            True
            sage: WordMorphism('0->01,1->1').is_growing()
            False

        The domain needs to be equal to the codomain::

            sage: WordMorphism('0->01,1->0,2->1',codomain=Words('012')).is_growing()
            True

        Test of erasing morphisms::

            sage: WordMorphism('0->01,1->').is_growing('0')
            False
            sage: m = WordMorphism('a->bc,b->bcc,c->',codomain=Words('abc'))
            sage: m.is_growing('a')
            False
            sage: m.is_growing('b')
            False
            sage: m.is_growing('c')
            False

        REFERENCES:

        ..  [CassNic10] Cassaigne J., Nicolas F. Factor complexity.
            Combinatorics, automata and number theory, 163--247, Encyclopedia
            Math. Appl., 135, Cambridge Univ. Press, Cambridge, 2010.
        """
        if self.is_primitive():
            return True
        if letter is None:
            I = range(self.domain().alphabet().cardinality())
        else:
            if not letter in self.domain().alphabet():
                raise TypeError("letter (=%s) is not in the domain of self"%letter)
            I = [self.domain().alphabet().rank(letter)]

        last_coef = 0
        coefs = self.incidence_matrix().charpoly().coefficients(sparse=False)
        while coefs[last_coef] == 0:
            last_coef += 1
        V = self.abelian_rotation_subspace() + (self.incidence_matrix()**last_coef).right_kernel().change_ring(QQ)
        basis = V.ambient_vector_space().basis()

        return not any(basis[i] in V for i in I)

    def growing_letters(self):
        r"""
        Returns the list of growing letters.

        See :meth:`.is_growing` for more information.

        EXAMPLES::

            sage: WordMorphism('0->01,1->10').growing_letters()
            ['0', '1']
            sage: WordMorphism('0->01,1->1').growing_letters()
            ['0']
            sage: WordMorphism('0->01,1->0,2->1',codomain=Words('012')).growing_letters()
            ['0', '1', '2']
        """
        if self.is_primitive():
            return self.domain().alphabet().list()
        last_coef = 0
        coefs = self.incidence_matrix().charpoly().coefficients(sparse=False)
        while coefs[last_coef] == 0:
            last_coef += 1
        V = self.abelian_rotation_subspace() + (self.incidence_matrix()**last_coef).right_kernel().change_ring(QQ)
        basis = V.ambient_vector_space().basis()
        A = self.domain().alphabet()

        return list(A.unrank(i) for i in range(A.cardinality()) if basis[i] not in V)

    def abelian_rotation_subspace(self):
        r"""
        Returns the subspace on which the incidence matrix of ``self`` acts by
        roots of unity.

        EXAMPLES::

            sage: WordMorphism('0->1,1->0').abelian_rotation_subspace()
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
            sage: WordMorphism('0->01,1->10').abelian_rotation_subspace()
            Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []
            sage: WordMorphism('0->01,1->1').abelian_rotation_subspace()
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]
            sage: WordMorphism('1->122,2->211').abelian_rotation_subspace()
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -1]
            sage: WordMorphism('0->1,1->102,2->3,3->4,4->2').abelian_rotation_subspace()
            Vector space of degree 5 and dimension 3 over Rational Field
            Basis matrix:
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]

        The domain needs to be equal to the codomain::

            sage: WordMorphism('0->1,1->',codomain=Words('01')).abelian_rotation_subspace()
            Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        if not self.domain() == self.codomain():
            raise TypeError("self (=%s) is not an endomorphism"%self)

        if self.domain().alphabet().cardinality() == Infinity:
            raise ValueError("the alphabet is infinite")

        M = self.incidence_matrix()
        p = M.charpoly().factor()
        basis = []
        for factor in p:
            if factor[0].is_cyclotomic():
                basis.extend((factor[0])(M).right_kernel().basis())

        return M._column_ambient_module().change_ring(QQ).subspace(basis)
