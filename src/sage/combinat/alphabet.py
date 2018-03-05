class Letter:


    def __init__(self, value, is_even):
        self.value = value
        self.is_even = is_even

    @property
    def is_odd(self):
        return not self.is_even

    def __eq__(self, other):
        return self.value == other.value and self.is_even == other.is_even

    def __lt__(self, other):
        # value tie
        if self.value == other.value:
            # both of same evenness
            if self.is_even == other.is_even:
                return False
            # different evenness
            elif self.is_even and (not other.is_even):
                return False
            elif (not self.is_even) and other.is_even:
                return True
            else:
                raise ValueError
        # no value tie
        else:
            return self.value < other.value

    def __le__(self, other):
        return self < other or self == other

    def __gt__(self, other):
        return not self <= other

    def __ge__(self, other):
        return not self < other

    def __str__(self):
        if self.is_even:
            s = "even({})"
        else:
            s = "odd({})"
        return s.format(self.value)

    def __repr__(self):
        return str(self)

def even(n):
    return Letter(n, True)

def odd(n):
    return Letter(n, False)


class BiLetter:


    def __init__(self, a, b):
        self.value = (a, b)

    def __repr__(self):
        return '{}'.format(self.value)

    def __lt__(self, other):
        return self.value[1] < other.value[1] or (self.value[1] == other.value[1] and self.value[0] < other.value[0])

    def __le__(self, other):
        return self < other or self == other

    def __gt__(self, other):
        return not self <= other

    def __ge__(self, other):
        return not self < other

    def is_mixed(self):
        return self.value[0].is_even != self.value[1].is_even

    def inverse(self):
        (a, b) = self.value
        return BiLetter(b, a)


from collections import Counter
class BiWord:


    def __init__(self, obj1, obj2=None):
        if obj2 is None:
            # it must be a list of biletters or tuples
            if not obj1:
                # it's empty
                self.value = []
            elif obj1 and isinstance(obj1[0], tuple):
                # it's a list of tuples
                self.value = [BiLetter(x, y) for x,y in obj1]
            elif obj1 and isinstance(obj1[0], BiLetter):
                # it's a list of biletters
                self.value = obj1
            elif isinstance(obj1, list) and isinstance(obj2, list):
                # it's two separate words
                pairs = zip(obj1, obj2)
                self.value = [BiLetter(x, y) for x,y in pairs]
            else:
                raise ValueError('unsupported input format')

    def __repr__(self):
        return 'BiWord({})'.format(self.value)

    def __len__(self):
        return len(self.to_list())

    def is_restricted(self):
        mixed_biletter_lis = [e for e in self.value if e.is_mixed()]
        mixed_biletter_set = set(mixed_biletter_lis)
        return len(mixed_biletter_lis) == len(mixed_biletter_set)

    def is_ordered(self):
        """
        TESTS:

            sage: BiWord([(even(1), odd(3)), (even(1), even(3))]).is_ordered()
            True
            sage: BiWord([(odd(1), even(3)), (even(1), even(3))]).is_ordered()
            True
            sage: BiWord([(even(1), even(3)), (even(1), even(3))]).is_ordered()
            True
            sage: BiWord([(even(1), even(3)), (even(1), odd(3))]).is_ordered()
            False
            sage: BiWord([(odd(5), even(2)), (odd(4), even(2))]).is_ordered()
            False
        """
        return self.value == sorted(self.value)

    def lword(self):
        return [b.value[0] for b in self.value]

    def lcon(self):
        return Counter(self.lword())

    def rword(self):
        return [b.value[1] for b in self.value]

    def rcon(self):
        return Counter(self.lword())

    def upside_down(self):
        ud = [b.inverse() for b in self.value]
        return BiWord(ud)

    def to_list(self):
        return self.value

    def sorted(self):
        return BiWord(sorted(self.to_list()))

    def inverse(self):
        """
        EXAMPLES:

            sage: BiWord([(odd(5), even(2)), (odd(4), even(2))]).inverse()
            BiWord([(even(2), odd(4)), (even(2), odd(5))])

        TESTS:

            sage: BiWord([(odd(1), even(2)), (odd(2), even(3))]).inverse()
            BiWord([(even(2), odd(1)), (even(3), odd(2))])
        """
        return self.upside_down().sorted()

    def _rsk_iter(self):
        return [(b.value[1], b.value[0]) for b in self.value]
