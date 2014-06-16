"""
Tools to combine dependencies
"""


class CommutingSha1Accumulator(object):
    """
    Accumulate SHA1 checksums disregarding their order.

    EXAMPLES::

        >>> from sage_pkg.package.base import CommutingSha1Accumulator
        >>> acc = CommutingSha1Accumulator();  str(acc)
        '0000000000000000000000000000000000000000'
        >>> acc += 'bbd8ea60604b02bf7e33c8c3ca86c2b6ac032b3e'
        >>> acc += '2cdba115ff4cb73bb14bb02881e31f3dfe871ec6'
        >>> str(acc)
        'e8b48b765f97b9fb2f7f78ec4c69e1f4aa8a4a04'
        >>> rev = CommutingSha1Accumulator()
        >>> rev += '2cdba115ff4cb73bb14bb02881e31f3dfe871ec6'
        >>> rev += 'bbd8ea60604b02bf7e33c8c3ca86c2b6ac032b3e'
        >>> str(acc) == str(rev)
        True
    """

    def __init__(self):
        self.value = 0
    
    def __iadd__(self, sha1):
        if not len(sha1) == 40:
            raise ValueError('sha1 must be 40 digits')
        self.value += int(sha1, 16)
        return self

    def __repr__(self):
        return  '{0:040x}'.format(self.value)[-40:]

