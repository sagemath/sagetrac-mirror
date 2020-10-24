r"""
Rolfsen database of knots with at most 10 crossings.

Each entry is indexed by a pair `(n, k)` where `n` is the crossing number
and `k` is the index in the table.

Knots are described by words in braids groups, more precisely by pairs
`(m, w)` where `m` is the number of strands and `w` is a list of
indices of generators (negative for the inverses).

This data has been obtained from the [KnotAtlas]_ ,
with the kind permission of Dror Bar-Natan.

All knots in the Rolfsen table are depicted nicely in the following page:

http://katlas.math.toronto.edu/wiki/The_Rolfsen_Knot_Table

.. NOTE::

    Some of the knots are not represented by a diagram with the
    minimal number of crossings.

.. TODO::

    - Extend the table to knots with more crossings.

    - Do something similar for links using
      http://katlas.org/wiki/The_Thistlethwaite_Link_Table
"""


small_knots_table = {
    (0, 1): (1, []),
    (3, 1): (2, [-1,-1,-1]),
    (4, 1): (3, [-1,2,-1,2]),
    (5, 1): (2, [-1,-1,-1,-1,-1]),
    (5, 2): (3, [-1,-1,-1,-2,1,-2]),
    (6, 1): (4, [-1,-1,-2,1,3,-2,3]),
    (6, 2): (3, [-1,-1,-1,2,-1,2]),
    (6, 3): (3, [-1,-1,2,-1,2,2]),
    (7, 1): (2, [-1,-1,-1,-1,-1,-1,-1]),
    (7, 2): (4, [-1,-1,-1,-2,1,-2,-3,2,-3]),
    (7, 3): (3, [1,1,1,1,1,2,-1,2]),
    (7, 4): (4, [1,1,2,-1,2,2,3,-2,3]),
    (7, 5): (3, [-1,-1,-1,-1,-2,1,-2,-2]),
    (7, 6): (4, [-1,-1,2,-1,-3,2,-3]),
    (7, 7): (4, [1,-2,1,-2,3,-2,3]),
    (8, 1): (5, [-1,-1,-2,1,-2,-3,2,4,-3,4]),
    (8, 2): (3, [-1,-1,-1,-1,-1,2,-1,2]),
    (8, 3): (5, [-1,-1,-2,1,3,-2,3,4,-3,4]),
    (8, 4): (4, [-1,-1,-1,2,-1,2,3,-2,3]),
    (8, 5): (3, [1,1,1,-2,1,1,1,-2]),
    (8, 6): (4, [-1,-1,-1,-1,-2,1,3,-2,3]),
    (8, 7): (3, [1,1,1,1,-2,1,-2,-2]),
    (8, 8): (4, [1,1,1,2,-1,-3,2,-3,-3]),
    (8, 9): (3, [-1,-1,-1,2,-1,2,2,2]),
    (8, 10): (3, [1,1,1,-2,1,1,-2,-2]),
    (8, 11): (4, [-1,-1,-2,1,-2,-2,3,-2,3]),
    (8, 12): (5, [-1,2,-1,-3,2,4,-3,4]),
    (8, 13): (4, [-1,-1,2,-1,2,2,3,-2,3]),
    (8, 14): (4, [-1,-1,-1,-2,1,-2,3,-2,3]),
    (8, 15): (4, [-1,-1,2,-1,-3,-2,-2,-2,-3]),
    (8, 16): (3, [-1,-1,2,-1,-1,2,-1,2]),
    (8, 17): (3, [-1,-1,2,-1,2,-1,2,2]),
    (8, 18): (3, [-1,2,-1,2,-1,2,-1,2]),
    (8, 19): (3, [1,1,1,2,1,1,1,2]),
    (8, 20): (3, [1,1,1,-2,-1,-1,-1,-2]),
    (8, 21): (3, [-1,-1,-1,-2,1,1,-2,-2]),
    (9, 1): (2, [-1,-1,-1,-1,-1,-1,-1,-1,-1]),
    (9, 2): (5, [-1,-1,-1,-2,1,-2,-3,2,-3,-4,3,-4]),
    (9, 3): (3, [1,1,1,1,1,1,1,2,-1,2]),
    (9, 4): (4, [-1,-1,-1,-1,-1,-2,1,-2,-3,2,-3]),
    (9, 5): (5, [1,1,2,-1,2,2,3,-2,3,4,-3,4]),
    (9, 6): (3, [-1,-1,-1,-1,-1,-1,-2,1,-2,-2]),
    (9, 7): (4, [-1,-1,-1,-1,-2,1,-2,-3,2,-3,-3]),
    (9, 8): (5, [-1,-1,2,-1,2,3,-2,-4,3,-4]),
    (9, 9): (3, [-1,-1,-1,-1,-1,-2,1,-2,-2,-2]),
    (9, 10): (4, [1,1,2,-1,2,2,2,2,3,-2,3]),
    (9, 11): (4, [1,1,1,1,-2,1,3,-2,3]),
    (9, 12): (5, [-1,-1,2,-1,-3,2,-3,-4,3,-4]),
    (9, 13): (4, [1,1,1,1,2,-1,2,2,3,-2,3]),
    (9, 14): (5, [1,1,2,-1,-3,2,-3,4,-3,4]),
    (9, 15): (5, [1,1,1,2,-1,-3,2,4,-3,4]),
    (9, 16): (3, [1,1,1,1,2,2,-1,2,2,2]),
    (9, 17): (4, [1,-2,1,-2,-2,-2,3,-2,3]),
    (9, 18): (4, [-1,-1,-1,-2,1,-2,-2,-2,-3,2,-3]),
    (9, 19): (5, [1,-2,1,-2,-2,-3,2,4,-3,4]),
    (9, 20): (4, [-1,-1,-1,2,-1,-3,2,-3,-3]),
    (9, 21): (5, [1,1,2,-1,2,-3,2,4,-3,4]),
    (9, 22): (4, [-1,2,-1,2,-3,2,2,2,-3]),
    (9, 23): (4, [-1,-1,-1,-2,1,-2,-2,-3,2,-3,-3]),
    (9, 24): (4, [-1,-1,2,-1,-3,2,2,2,-3]),
    (9, 25): (5, [-1,-1,2,-1,-3,-2,-2,4,-3,4]),
    (9, 26): (4, [1,1,1,-2,1,-2,3,-2,3]),
    (9, 27): (4, [-1,-1,2,-1,2,2,-3,2,-3]),
    (9, 28): (4, [-1,-1,2,-1,-3,2,2,-3,-3]),
    (9, 29): (4, [1,-2,-2,3,-2,1,-2,3,-2]),
    (9, 30): (4, [-1,-1,2,2,-1,2,-3,2,-3]),
    (9, 31): (4, [-1,-1,2,-1,2,-3,2,-3,-3]),
    (9, 32): (4, [1,1,-2,1,-2,1,3,-2,3]),
    (9, 33): (4, [-1,2,-1,2,2,-1,-3,2,-3]),
    (9, 34): (4, [-1,2,-1,2,-3,2,-1,2,-3]),
    (9, 35): (5, [-1,-1,-2,1,-2,-2,-3,2,2,-4,3,-2,-4,-3]),
    (9, 36): (4, [1,1,1,-2,1,1,3,-2,3]),
    (9, 37): (5, [-1,-1,2,-1,-3,2,1,4,-3,2,-3,4]),
    (9, 38): (4, [-1,-1,-2,-2,3,-2,1,-2,-3,-3,-2]),
    (9, 39): (5, [1,1,2,-1,-3,-2,1,4,3,-2,3,4]),
    (9, 40): (4, [-1,2,-1,-3,2,-1,-3,2,-3]),
    (9, 41): (5, [-1,-1,-2,1,3,2,2,-4,-3,2,-3,-4]),
    (9, 42): (4, [1,1,1,-2,-1,-1,3,-2,3]),
    (9, 43): (4, [1,1,1,2,1,1,-3,2,-3]),
    (9, 44): (4, [-1,-1,-1,-2,1,1,3,-2,3]),
    (9, 45): (4, [-1,-1,-2,1,-2,-1,-3,2,-3]),
    (9, 46): (4, [-1,2,-1,2,-3,-2,1,-2,-3]),
    (9, 47): (4, [-1,2,-1,2,3,2,-1,2,3]),
    (9, 48): (4, [1,1,2,-1,2,1,-3,2,-1,2,-3]),
    (9, 49): (4, [1,1,2,1,1,-3,2,-1,2,3,3]),
    (10, 1): (6, [-1,-1,-2,1,-2,-3,2,-3,-4,3,5,-4,5]),
    (10, 2): (3, [-1,-1,-1,-1,-1,-1,-1,2,-1,2]),
    (10, 3): (6, [-1,-1,-2,1,-2,-3,2,4,-3,4,5,-4,5]),
    (10, 4): (5, [-1,-1,-1,2,-1,2,3,-2,3,4,-3,4]),
    (10, 5): (3, [1,1,1,1,1,1,-2,1,-2,-2]),
    (10, 6): (4, [-1,-1,-1,-1,-1,-1,-2,1,3,-2,3]),
    (10, 7): (5, [-1,-1,-2,1,-2,-3,2,-3,-3,4,-3,4]),
    (10, 8): (4, [-1,-1,-1,-1,-1,2,-1,2,3,-2,3]),
    (10, 9): (3, [1,1,1,1,1,-2,1,-2,-2,-2]),
    (10, 10): (5, [-1,-1,2,-1,2,2,3,-2,3,4,-3,4]),
    (10, 11): (5, [-1,-1,-1,-1,-2,1,3,-2,3,4,-3,4]),
    (10, 12): (4, [1,1,1,1,1,2,-1,-3,2,-3,-3]),
    (10, 13): (6, [-1,-1,-2,1,3,-2,-4,3,5,-4,5]),
    (10, 14): (4, [-1,-1,-1,-1,-1,-2,1,-2,3,-2,3]),
    (10, 15): (4, [1,1,1,1,-2,1,-2,-3,2,-3,-3]),
    (10, 16): (5, [1,1,2,-1,2,2,-3,2,-3,-4,3,-4]),
    (10, 17): (3, [-1,-1,-1,-1,2,-1,2,2,2,2]),
    (10, 18): (5, [-1,-1,-1,-2,1,-2,3,-2,3,4,-3,4]),
    (10, 19): (4, [-1,-1,-1,-1,2,-1,2,2,3,-2,3]),
    (10, 20): (5, [-1,-1,-1,-1,-2,1,-2,-3,2,4,-3,4]),
    (10, 21): (4, [-1,-1,-2,1,-2,-2,-2,-2,3,-2,3]),
    (10, 22): (4, [1,1,1,1,2,-1,-3,2,-3,-3,-3]),
    (10, 23): (4, [-1,-1,2,-1,2,2,2,2,3,-2,3]),
    (10, 24): (5, [-1,-1,-2,1,-2,-2,-2,-3,2,4,-3,4]),
    (10, 25): (4, [-1,-1,-1,-1,-2,1,-2,-2,3,-2,3]),
    (10, 26): (4, [-1,-1,-1,2,-1,2,2,2,3,-2,3]),
    (10, 27): (4, [-1,-1,-1,-1,-2,1,-2,3,-2,3,3]),
    (10, 28): (5, [1,1,2,-1,2,2,3,-2,-4,3,-4,-4]),
    (10, 29): (5, [-1,-1,-1,2,-1,-3,2,4,-3,4]),
    (10, 30): (5, [-1,-1,-2,1,-2,-2,-3,2,-3,4,-3,4]),
    (10, 31): (5, [-1,-1,-1,-2,1,3,-2,3,3,4,-3,4]),
    (10, 32): (4, [1,1,1,-2,1,-2,-2,-3,2,-3,-3]),
    (10, 33): (5, [-1,-1,-2,1,-2,3,-2,3,3,4,-3,4]),
    (10, 34): (5, [1,1,1,2,-1,2,3,-2,-4,3,-4,-4]),
    (10, 35): (6, [-1,2,-1,2,3,-2,-4,3,5,-4,5]),
    (10, 36): (5, [-1,-1,-1,-2,1,-2,-3,2,-3,4,-3,4]),
    (10, 37): (5, [-1,-1,-1,-2,1,3,-2,3,4,-3,4,4]),
    (10, 38): (5, [-1,-1,-1,-2,1,-2,-2,-3,2,4,-3,4]),
    (10, 39): (4, [-1,-1,-1,-2,1,-2,-2,-2,3,-2,3]),
    (10, 40): (4, [1,1,1,2,-1,2,2,-3,2,-3,-3]),
    (10, 41): (5, [1,-2,1,-2,-2,3,-2,-4,3,-4]),
    (10, 42): (5, [-1,-1,2,-1,2,-3,2,4,-3,4]),
    (10, 43): (5, [-1,-1,2,-1,-3,2,4,-3,4,4]),
    (10, 44): (5, [-1,-1,2,-1,-3,2,-3,4,-3,4]),
    (10, 45): (5, [-1,2,-1,2,-3,2,-3,4,-3,4]),
    (10, 46): (3, [1,1,1,1,1,-2,1,1,1,-2]),
    (10, 47): (3, [1,1,1,1,1,-2,1,1,-2,-2]),
    (10, 48): (3, [-1,-1,-1,-1,2,2,-1,2,2,2]),
    (10, 49): (4, [-1,-1,-1,-1,2,-1,-3,-2,-2,-2,-3]),
    (10, 50): (4, [1,1,2,-1,2,2,-3,2,2,2,-3]),
    (10, 51): (4, [1,1,2,-1,2,2,-3,2,2,-3,-3]),
    (10, 52): (4, [1,1,1,-2,1,1,-2,-2,-3,2,-3]),
    (10, 53): (5, [-1,-1,-2,1,-2,3,-2,-4,-3,-3,-3,-4]),
    (10, 54): (4, [1,1,1,-2,1,1,-2,-3,2,-3,-3]),
    (10, 55): (5, [-1,-1,-1,-2,1,3,-2,-4,-3,-3,-3,-4]),
    (10, 56): (4, [1,1,1,2,-1,2,-3,2,2,2,-3]),
    (10, 57): (4, [1,1,1,2,-1,2,-3,2,2,-3,-3]),
    (10, 58): (6, [1,-2,1,3,-2,-4,-3,-3,5,-4,5]),
    (10, 59): (5, [-1,2,-1,2,-3,2,2,4,-3,4]),
    (10, 60): (5, [-1,2,-1,2,2,-3,2,-3,-2,-4,3,-4]),
    (10, 61): (4, [1,1,1,-2,1,1,1,-2,-3,2,-3]),
    (10, 62): (3, [1,1,1,1,-2,1,1,1,-2,-2]),
    (10, 63): (5, [-1,-1,2,-1,-3,-2,-2,-2,-3,-4,3,-4]),
    (10, 64): (3, [1,1,1,-2,1,1,1,-2,-2,-2]),
    (10, 65): (4, [1,1,2,-1,2,-3,2,2,2,-3,-3]),
    (10, 66): (4, [-1,-1,-1,2,-1,-3,-2,-2,-2,-3,-3]),
    (10, 67): (5, [-1,-1,-1,-2,1,-2,-3,2,2,4,-3,-2,4,-3]),
    (10, 68): (5, [1,1,-2,1,-2,-2,-3,2,2,-4,3,-2,-4,-3]),
    (10, 69): (5, [1,1,2,-1,-3,2,1,4,-3,2,-3,4]),
    (10, 70): (5, [-1,2,-1,-3,2,2,2,4,-3,4]),
    (10, 71): (5, [-1,-1,2,-1,-3,2,2,4,-3,4]),
    (10, 72): (4, [1,1,1,1,2,2,-1,2,-3,2,-3]),
    (10, 73): (5, [-1,-1,-2,1,-2,-1,3,-2,3,-4,3,-4]),
    (10, 74): (5, [-1,-1,-2,1,-2,-2,-3,2,2,4,-3,-2,4,-3]),
    (10, 75): (5, [1,-2,1,-2,3,-2,-2,4,-3,2,4,3]),
    (10, 76): (4, [1,1,1,1,2,-1,-3,2,2,2,-3]),
    (10, 77): (4, [1,1,1,1,2,-1,-3,2,2,-3,-3]),
    (10, 78): (5, [-1,-1,-2,1,-2,-1,3,-2,-4,3,-4,-4]),
    (10, 79): (3, [-1,-1,-1,2,2,-1,-1,2,2,2]),
    (10, 80): (4, [-1,-1,-1,2,-1,-1,-3,-2,-2,-2,-3]),
    (10, 81): (5, [1,1,-2,1,3,2,2,-4,-3,-3,-3,-4]),
    (10, 82): (3, [-1,-1,-1,-1,2,-1,2,-1,2,2]),
    (10, 83): (4, [1,1,2,-1,2,-3,2,2,-3,2,-3]),
    (10, 84): (4, [1,1,1,2,-1,-3,2,2,-3,2,-3]),
    (10, 85): (3, [-1,-1,-1,-1,2,-1,-1,2,-1,2]),
    (10, 86): (4, [-1,-1,2,-1,2,-1,2,2,3,-2,3]),
    (10, 87): (4, [1,1,1,2,-1,-3,2,-3,2,-3,-3]),
    (10, 88): (5, [-1,2,-1,-3,2,-3,2,4,-3,4]),
    (10, 89): (5, [-1,2,-1,2,3,-2,-1,-4,-3,2,-3,-4]),
    (10, 90): (4, [-1,-1,2,-1,2,3,-2,-1,3,2,2]),
    (10, 91): (3, [-1,-1,-1,2,-1,2,2,-1,2,2]),
    (10, 92): (4, [1,1,1,2,2,-3,2,-1,2,-3,2]),
    (10, 93): (4, [-1,-1,2,-1,-1,2,-1,2,3,-2,3]),
    (10, 94): (3, [1,1,1,-2,1,1,-2,-2,1,-2]),
    (10, 95): (4, [-1,-1,2,2,-3,2,-1,2,3,3,2]),
    (10, 96): (5, [-1,2,1,-3,2,1,-3,4,-3,2,-3,4]),
    (10, 97): (5, [1,1,2,-1,2,1,-3,2,-1,2,3,-4,3,-4]),
    (10, 98): (4, [-1,-1,-2,-2,3,-2,1,-2,-2,3,-2]),
    (10, 99): (3, [-1,-1,2,-1,-1,2,2,-1,2,2]),
    (10, 100): (3, [-1,-1,-1,2,-1,-1,2,-1,-1,2]),
    (10, 101): (5, [1,1,1,2,-1,3,-2,1,3,2,2,4,-3,4]),
    (10, 102): (4, [-1,-1,2,-1,-3,2,-1,2,2,3,3]),
    (10, 103): (4, [-1,-1,-2,1,3,-2,-2,3,-2,-2,3]),
    (10, 104): (3, [-1,-1,-1,2,2,-1,2,-1,2,2]),
    (10, 105): (5, [1,1,-2,1,3,2,2,-4,-3,2,-3,-4]),
    (10, 106): (3, [1,1,1,-2,1,-2,1,1,-2,-2]),
    (10, 107): (5, [-1,-1,2,-1,3,2,2,-4,3,-2,3,-4]),
    (10, 108): (4, [1,1,-2,1,1,3,-2,1,-2,-3,-3]),
    (10, 109): (3, [-1,-1,2,-1,2,2,-1,-1,2,2]),
    (10, 110): (5, [-1,2,-1,-3,-2,-2,-2,4,3,-2,3,4]),
    (10, 111): (4, [1,1,2,2,-3,2,2,-1,2,-3,2]),
    (10, 112): (3, [-1,-1,-1,2,-1,2,-1,2,-1,2]),
    (10, 113): (4, [1,1,1,2,-3,2,-1,2,-3,2,-3]),
    (10, 114): (4, [-1,-1,-2,1,3,-2,3,-2,3,-2,3]),
    (10, 115): (5, [1,-2,1,3,2,2,-4,-3,2,-3,-3,-4]),
    (10, 116): (3, [-1,-1,2,-1,-1,2,-1,2,-1,2]),
    (10, 117): (4, [1,1,2,2,-3,2,-1,2,-3,2,-3]),
    (10, 118): (3, [1,1,-2,1,-2,1,-2,-2,1,-2]),
    (10, 119): (4, [-1,-1,2,-1,-3,2,-1,2,3,3,2]),
    (10, 120): (5, [-1,-1,-2,1,3,2,-1,-4,-3,-2,-2,-3,-3,-4]),
    (10, 121): (4, [-1,-1,-2,3,-2,1,-2,3,-2,3,-2]),
    (10, 122): (4, [1,1,2,-3,2,-1,-3,2,-3,2,-3]),
    (10, 123): (3, [-1,2,-1,2,-1,2,-1,2,-1,2]),
    (10, 124): (3, [1,1,1,1,1,2,1,1,1,2]),
    (10, 125): (3, [1,1,1,1,1,-2,-1,-1,-1,-2]),
    (10, 126): (3, [-1,-1,-1,-1,-1,-2,1,1,1,-2]),
    (10, 127): (3, [-1,-1,-1,-1,-1,-2,1,1,-2,-2]),
    (10, 128): (4, [1,1,1,2,1,1,2,2,3,-2,3]),
    (10, 129): (4, [1,1,1,-2,-1,-1,3,-2,-1,3,-2]),
    (10, 130): (4, [1,1,1,-2,-1,-1,-2,-2,-3,2,-3]),
    (10, 131): (4, [-1,-1,-1,-2,1,1,-2,-2,-3,2,-3]),
    (10, 132): (4, [1,1,1,-2,-1,-1,-2,-3,2,-3,-3]),
    (10, 133): (4, [-1,-1,-1,-2,1,1,-2,-3,2,-3,-3]),
    (10, 134): (4, [1,1,1,2,1,1,2,3,-2,3,3]),
    (10, 135): (4, [1,1,1,2,-1,2,-3,-2,-2,-2,-3]),
    (10, 136): (5, [1,-2,1,-2,-3,2,2,4,-3,4]),
    (10, 137): (5, [-1,2,-1,2,-3,-2,-2,4,-3,4]),
    (10, 138): (5, [-1,2,-1,2,3,2,2,-4,3,-4]),
    (10, 139): (3, [1,1,1,1,2,1,1,1,2,2]),
    (10, 140): (4, [1,1,1,-2,-1,-1,-1,-2,-3,2,-3]),
    (10, 141): (3, [1,1,1,1,-2,-1,-1,-1,-2,-2]),
    (10, 142): (4, [1,1,1,2,1,1,1,2,3,-2,3]),
    (10, 143): (3, [-1,-1,-1,-1,-2,1,1,1,-2,-2]),
    (10, 144): (4, [-1,-1,-2,1,-2,-1,3,-2,-1,3,2]),
    (10, 145): (4, [-1,-1,-2,1,-2,-1,-3,-2,1,-2,-3]),
    (10, 146): (4, [-1,-1,2,-1,2,1,-3,2,-1,2,-3]),
    (10, 147): (4, [1,1,1,-2,1,-2,-3,2,-1,2,-3]),
    (10, 148): (3, [-1,-1,-1,-1,-2,1,1,-2,1,-2]),
    (10, 149): (3, [-1,-1,-1,-1,-2,1,-2,1,-2,-2]),
    (10, 150): (4, [1,1,1,-2,1,1,3,-2,-1,3,2]),
    (10, 151): (4, [1,1,1,2,-1,-1,3,-2,1,3,-2]),
    (10, 152): (3, [-1,-1,-1,-2,-2,-1,-1,-2,-2,-2]),
    (10, 153): (4, [-1,-1,-1,-2,-1,-1,3,2,2,2,3]),
    (10, 154): (4, [1,1,2,-1,2,1,3,2,2,2,3]),
    (10, 155): (3, [1,1,1,2,-1,-1,2,-1,-1,2]),
    (10, 156): (4, [-1,-1,-1,2,1,1,-3,-2,1,-2,-3]),
    (10, 157): (3, [1,1,1,2,2,-1,2,-1,2,2]),
    (10, 158): (4, [-1,-1,-1,-2,1,1,3,2,-1,2,3]),
    (10, 159): (3, [-1,-1,-1,-2,1,-2,1,1,-2,-2]),
    (10, 160): (4, [1,1,1,2,1,1,-3,2,-1,2,-3]),
    (10, 161): (3, [-1,-1,-1,-2,1,-2,-1,-1,-2,-2]),
    (10, 162): (4, [-1,-1,-2,1,1,-2,-2,-1,3,-2,3]),
    (10, 163): (4, [1,1,-2,-1,-1,3,2,-1,2,2,3]),
    (10, 164): (4, [1,1,-2,1,-2,-2,-3,2,-1,2,-3]),
    (10, 165): (4, [1,1,2,-1,-3,2,-1,2,3,3,2])
}
