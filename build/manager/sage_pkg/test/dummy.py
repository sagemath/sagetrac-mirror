"""
Dummy test
"""

import unittest


class DummyTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testDummy(self):
        self.assertTrue(True)



if __name__ == '__main__':
    unittest.main()
