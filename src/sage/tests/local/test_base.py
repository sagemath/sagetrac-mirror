"""
Base class for tests to be run on all files in ``SAGE_LOCAL``
"""



class TestException(Exception):
    pass


class LocalFileTest(object):

    @classmethod
    def is_magic(cls, head):
        """
        Test whether the test is applicable to a file

        To be overridden in a derived class

        INPUT:

        - ``head`` -- bytes. The first 512 bytes of the file.

        OUTPUT:
        
        Boolean.

        EXAMPLES::

            sage: from sage.tests.local.test_base import LocalFileTest
            sage: LocalFileTest.is_magic('foo')
            Traceback (most recent call last):
            ...
            NotImplementedError: to be implemented in derived class
        """
        raise NotImplementedError('to be implemented in derived class')

    def __init__(self, filename):
        """
        Base class for wrapper around a file in ``SAGE_LOCAL``

        EXAMPLES::

            sage: from sage.tests.local.test_base import LocalFileTest
            sage: LocalFileTest('foo').filename
            'foo'
        """
        self.filename = filename

    def perform_tests(self):
        """
        Perform tests

        To be overridden in a derived class

        RAISES:

        :class:`TestException` if a test fails.
        
        EXAMPLES::

            sage: from sage.tests.local.test_base import LocalFileTest
            sage: LocalFileTest('foo').perform_tests()
            Traceback (most recent call last):
            ...
            NotImplementedError: to be implemented in derived class
        """
        raise NotImplementedError('to be implemented in derived class')
