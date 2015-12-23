"""
Base class for tests to be run on all files in ``SAGE_LOCAL``
"""



class TestException(Exception):
    """
    A test failed that was supposed to succeed

    This is raised if something is wrong with how we built Sage, and a
    problem with one of Sage's files was detected.
    """
    pass


class FileMagicException(Exception):
    """
    File header does not match the expected file type
    """


class LocalFileTest(object):

    def __init__(self, filename, head):
        """
        Base class for wrapper around a file in ``SAGE_LOCAL``

        INPUT:

        - ``filename`` -- string. The filename.
        
        - ``head`` -- bytes. The first 512 bytes of the file.

        EXAMPLES::

            sage: from sage.tests.local.test_base import LocalFileTest
            sage: LocalFileTest('foo', 'magic').filename
            'foo'
        """
        self.filename = filename
        self.head = head

    def perform_tests(self):
        """
        Perform tests

        To be overridden in a derived class

        RAISES:

        :class:`TestException` if a test fails.
        
        EXAMPLES::

            sage: from sage.tests.local.test_base import LocalFileTest
            sage: LocalFileTest('foo', 'magic').perform_tests()
            Traceback (most recent call last):
            ...
            NotImplementedError: to be implemented in derived class
        """
        raise NotImplementedError('to be implemented in derived class')

    @classmethod
    def load(cls, filename):
        """
        Initialize just from a file name for testing purposes

        INPUT:

        - ``filename`` -- string. The filename.
        
        OUTPUT:

        :class:`LocalFileTest` instance

        EXAMPLES::

            sage: from sage.tests.local.test_base import LocalFileTest
            sage: LocalFileTest.load(os.path.join(SAGE_LOCAL, 'bin', 'patch'))
            <sage.tests.local.test_base.LocalFileTest object at 0x...>
        """
        with open(filename, 'rb') as f:
            head = f.read(512)
        return cls(filename, head)
