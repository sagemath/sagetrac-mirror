

def is_valid_sha1(sha1):
    """
    Test whether sha1 is a ascii sha1

    EXAMPLES::

        >>> is_valid_sha1('e282e1537ae6614980b3a8a46e2e3a99ecff5096')
        True
    """
    return len(sha1) == 40 and int(sha1, 16) >= 0
