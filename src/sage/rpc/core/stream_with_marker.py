"""
Stream With Embedded End Marker

EXAMPLES::

    sage: from sage.rpc.core.stream_with_marker import StreamWithMarker
    sage: s = StreamWithMarker()
    sage: s.write('foo')
    sage: s.read()
    'foo'

Now we set an end marker into the stream::

    sage: s.stop_at('bar')
    sage: s.write('123b')
    sage: s.read()
    '12'
    sage: s.write('artender')
    sage: s.is_stopped()
    False
    sage: s.read()
    '3'
    sage: s.is_stopped()
    True

One way to continue is to remove the end marker. The output continues
past the end marker::

    sage: s.stop_at(None)
    sage: s.read()
    'tender'

It is possible to have empty reads under certain circumstances::

    sage: s.stop_at('bar')
    sage: s.write('123ba')
    sage: s.read()
    '123'
    sage: s.write('rtender')
    sage: s.is_stopped()
    False
    sage: s.read()
    ''
    sage: s.is_stopped()
    True
    sage: s.read()
    Traceback (most recent call last):
    ...
    StreamStoppedException: end marker encountered
"""

class StreamStoppedException(Exception):
    pass


class StreamWithMarker(object):

    def __init__(self):
        """
        Stream with embedded end marker.

        If the end marker is found, the stream stops (but skips over
        the end marker). All further reads raise an exception until
        the stream is restarted.
        """
        self._buf = ''
        self._stopped = False
        self._end_marker = None

    def write(self, data):
        self._buf += data

    def read(self):
        """
        Read from the stream.

        The end marker is not returned if it is present in the input
        stream!
        
        A :class:`StreamStoppedException` is raised if the stream is
        currently stopped, see :meth:`is_stopped`.
        """
        if self._stopped:
            raise StreamStoppedException('end marker encountered')
        buf = self._buf
        marker = self._end_marker
        if marker is None:   # not currently searching
            self._buf = ''
            return buf
        pos = buf.find(marker)
        if pos >= 0:         # found
            self._stopped = True
            self._buf = buf[pos+len(marker):]
            return buf[0:pos]
        else:                # not found, keep len(marker)-1 chars in buffer for later
            last_safe_pos = len(buf) - len(marker) + 1
            if last_safe_pos <= 0:
                return ''
            else:
                self._buf = buf[last_safe_pos:]
                return buf[0:last_safe_pos]
            
    def stop_at(self, end_marker=None):
        """
        Process until an ``end_marker`` is encountered.

        This also continues the stream if it is currently stopped.

        INPUT:

        - ``end_marker`` -- string or ``None``. End marker to search
          for and then stop output. If ``None``, no end marker is
          searched for.
        """
        self._end_marker = end_marker
        self._stopped = False

    def is_stopped(self):
        """
        Return whether the stream is currently stopped.

        OUTPUT:

        Boolean. If ``True``, it is guaranteed that the next
        :meth:`read` will not raise a :class:`StreamStoppedException`.
        """
        return self._stopped

    def restart(self):
        """
        Restart a stopped stream.
        """
        self._stopped = False
