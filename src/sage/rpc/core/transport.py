"""
Transport

This module abstracts the underlying socket IO mechanism. Think either
network or Unix sockets. However, note that the interface is about
sending complete messages ``read``/``write``, not streaming
``send``/``recv``.

Since a single method in socket communication code is generally not
that useful, there are no really meaningful doctests in this
file. That everything is working together is tested in
:mod:`sage.rpc.core.transport_test`.

There are, in fact, two classes. They differ in how they establish a
connection:

* :class:`Transport` will make a connection to an already-listening
  port. You need to call :meth:`~Transport.connect` to establish the
  connection.

* :class:`TransportListen` will bind to a port and listen. You need to
  call :meth:`~TransportListen.accept` to establish the connection. By
  passing port 0, one can bind to a random unused port.

Both will block until the connection is established.
"""

import socket
import select
import re
import errno
import json

DNS_REGEX_PATTERN = '[a-zA-Z0-9][a-zA-Z0-9\-\.]*'
URI_RE = re.compile(r'(?P<scheme>[A-Za-z]+)://(?P<interface>'
                    + DNS_REGEX_PATTERN + 
                    r'):(?P<port>[0-9]+)')
SEPARATOR = '\0'


class TransportError(Exception):
    pass


class TransportNotConnected(TransportError):
    pass


class TransportBase(object):

    BUFFER_SIZE = 4096

    def _init_buffers(self):
        """
        Initialize the internal buffers

        TESTS::

        
        """
        self._read_buf = ''
        self._read_pos = -1
        self._write_buf = ''
    
    def _init_uri(self, uri):
        """
        Parse the uri

        TESTS::

            sage: from sage.rpc.core.transport import URI_RE
            sage: URI_RE.search('tcp://localhost:123').groups()
            ('tcp', 'localhost', '123')
        """
        m = URI_RE.search(uri)
        self._scheme = m.group('scheme').lower()
        self._interface = m.group('interface').lower()
        self._port = int(m.group('port'))

    def interface(self):
        """
        Return the interface name.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.rpc.core.transport import TransportListen
            sage: transport = TransportListen('tcp://localhost:0')
            sage: transport.interface()
            'localhost'
            sage: transport.port()   # random output
            48647
            sage: transport.close()
        """
        return self._interface

    def is_connected(self):
        """
        Return whether the transport is connected to another endpoint.

        OUTPUT:

        Boolean
        
        EXAMPLES::
        
            sage: from sage.rpc.core.transport import Transport
            sage: transport = Transport('tcp://localhost:80')
            sage: transport.is_connected()
            False
        """
        return hasattr(self, '_sock')

    def _socket(self):
        try:
            return self._sock
        except AttributeError:
            raise TransportNotConnected('transport is not connected')

    def port(self):
        """
        Return the port.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: from sage.rpc.core.transport import TransportListen
            sage: transport = TransportListen('tcp://localhost:0')
            sage: transport.interface()
            'localhost'
            sage: transport.port()   # random output
            48647
            sage: transport.close()
        """
        return self._port

    def fileno(self):
        """
        Return the FD of the socket.

        This is important when integrating into external event loops.

        OUTPUT:

        A small integer.

        EXAMPLES::

            sage: from sage.rpc.core.transport import Transport
            sage: transport = Transport('tcp://localhost:0')
            sage: transport.fileno()            
            Traceback (most recent call last):
            ...
            TransportNotConnected: transport is not connected
        """
        try:
            return self._socket().fileno()
        except socket.error as err:
            raise TransportError(str(err))

    def nonblocking_read(self):
        """
        Perform a non-blocking read

        There is no guarantee that anything is read, or that a whole
        message is read.

        This should only be called from an idle loop. If you want to
        write a message to the transport, use :meth:`write`.

        OUTPUT:

        Integer. The number of bytes read.
        """
        try:
            data = self._socket().recv(self.BUFFER_SIZE)
        except socket.error as e:
            if e.errno == errno.EAGAIN or e.errno == errno.EWOULDBLOCK:
                # there happened to be nothing to read, fine
                return
            raise TransportError(str(e))
        if data == '':
            raise TransportError('Remote end closed connection')
        if self._read_pos == -1:
            pos = data.find(SEPARATOR)
            if pos != -1:
                self._read_pos = len(self._read_buf) + pos
        self._read_buf += data
        return len(data)

    def nonblocking_write(self):
        """
        Perform a non-blocking write

        There is no guarantee that anything is written, or that a
        whole message is written.

        This should only be called from an idle loop.
        If you want to read a message from the
        transport, use :meth:`read`.

        OUTPUT:

        Integer. The number of bytes written.
        """
        buf = self._write_buf
        try:
            sent = self._socket().send(buf)
        except socket.error as e:
            if e.errno == errno.EAGAIN or e.errno == errno.EWOULDBLOCK:
                return 0
            raise TransportError(str(e))
        buf = buf[sent:]
        self._write_buf = buf
        return sent
        
    def can_read(self, try_read=True):
        """
        Return whether reading would block
       
        INPUT:

        - ``try_read`` -- boolean (default: ``True``). Whether to
          attempt to do a non-blocking read. Can be turned off for
          debugging purposes.

        OUTPUT:

        Boolean. If ``True``, then :meth:`read` will return a complete
        message.
        """
        if try_read:
            self.nonblocking_read()
        return self._read_pos != -1

    def read(self):
        """
        Read a complete message from the transport

        This method will block until a complete message has been
        received. Query :meth:`can_read` to avoid blocking.

        OUTPUT:

        A dictionary.

        EXAMPLES::

        
        """
        while self._read_pos == -1:
            fd = self.fileno()
            rlist, wlist, xlist = select.select([fd], [], [fd])
            self.nonblocking_read()
        pos = self._read_pos
        buf = self._read_buf
        assert buf[pos] == SEPARATOR
        data = buf[:pos]
        self._read_buf = buf[pos+1:]
        self._read_pos = self._read_buf.find(SEPARATOR)
        return json.loads(data)

    def write(self, data):
        """
        Write to the transport

        This buffers the output, but does not necessarily send it all!
        To ensure that the message has been sent, you must either call
        
        * :meth:`flush`, or

        * :meth:`is_written` until that returns ``True``

        INPUT:
        
        - ``data`` -- string. Must not contain NULL characters.

        OUTPUT:
        
        Boolean. Whether everything was written. If ``False``,
        """
        self._write_buf = self._write_buf + json.dumps(data) + '\0'
        self.nonblocking_write()

    def is_written(self, try_write=True):
        """
        Return whether all messages have been written.

        INPUT:

        - ``try_write`` -- boolean (default: ``True``). Whether to
          attempt to do a non-blocking write. Can be turned off for
          debugging purposes.
        """
        if try_write and len(self._write_buf) > 0:
            self.nonblocking_write()
        return len(self._write_buf) == 0

    def flush(self):
        while not self.is_written():
            select.select([], [self], [self])

    def close(self):
        if self.is_connected():
            self.flush()
            s = self._socket()
            try:
                s.shutdown(socket.SHUT_RDWR)
                s.close()
            except socket.error as e:
                if e.errno == errno.EBADF:
                    pass   # the remote side called close()
                else:
                    raise TransportError(str(e))



class Transport(TransportBase):

    def __init__(self, uri):
        """
        Transport

        This class implements the client part, it must be able to
        connect to the server or the constructor will raise a
        :class:`TransportError` exception.

        INPUT:

        - ``uri`` -- string. A valid URI, including port number. Only
          the ``'tcp://'`` scheme is currenly supported.

        EXAMPLES::
        
            sage: from sage.rpc.core.transport import Transport
            sage: t = Transport('tcp://www.sagemath.org:80');  t
            <sage.rpc.core.transport.Transport object at 0x...>
            sage: t.connect()    # optional - internet
        """
        self._init_uri(uri)
        self._init_buffers()

    def connect(self):
        """
        Connect to the listening port.

        EXAMPLES::
        
            sage: from sage.rpc.core.transport import Transport
            sage: t = Transport('tcp://www.sagemath.org:80');  t
            <sage.rpc.core.transport.Transport object at 0x...>
            sage: t.is_connected()
            False
        """
        if self._scheme == 'tcp':
            s = socket.socket()
            connect_arg = (self._interface, self._port)
        else:
            raise NotImplementedError()
        try:
            s.connect(connect_arg)
        except socket.error as err:
            s.close()
            raise TransportError(err)
        s.setblocking(0)
        self._sock = s
        self._init_buffers()



class TransportListen(TransportBase):

    def __init__(self, uri):
        """
        Bind a port and listen for the connection

        INPUT:
        
        See :meth:`Transport.__init__`.
        """
        self._init_uri(uri)
        self._init_buffers()
        listen = socket.socket()
        listen.bind((self._interface, self._port))
        if self._port == 0:
            self._port = listen.getsockname()[1]
        listen.listen(1)
        self._listen = listen

    def accept(self):
        s, self._remote = self._listen.accept()
        s.setblocking(0)
        self._sock = s
        self._close_listen()

    def _close_listen(self):
        if hasattr(self, '_listen'):
            self._listen.shutdown(socket.SHUT_RDWR)
            self._listen.close()
            del self._listen
        
    def remote_address(self):
        return self._remote

    def close(self):
        self._close_listen()
        super(TransportListen, self).close()
