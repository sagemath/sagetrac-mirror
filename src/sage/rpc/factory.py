"""
Factory Object to Create Sage Client/Server Instances

This defines a convenience object `sage_remote` in the global
namespace for ease of doctesting.

The simplest example is the ping test::

    sage: client, server = sage_remote.test_instance()
    sage: client.ping()
    sage: sage_remote.test_idle()   # simulate idle loop
    pong #0 (...ms)

The actual rpc call is under the ``rpc`` attribute::

    sage: client.rpc
    Available RPC calls:
    * util.ping
    * util.pong
    * util.quit
    * util.remote_is_ready
    sage: client.rpc.util.ping(int(0), 'count')
    sage: sage_remote.test_idle()
    pong #count (...ms)

Another feature is that server logs are sent and displayed at the
client::

    sage: server.log.error('foobar')
    sage: sage_remote.test_idle()
    [RPC] {server} ERROR: foobar
"""

import os

sage_remote = None


class SageRemoteFactory(object):

    def __init__(self):
        """
        Factory object to create frontend/backend instances.

        EXAMPLES::

            sage: sage_remote
            <sage.rpc.factory.SageRemoteFactory object at 0x...>
        """
        assert sage_remote is None   # singleton class 
        import weakref
        self._idle = weakref.WeakSet()
        from sage.doctest import DOCTEST_MODE
        if not DOCTEST_MODE:
            self.set_debug(True)

    def _add_idle(self, obj):
        self._idle.add(obj)
        from sage.misc.inputhook import sage_inputhook
        sage_inputhook.add(self._inputhook)
    
    def _inputhook(self):
        """
        Main loop idle handler.
        
        Bidirectional RPC is 100% asynchronous, so ideally the
        different endpoints should run in separate processes or
        threads. However, on the Sage command line we have only a
        single thread. Hence we fake the asynchronous operation by
        hooking in the idle loop (when Sage is sitting at the input
        prompt).
        """
        from sage.rpc.core.transport import TransportError
        rlist = []
        wlist = []
        xlist = []
        for obj in list(self._idle):
            try:
                rl, wl, xl = obj.select_args()
            except TransportError:
                self._idle.discard(obj)
                continue
            rlist += rl
            wlist += wl
            xlist += xl
        import select
        rlist, wlist, xlist = select.select(rlist, wlist, xlist, 0)
        if rlist == [] and wlist == [] and xlist == []:
            return True   # timeout
        from sage.libs.readline import interleaved_output
        with interleaved_output():
            for obj in list(self._idle):
                try:
                    obj.select_handle(rlist, wlist, xlist)
                except TransportError:
                    self._idle.discard(obj)
                    obj.close()
        return len(self._idle) > 0   # remove idle handler when empty

    def set_debug(self, value=True):
        """
        Set the logging level.
        
        INPUT:

        - ``value`` -- integer logging level or boolean. If ``True``,
          debug logging will be enabled.
        """
        import logging
        if value is True:
            value = logging.DEBUG
        elif value is False:
            value = logging.WARNING  # default level
        from sage.rpc.core.logging_origin import logger
        logger.setLevel(value)

    def random_cookie(self, length=30):
        """
        Return a new random string.

        OUTPUT:

        A random string.

        EXAMPLES::

            sage: sage_remote.random_cookie()    # random output
            'FayJpeGUjD7wg0tSqQGEpzupkWX1km'
        """
        import string
        import random
        return ''.join(random.choice(string.ascii_letters + string.digits) 
                       for x in range(length))

    def random_port(self):
        """
        Return a random port that is not currently being listened on.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: sage_remote.random_port()
            8214
        """
        from sage.misc.prandom import randint
        return randint(1025, 65536)

    def start_server(self, cookie, port, interface='localhost'):
        """
        Run server in the current process.

        INPUT:

        - ``cookie`` -- string. Magic string that authenticates the
          client connection.
        
        - ``port`` -- integer. TCP port to listen on.
        
        - ``interface`` -- string (default: ``'localhost'``). The
          interface to listen on.

        .. WARNING::

            If ``interface=''`` then remote computers can connect and
            execute arbitrary code. It is your responsibility to
            secure this, for example by firewalling the port or
            running the server in an unprivileged Unix account or
            virtual machine.

        This method only returns after the server closed down, for
        example in response to the client requesting a server
        shutdown.
        """
        from sage.rpc.core.transport import Transport
        uri = 'tcp://{0}:{1}'.format(interface, port)
        transport = Transport(uri)
        transport.connect()
        from sage.rpc.core.server_base import ServerBase
        server = ServerBase(transport, cookie)
        while True:
            server.loop()

    def spawn_server(self, cookie, port, interface):
        """
        Create a server in a new process.

        See :meth:`start_server` for a description of the input.

        OUTPUT:

        The subprocess running the server.
        """
        env = dict(os.environ)
        env['COOKIE'] = cookie
        cmd = ['sage', '-c',
               'sage_remote.start_server({0}, {1}, "{2}")'
               .format('os.environ["COOKIE"]', port, interface)]
        import sys
        from subprocess import Popen
        return Popen(cmd, env=env, stdout=sys.stdout, stderr=sys.stderr)

    def new_client_base(self, interface='localhost'):
        """
        Construct a new client instance

        INPUT:

        - ``cookie`` -- string. Secret that authenticates the
          connection.

        - ``server_uri`` -- string. URI of the backend server. The
          server to connect to must already be running, or else FIXME

        OUTPUT:

        A new :class:`~sage.rpc.core.client_base.ClientBase` instance.
        
        EXAMPLES::

            sage: client = sage_remote.new_client_base()
            sage: client.ping()
            sage: client.loop()
            True
            sage: client.loop()
            pong #0 (...ms)
            True
        """
        cookie = self.random_cookie()
        from sage.rpc.core.transport import TransportListen
        uri = 'tcp://localhost:0'.format(interface)
        transport = TransportListen(uri)
        proc = self.spawn_server(cookie, transport.port(), interface)
        transport.accept()
        from sage.rpc.core.client_base import ClientBase
        client = ClientBase(transport, cookie)
        client.wait_for_initialization()
        self._add_idle(client)
        return client

    def spawn_monitor(self, cookie, port, interface):
        """
        Create a server in a new process.

        See :meth:`start_server` for a description of the input.

        OUTPUT:

        The subprocess running the server.
        """
        env = dict(os.environ)
        env['COOKIE'] = cookie
        cmd = ['sage', '-python', '-c',
               'from sage.rpc.core.monitor import start_monitor; '
               'start_monitor({0}, "{1}")'.format(port, interface)]
        import sys
        from subprocess import Popen
        return Popen(cmd, env=env, stdout=sys.stdout, stderr=sys.stderr)

    def new_client(self, interface='localhost'):
        """
        Construct a new client instance.

        INPUT:

        - ``cookie`` -- string. Secret that authenticates the
          connection.

        - ``server_uri`` -- string. URI of the backend server. The
          server to connect to must already be running, or else FIXME

        OUTPUT:

        A new :class:`~sage.rpc.core.monitor.MonitorClient` instance.

        EXAMPLES::

            sage: client = sage_remote.new_client()
            sage: client.ping()
            sage: client.wait('util.pong')
            pong #0 (...ms)
            <sage.rpc.core.monitor.MonitorClient object at 0x...>
        """
        cookie = self.random_cookie()
        from sage.rpc.core.transport import TransportListen
        uri = 'tcp://localhost:0'.format(interface)
        transport = TransportListen(uri)
        proc = self.spawn_monitor(cookie, transport.port(), interface)
        transport.accept()
        from sage.rpc.core.monitor import MonitorClient
        client = MonitorClient(transport, cookie)
        client.wait_for_initialization()
        self._add_idle(client)        
        return client

    def _new_monitor_client(self, monitor_client_class, interface='localhost'):
        """
        Construct a new client instance.

        This method is for internal use only.

        INPUT:

        - ``monitor_client_class`` -- subclass of
          :class:`~sage.rpc.core.monitor.MonitorClient`.

        - ``interface`` -- string, defaults to
          ``'localhost'``. Network interface to listen on.

        OUTPUT:

        A new ``monitor_client_class`` instance.

        EXAMPLES::

            sage: client = sage_remote.new_client()
            sage: client.ping()
            sage: client.wait('util.pong')
            pong #0 (...ms)
            <sage.rpc.core.monitor.MonitorClient object at 0x...>
        """
        cookie = self.random_cookie()
        from sage.rpc.core.transport import TransportListen
        uri = 'tcp://localhost:0'.format(interface)
        transport = TransportListen(uri)
        proc = self.spawn_monitor(cookie, transport.port(), interface)
        transport.accept()
        client = monitor_client_class(transport, cookie)
        client.wait_for_initialization()
        self._add_idle(client)        
        return client

    def new_client(self, interface='localhost'):
        """
        Construct a new client instance.

        INPUT:

        - ``interface`` -- string, defaults to
          ``'localhost'``. Network interface to listen on.

        OUTPUT:

        A new :class:`~sage.rpc.core.monitor.MonitorClient` instance.

        EXAMPLES::

            sage: client = sage_remote.new_client()
            sage: client.ping()
            sage: client.wait('util.pong')
            pong #0 (...ms)
            <sage.rpc.core.monitor.MonitorClient object at 0x...>
        """
        from sage.rpc.core.monitor import MonitorClient
        return self._new_monitor_client(MonitorClient, interface)

    def new_blocking(self, interface='localhost'):
        """
        Construct a new blocking interface to another Sage session.

        INPUT:

        - ``interface`` -- string, defaults to
          ``'localhost'``. Network interface to listen on.

        OUTPUT:

        A new :class:`~sage.rpc.sage_interface.BlockingSage` instance.

        EXAMPLES::

            sage: s = sage_remote.test_blocking()  # indirect doctest
            sage: s.ping()
            pong #0 (...ms)
        """
        from sage.rpc.sage_client import SageClient
        from sage.rpc.sage_interface import BlockingSageInterface
        client = self._new_monitor_client(SageClient, interface)
        return BlockingSageInterface(client)

    def test_blocking(self):
        try:
            return self._test_blocking
        except AttributeError:
            pass
        self._test_blocking = self.new_blocking()
        return self._test_blocking

    def test_instance(self):
        """
        Create a test client/server combination.

        The instance is cached, so subsequent calls return the same
        instance.
        """
        try:
            return self._test_instance
        except AttributeError:
            pass
        # set up the transport 1: listen on port
        from sage.rpc.core.transport import Transport, TransportListen
        client_uri = 'tcp://localhost:0'
        listen = TransportListen(client_uri)
        # set up transport part 2: connect to listening port
        server_uri = 'tcp://{0}:{1}'.format(listen.interface(), listen.port())
        transport = Transport(server_uri)
        transport.connect()
        listen.accept()
        # Create client/server pair
        from sage.rpc.core.client_base import ClientBase
        from sage.rpc.core.server_base import ServerBase
        cookie = self.random_cookie()
        client = ClientBase(listen, cookie)
        server = ServerBase(transport, cookie)
        for i in range(5):
            server.loop(0)
            client.loop(0)
        self._test_instance = (client, server)
        self._add_idle(client)
        self._add_idle(server)
        return self._test_instance

    def test_idle(self):
        """
        """
        client, server = self.test_instance()
        import time
        for i in range(10):
            self._inputhook()
            time.sleep(0.01)

sage_remote = SageRemoteFactory()
