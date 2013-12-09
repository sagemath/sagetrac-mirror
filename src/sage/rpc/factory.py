"""
Factory Object to Create Front/Backend Instances

This defines a convenience object `sage_remote` in the global
namespace for ease of doctesting.
"""

import os

class SageRemoteFactory(object):
    """
    Factory object to create frontend/backend instances.
    """

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

            sage: server.random_port()
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
        print('server started')
        from sage.rpc.client.transport import Transport
        uri = 'tcp://{0}:{1}'.format(interface, port)
        transport = Transport(uri)
        transport.connect()
        from sage.rpc.base import ServerBase
        server = ServerBase(transport, cookie)
        while True:
            print 'server main loop'
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

    def new_client(self, interface='localhost'):
        """
        Construct a new client instance.

        INPUT:

        - ``cookie`` -- string. Secret that authenticates the
          connection.

        - ``server_uri`` -- string. URI of the backend server. The
          server to connect to must already be running, or else FIXME

        OUTPUT:

        A new :class:`~sage.rpc.client.base.Client` instance.
        """
        cookie = self.random_cookie()
        from sage.rpc.client.transport import TransportListen
        uri = 'tcp://localhost:0'.format(interface)
        transport = TransportListen(uri)
        proc = self.spawn_server(cookie, transport.port(), interface)
        transport.accept()
        from sage.rpc.client.base import ClientBase
        return ClientBase(transport, cookie)

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
        from sage.rpc.client.transport import Transport, TransportListen
        client_uri = 'tcp://localhost:0'
        listen = TransportListen(client_uri)
        # set up transport part 2: connect to listening port
        server_uri = 'tcp://{0}:{1}'.format(listen.interface(), listen.port())
        transport = Transport(server_uri)
        transport.connect()
        listen.accept()
        # Create client/server pair
        from sage.rpc.client.base import ClientBase
        from sage.rpc.base import ServerBase
        cookie = self.random_cookie()
        client = ClientBase(listen, cookie)
        server = ServerBase(transport, cookie)
        for i in range(5):
            server.loop(0)
            client.loop(0)
        self._test_instance = (client, server)
        return self._test_instance

    def _test_idle_loop(self):
        client, server = self._test_instance
        client.loop(0)
        server.loop(0)
        return True
            





sage_remote = SageRemoteFactory()
