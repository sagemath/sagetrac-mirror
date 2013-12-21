"""
The Sage Compute Server

This is a separate process that evaluates code and writes the result
to stdout/stderr. To be useful, you have to capture the output
somehow. Usually, it is running under the control of
:mod:`~sage.rpc.core.monitor`.

EXAMPLES::

    
"""

import os
import sys
import time
from sage.rpc.core.server_base import ServerBase
from sage.rpc.core.transport import TransportError


def start_server(port, interface):
    cookie = os.environ['COOKIE']
    from sage.rpc.core.transport import Transport
    uri = 'tcp://{0}:{1}'.format(interface, port)
    transport = Transport(uri)
    transport.connect()
    server = SageComputeServer(transport, cookie)
    while True:
        # print 'compute loop'
        server.loop()


class SageComputeServer(ServerBase):

    def api_version(self):
        return 'compute v1'

    def construct_rpc_table(self):
        rpc = super(ServerBase, self).construct_rpc_table()
        rpc['compute.sage_eval'] = self.rpc_sage_eval_init
        return rpc

    def __init__(self, transport, cookie):
        super(SageComputeServer, self).__init__(transport, cookie)
        self._init_shell()

    def _init_shell(self):
        from sage.misc.interpreter import get_test_shell
        self._shell = get_test_shell()

    def init_receive_extra(self, extra):
        """
        Receive additional data during initialization.

        Override this method for an easy way to get additional
        initialization data. It is not permitted to make RPC calls at
        this point in the negotiation.

        INPUT:

        - ``extra`` -- the output of :meth:`client_base.send_extra` on
          the other end of the connection.
        """
        self.end_marker = extra['end_marker']
        
    def rpc_sage_eval_init(self, code_string, label):
        t0_cpu = time.clock()
        t0_wall = time.time()
        self._shell.run_cell(code_string)
        t1_cpu = time.clock()
        t1_wall = time.time()
        sys.stdout.flush()
        sys.stderr.flush()
        # This kills the process... why?
        #os.fsync(sys.stdout.fileno())
        #os.fsync(sys.stderr.fileno())
        self.rpc.sage_eval.finished(t1_cpu - t0_cpu, t1_wall - t0_wall, label)
