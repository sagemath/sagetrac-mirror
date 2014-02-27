"""
The Sage Compute Server

This is a separate process that evaluates code and writes the result
to stdout/stderr. To be useful, you have to capture the output
somehow. Usually, it is running under the control of
:mod:`~sage.rpc.core.monitor`.

EXAMPLES::

    sage: from sage.rpc.factory import sage_remote
    sage: c = sage_remote.new_client().wait_for_initialization()
    sage: c.sage_eval('1+1')
    sage: c.loop()
    sage: c.loop()
    sage: c.loop()
"""

from __future__ import print_function

import os
import sys
import time
from sage.misc.cstdio_buffer import linebuffer_output
from sage.rpc.core.server_base import ServerBase
from sage.rpc.core.transport import TransportError



def start_server(port, interface):
    """
    Start the compute server

    INPUT:

    - ``port`` -- integer. The port to connect to.

    - ``interface`` -- string. The interface to listen on.

    OUTPUT:

    This function never returns. Instead it keeps waiting for new
    compute tasks, and processing them as they come in.

    EXAMPLES::

        sage: from sage.rpc.compute_server import start_server
    """
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
        rpc['print'] = self.rpc_print
        rpc['print_end_marker'] = self.rpc_print_end_marker
        rpc['sage_eval'] = self.rpc_sage_eval_init
        rpc['code_completion.start'] = self.rpc_code_completion_start
        return rpc

    def __init__(self, transport, cookie):
        super(SageComputeServer, self).__init__(transport, cookie)
        self._init_shell()

    def _init_shell(self):
        """
        Setup the Sage/IPython shell
        """
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
        linebuffer_output()
        self._shell.run_cell(code_string)
        t1_cpu = time.clock()
        t1_wall = time.time()
        sys.stdout.flush()
        sys.stderr.flush()
        self.rpc.sage_eval.finished(t1_cpu - t0_cpu, t1_wall - t0_wall, label)

    def rpc_print(self, value):
        print(value)
        sys.stdout.flush()        

    def rpc_print_end_marker(self):
        sys.stdout.write(self.end_marker)
        sys.stderr.write(self.end_marker)
        sys.stdout.flush()        
        sys.stderr.flush()

    def rpc_code_completion_start(self, line, position, label):
        basestr, completions = self._shell.complete(None, line, position)
        self.rpc.code_completion.finished(basestr, completions, label)
