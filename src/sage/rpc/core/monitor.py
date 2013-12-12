"""
The Monitor Process

The monitor sits between the client (your code trying to perform a
remote computation) and the sage compute server. The monitor is a
separate process so that it does not block while a computation is
running::

    MonitorClient <-> [MonitorServer, ComputeClient] <-> Sage Compute Server
"""

import os
import sys

from sage.rpc.core.transport import Transport, TransportListen
from sage.rpc.core.client_base import ClientBase
from sage.rpc.core.server_base import ServerBase



class MonitorClient(ClientBase):
    
    def api_version(self):
        return 'monitor v1'

    def construct_rpc_table(self):
        rpc = super(MonitorClient, self).construct_rpc_table()
        rpc['sage_eval.result'] = self._impl_sage_eval_result
        return rpc        

    def __init__(self, transport, cookie):
        """
        This is the end-user visible client.

        Your code should use it to make queries to the compute server
        via the monitor, and in return receive callbacks.
        """
        super(MonitorClient, self).__init__(transport, cookie)
        self.log.info('MonitorClient started')

    def sage_eval(self, code_string):
        self.rpc.sage_eval.start(code_string)

    def _impl_sage_eval_result(self, cpu_time, wall_time):
        """
        RPC callback when evaluation is finished
        """
        print('Evaluation finished in cpu={0}ms, wall={1}ms'
              .format(int(1000*cpu_time), int(1000*wall_time)))


class MonitorServer(ServerBase):
    
    def api_version(self):
        return 'monitor v1'

    def construct_rpc_table(self):
        rpc = super(MonitorServer, self).construct_rpc_table()
        rpc['sage_eval.start'] = self.monitor._impl_sage_eval_start
        rpc['util.ping'] = self.monitor._impl_ping
        rpc['util.quit'] = self.monitor._impl_quit
        return rpc        

    def __init__(self, monitor, transport, cookie):
        """
        The server half of the monitor
        
        This part connects to the :class:`MonitorClient`.
        """
        self.monitor = monitor
        super(MonitorServer, self).__init__(transport, cookie)
        self.log.info('MonitorServer started')


class ComputeClient(ClientBase):
    
    def api_version(self):
        return 'compute v1'

    def construct_rpc_table(self):
        rpc = super(ComputeClient, self).construct_rpc_table()
        rpc['sage_eval.finished'] = self.monitor._impl_sage_eval_finished
        rpc['util.pong'] = self.monitor._impl_pong
        return rpc

    def __init__(self, monitor, transport, cookie):
        """
        The client half of the monitor.

        This part connects to the compute server.
        """
        self.monitor = monitor
        super(ComputeClient, self).__init__(transport, cookie)
        from sage.rpc.core.server_base import RemoteProcedureLogger
        self.log = RemoteProcedureLogger(monitor.server, origin='compute')
        self.log.info('ComputeClient started')

 
class Monitor(object):

    def __init__(self, process, server_transport, client_transport, cookie):
        """
        The monitor
        """
        self.server = MonitorServer(self, server_transport, cookie)
        self.client = ComputeClient(self, client_transport, cookie)
        self.process = process
        self.log = self.server.log
        
    def _impl_sage_eval_start(self, code_string):
        self.client.rpc.compute.sage_eval(code_string)

    def _impl_sage_eval_finished(self, cpu_time, wall_time):
        self.server.rpc.sage_eval.result(cpu_time, wall_time)

    def loop(self):
        from sage.rpc.core.transport import TransportError
        try:
            srl, swl, sxl = self.server.select_args()
            crl, cwl, cxl = self.client.select_args()
        except TransportError:
            self.quit()
            sys.exit(0)
        rlist = srl + crl
        wlist = swl + cwl
        xlist = sxl + cxl
        import select
        rlist, wlist, xlist = select.select(rlist, wlist, xlist)
        if rlist == [] and wlist == [] and xlist == []:
            return True   # timeout
        self.client.select_handle(rlist, wlist, xlist)
        self.server.select_handle(rlist, wlist, xlist)

    def _impl_ping(self, count, time):
        self.log.debug('ping #%s', count)
        self.client.rpc.util.ping(count, time)

    def _impl_pong(self, count, time):
        self.log.debug('pong #%s', count)
        self.server.rpc.util.pong(count, time)

    def _impl_quit(self):
        self.quit()

    def quit(self):
        print('monitor quit')
        self.log.debug('quit')
        self.client.rpc.util.quit()
        self.client._transport.flush()
        self.process.wait()
        self.client.close()

        # todo: send a notification to the Client that we are about to close
        self.server.close()
        self.server._transport.flush()
        sys.exit(0)
        



class MonitoredProcess(object):

    def __init__(self, argv):
        """
        A monitored external process
        
        INPUT:

        - ``argv`` -- list of string.
        """
        self._init_listen()
        self._init_process(argv)
        self.transport().accept()

    def _init_listen(self):
        client_uri = 'tcp://localhost:0'
        self._transport = TransportListen(client_uri)

    def _init_process(self, argv):
        t = self.transport()
        self._cmd = [cmd.format(port=t.port(), interface=t.interface())
                     for cmd in argv]
        from subprocess import Popen, PIPE
        self._proc = Popen(self._cmd, stdout=sys.stdout, stderr=sys.stdout)

    def transport(self):
        return self._transport

    def port(self):
        return self.transport().port()

    def close(self):
        self._proc.kill()
        
    def wait(self):
        self._proc.wait()


def start_monitor(port, interface):
    print 'starting monitor process'

    # set up the transport 1: connecting to the client
    uri = 'tcp://{0}:{1}'.format(interface, port)
    transport = Transport(uri)
    transport.connect()

    # set up transport part 2: listen for and launch the compute server
    cmd = ['sage', '-c', 
           'from sage.rpc.compute_server import start_server; '
           'start_server({port}, "{interface}")']
    process = MonitoredProcess(cmd)

    # start up
    cookie = os.environ['COOKIE']
    monitor = Monitor(process, transport, process.transport(), cookie)
    while True:
        # print 'monitor loop'
        monitor.loop()
