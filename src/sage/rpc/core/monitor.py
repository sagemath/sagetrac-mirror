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

from sage.rpc.core.transport import Transport, TransportListen, TransportError
from sage.rpc.core.client_base import ClientBase
from sage.rpc.core.server_base import ServerBase



class MonitorClient(ClientBase):
    
    def api_version(self):
        return 'monitor v1'

    def construct_rpc_table(self):
        rpc = super(MonitorClient, self).construct_rpc_table()
        rpc['sage_eval.result'] = self._impl_sage_eval_result
        rpc['sage_eval.stdin'] = self._impl_sage_eval_stdin
        rpc['sage_eval.stdout'] = self._impl_sage_eval_stdout
        rpc['sage_eval.stderr'] = self._impl_sage_eval_stderr
        rpc['sage_eval.crash'] = self._impl_sage_eval_crash
        return rpc        

    def __init__(self, transport, cookie):
        """
        This is the end-user visible client.

        Your code should use it to make queries to the compute server
        via the monitor, and in return receive callbacks.
        """
        super(MonitorClient, self).__init__(transport, cookie)
        self._label_counter = 0
        self.log.info('MonitorClient started')

    def random_label(self):
        """
        Return a new random label.

        OUTPUT:

        A JSON-serializable value.
        """
        self._label_counter += 1
        return self._label_counter

    def sage_eval(self, code_string, label=None):
        """
        Initiate the evaluation of code
        
        INPUT:

        - ``code_string`` -- string. The code to evaluate.

        - ``label`` -- anything JSON-serializable. A label to attach
          to the compute request. Can be ``None``, in which case one
          is randomly chosen for you. See also :meth:`random_label`.

        OUTPUT:
    
        The label of the compute request. Equals ``label`` if one was
        specified.
        """
        if label is None:
            label = self.random_label()
        self.rpc.sage_eval.start(code_string, label)
        return label

    def _impl_sage_eval_result(self, cpu_time, wall_time, label):
        """
        RPC callback when evaluation is finished
        """
        print('Evaluation finished in cpu={0}ms, wall={1}ms'
              .format(int(1000*cpu_time), int(1000*wall_time)))

    def _impl_sage_eval_stdin(self, label):
        """
        RPC callback when evaluation requests stdin
        """
        self.log.debug('stdin')
        print('Input requested')

    def _impl_sage_eval_stdout(self, stdout, label):
        """
        RPC callback when evaluation produces stdout
        """
        self.log.debug('stdout %s', stdout.strip())
        print('STDOUT ' + stdout)

    def _impl_sage_eval_stderr(self, stderr, label):
        """
        RPC callback when evaluation produces stderr
        """
        self.log.debug('stderr %s', stderr.strip())
        print('STDERR ' + stderr)

    def _impl_sage_eval_crash(self):
        """
        RPC callback when the compute server crashed
        """
        print('Compute server crashed.')


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
        self.current_label = None
        
    def _impl_sage_eval_start(self, code_string, label):
        self.current_label = label
        self.client.rpc.compute.sage_eval(code_string, label)

    def _impl_sage_eval_finished(self, cpu_time, wall_time, label):
        self.server.rpc.sage_eval.result(cpu_time, wall_time, label)
        self.current_label = None

    def loop(self):
        try:
            crl, cwl, cxl = self.client.select_args()
        except TransportError:
            self.handle_client_disconnect()
        try:
            srl, swl, sxl = self.server.select_args()
        except TransportError:
            self.handle_server_disconnect()
        prl, pwl, pxl = self.process.select_args()
        rlist = srl + crl + prl
        wlist = swl + cwl + pwl
        xlist = sxl + cxl + pxl
        import select
        rlist, wlist, xlist = select.select(rlist, wlist, xlist)
        if rlist == [] and wlist == [] and xlist == []:
            return True   # timeout
        try:
            self.client.select_handle(rlist, wlist, xlist)
        except TransportError:
            self.handle_client_disconnect()
        try:
            self.server.select_handle(rlist, wlist, xlist)
        except TransportError:
            self.handle_server_disconnect()
        self.process.select_handle(self, rlist, wlist, xlist)

    def handle_server_disconnect(self):
        """
        React to the server disconnecting from the user.
        
        Compare with :meth:`handle_client_disconnect`.
        """
        print('Monitor: disconnected, killing subprocess and exiting.')
        try:
            self.client.rpc.util.quit()
            self.client._transport.flush()
            self.process.wait()
            self.client.close()
        except TransportError:
            pass
        sys.exit(0)                    

    def handle_client_disconnect(self):
        """
        React to the client disconnecting from the compute server

        Compare with :meth:`handle_server_disconnect`.
        """
        print('Monitor: Compute server died.')
        self.server.rpc.sage_eval.crash()
        self.server._transport.flush()
        # todo: restart?
        self.process.wait()
        sys.exit(0)                    
        
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
        
    def need_stdin(self):
        """
        Called by the monitored process if it requests user input (overriding 
        """
        print('need stdin')
        self.server.rpc.sage_eval.stdin(self.current_label)

    def got_stdout(self, stdout):
        """
        Called by the monitored process if it produces stderr
        """
        self.server.rpc.sage_eval.stdout(stdout, self.current_label)

    def got_stderr(self, stderr):
        """
        Called by the monitored process if it produces stderr
        """
        self.server.rpc.sage_eval.stdin(stderr, self.current_label)


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
        self._proc = Popen(self._cmd, stdin=None, stdout=PIPE, stderr=PIPE)
        import fcntl
        for fd in [self._proc.stdout.fileno(), self._proc.stderr.fileno()]:
            fl = fcntl.fcntl(fd, fcntl.F_GETFL)
            fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)

    def select_args(self):
        proc = self._proc
        fd_out = proc.stdout.fileno()
        fd_err = proc.stderr.fileno()
        return ([fd_out, fd_err], [], [fd_out, fd_err])

    def select_handle(self, monitor, rlist, wlist, xlist):
        stdout, stderr = self._proc.stdout, self._proc.stderr
        fd_out = stdout.fileno()
        fd_err = stderr.fileno()
        if fd_out in rlist:
            monitor.got_stdout(stdout.read())
        if fd_err in rlist:
            monitor.got_stderr(stderr.read())

    def transport(self):
        return self._transport

    def port(self):
        return self.transport().port()

    def close(self):
        self._proc.kill()
        
    def wait(self):
        self._proc.wait()


def start_monitor(port, interface):
    print('starting monitor process')

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
        # print('monitor loop')
        monitor.loop()
