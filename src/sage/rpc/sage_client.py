"""
Client for the Compute Server

See :mod:`sage.rpc.sage_interface` for an example
"""


from sage.rpc.core.transport import TransportListen
from sage.rpc.core.monitor import MonitorClient
from sage.rpc.core.decorator import remote_callable




class SageClient(MonitorClient):
    
    def __init__(self, *args, **kwds):
        super(SageClient, self).__init__(*args, **kwds)
        self.interface = None

    @remote_callable('sage_eval.stdin')
    def _impl_sage_eval_stdin(self, label):
        """
        RPC callback when evaluation requests stdin
        """
        self.interface._impl_sage_eval_stdin(label)

    @remote_callable('sage_eval.stdout')
    def _impl_sage_eval_stdout(self, stdout, label):
        """
        RPC callback when evaluation produces stdout
        """
        self.interface._impl_sage_eval_stdout(stdout, label)

    @remote_callable('sage_eval.stderr')
    def _impl_sage_eval_stderr(self, stderr, label):
        """
        RPC callback when evaluation produces stderr
        """
        self.interface._impl_sage_eval_stderr(stderr, label)

    @remote_callable('sage_eval.result')
    def _impl_sage_eval_result(self, cpu_time, wall_time, label):
        """
        RPC callback when evaluation is finished
        """
        self.interface._impl_sage_eval_result(cpu_time, wall_time, label)

    @remote_callable('sage_eval.crash')
    def _impl_sage_eval_crash(self, label):
        """
        RPC callback when the compute server crashed
        """
        self.interface._impl_sage_eval_crash(label)

    @remote_callable('code_completion.finished')
    def _impl_code_completion_finished(self, base, completions, label):
        self.interface._impl_code_completion_finished(base, completions, label)

    @remote_callable('util.pong')
    def _impl_pong(self, count, start_time):
        self.interface._impl_pong(count, start_time)
