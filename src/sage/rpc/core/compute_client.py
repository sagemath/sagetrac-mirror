"""
Client for the Sage Compute Service

The corresponding server implementation is at
:mod:`sage.rpc.compute_sever'. The client/server pair negotiate a
random string that can be used to mark the end of output, but the
:class:`ComputeClient` does not capture the server's output (and, in
particular, has no mechanism to check that this end marker appers in
the server stdout). See :mod:`sage.rpc.core.monitor.Monitor` for an
example that ties everything together.
"""

from sage.rpc.core.client_base import ClientBase



class ComputeClient(ClientBase):
    
    def api_version(self):
        return 'compute v1'

    def construct_rpc_table(self):
        rpc = super(ComputeClient, self).construct_rpc_table()
        rpc['sage_eval.finished'] = self._impl_sage_eval_finished
        return rpc

    def __init__(self, transport, cookie):
        """
        The client half of the compute server.
        """
        self._init_end_marker()
        super(ComputeClient, self).__init__(transport, cookie)

    def _init_end_marker(self):
        """
        To mark the end of stdout/stderr streams on the compute server

        Note that the socket and stdout streams are independent and
        can arrive in random order.
        """
        import string, random
        self._end_marker_len = 78
        self._end_marker = ''.join(random.choice(string.ascii_letters + string.digits)
                                   for x in range(self._end_marker_len))

    def init_send_extra(self):
        """
        Send additional data during initialization.

        OUTPUT:

        Anything that is a JSON serializeable dictionary value.
        """
        return {'end_marker': self._end_marker}

    @property
    def end_marker(self):
        return self._end_marker

    def _impl_sage_eval_finished(self, cpu_time, wall_time, label):
        raise NotImplementedError


