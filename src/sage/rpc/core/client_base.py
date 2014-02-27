"""
Base Class for Sage Clients


"""

from sage.rpc.core.common import (
    RemoteProcedureCaller,
    RemoteProcedureException,
    TYPE_INIT_CONNECTION,
    TYPE_INIT_REPLY,
)

from sage.rpc.core.logging_origin import logger_client, logger_server



class ClientException(RemoteProcedureException):
    pass


class ClientBase(RemoteProcedureCaller):
    
    def api_version(self):
        return 'base class v1'
q
    def construct_rpc_table(self):
        rpc = super(ClientBase, self).construct_rpc_table()
        logger = self.log
        rpc['log.debug']     = logger.debug
        rpc['log.info']      = logger.info
        rpc['log.warn']      = logger.warning
        rpc['log.warning']   = logger.warning
        rpc['log.error']     = logger.error
        rpc['log.exception'] = logger.exception
        rpc['log.critical']  = logger.critical
        rpc['log.fatal']     = logger.critical
        return rpc

    def __init__(self, transport, cookie, logger=None, remote_logger=None):
        self.log = logger_client if logger is None else logger
        self.remote_log = logger_server if remote_logger is None else remote_logger
        super(ClientBase, self).__init__(transport, cookie)

    def init_receive_extra(self, extra):
        """
        Receive additional data during initialization.

        Override this method for an easy way to get additional
        initialization data. It is not permitted to make RPC calls at
        this point in the negotiation.

        INPUT:

        - ``extra`` -- the output of :meth:`server_base.send_extra` on
          the other end of the connection.
        """
        pass

    def init_send_extra(self):
        """
        Send additional data during initialization.

        Override this method for an easy way to pass additional
        initialization data. It is not permitted to make RPC calls at
        this point in the negotiation.

        OUTPUT:

        Anything that is a JSON serializeable dictionary value.
        """
        return None

    def handle_init_connection(self, request):
        self.log.debug('initializing connection %s', request)
        if request.get('type') != TYPE_INIT_CONNECTION:
            error = 'invalid message type'
            self.log.critical(error)
            raise ClientException(error)
        if request.get('cookie') != self._cookie:
            error = 'authentication cookie mismatch'
            self.log.critical(error)
            raise ClientException(error)
        if not self.check_api_version(request.get('api_version')):
            error = 'wrong api version'
            self.log.critical(error)
            raise ClientException(error)
        try:
            rpc_table = request['rpc_table']
        except KeyError:
            error = 'no table of rpc calls'
            self.log.critical(error)
            raise ClientException(error)
        self.init_receive_extra(request.get('extra', None))
        self._init_remote_rpc(rpc_table)
        self._transport.write({
            'type': TYPE_INIT_REPLY,
            'rpc_table': list(self._rpc.keys()),
            'ok': True,
            'extra': self.init_send_extra(),
        })
        self.rpc.util.remote_is_ready()   # the first RPC call

