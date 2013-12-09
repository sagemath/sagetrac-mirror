"""
Base class for Sage Server
"""

from sage.rpc.client.transport import Transport
from sage.rpc.client.common import (
    RemoteProcedureCaller,
    RemoteProcedureException,
    TYPE_INIT_CONNECTION,
    TYPE_INIT_REPLY,
    log_client, log_server,
)


class ServerException(RemoteProcedureException):
    pass



class ServerBase(RemoteProcedureCaller):
    
    _api_version = 'abstract_base_class'

    def __init__(self, transport, cookie, rpc={}):
        super(ServerBase, self).__init__(transport, cookie, rpc)
        self.log = log_server
        self._init_rpc_logger(log_client)
        self._init_client()

    def _init_client(self):
        self._transport.write({
            'type': TYPE_INIT_CONNECTION,
            'cookie': self._cookie,
            'api_version': self._api_version})

    def handle_init_reply(self, reply):
        if self._initialized:
            raise ServerException('multiple initialization replies')
        msg_type = reply['type']
        ok = reply.get('ok', False)
        if not(msg_type == TYPE_INIT_REPLY and ok):
            raise ServerException('client refused connection')
        self._initialized = True


