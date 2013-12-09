"""
Base Class for Sage Clients


"""


from sage.rpc.client.common import (
    RemoteProcedureCaller,
    RemoteProcedureException,
    TYPE_INIT_CONNECTION,
    TYPE_INIT_REPLY,
    log_client, log_server,
)




class ClientException(RemoteProcedureException):
    pass


class ClientBase(RemoteProcedureCaller):
    
    _api_version = 'abstract_base_class'

    def __init__(self, transport, cookie, rpc={}):
        super(ClientBase, self).__init__(transport, cookie, rpc)
        self.log = log_client
        self._init_rpc_logger(log_server)

    def handle_init_connection(self, request):
        self._initialized = True
        if request.get('type') != TYPE_INIT_CONNECTION:
            raise ClientException('invalid message type')
        if request.get('cookie') != self._cookie:
            raise ClientException('authentication cookie mismatch')
        if  request.get('api_version') != self._api_version:
            raise ClientException('wrong api version')
        self._transport.write({
            'type': TYPE_INIT_REPLY,
            'ok': True})
