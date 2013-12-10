"""
Base class for Sage Server
"""

import logging

from sage.rpc.core.transport import Transport
from sage.rpc.core.common import (
    RemoteProcedureCaller,
    RemoteProcedureException,
    TYPE_INIT_CONNECTION,
    TYPE_INIT_REPLY,
    TYPE_REMOTE_LOG,
)



class RemoteProcedureLogger(object):
    
    def __init__(self, rpc_caller):
        self._rpc_caller = rpc_caller

    def log(self, level, msg, *args, **kwds):
        msg = str(msg)
        args = map(str, args)
        kwds = dict((k,str(v)) for k,v in kwds.iteritems())
        self._rpc_caller._transport.write({
            'type': TYPE_REMOTE_LOG,
            'level': level,
            'msg': msg,
            'args': args,
            'kwds': kwds})

    def debug(self, *args, **kwds):
        self.log(logging.DEBUG, *args, **kwds)

    def info(self, *args, **kwds):
        self.log(logging.INFO, *args, **kwds)

    def warning(self, *args, **kwds):
        self.log(logging.WARNING, *args, **kwds)
    
    warn = warning

    def error(self, *args, **kwds):
        self.log(logging.ERROR, *args, **kwds)

    exception = error

    def critical(self, *args, **kwds):
        self.log(logging.CRITICAL, *args, **kwds)

    fatal = critical



class ServerException(RemoteProcedureException):
    pass



class ServerBase(RemoteProcedureCaller):
    
    _api_version = 'abstract_base_class'

    def __init__(self, transport, cookie, rpc={}):
        self.log = RemoteProcedureLogger(self)
        super(ServerBase, self).__init__(transport, cookie, rpc)
        self._init_client()

    def _init_client(self):
        self._transport.write({
            'type': TYPE_INIT_CONNECTION,
            'cookie': self._cookie,
            'api_version': self._api_version,
            'rpc_table': self._rpc.keys()})

    def handle_init_reply(self, reply):
        self.log.debug('initialization reply %s', reply)
        if self._initialized:
            raise ServerException('multiple initialization replies')
        msg_type = reply['type']
        ok = reply.get('ok', False)
        if not(msg_type == TYPE_INIT_REPLY and ok):
            raise ServerException('client refused connection')
        self._initialized = True



