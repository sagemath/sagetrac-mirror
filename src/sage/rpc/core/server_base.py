"""
Base class for Sage Server
"""

import logging

from sage.rpc.core.logging_origin import logger_client, logger_server
from sage.rpc.core.transport import Transport
from sage.rpc.core.common import (
    RemoteProcedureCaller,
    RemoteProcedureException,
    TYPE_INIT_CONNECTION,
    TYPE_INIT_REPLY,
    TYPE_REMOTE_LOG,
)



class RemoteProcedureLogger(object):
    
    def __init__(self, rpc_caller, origin='server'):
        self._extra = {'origin': origin}
        self._rpc_caller = rpc_caller

    def log(self, level, msg, *args, **kwds):
        msg = str(msg)
        args = map(str, args)
        extra = dict(self._extra)
        extra.update(kwds.get('extra', {}))
        kwds = dict((k,str(v)) for k,v in kwds.iteritems() if k != 'extra')
        kwds['extra'] = extra
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
    
    def api_version(self):
        return 'base class v1'

    def construct_rpc_table(self):
        rpc = super(ServerBase, self).construct_rpc_table()
        return rpc

    def __init__(self, transport, cookie):
        self.log = RemoteProcedureLogger(self)
        super(ServerBase, self).__init__(transport, cookie)
        self._init_client()

    def _init_client(self):
        self._transport.write({
            'type': TYPE_INIT_CONNECTION,
            'cookie': self._cookie,
            'api_version': self.api_version(),
            'rpc_table': self._rpc.keys()})

    def handle_init_reply(self, reply):
        self.log.debug('initialization reply %s', reply)
        msg_type = reply['type']
        ok = reply.get('ok', False)
        if not(msg_type == TYPE_INIT_REPLY and ok):
            error = 'client refused connection'
            self.log.critical(error)
            raise ServerException(error)
        try:
            rpc_table = reply['rpc_table']
        except KeyError:
            error = 'no table of rpc calls'
            self.log.critical(error)
            raise ServerException(error)
        self._init_remote_rpc(rpc_table)
        self.rpc.util.remote_is_ready()   # the first RPC call



