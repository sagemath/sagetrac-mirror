"""
Base Class for Sage Clients


"""

import logging

from sage.rpc.core.common import (
    RemoteProcedureCaller,
    RemoteProcedureException,
    TYPE_INIT_CONNECTION,
    TYPE_INIT_REPLY,
)


def make_logger(name):
    try:
        from sage.doctest import DOCTEST_MODE
    except ImportError:
        DOCTEST_MODE = False
    logger = logging.getLogger(name)
    if len(logger.handlers) == 0:
        if DOCTEST_MODE:
            formatter = logging.Formatter('[%(name)s] %(levelname)s: %(message)s')
        else:
            formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s: %(message)s',
                                          datefmt='%H:%M:%S')
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger

class ClientLoggerAdapter(logging.LoggerAdapter):
    
    def process(self, msg, kwargs):
        return '{%s} %s' % (self.extra['origin'], msg), kwargs

logger = make_logger('RPC')    
logger_client = ClientLoggerAdapter(logger, {'origin': 'client'})
logger_server = ClientLoggerAdapter(logger, {'origin': 'server'})


class ClientException(RemoteProcedureException):
    pass


class ClientBase(RemoteProcedureCaller):
    
    _api_version = 'abstract_base_class'

    def __init__(self, transport, cookie, rpc={}, logger=None, remote_logger=None):
        self.log = logger_client if logger is None else logger
        self.remote_log = logger_server if remote_logger is None else remote_logger
        super(ClientBase, self).__init__(transport, cookie, rpc)
        self._init_rpc_logger(logger_server)

    def _init_rpc_logger(self, logger):
        self._rpc['log.debug']     = logger.debug
        self._rpc['log.info']      = logger.info
        self._rpc['log.warn']      = logger.warning
        self._rpc['log.warning']   = logger.warning
        self._rpc['log.error']     = logger.error
        self._rpc['log.exception'] = logger.exception
        self._rpc['log.critical']  = logger.critical
        self._rpc['log.fatal']     = logger.critical

    def handle_init_connection(self, request):
        self.log.debug('initializing connection %s', request)
        self._initialized = True
        if request.get('type') != TYPE_INIT_CONNECTION:
            raise ClientException('invalid message type')
        if request.get('cookie') != self._cookie:
            raise ClientException('authentication cookie mismatch')
        if  request.get('api_version') != self._api_version:
            raise ClientException('wrong api version')
        try:
            rpc_table = request['rpc_table']
        except KeyError:
            raise ClientException('no table of rpc calls')
        self._make_rpc_proxy(rpc_table)
        self._transport.write({
            'type': TYPE_INIT_REPLY,
            'rpc_table': self._rpc.keys(),
            'ok': True})

