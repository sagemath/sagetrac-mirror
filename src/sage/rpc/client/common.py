"""
RPC handling that is common to both Server and Client
"""


from sage.rpc.client.transport import TransportError


def make_logger(name):
    import logging
    logger = logging.getLogger(name)
    if len(logger.handlers) == 0:
        formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s: %(message)s',
                                      datefmt='%Y-%m-%d %I:%M:%S %p')
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)
    return logger
    
log_client = make_logger('client')
log_server = make_logger('server')


TYPE_INIT_CONNECTION = 'init_connection'
TYPE_INIT_REPLY = 'init_reply'
TYPE_RPC_METHOD_CALL = 'rpc_method_call'


class RemoteProcedureException(TransportError):
    pass


class RemoteProcedureUtils(object):
    
    def __init__(self, rpc_caller):
        self._rpc_caller = rpc_caller

    def ping(self, count):
        self._rpc_caller.log.debug('ping #%s', count)
        self._rpc_caller('util.pong', count)

    def pong(self, count):
        self._rpc_caller.log.debug('pong #%s', count)
        print('pong #{0}'.format(count))
    
    def quit(self):
        self._rpc_caller.log.debug('received quit command', count)
        import sys
        sys.exit(0)


class RemoteProcedureCaller(object):

    def __init__(self, transport, cookie, rpc):
        self._rpc = dict(rpc)
        utils = RemoteProcedureUtils(self)
        self._rpc['util.ping'] = utils.ping
        self._rpc['util.pong'] = utils.pong
        self._rpc['util.quit'] = utils.quit
        self._rpc_count = 0
        self._ping_count = 0
        self._transport = transport
        self._cookie = cookie
        self._initialized = False   # RPC calls are only allowed after this is True

    def _init_rpc_logger(self, logger):
        self._rpc_logger = logger
        self._rpc['log.debug']    = logger.debug
        self._rpc['log.info']     = logger.info
        self._rpc['log.warn']     = logger.warn
        self._rpc['log.error']    = logger.error
        self._rpc['log.critical'] = logger.critical

    def close(self):
        """
        Close down
        """
        self._transport.close()
        self._rpc_callee._rpc_caller = None
        
    def call(self, cmd, *args, **kwds):
        call = {'type': TYPE_RPC_METHOD_CALL,
                'cmd': cmd,
                'args': args,
                'kwds': kwds,
                'id': self._rpc_count}
        self._transport.write(call)
        self._rpc_count += 1

    __call__ = call

    def ping(self):
        self('util.ping', self._ping_count)
        self._ping_count += 1

    def can_handle(self):
        return self._transport.can_read()

    def loop(self, timeout=None):
        """
        One iteration of the main loop

        This method exist mostly for pedagogical reasons, in real
        applications the select call should be integrated into the app
        main loop.

        OUTPUT:

        Boolean. Whether anything was read or written during the loop.
        """
        import select
        rlist, wlist, xlist = self.select_args()
        rlist, wlist, xlist = select.select(rlist, wlist, xlist, timeout)
        if rlist == [] and wlist == [] and xlist == []:
            return False
        return self.select_handle(rlist, wlist, xlist)

    def select_args(self):
        """
        Arguments for a select call
        """
        transport = self._transport
        is_written = transport.is_written(False)
        fd = transport.fileno()
        return ([fd], [] if is_written else [fd], [fd])

    def select_handle(self, rlist, wlist, xlist):
        """
        Handle the result of a select call

        OUTPUT:

        Boolean. Whether anything was read or written.
        """
        transport = self._transport
        fd = transport.fileno()
        r = w = 0
        if fd in rlist:
            r = transport.nonblocking_read()
            while transport.can_read():
                self.handle()
        if fd in wlist:
            w = transport.nonblocking_write()
        if fd in xlist:
            # can this happen?
            raise RemoteProcedureException('socket in exception state')
        # print 'select', r,w, rlist, wlist, xlist
        if r == w == 0:
            self.log.debug('Nothing read or written: %s, %s, %s', rlist, wlist, xlist)
        return not(r == w == 0)

    def handle(self):
        """
        Read and deal with request.

        This is the central entry point for handling incomming
        requests.
        """
        msg = self._transport.read()
        self.log.debug('handling message %s', msg)
        return self.handle_msg(msg)

    def handle_msg(self, msg):
        msg_type = msg['type']
        if msg_type == TYPE_RPC_METHOD_CALL:   # most common
            return self.handle_rpc(msg)
        elif msg_type == TYPE_INIT_CONNECTION:
            return self.handle_init_connection(msg)
        elif msg_type == TYPE_INIT_REPLY:
            return self.handle_init_reply(msg)
        raise RemoteProcedureException('rpc type {0} does not exist'.format(msg_type))

    def handle_init_connection(self, msg):
        raise RemoteProcedureException('initialization not defined')
        self._initialized = True # Derived clasess must set this to True during initialization

    def handle_init_reply(self, msg):
        raise RemoteProcedureException('initialization reply not defined')
        self._initialized = True # Derived clasess must set this to True during initialization

    def handle_rpc(self, msg):
        """
        Handle a message that was read in :meth:`handle`.

        You can override this method to implement custom handlers.
        """
        if not self._initialized:
            raise RemoteProcedureException('RPC call before initialization')
        cmd = msg.get('cmd', None)
        args = msg.get('args', ())
        kwds = msg.get('kwds', {})
        return self.handle_method_call(cmd ,args, kwds)
        
    def handle_method_call(self, cmd, args, kwds):
        """
        Handle a rpc method call.
        """
        try:
            method = self._rpc[cmd]
        except KeyError:
            raise RemoteProcedureException('no RPC command {0} registered'.format(cmd))
        return method(*args, **kwds)
