



import logging


def make_logger(name):
    try:
        from sage.doctest import DOCTEST_MODE
    except ImportError:
        DOCTEST_MODE = False
    logger = logging.getLogger(name)
    if len(logger.handlers) == 0:
        if DOCTEST_MODE:
            formatter = logging.Formatter('[%(name)s] {%(origin)s} %(levelname)s: %(message)s')
        else:
            formatter = logging.Formatter('%(asctime)s [%(name)s] {%(origin)s} %(levelname)s: %(message)s',
                                          datefmt='%H:%M:%S')
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger

class OriginLoggerAdapter(logging.LoggerAdapter):
    
    def __init__(self, logger, origin):
        super(OriginLoggerAdapter, self).__init__(logger, extra={'origin': origin})
        self._origin = origin

    def process(self, msg, kwargs):
        extra = kwargs.get('extra', {})
        extra.setdefault('origin', self._origin)
        kwargs['extra'] = extra
        return msg, kwargs

logger = make_logger('RPC')    
logger_client = OriginLoggerAdapter(logger, 'client')
logger_server = OriginLoggerAdapter(logger, 'server')
