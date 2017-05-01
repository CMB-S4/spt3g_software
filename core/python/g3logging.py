from spt3g.core import G3Logger, G3LogLevel
import traceback
from functools import reduce

def _make_logger(log_level):
    def logger(*args, **kwargs):
        '''
        Behaves like print but with logging behavior as described in logging.rst

        The only kwarg that it recognizes is "unit", for setting the logging unit.
        '''
        message = reduce(lambda a,b: '%s %s' % ( str(a), str(b)), args)
        unit = kwargs.get('unit')
        unit = unit if unit != None else 'Python'
        tb = traceback.extract_stack(limit=2)[0]
        G3Logger.global_logger.log(log_level, unit, tb[0], tb[1],
                       tb[2], message)
    return logger

log_trace = _make_logger(G3LogLevel.LOG_TRACE)
log_debug = _make_logger(G3LogLevel.LOG_DEBUG)
log_info = _make_logger(G3LogLevel.LOG_INFO)
log_notice = _make_logger(G3LogLevel.LOG_NOTICE)
log_warn = _make_logger(G3LogLevel.LOG_WARN)
log_error = _make_logger(G3LogLevel.LOG_ERROR)

def log_fatal(*args, **kwargs):
    '''
    Behaves like print but with logging behavior as described in logging.rst
    
    The only kwarg that it recognizes is "unit", for setting the logging unit.
    '''
    message = reduce(lambda a,b: '%s %s' % ( str(a), str(b)), args)
    unit = kwargs.get('unit')
    unit = unit if unit != None else 'Python'    
    tb = traceback.extract_stack(limit=2)[0]
    G3Logger.global_logger.log(G3LogLevel.LOG_FATAL, unit, tb[0], tb[1],
        tb[2], message)
    raise RuntimeError(message + " (in " + tb[2] + ")")

def set_log_level(level, unit=None):
    '''
    Set log level to the requested level. If unit is not None, set the
    log level for the given logging unit only.
    
    Example: core.set_log_level(core.G3LogLevel.LOG_DEBUG, 'GCPMuxDataDecoder')
    '''

    if unit is not None:
        G3Logger.global_logger.set_level_for_unit(unit, level)
    else:
        G3Logger.global_logger.set_level(level)

