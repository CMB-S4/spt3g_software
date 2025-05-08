from . import G3Logger, G3LogLevel, usefulfunc
import traceback
from functools import reduce

__all__ = [
    "log_trace",
    "log_debug",
    "log_info",
    "log_notice",
    "log_warn",
    "log_error",
    "log_fatal",
    "set_log_level",
]


def log(log_level, unit, *args):
    message = reduce(lambda a,b: '%s %s' % ( str(a), str(b)), args)
    tb = traceback.extract_stack(limit=3)[0]
    G3Logger.global_logger.log(log_level, unit, tb[0], tb[1], tb[2], message)
    if log_level == G3LogLevel.LOG_FATAL:
        raise RuntimeError(message + " (in " + tb[2] + ")")


@usefulfunc
def log_trace(*args, unit="Python"):
    """
    Log a message for the given unit at level LOG_TRACE.
    """
    return log(G3LogLevel.LOG_TRACE, unit, *args)

@usefulfunc
def log_debug(*args, unit="Python"):
    """
    Log a message for the given unit at level LOG_DEBUG.
    """
    return log(G3LogLevel.LOG_DEBUG, unit, *args)

@usefulfunc
def log_info(*args, unit="Python"):
    """
    Log a message for the given unit at level LOG_INFO.
    """
    return log(G3LogLevel.LOG_INFO, unit, *args)

@usefulfunc
def log_notice(*args, unit="Python"):
    """
    Log a message for the given unit at level LOG_NOTICE.
    """
    return log(G3LogLevel.LOG_NOTICE, unit, *args)

@usefulfunc
def log_warn(*args, unit="Python"):
    """
    Log a message for the given unit at level LOG_WARN.
    """
    return log(G3LogLevel.LOG_WARN, unit, *args)

@usefulfunc
def log_error(*args, unit="Python"):
    """
    Log a message for the given unit at level LOG_ERROR.
    """
    return log(G3LogLevel.LOG_ERROR, unit, *args)

@usefulfunc
def log_fatal(*args, unit="Python"):
    """
    Log a message for the given unit at level LOG_FATAL,
    and raise a RuntimeError.
    """
    message, tb = log(G3LogLevel.LOG_FATAL, unit, *args)

@usefulfunc
def set_log_level(level, unit=None):
    '''
    Set log level to the requested level. If unit is not None, set the
    log level for the given logging unit only.
    
    Example
    -------
    ::

        core.set_log_level(core.G3LogLevel.LOG_DEBUG, 'GCPMuxDataDecoder')
    '''

    if unit is not None:
        G3Logger.global_logger.set_level_for_unit(unit, level)
    else:
        G3Logger.global_logger.set_level(level)


import atexit
def fix_logging_crash():
    # Unload any python loggers at exit to prevent Py_DECREF() after
    # interpreter destruction
    G3Logger.global_logger = None
atexit.register(fix_logging_crash)
