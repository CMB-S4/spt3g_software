try:
    from .._libcore import *
except ImportError:
    from .load_pybindings import load_pybindings
    load_pybindings(__name__, __path__)

from .g3logging import log_trace, log_debug, log_info, log_notice, log_warn, log_error, log_fatal, set_log_level

import atexit
def fix_logging_crash():
    # Unload any python loggers at exit to prevent Py_DECREF() after
    # interpreter destruction
    G3Logger.global_logger = None
atexit.register(fix_logging_crash)
del fix_logging_crash

from .fileio import *
from .modconstruct import pipesegment, indexmod, pipesegment_nodoc
from .funcconstruct import usefulfunc
try:
	from .parse_pipeline_graph import plot_frame_processing_info
except ImportError:
	pass
try:
	from .multiprocess import Subproc
except ImportError:
	pass
from .util import *

from .docparser import *
from .dataextensions import *
from .frameextensions import *
from .timestreamextensions import *

from .g3decorators import cache_frame_data, scan_func_cache_data
