from spt3g.core.load_pybindings import load_pybindings
load_pybindings(__name__, __path__)

from spt3g.core.g3logging import log_trace, log_debug, log_info, log_notice, log_warn, log_error, log_fatal, set_log_level

import atexit
def fix_logging_crash():
    # Unload any python loggers at exit to prevent Py_DECREF() after
    # interpreter destruction
    G3Logger.global_logger = None
atexit.register(fix_logging_crash)
del fix_logging_crash

from spt3g.core.fileio import *
from spt3g.core.modconstruct import pipesegment, indexmod, pipesegment_nodoc
from spt3g.core.funcconstruct import usefulfunc
try:
	from spt3g.core.parse_pipeline_graph import plot_frame_processing_info
except ImportError:
	pass
try:
	from spt3g.core.multiprocess import Subproc
except ImportError:
	pass
from spt3g.core.util import *

from spt3g.core.docparser import *
from spt3g.core.dataextensions import *
from spt3g.core.frameextensions import *
from spt3g.core.timestreamextensions import *

from spt3g.core.g3decorators import cache_frame_data, scan_func_cache_data
