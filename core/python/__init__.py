# We need to share certain global variables (cereal and pybind11
# type registries between shared objects, so make sure the RTLD
# is configured to allow that to happen. This can sometimes work
# without this on glibc due to a glibc bug, but it will fail on
# other libc implementations
import os, sys
sys.setdlopenflags(sys.getdlopenflags() | os.RTLD_GLOBAL)

from .._libcore import *

from .modconstruct import usefulfunc, pipesegment, indexmod
from .g3logging import *

from .fileio import *
try:
	from .parse_pipeline_graph import plot_frame_processing_info
except ImportError:
	pass
try:
	from .multiprocess import Subproc
except ImportError:
	pass
from .util import *

from .apidoc import *
from .timestreamextensions import *
from .quatextensions import *

from .g3decorators import cache_frame_data, scan_func_cache_data
