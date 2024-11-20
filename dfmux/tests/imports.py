import spt3g._libdfmux
for k in dir(spt3g._libdfmux):
    print(k, getattr(spt3g._libdfmux, k))
from spt3g._libdfmux import DfMuxSample
import spt3g.dfmux
from spt3g import dfmux
from spt3g.dfmux import DfMuxSample
