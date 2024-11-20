import spt3g._libgcp
for k in dir(spt3g._libgcp):
    print(k, getattr(spt3g._libgcp, k))
from spt3g._libgcp import ACUStatus
import spt3g.gcp
from spt3g import gcp
from spt3g.gcp import ACUStatus
