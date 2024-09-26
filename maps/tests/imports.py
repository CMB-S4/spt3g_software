import spt3g._libmaps
for k in dir(spt3g._libmaps):
    print(k, getattr(spt3g._libmaps, k))
from spt3g._libmaps import FlatSkyMap
import spt3g.maps
from spt3g import maps
from spt3g.maps import FlatSkyMap
