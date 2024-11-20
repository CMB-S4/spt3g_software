import spt3g
import spt3g._libcore
for k in dir(spt3g._libcore):
    print(k, getattr(spt3g._libcore, k))
from spt3g._libcore import BoolVector
import spt3g.core
from spt3g import core
from spt3g.core import BoolVector
