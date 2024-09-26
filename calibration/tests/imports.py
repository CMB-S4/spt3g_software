import spt3g._libcalibration
for k in dir(spt3g._libcalibration):
    print(k, getattr(spt3g._libcalibration, k))
from spt3g._libcalibration import BolometerProperties
import spt3g.calibration
from spt3g import calibration
from spt3g.calibration import BolometerProperties
