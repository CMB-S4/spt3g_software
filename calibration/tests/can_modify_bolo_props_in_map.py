#!/usr/bin/env python
from spt3g import core, calibration

# Test that parts of BolometerPropertiesMap are appropriately settable

a = calibration.BolometerProperties()
a.pol_angle = 3
a.pol_angle += 3
assert(a.pol_angle == 6)
b = calibration.BolometerPropertiesMap()
b['test'] = a
assert(b['test'].pol_angle == 6)
b['test'].pol_angle += 6
print(b['test'].pol_angle)
assert(b['test'].pol_angle == 12)

