#!/usr/bin/env python

from spt3g import core
import sys

b = core.G3MapDouble()
b['gh'] = 5

a = core.G3MapDouble(b)
print(a['gh'])
if a['gh'] != b['gh']:
	sys.exit(1)

sys.exit(0)

