#!/usr/bin/env python

from spt3g import core
import time, sys

# Check round trips
assert(str(core.G3Time('20170101_000000')) == '01-Jan-2017:00:00:00.000000000')
assert(str(core.G3Time('01-Jan-2017:00:00:00')) == '01-Jan-2017:00:00:00.000000000')
assert(str(core.G3Time('170101_000000')) == '01-Jan-2017:00:00:00.000000000')
assert(str(core.G3Time('170101 00:00:00')) == '01-Jan-2017:00:00:00.000000000')
assert(str(core.G3Time('2017-01-01T00:00:00+0000')) == '01-Jan-2017:00:00:00.000000000')

# Check time zones for ISO dates
assert(str(core.G3Time('2017-01-01T00:00:00-0300')) == '01-Jan-2017:03:00:00.000000000')
t = core.G3Time.Now()
assert(core.G3Time(t.isoformat()) == t)

# Check that internal values are right. String identities above check for
# consistency, so we only need one case here.
assert(core.G3Time('20170101_000000').time == 148322880000000000)

