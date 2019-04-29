#!/usr/bin/env python

# Test that various styles of return values do the right thing in C++.
# Empty frames can confuse Boost, which will trigger exceptions on the C++
# side (hence the lack of assert(), etc. here)

from spt3g import core

class Source(object):
    def __init__(self):
        self.trig = True
    def __call__(self, frame):
        if self.trig:
            self.trig = False
            f = core.G3Frame(core.G3FrameType.Timepoint)
            f['Foo'] = 5
            return f
        return []

class Returner(object):
    def __init__(self):
        return
    def __call__(self, frame):
        print(1, frame)
        return 


class Returner2(object):
    def __init__(self):
        return
    def __call__(self, frame):
        print(2, frame.__class__, frame)
        return frame

p = core.G3Pipeline()
p.Add(Source)
p.Add(Returner)
p.Add(Returner2)
p.Run()
