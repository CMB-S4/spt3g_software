#!/usr/bin/env python

from spt3g.core import G3Pipeline, G3Module
import unittest

class MySubModule(G3Module):
    def __init__(self, myarg=None):
        super(MySubModule, self).__init__()

    def Process(self, frame):
        pass

class MyModule:
    def __init__(self, myarg=None):
        pass
    def __call__(self, frame):
        pass

class TestPipelineModules(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.pipe = G3Pipeline()

    def test_add_g3module(self):
        self.pipe.Add(MySubModule)
        self.pipe.Add(MySubModule, "argval")

    def test_add_generic_module(self):
        self.pipe.Add(MyModule)
        self.pipe.Add(MyModule, "argval")

if __name__ == '__main__':
    unittest.main()
