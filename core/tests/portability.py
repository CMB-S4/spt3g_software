#!/usr/bin/env python

from spt3g import core
import os, sys

# Test that we can read files written on a variety of platforms. Pass a path
# to generate test data for whatever this platform is.

testpath = os.path.join(os.environ['SPT3G_SOFTWARE_PATH'], 'core/tests/portability')

# Test data. Exercise some complicated things (STL bits) that we don't
# necessarily have control over, mapping a few primitive types.
f = core.G3Frame()
f['Five'] = 5
v = core.G3VectorDouble([2.6, 7.2])
f['Vec'] = v
v = core.G3VectorInt([17, 42, 87])
f['VecInt'] = v
m = core.G3MapDouble()
m['Six'] = 6
m['GoingOnSixteen'] = 15.9
f['Map'] = m

if len(sys.argv) > 1:
	core.G3Writer(sys.argv[1])(f)
	sys.exit(0)

# For now, we test files from big-endian (PPC64) and little-endian (amd64)
# 64-bit systems. Should include some 32-bit ones.

for test in ['test-be.g3', 'test-le.g3']:
	print(test)
	testdata = core.G3Reader(os.path.join(testpath, test))(None)[0]

	assert(testdata['Five'] == f['Five'])
	assert(len(testdata['Vec']) == len(f['Vec']))
	for i in range(len(testdata['Vec'])):
		assert(testdata['Vec'][i] == f['Vec'][i])
	assert(len(testdata['VecInt']) == len(f['VecInt']))
	for i in range(len(testdata['VecInt'])):
		assert(testdata['VecInt'][i] == f['VecInt'][i])
	assert(len(testdata['Map']) == len(f['Map']))
	assert(testdata['Map'].keys() == f['Map'].keys())
	for i in testdata['Map'].keys():
		assert(testdata['Map'][i] == f['Map'][i])
	

