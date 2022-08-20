#!/usr/bin/env python

from spt3g import core
import numpy
import pickle

f64 = numpy.ones(20,dtype='float64')
f32 = numpy.ones(20,dtype='float32')
i32 = numpy.ones(20,dtype='int32')
i64 = numpy.ones(20,dtype='int64')
i16 = numpy.ones(20,dtype='int16') # This one shouldn't work directly

for d in [f64, f32, i32, i64, i16]:
	print(d.dtype)
	try:
		t = core.G3Timestream(d)
	except TypeError:
		if d.dtype == 'int16': # Meant to fail
			continue
		raise
	else:
		if d.dtype == 'int16': # Meant to fail
			raise TypeError('int16 should not work')


	assert(len(t) == 20)
	assert(t[12] == 1.0)
	t[12] = 5
	assert(t[12] == 5.0)
	assert(t[13] == 1.0)
	d2 = numpy.asarray(t)
	assert(d2.dtype == d.dtype)
	assert(not (d2 == d).all())
	
	# Test shared buffer, both ways
	t[12] = 1
	assert((d2 == d).all())
	d2[12] = 3
	assert(t[12] == 3)

	# Test serialization
	buf = pickle.dumps(t)
	t2 = pickle.loads(buf)
	assert(numpy.asarray(t2).dtype == d.dtype) # Datatype preserved
	assert((numpy.asarray(t2) == d2).all())

	# Test copy constructor
	t2 = core.G3Timestream(t)
	assert(numpy.asarray(t2).dtype == d.dtype) # Datatype preserved
	assert((numpy.asarray(t2) == d2).all())

