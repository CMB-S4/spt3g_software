#!/usr/bin/env python

import pickle, numpy, sys
from spt3g import core

dat = numpy.floor(numpy.random.normal(size=1200, scale=20, loc=0))
ts = core.G3Timestream(dat)
ts.units = core.G3TimestreamUnits.Counts
ts.SetFLACCompression(5)
ts[5] = numpy.nan
dat[5] = numpy.nan

pickle.dump(ts, open('tsdump.pkl', 'wb'))

ts_rehyd = pickle.load(open('tsdump.pkl', 'rb'))
print('Original')
print(ts.units)
print(ts[4], ts[5], ts[6])
print('Rehydrated')
print(ts_rehyd.units)
print(ts_rehyd[4], ts_rehyd[5], ts_rehyd[6])

if ts_rehyd.units != ts.units:
	print('Units do not match')
	sys.exit(1)

if numpy.isfinite(ts_rehyd[5]):
	print('Element 5 finite!')
	sys.exit(1)

if (numpy.asarray(ts_rehyd)[numpy.isfinite(ts_rehyd)] != dat[numpy.isfinite(dat)]).any():
	print('Elements not preserved!')
	sys.exit(1)

