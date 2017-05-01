import numpy
from spt3g import core, dfmux

pipe = core.G3Pipeline()
builder = dfmux.DfMuxBuilder(1)
collector = dfmux.DfMuxCollector("192.168.1.4", builder)

pipe.Add(builder)

startpoint = {}
nframes = 0

def dump(frame):
	global nframes
	#print('Frame: %s' % str(frame['EventHeader']))
	nframes += 1
	for i in frame['DfMux']:
		ip = i.key()
		data = i.data()
		for module in data:
			if module.key() not in startpoint:
				startpoint[module.key()] = {}
			samps = numpy.abs(module.data())
			for chan in range(0,16):
				if chan not in startpoint[module.key()]:
					startpoint[module.key()][chan] = 0
				if nframes < 200:
					startpoint[module.key()][chan] = max(startpoint[module.key()][chan], max(samps[chan*2], samps[chan*2 + 1]))
				rail = startpoint[module.key()][chan]*3
				if samps[chan*2] > rail or samps[chan*2 + 1] > rail:
					print('Channel %d hosed (I %d, Q %d)' % (chan, samps[chan*2], samps[chan*2 + 1]))


pipe.Add(dump)
#pipe.Add(core.G3Writer('dfmux.g3'))

collector.Start()
pipe.Run()
