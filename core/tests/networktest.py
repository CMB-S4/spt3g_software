#!/usr/bin/env python
from __future__ import print_function
import random, os, sys, numpy, time
from spt3g import core

port = random.randint(10000, 60000)

frames = []
for i in range(0, 20):
	f = core.G3Frame()
	f['Sequence'] = i
	f['Data'] = core.G3Timestream(numpy.zeros(100000))
	frames.append(f)

print('Port: ', port)
child = os.fork()

if child != 0:
	# Parent
	print('Parent')	

	send = core.G3NetworkSender(hostname='*', port=port)
	time.sleep(1) # XXX: how to signal that the remote end is ready?
	print('Sending')

	for f in frames:
		send(f)
	send(core.G3Frame(core.G3FrameType.EndProcessing))

	pid, status = os.wait()
	print('Child Status: ', status)
	if status == 0:
		print('OK')
	sys.exit(status)
else:
	# Child
	print('Child')

	recv = core.G3Reader(filename='tcp://localhost:%d' % port)
	rframes = []
	for k in range(len(frames)):
		chunk = recv(None)
		print(chunk[0])
		rframes += chunk
	assert(len(rframes)) == len(frames)
	for i in range(len(rframes)):
		assert(rframes[i]['Sequence'] == i)

	sys.exit(0)
