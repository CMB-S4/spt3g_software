#!/usr/bin/env python

from spt3g import core, dfmux
import socket, argparse

parser = argparse.ArgumentParser(description='Record dfmux data to disk')
parser.add_argument('hardware_map', metavar='/path/to/hwm.yaml', help='Path to hardware map YAML file')
parser.add_argument('output', metavar='output', help='Base for output file names')

parser.add_argument('-v', dest='verbose', action='store_true', help='Verbose mode (print all frames)')
parser.add_argument('-a', dest='align', action='store_true', help='Align sampling. This has to happen once, but will break any existing DAN loops when run.')
parser.add_argument('-s', dest='system_time', action='store_true', help='Replace board time with system time when data received. Useful if your board is in IRIG_TEST mode.')
parser.add_argument('--max_file_size', dest='max_file_size', help='Maximum file size in MB')
args = parser.parse_args()

# Import pydfmux later since it can take a while
import pydfmux

print('Initializing hardware map and boards')
hwm = pydfmux.load_session(open(args.hardware_map, 'r'))['hardware_map']
if hwm.query(pydfmux.core.dfmux.IceCrate).count() > 0:
        hwm.query(pydfmux.core.dfmux.IceCrate).resolve()

if args.align:
	print('Aligning board sampling, this will break any existing DAN loops!')
	hwm.query(pydfmux.core.dfmux.IceBoard).set_fir_stage(6)
	hwm.query(pydfmux.core.dfmux.IceBoard).align_sampling()

print('Beginning data acquisition')
# Set up DfMux consumer
pipe = core.G3Pipeline()
builder = dfmux.DfMuxBuilder([int(board.serial) for board in hwm.query(pydfmux.IceBoard)])

# Get the local IP to connect to the boards by opening a test connection.
testsock = socket.create_connection(('iceboard' + str(hwm.query(pydfmux.core.dfmux.IceBoard).first().serial) + '.local', 80))
local_ip = testsock.getsockname()[0]
testsock.close()

collector = dfmux.DfMuxCollector(local_ip, builder)
pipe.Add(builder)

# Insert current hardware map into data stream. This is critical to get the
# board ID -> IP mapping needed to do anything useful with the data
pipe.Add(dfmux.PyDfMuxHardwareMapInjector, pydfmux_hwm=hwm)

if args.system_time:
	def sub_system_time(frame):
		if frame.type != core.G3FrameType.Timepoint:
			return
		del frame['EventHeader']
		frame['EventHeader'] = core.G3Time.Now()
	pipe.Add(sub_system_time)

pipe.Add(dfmux.PeriodicHousekeepingCollector)
pipe.Add(dfmux.HousekeepingConsumer, subprocess=True)

if args.verbose:
	pipe.Add(core.Dump)

pipe.Add(core.G3Writer, filename=args.output)

collector.Start()
pipe.Run()

