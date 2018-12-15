#!/usr/bin/env python

from spt3g import core, dfmux
import socket, argparse, fnmatch

parser = argparse.ArgumentParser(description='Record dfmux data to a NetCDF file based on connection to a G3NetworkSender tee in a running DAQ script', prog='ledgerman')
parser.add_argument('output', metavar='output.nc', help='Path to output NetCDF file')
parser.add_argument('server', metavar='host.name', help='Server running data acquisition')
parser.add_argument('port', metavar='port', default=5353, help='Port of DAQ server')

parser.add_argument('-m', dest='hardware_map', default=None,
                    help='Path to hardware map YAML file. Uses the network version if unspecified.')

parser.add_argument('-v', dest='verbose', action='store_true',
                    help='Verbose mode (print all frames)')
parser.add_argument('-P', dest='physnames', action='store_true',
                    help='Use physical bolometer names rather than channel names')
parser.add_argument('-p', dest='pathstring', action='store', default=[], nargs='+',
                    help='Only record channels that match channel path string,'
                    ' e.g. with a local pydfmux hardware map, 005/5/2/3/*'
                    ' saves data for all bolometers on crate 005, slot 5, '
                    'mezzanine 2, module 3. When used with a remote hardware '
                    'map, applied as a glob to the stored bolometer names.')
parser.add_argument('-b', dest='state', action='store', default=[], nargs='+',
                    help='Only record bolometers that have a state that'
                    ' matches the supplied string(s), e.g. overbiased tuned')
args = parser.parse_args()


if args.hardware_map is not None:
    # Import pydfmux later since it can take a while
    import pydfmux

    core.log_notice('Initializing hardware map and boards', unit='TeeLedgerman')
    hwm = pydfmux.load_session(open(args.hardware_map, 'r'))['hardware_map']
    if hwm.query(pydfmux.IceCrate).count() > 0:
        hwm.query(pydfmux.IceCrate).resolve()

    # Retrieve pydfmux state
    # Make sure the hardware map is consistent with what's on the IceBoards
    if args.state:
        hwm.query(pydfmux.Bolometer).load_bolo_states()

core.log_notice('Beginning data acquisition', unit='TeeLedgerman')

# Set up DfMux consumer
pipe = core.G3Pipeline()
pipe.Add(core.G3Reader, filename='tcp://' + args.server + ':' + str(args.port))

if args.hardware_map is not None:
    # Insert current hardware map, replacing whatever we got from the server,
    # if user requests it
    pipe.Add(lambda fr: fr.type != core.G3FrameType.Wiring)
    pipe.Add(dfmux.PyDfMuxHardwareMapInjector, pydfmux_hwm=hwm)

    if args.physnames:
        # Need bolo properties below
        pipe.Add(dfmux.PyDfMuxBolometerPropertiesInjector, pydfmux_hwm=hwm)

if args.physnames:
    from spt3g import calibration
    # Swap the normal bolometer names for their physical_name equivalents,
    # which can make laboratory tasks simpler -- we don't care about
    # long-term data archiving here.
    class SwapPhysLogicalNames(object):
        def __init__(self):
            self.calframe = None
            self.wiringframe = None
            self.sent = False

        def __call__(self, frame):
            if self.sent:
                return
            if frame.type == core.G3FrameType.Wiring:
                self.wiringframe = frame
                return []
            elif frame.type == core.G3FrameType.Calibration:
                self.calframe = frame
                return []
            elif frame.type != core.G3FrameType.Timepoint:
                # Just pass this along...
                return

            if self.calframe is None or self.wiringframe is None:
                raise Exception('Data before wiring and cal frames!')
            w = self.wiringframe['WiringMap']
            del self.wiringframe['WiringMap']
            c = self.calframe['NominalBolometerProperties']

            new_wiring = dfmux.DfMuxWiringMap()
            for k, d in w.iteritems():
                new_wiring[c[k].physical_name] = d
            self.wiringframe['WiringMap'] = new_wiring

            self.sent = True
            return [self.wiringframe, frame]
    pipe.Add(SwapPhysLogicalNames)

if len(args.pathstring) > 0:
    def trim_wiring_map(fr):
        if fr.type != core.G3FrameType.Wiring:
            return

        fr['OriginalWiringMap'] = fr['WiringMap']
        del fr['WiringMap']
        keys_to_keep = set()
        for item in args.pathstring:
            keys_to_keep |= set(fnmatch.filter(fr['OriginalWiringMap'].keys(), item))
        new_wiring = dfmux.DfMuxWiringMap()
        for k in keys_to_keep:
            new_wiring[k] = fr['OriginalWiringMap'][k]
        fr['WiringMap'] = new_wiring
        core.log_notice('Taking data from %d/%d detectors' % (len(keys_to_keep), len(fr['OriginalWiringMap'])), unit='TeeLedgerman')
    pipe.Add(trim_wiring_map)

if args.verbose:
    pipe.Add(core.Dump)

pipe.Add(dfmux.NetCDFDump, filename=args.output)
pipe.Run()

