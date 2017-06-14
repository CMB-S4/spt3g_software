from spt3g import core
from spt3g.dfmux import DfMuxWiringMap, DfMuxChannelMapping
import struct, socket

'''
Code to generate wiring map objects (and wiring frames in general) from
control system hardware maps.
'''

@core.indexmod
class GenerateFakeHardwareMap(object):
    '''Inserts a fake hardware map into the data stream. Takes a wiring frame as an argument to the constructor, which it will inject before any other frames.'''
    def __init__(self, frame):
        self.ran = False
        self.insert = frame
    def __call__(self, frame):
        if self.ran:
            return frame
        else:
            self.ran = True
            return [self.insert, frame]

@core.indexmod
class PyDfMuxWiringMapInjector(object):
    '''
    Insert a wiring map derived from a pydfmux hardware map into the data
    stream ahead of what would otherwise be the first frame.
   
    Optionally filter for detectors described by the mask in <pathstring>
    (see pydfmux documentation for hwm.channel_maps_from_pstring()) and
    detectors in one of the states identified by the state argument.
    '''
    def __init__(self, pydfmux_hwm, pathstring=None, state=[]):
        self.ran = False
        self.hwm = pydfmux_hwm
        self.pathstring = pathstring
        self.state = state
    def __call__(self, frame):
        if self.ran:
            return frame

        from pydfmux.core import dfmux as pydfmux
        hwmf = core.G3Frame(core.G3FrameType.Wiring)
        hwm = DfMuxWiringMap()
        
        if self.pathstring:
            chan_map_query = self.hwm.channel_maps_from_pstring(self.pathstring)
        else:
            chan_map_query = self.hwm.query(pydfmux.ChannelMapping)

        if self.state:
            for bolo in self.hwm.query(pydfmux.Bolometer):
                if bolo.readout_channel:
                    bolo.state = bolo.retrieve_bolo_state().state
            self.hwm.commit()
            chan_map_query = chan_map_query.join(pydfmux.ChannelMapping, pydfmux.Bolometer).filter(pydfmux.Bolometer.state._in(self.state))

        for bolo in chan_map_query:
            mapping = DfMuxChannelMapping()
            mapping.board_ip = struct.unpack("i", socket.inet_aton(socket.gethostbyname('iceboard' + str(bolo.iceboard.serial) + '.local')))[0]
            mapping.board_serial = int(bolo.iceboard.serial)
            mapping.board_slot = bolo.iceboard.slot if bolo.iceboard.slot else -1
            mapping.crate_serial = int(bolo.iceboard.crate.serial) if bolo.iceboard.slot else -1
            mapping.module = bolo.readout_channel.module.module - 1 # pydfmux HWMs use 1-indexing of modules, while FPGA uses 0-indexing
            if bolo.readout_channel.mezzanine.mezzanine == 2:
                mapping.module += 4
            mapping.channel = bolo.readout_channel.channel - 1 # pydfmux HWMs use 1-indexing of channels, while FPGA uses 0-indexing
            hwm[str(bolo.bolometer.global_name)] = mapping
        hwmf['WiringMap'] = hwm
        hwmf['ReadoutSystem'] = 'ICE'
        self.ran = True

        return [hwmf, frame]

# Compatibility
PyDfMuxHardwareMapInjector = PyDfMuxWiringMapInjector

@core.indexmod
def PyDfMuxBolometerPropertiesInjector(frame, pydfmux_hwm=None, angle_per_mm = 4.186*core.G3Units.deg/1000):
    '''
    Insert a calibration frame following any wiring frame containing a
    BolometerPropertiesMap named "NominalBolometerProperties" that has
    the properties of each bolometer as defined by the given pydfmux 
    hardware map.
    '''

    if frame.type != core.G3FrameType.Wiring:
        return

    from spt3g import calibration
    from pydfmux.core import dfmux as pydfmux

    cal = core.G3Frame(core.G3FrameType.Calibration)
    bpm = calibration.BolometerPropertiesMap()

    for bolo in pydfmux_hwm.query(pydfmux.Bolometer):
        bp = calibration.BolometerProperties()
        if hasattr(bolo, 'physical_name'):
            bp.physical_name = str(bolo.wafer.name) + '_' + str(bolo.physical_name)
        if hasattr(bolo, 'x_mm') and bolo.x_mm is not None:
            bp.x_offset = bolo.x_mm * angle_per_mm
        if hasattr(bolo, 'y_mm') and bolo.y_mm is not None:
            bp.y_offset = bolo.y_mm * angle_per_mm
        if bolo.wafer is not None:
            bp.wafer_id = str(bolo.wafer.name)
        if hasattr(bolo, 'pixel'):
            bp.pixel_id = bolo.pixel
        if hasattr(bolo, 'observing_band') and bolo.observing_band is not None:
            bp.band = float(bolo.observing_band)*core.G3Units.GHz
        if hasattr(bolo, 'polarization_angle') and bolo.polarization_angle is not None:
            bp.pol_angle = float(bolo.polarization_angle)*core.G3Units.deg

        bpm[str(bolo.global_name)] = bp

    cal['NominalBolometerProperties'] = bpm
 
    return [frame, cal]
    
core.indexmod
class PyDfMuxWiringMapInjectorAllChannels(object):
    '''
    Insert a wiring map derived from a pydfmux hardware map, but
    ignoring the channel mappings and adding all the channels on each defined
    module from the HWM.
    '''
    def __init__(self, pydfmux_hwm, channels=64):
        '''
        Use readout modules from provided pydfmux_hwm and add channels to it
        up to the index specified in channels.
        '''
        self.ran = False
        self.hwm = pydfmux_hwm
        self.channels = channels
    def __call__(self, frame):
        if self.ran:
            return frame

        from pydfmux.core import dfmux as pydfmux
        hwmf = core.G3Frame(core.G3FrameType.Wiring)
        hwm = DfMuxWiringMap()
        for mod in self.hwm.query(pydfmux.ReadoutModule):
            for bolo in range(self.channels):
                mapping = DfMuxChannelMapping()
                mapping.board_ip = struct.unpack("i", socket.inet_aton(socket.gethostbyname('iceboard' + str(mod.iceboard.serial) + '.local')))[0]
                mapping.board_serial = int(mod.iceboard.serial)
                mapping.board_slot = mod.iceboard.slot if mod.iceboard.slot else -1
                mapping.crate_serial = int(mod.iceboard.crate.serial) if mod.iceboard.slot else -1
                mapping.module = mod.module - 1 # pydfmux HWMs use 1-indexing of modules, while FPGA uses 0-indexing
                if mod.mezzanine.mezzanine == 2:
                    mapping.module += 4
                mapping.channel = bolo
                hwm[str(mapping)] = mapping
        hwmf['WiringMap'] = hwm
        hwmf['ReadoutSystem'] = 'ICE'
        self.ran = True

        return [hwmf, frame]

# Compatibility
PyDfMuxHardwareMapInjectorAllChannels = PyDfMuxWiringMapInjectorAllChannels

@core.indexmod
class DfmlHardwareMapInjector(object):
    def __init__(self, dfml_hwm):
        self.ran = False
        self.dfml_hwm = dfml_hwm

        self.hwm = DfMuxWiringMap()
        for ch in self.dfml_hwm.channels:
            bolo_id = str(self.dfml_hwm(ch, 'channel_id'))
            if self.dfml_hwm(ch, 'bolo_id') == '':
                # Check bolo_id, which is set iff this is a channel we care about
                continue
            mapping = DfMuxChannelMapping()
            mapping.board_ip = struct.unpack("i", socket.inet_aton( self.dfml_hwm(ch, 'dfmux_ip') ))[0]
            mapping.board_serial = int(self.dfml_hwm(ch, 'dfmux_id'))
            if self.dfml_hwm(ch, 'crate_slot') is not None:
                mapping.board_slot = int(self.dfml_hwm(ch, 'crate_slot'))
            mapping.crate_serial = -1
            mapping.module = self.dfml_hwm(ch, 'module') - 1
            mapping.channel = self.dfml_hwm(ch, 'chan_num') - 1
            self.hwm[bolo_id] = mapping
    def __call__(self, frame):
        if self.ran:
            return frame
        hwmf = core.G3Frame(core.G3FrameType.Wiring)
        hwmf['WiringMap'] = self.hwm
        hwmf['ReadoutSystem'] = 'DfMux'
        self.ran = True
        return [hwmf, frame]

