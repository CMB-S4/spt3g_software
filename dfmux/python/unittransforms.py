import numpy
from spt3g import core
from .Housekeeping import HousekeepingForBolo

# Transfer functions for 3G and SPTpol boards
from .IceboardConversions import convert_TF, convert_adc_samples
from .wtl_ConvertUtils import convert_squid, convert_demod, convert_mb

def counts_to_rms_amps(wiringmap, hkmap, bolo, system, tf=None):
    boardhk, mezzhk, modhk, chanhk = HousekeepingForBolo(hkmap, wiringmap, bolo, True)

    if system == 'ICE':
        if chanhk.dan_streaming_enable:
            tf_I = convert_TF(modhk.nuller_gain, target='nuller', custom_TF=tf, unit='RAW', frequency=chanhk.carrier_frequency/core.G3Units.Hz)
        else:
            tf_V = convert_adc_samples('Streamer')
            tf_I = tf_V/modhk.squid_transimpedance
        if chanhk.carrier_frequency != 0:
            tf_I /= numpy.sqrt(2) # Use RMS amps
    elif system == 'DfMux':
        if chanhk.dan_streaming_enable:
            squid_converter = convert_squid(units='Normalized')
            tf_I = 1. / (2.**24) * squid_converter.theo_AnToIbolo(1., Gnuller=modhk.nuller_gain, eff=1.)
        else:
            dfmux_converter = convert_demod('gcp')
            resistor='10K' # Assume this for DAN Off
            # Handle the three exceptions for SPTpol
            if bolo in ['Sq3SBpol23Ch7', 'Sq3SBpol23Ch3', 'Sq3SBpol23Ch5']:
                resistor='5K'
            # XXX above list is wrong for 100d data
            if boardhk.timestamp.mjd < 56293:
                core.log_fatal('Set correct resistor values for 2012 SPTpol data', unit='dfmux.unittransforms')
            tf_I = dfmux_converter.theo_AdcToIbolo(1., modhk.demod_gain, resistor, f_carrier_equals_f_demod=True)
            tf_I *= 256 # 3G software shifts 8 bits right compared to SPTpol, so account for difference in meaning of "counts"
    else:
        core.log_fatal('Unknown transfer function profile %s' % system, unit='dfmux.unittransforms')

    return tf_I * core.G3Units.amp

def bolo_bias_voltage_rms(wiringmap, hkmap, bolo, system, tf=None):
    boardhk, mezzhk, modhk, chanhk = HousekeepingForBolo(hkmap, wiringmap, bolo, True)

    if system == 'ICE':
        tf_V = convert_TF(modhk.carrier_gain, target='carrier', custom_TF=tf, unit='NORMALIZED', frequency=chanhk.carrier_frequency/core.G3Units.Hz)
        volts = tf_V * chanhk.carrier_amplitude * core.G3Units.V
        if chanhk.carrier_frequency != 0:
            volts /= numpy.sqrt(2) # Use RMS amps
    elif system == 'DfMux':
        squid_converter = convert_squid(units='Normalized')
        volts = squid_converter.theo_AcToVbias(chanhk.carrier_amplitude, modhk.carrier_gain) * core.G3Units.V
    else:
        core.log_fatal('Unknown transfer function profile %s' % system, unit='dfmux.unittransforms')

    return volts

@core.usefulfunc
def get_timestream_unit_conversion(from_units, to_units, bolo, wiringmap=None, hkmap=None, system=None, tf=None):
    '''
    Return the scalar conversion factor to move timestream data from a
    given system of units (Power, Resistance, Current, Counts) to another one.
    Requires a wiring map and recent housekeeping data.
    Returned quantities are RMS for currents and time-averaged for power.

    Note that this does not handle conversions to on-sky quantities (e.g. K_cmb)
    '''

    # Rather than special-casing the full mesh of conversions,
    # first calculate a conversion to watts from whatever units
    # we have, then a conversion from watts

    if wiringmap is None:
        core.log_fatal('Try to convert data with no wiring map', unit='dfmux.unittransforms')
    if hkmap is None:
        core.log_fatal('Try to convert data with no housekeeping', unit='dfmux.unittransforms')

    if from_units == to_units:
        return 1.

    if system not in ["ICE", "DfMux"]:
        core.log_fatal(f"Cannot convert data for {system} readout system", unit="dfmux.unittransforms")

    I = counts_to_rms_amps(wiringmap, hkmap, bolo, system, tf)
    V = bolo_bias_voltage_rms(wiringmap, hkmap, bolo, system, tf)

    if from_units == core.G3TimestreamUnits.Resistance:
        if to_units == core.G3TimestreamUnits.Counts:
            return V / I
        if to_units == core.G3TimestreamUnits.Current:
            return V
        if to_units == core.G3TimestreamUnits.Power:
            return V * V
    if to_units == core.G3TimestreamUnits.Resistance:
        if from_units == core.G3TimestreamUnits.Counts:
            return V / I
        if from_units == core.G3TimestreamUnits.Current:
            return V
        if from_units == core.G3TimestreamUnits.Power:
            return V * V

    # First, convert to watts
    if from_units == core.G3TimestreamUnits.Counts:
        to_watts = I * V
    elif from_units == core.G3TimestreamUnits.Current:
        to_watts = V
    elif from_units == core.G3TimestreamUnits.Power:
        to_watts = 1.
    else:
        core.log_fatal('Unknown units scheme %s in timestream for bolo %s' % (repr(from_units), bolo), unit='dfmux.unittransforms')

    # Now the conversion from watts
    if to_units == core.G3TimestreamUnits.Counts:
        from_watts = I * V
        if from_watts != 0:
            from_watts = 1./from_watts
    elif to_units == core.G3TimestreamUnits.Current:
        from_watts = V
        if from_watts != 0:
            from_watts = 1./from_watts
    elif to_units == core.G3TimestreamUnits.Power:
        from_watts = 1.
    else:
        core.log_fatal('Trying to convert timestreams to unknown units %s' % repr(to_units), unit='dfmux.unittransforms')

    # Multiply them together
    return to_watts * from_watts

@core.indexmod
class ConvertTimestreamUnits(object):
    '''
    Changes timestream units from one set of units (e.g. ADC counts) to
    another (e.g. Power). Currents and power are time averaged quantities
    (i.e. currents give RMS values).

    Note that this does not handle conversions to on-sky quantities (e.g. K_cmb)
    '''

    def __init__(self, Input='RawTimestreams', Output='CalTimestreams', Units=core.G3TimestreamUnits.Power, SkipUncalibratable=False, KeepConversionsForObservation=True):
        '''
        Copy data in the timestream map Input to the timestream map Output,
        converting the units from whatever they were to those specified by
        Units.

        If SkipUncalibratable is true, copy timestreams for which the
        unit conversions could not be evaluated into the output timestream
        in their original units. If false, throws an exception if this occurs.

        If KeepConversionsForObservation is True (default), conversion factors
        will only be evaluated once per observation. Note that this will cause
        the wrong conversions to be applied if any parameters of the mux system
        are modified during the observation but will substantially increase
        performance.
        '''
        self.units = Units
        self.input = Input
        self.output = Output

        self.wiringmap = None
        self.system = None
        self.default_tf = None
        self.hkmap = None
        self.convfactors = {}
        self.skiperrors = SkipUncalibratable
        self.keepconversions = KeepConversionsForObservation

    def __call__(self, frame):
        if 'DfMuxTransferFunction' in frame:
            self.default_tf = frame['DfMuxTransferFunction']
        if frame.type == core.G3FrameType.Wiring:
            self.wiringmap = frame['WiringMap']
            self.system = frame['ReadoutSystem']
            self.convfactors = {}
            return
        if frame.type == core.G3FrameType.Housekeeping:
            self.hkmap = frame['DfMuxHousekeeping']
            self.convfactors = {}
            return
        if self.keepconversions and frame.type == core.G3FrameType.Observation:
            self.convfactors = {}
            self.hkmap = None
            return

        if frame.type != core.G3FrameType.Scan:
            return

        # Housekeeping data can also be in Scan frames
        if 'DfMuxHousekeeping' in frame:
            if self.hkmap is None or (frame['DfMuxHousekeeping'].values()[0].timestamp != self.hkmap.values()[0].timestamp and not self.keepconversions):
                # XXX: finer-grained check for same values?
                self.hkmap = frame['DfMuxHousekeeping']
                self.convfactors = {}

        # Now the meat
        newts = core.G3TimestreamMap()
        oldts = frame[self.input]
        for bolo,ts in oldts.items():
            if ts.units not in self.convfactors:
                self.convfactors[ts.units] = {}
            if ts.units == self.units:
                newts[bolo] = ts
                continue
            elif bolo in self.convfactors[ts.units]:
                convfactor = self.convfactors[ts.units][bolo]
            else:
                tf = self.default_tf # XXX: might get from HK data
                try:
                    convfactor = get_timestream_unit_conversion(ts.units, self.units, bolo, wiringmap=self.wiringmap, hkmap=self.hkmap, system=self.system, tf=tf)
                except KeyError:
                    if not self.skiperrors:
                        raise
                    newts[bolo] = ts
                    continue
                self.convfactors[ts.units][bolo] = convfactor # And cache it

            # Convert timestream and store results
            convts = core.G3Timestream(ts)
            if core.G3TimestreamUnits.Resistance in [self.units, ts.units] and self.units != ts.units:
                convts = 1. / convts
            convts.units = self.units
            convts *= convfactor
            if convts.units != core.G3TimestreamUnits.Counts:
                convts.SetFLACCompression(False)
            newts[bolo] = convts
        frame[self.output] = newts

