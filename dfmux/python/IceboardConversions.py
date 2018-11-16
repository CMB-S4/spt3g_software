"""
Author : Joshua Montgomery <Joshua.Montgomery@mail.mcgill.ca>
Date   : June 16, 2014

This module contains algorithms that measure and calculate quantities related
to the readout, such as SQUID Transimpedance and detector Resistance, Current
and Voltage. It also contains the transfer functions for the ADC and Synthesizer
chains, the SQUID Bias DACs, and the streamers.


for more information, consult the doc-strings
"""

# NOTE: This file copied from pydfmux. Please ensure they stay in sync!

import numpy as np
import json
import functools
import os
from scipy import interpolate

# Figure out pydfmux directory and names and fine the json file containing
# the frequency dependent transfer function
TFpath = os.path.dirname(os.path.realpath(__file__))
TFDicts = {}

def frequency_correction(frequency, gain, target='nuller', custom_TF=None):
    """Returns fractional correction to the DC transfer function for a requested
    frequency and analog gain (defaults to baseline gain if none is provided)

    When the transferfunction_utils module is imported the pickle file
    containing the splines with the frequency corrections to the synthesizer
    transfer functions are loaded. This just selects the correct spline from that
    dictionary, evaluates it at the given frequency and returns the correction
    as a float.

    Parameters
    ----------
    frequency : float
        The frequency the correction is requested for (in HZ)

    gain : integer [0, 15]
        The analog gain of the system. The bandwidth changes slightly as
        a function of gain.

    target : string ['nuller', 'carrier'], default='nuller'

    custom_TF : string, default=None
        Optional, allows using a custom frequency correction.  Looks for a file
        in this directory with the argument name

    Returns
    -------
    float
        fractional correction to the transfer function at the provided frequency"""


    try: ## Assume frequency is an iterable
        warn = max(frequency)>9e6
    except TypeError: ## Only given one number
        warn = frequency>9e6

    if custom_TF is None:
        raise ValueError('custom_TF argument is required')

    try:
        tf_model = TFDicts[custom_TF]
    except KeyError:
        TFfile = open('{0}/{1}.txt'.format(TFpath, custom_TF), 'r')
        tf_model = json.load(TFfile)
        TFfile.close()
        TFDicts[custom_TF] = tf_model

    if 'type' in tf_model and tf_model['type'] == 'spline_interpolation':
        tf_funct = lambda x: interpolate.splev(x, tf_model[target][str(gain)])
    else:
        tf_funct = functools.partial(np.polyval, tf_model['Gain_{0}'.format(int(gain))][target])
    return tf_funct(frequency)

def gsetting_to_R(G, oldmezz=False):
    """This function converts an analog gain setting into the corresponding effective
    resistance set in the switchable resistor.

    Parameters
    ----------
    G : Integer [0, 15]
        Mezzanine gain setting

    oldmezz : Boolean, default=False
        For use with the Rev2 mezzanines (Not part of any PB2 or SPT3G allocation)

    Returns
    -------
    float
        Effective resistance of the switchable resitive network for the gain
        stage.
    """

    if G not in range(16):
        raise(Exception("Mezzanine Gain settings are integers in set [0, 15]"))

    if oldmezz:
        return 100+(15-G)*15.
    else:
        Rs = [300.0000, 212.0000, 174.2857, 153.3333, 140.0000, 130.7692,
              124.0000, 118.8235, 114.7368, 111.4286, 108.6957, 106.4000,
              104.4444, 102.7586, 101.2903, 100.0000]
        return Rs[G]

def convert_TF(gain, target, unit='NORMALIZED', frequency=False, custom_TF=None):
    """Analytic Transfer Function.

     Parameters
     ----------
    gain : Integer [0, 15]
        Analog gain setting

    target : string
        Options are contained as attributes of DFMUX boards [d.TARGET.NULLER or
        d.TARGET.CARRIER], or alternatively, 'nuller', or 'carrier'

    units : string ['NORMALIZED', 'RAW'], default='NORMALIZED'
        'NORMALIZED' for programmed or 'RAW' for dan streamer

    frequency : float, default=False
        Optional, will apply a frequency dependent correction if included.

    custom_TF : string, default=None
        Optional, allows using a custom frequency correction.  Looks for a file
        in this directory with the argument name

    Note
    ----
    All amplitudes in the DfMUX system are expressed as Peak Amplitudes.
    This function converts those quantities into real units (Volts and Amps),
    which are still expressed in Peak Amplitude. To convert to RMS, divide
    the resulting values by sqrt(2).

    Returns
    -------
    float
        If target='carrier' the return is Volts Peak Amplitude across Bolometer
        If target='nuller'  the return is Amps Peak Amplitude across the SQUID
    """

    if unit.upper() not in ['RAW', 'NORMALIZED']:
        raise(NameError('Acceptable units are "RAW" or "NORMALIZED"'))
    elif unit.upper()=='RAW':
        prefactor = 1/2.**23
    else:
        prefactor = 1

    if frequency:
        freq_corr = frequency_correction(frequency, gain, target, custom_TF)
    else:
        freq_corr = 1

    if target.lower()=='carrier':
        # The 10 in front is the range of the DAC in mA. the 1e-3 at the end
        # (with the 0.03 mOhm bias resistor) converts this to Volts
        TF_corr = 10*(300*(200/gsetting_to_R(gain))*(1/180.)*(.03))*1e-3
    elif target.lower()=='nuller':
        # The 10 in front is the range of the DAC in mA. the 1e-3 at the end
        # converts this to Amps
        TF_corr = 10*(300. * (200./gsetting_to_R(gain))*(96.77/(100+96.77)) / (750.*4)) * 1e-3
    else:
        raise NameError('Target must be NULLER or CARRIER')
    return prefactor*freq_corr*TF_corr

def convert_adc_samples(unit='ADC Counts'):
    """
    Provides the conversion factor to go from raw ADC samples or Streamer counts
    to a voltage output at the SQUID in Volts.

    Parameters
    ----------
    unit : string ['ADC Counts', 'Streamer'], default='ADC Counts'
        Starting unit. 'ADC Counts' to convert from adc units such as those
        from get_fast_samples. 'Streamer' to convert from the units used in the
        demodulator streamer channels (ledgerman, parser, get_samples, etc)

    Returns
    -------
    float
        Value in volts at the SQUID Output.
    """

    if not unit.upper()=='ADC COUNTS':
        prefactor = 1/128. # This came from a demod streamer -- between the
                           # 16-vs-24 bits scaling and the presence of the
                           # mixer, this is a factor of 128 scaling
    else:
        prefactor = 1.

    gain_1st_stage = 16 # On SQCB (topology: 1+G1/G2)
    gain_2nd_stage = 4  # On SQCB (inverting topology: G1/G2)
    gain_3rd_stage = 4  # On Mezz, with 50ohm on SQCB and 50ohn on Mezz (inverting topology: G1/G2)
    gain_4th_stage = 0.5 # Passive network between last amplifier and ADC chip
    adc_bits = 16 # SQUID has P-P voltage of 2 Volts. AROCs are in Peak Amplitudes

    vsq_per_aroc = ( (1./ (gain_1st_stage*gain_2nd_stage*gain_3rd_stage*gain_4th_stage))
                    /(2**(adc_bits-1)) )

    return vsq_per_aroc*prefactor

