# -*- coding: utf-8 -*-
'''
'''
import unittest
import math


# From PyWTL. In the unlikely event PyWTL is ever updated, please synchronize.

__all__ = ['convert_mb', 'convert_mezz', 'convert_squid', 'convert_demod']

class wtl_convert:
    '''
       Equations for basic conversions to be inheritted by other derives classes.
       All methods in this class are static
       (there are no data members of the class), meaning there is never any reason to instantiate
       an instance of this class. You can always safely use  wtl_convert.method().
    '''
    #_v1 = 3.3
    global _Vadc, _drHK
    _Vadc = 3.3
    _drHK = 2**12

    def _v3_voltage_divider(self, r1, r2, r3, v1, v2):
        v3 = ((r3*v1)/r1) - (((r1*r2+r2*r3+r1*r3)*v2)/(r1*r2))
        return -v3

    def _v1_voltage_divider(self, r1, r2, v2):
        return v2 * (1 + (r1/r2))

    def _count_to_voltage(self, adc_count):
        return adc_count * (_Vadc/_drHK)

    def _voltage_to_count(self, voltage):
        return voltage * (_drHK/_Vadc)

    def temperature(self, adc_count):
        if(adc_count >= 8192):
            adc_count = adc_count - 16384;
        return adc_count / 32.0

class convert_mezz(wtl_convert):
    '''
        Handles conversions for the mezzanine board. Inherits from wtl_convert.

    '''
    def neg_voltage(self, adc_count_pos, adc_count_neg):
        return self._v3_voltage_divider(10000, 10000, 100000,
                                        self.pos_voltage(adc_count_pos),
                                        self._count_to_voltage(adc_count_neg))

    def pos_voltage(self, adc_count):
        return self._v1_voltage_divider(20000, 10000, self._count_to_voltage(adc_count))

    def demod_voltage_se(self, adc_count_pos, adc_count_demod):
        return self._v3_voltage_divider(10000, 5000, 20050,
                                        self.pos_voltage(adc_count_pos),
                                        self._count_to_voltage(adc_count_demod))

    def demod_voltage_de(self, adc_count_pos, adc_count_demod_pos, adc_count_demod_neg):
        vp = self._v3_voltage_divider(10000, 10000, 5050,
                                      self.pos_voltage(adc_count_pos),
                                      self._count_to_voltage(adc_count_demod_pos))
        vn = self._v3_voltage_divider(10000, 10000, 5050,
                                      self.pos_voltage(adc_count_pos),
                                      self._count_to_voltage(adc_count_demod_neg))
        return vp-vn




class convert_mb(wtl_convert):
    '''
    '''
    global _v1
    _v1 = 3.3

    def __init__(self, units='Normalized_16ch'):
        if units=='Normalized':
            '''
            This converts a 24 bit number to Normalized units. Useful for e.g.
            GCP housekeeping data.
            TdH - 19Apr2012
            '''
            self.human2amp_fact = 2.**23 - 1.
        elif units=='Normalized_16ch':
            # I'm not sure what this is, but I am leaving it in place.
            self.human2amp_fact = 2.**16 - 1.
        else:
            raise ValueError('convert_mb :: Unknown argument %s' % str(units))

    def neg_voltage(self, adc_count):
        return self._v3_voltage_divider(10000, 10000, 100000,
                                        _v1,
                                        self._count_to_voltage(adc_count))

    def pos_voltage(self, adc_count):
        return self._v1_voltage_divider(30100, 10000, self._count_to_voltage(adc_count))

    def freq_to_human(self, raw_frequency):
        return (raw_frequency / (2.**32-1)) * 25.e6

    def human_to_freq(self, f,round=3):
        return self.round_and_convert_frequencyHz_to_raw(f,round)

    def phase_to_human(self, raw_phase, rads=False):
        if rads:
            return (raw_phase / (2.**32-1)) * 2 * math.pi
        return (raw_phase / (2.**32-1)) * 360.

    def phase_to_raw(self, human_phase, rads=False):
        if rads:
            return human_phase*(2.**32-1)/(2*math.pi)
        else:
            human_phase%360 #Make phase between 0 and 360.
            return (human_phase*2.**32-1)/360.

    def amp_to_human(self, raw_amplitude):
        return raw_amplitude / self.human2amp_fact

    def human_to_amp(self, human_amplitude):
        return int(human_amplitude * self.human2amp_fact)

        def amp_to_raw(self, human_amp):
                return human_amp * (2.**16-1)


    def round_and_convert_frequencyHz_to_raw( self, freqHz, round=3 ):
        '''
            helper routine for set_frequency
            Funny things happen when you ask a digital sine wave generator to operate at integer multiples of the clock.
            Our clock is 25 MHz. To ensure this never happens, we ONLY allow you to set frequencies that are
            multiples of 3 Hz.
            THIS FUNCTION IS NOT SUFFICIENT. Still have a problem if the difference from a fraction of the
            clock is, itself, a fraction of the clock as the samples repeat each nth time and form noise
            spikes. New function pending.
        '''
        intFreqHz = integer_round_off(  ( integer_round_off(freqHz/round) ) * round  )
        freqRaw = integer_round_off( intFreqHz / 25000000.0 * (2**32 - 1) )
        return freqRaw

    def freqHz_to_raw(self,freq,rounding=3):
        rounded = (int(freq)/rounding)*rounding
# Work around for python bug - casting floats of a certain magnitude to integers -> 0
# Raw frequencies starting at ~6.2MHz (raw = 2**30) set it off.
# Float to long and long to int is safe.
        raw = int(int(round((rounded/25e6)*(2**32-1))))
        return raw


class convert_squid:
    '''
        Class for conversion from settings to physical units at the SQUID.
    '''
    global _dr, _Vdac
    _dr = 2**14
    _Vdac = 3.0

    def __init__(self, units=None):
      if units==None:
        print("")
        print("convert_squid was not called with the 'units' argument.")
        print("Please specify whether you're using dfmux.NORMALIZED or")
        print("dfmux.NORMALIZED_16CH units. e.g.")
        print("  convert_squid(units=b.NORMALIZED)")
        print("For those running the old (classic / pre-DAN) firmware")
        print("I'll assume you're using NORMALIZED_16CH for now,")
        print("but please replace")
        print("  convert_squid()")
        print("by")
        print("  convert_squid('Normalized_16ch')")
        print("in your code and don't worry about this.")
        print("")
        self.units = 'Normalized_16ch'
      else:
        self.units = units
      if self.units=='Normalized':
        self.amplitude_to_normalized16ch_conversion = 16.
      elif self.units=='Normalized_16ch':
        self.amplitude_to_normalized16ch_conversion = 1.
      else:
        raise Exception("convert_squid was called with an erroneous 'units' argument. Please provide dfmux.NORMALIZED or dfmux.NORMALIZED_16CH")

    def toVoltageSquidDACOffset(self, value):
        '''
            Converts raw output of squid control board DAC Offset field to volts at the DAC. This
            field sets the zero point of the DAC.
        '''
        return 6 * (value * _Vdac) / _dr

    def toVoltageSquidChannel(self, value):
        '''
            Converts raw output of squid control board channel to volts at the DAC. This applies
            to all squid board channels. ie: squid biases, flux biases, heaters and stage1_offset.
        '''
        return ((7 * _Vdac * value) / _dr) - 8

    def toRawVoltage(self, voltage):
        '''
            Converts voltage at DAC to raw (machine number) for setting on the squid control board
        '''
        return math.floor((((voltage+8.0)/7.0) * _dr)/_Vdac)

    def toVsquid(self, value, type='default'):
        ''' converts Voffset voltage to voltage at input of the 1st stage
            amplifier.  This is also the output voltage of the Squid

            type is made for custom changes.
                default : typical DfMUX and SQUID controller board
                EBEX : Rx36 and Rx43 are changed from 5K to 2.5K in EBEX, the
                       transfer function is 20/15K instead of 20/20K
        '''
        r1 = 20.
        r2 = 20.06e3
        if type=='EBEX':
            value = value*4./3.
        return value*r1/(r1+r2)

    def toFluxQuanta(self, value):
        ''' converts DAC units to DC flux applied to the SQUID
            in units of the flux quanta
        '''
        R = 50.28e3
        M = 26e-6     # mutual inductance of the squid 26uA/phi_o
        return value/R/M     # flux to dac units

    def fluxQuantaToVFluxBias(self, value):
        '''
            converts from fraction of flux quanta to the voltage difference
            at the flux bias coil. (Inverse of toFluxQuanta)
        '''
        R = 50.28e3
        M = 26e-6
        return value*R*M

    def toFluxCurrent(self, value):
        '''
        Converts V_flux_bias at the DAC to current at the flux bias inductor in
        ** Amps **
        '''
        R = 50.28e3
        return value/R

    def toSquidCurrent(self, value):
        '''
        Converts V_SQUID_bias at the DAC to current at the SQUID in
        ** Amps **
        '''
        R = 50.28e3
        return value/R

    def nullerAmpToSquidCurrent(self, amp, gain, frequency=500000):
        '''
        Convert a nuller amplitude setting to the current at the SQUID coil. Based
        on measurements done at McGill, see the following wiki page for details:
        http://kingspeak.physics.mcgill.ca/twiki/bin/view/DigitalFMux/T_null

        Notes:
            - Based on 16 channel non-DAN/DAF firmware "Antlion"

        Input units:
        amp - Nuller fractional amplitude
        gain - nuller gain number (0,1,2,3)
        frequency - Hertz

                Output units:
                current - Arms
        '''
        amp = amp * self.amplitude_to_normalized16ch_conversion
        Tgain = [1.50, 3.39, 10.4, 29.7]
        return Tgain[gain]/(1. + (float(frequency)/2.365E6)**3)**(0.5) * amp * 10**(-6)


    def theo_AnToIbolo(self, An, Gnuller, N=16, eff=0.9): # Kam Arnolds's first cut measurement of efficiency`
        '''
        Theoretical convertion from fractionnal nuller amplitude into
            A_rms at the SQUID coil from schematics
            Francois Aubin <francois.aubin@mail.mcgill.ca> July 2009

        An is the nuller fractionnal amplitude (either a pylab array or
            a scalar)
        Gnuller is the nuller gain (0-3)
        N is the multiplexing factor (8 or 16)
        returns the A_rms through the SQUID coil
        NOTE: The input (An) is an *amplitude*, the output is an *RMS* (7Nov2011 - TdH)
        '''
        '''
        if not eff==1:
            print '\n##############################################################'
            print '# Warning : An efficiency of '+str(eff)+' is used in theo_AnToIbolo.  #'
            print '# So this function does dot reflect the theoretical transfer #'
            print '# function as it is supposed.                                #'
            print '##############################################################\n'
            '''
        An = An * self.amplitude_to_normalized16ch_conversion

        #variable resistor on demod chain (mezz rev2 and 3)
        R=[2000.0, 820.0, 200.0, 0.0]

        #transfer function
        I = An * 3.3179e-3/(100.0+R[Gnuller])

        #multiplexing factor
        if N==8:
            I*=2.0
        elif N==16:
            pass
        else:
            return -1

        return I*eff


    def theo_IboloToAn(self, I, Gnuller, N=16, eff=0.9):
        '''
        Theoretical convertion from A_rms at the SQUID coil  to
            fractionnal nuller amplitude from schematics
            Francois Aubin <francois.aubin@mail.mcgill.ca> July 2009

        I is the A_rms through the SQUID coil
        Gnuller is the nuller gain (0-3)
        N is the multiplexing factor (8 or 16)
        returns the nuller fractionnal amplitude
        '''
        #variable resistor on demod chain (mezz rev2 and 3)
        R=[2000.0, 820.0, 200.0, 0.0]

        #transfer function
        An = I * (100.0+R[Gnuller])/3.3179e-3

        #multiplexing factor
        if N==8:
            An/=2.0
        elif N==16:
            pass
        else:
            return -1

        return An / self.amplitude_to_normalized16ch_conversion/eff


    def theo_AcToVbias(self,Ac,Gcarrier,N=16, correct_for_efficiency=False):
        '''
        Theoretical conversion from fractional carrier amplitude to
            V_rms across bias resistor. Calculated using schematics
            Francois Aubin <francois.aubin@mail.mcgill.ca> July 2009

        Arguments:

          Ac - is the carrier fractional amplitude and should be in the units
               provided to convert_squid. Note this can be a single value, or
               and array.

          Gcarrier - the carrier gain setting (0-3)

          N  - This is a hack which you should ignore unless you are running eight
               channel firmware, in which case you should provide convert_squid
               with 'Normalized_16ch' and set N=8, here.

          correct_for_efficiency - Set this to True if you want some kind of
                                   correction to be applied. I have no idea where
                                   this correction came from so I wouldn't do it
                                   (TdH - 19Apr2012)

        '''
        Ac = Ac * self.amplitude_to_normalized16ch_conversion
        #variable resistor on carrier chain (mezz rev2 and 3)
        R = [2000.0, 820.0, 200.0, 0.0]
        Eff = [0.9030, 0.8876, 0.8718, 0.8439]

        #transfer functions
        V = Ac * 1.5972e-3 / (100.0+R[Gcarrier])

        #multiplexing factor
        if N==8:
            V*=2.0
        elif N==16:
            pass
        else:
            return -1

        if correct_for_efficiency:
          return V * Eff[Gcarrier]
        else:
          return V


    def theo_VbiasToAc(self,V,Gcarrier,N=16, correct_for_efficiency=False):
        '''
        Theoretical convertion from V_rms across bias resistor to
            fractionnal carrier amplitude from schematics
        Francois Aubin <francois.aubin@mail.mcgill.ca> July 2009

        V is the V_rms across the bias resistor
        Gcarrier is the carrier gain (0-3)
        N is the multiplexing factor (8 or 16)
        returns the carrier fractionnal amplitude
        '''
        #variable resistor for carrier gain
        R=[2000.0, 820.0, 200.0, 0.0]
        Eff = [0.9030, 0.8876, 0.8718, 0.8439]

        #transfer function
        Ac = V * (100.0+R[Gcarrier]) / 1.5972e-3

        #multiplexing factor
        if N==8:
            Ac/=2.0
        elif N==16:
            pass
        else:
            return -1

        Ac = Ac / self.amplitude_to_normalized16ch_conversion

        if correct_for_efficiency:
          return Ac / Eff[Gcarrier]
        else:
          return Ac


    def DDStoVbias(self, Carrier_amplitude=1, Carrier_gain=0, firmware_version='16ch'):
        '''Converts the carrier gain settings to the rms voltage bias
        across a bolometer.  This function does not take strays
        into account.

        This is a measured transfer function
        Function written by Hannes in ~2008
        '''
        Carrier_amplitude = Carrier_amplitude * self.amplitude_to_normalized16ch_conversion

        if Carrier_gain == 0:
          g = 2.48
        elif Carrier_gain == 1:
          g = 5.56
        elif Carrier_gain == 2:
          g = 16.77
        elif Carrier_gain == 3:
          g = 48.47
        else:
          raise Exception("wtl_ConvertUtils.convert_squid.DDStoVbias given invalid input. Carrier gain must be either 0,1,2,3. Got %s"%(str(Carrier_gain)))

        if firmware_version=='16ch':
          z = 0.5
        elif firmware_version=='8ch':
          z=1
        else:
          raise Exception ("wtl_ConvertUtils.convert_squid.DDStoVbias given invalid input. firmware_version must be '8ch' or '16ch'. Got %s"%(str(firmware_version)))
          return False

        g0 = 2.48
        Rb = .03
        R_current = 210 #includes dewar stray resistance

        Vbias = 9.939e-3*Rb/R_current*g/g0*Carrier_amplitude*z
        return Vbias

class convert_demod:
    def __init__(self, dfmuxApi=None):
      """
          Demod / ADC conversion utilities

          Arguments:
            dfmuxApi - This argument allows the user to choose which API is used and the bit width
                       that is used. The options are as follows:
                        - True
                            The standard "dfmux" API. (24 bit numbers)
                        - False
                            The classic API (wtl_FpgaMotherboard, wtl_DataStreamer, wtl_SquidControl)
                        - 'parser'
                            Data taken from the parser. Like "True" except the numbers are 32 bits wide.
                        - 'dfmux'
                            Equivalent to True
                        - 'classic'
                            Equivalent to False
                        - 'gcp'
                            Equivalent to 'parser'
                        - None
                            Prints a warning message, otherwise identical to 'classic' or False
                        - <integer>
                            Provide the LSB bit shift relative to a 16 bit output.
                            e.g. providing 8 is equivalent to True or 'dfmux'
                            providing 16 is equivalent to 'parser' or 'gcp'
      """
      if dfmuxApi==None:
        print("""
              convert_demod was not called with the dfmuxApi flag.
              Please specify whether you're using the old API (wtl_DataStreamer) with:
              convert_demod(dfmuxApi=False)
              or the new API (libdfmux) with:
              convert_demod(dfmuxApi=True)
              or the parsed (e.g. 32 bit GCP or parser) data with:
              convert_demod(dfmuxApi='parser')

              For backwards compatibility, I'll assume the old API (wtl_DataStreamer)
              and not correct for DAN transfer function.
              Please change your code to explicitly specify the dfmuxApi flag, though.
              """)
        self._to_classic_conversion = 1.
      elif dfmuxApi==True or dfmuxApi=='dfmux':
        self._to_classic_conversion = 1./(2.**(8)) * 0.47685/0.5
      elif dfmuxApi==False or dfmuxApi=='classic':
        self._to_classic_conversion = 1.
      elif dfmuxApi=='parser' or dfmuxApi=='gcp':
        self._to_classic_conversion = 1./(2.**(16)) * 0.47685/0.5
      else:
        self._to_classic_conversion = 1./(2.**(dfmuxApi)) * 0.47685/0.5

    def _rfb_to_float(self, R_FB):
      if R_FB == 'openloop':
        feedback_resistor_ohms = 30600. # this is not actually a feedback resistor,
                                          # but rather a hack to give the correct transfer function
      elif R_FB == '10K':
        feedback_resistor_ohms = 10000.
      elif R_FB == '5K':
        feedback_resistor_ohms = 5000.
      elif R_FB == '3.3K':
        feedback_resistor_ohms = 3333.
      elif (R_FB > 500) and (R_FB < 900000): # numerical value in ohms
        feedback_resistor_ohms = float(R_FB)
      else:
          raise Exception('convert_demod._rfb_to_float :: unable to parse R_FB, which was given as '+str(R_FB))
      return feedback_resistor_ohms        
        
    def convert_rawdump(self, raw_adc_data, Gdemod, R_FB):
        '''
        Convert a rawdump (pre-demodulation timestream from the
        25 MHz, 14 bit ADC) to physical current at the SQUID input
        coil.
        7Nov2011 - TdH
        '''
        Y = [10000.,2000.,200.,0.][Gdemod] + 10.
        current_A = raw_adc_data * (1.
                                    / 2.**14 # 14 bit ADC
                                    * 2. # volts at the ADC
                                    * (50. + 100. + 50.)/100. # voltage divider before ADC
                                    * (50. + Y + 50)/10000. # gain of second stage amplifier
                                    * (50. + 50. + 50. + 50.)/(200. + 200.) # gain of first stage amplifier
                                    * (100.+82.5+1./(1./10. + 1./121.))/(500. + 500.) * (1./self._rfb_to_float(R_FB)) # voltage to current at SQUID coil
                                    * 1.13 # Franky's efficiency factor
                                    )
        return current_A

    def convert_to_volts_dc(self, raw_data, Gdemod=0, dc_bypass_resistors=5000.):
        '''
        Take data from a board that has been modified to read at DC. This assumes you are
        looking at an I channel initialized at zero frequency.

        The defaults (Gdemod=0, dc_bypass_resistors=5000. are for the sptpol calibrator readout
        board.
        '''
        data = raw_data * self._to_classic_conversion

        Y = [10000.,2000.,200.,0.][Gdemod] + 10.
        V_input = data * (1.
                          / 2.**16 # 16 bit numbers in the classic API
                          * 2. # volts at the ADC
                          * (50. + 100. + 50.)/100. # voltage divider before ADC
                          * (50. + Y + 50)/10000. # gain of second stage amplifier
                          * (50. + 50. + 2.*dc_bypass_resistors)/(200. + 200.) # gain of first stage amplifier
                          )

        return V_input

    def theo_AdcToIbolo(self,ADC,Gdemod,R_FB,f_carrier_equals_f_demod,L=10000.):
        '''
        Theoretical convertion ADC to A_rms at bolometers from schematics
        Francois Aubin <francois.aubin@mail.mcgill.ca> July 2009

        As of August 2011, this function is know to be off.  More measurements
        need to be taken more precisely, but at the moment, a factor of 1.13
        needs to be multiplied to the output of this function to get the
        \'mesured\' transfer function (current is actually higher by 13% then
        what this function reports.

        ADC is the ADC count number (either a pylab array or a scalar
        Gdemod is the demod gain (0-3)
        R_FB is the feedback resistor (10000., 5000. or 3333.)
        f_carrier_equals_f_demod is a boolean : True if f_c==f_d
                                                False if not f_c==f_d
        L is the loopgain, default assumes high loopgain (0.01% effect)
        returns the current in the bolometer in A_rms
        '''
        #variable resistor on demod chain (mezz rev2 and 3)
        ADC *= self._to_classic_conversion
        R=[10000.0, 2000.0, 200.0, 0.0]
          
        G_DF=1.9074 #1.6678=gain of square mixer, 1.9074=gain of semi-sine mixer
        #transfer function
        I = ADC * ((100.0+R[Gdemod]) * 2.3405e-9) / (self._rfb_to_float(R_FB) * G_DF) *abs((1-L)/L)

        #RMS and peak are different only if fc==fd
        if f_carrier_equals_f_demod:
            I=I/2**0.5

        return I


    def theo_IboloToAdc(self,I,Gdemod,R_FB,f_carrier_equals_f_demod,L=10000.):
        '''
        Theoretical convertion A_rms at bolometers to ADC from schematics
        Francois Aubin <francois.aubin@mail.mcgill.ca> July 2009

        I is the bolometer current in A_rms (either a pylab array or a
            scalar)
        Gdemod is the demod gain (0-3)
        R_FB is the feedback resistor (10000., 5000. or 3333.)
        f_carrier_equals_f_demod is a boolean : True if f_c==f_d
                                                False if not f_c==f_d
        L is the loopgain, default assumes high loopgain (0.01% effect)
        f_carrier_equals_f_demod is a boolean : True if f_c==f_d
                                                    False if not f_c==f_d
        returns the ADC count number
        '''
        #variable resistor on demod chain (mezz rev2 and 3)
        R=[10000.0, 2000.0, 200.0, 0.0]

        G_DF=1.90742 #1.6678=gain of square mixer, 1.9074=gain of semi-sine mixer
        #transfer function
        ADC = I * (self._rfb_to_float(R_FB) * G_DF) / ((100.0+R[Gdemod]) * 2.3405e-9) /abs((1-L)/L)

        #RMS and peak are different only if fc==fd
        if f_carrier_equals_f_demod:
            ADC=ADC*2**0.5

        ADC /= self._to_classic_conversion

        return ADC


    def theo_AdcToVopenloop(self,ADC,fd,Gdemod,f_carrier_equals_f_demod):
        '''
        Converts ADC counts to voltage at the input of the first stage
            amplifier of the SQUID controller board from schematics
        Francois Aubin <francois.aubin@mail.mcgill.ca> July 2009

        x is the ADC count number
        fd is the demodulated frequency
        Gdemod is the demod gain (0-3)
        f_carrier_equals_f_demod is a boolean : True if f_c==f_d
                                                    False if not f_c==f_d
        returns V_rms at the input of the 1st stage amp of the SQUID ctrl
        '''
        ADC *= self._to_classic_conversion
        #variable resistor in 2nd stage ampifier of demod
        R=[10000.0, 2000.0, 200.0, 0.0]
        Y=R[Gdemod]

                #ADC counts
        V=ADC
        #Digital filter gains
        V=V/1.9074
        #ADC to volt conversion
        V=V*2./16384.
        #voltage divider at input of ADC
        V=V*2.
                #2nd stage demod amplifier gain
        V=V*((100.0+Y)/10000.0)
        #1st stage demod amplifier gain
        V=V/2.
        #2nd stage SQUID controller amplifier gain
        V=V*(100.+82.5+(1./(1./10.+1./121.)))/(500.+500.)
                #current divider
        V=V*(121.+10.)/10.
        #SQUID ctrl 1st stage amplifier gain
        V=V*20./69820.
        #1MHz rolloff of the 1st dtage amplifier of the SQUID controller board
        V=V*(1+(fd/1e6)**2)**0.5

                #if demod and carrier frequencies are the same, rms is not == p
        if f_carrier_equals_f_demod:
            V=V/2**0.5

        return V


    def theo_VopenloopToAdc(self,V,fd,Gdemod,f_carrier_equals_f_demod):
        '''
        Converts voltage at the input of the first stage amplifier of the
            SQUID controller board to ADC counts from schematics
        Francois Aubin <francois.aubin@mail.mcgill.ca> July 2009

        V is the voltage at the input of the first stage amplifier of the
            SQUID controller board
            fd is the demodulated frequency
        Gdemod is the demod gain (0-3)
        f_carrier_equals_f_demod is a boolean : True if f_c==f_d
                                                    False if not f_c==f_d
        returns ADC counts
        '''
            #variable resistor in 2nd stage ampifier of demod
        R=[10000.0, 2000.0, 200.0, 0.0]
        Y=R[Gdemod]

                #ADC counts
        ADC=V
        #1MHz rolloff of the 1st dtage amplifier of the SQUID controller board
        ADC=ADC/(1+(fd/1e6)**2)**0.5
        #SQUID ctrl 1st stage amplifier gain
        ADC=ADC*69820./20.
                #current divider
        ADC=ADC*10./(121.+10.)
        #2nd stage SQUID controller amplifier gain
        ADC=ADC*(500.+500.)/(100.+82.5+(1./(1./10.+1./121.)))
        #1st stage demod amplifier gain
        ADC=ADC*2.
        #2nd stage demod amplifier gain
        ADC=ADC*10000./(100.0+Y)
        #voltage divider at input of ADC
        ADC=ADC/2.
        #volt to ADC conversion
        ADC=ADC*16384./2.
        #Digital filter gains
        ADC=ADC*1.9074

                #if demod and carrier frequencies are the same, rms is not = p
        if f_carrier_equals_f_demod:
            ADC=ADC*2**0.5

        ADC /= self._to_classic_conversion

        return ADC

    #######################################################################################
    # Converts ADC counts_DC to current_rms at bolometers and includes transfer functions #
    #######################################################################################
    def ADCtoItxfr(self,x,Gdemod,R_FB,freq):
        x *= self._to_classic_conversion
        R=[10000.0, 2000.0, 200.0, 0.0]
        m=[-0.2696e-7, -0.2691e-7, -0.5229e-7, -1.8598e-7]
        b=[1.0054, 1.0054, 1.0105, 1.0372]

        G_DF=1.90742 #1.6678=gain of square mixer, 1.9074=gain of semi-sine mixer
        rt2 = 1.41421356
        I = x * ((100.0+R[Gdemod]) * 2.3405e-9) / (self._rfb_to_float(R_FB) * G_DF * rt2)
        if freq >= 200000:
            Itxfr = I / (m[Gdemod] * freq + b[Gdemod])
        else:
            Itxfr = I

        return Itxfr

    #########################################################################################
    # Converts current_rms at bolometers into ADC counts_DC and includes transfer functions #
    #########################################################################################
    def ItoADCtxfr(self,I,Gdemod,R_FB,freq):

        R=[10000.0, 2000.0, 200.0, 0.0]
        m=[-0.2696e-7, -0.2691e-7, -0.5229e-7, -1.8598e-7]
        b=[1.0054, 1.0054, 1.0105, 1.0372]

        G_DF=1.90742 #1.6678=gain of square mixer, 1.9074=gain of semi-sine mixer
        rt2 = 1.41421356
        x = I * (self._rfb_to_float(R_FB) * G_DF * rt2) / ((100.0+R[Gdemod]) * 2.3405e-9)
        if freq >= 200000:
            xtxfr = x * (m[Gdemod] * freq + b[Gdemod])
        else:
            xtxfr = x
        x /= self._to_classic_conversion
        return xtxfr

##
##
#######################################################
##
##  T E S T  C A S E S
##
#######################################################
##
##
class ConvertTestCase (unittest.TestCase):

    def setUp(self):
        #print "Setting up WtlTest cases"
        self.module_mb = convert_mb()
        self.module_mezz = convert_mezz()

    def tearDown(self):
        #print "Cleaning up WtlTest cases"
        self.module_mb = None
        self.module_mezz = None

    def test_01_Convert_Temperature(self):
        '''
            Test the temperature conversion. These examples are found in
            Table 5 of the ADT7301 datasheet.
        '''

        trials = {
            0x3b00 : -40.0,
            0x3c40 : -30.0,
            0x3ce0 : -25.0,
            0x3ec0 : -10.0,
            0x3fff : -0.03125,
            0x0000 : 0,
            0x0001 : 0.03125,
            0x0140 : +10,
            0x0320 : +25,
            0x0640 : +50,
            0x0960 : +75,
            0x0c80 : +100,
            0x0fa0 : +125,
            0x12c0 : +150
        }

        for value, expected in trials.items():
            calculated = self.module_mb.temperature(value)
            self.assertEqual(calculated, expected,
                'error converting 0x%08x to mb temperature %f (calculated %f)' %
                (value, expected, calculated))


    def test_02_Convert_MB_neg(self):
        ''' 2: Test the MB Negative conversion
        '''
        self.assertEqual(self.module_mb.neg_voltage(1520), -7.283203125, 'error converting mb neg voltage')
        self.assertEqual(self.module_mb.neg_voltage(2520), 9.6357421875, 'error converting mb neg voltage')

    def test_03_Convert_MB_pos(self):
        ''' 2: Test the MB Positive conversion
        '''
        self.assertEqual(self.module_mb.pos_voltage(1024), 3.2999999999999998, 'error converting mb pos voltage')
        self.assertEqual(self.module_mb.pos_voltage(4096), 13.199999999999999, 'error converting mb pos voltage')

    def test_04_Convert_Mezz_neg(self):
        ''' 2: Test the Mezz Negative conversion
        '''
        self.assertEqual(self.module_mezz.neg_voltage(1768, 2020), -8.5561523437499929, 'error converting mezz neg voltage')
        self.assertEqual(self.module_mezz.neg_voltage(1224, 1520), -3.8671874999999929, 'error converting mezz neg voltage')

    def test_05_Convert_Mezz_pos(self):
        ''' 2: Test the Mezz Positive conversion
        '''
        self.assertEqual(self.module_mezz.pos_voltage(1024), 2.4749999999999996, 'error converting mezz pos voltage')
        self.assertEqual(self.module_mezz.pos_voltage(4096), 9.8999999999999986, 'error converting mezz pos voltage')


def suite():
    suite = unittest.makeSuite(ConvertTestCase,'test')
    return suite

if __name__ == '__main__':
    unittest.main()
