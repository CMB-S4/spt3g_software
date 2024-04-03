import numpy as np
import copy
from spt3g import core
from spt3g.gcp import ACUStatus, ACUState, TrackerStatus, TrackerState, TrackerPointing, CalFile

@core.usefulfunc
def UnitValue(caldict_entry):
    '''Turn unit name into floating point unit value'''

    if 'UnitValue' in caldict_entry:
        return caldict_entry['UnitValue']

    try: 
        uname = caldict_entry['UnitName']
        if uname and uname != 'None':
            try:
                if '/' in uname:
                    unames = list(filter(None,uname.split('/')))
                    uvalue1 = getattr(core.G3Units, 
                                      list(filter(None,unames[0].split(' ')))[0])
                    uvalue2 = getattr(core.G3Units, 
                                      list(filter(None,unames[1].split(' ')))[0])
                    uvalue = uvalue1 / uvalue2
                else:
                    uvalue = getattr(core.G3Units, uname)
            except AttributeError:
                uvalue = 1.
                core.log_warn('No entry in G3Units for ' + uname + '. Setting UnitValue to 1.0\n')
        else:
            uvalue = 1.
    except KeyError:
        uvalue = 1.

    caldict_entry['UnitValue'] = uvalue

    return uvalue


@core.usefulfunc
def CalibrateValue(data, caldict_entry):
    '''Apply gain / offset units from G3 cal file to register'''

    uvalue = UnitValue(caldict_entry)

    # if a register has units, it can't be an int anymore.  well, actually,
    # it can't be an int if we're adding floats to it or multiplying it by
    # floats either, so convert everything that has an entry in the cal file
    # to float/double.
    if 'OutputType' not in caldict_entry:
        g3type = type(data)
        if g3type == core.G3VectorInt:
            g3type = core.G3VectorDouble
        elif g3type == core.G3MapInt:
            g3type = core.G3MapDouble
        elif g3type == core.G3Int:
            g3type = core.G3Double
        caldict_entry['OutputType'] = g3type
    g3type = caldict_entry['OutputType']

    # make a copy
    if np.size(data) == 1:
        data = data.value
    data2 = np.asarray(data, dtype='float64')

    # calibrate units
    offset, recip = caldict_entry['Offset'], caldict_entry['ReciprocalFactor']
    if offset != 0:
        data2 += float(offset)
    data2 *= float(uvalue / recip)
    if not data2.shape:
        data2 = data2.tolist()

    return g3type(data2)


@core.indexmod
class CalibrateFrame:
    '''Apply gain / offset / units from G3 cal file'''

    def __init__(self, calibration_file=None):
        cf = CalFile.CalFileReader()
        self.cal = cf.readCalFile(calibration_file)

    def __call__(self, frame):

        if frame.type != core.G3FrameType.GcpSlow or 'Calibrated' in frame:
            return

        for board in frame.keys():
            cboard = frame.pop(board)
            if board not in self.cal:
                continue
            bcal = self.cal[board]
            for rmap, crmap in cboard.items():
                if rmap not in bcal:
                    continue
                rcal = bcal[rmap]
                for reg, creg in crmap.items():
                    if reg not in rcal:
                        continue
                    rcd = rcal[reg]
                    rsize = np.size(creg)
                    rshape = np.shape(creg)
                    if rsize > 1 and len(rshape) > 1:
                        for i in range(rshape[0]):
                            try:
                                rcdi = rcd[i]
                            except KeyError:
                                rcdi = rcd
                            creg[i] = CalibrateValue(creg[i], rcdi)
                    else:
                        try:
                            rcdi = rcd[0]
                        except KeyError:
                            rcdi = rcd
                        crmap[reg] = CalibrateValue(creg, rcdi)
            frame[board] = cboard

        frame['Calibrated'] = True


@core.indexmod
def UnpackACUData(f):
    '''Extracts ACU status information to ACUStatus key in frame'''

    if f.type != core.G3FrameType.GcpSlow:
        return

    board = f['antenna0']['acu']

    a = ACUStatus()
    a.time = f['antenna0']['frame']['utc']
    a.az_pos = board['az_pos'].value
    a.el_pos = board['el_pos'].value
    a.az_rate = board['az_rate'].value
    a.el_rate = board['el_rate'].value

    # 'new_*' registers not actually filled by GCP; ignore them

    a.px_checksum_error_count = board['px_checksum_error_count'].value
    a.px_resync_count = board['px_resync_count'].value
    a.px_resync_timeout_count = board['px_resync_timeout_count'].value
    a.px_resyncing = board['px_resyncing'].value
    a.px_timeout_count = board['px_timeout_count'].value
    a.restart_count = board['restart_count'].value

    a.state = ACUState(board['state'].value)
    a.status = board['acu_status'].value
    try:
        a.error = board['acu_error'].value
    except KeyError:
        # This register was some time in early 2018.  In order to read
        # older data, just set the error code to 0.
        a.error = 0

    f['ACUStatus'] = a


@core.indexmod
def UnpackTrackerMinimal(f, rewrite_source_from_feature_bits=True):
    '''
    Construct SourceName and ObservationId keys from frame.

    If rewrite_source_from_feature_bits is True (the default), will try to
    rewrite source names if DecryptFeatureBit() has been run and either
    "elnod", "calibrator", or "noise" is present in the feature bit list
    to that value.
    '''

    if f.type != core.G3FrameType.GcpSlow:
        return

    # Grab the GCP source name. If it is "current", fill in something more
    # helpful from the feature bits if possible.
    source = f['antenna0']['tracker']['source'].value
    if rewrite_source_from_feature_bits and 'GCPFeatureBits' in f:
        if 'elnod' in f['GCPFeatureBits']:
            source = 'elnod'
        if 'calibrator' in f['GCPFeatureBits'] and 'source_scan' not in f['GCPFeatureBits']:
            source = 'calibrator'
        if 'noise' in f['GCPFeatureBits']:
            source = 'noise'
        if 'debug' in f['GCPFeatureBits']:
            source = 'debug-forced-scanify'
        if 'every_pixel_on_src' in f['GCPFeatureBits']:
            source = source + '-pixelraster' # NB: Do NOT use in-place +=
    f['SourceName'] = source

    # And observation ID, if present
    if 'obs_id' in f['antenna0']['tracker']:
        f['ObservationID'] = f['antenna0']['tracker']['obs_id']


@core.indexmod
def UnpackTrackerData(f, rewrite_source_from_feature_bits=True):
    '''
    Extracts tracker status information to frame into the TrackerStatus key,
    along with the observation processing handled by UnpackTrackerMinimal.

    If rewrite_source_from_feature_bits is True (the default), will try to
    rewrite source names if DecryptFeatureBit() has been run and either
    "elnod", "calibrator", or "noise" is present in the feature bit list
    to that value.
    '''

    if f.type != core.G3FrameType.GcpSlow:
        return

    UnpackTrackerMinimal(f, rewrite_source_from_feature_bits)

    board = f['antenna0']['tracker']

    t = TrackerStatus()
    # List comprehensions are due to funny business with G3VectorFrameObject
    t.time = board['utc'][0]

    # Measured values
    t.az_pos = np.asarray(board['actual'][0])
    t.el_pos = np.asarray(board['actual'][1])
    # XXX units for rates seem to be wrong. I think this is in encoder counts
    t.az_rate = np.asarray(board['actual_rates'][0], dtype = float)
    t.el_rate = np.asarray(board['actual_rates'][1], dtype = float)
    
    # Expected values
    t.az_command = np.asarray(board['expected'][0])
    t.el_command = np.asarray(board['expected'][1])
    t.az_rate_command = np.asarray(board['expected_rates'][0], dtype = float)
    t.el_rate_command = np.asarray(board['expected_rates'][1], dtype = float)

    # Status params
    if isinstance(board['state'][0], core.G3String):
        # If state is all zero (LACKING), for example due to an ACU glitch,
        # the ARC reader may decide that the 8-bit array field is a string.
        # Treat it as one.
        t.state = [TrackerState(0) for s in board['inControl'][0]]
    else:
        t.state = [TrackerState(s) for s in board['state'][0]]
    t.acu_seq = board['acu_seq'][0]
    t.in_control = core.BoolVector(board['inControl'][0])
    t.in_control_int = core.IntVector(board['inControl'][0])
    t.scan_flag = core.BoolVector(board['scan_flag'][0])
    
    t.lst = np.asarray(board['lst'][0], dtype=float)

    t.source_acquired = np.asarray(board['off_source'][0])
    t.source_acquired_threshold = np.asarray(board['source_acquired_threshold'])
    t.tracker_mode = np.asarray(board['mode'][0])
    t.tracker_lacking = np.asarray(board['lacking'][0])
    t.time_status = np.asarray(board['time_status'][0])
    try:
        t.schedule_name = np.asarray(board['schedule_name'].value)
    except AttributeError:
        t.schedule_name = np.asarray(''.join([chr(x) for x in board['schedule_name']]))

    f['TrackerStatus'] = t


@core.indexmod
def UnpackTrackerPointingData(f):
    '''
    Extracts tracker registers relevant to online and offline pointing.
    Calibration values (offsets and multiplicative constants) are from
    gcp/control/conf/spt/cal.
    '''

    if f.type != core.G3FrameType.GcpSlow:
        return

    board = f['antenna0']['tracker']

    t = TrackerPointing()
    t.time = board['utc'][0]
    t.scu_temp = np.asarray(f['antenna0']['scu']['temp'])
    t.features = core.IntVector([f['array']['frame']['features'].value])

    t.encoder_off_x = np.asarray([board['encoder_off'][0]], dtype=np.double)
    t.encoder_off_y = np.asarray([board['encoder_off'][1]], dtype=np.double)
    
    t.low_limit_az = np.asarray([board['az_limits'][0]], dtype=np.double)
    t.high_limit_az = np.asarray([board['az_limits'][1]], dtype=np.double)
    t.low_limit_el = np.asarray([board['el_limits'][0]], dtype=np.double)
    t.high_limit_el = np.asarray([board['el_limits'][1]], dtype=np.double)

    t.tilts_x = np.asarray(board['tilt_xy_avg'][0], dtype=np.double)
    t.tilts_y = np.asarray(board['tilt_xy_avg'][1], dtype=np.double)
    t.refraction = np.asarray(board['refraction'][2], dtype=np.double)

    t.horiz_mount_x = np.asarray(board['horiz_mount'][0])
    t.horiz_mount_y = np.asarray(board['horiz_mount'][1])
    t.horiz_off_x = np.asarray(board['horiz_off'][0])
    t.horiz_off_y = np.asarray(board['horiz_off'][1])

    t.scan_off_x = np.asarray(board['scan_off'][0])
    t.scan_off_y = np.asarray(board['scan_off'][1])
    t.sky_off_x = np.asarray(board['sky_xy_off'][0])
    t.sky_off_y = np.asarray(board['sky_xy_off'][1])
    t.equat_off_x = np.asarray(board['equat_off'][0])
    t.equat_off_y = np.asarray(board['equat_off'][1])

    t.equat_geoc_ra = np.asarray(board['equat_geoc'][0])
    t.equat_geoc_dec = np.asarray(board['equat_geoc'][1])
    t.horiz_topo_az = np.asarray(board['horiz_topo'][0])
    t.horiz_topo_el = np.asarray(board['horiz_topo'][1])

    t.error_az = np.asarray(board['errors'][0])
    t.error_el = np.asarray(board['errors'][1])

    t.linsens_avg_l1 = np.asarray(board['linear_sensor_avg'][0])
    t.linsens_avg_l2 = np.asarray(board['linear_sensor_avg'][1])
    t.linsens_avg_r1 = np.asarray(board['linear_sensor_avg'][2])
    t.linsens_avg_r2 = np.asarray(board['linear_sensor_avg'][3])
    
    t.telescope_temp = np.asarray([f['array']['weather']['airTemperature'].value])
    t.telescope_pressure = np.asarray([f['array']['weather']['pressure'].value])

    f['TrackerPointing'] = t

    p = core.G3MapVectorDouble()
    pc = core.G3MapVectorDouble()

    # core model parameters
    for k in ["tilts", "flexure", "fixedCollimation"]:
        p[k] = np.asarray(board[k], dtype=np.double)
        kadj = k + "Adjust"
        if kadj in board.keys():
            v = np.asarray(board[kadj], dtype=np.double)
            if np.any(v):
                pc[k] = v

    p['time'] = np.asarray(t.time, dtype=np.double)
    if 'linear_sensor_enabled' in board.keys():
        p['linsensEnabled'] = np.asarray(board['linear_sensor_enabled'][0], dtype=np.double)
        if 'linear_sensor_coeff' in board.keys():
            p['linsensCoeffs'] = np.asarray(board['linear_sensor_coeff'], dtype=np.double)
    else:
        p['linsensEnabled'] = np.zeros_like(p['time'], dtype=np.double)

    f['OnlinePointingModel'] = p
    if len(pc.keys()):
        pc['time'] = p['time']
        f['OnlinePointingModelCorrection'] = pc


@core.indexmod
def UpdateLinearSensorDeltas(f):

    if 'TrackerPointing' not in f:
        return

    if 'LinearSensorDeltas' in f:
        return

    # Yoke dimensions in mm.
    Rs = 1652.0 * core.G3Units.mm
    Rh = 3556.0 * core.G3Units.mm
    Ry = 6782.0 * core.G3Units.mm

    l1 = f["TrackerPointing"].linsens_avg_l1
    l2 = f["TrackerPointing"].linsens_avg_l2
    r1 = f["TrackerPointing"].linsens_avg_r1
    r2 = f["TrackerPointing"].linsens_avg_r2

    # Calculate corrections in radians.
    p = core.G3MapVectorDouble()
    p['time'] = np.asarray(f['TrackerPointing'].time, dtype=np.double)
    p["delta_az"] = (Rh / (Ry * Rs)) * (l1 - l2 - r1 + r2) * core.G3Units.rad
    p["delta_el"] = (1. / (2. * Rs)) * (l2 - l1 + r2 - r1) * core.G3Units.rad
    p["delta_et"] = (1. / (2. * Ry)) * (l1 + l2 - r1 - r2) * core.G3Units.rad

    f['LinearSensorDeltas'] = p


@core.indexmod
def DecryptFeatureBit(f):
    '''
    Unpacks the GCP feature flags
    '''

    if f.type != core.G3FrameType.GcpSlow:
        return

    flag_array = core.G3VectorString()
    feature_bit = f['array']['frame']['features'].value

    flags = ['analyze', 'source_scan', 'cabin_shutter', 'elnod', 'pol_cal',
             'calibrator', 'every_pixel_on_src', 'skydip', 'optical', 'noise',
             'trail', 'el_scan', None, None, None, None, None, None, None,
             'debug']
             # Sorry... NDH

    for i in enumerate(flags):
        if feature_bit & (1 << i[0]):
            if i[1] is None:
                core.log_error('Got an unused feature bit: {:d}'.format(i[0]))
            flag_array.append(i[1])

    f['GCPFeatureBits'] = flag_array


@core.indexmod
def AddBenchData(f):
    '''
    Add the optical bench positions to the frame.
    '''
    if f.type != core.G3FrameType.GcpSlow:
        return
    bench_axes = ['y1', 'y2', 'y3', 'x4', 'x5', 'z6']

    board = f['antenna0']['scu']

    # As of 2017-08-03, SCU time is not trustworthy
    # time = board['benchSampleTime'][0]
    # For now, do this bit of evil
    time = f['antenna0']['tracker']['utc'][0]
    start = time[0]
    stop = time[-1]

    benchcom = core.G3TimestreamMap()
    benchpos = core.G3TimestreamMap()
    benchzero = core.G3TimestreamMap()
    benchoff = core.G3TimestreamMap()
    bencherr = core.G3TimestreamMap()
    bench_info = core.G3TimestreamMap()
    for i, key in enumerate(bench_axes):
        benchcom[key] = core.G3Timestream(board['benchExpected'][i])
        benchcom[key].start = start
        benchcom[key].stop = stop

        benchpos[key] = core.G3Timestream(board['benchActual'][i])
        benchpos[key].start = start
        benchpos[key].stop = stop

        benchzero[key] = core.G3Timestream(board['benchZeros'][i])
        benchzero[key].start = start
        benchzero[key].stop = stop

        benchoff[key] = core.G3Timestream(board['benchOffsets'][i])
        benchoff[key].start = start
        benchoff[key].stop = stop

        bencherr[key] = core.G3Timestream(board['benchErrors'][i])
        bencherr[key].start = start
        bencherr[key].stop = stop

    info_items = ['benchFocus', 'benchDeadBand', 'benchAcquiredThreshold',
                  'benchPrimaryState', 'benchSecondaryState', 
                  'benchFault', 'timeLocked']
    bench_info = core.G3TimestreamMap()
    for i, key in enumerate(info_items):
        bench_info[key] = core.G3Timestream(board[key][0])
        bench_info[key].start = start
        bench_info[key].stop = stop

    f['BenchPosition'] = benchpos
    f['BenchCommandedPosition'] = benchcom
    f['BenchZeros'] = benchzero
    f['BenchOffsets'] = benchoff
    f['BenchErrors'] = bencherr
    f['BenchInfo'] = bench_info
    f['BenchSampleTime'] = board['benchSampleTime'][0]
    
@core.indexmod
def UnpackCryoData(f):
    '''
    Extracts cryo information into CryoStatus key
    '''

    if f.type != core.G3FrameType.GcpSlow:
        return

    if 'cryo' not in f['array']:
        return

    board = f['array']['cryo']
    temps = board['temperature']
    heaters = board['heater_dac']

    t = core.G3MapDouble()
    t['cryo_is_valid'] = board['cryoIsValid'][0]

    # Measured values
    # He10
    t['uc_head'] = temps[0][0]
    t['ic_head'] = temps[0][1]
    t['he4_head'] = temps[0][2]
    t['he4_fb'] = temps[0][3]
    t['he4_pump'] = temps[0][4]
    t['ic_pump'] = temps[0][5]
    t['uc_pump'] = temps[0][6]
    t['he4_sw'] = temps[0][7]
    t['ic_sw'] = temps[0][8]
    t['uc_sw'] = temps[0][9]
    t['uc_stage'] = temps[0][10]
    t['lc_tower'] = temps[0][11]
    t['ic_stage'] = temps[0][12]
    t['t4k_head'] = temps[0][13]
    t['t4k_squid_strap'] = temps[0][14]
    t['t50k_head'] = temps[0][15]

    # Optics
    t['b1_50k_wbp_near'] = temps[1][0]
    t['b2_50k_wbp_far'] = temps[1][1]
    t['b3_50k_diving_board'] = temps[1][2]
    t['b4_50k_top_bot_ptc'] = temps[1][3]
    t['y1_50k_head'] = temps[1][4]
    t['y2_50k_window_strap_near'] = temps[1][5]
    t['y3_50k_tube_strap_near'] = temps[1][6]
    t['y4_50k_tube'] = temps[1][7]
    t['g1_4k_head'] = temps[1][8]
    t['g2_4k_strap'] = temps[1][9]
    t['g3_4k_lens_tab'] = temps[1][10]
    t['g4_4k_lens_tab_far'] = temps[1][11]
    t['r1_4k_top_top_ptc'] = temps[1][12]
    t['r2_50k_midop_bot_ptc'] = temps[1][13]
    t['r3_4k_lyot_flange'] = temps[1][14]
    t['r4_4k_lyot'] = temps[1][15]

    # Receiver
    t['t4k_plate_far'] = temps[2][0]
    t['t4k_strap_optics'] = temps[2][1]
    t['t4k_plate_mid'] = temps[2][2]
    t['t4k_plate_top'] = temps[2][3]
    t['t4k_plate_ptc'] = temps[2][4]
    t['t50k_harness_middle'] = temps[2][6]
    t['t50k_strap'] = temps[2][7]
    t['squid_wh1_sl1'] = temps[2][8]
    t['squid_wh5_sl1'] = temps[2][9]
    t['squid_wh3_sl7'] = temps[2][10]
    t['cal_filament'] = temps[2][11]
    t['cal_ambient1'] = temps[2][12]
    t['cal_ambient2'] = temps[2][13]
    t['cal_ambient3'] = temps[2][14]

    # Heaters
    t['heat_he4_pump'] = heaters[0][3]
    t['heat_ic_pump'] = heaters[0][4]
    t['heat_uc_pump'] = heaters[0][5]
    t['heat_he4_sw'] = heaters[0][0]
    t['heat_ic_sw'] = heaters[0][1]
    t['heat_uc_sw'] = heaters[0][2]

    f['CryoStatus'] = t
    f['CryoStatusTime'] = board['utc']


@core.indexmod
def UnpackPTData(f):
    '''Extracts pulse tube status information to PTStatus key 
    in frame'''

    if f.type != core.G3FrameType.GcpSlow:
        return

    if 'pt415' not in f['array']:
        return

    board = f['array']['pt415']

    p = core.G3MapDouble()

    p['optics_lowp'] = board['pressure_low'][0]
    p['min_optics_lowp'] = board['min_pressure_low'][0]
    p['max_optics_lowp'] = board['max_pressure_low'][0]
    p['optics_highp'] = board['pressure_high'][0]
    p['min_optics_highp'] = board['min_pressure_high'][0]
    p['max_optics_highp'] = board['max_pressure_high'][0]
    p['optics_tempoil'] = board['temp_oil'][0]
    p['min_optics_tempoil'] = board['min_temp_oil'][0]
    p['max_optics_tempoil'] = board['max_temp_oil'][0]

    p['receiver_lowp'] = board['pressure_low'][1]
    p['min_receiver_lowp'] = board['min_pressure_low'][1]
    p['max_receiver_lowp'] = board['max_pressure_low'][1]
    p['receiver_highp'] = board['pressure_high'][1]
    p['min_receiver_highp'] = board['min_pressure_high'][1]
    p['max_receiver_highp'] = board['max_pressure_high'][1]
    p['receiver_tempoil'] = board['temp_oil'][1]
    p['min_receiver_tempoil'] = board['min_temp_oil'][1]
    p['max_receiver_tempoil'] = board['max_temp_oil'][1]

    p['optics_is_valid'] = board['deviceIsValid'][0]
    p['receiver_is_valid'] = board['deviceIsValid'][1]

    f['PTStatus'] = p
    f['PTStatusTime'] = board['utc']


@core.indexmod
def UnpackMuxData(f):
    '''
    Add the DFMux data to the frame.
    '''
    if f.type != core.G3FrameType.GcpSlow:
        return

    try:
        mux = f['array']['muxHousekeeping']
        boards = mux['boardname']
    except KeyError:
        return

    fpga_temp = core.G3MapDouble()
    board_name = core.G3MapString()

    for i, bn in enumerate(boards):
        bn = str(bn).replace('"', '') # get rid of extra quotes in board name
        if bn != "":
            board_name[str(i)] = bn
            fpga_temp[str(i)] = mux['MB_TEMPERATURE_FPGA_DIE'][i]
    f['MuxFPGATemp'] = fpga_temp
    f['MuxBoardName'] = board_name
    f['MuxTime'] = mux['utc']


@core.indexmod
def UnpackWeatherData(f):
    '''Extracts weather status information to Weather key 
    in frame'''

    if f.type != core.G3FrameType.GcpSlow:
        return

    if 'weather' not in f['array']:
        return

    board = f['array']['weather']

    t = core.G3MapDouble()
    t['telescope_temp'] = board['airTemperature'].value
    t['telescope_pressure'] = board['pressure'].value
    t['inside_dsl_temp'] = board['internalTemperature'].value
    t['wind_speed'] = board['windSpeed'].value
    t['wind_direction'] = board['windDirection'].value
    t['battery'] = board['battery'].value
    t['rel_humidity'] = board['relativeHumidity'].value
    t['power'] = board['power'].value
    t['tau'] = f['array']['tipper']['tau'].value
    t['tatm'] = f['array']['tipper']['tatm'].value

    f['Weather'] = t
    f['WeatherTime'] = board['utc']
    f['TipperTime'] = f['array']['tipper']['utc']


@core.pipesegment
def ARCExtract(pipe):
    '''Extract GCP registers into native objects'''
    pipe.Add(CalibrateFrame)
    pipe.Add(UnpackACUData)
    pipe.Add(UnpackTrackerPointingData)
    pipe.Add(UpdateLinearSensorDeltas)
    pipe.Add(DecryptFeatureBit)
    pipe.Add(UnpackTrackerData)
    pipe.Add(AddBenchData)
    pipe.Add(UnpackCryoData)
    pipe.Add(UnpackPTData)
    pipe.Add(UnpackMuxData)
    pipe.Add(UnpackWeatherData)


@core.pipesegment
def ARCExtractMinimal(pipe):
    '''
    Extract bare minimum GCP registers into native objects.

    Includes only GCPFeatureBits, SourceName and ObservationID keys.
    Use ARCExtract to calibrate and extract the complete frame.
    '''
    pipe.Add(DecryptFeatureBit)
    pipe.Add(UnpackTrackerMinimal)

# Need tool for tilt meter next
