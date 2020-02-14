import numpy, copy
from spt3g import core
from spt3g.gcp import ACUStatus, ACUState, TrackerStatus, TrackerState, TrackerPointing, CalFile

@core.indexmod
def UnitValue(caldict_entry):
    '''Turn unit name into floating point unit value'''

    try: 
        uname = caldict_entry['UnitName']
        if uname is not 'None':
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

    return uvalue


@core.indexmod
def CalibrateFrame(f, calibration_file=None):
    '''Apply gain / offset / units from G3 cal file'''
    
    if f.type != core.G3FrameType.GcpSlow:
        return

    try:
        if f['Calibrated'] == True:
            print('Already calibrated!\n')
            return
    except KeyError:
        f['Calibrated'] = True

    cf = CalFile.CalFileReader()
    cd = cf.readCalFile(calibration_file)

    for board in f.keys():
        if board == 'Calibrated':
            continue
        cboard = copy.deepcopy(f[board])
        for rmap in cboard.keys():
            for reg in cboard[rmap].keys():
                try: 
                    rcd = cd[board][rmap][reg]
                except KeyError:
                    continue
                rsize = numpy.size(cboard[rmap][reg])
                if rsize > 1:
                    rshape = numpy.shape(cboard[rmap][reg])
                    if len(rshape) > 1:
                        for i in range(rshape[0]):
                            try:
                                rcdi = rcd[i]
                            except KeyError:
                                rcdi = rcd
                            uvalue = UnitValue(rcdi)
                            datatemp = numpy.asarray(cboard[rmap][reg][i])
                            datatemp2 = datatemp.copy()
                            # if a register has units, it can't be an
                            # int anymore.
                            # well, actually, it can't be an int if
                            # we're adding floats to it or multiplying
                            # it by floats either, so convert
                            # everything that has an entry in the cal
                            # file to float/double.
                            datatemp2 = numpy.asarray(datatemp2,dtype='float64')
                            thisdtype = datatemp2.dtype
                            datatemp2 += \
                                numpy.array(rcdi['Offset'],dtype=thisdtype)
                            datatemp2 *= numpy.array(uvalue / 
                                                     rcdi['ReciprocalFactor'],
                                                     dtype=thisdtype)
                            if type(cboard[rmap][reg][i]) \
                                    is core.G3VectorInt:
                                regitemp = core.G3VectorDouble(datatemp2)
                            elif type(cboard[rmap][reg][i]) \
                                    is core.G3MapInt:
                                regitemp = core.G3MapDouble(datatemp2)
                            elif type(cboard[rmap][reg][i]) \
                                    is core.G3Int:
                                regitemp = core.G3Double(datatemp2)
                            else:
                                regitemp = \
                                    (type(cboard[rmap][reg][i]))(datatemp2)
                            cboard[rmap][reg][i] = regitemp
                    else:
                        try:
                            rcdi = rcd[0]
                        except KeyError:
                            rcdi = rcd
                        uvalue = UnitValue(rcdi)
                        datatemp = numpy.asarray(cboard[rmap][reg])
                        datatemp2 = datatemp.copy()
                        # if a register has units, it can't be an
                        # int anymore. well, actually (see above)...
                        datatemp2 = numpy.asarray(datatemp2,dtype='float64')
                        thisdtype = datatemp2.dtype
                        datatemp2 += \
                            numpy.array(rcdi['Offset'],dtype=thisdtype)
                        datatemp2 *= numpy.array(uvalue / rcdi['ReciprocalFactor'],dtype=thisdtype)
                        if type(cboard[rmap][reg]) \
                                is core.G3VectorInt:
                            regtemp = core.G3VectorDouble(datatemp2)
                        elif type(cboard[rmap][reg]) \
                                is core.G3MapInt:
                            regtemp = core.G3MapDouble(datatemp2)
                        elif type(cboard[rmap][reg]) \
                                is core.G3Int:
                            regtemp = core.G3Double(datatemp2)
                        else:
                            regtemp = \
                                (type(cboard[rmap][reg]))(datatemp2)
                        cboard[rmap][reg] = regtemp
                else:
                    try:
                        rcdi = rcd[0]
                    except KeyError:
                        rcdi = rcd
                    uvalue = UnitValue(rcdi)
                    datatemp = cboard[rmap][reg].value
                    datatemp2 = datatemp
                    # if a register has units, it can't be an
                    # int anymore. well, actually (see above)...
                    datatemp2 = numpy.float(datatemp2)
                    datatemp2 = datatemp2 + rcdi['Offset']
                    datatemp2 *= uvalue / rcdi['ReciprocalFactor']
                    if type(cboard[rmap][reg]) \
                            is core.G3VectorInt:
                        regtemp = core.G3VectorDouble(datatemp2)
                    elif type(cboard[rmap][reg]) \
                            is core.G3MapInt:
                        regtemp = core.G3MapDouble(datatemp2)
                    elif type(cboard[rmap][reg]) \
                            is core.G3Int:
                        regtemp = core.G3Double(datatemp2)
                    else:
                        regtemp = \
                            (type(cboard[rmap][reg]))(datatemp2)
                    cboard[rmap][reg] = regtemp
        del f[board]
        f[board] = cboard

@core.indexmod
def UnpackACUData(f):
    '''Extracts ACU status information to ACUStatus key in frame'''

    if f.type != core.G3FrameType.GcpSlow:
        return

    a = ACUStatus()
    a.time = f['antenna0']['frame']['utc']
    a.az_pos = f['antenna0']['acu']['az_pos'].value
    a.el_pos = f['antenna0']['acu']['el_pos'].value
    a.az_rate = f['antenna0']['acu']['az_rate'].value
    a.el_rate = f['antenna0']['acu']['el_rate'].value

    # 'new_*' registers not actually filled by GCP; ignore them

    a.px_checksum_error_count = f['antenna0']['acu']['px_checksum_error_count'].value
    a.px_resync_count = f['antenna0']['acu']['px_resync_count'].value
    a.px_resync_timeout_count = f['antenna0']['acu']['px_resync_timeout_count'].value
    a.px_resyncing = f['antenna0']['acu']['px_resyncing'].value
    a.px_timeout_count = f['antenna0']['acu']['px_timeout_count'].value
    a.restart_count = f['antenna0']['acu']['restart_count'].value

    a.state = ACUState(f['antenna0']['acu']['state'].value)
    a.status = f['antenna0']['acu']['acu_status'].value
    a.error = f['antenna0']['acu']['acu_error'].value

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
        if 'calibrator' in f['GCPFeatureBits']:
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

    t = TrackerStatus()
    # List comprehensions are due to funny business with G3VectorFrameObject
    t.time = [tm for tm in f['antenna0']['tracker']['utc'][0]]

    # Measured values
    t.az_pos = numpy.asarray(f['antenna0']['tracker']['actual'][0])
    t.el_pos = numpy.asarray(f['antenna0']['tracker']['actual'][1])
    # XXX units for rates seem to be wrong. I think this is in encoder counts
    t.az_rate = numpy.asarray(f['antenna0']['tracker']['actual_rates'][0],
                              dtype = float)
    t.el_rate = numpy.asarray(f['antenna0']['tracker']['actual_rates'][1],
                              dtype = float)
    
    # Expected values
    t.az_command = numpy.asarray(f['antenna0']['tracker']['expected'][0])
    t.el_command = numpy.asarray(f['antenna0']['tracker']['expected'][1])
    t.az_rate_command = numpy.asarray(f['antenna0']['tracker']['expected_rates'][0], dtype = float)
    t.el_rate_command = numpy.asarray(f['antenna0']['tracker']['expected_rates'][1], dtype = float)

    # Status params
    if isinstance(f['antenna0']['tracker']['state'][0], core.G3String):
        # If state is all zero (LACKING), for example due to an ACU glitch,
        # the ARC reader may decide that the 8-bit array field is a string.
        # Treat it as one.
        t.state = [TrackerState(0) for s in f['antenna0']['tracker']['inControl'][0]]
    else:
        t.state = [TrackerState(s) for s in f['antenna0']['tracker']['state'][0]]
    t.acu_seq = f['antenna0']['tracker']['acu_seq'][0]
    t.in_control = core.BoolVector(f['antenna0']['tracker']['inControl'][0])
    t.in_control_int = core.IntVector(f['antenna0']['tracker']['inControl'][0])
    t.scan_flag = core.BoolVector(f['antenna0']['tracker']['scan_flag'][0])
    
    t.lst = numpy.asarray(f['antenna0']['tracker']['lst'][0], dtype=float)

    t.source_acquired = numpy.asarray(f['antenna0']['tracker']['off_source'][0])
    t.source_acquired_threshold = numpy.asarray(f['antenna0']['tracker']['source_acquired_threshold'])
    t.tracker_mode = numpy.asarray(f['antenna0']['tracker']['mode'][0])
    t.tracker_lacking = numpy.asarray(f['antenna0']['tracker']['lacking'][0])
    t.time_status = numpy.asarray(f['antenna0']['tracker']['time_status'][0])
    try:
        t.schedule_name = numpy.asarray(f['antenna0']['tracker']['schedule_name'].value)
    except AttributeError:
        t.schedule_name = numpy.asarray(''.join([chr(x) for x in f['antenna0']['tracker']['schedule_name']]))

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

    t = TrackerPointing()
    t.time = [tm for tm in f['antenna0']['tracker']['utc'][0]]
    t.scu_temp = numpy.asarray(f['antenna0']['scu']['temp'])
    t.features = core.IntVector([f['array']['frame']['features'].value])

    t.encoder_off_x = numpy.asarray([f['antenna0']['tracker']['encoder_off'][0]], dtype=numpy.double)
    t.encoder_off_y = numpy.asarray([f['antenna0']['tracker']['encoder_off'][1]], dtype=numpy.double)
    
    t.low_limit_az = numpy.asarray([f['antenna0']['tracker']['az_limits'][0]], dtype=numpy.double)
    t.high_limit_az = numpy.asarray([f['antenna0']['tracker']['az_limits'][1]], dtype=numpy.double)
    t.low_limit_el = numpy.asarray([f['antenna0']['tracker']['el_limits'][0]], dtype=numpy.double)
    t.high_limit_el = numpy.asarray([f['antenna0']['tracker']['el_limits'][1]], dtype=numpy.double)

    t.tilts_x = numpy.asarray(f['antenna0']['tracker']['tilt_xy_avg'][0], dtype=numpy.double)
    t.tilts_y = numpy.asarray(f['antenna0']['tracker']['tilt_xy_avg'][1], dtype=numpy.double)
    t.refraction = numpy.asarray(f['antenna0']['tracker']['refraction'][2], dtype=numpy.double)

    t.horiz_mount_x = numpy.asarray(f['antenna0']['tracker']['horiz_mount'][0])
    t.horiz_mount_y = numpy.asarray(f['antenna0']['tracker']['horiz_mount'][1])
    t.horiz_off_x = numpy.asarray(f['antenna0']['tracker']['horiz_off'][0])
    t.horiz_off_y = numpy.asarray(f['antenna0']['tracker']['horiz_off'][1])

    t.scan_off_x = numpy.asarray(f['antenna0']['tracker']['scan_off'][0])
    t.scan_off_y = numpy.asarray(f['antenna0']['tracker']['scan_off'][1])
    t.sky_off_x = numpy.asarray(f['antenna0']['tracker']['sky_xy_off'][0])
    t.sky_off_y = numpy.asarray(f['antenna0']['tracker']['sky_xy_off'][1])
    t.equat_off_x = numpy.asarray(f['antenna0']['tracker']['equat_off'][0])
    t.equat_off_y = numpy.asarray(f['antenna0']['tracker']['equat_off'][1])

    t.equat_geoc_ra = numpy.asarray(f['antenna0']['tracker']['equat_geoc'][0])
    t.equat_geoc_dec = numpy.asarray(f['antenna0']['tracker']['equat_geoc'][1])
    t.horiz_topo_az = numpy.asarray(f['antenna0']['tracker']['horiz_topo'][0])
    t.horiz_topo_el = numpy.asarray(f['antenna0']['tracker']['horiz_topo'][1])

    t.error_az = numpy.asarray(f['antenna0']['tracker']['errors'][0])
    t.error_el = numpy.asarray(f['antenna0']['tracker']['errors'][1])

    t.linsens_avg_l1 = numpy.asarray(f['antenna0']['tracker']['linear_sensor_avg'][0])
    t.linsens_avg_l2 = numpy.asarray(f['antenna0']['tracker']['linear_sensor_avg'][1])
    t.linsens_avg_r1 = numpy.asarray(f['antenna0']['tracker']['linear_sensor_avg'][2])
    t.linsens_avg_r2 = numpy.asarray(f['antenna0']['tracker']['linear_sensor_avg'][3])
    
    t.telescope_temp = numpy.asarray([f['array']['weather']['airTemperature'].value])
    t.telescope_pressure = numpy.asarray([f['array']['weather']['pressure'].value])

    f['TrackerPointing'] = t

    p = core.G3MapVectorDouble()
    p['tilts'] = numpy.asarray(f['antenna0']['tracker']['tilts'], dtype=numpy.double)
    p['flexure'] = numpy.asarray(f['antenna0']['tracker']['flexure'], dtype=numpy.double)
    p['fixedCollimation'] = numpy.asarray(f['antenna0']['tracker']['fixedCollimation'], dtype=numpy.double)
    p['time'] = numpy.asarray(t.time, dtype=numpy.double)

    f['OnlinePointingModel'] = p

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

    benchcom = core.G3TimestreamMap()
    benchpos = core.G3TimestreamMap()
    benchzero = core.G3TimestreamMap()
    benchoff = core.G3TimestreamMap()
    bencherr = core.G3TimestreamMap()
    bench_info = core.G3TimestreamMap()
    for i, key in enumerate(bench_axes):
        # As of 2017-08-03, SCU time is not trustworthy
        # start = f['antenna0']['scu']['benchSampleTime'][0][0]
        # stop = f['antenna0']['scu']['benchSampleTime'][0][-1]
        # For now, do this bit of evil
        start = f['antenna0']['tracker']['utc'][0][0]
        stop = f['antenna0']['tracker']['utc'][0][-1]

        benchcom[key] = core.G3Timestream(f['antenna0']['scu']['benchExpected'][i])
        benchcom[key].start = start
        benchcom[key].stop = stop

        benchpos[key] = core.G3Timestream(f['antenna0']['scu']['benchActual'][i])
        benchpos[key].start = start
        benchpos[key].stop = stop

        benchzero[key] = core.G3Timestream(f['antenna0']['scu']['benchZeros'][i])
        benchzero[key].start = start
        benchzero[key].stop = stop

        benchoff[key] = core.G3Timestream(f['antenna0']['scu']['benchOffsets'][i])
        benchoff[key].start = start
        benchoff[key].stop = stop

        bencherr[key] = core.G3Timestream(f['antenna0']['scu']['benchErrors'][i])
        bencherr[key].start = start
        bencherr[key].stop = stop

    info_items = ['benchFocus', 'benchDeadBand', 'benchAcquiredThreshold',
                  'benchPrimaryState', 'benchSecondaryState', 
                  'benchFault', 'timeLocked']
    bench_info = core.G3TimestreamMap()
    for i, key in enumerate(info_items):
        start = f['antenna0']['tracker']['utc'][0][0]
        stop = f['antenna0']['tracker']['utc'][0][-1]

        bench_info[key] = core.G3Timestream(f['antenna0']['scu'][key][0])
        bench_info[key].start = start
        bench_info[key].stop = stop

    f['BenchPosition'] = benchpos
    f['BenchCommandedPosition'] = benchcom
    f['BenchZeros'] = benchzero
    f['BenchOffsets'] = benchoff
    f['BenchErrors'] = bencherr
    f['BenchInfo'] = bench_info
    f['BenchSampleTime'] = f['antenna0']['scu']['benchSampleTime'][0]
    
@core.indexmod
def UnpackCryoData(f):
    '''
    Extracts cryo information into CryoStatus key
    '''

    if f.type != core.G3FrameType.GcpSlow:
        return

    t = core.G3MapDouble()
    t.time = f['array']['cryo']['utc']
    t.cryo_is_valid = f['array']['cryo']['cryoIsValid'][0]

    # Measured values
    # He10
    t.uc_head = f['array']['cryo']['temperature'][0][0]
    t.ic_head = f['array']['cryo']['temperature'][0][1]
    t.he4_head = f['array']['cryo']['temperature'][0][2]
    t.he4_fb = f['array']['cryo']['temperature'][0][3]
    t.he4_pump = f['array']['cryo']['temperature'][0][4]
    t.ic_pump = f['array']['cryo']['temperature'][0][5]
    t.uc_pump = f['array']['cryo']['temperature'][0][6]
    t.he4_sw = f['array']['cryo']['temperature'][0][7]
    t.ic_sw = f['array']['cryo']['temperature'][0][8]
    t.uc_sw = f['array']['cryo']['temperature'][0][9]
    t.uc_stage = f['array']['cryo']['temperature'][0][10]
    t.lc_tower = f['array']['cryo']['temperature'][0][11]
    t.ic_stage = f['array']['cryo']['temperature'][0][12]
    t.t4k_head = f['array']['cryo']['temperature'][0][13]
    t.t4k_squid_strap = f['array']['cryo']['temperature'][0][14]
    t.t50k_head = f['array']['cryo']['temperature'][0][15]

    # Optics
    t.b1_50k_wbp_near = f['array']['cryo']['temperature'][1][0]
    t.b2_50k_wbp_far = f['array']['cryo']['temperature'][1][1]
    t.b3_50k_diving_board = f['array']['cryo']['temperature'][1][2]
    t.b4_50k_top_bot_ptc = f['array']['cryo']['temperature'][1][3]
    t.y1_50k_head = f['array']['cryo']['temperature'][1][4]
    t.y2_50k_window_strap_near = f['array']['cryo']['temperature'][1][5]
    t.y3_50k_tube_strap_near = f['array']['cryo']['temperature'][1][6]
    t.y4_50k_tube = f['array']['cryo']['temperature'][1][7]
    t.g1_4k_head = f['array']['cryo']['temperature'][1][8]
    t.g2_4k_strap = f['array']['cryo']['temperature'][1][9]
    t.g3_4k_lens_tab = f['array']['cryo']['temperature'][1][10]
    t.g4_4k_lens_tab_far = f['array']['cryo']['temperature'][1][11]
    t.r1_4k_top_top_ptc = f['array']['cryo']['temperature'][1][12]
    t.r2_50k_midop_bot_ptc = f['array']['cryo']['temperature'][1][13]
    t.r3_4k_lyot_flange = f['array']['cryo']['temperature'][1][14]
    t.r4_4k_lyot = f['array']['cryo']['temperature'][1][15]

    # Receiver
    t.t4k_plate_far = f['array']['cryo']['temperature'][2][0]
    t.t4k_strap_optics = f['array']['cryo']['temperature'][2][1]
    t.t4k_plate_mid = f['array']['cryo']['temperature'][2][2]
    t.t4k_plate_top = f['array']['cryo']['temperature'][2][3]
    t.t4k_plate_ptc = f['array']['cryo']['temperature'][2][4]
    t.t50k_harness_middle = f['array']['cryo']['temperature'][2][6]
    t.t50k_strap = f['array']['cryo']['temperature'][2][7]
    t.squid_wh1_sl1 = f['array']['cryo']['temperature'][2][8]
    t.squid_wh5_sl1 = f['array']['cryo']['temperature'][2][9]
    t.squid_wh3_sl7 = f['array']['cryo']['temperature'][2][10]
    t.cal_filament = f['array']['cryo']['temperature'][2][11]
    t.cal_ambient1 = f['array']['cryo']['temperature'][2][12]
    t.cal_ambient2 = f['array']['cryo']['temperature'][2][13]
    t.cal_ambient3 = f['array']['cryo']['temperature'][2][14]

    # Heaters
    t.heat_he4_pump = f['array']['cryo']['heater_dac'][0][3]
    t.heat_ic_pump = f['array']['cryo']['heater_dac'][0][4]
    t.heat_uc_pump = f['array']['cryo']['heater_dac'][0][5]
    t.heat_he4_sw = f['array']['cryo']['heater_dac'][0][0]
    t.heat_ic_sw = f['array']['cryo']['heater_dac'][0][1]
    t.heat_uc_sw= f['array']['cryo']['heater_dac'][0][2]

    f['CryoStatus'] = t


@core.indexmod
def UnpackPTData(f):
    '''Extracts pulse tube status information to PTStatus key 
    in frame'''

    if f.type != core.G3FrameType.GcpSlow:
        return

    p = core.G3MapDouble()

    p.time = f['array']['pt415']['utc']
    p.optics_lowp = f['array']['pt415']['pressure_low'][0]
    p.min_optics_lowp = f['array']['pt415']['min_pressure_low'][0]
    p.max_optics_lowp = f['array']['pt415']['max_pressure_low'][0]
    p.optics_highp = f['array']['pt415']['pressure_high'][0]
    p.min_optics_highp = f['array']['pt415']['min_pressure_high'][0]
    p.max_optics_highp = f['array']['pt415']['max_pressure_high'][0]
    p.optics_tempoil = f['array']['pt415']['temp_oil'][0]
    p.min_optics_tempoil = f['array']['pt415']['min_temp_oil'][0]
    p.max_optics_tempoil = f['array']['pt415']['max_temp_oil'][0]

    p.receiver_lowp = f['array']['pt415']['pressure_low'][1]
    p.min_receiver_lowp = f['array']['pt415']['min_pressure_low'][1]
    p.max_receiver_lowp = f['array']['pt415']['max_pressure_low'][1]
    p.receiver_highp = f['array']['pt415']['pressure_high'][1]
    p.min_receiver_highp = f['array']['pt415']['min_pressure_high'][1]
    p.max_receiver_highp = f['array']['pt415']['max_pressure_high'][1]
    p.receiver_tempoil = f['array']['pt415']['temp_oil'][1]
    p.min_receiver_tempoil = f['array']['pt415']['min_temp_oil'][1]
    p.max_receiver_tempoil = f['array']['pt415']['max_temp_oil'][1]

    p.optics_is_valid = f['array']['pt415']['deviceIsValid'][0]
    p.receiver_is_valid = f['array']['pt415']['deviceIsValid'][1]


    f['PTStatus'] = p

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
    fpga_temp.time = mux['utc']
    board_name.time = mux['utc']
    f['MuxFPGATemp'] = fpga_temp
    f['MuxBoardName'] = board_name

@core.indexmod
def UnpackWeatherData(f):
    '''Extracts weather status information to Weather key 
    in frame'''

    if f.type != core.G3FrameType.GcpSlow:
        return

    t = core.G3MapDouble()
    t.time = f['array']['weather']['utc']
    t.telescope_temp = f['array']['weather']['airTemperature'].value
    t.telescope_pressure = f['array']['weather']['pressure'].value
    t.inside_dsl_temp = f['array']['weather']['internalTemperature'].value
    t.wind_speed = f['array']['weather']['windSpeed'].value
    t.wind_direction = f['array']['weather']['windDirection'].value
    t.battery = f['array']['weather']['battery'].value
    t.rel_humidity = f['array']['weather']['relativeHumidity'].value
    t.power = f['array']['weather']['power'].value
    t.tau = f['array']['tipper']['tau'].value
    t.tatm = f['array']['tipper']['tatm'].value

    f['Weather'] = t

@core.pipesegment
def ARCExtract(pipe):
    '''Extract GCP registers into native objects'''
    pipe.Add(CalibrateFrame)
    pipe.Add(UnpackACUData)
    pipe.Add(UnpackTrackerPointingData)
    pipe.Add(DecryptFeatureBit)
    pipe.Add(UnpackTrackerData)
    pipe.Add(AddBenchData)
    pipe.Add(UnpackCryoData)
    pipe.Add(UnpackPTData)
    pipe.Add(UnpackMuxData)
    pipe.Add(UnpackWeatherData)

@core.pipesegment
def ARCExtractMinimal(pipe):
    '''Extract bare minimum GCP registers into native objects.

    Includes only GCPFeatureBits, SourceName and ObservationID keys.
    Use ARCExtract to calibrate and extract the complete frame.
    '''
    pipe.Add(DecryptFeatureBit)
    pipe.Add(UnpackTrackerMinimal)

# Need tool for tilt meter next
