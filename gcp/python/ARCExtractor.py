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
    a.acu_status = f['antenna0']['acu']['acu_status'].value

    f['ACUStatus'] = a

@core.indexmod
def UnpackTrackerData(f, rewrite_source_from_feature_bits=True):
    '''
    Extracts tracker status information to frame. If
    rewrite_source_from_feature_bits is True (the default), will try to
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
    t.scan_flag = core.BoolVector(f['antenna0']['tracker']['scan_flag'][0])

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
    t.tilts_x = numpy.asarray(f['antenna0']['tracker']['tilt_xy_avg'][0], dtype=numpy.double)
    t.tilts_y = numpy.asarray(f['antenna0']['tracker']['tilt_xy_avg'][1], dtype=numpy.double)
    t.refraction = numpy.asarray(f['antenna0']['tracker']['refraction'][2], dtype=numpy.double)

    t.horiz_mount_x = numpy.asarray(f['antenna0']['tracker']['horiz_mount'][0])
    t.horiz_mount_y = numpy.asarray(f['antenna0']['tracker']['horiz_mount'][1])
    t.horiz_off_x = numpy.asarray(f['antenna0']['tracker']['horiz_off'][0])
    t.horiz_off_y = numpy.asarray(f['antenna0']['tracker']['horiz_off'][1])

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

    f['BenchPosition'] = benchpos
    f['BenchCommandedPosition'] = benchcom
    f['BenchZeros'] = benchzero
    

@core.pipesegment
def ARCExtract(pipe):
    '''Extract GCP registers into native objects'''
    pipe.Add(CalibrateFrame)
    pipe.Add(UnpackACUData)
    pipe.Add(UnpackTrackerPointingData)
    pipe.Add(DecryptFeatureBit)
    pipe.Add(UnpackTrackerData)
    pipe.Add(AddBenchData)

# Need tool for tilt meter next
