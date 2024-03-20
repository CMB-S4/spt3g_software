import numpy as np
import re
import datetime as dt
from spt3g import core
from spt3g.core import G3Units as U
from . import ARCExtractor

def build_field_list(fr):
    """
    All the fields you can upload from gcp to InfluxDB
    """
    r = {
        # tracker status
        'az_actual': ['TrackerStatus', 'az_pos', U.deg],
        'el_actual': ['TrackerStatus', 'el_pos', U.deg],
        'az_rate_actual': ['TrackerStatus', 'az_rate', U.deg/U.sec],
        'el_rate_actual': ['TrackerStatus', 'el_rate', U.deg/U.sec],
        'az_command': ['TrackerStatus', 'az_command', U.deg],
        'el_command': ['TrackerStatus', 'el_command', U.deg],
        'az_rate_command': ['TrackerStatus', 'az_rate_command', U.deg/U.sec],
        'el_rate_command': ['TrackerStatus', 'el_rate_command', U.deg/U.sec],
        'tracker_state': ['TrackerStatus', 'state', None],
        'acu_seq': ['TrackerStatus', 'acu_seq', None],
        'in_control_int': ['TrackerStatus', 'in_control_int', None],
        'scan_flag': ['TrackerStatus', 'scan_flag', None],
        'lst': ['TrackerStatus', 'lst', U.hour],
        'source_acquired': ['TrackerStatus', 'source_acquired', None],
        'source_acquired_thresh': ['TrackerStatus', 'source_acquired_threshold', None],
        'tracker_mode': ['TrackerStatus', 'tracker_mode', None],
        'tracker_lacking': ['TrackerStatus', 'tracker_lacking', None],
        'time_status': ['TrackerStatus', 'time_status', None],
        'schedule': ['TrackerStatus', 'schedule_name', None],
        'raw_encoder_1': ['antenna0', 'tracker', 'raw_encoder', 0, U.deg],
        'raw_encoder_2': ['antenna0', 'tracker', 'raw_encoder', 1, U.deg],
        'drive_currents_el1': ['array', 'dc', 'currents', 0, U.volt],
        'drive_currents_el2': ['array', 'dc', 'currents', 1, U.volt],
        'drive_currents_el3': ['array', 'dc', 'currents', 2, U.volt],
        'drive_currents_el4': ['array', 'dc', 'currents', 3, U.volt],
        'drive_currents_az1': ['array', 'dc', 'currents', 4, U.volt],
        'drive_currents_az2': ['array', 'dc', 'currents', 5, U.volt],
        'drive_currents_az3': ['array', 'dc', 'currents', 6, U.volt],
        'drive_currents_az4': ['array', 'dc', 'currents', 7, U.volt],

        # tracker pointing
        'features': ['TrackerPointing', 'features', 1],
        'encoder_off_x': ['TrackerPointing', 'encoder_off_x', U.deg],
        'encoder_off_y': ['TrackerPointing', 'encoder_off_y', U.deg],
        'low_limit_az': ['TrackerPointing', 'low_limit_az', U.deg],
        'high_limit_az': ['TrackerPointing', 'high_limit_az', U.deg],
        'low_limit_el': ['TrackerPointing', 'low_limit_el', U.deg],
        'high_limit_el': ['TrackerPointing', 'high_limit_el', U.deg],
        'tilts_x': ['TrackerPointing', 'tilts_x', U.deg],
        'tilts_y': ['TrackerPointing', 'tilts_y', U.deg],
        'refraction': ['TrackerPointing', 'refraction', U.deg],
        'horiz_mount_x': ['TrackerPointing', 'horiz_mount_x', U.deg],
        'horiz_mount_y': ['TrackerPointing', 'horiz_mount_y', U.deg],
        'horiz_topo_az': ['TrackerPointing', 'horiz_topo_az', U.deg],
        'horiz_topo_el': ['TrackerPointing', 'horiz_topo_el', U.deg],
        'horiz_off_x': ['TrackerPointing', 'horiz_off_x', U.deg],
        'horiz_off_y': ['TrackerPointing', 'horiz_off_y', U.deg],
        'scan_off_x': ['TrackerPointing', 'scan_off_x', U.deg],
        'scan_off_y': ['TrackerPointing', 'scan_off_y', U.deg],
        'sky_off_x': ['TrackerPointing', 'sky_off_x', U.deg],
        'sky_off_y': ['TrackerPointing', 'sky_off_y', U.deg],
        'equat_off_x': ['TrackerPointing', 'equat_off_x', U.deg],
        'equat_off_y': ['TrackerPointing', 'equat_off_y', U.deg],
        'source_ra': ['TrackerPointing', 'equat_geoc_ra', U.rahr],
        'source_dec': ['TrackerPointing', 'equat_geoc_dec', U.deg],
        'error_az': ['TrackerPointing', 'error_az', U.deg],
        'error_el': ['TrackerPointing', 'error_el', U.deg],
        'linsens_avg_l1': ['TrackerPointing', 'linsens_avg_l1', U.mm],
        'linsens_avg_l2': ['TrackerPointing', 'linsens_avg_l2', U.mm],
        'linsens_avg_r1': ['TrackerPointing', 'linsens_avg_r1', U.mm],
        'linsens_avg_r2': ['TrackerPointing', 'linsens_avg_r2', U.mm],
        'linsens_daz': ['LinearSensorDeltas', 'delta_az', U.deg],
        'linsens_del': ['LinearSensorDeltas', 'delta_el', U.deg],
        'linsens_det': ['LinearSensorDeltas', 'delta_et', U.deg],

        # Weather
        'telescope_temp': ['Weather', 'telescope_temp', 'C'],
        'inside_dsl_temp': ['Weather', 'inside_dsl_temp', 'C'],
        'telescope_pressure': ['Weather', 'telescope_pressure', U.mb],
        'wind_speed': ['Weather', 'wind_speed', U.m / U.s],
        'wind_direction': ['Weather', 'wind_direction', U.deg],
        'battery': ['Weather', 'battery', U.mV],
        'rel_humidity': ['Weather', 'rel_humidity', None],
        'power': ['Weather', 'power', None],
        'tau': ['Weather', 'tau', None],
        'tatm': ['Weather', 'tatm', None],

        # Cryo -- units appear to just be in K. Don't recalibrate.
        # He10
        'uc_head': ['CryoStatus', 'uc_head', 1],
        'ic_head': ['CryoStatus', 'ic_head', 1],
        'he4_head': ['CryoStatus', 'he4_head', 1],
        'he4_fb': ['CryoStatus', 'he4_fb', 1],
        'he4_pump': ['CryoStatus', 'he4_pump', 1],
        'ic_pump': ['CryoStatus', 'ic_pump', 1],
        'uc_pump': ['CryoStatus', 'uc_pump', 1],
        'he4_sw': ['CryoStatus', 'he4_sw', 1],
        'ic_sw': ['CryoStatus', 'ic_sw', 1],
        'uc_sw': ['CryoStatus', 'uc_sw', 1],
        'uc_stage': ['CryoStatus', 'uc_stage', 1],
        'lc_tower': ['CryoStatus', 'lc_tower', 1],
        'ic_stage': ['CryoStatus', 'ic_stage', 1],
        '4k_head': ['CryoStatus', 't4k_head', 1],
        '4k_squid_strap': ['CryoStatus', 't4k_squid_strap', 1],
        '50k_head': ['CryoStatus', 't50k_head', 1],
        # Optics
        'b1_50k_wbp_near': ['CryoStatus', 'b1_50k_wbp_near', 1],
        'b2_50k_wbp_far': ['CryoStatus', 'b2_50k_wbp_far', 1],
        'b3_50k_diving_board': ['CryoStatus', 'b3_50k_diving_board', 1],
        'b4_50k_top_bot_ptc': ['CryoStatus', 'b4_50k_top_bot_ptc', 1],
        'y1_50k_head': ['CryoStatus', 'y1_50k_head', 1],
        'y2_50k_window_strap_near': ['CryoStatus', 'y2_50k_window_strap_near', 1],
        'y3_50k_tube_strap_near': ['CryoStatus', 'y3_50k_tube_strap_near', 1],
        'y4_50k_tube': ['CryoStatus', 'y4_50k_tube', 1],
        'g1_4k_head': ['CryoStatus', 'g1_4k_head', 1],
        'g2_4k_strap': ['CryoStatus', 'g2_4k_strap', 1],
        'g3_4k_lens_tab': ['CryoStatus', 'g3_4k_lens_tab', 1],
        'g4_4k_lens_tab_far': ['CryoStatus', 'g4_4k_lens_tab_far', 1],
        'r1_4k_top_top_ptc': ['CryoStatus', 'r1_4k_top_top_ptc', 1],
        'r2_50k_midop_bot_ptc': ['CryoStatus', 'r2_50k_midop_bot_ptc', 1],
        'r3_4k_lyot_flange': ['CryoStatus', 'r3_4k_lyot_flange', 1],
        'r4_4k_lyot': ['CryoStatus', 'r4_4k_lyot', 1],
        # Receiver
        '4k_plate_far': ['CryoStatus', 't4k_plate_far', 1],
        '4k_strap_optics': ['CryoStatus', 't4k_strap_optics', 1],
        '4k_plate_mid': ['CryoStatus', 't4k_plate_mid', 1],
        '4k_plate_top': ['CryoStatus', 't4k_plate_top', 1],
        '4k_plate_ptc': ['CryoStatus', 't4k_plate_ptc', 1],
        '50k_harness_middle': ['CryoStatus', 't50k_harness_middle', 1],
        '50k_strap': ['CryoStatus', 't50k_strap', 1],
        'squid_wh1_sl1': ['CryoStatus', 'squid_wh1_sl1', 1],
        'squid_wh5_sl1': ['CryoStatus', 'squid_wh5_sl1', 1],
        'squid_wh3_sl7': ['CryoStatus', 'squid_wh3_sl7', 1],
        'cal_filament': ['CryoStatus', 'cal_filament', 1],
        'cal_ambient1': ['CryoStatus', 'cal_ambient1', 1],
        'cal_ambient2': ['CryoStatus', 'cal_ambient2', 1],
        'cal_ambient3': ['CryoStatus', 'cal_ambient3', 1],
        # heaters
        'heat_he4_pump': ['CryoStatus', 'heat_he4_pump', 1],
        'heat_ic_pump': ['CryoStatus', 'heat_ic_pump', 1],
        'heat_uc_pump': ['CryoStatus', 'heat_uc_pump', 1],
        'heat_he4_sw': ['CryoStatus', 'heat_he4_sw', 1],
        'heat_ic_sw': ['CryoStatus', 'heat_ic_sw', 1],
        'heat_uc_sw': ['CryoStatus', 'heat_uc_sw', 1],
        # status bit
        'cryo_is_valid': ['CryoStatus', 'cryo_is_valid', None],

        # PT status
        'optics_low_p_now': ['PTStatus', 'optics_lowp', None],
        'optics_low_p_min': ['PTStatus', 'min_optics_lowp', None],
        'optics_low_p_max': ['PTStatus', 'max_optics_lowp', None],
        'optics_high_p_now': ['PTStatus', 'optics_highp', None],
        'optics_high_p_min': ['PTStatus', 'min_optics_highp', None],
        'optics_high_p_max': ['PTStatus', 'max_optics_highp', None],
        'optics_tempoil_now': ['PTStatus', 'optics_tempoil', None],
        'optics_tempoil_min': ['PTStatus', 'min_optics_tempoil', None],
        'optics_tempoil_max': ['PTStatus', 'max_optics_tempoil', None],

        'receiver_low_p_now': ['PTStatus', 'receiver_lowp', None],
        'receiver_low_p_min': ['PTStatus', 'min_receiver_lowp', None],
        'receiver_low_p_max': ['PTStatus', 'max_receiver_lowp', None],
        'receiver_high_p_now': ['PTStatus', 'receiver_highp', None],
        'receiver_high_p_min': ['PTStatus', 'min_receiver_highp', None],
        'receiver_high_p_max': ['PTStatus', 'max_receiver_highp', None],
        'receiver_tempoil_now': ['PTStatus', 'receiver_tempoil', None],
        'receiver_tempoil_min': ['PTStatus', 'min_receiver_tempoil', None],
        'receiver_tempoil_max': ['PTStatus', 'max_receiver_tempoil', None],

        'optics_is_valid': ['PTStatus', 'optics_is_valid', None],
        'receiver_is_valid': ['PTStatus', 'receiver_is_valid', None],

        # Online Pointing Model
        'tilts_hr_angle': ['OnlinePointingModel', 'tilts', 0, U.deg],
        'tilts_lat': ['OnlinePointingModel', 'tilts', 1, U.deg],
        'tilts_el': ['OnlinePointingModel', 'tilts', 2, U.deg],
        'flexure_sin': ['OnlinePointingModel', 'flexure', 0, U.deg],
        'flexure_cos': ['OnlinePointingModel', 'flexure', 1, U.deg],
        'fixed_collimation_x': ['OnlinePointingModel', 'fixedCollimation', 0, U.deg],
        'fixed_collimation_y': ['OnlinePointingModel', 'fixedCollimation', 1, U.deg],
        'tilts_hr_angle_adjust': ['OnlinePointingModelCorrection', 'tilts', 0, U.deg],
        'tilts_lat_adjust': ['OnlinePointingModelCorrection', 'tilts', 1, U.deg],
        'tilts_el_adjust': ['OnlinePointingModelCorrection', 'tilts', 2, U.deg],
        'flexure_sin_adjust': ['OnlinePointingModelCorrection', 'flexure', 0, U.deg],
        'flexure_cos_adjust': ['OnlinePointingModelCorrection', 'flexure', 1, U.deg],
        'fixed_collimation_x_adjust': ['OnlinePointingModelCorrection', 'fixedCollimation', 0, U.deg],
        'fixed_collimation_y_adjust': ['OnlinePointingModelCorrection', 'fixedCollimation', 1, U.deg],
        'linsens_coeff_az': ['OnlinePointingModel', 'linsensCoeffs', 0, None],
        'linsens_coeff_el': ['OnlinePointingModel', 'linsensCoeffs', 1, None],
        'linsens_coeff_et': ['OnlinePointingModel', 'linsensCoeffs', 2, None],
        'linsens_enabled': ['OnlinePointingModel', 'linsensEnabled', 0, None],

        # Other
        'obs_id': ['ObservationID', None],
        'source_name': ['SourceName', None],

        # ACUStatus
        'acu_state': ['ACUStatus', 'state', None],
        'acu_status': ['ACUStatus', 'status', None],
        'acu_error': ['ACUStatus', 'error', None],

        # Bench
        'bench_command_y1': ['BenchCommandedPosition', 'y1', U.mm],
        'bench_command_y2': ['BenchCommandedPosition', 'y2', U.mm],
        'bench_command_y3': ['BenchCommandedPosition', 'y3', U.mm],
        'bench_command_x4': ['BenchCommandedPosition', 'x4', U.mm],
        'bench_command_x5': ['BenchCommandedPosition', 'x5', U.mm],
        'bench_command_z6': ['BenchCommandedPosition', 'z6', U.mm],

        'bench_actual_y1': ['BenchPosition', 'y1', U.mm],
        'bench_actual_y2': ['BenchPosition', 'y2', U.mm],
        'bench_actual_y3': ['BenchPosition', 'y3', U.mm],
        'bench_actual_x4': ['BenchPosition', 'x4', U.mm],
        'bench_actual_x5': ['BenchPosition', 'x5', U.mm],
        'bench_actual_z6': ['BenchPosition', 'z6', U.mm],

        'bench_zero_y1': ['BenchZeros', 'y1', U.mm],
        'bench_zero_y2': ['BenchZeros', 'y2', U.mm],
        'bench_zero_y3': ['BenchZeros', 'y3', U.mm],
        'bench_zero_x4': ['BenchZeros', 'x4', U.mm],
        'bench_zero_x5': ['BenchZeros', 'x5', U.mm],
        'bench_zero_z6': ['BenchZeros', 'z6', U.mm],

        'bench_offset_y1': ['BenchOffsets', 'y1', U.mm],
        'bench_offset_y2': ['BenchOffsets', 'y2', U.mm],
        'bench_offset_y3': ['BenchOffsets', 'y3', U.mm],
        'bench_offset_x4': ['BenchOffsets', 'x4', U.mm],
        'bench_offset_x5': ['BenchOffsets', 'x5', U.mm],
        'bench_offset_z6': ['BenchOffsets', 'z6', U.mm],

        'bench_error_y1': ['BenchErrors', 'y1', U.mm],
        'bench_error_y2': ['BenchErrors', 'y2', U.mm],
        'bench_error_y3': ['BenchErrors', 'y3', U.mm],
        'bench_error_x4': ['BenchErrors', 'x4', U.mm],
        'bench_error_x5': ['BenchErrors', 'x5', U.mm],
        'bench_error_z6': ['BenchErrors', 'z6', U.mm],

        'bench_focus': ['BenchInfo', 'benchFocus', U.mm],
        'bench_dead_band': ['BenchInfo', 'benchDeadBand', U.mm],
        'bench_acquired_thresh': ['BenchInfo', 'benchAcquiredThreshold', U.mm],
        'bench_primary_state': ['BenchInfo', 'benchPrimaryState', None],
        'bench_secondary_state': ['BenchInfo', 'benchSecondaryState', None],
        'bench_fault': ['BenchInfo', 'benchFault', None],
        'bench_time_locked': ['BenchInfo', 'timeLocked', None],
    }

    # mux housekeeping
    for i in range(32):
        i = str(i)
        r['fpga_temp_ib{}'.format(i)] = ['MuxFPGATemp', i, None]
        r['name_ib{}'.format(i)] = ['MuxBoardName', i, None]

    # scu.temp - all temps documented given a name, others just a number
    scu_temps = {
        0: 'yoke_air',
        1: 'ctrl_room_air',
        2: 'glycol_supply',
        3: 'glycol_return',
        4: 'ctrl_room',
        20: 'secondary',
        21: 'icecrate',
        22: 'bench',
        23: 'attic',
        24: 'cabin',
        25: 'cryoboard',
    }
    for i in range(60):
        key = 't_scu_{}'.format(scu_temps.get(i, i))
        r[key] = ['TrackerPointing', 'scu_temp', i, 'C']

    return r


def make_lines(measurement, field, time, dat, tags=None):
    '''
    Return list of string lines to add to database
    '''
    if isinstance(dat[0], str):
        dat = ['"' + x + '"' for x in dat]
    if tags is None:
        fmt_str = measuremnt +' '+ field + '={val} {time}'
    else:
        fmt_str = measurement
        for tag, tag_val in tags.items():
            fmt_str += ',{}={}'.format(tag, tag_val)
        fmt_str += ' ' + field + '={val} {time}'
    # time in int because it's in nanoseconds for InfluxDB
    lines = [fmt_str.format(val=x, time=int(t0)) for x, t0 in zip(dat, time)]
    return lines


@core.indexmod
def WriteDB(fr, client, fields=None):
    '''
    Write points to the database for each field

    Arguments
    ---------
    client :
        InfluxDB client
    fields :
        Which gcp fields to add to database. See parse_field for options. If
        None, add all.
    '''
    from influxdb.exceptions import InfluxDBClientError
    from influxdb.exceptions import InfluxDBServerError

    if fr.type != core.G3FrameType.GcpSlow:
        return
    all_fields = build_field_list(fr)
    if fields is None:
        fields = all_fields.keys()
    dict_list = []
    for f in fields:
        field_dat = all_fields[f]
        if len(field_dat) == 5:
            # raw register
            stat = 'TrackerStatus'
            tmp, brd, attr, ind, unit = field_dat
            try:
                dat = fr[tmp][brd][attr][ind]
            except:
                continue
            try:
                if 'utc' in fr[tmp][brd].keys():
                    time = fr[tmp][brd]['utc'][0]
                else:
                    time = fr['antenna0']['tracker']['utc'][0][:len(dat)]
            except:
                continue
        elif len(field_dat) == 4:
            stat, attr, ind, unit = field_dat
            if stat not in fr:
                # OnlinePointingModelCorrection
                continue
            try:
                dat = getattr(fr[stat], attr)[ind]
                time = getattr(fr[stat], 'time')
            except AttributeError:
                # OnlinePointingModel
                try:
                    dat = fr[stat][attr][ind]
                except KeyError:
                    # OnlinePointingModelCorrection
                    continue
                time = fr[stat]['time']
        elif len(field_dat) == 3:
            stat, attr, unit = field_dat
            if stat not in fr:
                # Field only exists in live data stream
                continue
            try:
                dat = getattr(fr[stat], attr)
            except AttributeError:
                try:
                    dat = fr[stat][attr]
                except KeyError: # Field only exists in live data stream
                    continue
            if 'Bench' in stat: # funny time field for bench positions
                time = fr['BenchSampleTime']
            elif 'Mux' in stat:
                time = fr['MuxTime']
            elif stat in ['CryoStatus', 'Weather', 'PTStatus']:
                time = fr['{}Time'.format(stat)]
            else:
                try:
                    time = getattr(fr[stat], 'time')
                except AttributeError:
                    time = fr[stat]['time']
        elif len(field_dat) == 2:
            stat, unit = field_dat
            try:
                dat = fr[stat]
            except KeyError: #eg, no obsid
                core.log_warn('No key {}'.format(stat), unit='InfluxDB')
                continue
            try:
                time = getattr(fr[stat], 'time')
            except AttributeError as err:
                time = fr['antenna0']['tracker']['utc'][0]

        # InfluxDB wants time in nanoseconds since the UNIX epoch in UTC
        if isinstance(time, core.G3Time):
            time = np.asarray([time.time / U.nanosecond])
        elif isinstance(time, core.G3VectorTime):
            time = np.asarray(time) / U.nanosecond
        else:
            try:
                time = np.asarray(core.G3VectorTime(np.atleast_1d(time))) / U.nanosecond
            except Exception as e:
                core.log_error("Error converting time: {}".format(str(e)), unit="InfluxDB")
                continue
        # Avoid out-of-range errors in parsing time
        # time < 0 occurs when, e.g., the SCU loses the clock (e.g. during a brownout)
        time[time < 0] = 0
        if dat is None:
            core.log_warn('{} dat is None'.format(f), unit='InfluxDB')
            continue
        dat = np.atleast_1d(dat)
        try:
            dlen = len(dat)
        except TypeError:
            # sometimes source_name is a weird non-none value
            continue
        if unit is not None:
            if unit == 'C':
                zeropt_K = 273.15
                cal_dat = dat / U.K - zeropt_K
            else:
                cal_dat = dat / unit
        else:
            cal_dat = dat
        try:
            if np.any(np.isnan(cal_dat)):
                continue
        except TypeError:
            pass
        if 'heat' not in f:
            tag = f
        else:
            tag = f.replace('heat_', '')

        # for fields that have az/el components
        az_el_names = ['az', 'el', 'az', 'el', 'ra', 'dec', 'x', 'y',
                       'hr_angle', 'sin', 'cos', 'lat']
        tag2 = f
        for name in az_el_names:
            # require name_ at beginning or _name at end
            match1 = re.findall('^{}_'.format(name), f)
            match2 = re.findall('_{}$'.format(name), f)
            if len(match1):
                tag2 = f.replace(match1[0], '')
            if len(match2):
                tag2 = f.replace(match2[0], '')
        # also group source names
        if 'source' in f:
            tag2 = 'source'
            stat = 'TrackerPointing'
        if stat == 'PTStatus':
            groups = ['now', 'min', 'max']
            for g in groups:
                match = re.findall('_{}$'.format(g), f)
                if len(match):
                    tag2 = f.replace(match[0], '')
        # group bench positions
        # require bench_ at beginning
        match = re.findall('^bench', f)
        if len(match):
            tag2 = attr # y1, y2, etc
            stat = 'Bench'

        # group Mux properties
        if 'Mux' in stat:
            stat = 'muxHousekeeping'
            tag2 = 'ib'+f.split('ib')[-1]

        if 'linsens_coeff' in f:
            tag2 = 'linsens_coeff'

        if 'drive_currents_az' in f:
            tag2 = 'drive_currents_az'
        if 'drive_currents_el' in f:
            tag2 = 'drive_currents_el'

        dict_list += make_lines(
            measurement=stat,
            field=f,
            time=time,
            dat=cal_dat,
            tags={'label': tag, 'label2': tag2},
        )

    try:
        now = core.G3Time.Now()
        delay = float(now.time/U.nanosecond - time[-1])/1e9
        if delay > 5:
            core.log_info(
                '{} Delay: {} s'.format(now.isoformat(), delay), unit='InfluxDB'
            )
    except RuntimeError: # sometimes timestamp gets screwed up
        pass

    try:
        client.write_points(
            dict_list, batch_size=len(dict_list), protocol='line'
        )
    except (InfluxDBClientError, InfluxDBServerError) as v:
        core.log_error(
            'Error writing to database. {}'.format(v), unit='InfluxDB'
        )


@core.pipesegment
def UpdateDB(pipe, client, fields=None):
    '''
    Update InfluxDB with data in frame

    Arguments
    ---------
    client :
        InfluxDB client
    fields :
        Which gcp fields to add to database. See parse_field for options. If
        None, add all.
    '''
    pipe.Add(ARCExtractor.ARCExtract)
    pipe.Add(WriteDB, client=client, fields=fields)
