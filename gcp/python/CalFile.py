from spt3g import core
import os, time
import numpy as np

@core.indexmod
class CalFileReader(object):
    '''
    For now just reads a G3 calibration file and loads the contents into 
    a dictionary
    '''

    def __init__(self, calibration_file=None):
        self.default_calibration_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'G3RegCal.txt')

    def readCalFile(self, calfile=None):
        cal_dict = {}
        if calfile is None:
            self.calibration_file = self.default_calibration_file
        else:
            self.calibration_file = calfile
        try:
            f = open(self.calibration_file, 'r')
        except:
            core.log_warn('3G calibration file ' + 
                              self.calibration_file + ' not found.\n')
            return cal_dict
        ninvalid = 0
        for line in f:
            if line[0] != '#' and line[0] != '\n' and len(line) > 0:
                try:
                    line = line.replace('\n','')
                    fline = list(filter(None,line.split(' '))) #list() for py3
                    name = fline[0]
                    linedict = {}
                    linedict['Offset'] = float(fline[1])
                    linedict['ReciprocalFactor'] = float(fline[2])
                    try:
                        linedict['UnitName'] = fline[3]
                    except:
                        linedict['UnitName'] = 'None'
                    cal_dict[name] = linedict
                except:
                    ninvalid += 1
        if ninvalid > 0:
            if ninvalid == 1:
                core.log_warn('One invalid line found and ignored in calibration file.\n')
            else:
                core.log_warn(str(ninvalid) + ' invalid lines found and ignored in calibration file.\n')
        f.close()

        # now parse into dictionary structure NW uses in ARCExtractor
        cal_dict_orig = cal_dict.copy()
        cal_dict = {}
        origkeys = cal_dict_orig.keys()
        origkeys = sorted(origkeys)
        names = [key for key in origkeys]
        regmaps = [(list(filter(None,(tname).split('.'))))[0] for tname in names]
        uregmaps = np.unique(np.asarray(regmaps))
        for urm in uregmaps:
            cal_dict[urm] = {}
        regblocks = [(list(filter(None,(tname).split('.'))))[1] for tname in names]
        uregblocks,argurb = np.unique(np.asarray(regblocks),return_index=True)
        for j in np.arange(len(uregblocks)):
            cal_dict[regmaps[argurb[j]]][uregblocks[j]] = {}
        regs = [(list(filter(None,(tname).split('.'))))[2] for tname in names]
        for j in np.arange(len(regs)):
            reg = regs[j]
            if '[' in reg:
                reg_plus_indices = reg.split('[')
                realreg = reg_plus_indices[0]
                indices = ((reg_plus_indices[1]).split(']'))[0]
                if '-' in indices:
                    cal_dict[regmaps[j]][regblocks[j]][realreg] = {}
                    siindices = indices.split('-')
                    iindices = (np.array(siindices)).astype(int)
                    for i in np.arange(np.max(iindices)-np.min(iindices)+1) + \
                            np.min(iindices):
                        cal_dict[regmaps[j]][regblocks[j]]\
                            [realreg][i] = cal_dict_orig[names[j]]
                else:
                    if indices == '0':
                        cal_dict[regmaps[j]][regblocks[j]][realreg] = {}
                    cal_dict[regmaps[j]][regblocks[j]][realreg]\
                        [int(indices)] = cal_dict_orig[names[j]]
            else:
                cal_dict[regmaps[j]][regblocks[j]][reg] = \
                    cal_dict_orig[names[j]]

        return cal_dict

@core.usefulfunc
def create_g3_cal_file(path, read_from_gcp=True, extra_dict=None,
                       use_extra_info=True, gcp_cal_file=None):
    '''
    Create a G3 register calibration file. Usually reads in
    calibration and units from a GCP cal file then adds extra
    information that isn't in the GCP one. This extra information can
    be handed to the routine as a dict (note format below); otherwise,
    the hard-coded extra info will be used. (Feel free to edit this
    hard-coded info.)
    '''

    cal_dict = {}
    status = 0
    
    if read_from_gcp:
        try:
            default_calibration_file = \
                os.environ['GCP_DIR']+'/config/init/cal'
        except:
            default_calibration_file = \
                '/home/sptdat/gcproot/gcp/config/init/cal'
        if gcp_cal_file is None:
            calibration_file = default_calibration_file
        else:
            calibration_file = gcp_cal_file
        try:
            f = open(calibration_file)
        except:
            core.log_warn('GCP calibration file ' + 
                          calibration_file + ' not found.\n')
            f = []
        ninvalid = 0
        for line in f:
            if line[0] != '#' and line[0] != '\n' and len(line) > 0:
                line = line.replace('\n','')
                linedict = {}
                # cal file has spaces AND tabs (and comment characters), whee!
                info_and_comment = list(filter(None,line.split('#')))
                fline = list(filter(None,(info_and_comment[0]).split(' ')))
                if len(fline) == 1 and '\t' in fline[0]:
                    fline = list(filter(None,(info_and_comment[0]).split('\t')))
                name = (fline[0]).replace('*','0')
                try:
                    linedict['Offset'] = float(fline[1])
                    if '/' in fline[2]:
                        gainfacs = (fline[2]).split('/')
                        gainfac = float(gainfacs[1])/float(gainfacs[0])
                    else:
                        gainfac = 1./float(fline[2])
                    linedict['ReciprocalFactor'] = gainfac
                    linedict['UnitName'] = ''
                # try to figure out units. currently cal file 
                # only has 2 forms, but this is not robust.
                    if len(info_and_comment) > 1: 
                        comment = info_and_comment[1]
                        if '->' in comment:
                            comments = list(filter(None,(comment).split('->')))
                            unitname = \
                                list(filter(None,(comments[1]).split(' ')))[0]
                            unitname = list(filter(None,(unitname).split('/')))[0]
                        else:
                            comments = \
                                list(filter(None,(comment).split('display in ')))
                            unitname = \
                                list(filter(None,(comments[0]).split(' ')))[0]
                            unitname = unitname.replace('\n','')
                        linedict['UnitName'] = unitname
                # another hudge kludge
                        if 'rate' in name:
                            if '/' not in unitname:
                                unitname = unitname + '/s'
                                linedict['UnitName'] = unitname
                    cal_dict[name] = linedict
                except:
                    ninvalid += 1
        if ninvalid > 0:
            if ninvalid == 1:
                core.log_warn('One invalid line found and ignored in calibration file.\n')
            else:
                core.log_warn(str(ninvalid) + ' invalid lines found and ignored in calibration file.\n')
        f.close()

    # now add extra info
    if use_extra_info:
        if extra_dict is None:
            cal_dict['antenna0.acu.az_pos'] = \
                {'ReciprocalFactor':1.,'Offset':0,'UnitName':'deg'}
            cal_dict['antenna0.acu.el_pos'] = \
                {'ReciprocalFactor':1.,'Offset':0,'UnitName':'deg'}
            cal_dict['antenna0.acu.az_err'] = \
                {'ReciprocalFactor':1.,'Offset':0,'UnitName':'deg'}
            cal_dict['antenna0.acu.el_err'] = \
                {'ReciprocalFactor':1.,'Offset':0,'UnitName':'deg'}

            cal_dict['antenna0.scu.temp'] = \
                {'ReciprocalFactor':1.,'Offset':273.15,'UnitName':'K'}
            cal_dict['array.weather.airTemperature'] = \
                {'ReciprocalFactor':1.,'Offset':273.15,'UnitName':'K'}
        else:
            try:
                for key in extra_dict:
                    cal_dict[key] = {}
                    cal_dict[key]['Offset'] = extra_dict[key]['Offset']
                    cal_dict[key]['ReciprocalFactor'] = \
                        extra_dict[key]['ReciprocalFactor']
                    cal_dict[key]['UnitName'] = extra_dict[key]['UnitName']
            except:
                core.log_warn('extra_dict input format not recognized, ignoring')
                del cal_dict[key]

    # write cal file
    fcal = open(path,'w')
    fcal.write('#-----------------------------------------------------------------------\n')
    fcal.write('# This is an SPT-3G register calibration file. \n')
    fcal.write('#\n')
    fcal.write('# It was created on ' + time.ctime() + ' by user ' + os.getlogin() + '\n')
    fcal.write('#  with the spt3g_software function gcp.CalFile.create_g3_cal_file. \n')
    fcal.write('#\n')
    fcal.write('# Each line has the following format:\n')
    fcal.write('#\n')
    fcal.write('#   register offset 1/factor units\n')
    fcal.write('#\n')
    fcal.write('# (Note that factors are stored as reciprocals, because they are most \n')
    fcal.write('#  often one over a whole number.) \n')
    fcal.write('#\n')
    fcal.write('# The calibrated value of a register is given by:\n')
    fcal.write('#\n')
    fcal.write('#   reg_cal = factor * (offset + register_value) * core.G3Units.$UNIT$\n')
    fcal.write('#\n')
    fcal.write('#-----------------------------------------------------------------------\n')
    fcal.write('# Register      Offset       1/Factor      Units \n')
    fcal.write('#-----------------------------------------------------------------------\n')
    fcal.write('\n')
    ktemp = list(cal_dict.keys())
    ktemp.sort()
    for key in ktemp:
        kstr = (str(key)).ljust(42)
        gstr = '%e' % cal_dict[key]['ReciprocalFactor']
        ostr = '%10.4f' % cal_dict[key]['Offset']
        fcal.write(kstr + '  ' + ostr + '  ' + gstr + '  ' + cal_dict[key]['UnitName'] + '\n')
    fcal.close()

    return status

