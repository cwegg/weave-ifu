#!/usr/bin/env python

#
# Copyright (C) 2018 Cambridge Astronomical Survey Unit
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <https://www.gnu.org/licenses/>.
#


import xml.dom.minidom
import os
from astropy.io import fits
import numpy
import glob
import re
from operator import indexOf, countOf
from guidestar import GuideStar
from calibstar import CalibStar

class stage3:
    """
    Add things to the XML files
    
    This class provides the code to take a series of input XML files containing
    IFU data (mIFU or LIFU) and updates them according to the PROGTEMP, OBSTEMP
    values, in addition to adding guide star(s) and calibration (WD) targets.

    Parameters
    ----------
    input_loc : str
        Directory of the input XML files
    input_loc : str
        Directory of the output XML files to be written
    """

    def __init__(self,input_loc,output_loc,mifu_ncalibs=2):
        self.version = 0.6
        self.plate_scale = 17.8   ##  "/mm
        print 'Reading in XMLs from %s'%(input_loc)
        self.xmls = glob.glob(input_loc+'/*.xml')
        if len(self.xmls) == 0:
            raise SystemExit('No XML files found in %s'%(input_loc))
        self.output_loc = output_loc
        self.mifu_ncalibs = mifu_ncalibs
        return

    def go(self):
        a = 1
        for xml_file in self.xmls:
            try:
                dom = xml.dom.minidom.parse(xml_file)
            except xml.parsers.expat.ExpatError:
                print("File {0} would not parse".format(xml_file))
                raise SystemExit(0)

            output_file = self.output_loc+os.path.basename(xml_file)
            self.ingest_xml(dom)
            # 1. Update PROGTEMP-related quantities
            self.progtemp()
            # 2. Update OBSTEMP-related quantities
            self.obstemp()
            # 4. Generate calibs where required
            self.calibs(mifu_ncalibs=self.mifu_ncalibs)
            # 3. Generate guidestar(s)
            self.guidestars()

            
            self.write_xml(output_file)
            # tasks:


            
            



    def progtemp(self):
        progtemp_code = self.observation.getAttribute('progtemp')
        chained = False
        if '+' in progtemp_code:
            chained = True
        
        file = open('progtemp.dat','r')
        data = file.readlines()
        progtemp_data = {}
        progtemp_n = {}
        progtemp_o = {}
        progtemp_i = {}
        progtemp_forbidden = []
        progtemp_piforbidden = []
        for line in data:
            if (line[0] == '#') or (line in ['\n','\t']):
                continue
            _key = line.split(':')[0].strip()
            _data = line.split(_key+':')[1].strip()
            try:
                n_index = int(_key[0])
                progtemp_n[n_index] = _data
            except ValueError:
                if _data == '':
                    o_val = int(_key.replace('_',''))
                    progtemp_o[o_val] = {}
                elif _key == 'WEAVE_FORBIDDEN':
                    progtemp_forbidden = _data.split()
                elif _key == 'PI_FORBIDDEN':
                    progtemp_piforbidden = _data.split()
                elif _key[1] != '_':
                    assert int(_key[1]) == o_val
                    progtemp_o[o_val][_key] = _data.replace('{','').replace('}','').split()
                elif (_key[:4] == '____') and _key[-1] != '_':
                    progtemp_i[int(_key[-1])] = _data
                else:
                    import pdb
                    pdb.set_trace()

        # 0. Forbidden-ness checks
        core_progtemp = progtemp_code.split('.')[0]
        verboten_weave = False
        verboten_pi = False
        for forbidden in progtemp_forbidden:
            blank_indices = [a.start() for a in list(re.finditer('_', forbidden))]
            blanked_code = ''
            for i in xrange(len(core_progtemp)):
                if i in blank_indices:
                    blanked_code += '_'
                else:
                    blanked_code += core_progtemp[i]
            if blanked_code == forbidden:
                print 'WARNING: supplied progtemp code (%s) is forbidden for WEAVE science teams as per rule %s'%(progtemp_code,forbidden)
                verboten_weave = True

        for forbidden in progtemp_piforbidden:
            blank_indices = [a.start() for a in list(re.finditer('_', forbidden))]
            blanked_code = ''
            for i in xrange(len(core_progtemp)):
                if i in blank_indices:
                    blanked_code += '_'
                else:
                    blanked_code += core_progtemp[i]
            if blanked_code == forbidden:
                print 'WARNING: supplied progtemp code (%s) is forbidden for instrument as per rule %s'%(progtemp_code,forbidden)
                verboten_pi = True

        if verboten_pi:
            print 'ERROR: Supplied PROGTEMP cannot be used with WEAVE. This OB will be rejected by WASP'
            raise SystemExit(0)
        
        elif verboten_weave:
            print 'WARNING: Supplied PROGTEMP cannot be used by WEAVE science teams. This OB will rejected by WASP unless submitted by an open-time PI project'
                    
        # 1. N progtemp component
        xml_N = progtemp_n[int(progtemp_code[0])].split()
        set_xbinning = []
        for data in xml_N:
            element = self.root
            data_path = data.split(':')
            for _data in data_path:
                if not '=' in _data:
                    element = element.getElementsByTagName(_data)[0]
                else:
                    keyname = _data.split('=')[0]
                    value = _data.split('=')[1]
                    print 'Setting %s to %s in element %s'%(keyname,value,element.nodeName) 
                    element.setAttribute(keyname,value=value)


                    if ('_Arm' in element.nodeName) and (not element.nodeName in set_xbinning):
                        keyname = 'binning_X'
                        value = '1'
                        print 'Setting %s to %s in element %s'%(keyname,value,element.nodeName)
                        element.setAttribute(keyname,value=value)
                        set_xbinning.append(element.nodeName)               


        xml_O = progtemp_o[int(progtemp_code[1])]
        xml_R = int(progtemp_code[2])
        xml_B = int(progtemp_code[3])
        arm_locked = True
        orb_data = []
        if xml_R != xml_B:
            arm_locked = False
        
        if arm_locked:
            xml_ORB = '_'+progtemp_code[1:4]+'_'    
            try:
                xml_ORB = progtemp_o[int(progtemp_code[1])][xml_ORB]
            except KeyError:
                raise SystemExit('The ORB component (%s) of PROGTEMP was not found in valid PROGTEMP values'%(progtemp_code[1:4]))
            orb_data = [xml_ORB]
        else:
            # do red first
            xml_OR = '_'+progtemp_code[1:3]+'__'
            try:
                xml_OR = progtemp_o[int(progtemp_code[1])][xml_OR]
            except KeyError:
                raise SystemExit('The OR component (%s) of PROGTEMP was not found in valid PROGTEMP values'%('_'+progtemp_code[1:3]+'__'))

            # now blue
            xml_OB = '_'+progtemp_code[1]+'_'+progtemp_code[3]+'_'
            try:
                xml_OB = progtemp_o[int(progtemp_code[1])][xml_OB]
            except KeyError:
                raise SystemExit('The OB component (%s) of PROGTEMP was not found in valid PROGTEMP values'%('_'+progtemp_code[1]+'_'+progtemp_code[3]+'_'))
            orb_data = [xml_OR,xml_OB]

        # Now apply this information to the relevant <exposures> elements
        for exp_data in orb_data:
            nexp = int(exp_data[0])
            if arm_locked == False:
                nxp *= 2
            xml_exposures = self.root.getElementsByTagName('exposures')[0].getElementsByTagName('exposure')
            sci_exposures = []
            for exp in xml_exposures:
                if exp.getAttribute('type') == 'science':
                    sci_exposures.append(exp)
            try:
                assert len(sci_exposures) == nexp
            except:
                import pdb
                pdb.set_trace()
    
                raise SystemExit('Expecting %d science exposure elements, retrieve %d'%(nexp,len(sci_exposures)))
            for exp in sci_exposures:
                if exp_data[2].split('=')[-1] == exp.getAttribute('arm'):
                    print 'Setting %s to %s in element %s'%('exp_time',exp_data[1].split('=')[-1],exp.nodeName) 
                    exp.setAttribute('exp_time',value=exp_data[1].split('=')[-1])
        

        # 3. I progtemp component
        xml_I = progtemp_i[int(progtemp_code[4])].split()
        for data in xml_I:
            element = self.root
            data_path = data.split(':')
            for _data in data_path:
                if not '=' in _data:
                    element = element.getElementsByTagName(_data)[0]
                else:
                    keyname = _data.split('=')[0]
                    value = _data.split('=')[1]
                    print 'Setting %s to %s in element %s'%(keyname,value,element.nodeName) 
                    element.setAttribute(keyname,value=value)


        # 4. Set the chained boolean
        obs = self.root.getElementsByTagName('observation')[0]
        print 'Setting %s to %s in element %s'%('chained',str(chained),obs.nodeName) 
        obs.setAttribute('chained',value=str(chained))


    def obstemp(self):
        obstemp_code = self.observation.getAttribute('obstemp')
        file = open('obstemp.dat','r')
        data = file.readlines()
        obstemp_data = {}
        obstemp = 'STAMB'
        
        for line in data:
            if (line[0] == '#') or (line in ['\n','\t']):
                continue
            _key = line.split(':')[0].strip()
            _data = line.split(_key+':')[1].strip()

            if _data == '':
                stamb_key = _key.replace('_','')
                obstemp_data[stamb_key] = {}
                ii = indexOf(obstemp,stamb_key)
            else:
                assert _key[ii] != '_'
                obstemp_data[stamb_key][_key[ii]] = _data.split('#')[0].strip().split()

        for i in xrange(len(obstemp)):
            xml_code = obstemp_code[i]
            xml_data = obstemp_data[obstemp[i]][xml_code]
            for data in xml_data:
                element = self.root
                data_path = data.split(':')
                for _data in data_path:
                    if not '=' in _data:
                        element = element.getElementsByTagName(_data)[0]
                    else:
                        keyname = _data.split('=')[0]
                        value = _data.split('=')[1]
                        print 'Setting %s to %s in element %s'%(keyname,value,element.nodeName) 
                        element.setAttribute(keyname,value=value)

        return

    def calibs(self,mifu_ncalibs=2):
        obs_mode = self.observation.getAttribute('obs_type')
        if obs_mode == 'LIFU':
            # No calibration targets needed!
            return
        pa = 0.0
        field =  self.fields.getElementsByTagName('field')[0]
        # there can be only one field element at this stage...
        ra = float(field.getAttribute('RA_d'))
        dec = float(field.getAttribute('Dec_d'))
        cs = CalibStar(ra,dec,pa,'mIFU',annular=False,plot=False)
        calibs = cs.get_calib(print_xml=False)
        if calibs == None:
            raise SystemExit('Cannot proceed without valid calibration bundle(s)!')

        # need to choose mifu_ncalibs of these
        calib_selection = []
        calib_sector_selection = []
        field =  self.fields.getElementsByTagName('field')[0]
        
        # get angular distribution of calibs:
        calib_countdata = self.angular_dist(calibs,ra,dec)

        while len(calib_selection) != mifu_ncalibs:
            # get non-guide targets (ie bundle centres)
            all_targets = field.getElementsByTagName('target')
            targets = []
            for target in all_targets:
                if target.getAttribute('targuse') in ['C','T']:
                    targets.append(target)
                
            targ_angles = []

            counts = self.angular_dist(targets,ra,dec)
            _sectors = counts.keys()
            _sectors.sort()
            sectors = []
            for s in _sectors:
                if (s in calib_countdata.keys()) and (calib_countdata[s]['count'] > 0):
                    sectors.append(s)
            sectors = numpy.array(sectors)
            sector_counts = []
            sector_counts = numpy.array([counts[s]['count'] for s in sectors])
            if len(calib_sector_selection) > 0:
                sector_counts[3] += 1
            mins = numpy.where(sector_counts == min(sector_counts))[0]
            if len(mins) > 1:
                if len(calib_sector_selection) > 0:
                    # choose something not in already selected sector, ideally far away
                    mdist = []
                    mid = []
                    for m in mins:
                        if not m in calib_sector_selection:
                            mdist.append(sum([abs(m-c) for c in calib_sector_selection]))
                            mid.append(m)
                    if len(mdist) == 0:
                        # you just have to choose something...
                        mins = mins[0]
                    else:
                        mins = mid[indexOf(mdist,max(mdist))]
            else:
                mins = mins[0]
            min_sector = sectors[mins]
            # this is the sector the WD should be selected from
            calib_sector_selection.append(min_sector)
            added = False
            for calib_candidate in calib_countdata[min_sector]['selected']:
                if not calib_candidate in calib_selection:
                    calib_selection.append(calib_candidate)
                    calib_countdata[min_sector]['count'] -= 1
                    added = True
                    break
            if not added:
                raise SystemExit('Could not add viable calibration bundle to field')

        for targ_calib in calib_selection:
            targ0 = field.getElementsByTagName('target')[0]
            field.insertBefore(targ_calib,targ0)

            
    def angular_dist(self,data,ra0,dec0,da=36):
        counts = []
        targ_angles = []
        for target in data:
            # translated position:
            targ_ra = ra0 - float(target.getAttribute('targra'))
            targ_dec = dec0 - float(target.getAttribute('targdec'))
            # determine the azimuthal angles
            targ_ang = numpy.arctan2(targ_dec,targ_ra)
            targ_ang = (targ_ang*180.0) / numpy.pi
            if targ_ang < 0:
                targ_ang += 360.0
            targ_angles.append(targ_ang)


        # counts per sector
        # assume 10 sectors (2 arms per sector)

        count_data = {}
        
        for aa in numpy.arange(0,360,da):
            in_sector = numpy.array([(ang >= aa) and (ang < (aa + da)) for ang in targ_angles])
            count = countOf(in_sector,True)
            selected_targets = numpy.array(data)[in_sector]
            count_data[aa] = {'count':count,'selected':selected_targets}

        return count_data

    
    def guidestars(self):
        obs_mode = self.observation.getAttribute('obs_type')
        pa = 0.0
        max_guides = {}
        max_guides['LIFU'] = None
        max_guides['mIFU'] = int(self.configure.getAttribute('max_guide'))
        print 'WARNING: PA set to 0.0 for the moment'
        field =  self.fields.getElementsByTagName('field')[0]
        # there can be only one field element at this stage...
        ra = float(field.getAttribute('RA_d'))
        dec = float(field.getAttribute('Dec_d'))
        gs = GuideStar(ra,dec,pa,obs_mode,max_guides=max_guides[obs_mode])
        guides = gs.get_guide(as_xml=True)

        if guides == None:
            raise SystemExit('Cannot proceed without valid guidestar')
        if type(guides) != type([]):
            guides = [guides]
            
        for guide in guides[:min(len(guides),int(self.configure.getAttribute('max_guide')))]:
            # get the first <target> in the field, add this one before it
            targ0 = field.getElementsByTagName('target')[0]
            field.insertBefore(guide,targ0)

            
        # import numpy
        # index = numpy.arange(len(guides))
        # numpy.random.shuffle(index)
        # for i in index[:8]:
        #     g = guides[i]
        #     if 1:
        #         print gs.guides_filter['ANGLE'][i]
        #     print g.toxml()

    
            
    def ingest_xml(self,dom):
        self.dom = dom
        self.root = dom.childNodes[0]
        self.exposures = self.root.getElementsByTagName('exposures')[0]
        self.observation = self.root.getElementsByTagName('observation')[0]
        self.configure = dom.getElementsByTagName('configure')[0]
        self.surveys = dom.getElementsByTagName('surveys')[0]
        self.obsconstraints = dom.getElementsByTagName('obsconstraints')[0]
        self.dithering = dom.getElementsByTagName('dithering')[0]
        self.fields = dom.getElementsByTagName('fields')[0]
        self.base_target = self.fields.getElementsByTagName('target')[0]
        self.offset = self.observation.getElementsByTagName('offsets')[0]
        self.targets_base = self.fields.getElementsByTagName('target')
            


    def write_xml(self,filename):
        newxml = self.remove_empty_lines(self.dom.toprettyxml())
        finalxml = self.remove_xml_declaration(newxml)
        print 'Writing to %s'%(filename)
        with open(filename, 'w') as f:
            f.write(finalxml)    

    def remove_xml_declaration(self,xml_text):
        doc = xml.dom.minidom.parseString(xml_text)
        root = doc.documentElement
        xml_text_without_declaration = root.toxml(doc.encoding)
        return xml_text_without_declaration

    def remove_empty_lines(self,xml_text):
        reparsed = xml.dom.minidom.parseString(xml_text)
        return '\n'.join([line for line in reparsed.toprettyxml(indent=' ').split('\n') if line.strip()])

        

if __name__ == '__main__':

    import os


    
    input_xml_loc = './input/'
    
    try:
        assert os.path.isdir(input_xml_loc)
    except:
        print 'Please supply valid XML location'
        print 'Usage: ./stage3.py location'
        raise SystemExit(0)

    output_dir = '../stage4/input_tmp/'

    if 1:
        cmd = 'rm -Rf %s/*.xml'%(output_dir)
        os.system(cmd)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    
    stage3_ifu = stage3(input_xml_loc,output_dir)   ## add specifiers for mIFU grouping / overloading fields /etc behaviour here
    stage3_ifu.go()
    
    print('IFU XMLs written to: {0}'.format(output_dir))
