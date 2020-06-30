#!/usr/bin/env python3

#
# Copyright (C) 2020 Cambridge Astronomical Survey Unit
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


import glob
import os
import xml.dom.minidom
from operator import indexOf, countOf

import numpy as np

from workflow.utils.classes import GuideStars
from workflow.utils.classes import CalibStars
from workflow.utils.classes import OBXML


class stage3(OBXML):
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

            
    def _angular_dist(self,data,ra0,dec0,da=36,input_type='table'):
        counts = []
        targ_angles = []
        for target in data:
            # translated position:
            if input_type == 'table':
                targ_ra = ra0 - target['GAIA_RA']
                targ_dec = dec0 - target['GAIA_DEC']
            else:
                targ_ra = ra0 - float(target.getAttribute('targra'))
                targ_dec = dec0 - float(target.getAttribute('targdec'))
            # determine the azimuthal angles
            targ_ang = np.arctan2(targ_dec,targ_ra)
            targ_ang = (targ_ang*180.0) / np.pi
            if targ_ang < 0:
                targ_ang += 360.0
            targ_angles.append(targ_ang)


        # counts per sector
        # assume 10 sectors (2 arms per sector)

        count_data = {}
        
        for aa in np.arange(0,360,da):
            in_sector = np.array([(ang >= aa) and (ang < (aa + da)) for ang in targ_angles])
            count = countOf(in_sector,True)
            selected_targets = np.array(data)[in_sector]
            count_data[aa] = {'count':count,'selected':selected_targets}

        return count_data

    
    def _calibs_to_xml(self,calibs):
    
        xmls = []

        target_template = self.fields.getElementsByTagName('target')[2]
    
        for calib in calibs:

            xml_target = target_template.cloneNode(True)


            ## NB: this will FAIL when we move from OpR3b catalogues, due to new column names for errors etc

            xml_target.setAttribute('cname',str(calib['CNAME']))
            xml_target.setAttribute('targid',str(calib['TARGID']))
            xml_target.setAttribute('targra',str(calib['GAIA_RA']))
            xml_target.setAttribute('targdec',str(calib['GAIA_DEC']))
            xml_target.setAttribute('targpmra',str(calib['GAIA_PMRA']))
            xml_target.setAttribute('targpmdec',str(calib['GAIA_PMDEC']))
            xml_target.setAttribute('targprio',str(calib['TARGPRIO']))
            xml_target.setAttribute('targuse',str(calib['TARGUSE']))
            xml_target.setAttribute('targsrvy',str(calib['TARGSRVY']))
            xml_target.setAttribute('targname',str(calib['TARGNAME']))
            xml_target.setAttribute('targprog',str(calib['TARGPROG']))
            xml_target.setAttribute('targclass',str(calib['TARGCLASS']))
            xml_target.setAttribute('targcat',str(calib['TARGCAT']))
            xml_target.setAttribute('targepoch',str(calib['GAIA_EPOCH']))
            xml_target.setAttribute('targparal',str(calib['GAIA_PARAL']))
            xml_target.setAttribute('targprio',"10")

            xml_photom = xml_target.getElementsByTagName('photometry')[0]
            bands = ['g','r','i']
            for b in bands:
                xml_photom.setAttribute('mag_%s'%(b),"")
                xml_photom.setAttribute('emag_%s'%(b),"")
            xml_photom.setAttribute('mag_gg',str(calib['GAIA_MAG_GG']))
            xml_photom.setAttribute('emag_gg',str(calib['GAIA_EMAG_GG']))
            xml_photom.setAttribute('mag_bp',str(calib['GAIA_MAG_BP']))
            xml_photom.setAttribute('emag_bp',str(calib['GAIA_EMAG_BP']))
            xml_photom.setAttribute('mag_rp',str(calib['GAIA_MAG_RP']))
            xml_photom.setAttribute('emag_rp',str(calib['GAIA_EMAG_RP']))
            
            xmls.append(xml_target)

        return xmls

                        
    def _calibs(self,mifu_ncalibs=2):
        obs_mode = self.observation.getAttribute('obs_type')
        if obs_mode == 'LIFU':
            # No calibration targets needed!
            return
        pa = 0.0
        field =  self.fields.getElementsByTagName('field')[0]
        # there can be only one field element at this stage...
        ra = float(field.getAttribute('RA_d'))
        dec = float(field.getAttribute('Dec_d'))
        cs = CalibStars(ra,dec,pa,'mIFU',annular=False,plot=False)
        calibs_table = cs.get_calib()

        if calibs_table == None:
            raise SystemExit('Cannot proceed without valid calibration bundle(s)!')

        calibs = self._calibs_to_xml(calibs_table)

        # need to choose mifu_ncalibs of these
        calib_selection = []
        calib_sector_selection = []
        field =  self.fields.getElementsByTagName('field')[0]
        
        # get angular distribution of calibs:
        calib_countdata = self._angular_dist(calibs,ra,dec,input_type='xml')

        while len(calib_selection) != mifu_ncalibs:
            # get non-guide targets (ie bundle centres)
            all_targets = field.getElementsByTagName('target')
            targets = []
            for target in all_targets:
                if target.getAttribute('targuse') in ['C','T']:
                    targets.append(target)
                
            targ_angles = []

            counts = self._angular_dist(targets,ra,dec,input_type='xml')
            _sectors = list(counts.keys())
            _sectors.sort()
            sectors = []
            for s in _sectors:
                if (s in calib_countdata.keys()) and (calib_countdata[s]['count'] > 0):
                    sectors.append(s)
            sectors = np.array(sectors)
            sector_counts = []
            sector_counts = np.array([counts[s]['count'] for s in sectors])
            if len(calib_sector_selection) > 0:
                sector_counts[3] += 1
            mins = np.where(sector_counts == min(sector_counts))[0]
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


    def _guides_to_xml(self,guides):
    
        xmls = []

        target_template = self.fields.getElementsByTagName('target')[0]

        for guide in guides:

            xml_target = target_template.cloneNode(True)

            xml_target.setAttribute('cname',str(guide['CNAME']))
            xml_target.setAttribute('targid',str(guide['TARGID']))
            xml_target.setAttribute('targra',str(guide['GAIA_RA']))
            xml_target.setAttribute('targdec',str(guide['GAIA_DEC']))
            xml_target.setAttribute('targpmra',str(guide['GAIA_PMRA']))
            xml_target.setAttribute('targpmdec',str(guide['GAIA_PMDEC']))
            xml_target.setAttribute('targprio',str(guide['TARGPRIO']))
            xml_target.setAttribute('targuse',str(guide['TARGUSE']))
            xml_target.setAttribute('targsrvy',str(guide['TARGSRVY']))
            xml_target.setAttribute('targname',str(guide['TARGNAME']))
            xml_target.setAttribute('targprog',str(guide['TARGPROG']))
            xml_target.setAttribute('targclass',str(guide['TARGCLASS']))
            xml_target.setAttribute('targcat',str(guide['TARGCAT']))
            xml_target.setAttribute('targepoch',str(guide['GAIA_EPOCH']))
            xml_target.setAttribute('targparal',str(guide['GAIA_PARAL']))
            xml_target.setAttribute('targprio',"10")

            xml_photom = xml_target.getElementsByTagName('photometry')[0]
            bands = ['g','r','i']
            for b in bands:
                xml_photom.setAttribute('mag_%s'%(b),"")
                xml_photom.setAttribute('emag_%s'%(b),"")
            xml_photom.setAttribute('mag_gg',str(guide['GAIA_MAG_GG']))
            xml_photom.setAttribute('emag_gg',str(guide['GAIA_EMAG_GG']))
            xml_photom.setAttribute('mag_bp',str(guide['GAIA_MAG_BP']))
            xml_photom.setAttribute('emag_bp',str(guide['GAIA_EMAG_BP']))
            xml_photom.setAttribute('mag_rp',str(guide['GAIA_MAG_RP']))
            xml_photom.setAttribute('emag_rp',str(guide['GAIA_EMAG_RP']))
    
            xmls.append(xml_target)

        return xmls

    
    def _guidestars(self):

        obs_mode = self.observation.getAttribute('obs_type')
        max_guides = {}
        max_guides['LIFU'] = None
        max_guides['mIFU'] = int(self.configure.getAttribute('max_guide'))

        if obs_mode == 'mIFU':
            pa = 0.0
            
            if np.isnan(float(self.observation.getAttribute('pa'))):
                print('assertion error to be revisited')
            else:
                assert float(self.observation.getAttribute('pa')) == pa
            
        else:
            pa = 0.0
            if self.observation.getAttribute('pa') != '%%%':
                pa = float(self.observation.getAttribute('pa'))

        if 1:
            # testing override
            pa = 0.0
            print('WARNING: PA set to 0.0 for the moment')


        field =  self.fields.getElementsByTagName('field')[0]
        # there can be only one field element at this stage...
        ra = float(field.getAttribute('RA_d'))
        dec = float(field.getAttribute('Dec_d'))
        gs = GuideStars(ra,dec,pa,obs_mode,max_guides=max_guides[obs_mode])

        pa_actual = pa
        print('WARNING: guidestar search will not adopt new PA - implement this!')
        # guides,pa_actual = gs.get_guide(as_xml=True)
        guides_table = gs.get_guide(as_xml=True)

        guides = self._guides_to_xml(guides_table)
        if obs_mode == 'LIFU':   
            self.observation.setAttribute('pa',value=str(pa_actual))

        
        if guides == None:
            raise SystemExit('Cannot proceed without valid guidestar')
            
        for guide in guides[:min(len(guides),int(self.configure.getAttribute('max_guide')))]:
            # get the first <target> in the field, add this one before it
            targ0 = field.getElementsByTagName('target')[0]
            field.insertBefore(guide,targ0)

        # index = np.arange(len(guides))
        # np.random.shuffle(index)
        # for i in index[:8]:
        #     g = guides[i]
        #     if 1:
        #         print(gs.guides_filter['ANGLE'][i])
        #     print(g.toxml())

        
    def go(self,mifu_ncalibs=2):

        # 1. Generate calibs where required
        self._calibs(mifu_ncalibs=mifu_ncalibs)
        # 2. Generate guidestar(s)
        self._guidestars()


if __name__ == '__main__':

    
    input_xml_loc = './input/'
    
    try:
        assert os.path.isdir(input_xml_loc)
    except:
        print('Please supply valid XML location')
        print('Usage: ./stage3.py location')
        raise SystemExit(0)

    output_dir = '../stage3/output/'

    if 1:
        cmd = 'rm -Rf %s/*.xml'%(output_dir)
        os.system(cmd)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print('Reading in XMLs from %s'%(input_xml_loc))

    filename_list = glob.glob(input_xml_loc+'/*.xml')
    filename_list.sort()

    if len(filename_list) == 0:
        raise SystemExit('No XML files found in %s'%(input_loc))

    for filename in filename_list:
        output_file = output_dir+os.path.basename(filename)
        
        stage3_ifu = stage3(filename)
        stage3_ifu.go()
        stage3_ifu.write_xml(output_file)
    
    print('IFU XMLs written to: {0}'.format(output_dir))

