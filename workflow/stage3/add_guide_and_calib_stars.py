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
import logging
import os
import xml.dom.minidom
from operator import indexOf, countOf

import numpy as np

from workflow.utils.classes import GuideStars
from workflow.utils.classes import CalibStars
from workflow.utils.classes import OBXML


class stage3(OBXML):

    
    def _get_obsmode(self):

        obsmode = self.observation.getAttribute('obs_type')

        return obsmode

    
    def _get_central_ra_dec(self, obsmode):

        # Get the first field

        first_field = self.fields.getElementsByTagName('field')[0]

        # For non-LIFU observations, get the centre from the field element
        if obsmode != 'LIFU':

            central_ra = float(first_field.getAttribute('RA_d'))
            central_dec = float(first_field.getAttribute('Dec_d'))

        else:

            central_ra = None
            central_dec = None

            for target in first_field.getElementsByTagName('target'):
                
                targuse = target.getAttribute('targuse')

                if targuse == 'T':
                    central_ra = float(target.getAttribute('targra'))
                    central_dec = float(target.getAttribute('targdec'))

        return central_ra, central_dec

    
    def _get_pa(self):

        str_pa = self.observation.getAttribute('pa')

        if str_pa != '%%%':
            pa = float(str_pa)
        else:
            pa = np.nan

        return pa

    
    def _get_max_guide(self):

        max_guide = int(self.configure.getAttribute('max_guide'))

        return max_guide


    def _get_guide_stars(self):

        obsmode = self._get_obsmode()

        central_ra, central_dec = self._get_central_ra_dec(obsmode)
        pa = self._get_pa()
        max_guide = self._get_max_guide()

        guide_stars = GuideStars(central_ra, central_dec, pa, obsmode,
                                 max_guide=max_guide)

        actual_pa, full_guides_table = guide_stars.get_table()

        guides_table = full_guides_table[0:max_guide]

        return actual_pa, guides_table


    def _set_pa(self, pa):

        self._set_attribs(self.observation, {'pa': pa})

        
    def _add_target(self, field, target_attrib_dict, photometry_attrib_dict):

        for target in field.getElementsByTagName('target'):
            targuse = target.getAttribute('targuse')
            if targuse == 'T':
                first_science_target = target
                break

        new_target = first_science_target.cloneNode(True)

        self._set_attribs(new_target, target_attrib_dict)

        new_target_photometry = new_target.getElementsByTagName('photometry')[0]

        self._set_attribs(new_target_photometry, photometry_attrib_dict)

        field.insertBefore(new_target, first_science_target)


    def _set_guide_stars(self, actual_pa, guides_table):

        # Update the PA value if needed

        pa = self._get_pa()

        if actual_pa != pa:
            if not np.isnan(pa):
                logging.info('Requested value for PA has NOT been adopted')
            self._set_pa(actual_pa)

        # Add the guide stars to the XML file

        if len(guides_table) == 0:
            logging.error('There is not guide stars available')
            raise SystemExit(2)

        col_to_attrib_target_dict = {
            'CNAME': 'cname',
            'TARGCAT': 'targcat',
            'TARGCLASS': 'targclass',
            'GAIA_DEC': 'targdec',
            'GAIA_EPOCH': 'targepoch',
            'TARGID': 'targid',
            'TARGNAME': 'targname',
            'GAIA_PARAL': 'targparal',
            'GAIA_PMDEC': 'targpmdec',
            'GAIA_PMRA': 'targpmra',
            'TARGPRIO': 'targprio',
            'TARGPROG': 'targprog',
            'GAIA_RA': 'targra',
            'TARGSRVY': 'targsrvy',
            'TARGUSE': 'targuse'
        }
        
        col_to_attrib_photometry_dict = {
            'GAIA_EMAG_BP': 'emag_bp',
            'EMAG_G': 'emag_g',
            'GAIA_EMAG_GG': 'emag_gg',
            'EMAG_I': 'emag_i',
            'EMAG_R': 'emag_r',
            'GAIA_EMAG_RP': 'emag_rp',
            'GAIA_MAG_BP': 'mag_bp',
            'MAG_G': 'mag_g',
            'GAIA_MAG_GG': 'mag_gg',
            'MAG_I': 'mag_i',
            'MAG_R': 'mag_r',
            'GAIA_MAG_RP': 'mag_rp',
        }

        for guide in guides_table:
            
            target_attrib_dict = {
                col_to_attrib_target_dict[col]: guide[col]
                for col in col_to_attrib_target_dict.keys()
            }

            photometry_attrib_dict = {
                col_to_attrib_photometry_dict[col]: guide[col]
                for col in col_to_attrib_photometry_dict.keys()
            }

            for field in self.fields.getElementsByTagName('field'):
                self._add_target(field, target_attrib_dict,
                                 photometry_attrib_dict)


    def _add_guide_stars(self):

        actual_pa, guides_table = self._get_guide_stars()

        self._set_guide_stars(actual_pa, guides_table)

            
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

        target_template = self.fields.getElementsByTagName('target')[2]
    
        xmls = []
    
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

                        
    def _add_calib_stars(self,mifu_num_calibs=2):
        obsmode = self.observation.getAttribute('obs_type')
        if obsmode == 'LIFU':
            # No calibration targets needed!
            return
        field =  self.fields.getElementsByTagName('field')[0]
        all_targets = field.getElementsByTagName('target')
        targ0 = field.getElementsByTagName('target')[0]

        ra = float(field.getAttribute('RA_d'))
        dec = float(field.getAttribute('Dec_d'))

        pa = 0.0
        # there can be only one field element at this stage...
        cs = CalibStars(ra,dec,pa,'mIFU',annular=False,plot=False)
        calibs_table = cs.get_calib()

        if calibs_table == None:
            raise SystemExit('Cannot proceed without valid calibration bundle(s)!')

        calibs = self._calibs_to_xml(calibs_table)

        # need to choose mifu_num_calibs of these
        calib_selection = []
        calib_sector_selection = []
        
        # get angular distribution of calibs:
        calib_countdata = self._angular_dist(calibs,ra,dec,input_type='xml')

        while len(calib_selection) != mifu_num_calibs:
            # get non-guide targets (ie bundle centres)
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
            field.insertBefore(targ_calib,targ0)

        
    def add_guide_and_calib_stars(self,mifu_num_calibs=2):

        self._add_guide_stars()
        self._add_calib_stars(mifu_num_calibs=mifu_num_calibs)


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
        stage3_ifu.add_guide_and_calib_stars()
        stage3_ifu.write_xml(output_file)
    
    print('IFU XMLs written to: {0}'.format(output_dir))

