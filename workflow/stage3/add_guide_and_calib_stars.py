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


    def _add_table_as_targets(self, table):

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

        for row in table:
            
            target_attrib_dict = {
                col_to_attrib_target_dict[col]: row[col]
                for col in col_to_attrib_target_dict.keys()
            }

            photometry_attrib_dict = {
                col_to_attrib_photometry_dict[col]: row[col]
                for col in col_to_attrib_photometry_dict.keys()
            }

            for field in self.fields.getElementsByTagName('field'):
                self._add_target(field, target_attrib_dict,
                                 photometry_attrib_dict)


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

        self._add_table_as_targets(guides_table)

                
    def _add_guide_stars(self):

        actual_pa, guides_table = self._get_guide_stars()

        self._set_guide_stars(actual_pa, guides_table)

                        
    def _get_calib_stars(self, mifu_num_calibs=2):

        obsmode = self._get_obsmode()

        central_ra, central_dec = self._get_central_ra_dec(obsmode)
        pa = self._get_pa()

        calib_stars = CalibStars(central_ra, central_dec, pa, obsmode)

        full_calibs_table = calib_stars.get_table()

        calibs_table = full_calibs_table[0:mifu_num_calibs]

        return calibs_table


    def _set_calib_stars(self, calibs_table):

        # Add the guide stars to the XML file

        if len(calibs_table) == 0:
            logging.error('There is not guide stars available')
            raise SystemExit(2)

        self._add_table_as_targets(calibs_table)


    def _add_calib_stars(self, mifu_num_calibs=2):

        obsmode = self._get_obsmode()

        if obsmode == 'LIFU':
            # No calibration stars needed!
            return

        calibs_table = self._get_calib_stars(mifu_num_calibs=mifu_num_calibs)

        self._set_calib_stars(calibs_table)

        
    def add_guide_and_calib_stars(self, mifu_num_calibs=2):

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

