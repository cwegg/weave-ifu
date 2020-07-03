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


import logging
import xml.dom.minidom as _minidom

import numpy as _np

from ._calibstars import CalibStars as _CalibStars
from ._guidestars import GuideStars as _GuideStars


class OBXML:

    
    def __init__(self, filename):

        # Save the input filename
        
        self.filename = filename

        # Parse the XML template

        try:

            self.dom = _minidom.parse(self.filename)

        except:

            logging.error('File {} would not parse'.format(self.filename))

            raise SystemExit(1)

        # Ingest the XML template

        self._ingest_xml()

        # Set the order of the first science exposure as unknown

        self._set_first_science_order()


    def _ingest_xml(self):

        self.root = self.dom.childNodes[0]

        self.spectrograph = self.dom.getElementsByTagName('spectrograph')[0]
        self.exposures = self.dom.getElementsByTagName('exposures')[0]
        self.observation = self.dom.getElementsByTagName('observation')[0]
        self.configure = self.dom.getElementsByTagName('configure')[0]
        self.obsconstraints = self.dom.getElementsByTagName('obsconstraints')[0]
        self.dithering = self.dom.getElementsByTagName('dithering')[0]
        self.surveys = self.dom.getElementsByTagName('surveys')[0]
        self.fields = self.dom.getElementsByTagName('fields')[0]


    def _set_first_science_order(self):

        last_pre_science_order = 0

        for exposure in self.exposures.getElementsByTagName('exposure'):

            if exposure.getAttribute('type') != 'science':

                order = exposure.getAttribute('order')

                if (order != '%%%'):
                    last_pre_science_order = int(order)
                else:
                    break

            else:

                break

        self.first_science_order = last_pre_science_order + 1

        
    def get_datamver(self):

        datamver = self.root.getAttribute('datamver')

        return datamver


    def _remove_elem_list_by_name(self, elem_name_list):

        for elem_name in elem_name_list:
            for node in self.dom.getElementsByTagName(elem_name):
                parent = node.parentNode
                parent.removeChild(node)


    def remove_configure_outputs(self):

        # Remove some attributes from the observation element

        configure_attrib_list = ['configure_version', 'seed']

        for attrib in configure_attrib_list:
            if attrib in self.configure.attributes.keys():
                self.configure.removeAttribute(attrib)

        # Remove some attributes from the target elements

        target_attrib_list = [
            'automatic', 'configid', 'fibreid', 'ifu_pa', 'ifu_spaxel',
            'targx', 'targy'
        ]

        for node in self.dom.getElementsByTagName('target'):
            for attrib in target_attrib_list:
                if attrib in node.attributes.keys():
                    node.removeAttribute(attrib)

        # Remove some elements

        elem_list = [
            'optical_axis', 'distortion_coefficients',
            'telescope', 'focal_plane_map', 'hour_angle_limits', 'offsets'
        ]

        self._remove_elem_list_by_name(elem_list)

                
    def remove_non_used_elements(self, elem_list=['simulation',
                                                  'avoidance_list', 'group']):

        self._remove_elem_list_by_name(elem_list)


    def _get_num_decimals(self, str_value):

        if '.' in str_value:
            num_decimals = len(str_value) - str_value.find('.') - 1
        else:
            num_decimals = 0

        return num_decimals


    def _set_attribs(self, elem, attrib_dict):

        for key in attrib_dict.keys():

            if key in elem.attributes.keys():

                value = attrib_dict[key]

                try:
                    isnan_flag = _np.isnan(value)
                except:
                    isnan_flag = False

                if isnan_flag == False:
                    str_value = str(value)

                    # For some attributes, configure wants to recieve a
                    # precision better than 0.1 arcsec: We will get a string
                    # with trailing zeros if the default decimals are not enough
                    
                    if key in ['targra', 'targdec', 'RA_d', 'Dec_d']:
                        if self._get_num_decimals(str_value) < 5:
                            str_value = '{:.5f}'.format(value)
                else:
                    str_value = ''

                elem.setAttribute(key, value=str_value)

            else:

                raise KeyError


    def set_root_attrib(self, attrib_dict):

        self._set_attribs(self.root, attrib_dict)

        
    def set_spectrograph(self, binning_y=None, resolution=None,
                         red_vph=None, blue_vph=None):

        vph_dict = {'red': red_vph, 'blue': blue_vph}

        for colour in ['red', 'blue']:

            elem_name = '{}_Arm'.format(colour)

            elem = self.spectrograph.getElementsByTagName(elem_name)[0]

            vph = vph_dict[colour]

            attrib_dict = {}

            if binning_y is not None:
                attrib_dict['binning_Y'] = binning_y

            if resolution is not None:
                attrib_dict['resolution'] = resolution

            if vph_dict[colour] is not None:
                attrib_dict['VPH'] = vph_dict[colour]

            self._set_attribs(elem, attrib_dict)

                
    def set_exposures(self, num_science_exposures, science_exptime,
                      clean_comment=True):

        # Get the science exposure element from the template
        # (it must be the first one without assigned order)

        for exposure in self.exposures.getElementsByTagName('exposure'):
            if exposure.getAttribute('order') == '%%%':
                assert exposure.getAttribute('type') == 'science'
                science_exposure_template = exposure
                break

        # Add the requested amount of science exposures before the science
        # exposure from the template

        order = self.first_science_order - 1
        
        for i in range(num_science_exposures):

            order += 1

            science_exposure = science_exposure_template.cloneNode(True)

            attrib_dict = {
                'arm': 'both',
                'exp_time': science_exptime,
                'order': order
            }

            self._set_attribs(science_exposure, attrib_dict)
            
            self.exposures.insertBefore(science_exposure,
                                        science_exposure_template)

        # Remove the science exposure element from the template

        self.exposures.removeChild(science_exposure_template)

        # Assign the order value to the calibration exposures after the
        # science exposures

        previous_arm = 'both'
        current_arm = 'both'
        arm_flag = False

        for exposure in self.exposures.getElementsByTagName('exposure'):

            if exposure.getAttribute('order') == '%%%':

                # Guess whether the order should be increased or not

                previous_arm = current_arm
                current_arm = exposure.getAttribute('arm')

                if ((previous_arm == 'both') or (current_arm == 'both') or
                    (current_arm == previous_arm) or (arm_flag == True)):
                    order += 1
                    arm_flag = False
                else:
                    arm_flag = True

                # Set the order

                attrib_dict = {'order': order}

                self._set_attribs(exposure, attrib_dict)

        # Identify the last comment node in the exposures elements

        if clean_comment == True:

            comment_string = '[[ORDER undefined until exposure code known]]'

            for node in self.exposures.childNodes:
                if node.nodeType is _minidom.Node.COMMENT_NODE:
                    if comment_string in node.data:
                        self.exposures.removeChild(node)

    
    def set_observation(self, attrib_dict):

        self._set_attribs(self.observation, attrib_dict)


    def set_configure(self, attrib_dict):

        self._set_attribs(self.configure, attrib_dict)

        
    def set_obsconstraints(self, attrib_dict):

        self._set_attribs(self.obsconstraints, attrib_dict)

    
    def set_dithering(self, attrib_dict):

        self._set_attribs(self.dithering, attrib_dict)

        
    def set_surveys(self, targsrvy_list, max_fibres, priority='1.0'):

        # Get the survey element from the template (there should be only one)

        survey_template = self.surveys.getElementsByTagName('survey')[0]

        # Add the requested amount of survey elements before the template

        for targsrvy in targsrvy_list:

            survey = survey_template.cloneNode(True)

            attrib_dict = {
                'name': targsrvy,
                'priority': priority,
                'max_fibres': max_fibres
            }

            self._set_attribs(survey, attrib_dict)
            
            self.surveys.insertBefore(survey, survey_template)

        # Remove the survey element from the template

        self.surveys.removeChild(survey_template)


    def _get_mid_value(self, list):

        min_value = _np.min(list)
        max_value = _np.max(list)

        mid_value = (min_value + max_value) / 2

        return mid_value


    def _get_mid_ra(self, ra_list, dec_value, max_dist=3):

        min_ra_value = _np.min(ra_list)
        max_ra_value = _np.max(ra_list)

        diff_ra = (max_ra_value - min_ra_value) * _np.cos(dec_value)

        # If we are not in the edge, get the trivial solution

        if diff_ra < max_dist:

            mid_value = self._get_mid_value([min_ra_value, max_ra_value])

        # If it is in edge of the coordinate system, it is a bit tricky

        else:

            left_edge_ra_value = _np.min([ra for ra in ra_list if ra >= 180])
            right_edge_ra_value = _np.max([ra for ra in ra_list if ra < 180])

            edge_diff_ra = (((right_edge_ra_value + 360) - left_edge_ra_value) *
                            _np.cos(dec_value))

            if edge_diff_ra < max_dist:

                mid_value = self._get_mid_value([left_edge_ra_value,
                                                 (right_edge_ra_value + 360)])

                if mid_value >= 360:
                    midvalue -= 360

            else:

                raise ValueError('unexpected range of RA values')

        return mid_value

            
    def set_fields(self, obsmode, entry_group, targcat=None,
                   num_science_exposures=None):

        # Check the input parameters

        assert obsmode in ['LIFU', 'mIFU']

        if obsmode == 'LIFU':
            assert num_science_exposures % len(entry_group) == 0

        # A dictionary for mapping column names of the entries to attributes

        col_to_attrib_dict = {
            'TARGSRVY': 'targsrvy',
            'TARGPROG': 'targprog',
            'TARGID': 'targid',
            'TARGNAME': 'targname',
            'TARGPRIO': 'targprio',
            'GAIA_RA': 'targra',
            'GAIA_DEC': 'targdec',
            'GAIA_EPOCH': 'targepoch',
            'GAIA_PMRA': 'targpmra',
            'GAIA_PMDEC': 'targpmdec',
            'GAIA_PARAL': 'targparal'
        }

        # Get a empty field element to use it as template

        field_list = self.fields.getElementsByTagName('field')
        
        field_template = field_list[0].cloneNode(True)

        for node in field_template.getElementsByTagName('target'):
            field_template.removeChild(node)

        for node in field_template.childNodes:
            if node.nodeType is _minidom.Node.COMMENT_NODE:
                field_template.removeChild(node)

        # Get a target element and save it as template

        target_list = self.fields.getElementsByTagName('target')

        target_template = None

        for target in target_list:
            if target.getAttribute('targuse') == 'T':
                target_template = target.cloneNode(True)
                break

        assert target_template is not None

        # Clean the fields element

        for node in self.fields.getElementsByTagName('field'):
            self.fields.removeChild(node)

        # For LIFU:
        # - If it is a non-custom dither,
        #   one field with one target will be added
        # - If it is a custom dither,
        #   one field with one target will be added per science exposure
        #   (it is worth noting that in case of having 3 dither positions and
        #    6 exposures, we will repeat the dither pattern twice)

        if obsmode == 'LIFU':

            assert num_science_exposures % len(entry_group) == 0
            
            if len(entry_group) != 1:
                num_fields = num_science_exposures
            else:
                num_fields = 1

            for i in range(num_fields):

                entry = entry_group[i % len(entry_group)]

                # Create a new field for the entry

                order = self.first_science_order + i

                # field_ra = entry['GAIA_RA']
                # field_dec = entry['GAIA_DEC']

                field = field_template.cloneNode(True)

                field_attrib_dict = {
                    # 'Dec_d': field_dec,
                    # 'RA_d': field_ra,
                    'order': order
                    }

                self._set_attribs(field, field_attrib_dict)

                # Create a target for the entry

                target = target_template.cloneNode(True)

                target_attrib_dict = {}

                target_attrib_dict['targuse'] = 'T'
                target_attrib_dict['targcat'] = targcat

                for col in col_to_attrib_dict.keys():
                    target_attrib_dict[col_to_attrib_dict[col]] = entry[col]

                self._set_attribs(target, target_attrib_dict)

                # Add target to field

                field.appendChild(target)

                # Add field to fields

                self.fields.appendChild(field)

        # For mIFU: One field with sereral targets will be added

        elif obsmode == 'mIFU':

            # Create a new field

            order = self.first_science_order

            ra_list = [entry['GAIA_RA'] for entry in entry_group]
            dec_list = [entry['GAIA_DEC'] for entry in entry_group]

            field_dec = self._get_mid_value(dec_list)
            field_ra = self._get_mid_ra(ra_list, field_dec)

            field = field_template.cloneNode(True)

            field_attrib_dict = {
                'Dec_d': field_dec,
                'RA_d': field_ra,
                'order': order
                }

            self._set_attribs(field, field_attrib_dict)

            # Add a target per entry to the field

            for entry in entry_group:

                target = target_template.cloneNode(True)

                target_attrib_dict = {}

                target_attrib_dict['targuse'] = 'T'
                target_attrib_dict['targcat'] = targcat

                for col in col_to_attrib_dict.keys():
                    target_attrib_dict[col_to_attrib_dict[col]] = entry[col]

                self._set_attribs(target, target_attrib_dict)

                field.appendChild(target)

            # Add the field to fields

            self.fields.appendChild(field)


    def write_xml(self, filename):

        # Get a pretty indented XML with its proper XML declaration
    
        pretty_xml = self.dom.toprettyxml()

        parsed_xml = _minidom.parseString(pretty_xml)

        pretty_xml2 = parsed_xml.toprettyxml(indent='  ',
                                             encoding='utf-8').decode('utf-8')

        output_xml = '\n'.join([line for line in pretty_xml2.split('\n')
                                if line.strip()]) + '\n'
        
        # Write the XML text to a file

        logging.info('\tWriting to {}'.format(filename))

        with open(filename, 'w') as f:
            f.write(output_xml)

    
    def _get_obsmode(self):

        obsmode = self.observation.getAttribute('obs_type')

        return obsmode

    
    def _get_central_ra_dec(self, obsmode):

        # Get the first field

        first_field = self.fields.getElementsByTagName('field')[0]

        # For non-LIFU observations, get the centre from the field element,
        # for LIFU observations, get the centre from the central fibre
        # (which should be the only one with TARGUSE=T when this method is
        # called, or the first one)

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

                    break

        return central_ra, central_dec

    
    def _get_pa(self):

        str_pa = self.observation.getAttribute('pa')

        if str_pa != '%%%':
            pa = float(str_pa)
        else:
            pa = _np.nan

        return pa

    
    def _get_max_guide(self):

        max_guide = int(self.configure.getAttribute('max_guide'))

        return max_guide


    def _get_guide_stars(self, all_guide_stars=False):

        obsmode = self._get_obsmode()

        central_ra, central_dec = self._get_central_ra_dec(obsmode)
        pa = self._get_pa()
        max_guide = self._get_max_guide()

        guide_stars = _GuideStars(central_ra, central_dec, pa, obsmode,
                                  max_guide=max_guide)

        actual_pa, full_guides_table = guide_stars.get_table()

        if all_guide_stars is False:
            guides_table = full_guides_table[0:max_guide]
        else:
            guides_table = full_guides_table

        return actual_pa, guides_table


    def _set_pa(self, pa):

        self._set_attribs(self.observation, {'pa': pa})

        
    def _add_target(self, field, target_attrib_dict, photometry_attrib_dict,
                    assert_targuse=None):

        for target in field.getElementsByTagName('target'):
            targuse = target.getAttribute('targuse')
            if targuse == 'T':
                first_science_target = target
                break

        new_target = first_science_target.cloneNode(True)

        for key in new_target.attributes.keys():
            assert key in target_attrib_dict.keys()

        self._set_attribs(new_target, target_attrib_dict)

        if assert_targuse is not None:
            assert new_target.getAttribute('targuse') == assert_targuse

        new_target_photometry = new_target.getElementsByTagName('photometry')[0]

        for key in new_target_photometry.attributes.keys():
            assert key in photometry_attrib_dict.keys()

        self._set_attribs(new_target_photometry, photometry_attrib_dict)

        field.insertBefore(new_target, first_science_target)


    def _add_table_as_targets(self, table, assert_targuse=None):

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
                                 photometry_attrib_dict,
                                 assert_targuse=assert_targuse)


    def _set_guide_stars(self, actual_pa, guides_table):

        # Update the PA value if needed

        pa = self._get_pa()

        if actual_pa != pa:
            if not _np.isnan(pa):
                logging.info('Requested value for PA has NOT been adopted')
            self._set_pa(actual_pa)

        # Add the guide stars to the XML file

        if len(guides_table) == 0:
            logging.error('There is not guide stars available')
            raise SystemExit(2)

        self._add_table_as_targets(guides_table, assert_targuse='G')

                
    def _add_guide_stars(self, all_guide_stars=False):

        actual_pa, guides_table = self._get_guide_stars(
            all_guide_stars=all_guide_stars)

        self._set_guide_stars(actual_pa, guides_table)

                        
    def _get_calib_stars(self, mifu_num_calibs=2, all_calib_stars=False):

        obsmode = self._get_obsmode()

        central_ra, central_dec = self._get_central_ra_dec(obsmode)
        pa = self._get_pa()

        calib_stars = _CalibStars(central_ra, central_dec, pa, obsmode)

        full_calibs_table = calib_stars.get_table()

        if all_calib_stars is False:
            calibs_table = full_calibs_table[0:mifu_num_calibs]
        else:
            calibs_table = full_calibs_table

        return calibs_table


    def _set_calib_stars(self, calibs_table):

        # Add the guide stars to the XML file

        if len(calibs_table) == 0:
            logging.error('There is not guide stars available')
            raise SystemExit(2)

        self._add_table_as_targets(calibs_table, assert_targuse='C')


    def _add_calib_stars(self, mifu_num_calibs=2, all_calib_stars=False):

        obsmode = self._get_obsmode()

        if obsmode == 'LIFU':
            # No calibration stars needed!
            return

        calibs_table = self._get_calib_stars(mifu_num_calibs=mifu_num_calibs,
                                             all_calib_stars=all_calib_stars)

        self._set_calib_stars(calibs_table)

        
    def add_guide_and_calib_stars(self, mifu_num_calibs=2,
                                  all_guide_stars=False, all_calib_stars=False):

        self._add_guide_stars(all_guide_stars=all_guide_stars)
        self._add_calib_stars(mifu_num_calibs=mifu_num_calibs,
                              all_calib_stars=all_calib_stars)

