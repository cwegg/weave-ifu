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


import argparse
import glob
import logging
import os
import re
import xml.dom.minidom
from collections import OrderedDict

import numpy as np
from astropy.io import fits

from workflow.utils.get_progtemp_info import get_obsmode_from_progtemp
from workflow.utils.get_resources import get_blank_xml_template


class _OBXML:

    
    def __init__(self, xml_template):

        # Save the input filename
        
        self.xml_template = xml_template

        # Parse the XML template

        try:

            self.dom = xml.dom.minidom.parse(self.xml_template)

        except:

            logging.error('File {} would not parse'.format(self.xml_template))

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

        for attrib in configure_attrib_list:
            if attrib in self.configure.attributes.keys():
                self.configure.removeAttribute(attrib)

        # Remove some elements

        elem_list = [
            'optical_axis', 'distortion_coefficients',
            'telescope', 'focal_plane_map', 'hour_angle_limits', 'offsets'
        ]

        self._remove_elem_list_by_name(elem_list)

                
    def remove_non_used_elements(self, elem_list=['simulation',
                                                  'avoidance_list', 'group']):

        self._remove_elem_list_by_name(elem_list)


    def _set_attribs(self, elem, attrib_dict):

        for key in attrib_dict.keys():

            if key in elem.attributes.keys():

                value = str(attrib_dict[key])
                elem.setAttribute(key, value=str(value))

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
                try:
                    if comment_string in node.data:
                        self.exposures.removeChild(node)
                except:
                    pass

    
    def set_observation(self, attrib_dict):

        self._set_attribs(self.observation, attrib_dict)


    def set_configure(self, attrib_dict):

        self._set_attribs(self.configure, attrib_dict)

        
    def set_obsconstraints(self, obstemp):

        logging.warning('Obsconstraints TBDeveloped')

    
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

            
    def set_fields(self, obsmode, entry_group, targcat):
        logging.warning('TBR 3')
        
        #now generate field template
        fields_clone = self.fields.cloneNode(True)

        #clear it
        field_list = fields_clone.getElementsByTagName('field')
        field_clone = field_list[0].cloneNode(True)
        for field in field_list:
            fields_clone.removeChild(field)

        target_list = field_clone.getElementsByTagName('target')
        node_list = field_clone.childNodes

        while len(field_clone.childNodes) != 0:
            for target in field_clone.childNodes:
                field_clone.removeChild(target)

        #now remove the existing <fields> element
        self.observation.removeChild(self.fields)

        ##########

        base_target = self.fields.getElementsByTagName('target')[0]

        #assembly order:
        #1.target
        #2.field
        #3.fields
        order = self.first_science_order
        field = None
        col_names = entry_group[0].array.columns.names
        for i in range(len(entry_group)):
            row = entry_group[i]
            target = base_target.cloneNode(True)

            # remove things we shouldn't have as the xml gets passed into
            # Configure
            to_remove = ['configid', 'fibreid']
            for rem in to_remove:
                target.removeAttribute(rem)

            _row = {}
            for col in col_names:
                _row[col] = row[col]
                if (str(row[col]) == 'nan') and np.isnan(row[col]):
                    _row[col] = ''
            row = _row
             
            target.setAttribute('targra',value=str(row['GAIA_RA']))
            target.setAttribute('targdec',value=str(row['GAIA_DEC']))
            target.setAttribute('targpmra',value=str(row['GAIA_PMRA']))
            target.setAttribute('targpmdec',value=str(row['GAIA_PMDEC']))
            target.setAttribute('targepoch',value=str(row['GAIA_EPOCH']))
            target.setAttribute('targparal',value=str(row['GAIA_PARAL']))
            target.setAttribute('targid',value=str(row['TARGID']))
            target.setAttribute('targname',value=str(row['TARGNAME']))
            target.setAttribute('targprog',value=str(row['TARGPROG']))
            target.setAttribute('targcat',value=targcat)
            target.setAttribute('targprio',value=str(row['TARGPRIO']))
            target.setAttribute('targuse',value='T')
            target.setAttribute('targsrvy',value=str(row['TARGSRVY']))

            if (obsmode == 'LIFU'):
                #generate a <field> element to add this target to
                field = field_clone.cloneNode(True)
                field.setAttribute('order',value=str(order))
                field.setAttribute('RA_d',value=str(row['GAIA_RA']))
                field.setAttribute('Dec_d',value=str(row['GAIA_DEC']))
                field.appendChild(target)
                #now add this field to the fields element
                fields_clone.appendChild(field)

                order += 1

            elif (obsmode == 'mIFU'):
                #if this is the first of the targets, then generate the field,
                #otherwise use the existing <field>

                #remember - we require fixed IFU_DITHER patterns for mIFU, so
                #multiple <field> entries not needed
                if i == 0:
                    field = field_clone.cloneNode(True)
                    field.setAttribute('order',value=str(order))
                    field.setAttribute('RA_d',value=str(row['GAIA_RA']))
                    field.setAttribute('Dec_d',value=str(row['GAIA_DEC']))
                field.appendChild(target)

        if (obsmode == 'mIFU'):
            #now add this field to the fields element
            fields_clone.appendChild(field)
            

        #...and finally to observation
        self.observation.appendChild(fields_clone)

                
    def _remove_empty_lines(self, xml_text):

        parsed_xml = xml.dom.minidom.parseString(xml_text)

        pretty_xml = parsed_xml.toprettyxml(indent=' ')

        clean_xml_text = '\n'.join([line for line in pretty_xml.split('\n')
                                    if line.strip()])

        return clean_xml_text


    def _remove_xml_declaration(self, xml_text):

        parsed_xml = xml.dom.minidom.parseString(xml_text)

        root = parsed_xml.documentElement

        clean_xml_text = root.toxml(parsed_xml.encoding)

        return clean_xml_text


    def write_xml(self, filename, remove_empty_lines=True, declaration=False):

        pretty_xml = self.dom.toprettyxml()

        if remove_empty_lines is True:
            pretty_xml = self._remove_empty_lines(pretty_xml)

        if declaration is False:
            pretty_xml = self._remove_xml_declaration(pretty_xml)

        logging.info('\tWriting to {}'.format(filename))

        with open(filename, 'w') as f:
            f.write(pretty_xml)    


class _IFUDriverCat:
    """
    Convert the target-level data from an IFU driver catalogue to a set of XMLs.
    
    This class provides code to take an IFU driver catalogue containing
    information of the central fibres and convert it into XML targets that are
    written into a supplied XML.
    
    Parameters
    ----------
    filename : str
        A FITS file containing an IFU driver cat.
    targcat : str, optional
        Filename of the FITS catalogue which will be submitted to WASP (if None,
        it will be derived from the input filename).
    """

    
    def __init__(self, filename, targcat=None):

        # Save the filename

        self.filename = filename

        # Set the default prefix to be used for output XMLs

        self.default_prefix =os.path.splitext(
            os.path.basename(self.filename))[0].replace('-ifu_driver_cat', '')

        if self.default_prefix[-1] != '-':
            self.default_prefix = self.default_prefix + '-'

        # Set the targcat value

        if targcat is None:
            self.targcat = os.path.basename(self.filename).replace(
                '-ifu_driver_cat', '')
        else:
            self.targcat = targcat

        # Read and save the information into the IFU driver catalogue
        
        logging.info('Reading {}'.format(self.filename))

        with fits.open(self.filename) as hdu_list:

            self.datamver = hdu_list[0].header['DATAMVER']
            self.trimester = hdu_list[0].header['TRIMESTE']
            self.report_verbosity = hdu_list[0].header['VERBOSE']
            self.author = hdu_list[0].header['AUTHOR']
            self.cc_report = hdu_list[0].header['CCREPORT']

            self.data = hdu_list[1].data

        
    def _get_previous_files(self, output_dir='', prefix=None, suffix='-t'):

        if prefix is None:
            prefix = self.default_prefix

        glob_pattern = os.path.join(output_dir, '{}*{}.xml'.format(prefix,
                                                                   suffix))
        regex = '^{}[lm]ifu_[0-9]+{}.xml$'.format(prefix, suffix)

        candidate_list = glob.glob(glob_pattern)
        candidate_list.sort()

        filename_list = [filename for filename in candidate_list
                        if re.match(regex, os.path.basename(filename))]

        return filename_list

    
    def remove_xmls(self, output_dir='', prefix=None, suffix='-t'):

        prev_filename_list = self._get_previous_files(output_dir=output_dir,
                                                      prefix=prefix,
                                                      suffix=suffix)

        if len(prev_filename_list) > 0:

            logging.info('Removing previous files: {}'.format(
                prev_filename_list))

            for filename in prev_filename_list:
                os.remove(filename)


    def _get_output_path(self, obsmode, output_dir='', prefix=None,
                         suffix='-t'):

        # Set the prefix of the output file if it has not been provided

        if prefix is None:
            prefix = self.default_prefix

        # Choose the first filename which does not exist

        index = 0
        output_path = ''
        
        while (index < 1) or os.path.exists(output_path):

            index += 1

            output_basename = '{}{}_{:02d}{}.xml'.format(
                prefix, obsmode.lower(), index, suffix)

            output_path = os.path.join(output_dir, output_basename)

        return output_path

    
    def _process_ob(self, entry_group, xml_template, output_dir='',
                    prefix=None, suffix='-t', pass_datamver=False):

        # Get some information from the first entry of the group

        first_entry = entry_group[0]

        progtemp = first_entry['PROGTEMP']
        obstemp = first_entry['OBSTEMP']
        ifu_dither = first_entry['IFU_DITHER']
        ifu_pa_request = first_entry['IFU_PA_REQUEST']

        # Get some information from all the entries in the group

        targsrvy_list = list(set([entry['TARGSRVY'] for entry in entry_group]))
        targsrvy_list.sort()

        # Guess the OBSMODE from PROGTEMP

        obsmode = get_obsmode_from_progtemp(progtemp)

        assert obsmode in ['LIFU', 'mIFU']

        # Set the some paremeters which depends on OBSMODE:
        #  - Name of the observation
        #  - Plate
        #  - max_fibres

        if obsmode == 'LIFU':
            observation_name = '{}-{}'.format(first_entry['TARGNAME'],
                                              first_entry['TARGID'])
            plate = 'LIFU'
            max_calibration = 547
            max_guide = 1
            max_sky = 56
            num_sky_fibres = 603
            max_fibres = 603
        elif obsmode == 'mIFU':
            observation_name = first_entry['TARGNAME']
            plate = 'PLATE_B'
            max_calibration = 37
            max_guide = 8
            max_sky = 73
            num_sky_fibres = 740
            max_fibres = 740

        # Guess the number of dither positions

        if (obsmode == 'LIFU') and (len(entry_group) > 1):
            num_dither_positions = len(entry_group)
        else:
            num_dither_positions = ifu_dither

        # Guess the number of science exposures

        logging.warning('TBD (get it from progtemp)')
        num_science_exposures = num_dither_positions

        assert (num_dither_positions % num_science_exposures == 0)

        # Guess the time of each science exposure

        logging.warning('TBD (get it from progtemp)')
        science_exptime = 1200

        # Guess some parameters from progtemp

        logging.warning('TBD (get it from progtemp)')
        binning_y = 1
        resolution = 'low'
        red_vph = 'VPH1'
        blue_vph = 'VPH1'

        # Set the position angle (if possible)

        if obsmode == 'LIFU':

            if not np.isnan(ifu_pa_request):
                pa = ifu_pa_request
            else:
                pa = None

        elif obsmode == 'mIFU':

            pa = 0.0

        # Create an OB from the XML template

        ob_xml = _OBXML(xml_template)

        # Remove the elements and attributes which are configure outputs
                        
        ob_xml.remove_configure_outputs()

        # Remove the non-used elements
                        
        ob_xml.remove_non_used_elements()

        # Check DATAMVER of the IFU driver cat and the XML template

        xml_datamver = ob_xml.get_datamver()

        if self.datamver != xml_datamver:
            logging.critical(
               'DATAMVER mismatch ({} != {}): Stop unless you are sure!'.format(
                    self.datamver, xml_datamver))

            if pass_datamver == False:
                raise SystemExit(2)

        # Set the attributes of the root element
        
        root_attrib_dict = {
            'author': self.author,
            'cc_report': self.cc_report,
            'report_verbosity': self.report_verbosity
        }

        ob_xml.set_root_attrib(root_attrib_dict)

        # Set the contents of the spectrograph element

        ob_xml.set_spectrograph(binning_y=binning_y, resolution=resolution,
                                red_vph=red_vph, blue_vph=blue_vph)

        # Set the contents of the exposures element

        ob_xml.set_exposures(num_science_exposures, science_exptime)

        # Set the attributes of the observation element

        observation_attrib_dict = {
            'name': observation_name,
            'obstemp': obstemp,
            'obs_type': obsmode,
            'progtemp': progtemp,
            'trimester': self.trimester
        }

        if pa is not None:
            observation_attrib_dict['pa'] = pa

        ob_xml.set_observation(observation_attrib_dict)
            
        # Set the attributes of the dithering element

        configure_attrib_dict = {
            'max_calibration': max_calibration,
            'max_guide': max_guide,
            'max_sky': max_sky,
            'num_sky_fibres': num_sky_fibres,
            'plate': plate
        }

        ob_xml.set_configure(configure_attrib_dict)

        # Set the attributes of the obsconstraints element

        ob_xml.set_obsconstraints(obstemp)

        # Set the attributes of the dithering element

        dithering_attrib_dict = {
            'apply_dither': ifu_dither
        }

        ob_xml.set_dithering(dithering_attrib_dict)

        # Set the contents of the survey element

        ob_xml.set_surveys(targsrvy_list, max_fibres)

        # Set the contents of the fields element

        ob_xml.set_fields(obsmode, entry_group, self.targcat)

        # Write the OB XML to a file

        output_path = self._get_output_path(obsmode, output_dir=output_dir,
                                            prefix=prefix, suffix=suffix)

        ob_xml.write_xml(output_path)

        
    def _generate_lifu_xmls(self, lifu_entry_list, xml_template, output_dir='',
                            prefix='', suffix='-t', pass_datamver=False):

        # How do you group entries belonging to the same OB?
        # (for custom dithers)

        group_id = ('TARGID', 'TARGNAME', 'PROGTEMP', 'OBSTEMP')

        # Detect the non-custom dithers and save them into a list

        non_custom_dithers_list = []
        
        for lifu_entry in lifu_entry_list:

            if lifu_entry['IFU_DITHER'] != -1:

                non_custom_dithers_list.append(lifu_entry)

        # Generate the IFU XMLs of the non-custom dithers

        logging.info(
            'Processing {} non-custom dither OBs for LIFU'.format(
                len(non_custom_dithers_list)))

        for lifu_entry in non_custom_dithers_list:

            entry_group = [lifu_entry]

            self._process_ob(entry_group, xml_template, output_dir=output_dir,
                             prefix=prefix, suffix=suffix,
                             pass_datamver=pass_datamver)

        # Detect the custom dithers and save them into a dict

        custom_dithers_dict = OrderedDict()

        for lifu_entry in lifu_entry_list:

            if lifu_entry['IFU_DITHER'] == -1:

                key = tuple(lifu_entry[col] for col in group_id)

                if key not in custom_dithers_dict.keys():
                    custom_dithers_dict[key] = []
                
                custom_dithers_dict[key].append(lifu_entry)

        # Process the custom dither pointing groupings
        
        logging.info(
            'Processing {} custom dither OBs for LIFU'.format(
                len(custom_dithers_dict)))

        for entry_group in custom_dithers_dict.values():
            self._process_ob(entry_group, xml_template, output_dir=output_dir,
                             prefix=prefix, suffix=suffix,
                             pass_datamver=pass_datamver)

                
    def _generate_mifu_xmls(self, mifu_entry_list, xml_template,
                            mifu_mode='aufbau', mifu_num_calibs=2,
                            mifu_num_extra=0, output_dir='', prefix='',
                            suffix='-t', pass_datamver=False):

        # What happens if there are (e.g.) 100 bundles in a given key and
        # mifu_num_calibs is 2?
        #
        # User needs to decide 3 options:
        # - all: Include all
        # - aufbau: Aufbau principle, i.e. fill up to 18, leaving 2 for
        #           calibration bundles) then make a new XML
        # - equipartition: divide up the bundles equally:
        #                  100 / 18 = 5.5 --> 6 OBs
        #                  100 / 6 = 17 bundles per OB
        #                            (with remainder added to a final OB)
        #                  i.e. 5 OBs x 17 bundles + 1 OB x 15 bundles

        assert mifu_mode in ['aufbau', 'equipartition', 'all']
        assert (mifu_num_calibs >= 0) and (mifu_num_calibs < 20)

        # How do you group bundles belonging to the same field?

        group_id = ('TARGNAME', 'PROGTEMP', 'OBSTEMP', 'IFU_DITHER')

        # Group the mIFU entries per field

        fields_dict = OrderedDict()

        for mifu_entry in mifu_entry_list:

            key = tuple(mifu_entry[col] for col in group_id)

            if key not in fields_dict.keys():
                fields_dict[key] = []

            fields_dict[key].append(mifu_entry)

        # Group the targets in each field per OB

        ob_nested_list = []

        for key in fields_dict.keys():

            # Get the entries in the field

            field_entry_list = fields_dict[key]

            # Guess the number of entries in the field

            num_entries_in_field = len(field_entry_list)

            # Create the group of mIFU entries per OB depending on the mode

            if (mifu_mode == 'aufbau') or (mifu_mode == 'equipartition'):

                # Maximum number of target entries per OB

                max_entries_per_ob = 20 - mifu_num_calibs

                if mifu_mode == 'aufbau':

                    # All the OBs (except the last one) will contain the maximum
                    # number of target entries plus the amount of extra targets
                    
                    num_entries_per_ob = max_entries_per_ob + mifu_num_extra

                elif mifu_mode == 'equipartition':

                    # Compute the number of OBs

                    num_obs = int(np.ceil(num_entries_in_field /
                                          max_entries_per_ob))

                    # Divide the targets to be the same (except the last one,
                    # due to rounding effects)

                    num_entries_per_ob = int(np.ceil(num_entries_in_field /
                                                     num_obs))

                else:

                    raise ValueError

                # Save the targets of each OB using Python slices
                # (do not worry, if stop is higher than the lenght of the list,
                #  this is allowed and the trick)
                
                start = 0

                while start < num_entries_in_field:

                    stop = start + num_entries_per_ob

                    ob_nested_list.append(field_entry_list[start:stop])

                    start = stop
                    
            elif mifu_mode == 'all':

                ob_nested_list.append(field_entry_list)

            else:

                raise ValueError

        # Proccess OB grouping

        logging.info(
            'Processing {} mIFU fields into {} OBs according to {} mode'.format(
                len(fields_dict), len(ob_nested_list), mifu_mode))

        for entry_group in ob_nested_list:
            self._process_ob(entry_group, xml_template, output_dir=output_dir,
                             prefix=prefix, suffix=suffix,
                             pass_datamver=pass_datamver)

                    
    def generate_xmls(self, xml_template, mifu_mode='aufbau', mifu_num_calibs=2,
                      mifu_num_extra=0, output_dir='', prefix='', suffix='-t',
                      pass_datamver=False):

        # Classify the entries into LIFU and mIFU
        
        lifu_entry_list = []
        mifu_entry_list = []

        for i, entry in enumerate(self.data):

            progtemp = entry['PROGTEMP']

            try:
                obsmode = get_obsmode_from_progtemp(progtemp)
            except:
                obsmode = None

            if obsmode == 'LIFU':
                lifu_entry_list.append(entry)
            elif obsmode == 'mIFU':
                mifu_entry_list.append(entry)
            else:
                logging.warning('unexpected PROGTEMP in row {}: {}'.format(
                    i + 1, progtemp))

        # Generate the LIFU XMLs
                
        self._generate_lifu_xmls(lifu_entry_list, xml_template,
                                 output_dir=output_dir, prefix=prefix,
                                 suffix=suffix, pass_datamver=pass_datamver)

        # Generate the mIFU XMLs
                
        self._generate_mifu_xmls(mifu_entry_list, xml_template,
                                 mifu_mode=mifu_mode,
                                 mifu_num_calibs=mifu_num_calibs,
                                 mifu_num_extra=mifu_num_extra,
                                 output_dir=output_dir, prefix=prefix,
                                 suffix=suffix, pass_datamver=pass_datamver)


def create_xml_files(ifu_driver_cat_filename, output_dir, xml_template,
                     mifu_mode='aufbau', mifu_num_calibs=2, mifu_num_extra=0,
                     prefix=None, suffix='-t', pass_datamver=False,
                     overwrite=False):
    """
    Create XML files with targets from an IFU driver cat.
    
    Parameters
    ----------
    ifu_driver_cat_filename : str
        A FITS file containing an IFU driver cat.
    output_dir : str
        Name of the directory which will containe the output XML files.
    xml_template : str
        A blank XML template to be populated with the information of the OBs.
    mifu_mode : {'aufbau', 'equipartition', 'all'}, optional
        Grouping mode for mIFU targets.
    mifu_num_calibs : int, optional
        Number of mIFU calibration stars and sky bundles to be considered in
        the mIFU grouping.
    mifu_num_extra : int, optional
        Number of extra mIFU targets to be included in the OBs when 'aufbau'
        grouping mode is used.
    prefix : str, optional
        Prefix to be used in the output files (it will be derived from
        ifu_driver_cat_filename if None is provided).
    suffix : str, optional
        Suffix to be used in the output files.
    pass_datamver : bool, optional
        Continue even if DATAMVER mismatch is detected.
    overwrite : bool, optional
        Overwrite the output FITS file.
    """

    # Check that the input IFU driver cat exists and is a file

    assert os.path.isfile(ifu_driver_cat_filename)

    # Create an object with the IFU driver cat

    ifu_driver_cat = _IFUDriverCat(ifu_driver_cat_filename)

    # Remove the previous files if overwriting has been requested
    
    if overwrite == True:
        ifu_driver_cat.remove_xmls(output_dir=output_dir, prefix=prefix,
                                   suffix=suffix)

    # Create the XML files

    ifu_driver_cat.generate_xmls(xml_template, mifu_mode=mifu_mode,
                                 mifu_num_calibs=mifu_num_calibs,
                                 mifu_num_extra=mifu_num_extra,
                                 output_dir=output_dir,
                                 prefix=prefix, suffix=suffix,
                                 pass_datamver=pass_datamver)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
             description='Create XML files with targets from an IFU driver cat')

    parser.add_argument('ifu_driver_cat',
                        help="""a FITS file containing an IFU driver cat""")

    parser.add_argument('--xml_template', dest='xml_template',
                        default='aux/BlankXMLTemplate.xml',
                        help="""a blank XML template to be populated with the
                        information of the OBs""")

    parser.add_argument('--mifu_mode', dest='mifu_mode', default='aufbau',
                        choices=['aufbau', 'equipartition', 'all'],
                        help="""grouping mode for mIFU targets: 'aufbau'
                        follows the aufbau principle, 'equipartition' divides up
                        the bundles equally, 'all' includes all""")

    parser.add_argument('--mifu_num_calibs', dest='mifu_num_calibs', default=2,
                        choices=range(20), type=int,
                        help="""number of mIFU calibration stars and sky bundles
                        to be considered in the mIFU grouping""")

    parser.add_argument('--mifu_num_extra', dest='mifu_num_extra', default=0,
                        type=int,
                        help="""number of extra mIFU targets to be included in
                        the OBs when 'aufbau' grouping mode is used""")

    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name of the directory which will containe the
                        output XML files""")

    parser.add_argument('--prefix', dest='prefix', default=None,
                        help="""prefix to be used in the output files (it will
                        be derived from ifu_driver_cat if non provided)""")

    parser.add_argument('--pass_datamver', dest='pass_datamver',
                        action='store_true',
                        help='continue even if DATAMVER mismatch is detected')

    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='overwrite the output files')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}

    logging.basicConfig(level=level_dict[args.log_level])

    if not os.path.exists(args.output_dir):
        logging.info('Creating the output directory')
        os.mkdir(args.output_dir)

    xml_template_dir = os.path.dirname(args.xml_template)
    if not os.path.exists(xml_template_dir):
        logging.info('Creating the directory of the blank XML template')
        os.mkdir(xml_template_dir)

    if not os.path.exists(args.xml_template):
        logging.info('Downloading the blank XML template')
        get_blank_xml_template(file_path=args.xml_template)

    create_xml_files(args.ifu_driver_cat, args.output_dir, args.xml_template,
                     mifu_mode=args.mifu_mode,
                     mifu_num_calibs=args.mifu_num_calibs,
                     mifu_num_extra=args.mifu_num_extra,
                     prefix=args.prefix, pass_datamver=args.pass_datamver,
                     overwrite=args.overwrite)

