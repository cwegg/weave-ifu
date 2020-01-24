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


import xml.dom.minidom

import numpy as np


def _get_lookup():

    lookup = {}

    lookup['target:cname'] = 'CNAME'
    lookup['target:targsrvy'] = 'TARGSRVY'
    lookup['target:targprog'] = 'TARGPROG'
    lookup['target:targcat'] = 'TARGCAT'
    lookup['target:targid'] = 'TARGID'
    lookup['target:targname'] = 'TARGNAME'
    lookup['target:targprio'] = 'TARGPRIO'
    lookup['target:targuse'] = 'TARGUSE'
    lookup['target:targclass'] = 'TARGCLASS'
    lookup['observation:progtemp'] = 'PROGTEMP'
    lookup['obsconstraints:obstemp'] = 'OBSTEMP'
    lookup[''] = 'GAIA_ID'
    lookup[''] = 'GAIA_DR'
    lookup['target:targra'] = 'GAIA_RA'
    lookup['target:targdec'] = 'GAIA_DEC'
    lookup['target:targepoch'] = 'GAIA_EPOCH'
    lookup['target:targpmra'] = 'GAIA_PMRA'
    lookup[''] = 'GAIA_PMRA_ERR'
    lookup['target:targpmdec'] = 'GAIA_PMDEC'
    lookup[''] = 'GAIA_PMDEC_ERR'
    lookup['target:targparal'] = 'GAIA_PARAL'
    lookup[''] = 'GAIA_PARAL_ERR'
    lookup[''] = 'HEALPIX'
    lookup['target:ifu_spaxel'] = 'IFU_SPAXEL'
    lookup['observation:pa'] = 'IFU_PA'
    lookup['dithering:apply_dither'] = 'IFU_DITHER'
    lookup['photometry:mag_g'] = 'MAG_G'
    lookup['photometry:emag_g'] = 'MAG_G_ERR'
    lookup['photometry:mag_r'] = 'MAG_R'
    lookup['photometry:emag_r'] = 'MAG_R_ERR'
    lookup['photometry:mag_i'] = 'MAG_I'
    lookup['photometry:emag_i'] = 'MAG_I_ERR'
    lookup['photometry:mag_gg'] = 'GAIA_MAG_G'
    lookup['photometry:emag_gg'] = 'GAIA_MAG_G_ERR'
    lookup['photometry:mag_bp'] = 'GAIA_MAG_BP'
    lookup['photometry:emag_bp'] = 'GAIA_MAG_BP_ERR'
    lookup['photometry:mag_rp'] = 'GAIA_MAG_RP'
    lookup['photometry:emag_rp'] = 'GAIA_MAG_RP_ERR'

    lookup.pop('')

    return lookup


def _get_formats():

    formats = {}

    formats['target:targprio'] = float
    formats['target:targra'] = float
    formats['target:targdec'] = float
    formats['target:targepoch'] = float
    formats['target:targpmra'] = float
    formats[''] = float # 'GAIA_PMRA_ERR'
    formats['target:targpmdec'] = float
    formats[''] = float # 'GAIA_PMDEC_ERR'
    formats['target:targparal'] = float
    formats[''] = float # 'GAIA_PARAL_ERR'
    formats[''] = int # 'HEALPIX'
    formats['observation:pa'] = float
    formats['dithering:apply_dither'] = int
    formats['photometry:mag_g'] = float
    formats['photometry:emag_g'] = float
    formats['photometry:mag_r'] = float
    formats['photometry:emag_r'] = float
    formats['photometry:mag_i'] = float
    formats['photometry:emag_i'] = float
    formats['photometry:mag_gg'] = float
    formats['photometry:emag_gg'] = float
    formats['photometry:mag_bp'] = float
    formats['photometry:emag_bp'] = float
    formats['photometry:mag_rp'] = float
    formats['photometry:emag_rp'] = float

    formats.pop('')

    return formats


def _get_xml_data(xml_filename):

    xml_data = {}

    dom = xml.dom.minidom.parse(xml_filename)

    root = dom.childNodes[0]

    # programme = root.childNodes[3]
    observation = root.childNodes[5]

    obsconstraints = dom.getElementsByTagName('obsconstraints')[0]
    dithering = dom.getElementsByTagName('dithering')[0]
    target = dom.getElementsByTagName('target')

    # xml_data['dom'] = dom
    # xml_data['root'] = root
    # xml_data['programme'] = programme
    xml_data['observation'] = observation
    xml_data['obsconstraints'] = obsconstraints
    xml_data['dithering'] = dithering
    xml_data['target'] = target

    return xml_data


def get_spa_data_of_targets_from_xmls(xml_filename_list):
    """
    Get SPA data of the targets contained in a list of configure XML files.

    Parameters
    ----------
    xml_filename_list: list of str
        A list of configure XML files.

    Returns
    -------
    data_dict : dict
        A dictionary with the data. Its keys are the name of the SPA columns
        which are expected to be potencially in the configure XML files, while
        its values are lists contain these data.
    """

    # Get dictionaries with the lookup information and the formats

    lookup = _get_lookup()
    formats = _get_formats()

    # Create a dictionary with the desired keywords and empty lists in order to
    # get ready to store them the data

    data_dict = {lookup[key]: [] for key in lookup.keys()}

    # For each XML file
    
    for xml_filename in xml_filename_list:

        # Get the data from the XML file

        xml_data = _get_xml_data(xml_filename)

        # For each target in the XML

        for target in xml_data['target']:

            # Skip the target if it is a guide or sky fibre

            if str(target.getAttribute('targuse')) != 'T':
                continue

            # For each key in the lookup dictionary
            
            for key in lookup.keys():
            
                # Get the column, element and attribute name
                
                col_name = lookup[key]
                element_name, attribute_name = key.split(':')
                
                # It the name of the element is 'target', choose as the element
                # the current target; otherwise, get the element from the
                # xml_data
                
                if element_name == 'target':
                    element = target
                elif element_name == 'photometry':
                    element = target.getElementsByTagName('photometry')[0]
                else:
                    element = xml_data[element_name]
                
                # Get the raw value from the XML

                raw_value = str(element.getAttribute(attribute_name))

                if key in formats.keys():
                    
                    format_func = formats[key]
                    
                    if raw_value not in ['', '%%%']:
                        formatted_value = format_func(raw_value)
                    elif (raw_value == '') and (format_func == float):
                        formatted_value = np.nan
                    elif (raw_value == '%%%') and (format_func == float):
                        formatted_value = np.nan
                    else:
                        raise ValueError(
                                  'raw value {}, format {}'.format(
                                  raw_value, format_func))
                else:
                    formatted_value = raw_value

                data_dict[col_name].append(formatted_value)

    return data_dict

