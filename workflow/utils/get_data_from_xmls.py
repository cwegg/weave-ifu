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


def _get_value_from_xml_data(xml_data, target, key, formats):

    # Get the column, element and attribute name
    
    element_name, attribute_name = key.split(':')
    
    # It the name of the element is 'target', choose as the element
    # the current target; otherwise, get the element from the
    # xml_data
    
    if element_name == 'target':
        element = target
    elif element_name == 'photometry':
        photometry_element_list = target.getElementsByTagName('photometry')
        
        # If the photometry element exist, select it and continue
        
        if len(photometry_element_list) > 0:
            element = photometry_element_list[0]
        
        # If there is not a photometry element (like for the sky fibres), the
        # key should correspond to a float in the formats dictionary, so a NaN
        # will be returned
        
        else:
            
            assert (formats[key] == float)
            
            formatted_value = np.nan
            
            return formatted_value
        
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
                      'raw value {}, format {}'.format(raw_value, format_func))
    else:
        formatted_value = raw_value
    
    return formatted_value


def _get_spa_data_from_targuse(xml_filename_list, targuse_list=['T', 'S'],
                               replace_triple_percent=True):

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

            if str(target.getAttribute('targuse')) not in targuse_list:
                continue

            # For each key in the lookup dictionary
            
            for key in lookup.keys():
                
                # Get its column name and its value and save them
                
                col_name = lookup[key]
                
                value = _get_value_from_xml_data(xml_data, target, key, formats)
                
                # Replace the triple percent if requested
                
                if (replace_triple_percent is True) and (value == '%%%'):
                    value = ''

                # Save the value in the corresponding column
                
                data_dict[col_name].append(value)

    return data_dict


def get_spa_data_of_target_fibres_from_xmls(xml_filename_list,
                                            replace_triple_percent=True):
    """
    Get SPA data of the target fibres contained in a list of XML files.

    Parameters
    ----------
    xml_filename_list : list of str
        A list of configure XML files.
    replace_triple_percent : bool, optional
        It will replace the values of '%%%' by '' also for data of type string.

    Returns
    -------
    data_dict : dict
        A dictionary with the data. Its keys are the name of the SPA columns
        which are expected to be potencially in the configure XML files, while
        its values are lists containing these data.
    """
    data_dict = \
        _get_spa_data_from_targuse(xml_filename_list, targuse_list=['T'],
                                  replace_triple_percent=replace_triple_percent)

    return data_dict


def get_spa_data_of_target_and_sky_fibres_from_xmls(xml_filename_list,
                                                   replace_triple_percent=True):
    """
    Get SPA data of the target and sky fibres contained in a list of XML files.

    Parameters
    ----------
    xml_filename_list : list of str
        A list of configure XML files.
    replace_triple_percent : bool, optional
        It will replace the values of '%%%' by '' also for data of type string.

    Returns
    -------
    data_dict : dict
        A dictionary with the data. Its keys are the name of the SPA columns
        which are expected to be potencially in the configure XML files, while
        its values are lists containing these data.
    """
    data_dict = \
        _get_spa_data_from_targuse(xml_filename_list, targuse_list=['T', 'S'],
                                  replace_triple_percent=replace_triple_percent)

    return data_dict


def _get_attribute_in_simple_element_of_xml_file(xml_file, element_name,
                                                 attribute_name):
    
    dom = xml.dom.minidom.parse(xml_file)
    
    element = dom.getElementsByTagName(element_name)[0]
    
    attribute = element.getAttribute(attribute_name)
    
    return attribute


def _get_single_value_from_list(value_list):

    value_set = set(value_list)
    
    assert len(value_set) == 1
    
    value = value_list[0]
    
    return value


def get_trimester_from_xmls(input_xmls):
    """
    Get the trimester from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    trimester : str
        The trimester present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    trimester_list = [_get_attribute_in_simple_element_of_xml_file(
                          xml_file, 'observation', 'trimester')
                      for xml_file in xml_filename_list]
    
    trimester = _get_single_value_from_list(trimester_list)
    
    return trimester


def get_author_from_xmls(input_xmls):
    """
    Get the author from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    author : str
        The author present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    author_list = [_get_attribute_in_simple_element_of_xml_file(
                       xml_file, 'root', 'author')
                   for xml_file in xml_filename_list]
    
    author = _get_single_value_from_list(author_list)
    
    return author


def get_cc_report_from_xmls(input_xmls):
    """
    Get the cc_report from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    cc_report : str
        The cc_report present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    cc_report_list = [_get_attribute_in_simple_element_of_xml_file(xml_file,
                          'root', 'cc_report')
                      for xml_file in xml_filename_list]
    
    cc_report = _get_single_value_from_list(cc_report_list)
    
    return cc_report


def get_report_verbosity_from_xmls(input_xmls):
    """
    Get the report_verbosity from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    report_verbosity : str
        The report_verbosity present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    report_verbosity_list = [_get_attribute_in_simple_element_of_xml_file(
                                 xml_file, 'root', 'report_verbosity')
                      for xml_file in xml_filename_list]
    
    report_verbosity = _get_single_value_from_list(report_verbosity_list)
    
    report_verbosity = int(report_verbosity)
    print(report_verbosity)
    
    return report_verbosity


