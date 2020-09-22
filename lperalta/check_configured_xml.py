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
import logging
import xml.dom.minidom as minidom

import ipdb


def get_obsmode_and_fields(filename):

    dom = minidom.parse(filename)

    observation = dom.getElementsByTagName('observation')[0]

    obsmode = observation.getAttribute('obs_mode')

    fields = dom.getElementsByTagName('fields')[0]
    
    field_list = fields.getElementsByTagName('field')
    
    return obsmode, field_list


def check_attributes_of_targets(field_list):

    expected_attribs = ['cname', 'targcat', 'targclass', 'targdec', 'targepoch',
                        'targid', 'targname', 'targparal', 'targpmdec',
                        'targpmra', 'targprio', 'targprog', 'targra',
                        'targsrvy', 'targuse',
                        'automatic', 'configid', 'fibreid', 'ifu_pa',
                        'ifu_spaxel', 'targx', 'targy']

    for i, field in enumerate(field_list):
    
        for j, target in enumerate(field.getElementsByTagName('target')):
        
            # Get the attributes of the target
        
            attrib_list = list(target.attributes.keys())
            attrib_list.sort()
            
            # Guess the missing attributes
            
            missing_attribs = []
            
            for attrib in expected_attribs:
                if attrib not in attrib_list:
                    missing_attribs.append(attrib)
            
            if len(missing_attribs) > 0:
                logging.error(
                    'missing attributes in target {} of field {}: {}'.format(
                    j + 1, i + 1, missing_attribs))
            
            # Guess the unexpected attributes
            
            unexpected_attribs = []
            
            for attrib in attrib_list:
                if attrib not in expected_attribs:
                    unexpected_attribs.append(attrib)
            
            if len(unexpected_attribs) > 0:
                logging.error(
                    'unexpected attributes in target {} of field {}: {}'.format(
                    j + 1, i + 1, unexpected_attribs))


def check_attributes_of_targets_allowing_missing(obsmode, field_list):

    expected_attribs = ['cname', 'targcat', 'targclass', 'targdec', 'targepoch',
                        'targid', 'targname', 'targparal', 'targpmdec',
                        'targpmra', 'targprio', 'targprog', 'targra',
                        'targsrvy', 'targuse',
                        'automatic', 'configid', 'fibreid', 'ifu_pa',
                        'ifu_spaxel', 'targx', 'targy']
    
    lifu_guide_flag = False
    lifu_central_flag = False
    mifu_non_allocated_flag = False
    mifu_guide_flag = False
    mifu_central_flag = False

    for i, field in enumerate(field_list):
    
        for j, target in enumerate(field.getElementsByTagName('target')):
        
            # Get the attributes of the target
        
            attrib_list = list(target.attributes.keys())
            attrib_list.sort()
            
            # Guess the missing attributes
            
            missing_attribs = []
            
            for attrib in expected_attribs:
                if attrib not in attrib_list:
                    missing_attribs.append(attrib)
            
            if obsmode == 'LIFU':
                
                # In LIFU mode, automatic attribute is not added to guide stars
                
                if target.getAttribute('targuse') == 'G':
                
                    missing_attribs.remove('automatic')
                    
                    lifu_guide_flag = True
                
                # In LIFU mode, automatic attribute is not added to central
                # fibre
            
                if target.getAttribute('ifu_spaxel') == 'C14':
                
                    missing_attribs.remove('automatic')
                    
                    lifu_central_flag = True
            
            elif obsmode == 'mIFU':
            
                # In mIFU mode, non-allocated fibres will not contain some
                # attributes
                
                if 'fibreid' in missing_attribs:
                
                    missing_attribs.remove('automatic')
                    missing_attribs.remove('fibreid')
                    missing_attribs.remove('ifu_pa')
                    missing_attribs.remove('ifu_spaxel')
                    
                    mifu_non_allocated_flag = True
                    
                else:
                
                    # In mIFU mode, automatic, ifu_pa and ifu_spaxel are not
                    # added to guide stars
                    
                    if target.getAttribute('targuse') == 'G':
                    
                        missing_attribs.remove('automatic')
                        missing_attribs.remove('ifu_pa')
                        missing_attribs.remove('ifu_spaxel')
                        
                        mifu_guide_flag = True
                
                    # In mIFU mode, automatic is not added to central spaxel
                    # with allocated calib stars and targets
                    
                    if ((target.getAttribute('targuse') in ['C', 'T']) and
                        (target.getAttribute('ifu_spaxel')[-3:] == 'C04')):
                        
                        missing_attribs.remove('automatic')
                        
                        mifu_central_flag = True
            
            if len(missing_attribs) > 0:
                logging.error(
                    'missing attributes in target {} of field {}: {}'.format(
                    j + 1, i + 1, missing_attribs))
            
            # Guess the unexpected attributes
            
            unexpected_attribs = []
            
            for attrib in attrib_list:
                if attrib not in expected_attribs:
                    unexpected_attribs.append(attrib)
            
            if len(unexpected_attribs) > 0:
                logging.error(
                    'unexpected attributes in target {} of field {}: {}'.format(
                    j + 1, i + 1, unexpected_attribs))

    if lifu_guide_flag is True:
        logging.warning(
            'LIFU guide stars do not include the attribute ' +
            'automatic'
        )

    if lifu_central_flag is True:
        logging.warning(
            'LIFU central target does not include the attribute ' +
            'automatic'
        )

    if mifu_non_allocated_flag is True:
        logging.warning(
            'mIFU non-allocated targets do not include the attribute ' +
            'automatic, fibreid, ifu_pa and ifu_spaxel'
        )

    if mifu_guide_flag is True:
        logging.warning(
            'mIFU guide stars do not include the attribute ' +
            'automatic, ifu_pa and ifu_spaxel'
        )

    if mifu_central_flag is True:
        logging.warning(
            'mIFU central targets do not include the attribute ' +
            'automatic'
        )


def check_attributes_of_photometry(field_list):

    expected_attribs = [
        'mag_g', 'emag_g', 'mag_r', 'emag_r','mag_i', 'emag_i',
        'mag_gg', 'emag_gg', 'mag_bp', 'emag_bp', 'mag_rp', 'emag_rp']

    for i, field in enumerate(field_list):
    
        for j, target in enumerate(field.getElementsByTagName('target')):
        
            # Get the targuse of the target
        
            targuse = target.getAttribute('targuse')
            
            # Get the photometry elements
            
            photometry_list = target.getElementsByTagName('photometry')
            
            # Depending on targuse, we expect to contain a photometry element
            # or not
            
            if targuse in ['T', 'G', 'C']:
            
            
                if len(photometry_list) == 0:
                    logging.error(
                       'missing photometry in target ' +
                        '{} of field {}: {}'.format(j + 1, i + 1))
                elif len(photometry_list) > 1:
                    logging.error(
                        'unexpected photometry in target ' +
                        '{} of field {}: {}'.format(j + 1, i + 1))
            
                photometry = photometry_list[0]
            
                # Get the attributes of the target
            
                attrib_list = list(photometry.attributes.keys())
                attrib_list.sort()
                
                # Guess the missing attributes
                
                missing_attribs = []
                
                for attrib in expected_attribs:
                    if attrib not in attrib_list:
                        missing_attribs.append(attrib)
                
                if len(missing_attribs) > 0:
                    logging.error(
                        'missing attributes in photometry of target ' +
                        '{} of field {}: {}'.format(j + 1, i + 1,
                                                    missing_attribs))
                
                # Guess the unexpected attributes
                
                unexpected_attribs = []
                
                for attrib in attrib_list:
                    if attrib not in expected_attribs:
                        unexpected_attribs.append(attrib)
                
                if len(unexpected_attribs) > 0:
                    logging.error(
                        'unexpected attributes in photometry of target ' +
                        '{} of field {}: {}'.format(j + 1, i + 1,
                                                    unexpected_attribs))
            
            elif targuse in ['S', 'R']:
            
                if len(photometry_list) > 0:
                    logging.error(
                        'unexpected photometry in target ' +
                        '{} of field {}: {}'.format(j + 1, i + 1))
            
            else:
                logging.error(
                    'unexpected targuse in target {} of field {}: {}'.format(
                        j + 1, i + 1, targuse))

                
def check_inherited_values_in_mifu(field_list):

    inherited_attribs = [
        'targcat', 'targclass', 'targepoch', 'targid', 'targname', 'targparal',
        'targpmdec', 'targpmra', 'targprio', 'targprog', 'targsrvy', 'ifu_pa'
    ]

    # Get a dictionary with the central target information

    central_target_attrib_dict = {}

    for target in field_list[0].getElementsByTagName('target'):

        targuse = target.getAttribute('targuse')

        if (targuse != 'G') and ('fibreid' in target.attributes.keys()):

            ifu_spaxel = target.getAttribute('ifu_spaxel')

            prefix = ifu_spaxel[:-3]
            suffix = ifu_spaxel[-3:]

            if suffix == 'C04':

                central_target_attrib_dict[prefix] = {
                    key: target.getAttribute(key)
                    for key in target.attributes.keys()}

                if target.getAttribute('targuse') == 'T':
                    assert (central_target_attrib_dict[prefix]['targclass'] ==
                            '%%%')
                elif target.getAttribute('targuse') == 'S':
                    assert (central_target_attrib_dict[prefix]['cname'] ==
                            '%%%')

                    # assert (central_target_attrib_dict[prefix]['targsrvy'] ==
                    #         '')

                    assert (central_target_attrib_dict[prefix]['targprog'] ==
                            '')
                    assert (central_target_attrib_dict[prefix]['targcat'] ==
                            '')
                    assert (central_target_attrib_dict[prefix]['targid'] ==
                            '')
                    assert (central_target_attrib_dict[prefix]['targname'] ==
                            'auto generated mIFU target')
                    assert (central_target_attrib_dict[prefix]['targprio'] ==
                            '1')

                    assert (central_target_attrib_dict[prefix]['targclass'] ==
                            'UNKNOWN')
                    logging.warning(
                        'UNKNOWN != "" ({})'.format(
                            central_target_attrib_dict[prefix]['targclass']))

                    assert (central_target_attrib_dict[prefix]['targepoch'] ==
                            '')
                    assert (central_target_attrib_dict[prefix]['targpmra'] ==
                            '0')
                    assert (central_target_attrib_dict[prefix]['targpmdec'] ==
                            '0')
                    assert (central_target_attrib_dict[prefix]['targparal'] ==
                            '0')

    # Check the values of all targets

    targprio_flag = False
    automatic_flag = False

    for i, field in enumerate(field_list):

        for j, target in enumerate(field.getElementsByTagName('target')):

            if 'fibreid' in target.attributes.keys():

                targuse = target.getAttribute('targuse')

                if targuse != 'G':

                    ifu_spaxel = target.getAttribute('ifu_spaxel')

                    prefix = ifu_spaxel[:-3]
                    suffix = ifu_spaxel[-3:]

                    if targuse != 'C':
                        assert target.getAttribute('cname') == '%%%'

                    if suffix == 'C04':
                        if (i == 0) and (targuse in ['T', 'C']):
                            if target.getAttribute('automatic') != '0':
                                if target.getAttribute('automatic') == '':
                                    automatic_flag = True
                                else:
                                    raise ValueError
                        else:
                            if target.getAttribute('automatic') != '1':
                                if target.getAttribute('automatic') == '':
                                    automatic_flag = True
                                else:
                                    raise ValueError
                    else:
                        assert target.getAttribute('automatic') == '1'

                    for key in inherited_attribs:
                        try:
                            assert (target.getAttribute(key) ==
                                    central_target_attrib_dict[prefix][key])
                        except:
                            if key == 'targprio':
                                targprio_flag = True
                            else:
                                raise ValueError


                    if central_target_attrib_dict[prefix]['targuse'] == 'T':
                        assert targuse == 'T'
                    elif central_target_attrib_dict[prefix]['targuse'] == 'S':
                        assert targuse == 'S'
                    else:
                        assert targuse in ['C', 'R', 'S']

                else:

                    if 'automatic' in target.attributes.keys():
                        assert target.getAttribute('automatic') == '0'

                    assert target.getAttribute('ifu_spaxel') == ''

                int(target.getAttribute('fibreid'))
                
            int(target.getAttribute('configid'))
            float(target.getAttribute('targra'))
            float(target.getAttribute('targdec'))
            float(target.getAttribute('targx'))
            float(target.getAttribute('targy'))

    if targprio_flag == True:
        logging.error('targprio is not being inherited in mIFU bundles')

    if automatic_flag == True:
        logging.error('automatic is None in mIFU fibres')


def check_inherited_values(obsmode, field_list):

    if obsmode == 'LIFU':

        check_inherited_values_in_lifu(field_list)

    elif obsmode == 'mIFU':

        check_inherited_values_in_mifu(field_list)

    else:

        raise ValueError


def check_configured_xml(filename, allow_missing=True):

    obsmode, field_list = get_obsmode_and_fields(filename)
    
    if allow_missing is False:
        check_attributes_of_targets(field_list)
    else:
        check_attributes_of_targets_allowing_missing(obsmode, field_list)
    
    check_attributes_of_photometry(field_list)
    
    check_inherited_values(obsmode, field_list)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
             description='Check the contents of configured XML files')

    parser.add_argument('xml_file',
                        help="""a configured OB XML file""")

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}

    logging.basicConfig(level=level_dict[args.log_level])

    check_configured_xml(args.xml_file)
