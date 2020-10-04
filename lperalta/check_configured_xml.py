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


def get_info_from_xml(filename):

    dom = minidom.parse(filename)

    observation = dom.getElementsByTagName('observation')[0]
    obsmode = observation.getAttribute('obs_mode')

    survey = dom.getElementsByTagName('survey')[0]
    survey_name = survey.getAttribute('name')

    fields = dom.getElementsByTagName('fields')[0]
    field_list = fields.getElementsByTagName('field')

    try:
        offsets = dom.getElementsByTagName('offsets')[0]
        step_ra0 = float(offsets.getAttribute('offset_step_ra').split()[0])
        step_dec0 = float(offsets.getAttribute('offset_step_dec').split()[0])
        
        if (step_ra0 == 0.0) and (step_dec0 == 0.0):
            first_field_with_offset = False
        else:
            first_field_with_offset = True
    except:
        first_field_with_offset = None
    
    return obsmode, survey_name, field_list, first_field_with_offset


def check_attributes_of_targets(obsmode, field_list):

    expected_attribs = [
        'cname', 'targcat', 'targclass', 'targdec', 'targepoch', 'targid',
        'targname', 'targparal', 'targpmdec', 'targpmra', 'targprio',
        'targprog', 'targra', 'targsrvy', 'targuse',
        'automatic', 'configid', 'fibreid', 'ifu_pa', 'ifu_spaxel', 'targx',
        'targy']
    
    lifu_guide_automatic_flag = 0
    lifu_central_automatic_flag = 0
    mifu_non_allocated_attribs_flag = 0
    mifu_guide_attribs_flag = 0
    mifu_central_automatic_flag = 0

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
                    
                    lifu_guide_automatic_flag += 1
                
                # In LIFU mode, automatic attribute is not added to central
                # fibre
            
                if target.getAttribute('ifu_spaxel') == 'C14':
                
                    missing_attribs.remove('automatic')
                    
                    lifu_central_automatic_flag += 1
            
            elif obsmode == 'mIFU':
            
                # In mIFU mode, non-allocated fibres will not contain some
                # attributes
                
                if 'fibreid' in missing_attribs:
                
                    missing_attribs.remove('automatic')
                    missing_attribs.remove('fibreid')
                    missing_attribs.remove('ifu_pa')
                    missing_attribs.remove('ifu_spaxel')
                    
                    mifu_non_allocated_attribs_flag += 1
                    
                else:
                
                    # In mIFU mode, automatic, ifu_pa and ifu_spaxel are not
                    # added to guide stars
                    
                    if target.getAttribute('targuse') == 'G':
                    
                        missing_attribs.remove('automatic')
                        missing_attribs.remove('ifu_pa')
                        missing_attribs.remove('ifu_spaxel')
                        
                        mifu_guide_attribs_flag += 1
                
                    # In mIFU mode, automatic is not added to central spaxel
                    # with allocated calib stars and targets
                    
                    if ((target.getAttribute('targuse') in ['C', 'T', 'R']) and
                        (target.getAttribute('ifu_spaxel')[-3:] == 'C04')):
                        
                        missing_attribs.remove('automatic')
                        
                        mifu_central_automatic_flag += 1
            
            if len(missing_attribs) > 0:
            
                try:
                    ifu_spaxel = target.getAttribute('ifu_spaxel')
                except:
                    ifu_spaxel = None
            
                try:
                    targuse = target.getAttribute('targuse')
                except:
                    targuse = None
                
                logging.error(
                    'missing attributes in target {} of field {} '.format(
                        j + 1, i + 1) +
                    '(ifu_spaxel={}, targuse={}): {}'.format(
                        ifu_spaxel, targuse, missing_attribs))
            
            # Guess the unexpected attributes
            
            unexpected_attribs = []
            
            for attrib in attrib_list:
                if attrib not in expected_attribs:
                    unexpected_attribs.append(attrib)
            
            if len(unexpected_attribs) > 0:
                logging.error(
                    'unexpected attributes in target {} of field {}: {}'.format(
                    j + 1, i + 1, unexpected_attribs))

    if lifu_guide_automatic_flag > 0:
        logging.warning(
            'LIFU guide stars do not include the attribute ' +
            'automatic ({})'.format(lifu_guide_automatic_flag))

    if lifu_central_automatic_flag > 0:
        logging.warning(
            'LIFU central target does not include the attribute ' +
            'automatic ({})'.format(lifu_central_automatic_flag))

    if mifu_non_allocated_attribs_flag > 0:
        logging.warning(
            'mIFU non-allocated targets do not include the attributes ' +
            'automatic, fibreid, ifu_pa and ifu_spaxel ({})'.format(
                mifu_non_allocated_attribs_flag))

    if mifu_guide_attribs_flag > 0:
        logging.warning(
            'mIFU guide stars do not include the attributes ' +
            'automatic, ifu_pa and ifu_spaxel ({})'.format(
                mifu_guide_attribs_flag))

    if mifu_central_automatic_flag > 0:
        logging.warning(
            'mIFU central targets do not include the attribute ' +
            'automatic ({})'.format(mifu_central_automatic_flag))


def check_attributes_of_photometry(field_list):

    expected_attribs = [
        'emag_bp', 'emag_g', 'emag_gg', 'emag_i', 'emag_r', 'emag_rp',
        'mag_bp', 'mag_g', 'mag_gg', 'mag_i', 'mag_r', 'mag_rp']

    for i, field in enumerate(field_list):
    
        for j, target in enumerate(field.getElementsByTagName('target')):
        
            # Get the ifu_spaxel and targuse of the target
        
            ifu_spaxel = target.getAttribute('ifu_spaxel')
            targuse = target.getAttribute('targuse')
            
            # Get the photometry elements
            
            photometry_list = target.getElementsByTagName('photometry')
            
            # Depending on targuse, we expect to contain a photometry element
            # or not
            
            if targuse in ['T', 'G', 'C']:
            
                if len(photometry_list) == 0:
                
                    logging.error(
                       'missing photometry in target ' +
                        '{} of field {} (ifu_spaxel={}, targuse={})'.format(
                            j + 1, i + 1, ifu_spaxel, targuse))
                    
                    continue
                    
                elif len(photometry_list) > 1:
                    
                    logging.error(
                        'unexpected photometry in target ' +
                        '{} of field {} (ifu_spaxel={}, targuse={})'.format(
                            j + 1, i + 1, ifu_spaxel, targuse))
            
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
                        '{} of field {} (ifu_spaxel={}, targuse={}): {}'.format(
                            j + 1, i + 1, ifu_spaxel, targuse, missing_attribs))
                
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
                        '{} of field {} (ifu_spaxel={}, targuse={})'.format(
                            j + 1, i + 1, ifu_spaxel, targuse))
            
            else:
                logging.error(
                    'unexpected targuse in target ' +
                    '{} of field {} (ifu_spaxel={}): {}'.format(
                        j + 1, i + 1, ifu_spaxel, targuse))


def check_values_of_attributes_of_photometry(field_list):

    for i, field in enumerate(field_list):

        for j, target in enumerate(field.getElementsByTagName('target')):

            targuse = target.getAttribute('targuse')
            ifu_spaxel = target.getAttribute('ifu_spaxel')
        
            photometry_list = target.getElementsByTagName('photometry')
        
            if targuse in ['G', 'C']:
            
                if len(photometry_list) == 0:
                    
#                    logging.error(
#                       'missing photometry in target ' +
#                        '{} of field {} (ifu_spaxel={}, targuse={})'.format(
#                            j + 1, i + 1, ifu_spaxel, targuse))
                    
                    continue
            
                photometry = photometry_list[0]
                    
                for attrib in photometry.attributes.keys():
                    
                    value = photometry.getAttribute(attrib)
                    
                    if value != '':
                        float(value)
            
            elif targuse == 'T':
            
                photometry = photometry_list[0]
                    
                for attrib in photometry.attributes.keys():
                    
                    value = photometry.getAttribute(attrib)
                    
                    assert value == '%%%'
            
            elif targuse in ['S', 'R']:
            
                pass
                
#                if len(photometry_list) > 0:
#                    logging.error(
#                        'unexpected photometry in target ' +
#                        '{} of field {} (ifu_spaxel={}, targuse={})'.format(
#                            j + 1, i + 1, ifu_spaxel, targuse))
            
            else:
            
                raise ValueError


def check_values_of_attributes_of_lifu_targets(field_list):

    inherited_attribs = [
        'targcat', 'targepoch', 'targid', 'targname', 'targparal',
        'targpmdec', 'targpmra', 'targprio', 'targprog', 'targsrvy', 'ifu_pa']

    # Get a dictionary with the central target information

    for target in field_list[0].getElementsByTagName('target'):

        ifu_spaxel = target.getAttribute('ifu_spaxel')

        if ifu_spaxel == 'C14':

            central_target_attrib_dict = {
                key: target.getAttribute(key)
                for key in target.attributes.keys()}

    # Check the inherited values of all targets

    for i, field in enumerate(field_list):

        for j, target in enumerate(field.getElementsByTagName('target')):

            ifu_spaxel = target.getAttribute('ifu_spaxel')

            if ifu_spaxel != '':
            
                for key in inherited_attribs:
                
                    assert (target.getAttribute(key) ==
                            central_target_attrib_dict[key])
                
            else:
            
                assert (target.getAttribute('ifu_pa') ==
                        central_target_attrib_dict['ifu_pa'])

    # Check attributes which depend on the bundle

    for i, field in enumerate(field_list):

        for j, target in enumerate(field.getElementsByTagName('target')):

            ifu_spaxel = target.getAttribute('ifu_spaxel')

            if ifu_spaxel != '':
            
                assert len(ifu_spaxel) == 3
            
                if ifu_spaxel[0] not in ['R', 'S']:
                
                    assert (target.getAttribute('targuse') == 'T')
                    assert (target.getAttribute('targclass') ==
                            central_target_attrib_dict['targclass'])
                    
                else:
                
                    assert (target.getAttribute('targuse') == 'S')
                    assert (target.getAttribute('targclass') == '')
            
            else:
            
                assert (target.getAttribute('targuse') == 'G')

    # Check that the assigned values by configure are as expected

    for i, field in enumerate(field_list):

        for j, target in enumerate(field.getElementsByTagName('target')):
            
            float(target.getAttribute('targdec'))
            float(target.getAttribute('targra'))
            
            int(target.getAttribute('configid'))
            
            float(target.getAttribute('targx'))
            float(target.getAttribute('targy'))
        
            ifu_spaxel = target.getAttribute('ifu_spaxel')
            
            if ifu_spaxel != '':
                int(target.getAttribute('fibreid'))
            else:
                assert (target.getAttribute('fibreid') == '')
            
            if ifu_spaxel != '' and ifu_spaxel != 'C14':
                assert (target.getAttribute('automatic') == '1')
            elif ifu_spaxel == 'C14':
                if 'automatic' in target.attributes.keys():
                    assert (target.getAttribute('automatic') == '0')
            elif ifu_spaxel == '':
                if 'automatic' in target.attributes.keys():
                    assert (target.getAttribute('automatic') == '0')

    # Check that CNAMEs are empty unless in guide stars

    for i, field in enumerate(field_list):

        for j, target in enumerate(field.getElementsByTagName('target')):
        
            if target.getAttribute('targuse') != 'G':
                
                assert (target.getAttribute('cname') == '%%%')
            
            else:
                
                assert (target.getAttribute('cname')[:4] == 'WVE_')


def check_lifu_versus_pre_xml_file(field_list, pre_field_list,
                                   first_field_with_offset):

    for i, field in enumerate(field_list):
        
        if len(pre_field_list) == 1:
            ref_field = pre_field_list[0]
        else:
            ref_field = pre_field_list[i]
            
        for j, pre_target in enumerate(ref_field.getElementsByTagName('target')):
            
            targuse = pre_target.getAttribute('targuse')
            
            if targuse == 'G':
            
                cname = pre_target.getAttribute('cname')
                
                match_found = False
            
                for k, target in enumerate(field.getElementsByTagName('target')):
                    
                    if target.getAttribute('cname') == cname:
                    
                        for attrib in pre_target.attributes.keys():
                        
                            assert (pre_target.getAttribute(attrib) ==
                                    target.getAttribute(attrib))
                        
                        match_found = True
                        break
                
                if match_found == False:
                    logging.error('Match not-found for guiding star')
                
            elif targuse == 'T':
                
                match_found = False
            
                for k, target in enumerate(field.getElementsByTagName('target')):
                    
                    if target.getAttribute('ifu_spaxel') == 'C14':
                    
                        for attrib in pre_target.attributes.keys():
                            
                            if (attrib not in ['targra', 'targdec']):
                                assert (pre_target.getAttribute(attrib) ==
                                        target.getAttribute(attrib))
                            elif (i == 0) and (first_field_with_offset == False): 
                                assert (pre_target.getAttribute(attrib) ==
                                        target.getAttribute(attrib))
                        
                        match_found = True
                        break
                
                if match_found == False:
                    logging.error('Match not-found for central fibre')
            
            else:
            
                raise ValueError


def check_values_of_attributes_of_mifu_targets(field_list, survey_name=None):

    inherited_attribs = [
        'targcat', 'targepoch', 'targid', 'targname', 'targparal', 'targpmdec',
        'targpmra', 'targprio', 'targprog', 'targsrvy', 'ifu_pa']

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

                if target.getAttribute('targuse') in ['T', 'R']:
                    assert (central_target_attrib_dict[prefix]['cname'] ==
                            '%%%')
                    assert (central_target_attrib_dict[prefix]['targclass'] ==
                            '%%%')
                elif target.getAttribute('targuse') == 'S':
                    
                    assert (central_target_attrib_dict[prefix]['cname'] ==
                            '%%%')

                    if survey_name is not None:
                        assert (central_target_attrib_dict[prefix]['targsrvy'] ==
                                survey_name)

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
                            '')

                    assert (central_target_attrib_dict[prefix]['targepoch'] ==
                            '')
                    assert (central_target_attrib_dict[prefix]['targpmra'] ==
                            '0')
                    assert (central_target_attrib_dict[prefix]['targpmdec'] ==
                            '0')
                    assert (central_target_attrib_dict[prefix]['targparal'] ==
                            '0')

    # Check the inherited values of all targets

    for i, field in enumerate(field_list):

        for j, target in enumerate(field.getElementsByTagName('target')):

            if 'fibreid' in target.attributes.keys():
            
                ifu_spaxel = target.getAttribute('ifu_spaxel')
                targuse = target.getAttribute('targuse')
                targclass = target.getAttribute('targclass')

                if len(ifu_spaxel) > 0:
                    prefix = ifu_spaxel[:-3]
                    suffix = ifu_spaxel[-3:]
                    
                    central_targuse = central_target_attrib_dict[prefix]['targuse']
                else:
                    prefix = ''
                    suffix = ''
                
                    central_targuse = 'G'

                if central_targuse in ['T', 'S', 'C']:

                    for key in inherited_attribs:
                        assert (target.getAttribute(key) ==
                                central_target_attrib_dict[prefix][key])

                    if central_targuse in ['T', 'S']:
                        assert targuse == central_targuse
                    elif central_targuse == 'C':
                        # ***Improve this assert
                        if suffix == 'C04':
                            if i == 0:
                                assert targuse == 'C'
                            else:
                                assert targuse == 'R'
                        elif suffix in ['C03', '03C', '03D', 'C05', '05D', '05C']:
                            if targuse != 'R':
                                logging.error(
                                    'unexpected targuse in target ' +
                                    '{} of field {} (ifu_spaxel={}, targuse={})'.format(
                                        j + 1, i + 1, ifu_spaxel, targuse))
                        else:
                            assert targuse == 'S'
                    
                    if targuse == 'T':
                        assert targclass == '%%%'
                    elif targuse == 'S':
                        assert targclass == ''

                elif central_targuse == 'G':

                    pass

                else:
                    
                    raise ValueError

    # Check that the assigned values by configure are as expected

    for i, field in enumerate(field_list):

        for j, target in enumerate(field.getElementsByTagName('target')):

            if 'fibreid' in target.attributes.keys():
            
                float(target.getAttribute('targdec'))
                float(target.getAttribute('targra'))
                
                int(target.getAttribute('fibreid'))
                int(target.getAttribute('configid'))
                
                float(target.getAttribute('targx'))
                float(target.getAttribute('targy'))
            
                targuse = target.getAttribute('targuse')
                ifu_spaxel = target.getAttribute('ifu_spaxel')
                
                if targuse != 'G':
                    assert len(ifu_spaxel) == 6
                else:
                    assert (ifu_spaxel == '')
                
                if ifu_spaxel != '' and ifu_spaxel[-3:] != 'C04':
                    assert (target.getAttribute('automatic') == '1')
                elif (ifu_spaxel[-3:] == 'C04') and (targuse != 'S'):
                    if 'automatic' in target.attributes.keys():
                        assert (target.getAttribute('automatic') == '0')
                elif (ifu_spaxel[-3:] == 'C04') and (targuse == 'S'):
                    if 'automatic' in target.attributes.keys():
                        assert (target.getAttribute('automatic') == '1')
                elif ifu_spaxel == '':
                    if 'automatic' in target.attributes.keys():
                        assert (target.getAttribute('automatic') == '0')

    # Check that CNAMEs are empty unless for guide stars and central fibre of
    # the calib stars

    for i, field in enumerate(field_list):

        for j, target in enumerate(field.getElementsByTagName('target')):

            if 'fibreid' in target.attributes.keys():
            
                targuse = target.getAttribute('targuse')
                ifu_spaxel = target.getAttribute('ifu_spaxel')
        
                if targuse == 'G':
                    
                    assert (target.getAttribute('cname')[:4] == 'WVE_')
            
                elif ((targuse == 'C') and (ifu_spaxel[-3:] == 'C04')):
                    
                    assert (target.getAttribute('cname')[:4] == 'WVE_')
                
                else:
                    
                    if (target.getAttribute('cname') != '%%%'):
                        logging.error(
                            'unexpected cname in target ' +
                            '{} of field {} (ifu_spaxel={}, targuse={})'.format(
                                j + 1, i + 1, ifu_spaxel, targuse))


def check_mifu_versus_pre_xml_file(field_list, pre_field_list):

    for i, field in enumerate(field_list):
        
        ref_field = pre_field_list[0]
            
        for j, pre_target in enumerate(ref_field.getElementsByTagName('target')):
            
            targuse = pre_target.getAttribute('targuse')
            
            if targuse == 'G':
            
                cname = pre_target.getAttribute('cname')
                
                match_found = False
            
                for k, target in enumerate(field.getElementsByTagName('target')):
                    
                    if target.getAttribute('cname') == cname:
                    
                        for attrib in pre_target.attributes.keys():
                        
                            assert (pre_target.getAttribute(attrib) ==
                                    target.getAttribute(attrib))
                        
                        match_found = True
                        break
                
                if match_found == False:
                    logging.error('Match not-found for guiding star')
                
            elif targuse == 'C':
            
                if i == 0:
                
                    cname = pre_target.getAttribute('cname')
                    
                    match_found = False
                
                    for k, target in enumerate(field.getElementsByTagName('target')):
                        
                        if target.getAttribute('cname') == cname:
                        
                            for attrib in pre_target.attributes.keys():
                            
                                assert (pre_target.getAttribute(attrib) ==
                                        target.getAttribute(attrib))
                            
                            match_found = True
                            break
                    
                    if match_found == False:
                        logging.error('Match not-found for calib star')
                
                else:
                
                    targname = pre_target.getAttribute('targname')
                    targid = pre_target.getAttribute('targid')
                    
                    match_found = False
                
                    for k, target in enumerate(field.getElementsByTagName('target')):
                        
                        if ((target.getAttribute('targname') == targname) and
                            (target.getAttribute('targid') == targid) and
                            (target.getAttribute('ifu_spaxel')[-3:] == 'C04')):
                        
                            for attrib in pre_target.attributes.keys():
                            
                                if attrib not in ['cname', 'targdec', 'targra', 'targuse']:
                                    assert (pre_target.getAttribute(attrib) ==
                                            target.getAttribute(attrib))
                                elif attrib == 'cname':
                                    ifu_spaxel = target.getAttribute('ifu_spaxel')
                                    targuse = target.getAttribute('targuse')
                                    if target.getAttribute(attrib) != '%%%':
                                        logging.error(
                                            'unexpected cname in target ' +
                                            '{} of field {} (ifu_spaxel={}, targuse={})'.format(
                                                k + 1, i + 1, ifu_spaxel, targuse))
                                elif attrib == 'targuse':
                                    assert (target.getAttribute(attrib) != 'C')
                            
                            match_found = True
                            break
                    
                    if match_found == False:
                        logging.error('Match not-found for calib star')
                
            elif targuse == 'T':
            
                if i == 0:
                
                    targname = pre_target.getAttribute('targname')
                    targid = pre_target.getAttribute('targid')
                    
                    match_found = False
                
                    for k, target in enumerate(field.getElementsByTagName('target')):
                        
                        if ((target.getAttribute('targname') == targname) and
                            (target.getAttribute('targid') == targid) and
                            ((target.getAttribute('ifu_spaxel')[-3:] == 'C04') or
                             (target.getAttribute('ifu_spaxel') == ''))):
                        
                            for attrib in pre_target.attributes.keys():
                                
                                if attrib not in ['targdec', 'targra']:
                                    assert (pre_target.getAttribute(attrib) ==
                                            target.getAttribute(attrib))
                                else:
                                    if (pre_target.getAttribute(attrib) !=
                                        target.getAttribute(attrib)):
                                        ifu_spaxel = target.getAttribute('ifu_spaxel')
                                        targuse = target.getAttribute('targuse')
                                        logging.error(
                                            'unexpected change in attribe in target ' +
                                            '{} of field {} (ifu_spaxel={}, targuse={}): {}'.format(
                                                k + 1, i + 1, ifu_spaxel, targuse, attrib))
                            
                            match_found = True
                            break
                    
                    if match_found == False:
                        logging.error(
                            'Match not-found for target ' +
                            'in field {} (targname={}, targid={})'.format(
                                i + 1, targname, targid))
                
                else:
                
                    targname = pre_target.getAttribute('targname')
                    targid = pre_target.getAttribute('targid')
                    
                    match_found = False
                
                    for k, target in enumerate(field.getElementsByTagName('target')):
                        
                        if ((target.getAttribute('targname') == targname) and
                            (target.getAttribute('targid') == targid) and
                            ((target.getAttribute('ifu_spaxel')[-3:] == 'C04') or
                             (target.getAttribute('ifu_spaxel') == ''))):
                        
                            for attrib in pre_target.attributes.keys():
                            
                                if attrib not in ['cname', 'targdec', 'targra']:
                                    assert (pre_target.getAttribute(attrib) ==
                                            target.getAttribute(attrib))
                                elif attrib == 'cname':
                                    ifu_spaxel = target.getAttribute('ifu_spaxel')
                                    targuse = target.getAttribute('targuse')
                                    if target.getAttribute(attrib) != '%%%':
                                        logging.error(
                                            'unexpected cname in target ' +
                                            '{} of field {} (ifu_spaxel={}, targuse={})'.format(
                                                k + 1, i + 1, ifu_spaxel, targuse))
                            
                            match_found = True
                            break
                    
                    if match_found == False:
                        logging.error(
                            'Match not-found for target ' +
                            'in field {} (targname={}, targid={})'.format(
                                i + 1, targname, targid))
            
            else:
            
                raise ValueError


def check_configured_xml(filename, pre_xml_file=None):

    obsmode, survey_name, field_list, first_field_with_offset = (
        get_info_from_xml(filename))
        
    if pre_xml_file is not None:
        pre_obsmode, pre_survey_name, pre_field_list, pre_first_field_with_offset = (
            get_info_from_xml(pre_xml_file))
    
    check_attributes_of_targets(obsmode, field_list)
    
    check_attributes_of_photometry(field_list)
    
    check_values_of_attributes_of_photometry(field_list)

    if obsmode == 'LIFU':

        check_values_of_attributes_of_lifu_targets(field_list)
        
        if pre_xml_file is not None:
        
            check_lifu_versus_pre_xml_file(
                field_list, pre_field_list, first_field_with_offset)

    elif obsmode == 'mIFU':

        check_values_of_attributes_of_mifu_targets(field_list,
                                                   survey_name=survey_name)
        
        if pre_xml_file is not None:
        
            check_mifu_versus_pre_xml_file(field_list, pre_field_list)

    else:

        raise ValueError


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
             description='Check the contents of configured XML files')

    parser.add_argument('xml_file', help='a configured OB XML file')

    parser.add_argument('--pre', dest='pre_xml_file', default=None,
                        help='the pre-configured OB XML file')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    check_configured_xml(args.xml_file, pre_xml_file=args.pre_xml_file)

