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
import os
import re
from collections import OrderedDict

import numpy as np
from astropy.io import fits

from workflow.utils import check_equal_headers
from workflow.utils.get_progtemp_info import get_obsmode_from_progtemp


def _check_versus_template(cat_filename, template):

    ignore_values = ['TRIMESTE', 'DATETIME', 'VERBOSE', 'AUTHOR', 'CCREPORT']

    result = check_equal_headers(cat_filename, template,
                                 ignore_values=ignore_values)

    return result


def _check_trimester(trimester):

    if re.match('^[0-9]{4}[AB][12]$', trimester):
        result = True
    else:
        result = False

    return result


def _check_verbose(verbose):

    result = (verbose in [0, 1])

    return result


def _check_email(string):

    arroba_index = string.find('@')

    if arroba_index != -1:

        substring = string[arroba_index + 1:]
        dot_index = substring.find('.')

        if dot_index != -1:
            result = True
        else:
            result = False

    else:
        result = False

    return result


def _check_author(author):

    result = _check_email(author)

    return result


def _check_ccreport(ccreport):

    result = True

    for string in ccreport.split(','):
        if not _check_email(string):
           result = False

    return result


def _check_targsrvy(array):

    icd_030_targsrvy_set = [
        'GA-LRDISC', 'GA-LRHIGHLAT', 'GA-HR', 'GA-OC', 'GA-CALIB', 'STEPS',
        'SCIP-AC', 'SCIP-CYG', 'SCIP-LR', 'WA', 'WC', 'WL-WIDE', 'WL-MID',
        'WL-DEEP', 'WQ', 'ASTRO-CALIB', 'WD', 'GS', 'ING-SYSCAT'
    ]

    result = True

    targsrvy_set = set(array)

    for value in targsrvy_set:
        if value not in icd_030_targsrvy_set:
            logging.error('unexpected TARGSRVY value: {}'.format(value))
            result = False

    if len(targsrvy_set) > 1:
        logging.warning(
            'there is more than one TARGSRV value in the catalogue, ' +
            'are you sure?'
        )

    return result


def _check_non_empty_string(array):

    result = True

    set_of_values = set(array)

    if '' in set_of_values:
        result = False

    return result


def _check_targid(array):

    result = _check_non_empty_string(array)

    if result is False:
        logging.error('TARGID value should not be populated by empty strings')

    return result


def _check_targname(array):

    result = _check_non_empty_string(array)

    if result is False:
        logging.error('TARGNAME value should not be populated by empty strings')

    return result


def _check_non_null_value(array):

    result = True

    for value in array:
        if np.isnan(value):
            result = False
            break

    return result


def _get_td(hdu, colname, minmax):

    assert minmax in ['min', 'max']

    if colname in hdu.data.columns.names:

        index = hdu.data.columns.names.index(colname)

        kwd = 'TL{}{}'.format(minmax.upper(), index + 1)

        if kwd in hdu.header:
            value = hdu.header[kwd]
        else:
            value = None

    else:
        value = None

    return value


def _check_in_td_range(hdu, colname):

    result = True

    min_value = _get_td(hdu, colname, 'min')
    max_value = _get_td(hdu, colname, 'max')
    
    if min_value is None:
        logging.warning(
            'TDMIN value is not available for column {}'.format(colname))
    
    if max_value is None:
        logging.warning(
            'TDMAX value is not available for column {}'.format(colname))

    for value in hdu.data[colname]:

        if not np.isnan(value):

            if min_value is not None:
                if not (value >= min_value):
                    result = False
                    break

            if max_value is not None:
                if not (value <= max_value):
                    result = False
                    break

    return result


def _check_targprio(hdu):

    result = True

    null_result = _check_non_null_value(hdu.data['TARGPRIO'])

    if null_result is False:
        logging.error('TARGPRIO values contain NaN')
        result = False

    range_result = _check_in_td_range(hdu, 'TARGPRIO')

    if range_result is False:
        logging.error('TARGPRIO values out of range')
        result = False

    return result


def _check_progtemp(array):

    progtemp_regex = '^[4-9][01239][0-9][0-9][12349](\.[0-9]+\+?)?$'

    result = True

    progtemp_set = set(array)

    for value in progtemp_set:

        value_result = True

        if not re.match(progtemp_regex, value):
            value_result = False

        try:
            o_char = value[1]

            r_char = value[2]
            b_char = value[3]

            if o_char == '0':
                assert r_char not in '12357'
                assert b_char not in '12357'
            elif o_char == '1':
                assert r_char not in '1'
                assert b_char not in '1'
            elif o_char == '2':
                assert r_char not in '15'
                assert b_char not in '15'
        except:
            value_result = False

        if value_result is False:
            logging.error('unexpected PROGTEMP value: {}'.format(value))
            result = False

    return result


def _check_obstemp(array):

    obstemp_regex = '^[A-X][A-E][A-F][A-E][A-F]$'

    result = True

    obstemp_set = set(array)

    for value in obstemp_set:
        if not re.match(obstemp_regex, value):
            logging.error('unexpected OBSTEMP value: {}'.format(value))
            result = False

    return result


def _check_gaia_ra(hdu):

    result = True

    null_result = _check_non_null_value(hdu.data['GAIA_RA'])

    if null_result is False:
        logging.error('GAIA_RA values contain NaN')
        result = False

    range_result = _check_in_td_range(hdu, 'GAIA_RA')

    if range_result is False:
        logging.error('GAIA_RA values out of range')
        result = False

    return result


def _check_gaia_dec(hdu):

    result = True

    null_result = _check_non_null_value(hdu.data['GAIA_DEC'])

    if null_result is False:
        logging.error('GAIA_DEC values contain NaN')
        result = False

    range_result = _check_in_td_range(hdu, 'GAIA_DEC')

    if range_result is False:
        logging.error('GAIA_DEC values out of range')
        result = False

    return result


def _check_gaia_epoch(hdu):

    result = _check_in_td_range(hdu, 'GAIA_EPOCH')

    if result is False:
        logging.error('GAIA_EPOCH values out of range')

    return result


def _check_gaia_pm(hdu):

    result = True

    pmra_array = hdu.data['GAIA_PMRA']
    pmdec_array = hdu.data['GAIA_PMDEC']
    epoch_array = hdu.data['GAIA_EPOCH']

    for i, (pmra, pmdec, epoch) in enumerate(
            zip(pmra_array, pmdec_array, epoch_array)):

        if (not np.isnan(pmra)) or (not np.isnan(pmdec)):

            if np.isnan(pmra) or np.isnan(pmdec) or np.isnan(epoch):
                logging.error(
                    'unexpected GAIA_PMRA/GAIA_PMDEC/GAIA_EPOCH values in ' +
                    'row {}'.format(i + 1))
                result = False
    
    return result


def _check_gaia_paral(hdu):

    result = True

    paral_array = hdu.data['GAIA_PARAL']
    epoch_array = hdu.data['GAIA_EPOCH']

    for i, (paral, epoch) in enumerate(zip(paral_array, epoch_array)):

        if not np.isnan(paral):

            if np.isnan(epoch):
                logging.error(
                    'unexpected GAIA_PARAL/GAIA_EPOCH values in row {}'.format(
                        i + 1))
                result = False
    
    return result


def _check_ifu_pa_request(hdu):

    result = True

    range_result = _check_in_td_range(hdu, 'IFU_PA_REQUEST')

    if range_result is False:
        logging.error('IFU_PA_REQUEST values out of range')
        result = False

    ifu_pa_request_array = hdu.data['IFU_PA_REQUEST']
    progtemp_array = hdu.data['PROGTEMP']

    for i, (ifu_pa_request, progtemp) in enumerate(
            zip(ifu_pa_request_array, progtemp_array)):

        if (not np.isnan(ifu_pa_request)):

            obsmode = get_obsmode_from_progtemp(progtemp)

            if obsmode != 'LIFU':

                logging.error(
                    'unexpected non-null IFU_PA_REQUEST for non-LIFU PROGTEMP' +
                    ' in row {}'.format(i + 1))
                result = False

    return result


def _check_ifu_dither(hdu):

    result = True

    range_result = _check_in_td_range(hdu, 'IFU_DITHER')

    if range_result is False:
        logging.error('IFU_DITHER values out of range')
        result = False

    ifu_dither_array = hdu.data['IFU_DITHER']
    progtemp_array = hdu.data['PROGTEMP']

    for i, (ifu_dither, progtemp) in enumerate(
            zip(ifu_dither_array, progtemp_array)):


        obsmode = get_obsmode_from_progtemp(progtemp)

        if (obsmode == 'LIFU'):
            if (ifu_dither not in [-1, 0, 3, 5]):
                logging.error(
                    'unexpected IFU_DITHER for LIFU PROGTEMP in row {}'.format(
                        i + 1))
                result = False
        elif (obsmode == 'mIFU'):
            if (ifu_dither not in [0, 3, 5]):
                logging.error(
                    'unexpected IFU_DITHER for mIFU PROGTEMP in row {}'.format(
                        i + 1))
                result = False
        elif (ifu_dither != 0):
            logging.error(
                'unexpected IFU_DITHER and PROGTEMP in row {}'.format(i + 1))
            result = False

    return result

def _nan_to_none(value):

    if np.isnan(value):
        result = None
    else:
        result = value

    return result


def _check_repeated_objects(hdu, offset_tol_arcsec=10.0):

    result = True

    # Create a dictionary with TARGNAME/TARGID as keys and empty lists as values

    targname_array = hdu.data['TARGNAME']
    targid_array = hdu.data['TARGID']

    object_dict = {(targname, targid): []
                   for targname, targid in zip(targname_array, targid_array)}

    # Populate the lists with the information of the dithering pattern and
    # coordinates
    
    for i in range(len(hdu.data)):

        targname = hdu.data['TARGNAME'][i]
        targid = hdu.data['TARGID'][i]

        ifu_dither = hdu.data['IFU_DITHER'][i]

        gaia_ra = hdu.data['GAIA_RA'][i]
        gaia_dec = hdu.data['GAIA_DEC'][i]
        gaia_epoch = hdu.data['GAIA_EPOCH'][i]
        gaia_pmra = hdu.data['GAIA_PMRA'][i]
        gaia_pmdec = hdu.data['GAIA_PMDEC'][i]
        gaia_paral = hdu.data['GAIA_PARAL'][i]

        object_dict[(targname, targid)].append(
            (ifu_dither,
             (gaia_ra, gaia_dec, _nan_to_none(gaia_epoch),
              _nan_to_none(gaia_pmra), _nan_to_none(gaia_pmdec),
              _nan_to_none(gaia_paral)))
        )

    # Check the uniqueness of the coordinates

    for key in object_dict.keys():

        targname, targid = key

        list_in_value = object_dict[key]

        ifu_dither_list = [elem[0] for elem in list_in_value]
        coord_list = [elem[1] for elem in list_in_value]

        # Case 1: We do not have custom dithers
        if (-1 not in ifu_dither_list):
            if len(set(coord_list)) != 1:
                logging.error(
                    'unexpected change of coordinates in non-custom dither {}:{}'.format(
                        targname, targid))
                result = False
        else:

            # If we have a mix of non-custom and custom dithers, take the non-
            # custom as reference. If we have only custom dithers, take the
            # first position as reference
            
            non_custom_coord_list = [elem[1] for elem in list_in_value
                                     if elem[0] != -1]

            non_custom_coord_set = set(non_custom_coord_list)

            # All the non-custom dithers should have the same coordinates
            
            if len(non_custom_coord_set) > 1:
                logging.error(
                    'unexpected change of coordinates for {}:{}'.format(
                        targname, targid))
                result = False
                continue
            elif len(non_custom_coord_set) == 1:
                ref_coord = non_custom_coord_list[0]
            else:
                ref_coord = coord_list[0]

            # Check that the offsets are not far from the reference
            
            for ifu_dither, coord in list_in_value:
                
                if ifu_dither == -1:

                    ra_offset = np.abs((coord[0] - ref_coord[0]) *
                                       np.cos(np.deg2rad(ref_coord[1]))) * 3600
                    dec_offset = np.abs(coord[1] - ref_coord[1]) * 3600

                    if ((ra_offset > offset_tol_arcsec) or
                        (dec_offset > offset_tol_arcsec)):
                        logging.error(
                            'too big offset for custom dither in {}:{}'.format(
                                targname, targid))
                        result = False
                        continue

                    for i in range(2, 6):
                        if ((coord[i] != ref_coord[i]) and not
                            (np.isnan(coord[i]) and np.isnan(ref_coord[i]))):
                            logging.error(
                                'unexpected change of GAIA_EPOCH/GAIA_PMRA/'
                                'GAIA_PMDEC/GAIA_PARAL for {}:{}'.format(
                                    targname, targid))
                            result = False
                            continue

    return result


def _progtemp_to_num_exp(progtemp):

    o_char = progtemp[1]
    r_char = progtemp[2]
    b_char = progtemp[3]

    convert_dict = {
        '0': {'0': 1, '1': None, '2': None, '3': None, '4': 2,
              '5': None, '6': 3, '7': None, '8': 4, '9': 5},
        '1': {'0': 1, '1': None, '2': 2, '3': 3, '4': 4,
              '5': 5, '6': 6, '7': 7, '8': 8, '9': 9},
        '2': {'0': 1, '1': None, '2': 3, '3': 5, '4': 6,
              '5': None, '6': 9, '7': 10, '8': 12, '9': 15},
        '3': {'0': 1, '1': 2, '2': 4, '3': 6, '4': 8,
              '5': 10, '6': 12, '7': 14, '8': 16, '9': 20}
    }

    if o_char in convert_dict.keys():
        num_exp_r = convert_dict[o_char][r_char]
        num_exp_b = convert_dict[o_char][b_char]
    else:
        num_exp_r = None
        num_exp_b = None

    result = (num_exp_r, num_exp_b)

    return result


def _check_progtemp_ifu_dither_compatibility(hdu):

    result = True

    progtemp_array= hdu.data['PROGTEMP']
    ifu_dither_array = hdu.data['IFU_DITHER']

    # Check the non-custom dithers

    for i, (progtemp, ifu_dither) in enumerate(
            zip(progtemp_array, ifu_dither_array)):

        if ifu_dither > 1:

            num_dither = ifu_dither
            num_exp_r, num_exp_b = _progtemp_to_num_exp(progtemp)
 
            if ((num_exp_r % num_dither != 0) or 
                (num_exp_b % num_dither != 0)):
                logging.error(
                    'PROGTEMP and IFU_DITHER are not compatible in row {}'.format(
                        i + 1))
                result = False

    # Check the custom dithers

    targname_array = hdu.data['TARGNAME']
    targid_array = hdu.data['TARGID']
    obstemp_array = hdu.data['OBSTEMP']

    custom_dict = OrderedDict()

    for i, (ifu_dither, targname, targid, progtemp, obstemp) in enumerate(
            zip(ifu_dither_array,
                targname_array, targid_array, progtemp_array, obstemp_array)):

        if ifu_dither == -1:

            key = (targname, targid, progtemp, obstemp)

            if key in custom_dict:
                custom_dict[key] += 1
            else:
                custom_dict[key] = 1

    for key in custom_dict.keys():

        targname, targid, progtemp, obstemp = key

        num_dither = custom_dict[key]
        num_exp_r, num_exp_b = _progtemp_to_num_exp(progtemp)

        if ((num_exp_r % num_dither != 0) or 
            (num_exp_b % num_dither != 0)):
            logging.error(
                'PROGTEMP and custom dither {}:{}:{}:{} are not compatible'.format(
                    targname, targid, progtemp, obstemp)
            )
            result = False

    return result


def check_ifu_driver_cat(cat_filename, template=None, check_vs_template=True):
    """
    Check the contents of an IFU driver catalogue.

    Parameters
    ----------
    cat_filename : str
        A FITS file with an IFU driver catalogue.
    template : str, optional
        A FITS file containing an IFU driver template.
    check_vs_template: bool, optional
        An option to indicate whether the file should be checked against a
        template or not.

    Returns
    ----------
    result : bool
        True if the file passes all the checks, otherwise False.
    """
    
    result = True

    # Check the file versus a template (if requested)
    
    if check_vs_template is True:

        if template is not None:

            if not _check_versus_template(cat_filename, template):
                logging.error(
                    '{} does not match the template {}'.format(
                        cat_filename, template))
                result = False

        else:
            logging.error(
                'a template for making the checks has not been provided')
            result = False

    # Open the catalogue file

    with fits.open(cat_filename) as hdu_list:

        # Check the values of its primary header

        trimester = hdu_list[0].header['TRIMESTE']
        if not _check_trimester(trimester):
            logging.error(
                'invalid value in the TRIMESTE keyword: {}'.format(trimester))
            result = False

        verbose = hdu_list[0].header['VERBOSE']
        if not _check_verbose(verbose):
            logging.error(
                'invalid value in the VERBOSE keyword: {}'.format(verbose))
            result = False

        author = hdu_list[0].header['AUTHOR']
        if not _check_author(author):
            logging.error(
                'invalid email address in the AUTHOR keyword: {}'.format(
                    author))
            result = False

        ccreport = hdu_list[0].header['CCREPORT']
        if not _check_ccreport(ccreport):
            logging.error(
                'invalid email addresses in the CCREPORT keyword: {}'.format(
                    ccreport))
            result = False

        # Check the data in the table of its first extension

        hdu1 = hdu_list[1]
        data = hdu_list[1].data

        if not _check_targsrvy(data['TARGSRVY']):
            result = False

        # Nothing to check in TARGPROG

        if not _check_targid(data['TARGID']):
            result = False

        if not _check_targname(data['TARGNAME']):
            result = False

        if not _check_targprio(hdu1):
            result = False

        if not _check_progtemp(data['PROGTEMP']):
            result = False

        if not _check_obstemp(data['OBSTEMP']):
            result = False

        if not _check_gaia_ra(hdu1):
            result = False

        if not _check_gaia_dec(hdu1):
            result = False

        if not _check_gaia_epoch(hdu1):
            result = False

        if not _check_gaia_pm(hdu1):
            result = False

        if not _check_gaia_paral(hdu1):
            result = False

        if not _check_ifu_pa_request(hdu1):
            result = False

        if not _check_ifu_dither(hdu1):
            result = False

        if not _check_repeated_objects(hdu1):
            result = False

        if not _check_progtemp_ifu_dither_compatibility(hdu1):
            result = False

    return result

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description='Check the contents of an IFU driver catalogue')
    
    parser.add_argument('catalogue',
                        help='a FITS file with an IFU driver catalogue')
    
    parser.add_argument('--template',
                        default=os.path.join('aux', 'ifu_driver_template.fits'),
                        help='a FITS file containing an IFU driver template')

    parser.add_argument('--no_check_vs_template', dest='check_vs_template',
                        action='store_false',
                        help=
                        'skip the check of the catalogue versus the template')
    
    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')
    
    args = parser.parse_args()
    
    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}
    
    logging.basicConfig(level=level_dict[args.log_level])
    
    check_ifu_driver_cat(args.catalogue, template=args.template,
                         check_vs_template=args.check_vs_template)

