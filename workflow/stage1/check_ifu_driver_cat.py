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
from workflow.utils.get_obstemp_info import (get_obstemp_info,
                                             get_obstemp_dict)
from workflow.utils.get_progtemp_info import (get_progtemp_dict,
                                              get_progtemp_info,
                                              get_obsmode_from_progtemp)
from workflow.utils.get_resources import (get_progtemp_file,
                                          get_obstemp_file)


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


def _check_targsrvy(array, trimester, cat_filename):

    icd_030_targsrvy_set = [
        'GA-LRDISC', 'GA-LRHIGHLAT', 'GA-HR', 'GA-OC', 'GA-CALIB', 'STEPS',
        'SCIP-AC', 'SCIP-CYG', 'SCIP-LR', 'WA', 'WC', 'WL-WIDE', 'WL-MID',
        'WL-DEEP', 'WQ', 'ASTRO-CALIB', 'WD', 'GS', 'ING-SYSCAT'
    ]

    sv_ot_regex = '^W[SV]([0-9]{4}[AB][12])-[0-9]{3}$'

    result = True

    targsrvy_set = set(array)

    for targsrvy_value in targsrvy_set:
        if targsrvy_value not in icd_030_targsrvy_set:

            if not re.match(sv_ot_regex, value):
                logging.error('unexpected TARGSRVY value: {}'.format(
                    targsrvy_value))
                result = False
            else:
                match = re.match(sv_ot_regex, targsrvy_value)

                trimester_in_sv_ot, = match.groups()

                if trimester_in_sv_ot != trimester:
                    logging.error(
                        'unexpected TARGSRVY value for TRIMESTER {}: {}'.format(
                            trimester, targsrvy_value))
                    result = False

    if len(targsrvy_set) > 1:
        logging.error(
            'there is more than one TARGSRV value in the catalogue')
        result = False

    if result == True:

        targsrvy = list(targsrvy_set)[0]
        
        if re.match(sv_ot_regex, targsrvy_value):

            expected_filename = '{}-ifu_driver_cat.fits'.format(targsrvy)

        else:

            expected_filename = '{}_{}-ifu_driver_cat.fits'.format(
                targsrvy, trimester)

        if expected_filename != os.path.basename(cat_filename):
            logging.error(
                'unexpected filename for TARGSRVY/TRIMESTE {}/{}: {}'.format(
                    targsrvy, trimester, expected_filename))
        

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


def _get_tl(hdu, colname, minmax):

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


def _check_in_tl_range(hdu, colname):

    result = True

    min_value = _get_tl(hdu, colname, 'min')
    max_value = _get_tl(hdu, colname, 'max')
    
    if min_value is None:
        logging.warning(
            'TLMIN value is not available for column {}'.format(colname))
    
    if max_value is None:
        logging.warning(
            'TLMAX value is not available for column {}'.format(colname))

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


def _check_targprio(hdu, progtemp_dict):

    result = True

    null_result = _check_non_null_value(hdu.data['TARGPRIO'])

    if null_result is False:
        logging.error('TARGPRIO values contain NaN')
        result = False

    range_result = _check_in_tl_range(hdu, 'TARGPRIO')

    if range_result is False:
        logging.error('TARGPRIO values out of range')
        result = False

    for targprio, progtemp in zip(hdu.data['TARGPRIO'], hdu.data['PROGTEMP']):

        obsmode = get_obsmode_from_progtemp(progtemp,
                                            progtemp_dict=progtemp_dict)

        if obsmode == 'LIFU':
            if targprio != 10.0:
                logging.error('TARGPRIO value must be 10.0 for LIFU mode')
                result = False

    return result


def _check_progtemp(array, progtemp_dict, forbidden_dict):

    result = True

    progtemp_set = set(array)

    for progtemp_value in progtemp_set:

        try:
            spectrograph_dict = get_progtemp_info(progtemp_value,
                                                  progtemp_dict=progtemp_dict)
        except:
            logging.error('unexpected PROGTEMP value: {}'.format(value))
            result = False

        for survey_type in forbidden_dict.keys():
            for regex in forbidden_dict[survey_type]:
                if re.match(regex, progtemp_value):
                    logging.warning(
                        'PROGTEMP value forbidden for {} surveys: {}'.format(
                            survey_type.upper(), progtemp_value))

    return result


def _check_obstemp(array, obstemp_dict):

    result = True

    obstemp_set = set(array)

    for obstemp_value in obstemp_set:

        try:
            obsconstraints_dict = get_obstemp_info(obstemp_value,
                                                   obstemp_dict=obstemp_dict)
        except:
            logging.error('unexpected OBSTEMP value: {}'.format(value))
            result = False

    return result


def _check_gaia_ra(hdu):

    result = True

    null_result = _check_non_null_value(hdu.data['GAIA_RA'])

    if null_result is False:
        logging.error('GAIA_RA values contain NaN')
        result = False

    range_result = _check_in_tl_range(hdu, 'GAIA_RA')

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

    range_result = _check_in_tl_range(hdu, 'GAIA_DEC')

    if range_result is False:
        logging.error('GAIA_DEC values out of range')
        result = False

    return result


def _check_gaia_epoch(hdu):

    result = True

    null_result = _check_non_null_value(hdu.data['GAIA_EPOCH'])

    if null_result is False:
        logging.error('GAIA_EPOCH values contain NaN')
        result = False

    range_result = _check_in_tl_range(hdu, 'GAIA_EPOCH')

    if range_result is False:
        logging.error('GAIA_EPOCH values out of range')
        result = False

    return result


def _check_gaia_pm(hdu):

    result = True

    pmra_null_result = _check_non_null_value(hdu.data['GAIA_PMRA'])

    if pmra_null_result is False:
        logging.error('GAIA_PMRA values contain NaN')
        result = False

    pmdec_null_result = _check_non_null_value(hdu.data['GAIA_PMDEC'])

    if pmdec_null_result is False:
        logging.error('GAIA_PMDEC values contain NaN')
        result = False
    
    return result


def _check_gaia_paral(hdu):

    result = True

    null_result = _check_non_null_value(hdu.data['GAIA_PARAL'])

    if null_result is False:
        logging.error('GAIA_PARAL values contain NaN')
        result = False
    
    return result


def _check_ifu_pa_request(hdu, progtemp_dict):

    result = True

    range_result = _check_in_tl_range(hdu, 'IFU_PA_REQUEST')

    if range_result is False:
        logging.error('IFU_PA_REQUEST values out of range')
        result = False

    ifu_pa_request_array = hdu.data['IFU_PA_REQUEST']
    progtemp_array = hdu.data['PROGTEMP']

    for i, (ifu_pa_request, progtemp) in enumerate(
            zip(ifu_pa_request_array, progtemp_array)):

        if (not np.isnan(ifu_pa_request)):

            obsmode = get_obsmode_from_progtemp(progtemp, progtemp_dict)

            if obsmode != 'LIFU':

                logging.error(
                    'unexpected non-null IFU_PA_REQUEST for non-LIFU PROGTEMP' +
                    ' in row {}'.format(i + 1))
                result = False

    return result


def _check_ifu_dither(hdu, progtemp_dict):

    result = True

    range_result = _check_in_tl_range(hdu, 'IFU_DITHER')

    if range_result is False:
        logging.error('IFU_DITHER values out of range')
        result = False

    ifu_dither_array = hdu.data['IFU_DITHER']
    progtemp_array = hdu.data['PROGTEMP']

    for i, (ifu_dither, progtemp) in enumerate(
            zip(ifu_dither_array, progtemp_array)):


        obsmode = get_obsmode_from_progtemp(progtemp, progtemp_dict)

        if (obsmode == 'LIFU'):
            if (ifu_dither not in [-1, 0, -3, 3, 4, 5, 6]):
                logging.error(
                    'unexpected IFU_DITHER for LIFU PROGTEMP in row {}'.format(
                        i + 1))
                result = False
        elif (obsmode == 'mIFU'):
            if (ifu_dither not in [0, -3, 3, 4, 5, 6]):
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


def _check_progtemp_ifu_dither_compatibility(hdu, progtemp_dict):

    result = True

    progtemp_array= hdu.data['PROGTEMP']
    ifu_dither_array = hdu.data['IFU_DITHER']

    # Check the non-custom dithers

    for i, (progtemp, ifu_dither) in enumerate(
            zip(progtemp_array, ifu_dither_array)):

        if np.abs(ifu_dither) != 1:

            num_dither = np.abs(ifu_dither)

            spectrograph_dict = get_progtemp_info(progtemp,
                                                  progtemp_dict=progtemp_dict)

            num_exp_r = spectrograph_dict['red_num_exposures']
            num_exp_b = spectrograph_dict['blue_num_exposures']
 
            if ((num_exp_r % num_dither != 0) or 
                (num_exp_b % num_dither != 0)):
                logging.error(
                    'PROGTEMP and IFU_DITHER are not compatible in row {}'.format(
                        i + 1))
                result = False
 
            if not ((num_exp_r % num_exp_b == 0) or
                    (num_exp_b % num_exp_r == 0)):
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

        spectrograph_dict = get_progtemp_info(progtemp,
                                              progtemp_dict=progtemp_dict)

        num_exp_r = spectrograph_dict['red_num_exposures']
        num_exp_b = spectrograph_dict['blue_num_exposures']

        if ((num_exp_r % num_dither != 0) or 
            (num_exp_b % num_dither != 0)):
            logging.error(
                'PROGTEMP and custom dither {}:{}:{}:{} are not compatible'.format(
                    targname, targid, progtemp, obstemp)
            )
            result = False
 
        if not ((num_exp_r % num_exp_b == 0) or
                (num_exp_b % num_exp_r == 0)):
            logging.error(
                'PROGTEMP and custom dither {}:{}:{}:{} are not compatible'.format(
                    targname, targid, progtemp, obstemp)
            )
            result = False

    return result


def check_ifu_driver_cat(cat_filename, template=None, check_vs_template=True,
                         progtemp_file=None, obstemp_file=None):
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
    progtemp_file : str, optional
        A progtemp.dat file with the definition of PROGTEMP.
    obstemp_file : str, optional
        A obstemp.dat file with the definition of OBSTEMP.

    Returns
    -------
    result : bool
        True if the file passes all the checks, otherwise False.
    """
    
    result = True

    # Get the dictionaries to interpret PROGTEMP and OBSTEMP, and their
    # DATAMVER values

    progtemp_datamver, progtemp_dict, forbidden_dict = get_progtemp_dict(
        filename=progtemp_file, assert_orb=True)

    obstemp_datamver, obstemp_dict = get_obstemp_dict(
        filename=obstemp_file)

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

        # Check DATAMVER of the IFU driver cat, PROGTEMP file and OBSTEMP file
        # are consistent

        datamver = hdu_list[0].header['DATAMVER']

        if datamver != progtemp_datamver:
            logging.critical(
                'DATAMVER mismatch ({} != {}) for PROGTEMP file: '.format(
                    datamver, progtemp_datamver) +
                'Stop unless you are sure!')
            raise SystemExit(2)

        if datamver != obstemp_datamver:
            logging.critical(
                'DATAMVER mismatch ({} != {}) for OBSTEMP file: '.format(
                    datamver, obstemp_datamver) +
                'Stop unless you are sure!')
            raise SystemExit(2)

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

        if not _check_targsrvy(data['TARGSRVY'], trimester, cat_filename):
            result = False

        # Nothing to check in TARGPROG

        if not _check_targid(data['TARGID']):
            result = False

        if not _check_targname(data['TARGNAME']):
            result = False

        if not _check_targprio(hdu1, progtemp_dict):
            result = False

        if not _check_progtemp(data['PROGTEMP'], progtemp_dict, forbidden_dict):
            result = False

        if not _check_obstemp(data['OBSTEMP'], obstemp_dict):
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

        if not _check_ifu_pa_request(hdu1, progtemp_dict):
            result = False

        if not _check_ifu_dither(hdu1, progtemp_dict):
            result = False

        if not _check_repeated_objects(hdu1):
            result = False

        if not _check_progtemp_ifu_dither_compatibility(hdu1, progtemp_dict):
            result = False

    # Write a message

    if result == True:
        logging.info('{} has passed the checks'.format(cat_filename))
    else:
        logging.error('{} has NOT passed the checks'.format(cat_filename))

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

    parser.add_argument('--progtemp_file',
                        default=os.path.join('aux', 'progtemp.dat'),
                        help="""a progtemp.dat file with the definition of
                        PROGTEMP""")


    parser.add_argument('--obstemp_file',
                        default=os.path.join('aux', 'obstemp.dat'),
                        help="""a obstemp.dat file with the definition of
                        OBSTEMP""")
    
    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')
    
    args = parser.parse_args()
    
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.progtemp_file):
        logging.info('Downloading the progtemp file')
        get_progtemp_file(file_path=args.progtemp_file)

    if not os.path.exists(args.obstemp_file):
        logging.info('Downloading the obstemp file')
        get_obstemp_file(file_path=args.obstemp_file)
    
    check_ifu_driver_cat(args.catalogue, template=args.template,
                         check_vs_template=args.check_vs_template,
                         progtemp_file=args.progtemp_file,
                         obstemp_file=args.obstemp_file)

