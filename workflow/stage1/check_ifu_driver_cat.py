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
import re

from astropy.io import fits

from workflow.utils import check_equal_headers


def _check_versus_template(cat_filename, template):

    ignore_values = ['TRIMESTE', 'DATETIME', 'VERBOSE', 'AUTHOR', 'CCREPORT']

    result = check_equal_headers(cat_filename, template,
                                 ignore_values=ignore_values)

    return result


def _check_trimester(trimester):

    if re.match('[0-9]{4}[AB][12]', trimester):
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


def _check_targsrvy(cat_filename):

    # Avoid non empty

    raise NotImplementedError


def _check_targid(cat_filename):

    # Avoid non empty

    raise NotImplementedError


def _check_targname(cat_filename):

    # Avoid non empty

    raise NotImplementedError


def _check_targprio(cat_filename):

    # Avoid null
    # Check in range

    raise NotImplementedError


def _check_progtemp(cat_filename):

    raise NotImplementedError


def _check_obstemp(cat_filename):

    raise NotImplementedError


def _check_gaia_ra(cat_filename):

    # Avoid null
    # Check in range

    raise NotImplementedError


def _check_gaia_dec(cat_filename):

    # Avoid null
    # Check in range

    raise NotImplementedError


def _check_gaia_epoch(cat_filename):

    # Check in range

    raise NotImplementedError


def _check_gaia_paral(cat_filename):

    # Check in range

    raise NotImplementedError


def _check_gaia_ifu_request(cat_filename):

    # Check in range

    raise NotImplementedError


def _check_gaia_ifu_dither(cat_filename):

    # Check specific values

    raise NotImplementedError


def _check_consistency_progtemp_dither(cat_filename):

    raise NotImplementedError


def _check_dither_size(cat_filename):

    raise NotImplementedError


def _check_custom_dither(cat_filename):

    # Check the size depending on the mode

    # Check that some values are the same for different rows
    # (e.g. ifu_pa_request)

    raise NotImplementedError


def _check_locations_for_names(cat_filename):

    raise NotImplementedError


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

    with fits.open(cat_filename) as hdu_list:

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

    return result

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description='Check the contents of an IFU driver catalogue')
    
    parser.add_argument('catalogue',
                        help='a FITS file with an IFU driver catalogue')
    
    parser.add_argument('--template', default='./aux/ifu_driver_template.fits',
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

