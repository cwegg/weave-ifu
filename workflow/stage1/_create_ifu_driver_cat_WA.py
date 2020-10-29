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


import os

import numpy as np

from workflow.stage1 import create_ifu_driver_cat


def set_keywords_info():

    trimester = '2020A1'
    author = 'a@domain.com'
    report_verbosity = 1
    cc_report = 'b@domain.com,c@domain.com'

    return trimester, author, report_verbosity, cc_report


def get_data_dict():

    data_dict = {}

    data_dict['TARGSRVY'] = \
        ['WA', 'WA', 'WA', 'WA']

    data_dict['TARGPROG'] = \
        ['WA_HR', 'WA_LR', 'WA_LR', 'WA_BW_LR']

    data_dict['TARGID'] = \
        ['HSB_HR', 'HSB_LR', 'LSB_LR', 'BWT_LR']

    data_dict['TARGNAME'] = \
        ['PSZ1_nameHSB', 'PSZ1_nameHSB', 'PSZ1_nameLSB', 'PSZ1_nameBW']

    data_dict['TARGPRIO'] = \
        [10.0, 10.0, 10.0, 10.0, 1.0]

    data_dict['PROGTEMP'] = \
        ['61331', '41331', '41331', '41331']

    data_dict['OBSTEMP'] = \
        ['NBCEB', 'NBCEC', 'NBCEC', 'XDCED']

    data_dict['GAIA_RA'] = \
        [48.14007743, 48.14007743, 183.84864, 184.88371]

    data_dict['GAIA_DEC'] = \
        [39.31920434, 39.31920434, 51.349714, 49.815753]

    data_dict['GAIA_EPOCH'] = \
        [2015.5, 2015.5, 2015.5, 2015.5]

    data_dict['GAIA_PMRA'] = \
        [np.nan, np.nan, np.nan, np.nan]

    data_dict['GAIA_PMDEC'] = \
        [np.nan, np.nan, np.nan, np.nan]

    data_dict['GAIA_PARAL'] = \
        [np.nan, np.nan, np.nan, np.nan]

    data_dict['IFU_PA_REQUEST'] = \
        [0.0, 0.0, 0.0, 0.0]

    data_dict['IFU_DITHER'] = \
        [3, 3, 3, 3]

    return data_dict


if __name__ == '__main__':

    ############################################################################
    # Set the location of the template and the output file and directory

    ifu_driver_template = os.path.join('aux', 'ifu_driver_template.fits')

    output_dir = 'output'
    output_filename = os.path.join(output_dir, 'WA_2020A1-ifu_driver_cat.fits')

    ############################################################################
    # Create the output directory if it does not exist

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ############################################################################
    # Get a dictionary with the data which will populate the template

    # NOTE: See the above function get_data_dict to understand the structure of
    # the dictionary

    data_dict = get_data_dict()

    ############################################################################
    # Set the needed information to populate some keywords of the primary header

    trimester, author, report_verbosity, cc_report = set_keywords_info()

    ############################################################################
    # Create the IFU driver catalogue

    create_ifu_driver_cat(ifu_driver_template, data_dict, output_filename,
                          trimester, author, report_verbosity=report_verbosity,
                          cc_report=cc_report)

