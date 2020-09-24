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

from astropy.io import fits

from workflow.utils import check_equal_headers
from workflow.utils import populate_fits_table_template


def _get_data_dict(fits_filename):
    
    data_dict = {}
    
    with fits.open(fits_filename) as hdu_list:
        data = hdu_list[1].data
        
        for col in data.names:
            data_dict[col] = list(data[col])
    
    return data_dict


def _get_combo_data_dict(fits_filename1, fits_filename2):
    
    data_dict1 = _get_data_dict(fits_filename1)
    data_dict2 = _get_data_dict(fits_filename2)
    
    combo_data_dict = {}
    
    for key in data_dict1.keys():
        combo_data_dict[key] = data_dict1[key] + data_dict2[key]
    
    return combo_data_dict


def create_combo_fits_cat(mos_cat, ifu_cat, output_filename, overwrite=False):
    """
    Combine MOS and IFU catalogues to create a combo catalogue.
    
    Parameters
    ----------
    mos_cat : str
        Name of a FITS file containing a MOS catalogue.
    ifu_cat : list of str
        Name of a FITS file containing a IFU catalogue.
    output_filename : str
        The name of the output file which will be created.
    overwrite : bool, optional
        Overwrite the output FITS file.
    """
    
    equal_headers = check_equal_headers(mos_cat, ifu_cat,
                                        ignore_values=['DATETIME'])
    
    if not equal_headers:
        logging.error('equal headers are expected to create a combo FITS file')
        
        assert equal_headers
    
    combo_data_dict = _get_combo_data_dict(mos_cat, ifu_cat)
    
    populate_fits_table_template(mos_cat, combo_data_dict, output_filename,
                                 update_datetime=True, overwrite=overwrite)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create a combo catalogue')

    parser.add_argument('mos_cat',
                        help="""a FITS file containing a MOS catalogue""")

    parser.add_argument('ifu_cat',
                        help="""a FITS file containing a IFU catalogue""")

    parser.add_argument('--out', dest='output_filename',
                        default='output/combo_cat.fits',
                        help="""name for the output file which will contain the
                        combo catalogue""")

    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='overwrite the output file')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(os.path.dirname(args.output_filename)):
        logging.info('Creating the output directory')
        os.mkdir(os.path.dirname(args.output_filename))

    create_combo_fits_cat(args.mos_cat, args.ifu_cat, args.output_filename,
                          overwrite=args.overwrite)

