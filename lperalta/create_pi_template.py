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
import os.path

from workflow.utils import create_sub_template
from workflow.utils.get_resources import get_master_cat
        
    
def create_pi_template(catalogue_template, output_filename, col_list,
                       extname='CATALOGUE', update_datetime=True,
                       overwrite=False):
    """
    Create a template for the input IFU driver catalogues.

    Parameters
    ----------
    catalogue_template : str
        The name of the master catalogue template.
    output_filename : str
        The name of the output file for the new PI template.
    col_list : list of str
        A list containing all the columns which will be included in the new PI
        template.
    extname : str, optional
        Name to be used in the table extension.
    update_datetime : bool, optional
        Update DATETIME keyword in the new PI template.
    overwrite : bool, optional
        Overwrite the output FITS file containing the PI template.
    """
    
    create_sub_template(catalogue_template, output_filename, col_list,
                        extname=extname, inherit_primary_kwds=True,
                        update_datetime=update_datetime, overwrite=overwrite)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Create a template for a PI proposal')

    parser.add_argument('--in', dest='catalogue_template',
                        default='Master_CatalogueTemplate.fits',
                        help='name of the master catalogue template')

    parser.add_argument('--out', dest='pi_template',
                        default='pi_template.fits',
                        help='name for the output PI template')

    parser.add_argument('--cols', dest='cols_str',
                        default=('TARGSRVY,TARGPROG,TARGID,TARGNAME,TARGPRIO,' +
                                 'PROGTEMP,OBSTEMP,GAIA_RA,GAIA_DEC,' +
                                 'GAIA_EPOCH,GAIA_PMRA,GAIA_PMDEC,GAIA_PARAL,' +
                                 'IFU_PA,IFU_DITHER'),
                        help='name of the columns to be included in the table')

    parser.add_argument('--extname', dest='extname',
                        default='CATALOGUE',
                        help='name to be used in the table extension')

    parser.add_argument('--keep_datetime', dest='update_datetime',
                        action='store_false',
                        help='keep DATETIME keyword from the master template')

    parser.add_argument('--overwrite', dest='overwrite',
                        action='store_true',
                        help='overwrite the output file')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='The level for the logging messages')

    args = parser.parse_args()
    
    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}
    
    logging.basicConfig(level=level_dict[args.log_level])

    if not os.path.exists(args.catalogue_template):
        logging.info('Downloading the master catalogue template')
        get_master_cat(file_path=args.catalogue_template)

    col_list = args.cols_str.split(',')

    create_pi_template(args.catalogue_template, args.pi_template, col_list,
                       extname=args.extname,
                       update_datetime=args.update_datetime,
                       overwrite=args.overwrite)

