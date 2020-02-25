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
import collections
import logging
import os.path

from workflow.utils import create_sub_template
from workflow.utils.get_resources import get_master_cat


def create_ifu_driver_template(catalogue_template, output_filename,
                               update_datetime=False, overwrite=False):
    """
    Create a template for the IFU driver catalogues.

    Parameters
    ----------
    catalogue_template : str
        Any catalogue template containing the SPA columns.
    output_filename : str
        The name of the output file for the new IFU driver template.
    update_datetime : bool, optional
        Update DATETIME keyword in the new IFU driver template.
    overwrite : bool, optional
        Overwrite the output FITS file containing the IFU driver template.
    """

    # Set the list of columns which will be added to the table
    
    col_list = ['TARGSRVY', 'TARGPROG', 'TARGID', 'TARGNAME', 'TARGPRIO',
                'PROGTEMP', 'OBSTEMP', 'GAIA_RA', 'GAIA_DEC', 'GAIA_EPOCH',
                'GAIA_PMRA', 'GAIA_PMDEC', 'GAIA_PARAL', 'IFU_PA', 'IFU_DITHER']

    # Set the extension name

    extname = 'IFU DRIVER CATALOGUE'

    # Choose the keywords which will be inherited from the catalogue template

    inherited_kwds = ['DATAMVER', 'TRIMESTE', 'DATETIME']

    # Set the information of the new keywords that will be added

    new_primary_kwds = collections.OrderedDict()
    new_primary_kwds['VERBOSE'] = \
        (1, 'Attribute "report_verbosity" of the XML files')
    new_primary_kwds['AUTHOR'] = ('', None)
    new_primary_kwds['CCREPORT'] = ('', None)

    # Create a dictionary to renaming the column IFU_PA
    
    rename_col_dict = {'IFU_PA': 'IFU_PA_REQUEST'}

    # Create the sub-template with the above paremeters
    
    create_sub_template(catalogue_template, output_filename, col_list,
                        extname=extname, inherit_primary_kwds=False,
                        inherited_kwds=inherited_kwds,
                        new_primary_kwds=new_primary_kwds,
                        rename_col_dict=rename_col_dict,
                        update_datetime=update_datetime, overwrite=overwrite)

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Create a template for the IFU driver catalogues')

    parser.add_argument('--in', dest='catalogue_template',
                        default='aux/Master_CatalogueTemplate.fits',
                        help="""name of catalogue template containing the SPA
                        columns""")

    parser.add_argument('--out', dest='ifu_driver_template',
                        default='aux/ifu_driver_template.fits',
                        help="""name for the output file which will contain the
                        new template for the IFU driver catalogues""")

    parser.add_argument('--update_datetime', dest='update_datetime',
                        action='store_true',
                        help="""update DATETIME keyword from the catalogue
                        template""")

    parser.add_argument('--overwrite', dest='overwrite',
                        action='store_true',
                        help='overwrite the output file')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()
    
    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}
    
    logging.basicConfig(level=level_dict[args.log_level])

    if not os.path.exists(args.catalogue_template):
        logging.info('Downloading the master catalogue template')
        get_master_cat(file_path=args.catalogue_template)

    create_ifu_driver_template(args.catalogue_template,
                               args.ifu_driver_template,
                               update_datetime=args.update_datetime,
                               overwrite=args.overwrite)

