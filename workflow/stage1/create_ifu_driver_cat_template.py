#!/usr/bin/env python3

#
# Copyright (C) 2019 Cambridge Astronomical Survey Unit
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


import re
import logging

from astropy.io import fits


def _get_format_from_disp(col_disp):

    assert re.match('A[0-9]*', col_disp)

    number = col_disp[1:]
    col_format = number + 'A'

    return col_format

    
def create_ifu_driver_cat_template(catalogue_template, output_filename,
                                   col_list=['TARGID', 'TARGNAME', 'TARGPRIO',
                                             'PROGTEMP', 'OBSTEMP',
                                             'GAIA_RA', 'GAIA_DEC',
                                             'GAIA_EPOCH',
                                             'GAIA_PMRA', 'GAIA_PMDEC',
                                             'GAIA_PARAL',
                                             'IFU_PA', 'IFU_DITHER'],
                                    rename_dict={'IFU_PA': 'LIFU_PA_REQUEST'},
                                    fix_str_format=False):

    # Read the catalogue template

    template_hdulist = fits.open(catalogue_template)
    template_hdu = template_hdulist[1]

    # Create the column list

    column_list = []

    for i, col in enumerate(template_hdu.columns):

        if col.name in col_list:

            # Get the properties of the column

            if col.name in rename_dict.keys():
                col_name = rename_dict[col.name]
            else:
                col_name = col.name

            col_format = col.format

            col_disp = col.disp
            col_unit = col.unit
            col_null = col.null

            # Fix the format of the string columns if requested

            if ('A' in col_format) and (fix_str_format is True):
                col_format = _get_format_from_disp(col_disp)

                if col_format != col.format:
                    logging.warning(
                        'format of column {} updated from {} to {}'.format(
                            col_name, col.format, col_format))

            # Create the column and add it to the column list

            column = fits.Column(name=col_name, format=col_format,
                                 disp=col_disp, unit=col_unit, null=col_null)

            column_list.append(column)

    # Create the HDU from the column list

    coldefs = fits.ColDefs(column_list)
    hdu = fits.BinTableHDU.from_columns(coldefs)

    # Add other keywords related with each column which could not be added
    # in the defitions of the columns, and also copy the comments from the
    # original catalogue

    # ***

    # Give a name to the HDU and write it

    hdu.name = 'INPUT IFU DRIVER CATALOGUE'

    hdu.writeto(output_filename)


if __name__ == '__main__':

    catalogue_template = '../../test_data/stage0_base.fits'
    ifu_driver_template = './aux/ifu_driver_template.fits'

    create_ifu_driver_cat_template(catalogue_template, ifu_driver_template,
                                   fix_str_format=True)

