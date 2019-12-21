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


def _try_to_copy_comment(kwd_type, hdu, i, template_hdu, j):

    kwd_i = '{}{}'.format(kwd_type, i + 1)
    kwd_j = '{}{}'.format(kwd_type, j + 1)

    if kwd_i in hdu.header.keys():

        hdu.header.comments[kwd_i] = template_hdu.header.comments[kwd_j]


def _try_to_copy_keyword(kwd_type, hdu, i, template_hdu, j):

    kwd_i = '{}{}'.format(kwd_type, i + 1)
    kwd_j = '{}{}'.format(kwd_type, j + 1)

    if kwd_j in template_hdu.header.keys():

        # Copy the keyword before the next TTYPE kwd if possible, otherwise at
        # the end of the header

        next_ttype_kwd = 'TTYPE{}'.format(i + 2)

        if next_ttype_kwd in hdu.header.keys():
            hdu.header.insert(next_ttype_kwd,
                              (kwd_i, template_hdu.header[kwd_j]))
        else:
            hdu.header[kwd_i] = template_hdu.header[kwd_j]

        # Copy the comment of the keyword

        hdu.header.comments[kwd_i] = template_hdu.header.comments[kwd_j]
        
    
def create_ifu_driver_template(catalogue_template, output_filename,
                               col_list=['TARGID', 'TARGNAME', 'TARGPRIO',
                                         'PROGTEMP', 'OBSTEMP',
                                         'GAIA_RA', 'GAIA_DEC', 'GAIA_EPOCH',
                                         'GAIA_PMRA', 'GAIA_PMDEC',
                                         'GAIA_PARAL', 'IFU_PA', 'IFU_DITHER'],
                               rename_dict={'IFU_PA': 'IFU_PA_REQUEST'},
                               fix_str_format=False, overwrite=False):

    # Lists with the type of keywords that will be copied from the catalogue
    # template

    basic_kwd_type_list = ['TTYPE', 'TFORM', 'TDISP', 'TUNIT', 'TNULL']
    extra_kwd_type_list = ['TDMIN', 'TDMAX', 'TUCD', 'TPROP']

    # Read the catalogue template

    template_hdulist = fits.open(catalogue_template)
    template_hdu = template_hdulist[1]

    # Check that all the requested columns exist in the catalogue template and
    # that there are not repeated columns

    template_column_names = [col.name for col in template_hdu.columns]

    for col_name in col_list:
        assert col_name in template_column_names

    assert len(col_list) == len(set(col_list))

    # Create the column list and save the mapping between the columns of the
    # new IFU driver catalogue template and the input catalogue template

    column_list = []
    column_mapping = {}

    col_counter = 0

    for j, col in enumerate(template_hdu.columns):

        if col.name in col_list:

            # Save the mapping and increase the counter of created columns

            column_mapping[col_counter] = j
            col_counter += 1

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

    # Copy the comments of the original catalogue for the created keywords

    for i in range(len(col_list)):

        for kwd_type in basic_kwd_type_list:

            j = column_mapping[i]

            _try_to_copy_comment(kwd_type, hdu, i, template_hdu, j)

    # Add other keywords related with each column which could not be added
    # in the defitions of the columns

    for i in range(len(col_list)):

        for kwd_type in extra_kwd_type_list:

            j = column_mapping[i]

            _try_to_copy_keyword(kwd_type, hdu, i, template_hdu, j)
    
    # Give a name to the HDU and write it

    hdu.name = 'INPUT IFU DRIVER CATALOGUE'

    hdu.writeto(output_filename, overwrite=overwrite)


if __name__ == '__main__':

    catalogue_template = '../../test_data/stage0_base.fits'
    ifu_driver_template = './aux/ifu_driver_template.fits'

    create_ifu_driver_template(catalogue_template, ifu_driver_template,
                               fix_str_format=True)

