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
import re
import urllib.request

from astropy.io import fits


def _get_master_cat(file_path='Master_CatalogueTemplate.fits',
                    url=('http://casu.ast.cam.ac.uk/weave/data_model/cats/' +
                         'Master_CatalogueTemplate.fits')):

    urllib.request.urlretrieve(url, file_path)


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
                               col_list=['TARGSRVY', 'TARGPROG',
                                         'TARGID', 'TARGNAME', 'TARGPRIO',
                                         'PROGTEMP', 'OBSTEMP',
                                         'GAIA_RA', 'GAIA_DEC', 'GAIA_EPOCH',
                                         'GAIA_PMRA', 'GAIA_PMDEC',
                                         'GAIA_PARAL', 'IFU_PA', 'IFU_DITHER'],
                               rename_dict={'IFU_PA': 'IFU_PA_REQUEST'},
                               fix_str_format=False, overwrite=False):
    """
    Create a template for the input IFU driver catalogues.

    Parameters
    ----------
    catalogue_template : str
        Any catalogue template containing the SPA columns.
    output_filename : str
        The name of the output file for the new IFU driver template.
    col_list : list of str, optional
        A list containing all the columns which will be included in the new IFU
        driver template.
    rename_dict : dict, optional
        A dictionary used to rename the columns. Its keys should be the name of
        the columns in the input catalogue template, while its values should be
        the new name which will be used in the output IFU driver template.
    fix_str_format : bool, optional
        Activate an option for fixing faulty TFORMs of the string columns in
        the input catalogue template. If True, TFORMs will be fixed using the
        values available for TDISP in these columns.
    overwrite : bool, optional
        Overwrite the output FITS file containing the IFU driver template.
    """

    # Lists with the type of keywords that will be copied from the catalogue
    # template

    basic_kwd_type_list = ['TTYPE', 'TFORM', 'TDISP', 'TUNIT', 'TNULL']
    extra_kwd_type_list = ['TDMIN', 'TDMAX', 'TUCD', 'TPROP']

    # Read the catalogue template

    template_hdulist = fits.open(catalogue_template)
    template_primary_hdr = template_hdulist[0].header
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
    
    # Give a name to the HDU

    hdu.name = 'INPUT IFU DRIVER CATALOGUE'
    
    # Create the primary extension to contain some attributes of the XMLs
    
    primary_hdr = fits.Header()
    
    # Add a keyword with the WEAVE data model version to the primary header
    
    primary_hdr['DATAMVER'] = template_primary_hdr['DATAMVER']
    primary_hdr.comments['DATAMVER'] = template_primary_hdr.comments['DATAMVER']
    
    # Add some attributes of the XMLs to the primary header
    
    primary_hdr['VERBOSE'] = 1
    primary_hdr.comments['VERBOSE'] = \
       'Attribute "report_verbosity" of the XML files'
    
    primary_hdr['AUTHOR'] = ''
    
    primary_hdr['CCREPORT'] = ''
    
    primary_hdu = fits.PrimaryHDU(header=primary_hdr)
    
    # Create a HDU list and save it to a file
    
    hdulist = fits.HDUList([primary_hdu, hdu])

    hdulist.writeto(output_filename, overwrite=overwrite)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Create a template for the input IFU driver catalogues')

    parser.add_argument('--in', dest='catalogue_template',
                        default='aux/Master_CatalogueTemplate.fits',
                        help="""name of catalogue template containing the SPA
                        columns""")

    parser.add_argument('--out', dest='ifu_driver_template',
                        default='aux/ifu_driver_template.fits',
                        help="""name for the output file which will contain the
                        new template for the IFU driver catalogues""")

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
        _get_master_cat(file_path=args.catalogue_template)

    create_ifu_driver_template(args.catalogue_template,
                               args.ifu_driver_template,
                               overwrite=args.overwrite)

