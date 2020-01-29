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
from datetime import date
import logging
import os.path
import re
import urllib.request

from astropy.io import fits


def _get_master_cat(file_path='Master_CatalogueTemplate.fits',
                    url=('http://casu.ast.cam.ac.uk/weave/data_model/' +
                         'Master_CatalogueTemplate.fits')):

    urllib.request.urlretrieve(url, file_path)


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
        
    
def create_pi_template(catalogue_template, output_filename, col_list,
                       extname='CATALOGUE BINARY TABLE', primary_kwds={},
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
    primary_kwds: dict, optional
        A dictionary with the keywords and their values which should be added or
        updated in the primary header.
    overwrite : bool, optional
        Overwrite the output FITS file containing the PI template.
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

            col_name = col.name
            col_format = col.format
            col_disp = col.disp
            col_unit = col.unit
            col_null = col.null

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

    hdu.name = extname
    
    # Create the primary header copying its keywords from the template
    
    primary_hdr = fits.Header()
    
    basic_kwd_list =['SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'COMMENT']
    
    for kwd in template_primary_hdr.keys():
        if kwd not in basic_kwd_list:
            primary_hdr[kwd] = template_primary_hdr[kwd]
            primary_hdr.comments[kwd] = template_primary_hdr.comments[kwd]
    
    # Add/Overwrite the requested keywords
    
    for kwd in primary_kwds.keys():
        primary_hdr[kwd] = template_primary_hdr[kwd]

    # Create the primary HDU
    
    primary_hdu = fits.PrimaryHDU(header=primary_hdr)
    
    # Create a HDU list and save it to a file
    
    hdulist = fits.HDUList([primary_hdu, hdu])

    hdulist.writeto(output_filename, overwrite=overwrite)


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
                        help='name of the columns to be inclueded in the table')

    parser.add_argument('--extname', dest='extname',
                        default='CATALOGUE BINARY TABLE',
                        help='name to be used in the table extension')

    parser.add_argument('--keep_date', dest='keep_date',
                        action='store_true',
                        help='keep the date keyword from the master template')

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

    col_list = args.cols_str.split(',')

    if args.keep_date is True:
        primary_kwds = {}
    else:
        primary_kwds = {'DATE': date.today().strftime('%Y-%m-%d')}

    create_pi_template(args.catalogue_template, args.pi_template, col_list,
                       extname=args.extname, primary_kwds=primary_kwds,
                       overwrite=args.overwrite)

