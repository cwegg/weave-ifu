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

import numpy as np
from astropy.io import fits

from workflow.utils import populate_fits_table_template
from workflow.utils.get_data_from_xmls import get_spa_data_of_target_and_sky_fibres_from_xmls


def _get_col_null_dict_of_template(fits_template):
    
    # Get an empty dictionary for the result
    
    col_null_dict = {}
    
    # Open the template and get list of keywords in extension 1
    
    with fits.open(fits_template) as hdu_list:
    
        hdu = hdu_list[1]
        hdr = hdu.header
        kwd_list = list(hdr.keys())
        
        # For each TTYPE keyword
        
        for kwd in kwd_list:
            
            if re.match('TTYPE([0-9]+)', kwd):
                
                # Get the column number
                
                num_str = kwd[5:]
                
                # Get the values of TTYPE and TFORM
                
                ttype_value = hdr['TTYPE{}'.format(num_str)]
                tform_value = hdr['TFORM{}'.format(num_str)]
                
                # Get the name of the keyword which COULD containe the TNULL
                
                tnull_kwd = 'TNULL{}'.format(num_str)
                
                # Get the TNULL depending on the TFORM value:
                #  - For integers, get it from the TNULL keyword
                #  - For reals, it will be NaN
                #  - For strings, it will be an empty string
                #  - If fails, assign None
                
                try:
                    if (re.match('([1]?)I', tform_value) or
                        re.match('([1]?)J', tform_value) or
                        re.match('([1]?)K', tform_value)):
                        
                        assert tnull_kwd in kwd_list
                        tnull_value = hdr[tnull_kwd]
                        
                    elif (re.match('([1]?)E', tform_value) or
                          re.match('([1]?)D', tform_value)):
                        
                        assert tnull_kwd not in kwd_list
                        tnull_value = np.nan
                        
                    elif re.match('([0-9]*)A', tform_value):
                    
                        assert tnull_kwd not in kwd_list
                        tnull_value = ''
                        
                    else:
                    
                        raise ValueError
                except:
                    tnull_value = None
                
                # Save the column name and its TNULL value
                
                col_null_dict[ttype_value] = tnull_value
    
    return col_null_dict


def _get_col_len(data_dict):

    # Get the lenght from one random column in the dictionary

    random_key = list(data_dict.keys())[0]
    col_len = len(data_dict[random_key])
    
    # Assert that all the columns in the dictionary have the same lenght
    
    for key in data_dict.keys():
        assert len(data_dict[key]) == col_len
    
    return col_len


def _add_missing_cols_of_template(data_dict, fits_template):
    
    # Get lenght of each column contained in the data dictionary
    
    col_len = _get_col_len(data_dict)
    
    # Get a dictionary with the names of the columns in the template and their
    # null values
    
    col_null_dict = _get_col_null_dict_of_template(fits_template)
    
    assert None not in col_null_dict.values()
    
    # For each missing column, add a list of nulls to the dictionary    

    out_data_dict = data_dict.copy()

    for col_name in col_null_dict.keys():
        if col_name not in data_dict.keys():
            out_data_dict[col_name] = [col_null_dict[col_name]] * col_len
    
    return out_data_dict


def create_ifu_fits_cat(xml_files, fits_template, output_filename,
                        overwrite=False):
    """
    Create an IFU FITS catalogue.
    
    Parameters
    ----------
    xml_files : str or list of str
        A string with a pattern of the input XML files or list of strings with
        the filenames of the XML files.
    fits_template : list of str
        A FITS template with a primary HDU and a first extension with a table.
    output_filename : str
        The name of the output file which will be created.
    overwrite : bool, optional
        Overwrite the output FITS file.
    """
    
    # Get list of filenames depending on the type of input given for XML files
    
    if type(xml_files) == list:
        xml_filename_list = xml_files
    elif type(xml_files) == str:
        xml_filename_list = glob.glob(xml_files)
        xml_filename_list.sort()
    else:
        raise TypeError
    
    # Get a dictionary with the SPA data from targets and skies in the XML files
    
    spa_data_dict = get_spa_data_of_target_and_sky_fibres_from_xmls(
                         xml_filename_list)
    
    # Add missing columns of the template to the dictionary
    
    data_dict = _add_missing_cols_of_template(spa_data_dict, fits_template)
    
    # Populate the template
    
    populate_fits_table_template(fits_template, data_dict, output_filename,
                                 update_datetime=True, overwrite=overwrite)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create an IFU FITS catalogue')

    parser.add_argument('template',
                        help="""a FITS template with a primary HDU and a first
                        extension with a table""")

    parser.add_argument('xml_file', nargs='+', help='name of a XML file')

    parser.add_argument('--out', dest='output_filename',
                        default='output/cat-ifu_from_xmls.fits',
                        help="""name for the output file which will contain the
                        IFU FITS catalogue""")

    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='overwrite the output file')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='The level for the logging messages')

    args = parser.parse_args()

    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}

    logging.basicConfig(level=level_dict[args.log_level])

    if not os.path.exists(os.path.dirname(args.output_filename)):
        logging.info('Creating the output directory')
        os.mkdir(os.path.dirname(args.output_filename))

    create_ifu_fits_cat(args.xml_file, args.template, args.output_filename,
                        overwrite=args.overwrite)

