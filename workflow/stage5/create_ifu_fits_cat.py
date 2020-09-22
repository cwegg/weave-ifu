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
from workflow.utils.get_data_from_xmls import (
    get_datamver_from_xmls,
    get_author_from_xmls, get_cc_report_from_xmls, get_trimester_from_xmls,
    get_spa_data_of_target_random_and_sky_fibres_from_xmls)

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


def _add_sim_extension(filename):

    # Open the FITS file in append mode

    with fits.open(filename, mode='append') as hdu_list:

        # Get the number of rows in the first extension

        num_rows = len(hdu_list[1].data)

        # Create the column defitions and its corresponding list of UCDs

        sim_col_list = []
        sim_ucd_list = []

        sim_col_list.append(
            fits.Column(name='SIM_TEMPLATE', format='30A',
                        array=['' for i in range(num_rows)]))
        sim_ucd_list.append('phot.flux.density;meta.file')

        sim_col_list.append(
            fits.Column(name='SIM_MAG', format='D', unit='mag',
                        array=[np.nan for i in range(num_rows)]))
        sim_ucd_list.append('phot.mag')

        sim_col_list.append(
            fits.Column(name='SIM_FILTERID', format='10A',
                        array=['' for i in range(num_rows)]))
        sim_ucd_list.append('instr.filter')

        sim_col_list.append(
            fits.Column(name='SIM_FWHM', format='D', unit='arcsec',
                        array=[np.nan for i in range(num_rows)]))
        sim_ucd_list.append('instr.obsty.seeing')

        sim_col_list.append(
            fits.Column(name='SIM_VELOCITY', format='D', unit='km/s',
                        array=[np.nan for i in range(num_rows)]))
        sim_ucd_list.append('phys.veloc')

        sim_col_list.append(
            fits.Column(name='SIM_REDSHIFT', format='D',
                        array=[np.nan for i in range(num_rows)]))
        sim_ucd_list.append('src.redshift')

        # Create an HDU from the column definitions

        sim_hdu = fits.BinTableHDU.from_columns(sim_col_list)

        # Add TUCD keywords to the HDU in nice positions

        for i, ucd in enumerate(sim_ucd_list):

            ucd_kwd = 'TUCD{}'.format(i + 1)

            unit_kwd = 'TUNIT{}'.format(i + 1)
            form_kwd = 'TFORM{}'.format(i + 1)

            if unit_kwd in sim_hdu.header:
                after_kwd = unit_kwd
            else:
                after_kwd = form_kwd

            sim_hdu.header.set(ucd_kwd, ucd, after=after_kwd)

        # Add an extension name to the HDU

        sim_hdu.name = 'SIM'

        # Append the sim HDU

        hdu_list.append(sim_hdu)


def create_ifu_fits_cat(xml_files, fits_template, output_filename,
                        cat_nme1='', cat_nme2='', sim_ext=False,
                        overwrite=False):
    """
    Create an IFU FITS catalogue.
    
    Parameters
    ----------
    xml_files : str or list of str
        A string with a pattern of the input XML files or list of strings with
        the filenames of the XML files.
    fits_template : str
        A FITS template with a primary HDU and a first extension with a table.
    output_filename : str
        The name of the output file which will be created.
    cat_nme1 : str, optional
        Value for populating CAT_NME1 keyword of the output file.
    cat_nme2 : str, optional
        Value for populating CAT_NME2 keyword of the output file.
    sim_ext : bool, optional
        Add an extra extension to the output file for the information needed to
        generate the simulations.
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

    # Check that DATAMVER is consistent bewteen the XML files and the template

    xml_datamver = get_datamver_from_xmls(xml_filename_list)
    template_datamver = fits.getval(fits_template, 'DATAMVER')

    if template_datamver != xml_datamver:
        logging.critical(
            'DATAMVER mismatch ({} != {}) '.format(
                xml_datamver, template_datamver) +
            'between XML files and FITS template: ' +
            'Stop unless you are sure!')
        raise SystemExit(2)
    
    # Get the trimester, author and cc_report of the files (which should be the
    # same for all them)
    
    author = get_author_from_xmls(xml_filename_list)
    cc_report = get_cc_report_from_xmls(xml_filename_list)
    trimester = get_trimester_from_xmls(xml_filename_list)
    
    primary_kwds = {
        'CAT_NME1': cat_nme1,
        'CAT_NME2': cat_nme2,
        'CAT_MAIL': author,
        'CAT_CC': cc_report,
        'TRIMESTE': trimester
    }
    
    # Get a dictionary with the SPA data from targets and skies in the XML files
    
    spa_data_dict = get_spa_data_of_target_random_and_sky_fibres_from_xmls(
                        xml_filename_list)
    
    # Add missing columns of the template to the dictionary
    
    data_dict = _add_missing_cols_of_template(spa_data_dict, fits_template)
    
    # Populate the template
    
    populate_fits_table_template(fits_template, data_dict, output_filename,
                                 primary_kwds=primary_kwds,
                                 update_datetime=True, overwrite=overwrite)

    # Add an extra extension for the simulations if requested

    if sim_ext == True:
        _add_sim_extension(output_filename)


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

    parser.add_argument('--cat_nme1', default='',
                        help="""value for populating CAT_NME1 keyword of the
                        output file""")

    parser.add_argument('--cat_nme2', default='',
                        help="""value for populating CAT_NME1 keyword of the
                        output file""")

    parser.add_argument('--sim_ext',  action='store_true',
                        help="""add an extra extension to the output file for
                        the information needed to generate the simulations""")

    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite the output file')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}

    logging.basicConfig(level=level_dict[args.log_level])

    if not os.path.exists(os.path.dirname(args.output_filename)):
        logging.info('Creating the output directory')
        os.mkdir(os.path.dirname(args.output_filename))

    create_ifu_fits_cat(args.xml_file, args.template, args.output_filename,
                        cat_nme1=args.cat_nme1, cat_nme2=args.cat_nme2,
                        sim_ext=args.sim_ext, overwrite=args.overwrite)

