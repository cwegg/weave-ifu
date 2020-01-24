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
import glob

import numpy as np

from populate_fits_table_template import populate_fits_table_template


def _get_data_dict():
    
    data_dict = {}
    
    data_dict['TARGSRVY'] = []
    
    data_dict['TARGPROG'] = []
    
    data_dict['TARGID'] = []
    
    data_dict['TARGNAME'] = []
    
    data_dict['TARGPRIO'] = []
    
    data_dict['PROGTEMP'] = []
    
    data_dict['OBSTEMP'] = []
    
    data_dict['GAIA_RA'] = []
    
    data_dict['GAIA_DEC'] = []
    
    data_dict['GAIA_EPOCH'] = []
    
    data_dict['GAIA_PMRA'] = []
    
    data_dict['GAIA_PMDEC'] = []
    
    data_dict['GAIA_PARAL'] = []
    
    data_dict['IFU_PA_REQUEST'] = []
    
    data_dict['IFU_DITHER'] = []
    
    return data_dict


if __name__ == '__main__':
    
    ############################################################################
    # Set the location of the template and the output file and directory
    
    ifu_driver_template = './aux/ifu_driver_template.fits'
    
    output_dir = './output/'
    output_filename = output_dir + 'WC_IFU.fits'
    
    ############################################################################
    # Create the output directory if it does not exist
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    ############################################################################
    # Get a dictionary with the data which will populate the template
    
    # NOTE: See the above function _get_data_dict to understand the structure of
    # the dictionary
    
    data_dict = _get_data_dict()
    
    ############################################################################
    # Create a list with the information needed to populate the keywords of
    # of the primary header of the template
    
    kwd_value_list = [
        ('AUTHOR',   'jairomendezabreu@gmail.com'),
        ('CCREPORT', 'daniela.bettoni@oapd.inaf.it,jalfonso@iac.es'),
        ('VERBOSE',   1)
    ]
    
    ############################################################################
    # Populate the template with the provided data
    
    populate_fits_table_template(ifu_driver_template, data_dict,
                                 output_filename, kwd_value_list=kwd_value_list)

