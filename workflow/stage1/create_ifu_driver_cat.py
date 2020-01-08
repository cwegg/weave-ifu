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


import os
import glob

import numpy as np

from get_data_from_xmls import get_target_data_from_xmls
from populate_fits_table_template import populate_fits_table_template


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
    # Read a dictionary with the data which will populate the template
    
    # For this example, we cheat getting them from stage-5 XMLs. Obviously,
    # this will not be possible in the real SWG workflow, but we provide this
    # example to illustrate how the data should be provided to the function
    # "populate_fits_table_template"
    
    xml_filename_list = glob.glob('../stage5/input/*.xml')
    xml_filename_list.sort()
    
    data_dict = get_target_data_from_xmls(xml_filename_list)
    
    # Add a column with the IFU_PA_REQUEST values. It will be 0 for LIFU and
    # NULL for mIFU
    
    data_dict['IFU_PA_REQUEST'] = [0. if progtemp[0] in ['4', '5', '6']
                                   else np.nan
                                   for progtemp in data_dict['PROGTEMP']]
    
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

