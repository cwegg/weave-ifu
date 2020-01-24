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


from astropy.io import fits


def populate_fits_table_template(fits_template, data_dict, output_filename,
                                 overwrite=False, kwd_value_list=[]):
    
    # Read the FITS template
    
    template_hdulist = fits.open(fits_template)
    template_primary_hdu = template_hdulist[0]
    template_hdu = template_hdulist[1]
    
    # Check that all the columns are available in the dictionary with the data
    
    template_column_names = [col.name for col in template_hdu.columns]
    
    for col_name in template_column_names:
        assert col_name in data_dict.keys()
    
    # Create list of columns as described in the template populated with the
    # provided data
    
    column_list = []
    
    for col in template_hdu.columns:
        
        column = fits.Column(name=col.name, format=col.format,
                             disp=col.disp, unit=col.unit, null=col.null,
                             array=data_dict[col.name])
        
        column_list.append(column)
    
    # Create a HDU from the column list
    
    coldefs = fits.ColDefs(column_list)
    hdu = fits.BinTableHDU.from_columns(coldefs)
    
    # Copy the header from the template
    
    hdu.header = template_hdu.header
    
    # Create the primary extension and populate it with the provided information
    
    primary_hdu = template_primary_hdu
    
    for kwd, value in kwd_value_list:
        primary_hdu.header[kwd] = value
    
    # Create a HDU list and save it to a file
    
    hdulist = fits.HDUList([primary_hdu, hdu])
    
    hdulist.writeto(output_filename, overwrite=overwrite)

