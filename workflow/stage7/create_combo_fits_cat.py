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


from workflow.utils import check_equal_fits_structure
from workflow.utils import populate_fits_table_template

from astropy.table import Table, vstack


def _get_combo_table(fits_filename1, fits_filename2):
    
    table1 = Table.read(fits_filename1, format='fits')
    table2 = Table.read(fits_filename2, format='fits')
    
    combo_table = vstack([table1, table2])
    
    return combo_table


def _convert_table_to_data_dict(combo_table):
    
    data_dict = {}
    
    for col_name in combo_table.colnames:
        data_dict[col_name] = list(combo_table[col_name].data)
    
    return data_dict


def _get_combo_data_dict(fits_filename1, fits_filename2):
    
    combo_table = _get_combo_table(fits_filename1, fits_filename2)
    
    data_dict = _convert_table_to_data_dict(combo_table)
    
    return data_dict


def create_combo_fits_cat(mos_cat, ifu_cat, output_filename, overwrite=False):
    
    assert check_equal_fits_structure(mos_cat, ifu_cat)
    
    combo_data_dict = _get_combo_data_dict(mos_cat, ifu_cat)
    
    populate_fits_table_template(mos_cat, combo_data_dict, output_filename,
                                 update_datetime=True, overwrite=overwrite)

