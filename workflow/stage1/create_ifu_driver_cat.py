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


from workflow.utils import populate_fits_table_template


def create_ifu_driver_cat(ifu_driver_template, data_dict, output_filename,
                          trimester, author, report_verbosity=1, cc_report='',
                          overwrite=False):
    """
    Create a IFU driver catalogue using a template and the needed information.

    Parameters
    ----------
    ifu_driver_template : str
        A FITS file containing an IFU driver template.
    data_dict : dict
        A dictionary with the information needed to populate the columns of the
        IFU driver template. The keys of the dictionary should be the column
        names, while their values should be list or array-like objects with the
        information of each column.
    output_filename : str
        The name of the output file for the new IFU driver catalogue.
    trimester : str
        The trimester of the catalogue (e.g. 2020A1).
    author : str
        The email address of the author of the catalogue.
    report_verbosity : {0, 1}, optional
        The level of verbosity which will be inherited in the files to be
        submitted to WASP.
    cc_report : str, optional
        A comma separated list of email addresses to be CC'ed in WASP
        submisions.
    overwrite : bool, optional
        Overwrite the output FITS file containing the IFU driver catalogue.
    """

    assert report_verbosity in [0, 1]

    primary_kwds = {'TRIMESTE': trimester, 'VERBOSE': report_verbosity,
                    'AUTHOR': author, 'CCREPORT': cc_report}
    
    populate_fits_table_template(ifu_driver_template, data_dict,
                                 output_filename, primary_kwds=primary_kwds,
                                 update_datetime=True, overwrite=overwrite)

