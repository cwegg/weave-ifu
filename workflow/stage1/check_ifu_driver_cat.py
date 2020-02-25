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


def check_ifu_driver_cat(cat_filename):
    """
    Check the contents of an IFU driver catalogue.

    Parameters
    ----------
    cat_filename : str
        A FITS file with an IFU driver catalogue.

    Returns
    ----------
    result : bool
        True if the file passes all the checks, otherwise False.
    """
    
    result = True
    
    return result

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description='Check the contents of an IFU driver catalogue')
    
    parser.add_argument('catalogue',
                        help='a FITS file with an IFU driver catalogue')
    
    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')
    
    args = parser.parse_args()
    
    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}
    
    logging.basicConfig(level=level_dict[args.log_level])
    
    check_ifu_driver_cat(args.catalogue)

