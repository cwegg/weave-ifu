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


def get_obsmode_from_progtemp(progtemp):
    """
    Get the obsmode value from a given progtemp.

    Parameters
    ----------
    progtemp : str
        A progtemp value.
    
    Returns
    -------
    obsmode : {'MOS', 'LIFU', 'mIFU'}
        A progtemp value.
    """

    first_char = progtemp[0]

    if first_char in ['1', '2', '3']:
        obsmode = 'MOS'
    elif first_char in ['4', '5', '6']:
        obsmode = 'LIFU'
    elif first_char in ['7', '8', '9']:
        obsmode = 'mIFU'
    else:
        raise ValueError

    return obsmode

