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


from urllib.request import urlretrieve as _urlretrieve


def get_master_cat(file_path='Master_CatalogueTemplate.fits',
                   url=('http://casu.ast.cam.ac.uk/weave/data_model/cats/weave/' +
                        'Master_CatalogueTemplate.fits')):
    """
    Download the latest version of the master catalogue.

    Parameters
    ----------
    file_path : str, optional
        The path used to save the downloaded file.
    url : str, optional
        The URL with the location ot the file in Internet.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    _urlretrieve(url, file_path)
    
    return file_path


def get_aladin_jar(file_path='Aladin.jar',
                   url='https://aladin.u-strasbg.fr/java/Aladin.jar'):
    """
    Download the latest version of Aladin.

    Parameters
    ----------
    file_path : str, optional
        The path used to save the downloaded file.
    url : str, optional
        The URL with the location ot the file in Internet.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    _urlretrieve(url, file_path)
    
    return file_path

