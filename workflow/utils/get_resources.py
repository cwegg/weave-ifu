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


import os as _os
import urllib as _urllib


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
    
    _urllib.request.urlretrieve(url, file_path)
    
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
    
    _urllib.request.urlretrieve(url, file_path)
    
    return file_path


def get_blank_xml_template(file_path='BlankXMLTemplate.xml',
                           url=('http://casu.ast.cam.ac.uk/weave/data_model/' +
                                'misc/BlankXMLTemplate.xml')):
    """
    Download the latest version of the blank XML template.

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
    
    _urllib.request.urlretrieve(url, file_path)
    
    return file_path


def get_progtemp_file(file_path='progtemp.dat',
                      url=('http://casu.ast.cam.ac.uk/weave/data_model/misc/' +
                           'progtemp.dat')):
    """
    Download the latest version of progtemp.dat.

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
    
    _urllib.request.urlretrieve(url, file_path)
    
    return file_path


def get_obstemp_file(file_path='obstemp.dat',
                     url=('http://casu.ast.cam.ac.uk/weave/data_model/misc/' +
                          'obstemp.dat')):
    """
    Download the latest version of obstemp.dat.

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
    
    _urllib.request.urlretrieve(url, file_path)
    
    return file_path


def _get_srvy_cat(survey, healpix_index, directory='',
                  url_template='http://hornet.ast.cam.ac.uk/{}/getHPID'):
    """
    Download a catalogue of stars for the requested survey and HEALPix index.

    Parameters
    ----------
    survey : {'WD', 'WG'}
        The survey to be considered.
    healpix_index : int
        The desired HEALPix index.
    directory : str, optional
        The directory used to save the downloaded file.
    url_template : str, optional
        The URL with the location ot the catalogue in Internet.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """

    assert survey in ['WD', 'WG']

    filename = '{}_{}.fits'.format(survey, healpix_index)
    file_path = _os.path.join(directory, filename)

    url = url_template.format(survey.lower())
    data = _urllib.parse.urlencode({'id': healpix_index}).encode('ascii')

    with open(file_path, 'wb') as f:
        response = _urllib.request.urlopen(url, data=data)
        contents = response.read()
        f.write(contents)
    
    return file_path


def get_guide_cat(healpix_index, directory=''):
    """
    Download a catalogue of guide stars for the requested HEALPix index.

    Parameters
    ----------
    healpix_index : int
        The desired HEALPix index.
    directory : str, optional
        The directory used to save the downloaded file.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    file_path = _get_srvy_cat('WG', healpix_index, directory=directory)
    
    return file_path


def get_calib_cat(healpix_index, directory=''):
    """
    Download a catalogue of calib stars for the requested HEALPix index.

    Parameters
    ----------
    healpix_index : int
        The desired HEALPix index.
    directory : str, optional
        The directory used to save the downloaded file.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    file_path = _get_srvy_cat('WD', healpix_index, directory=directory)
    
    return file_path

