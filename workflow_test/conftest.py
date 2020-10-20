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


import glob
import os
import pathlib

import pytest

import workflow


@pytest.fixture(scope='session')
def aladin_jar(tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join('Aladin.jar'))
    
    workflow.utils.get_resources.get_aladin_jar(file_path=file_path)
    
    result = str(file_path)
    
    assert os.path.exists(result)
    
    return result


@pytest.fixture(scope='session')
def master_cat(tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join(
                        'Master_CatalogueTemplate.fits'))
    
    workflow.utils.get_resources.get_master_cat(file_path=file_path)
    
    result = str(file_path)
    
    assert os.path.exists(result)
    
    return result


def test_fitscheck_master_cat(master_cat):

    returncode = subprocess.call(['fitscheck', master_cat])
    
    assert returncode == 0


@pytest.fixture(scope='session')
def pkg_ifu_driver_template():

    pkg_file_path = str(pathlib.Path(workflow.__path__[0]) / 'stage1' / 'aux' /
                        'ifu_driver_template.fits')

    assert os.path.exists(pkg_file_path)
    
    return pkg_file_path


def test_fitscheck_pkg_ifu_driver_template(pkg_ifu_driver_template):

    returncode = subprocess.call(['fitscheck', pkg_ifu_driver_template])
    
    assert returncode == 0


@pytest.fixture(scope='session')
def pkg_ifu_driver_cat():

    pkg_file_path = str(pathlib.Path(workflow.__path__[0]) / 'stage2' /
                        'input' / 'WC_2020A1-ifu_driver_cat.fits')

    assert os.path.exists(pkg_file_path)
    
    return pkg_file_path


def test_fitscheck_pkg_ifu_driver_cat(pkg_ifu_driver_cat):

    returncode = subprocess.call(['fitscheck', pkg_ifu_driver_cat])
    
    assert returncode == 0


@pytest.fixture(scope='session')
def pkg_tgc_xml_files():
    
    xml_files_pattern = str(pathlib.Path(workflow.__path__[0]) / 'stage4' /
                            'input' / '*-tgc.xml')
    
    xml_filename_list = glob.glob(xml_files_pattern)
    xml_filename_list.sort()

    assert len(xml_filename_list) > 0
    
    for xml_filename in xml_filename_list:
        assert os.path.exists(xml_filename)
    
    return xml_filename_list

