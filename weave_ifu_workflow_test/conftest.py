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

import weave_ifu_workflow
import weave_utils

@pytest.fixture(scope='session')
def aladin_jar(tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join('Aladin.jar'))
    
    weave_utils.get_resources.get_aladin_jar(file_path=file_path)
    
    assert os.path.exists(file_path)
    
    return file_path


@pytest.fixture(scope='session')
def master_cat(tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join(
                        'Master_CatalogueTemplate.fits'))
    
    weave_utils.get_resources.get_master_cat(file_path=file_path)
    
    assert os.path.exists(file_path)
    
    return file_path


@pytest.fixture(scope='session')
def blank_xml_template(tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join(
                        'BlankXMLTemplate.xml'))
    
    weave_utils.get_resources.get_blank_xml_template(file_path=file_path)
    
    assert os.path.exists(file_path)
    
    return file_path


@pytest.fixture(scope='session')
def progtemp_file(tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join(
                        'progtemp.dat'))
    
    weave_utils.get_resources.get_progtemp_file(file_path=file_path)
    
    assert os.path.exists(file_path)
    
    return file_path


@pytest.fixture(scope='session')
def obstemp_file(tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join(
                        'obstemp.dat'))
    
    weave_utils.get_resources.get_obstemp_file(file_path=file_path)
    
    assert os.path.exists(file_path)
    
    return file_path


@pytest.fixture(scope='session')
def pkg_wc_cat():

    pkg_file_path = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage5' / 'aux' /
                        'WC_CatalogueTemplate.fits')
    
    assert os.path.exists(pkg_file_path)
    
    return pkg_file_path


@pytest.fixture(scope='session')
def pkg_ifu_driver_template():

    pkg_file_path = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage1' / 'aux' /
                        'ifu_driver_template.fits')

    assert os.path.exists(pkg_file_path)
    
    return pkg_file_path


@pytest.fixture(scope='session')
def pkg_ifu_driver_cat():

    pkg_file_path = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage2' /
                        'input' / 'WC_2020A1-ifu_driver_cat.fits')

    assert os.path.exists(pkg_file_path)
    
    return pkg_file_path


@pytest.fixture(scope='session')
def pkg_t_xml_files():
    
    xml_files_pattern = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage3' /
                            'input' / '*-t.xml')
    
    xml_filename_list = glob.glob(xml_files_pattern)
    xml_filename_list.sort()

    assert len(xml_filename_list) > 0
    
    for xml_filename in xml_filename_list:
        assert os.path.exists(xml_filename)
    
    return xml_filename_list


@pytest.fixture(scope='session')
def pkg_tgc_xml_files():
    
    xml_files_pattern = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage4' /
                            'input' / '*-tgc.xml')
    
    xml_filename_list = glob.glob(xml_files_pattern)
    xml_filename_list.sort()

    assert len(xml_filename_list) > 0
    
    for xml_filename in xml_filename_list:
        assert os.path.exists(xml_filename)
    
    return xml_filename_list


@pytest.fixture(scope='session')
def pkg_tgcs_xml_files():
    
    xml_files_pattern = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage5' /
                            'input' / '*-tgcs.xml')
    
    xml_filename_list = glob.glob(xml_files_pattern)
    xml_filename_list.sort()

    assert len(xml_filename_list) > 0
    
    for xml_filename in xml_filename_list:
        assert os.path.exists(xml_filename)
    
    return xml_filename_list


@pytest.fixture(scope='session')
def pkg_ifu_cat_from_xmls():

    pkg_file_path = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage6' /
                        'input' / 'WC_2020A1-ifu_from_xmls.fits')

    assert os.path.exists(pkg_file_path)
    
    return pkg_file_path


@pytest.fixture(scope='session')
def pkg_ifu_cat():

    pkg_file_path = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage7' /
                        'input' / 'WC_2020A1-ifu.fits')

    assert os.path.exists(pkg_file_path)
    
    return pkg_file_path


@pytest.fixture(scope='session')
def pkg_mos_cat():

    pkg_file_path = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage7' /
                        'input' / 'WC_2020A1-mos.fits')

    assert os.path.exists(pkg_file_path)
    
    return pkg_file_path


@pytest.fixture(scope='session')
def pkg_combo_cat():

    pkg_file_path = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage8' /
                        'input' / 'WC_2020A1.fits')

    assert os.path.exists(pkg_file_path)
    
    return pkg_file_path


@pytest.fixture(scope='session')
def pkg_wasp_cat():

    pkg_file_path = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage9' /
                        'input' / 'WC_2020A1.fits')

    assert os.path.exists(pkg_file_path)
    
    return pkg_file_path


@pytest.fixture(scope='session')
def pkg_copy_tgcs_xml_files():
    
    xml_files_pattern = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage9' /
                            'input' / '*-tgcs.xml')
    
    xml_filename_list = glob.glob(xml_files_pattern)
    xml_filename_list.sort()

    assert len(xml_filename_list) > 0
    
    for xml_filename in xml_filename_list:
        assert os.path.exists(xml_filename)
    
    return xml_filename_list


@pytest.fixture(scope='session')
def pkg_wasp_xml_files():
    
    xml_files_pattern = str(pathlib.Path(weave_ifu_workflow.__path__[0]) / 'stage10' /
                            'input' / '*.xml')
    
    xml_filename_list = glob.glob(xml_files_pattern)
    xml_filename_list.sort()

    assert len(xml_filename_list) > 0
    
    for xml_filename in xml_filename_list:
        assert os.path.exists(xml_filename)
    
    return xml_filename_list

