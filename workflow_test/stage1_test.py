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


import os
import pathlib
import subprocess

import pytest

import workflow


@pytest.fixture(scope='session')
def aladin_jar(tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join('Aladin.jar'))
    
    workflow.utils.get_resources.get_aladin_jar(file_path=file_path)
    
    return str(file_path)


def test_download_aladin(aladin_jar):

    assert os.path.exists(aladin_jar)


@pytest.fixture(scope='session')
def master_cat(tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join(
                        'Master_CatalogueTemplate.fits'))
    
    workflow.utils.get_resources.get_master_cat(file_path=file_path)
    
    return str(file_path)


def test_download_master_cat(master_cat):

    assert os.path.exists(master_cat)


@pytest.fixture(scope='session')
def ifu_driver_template(master_cat, tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join(
                        'ifu_driver_template.fits'))
    
    workflow.stage1.create_ifu_driver_template(master_cat, file_path)
    
    return file_path


def test_create_ifu_driver_template(ifu_driver_template):

    assert os.path.exists(ifu_driver_template)


@pytest.fixture(scope='session')
def pkg_ifu_driver_template():

    pkg_file_path = str(pathlib.Path(workflow.__path__[0]) / 'stage1' / 'aux' /
                        'ifu_driver_template.fits')
    
    return pkg_file_path


def test_pkg_ifu_driver_template(pkg_ifu_driver_template):

    assert os.path.exists(pkg_ifu_driver_template)


def test_compare_ifu_driver_template(ifu_driver_template,
                                     pkg_ifu_driver_template):

    returncode = subprocess.call(
                     ['fitsdiff', ifu_driver_template, pkg_ifu_driver_template])
    
    assert returncode == 0


@pytest.fixture(scope='session')
def ifu_driver_cat(ifu_driver_template, tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('output').join(
                        'ifu_driver_cat.fits'))

    data_dict = workflow.stage1._get_data_dict_for_example()

    trimester, author, report_verbosity, cc_report = \
        workflow.stage1._set_keywords_info_for_example()
    
    workflow.stage1.create_ifu_driver_cat(ifu_driver_template, data_dict,
                                          file_path, trimester, author,
                                          report_verbosity=report_verbosity,
                                          cc_report=cc_report)
    
    return file_path


def test_create_ifu_driver_cat(ifu_driver_cat):

    assert os.path.exists(ifu_driver_cat)


@pytest.fixture(scope='session')
def pkg_ifu_driver_cat():

    pkg_file_path = str(pathlib.Path(workflow.__path__[0]) / 'stage2' /
                        'input' / 'WC_2020A1-ifu_driver_cat.fits')
    
    return pkg_file_path


def test_pkg_ifu_driver_cat(pkg_ifu_driver_cat):

    assert os.path.exists(pkg_ifu_driver_cat)


def test_compare_ifu_driver_cat(ifu_driver_cat, pkg_ifu_driver_cat):

    returncode = subprocess.call(
                     ['fitsdiff', '-k', 'DATETIME',
                      ifu_driver_cat, pkg_ifu_driver_cat])
    
    assert returncode == 0


@pytest.fixture(scope='session')
def ifu_driver_cat_cheating(ifu_driver_template, tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('output').join(
                        'ifu_driver_cat-cheating.fits'))
    
    xml_files_pattern = str(pathlib.Path(workflow.__path__[0]) / 'stage4' /
                            'input' / '*.xml')

    data_dict = workflow.stage1._get_data_dict_for_cheating(xml_files_pattern)

    trimester, author, report_verbosity, cc_report = \
        workflow.stage1._set_keywords_info_for_cheating()
    
    workflow.stage1.create_ifu_driver_cat(ifu_driver_template, data_dict,
                                          file_path, trimester, author,
                                          report_verbosity=report_verbosity,
                                          cc_report=cc_report)
    
    return file_path


def test_create_ifu_driver_cat_cheating(ifu_driver_cat_cheating):

    assert os.path.exists(ifu_driver_cat_cheating)


def test_compare_ifu_driver_cat_cheating(ifu_driver_cat_cheating,
                                         pkg_ifu_driver_cat):

    returncode = subprocess.call(
                     ['fitsdiff', '-k', 'DATETIME',
                      ifu_driver_cat_cheating, pkg_ifu_driver_cat])
    
    assert returncode == 0

