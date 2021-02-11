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
import subprocess

import pytest

import weave_ifu_workflow


@pytest.fixture(scope='module')
def ifu_driver_template(master_cat, tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('aux').join(
                        'ifu_driver_template.fits'))
    
    weave_ifu_workflow.stage1.create_ifu_driver_template(master_cat, file_path)

    assert os.path.exists(file_path)
    
    return file_path


def test_fitscheck_ifu_driver_template(ifu_driver_template):

    returncode = subprocess.call(['fitscheck', ifu_driver_template])
    
    assert returncode == 0


def test_fitsdiff_ifu_driver_template(ifu_driver_template,
                                      pkg_ifu_driver_template):

    returncode = subprocess.call(
                     ['fitsdiff', '-k', 'CHECKSUM,DATASUM',
                      ifu_driver_template, pkg_ifu_driver_template])
    
    assert returncode == 0


@pytest.fixture(scope='module')
def ifu_driver_cat(ifu_driver_template, tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('output').join(
                        'ifu_driver_cat.fits'))

    data_dict = weave_ifu_workflow.stage1._get_data_dict_for_example()

    trimester, author, report_verbosity, cc_report = \
        weave_ifu_workflow.stage1._set_keywords_info_for_example()
    
    weave_ifu_workflow.stage1.create_ifu_driver_cat(ifu_driver_template, data_dict,
                                                    file_path, trimester, author,
                                                    report_verbosity=report_verbosity,
                                                    cc_report=cc_report)
    
    assert os.path.exists(file_path)
    
    return file_path


def test_fitscheck_ifu_driver_cat(ifu_driver_cat):

    returncode = subprocess.call(['fitscheck', ifu_driver_cat])
    
    assert returncode == 0


def test_fitsdiff_ifu_driver_cat(ifu_driver_cat, pkg_ifu_driver_cat):

    returncode = subprocess.call(
                     ['fitsdiff', '-k', 'CHECKSUM,DATASUM,DATETIME',
                      ifu_driver_cat, pkg_ifu_driver_cat])
    
    assert returncode == 0


@pytest.fixture(scope='module')
def ifu_driver_cat_cheating(ifu_driver_template, pkg_tgc_xml_files,
                            tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('output').join(
                        'ifu_driver_cat-cheating.fits'))

    data_dict = weave_ifu_workflow.stage1._get_data_dict_for_cheating(pkg_tgc_xml_files)

    trimester, author, report_verbosity, cc_report = \
        weave_ifu_workflow.stage1._get_keywords_info_for_cheating(pkg_tgc_xml_files)
    
    weave_ifu_workflow.stage1.create_ifu_driver_cat(ifu_driver_template, data_dict,
                                                    file_path, trimester, author,
                                                    report_verbosity=report_verbosity,
                                                    cc_report=cc_report)
    
    assert os.path.exists(file_path)
    
    return file_path


def test_fitscheck_ifu_driver_cat_cheating(ifu_driver_cat_cheating):

    returncode = subprocess.call(['fitscheck', ifu_driver_cat_cheating])
    
    assert returncode == 0


def test_fitsdiff_ifu_driver_cat_cheating(ifu_driver_cat_cheating,
                                          pkg_ifu_driver_cat):

    returncode = subprocess.call(
                     ['fitsdiff', '-k', 'CHECKSUM,DATASUM,DATETIME',
                      ifu_driver_cat_cheating, pkg_ifu_driver_cat])
    
    assert returncode == 0


def test_check_ifu_driver_cat(ifu_driver_cat, ifu_driver_template):
    
    assert weave_ifu_workflow.stage1.check_ifu_driver_cat(ifu_driver_cat,
                                                          template=ifu_driver_template,
                                                          check_vs_template=True)

