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


import shutil
import subprocess

import pytest

import workflow


@pytest.fixture(scope='module')
def ifu_cat(pkg_ifu_cat_from_xmls, tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('output').join(
                        'WC_2020A1-ifu.fits'))

    # This line should be replace by the functions developed for stage 6, once
    # they will have been developed. They should be avaialable as:
    # workflow.stage6.function_1
    # workflow.stage6.function_2
    # ...
    # workflow.stage6.function_n
    
    shutil.copyfile(pkg_ifu_cat_from_xmls, file_path)
    
    return file_path


def test_fitscheck_ifu_cat(ifu_cat):

    returncode = subprocess.call(['fitscheck', ifu_cat])
    
    assert returncode == 0


def test_fitsdiff_ifu_cat(ifu_cat, pkg_ifu_cat):

    returncode = subprocess.call(
                     ['fitsdiff', '-k', 'CHECKSUM,DATASUM,DATETIME',
                      ifu_cat, pkg_ifu_cat])
    
    assert returncode == 0

