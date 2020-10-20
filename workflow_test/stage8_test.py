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
import shutil
import subprocess

import pytest

import workflow


@pytest.fixture(scope='module')
def wasp_cat(pkg_combo_cat, tmpdir_factory):

    file_path = str(tmpdir_factory.mktemp('output').join(
                        os.path.basename(pkg_combo_cat)))

    # This will not be really a WASP file: we copy it from the previous stage,
    # and we will skip the comparison of the CNAME column in
    # test_fitsdiff_wasp_cat
    
    shutil.copyfile(pkg_combo_cat, file_path)
    
    return file_path


def test_fitscheck_wasp_cat(wasp_cat):

    returncode = subprocess.call(['fitscheck', wasp_cat])
    
    assert returncode == 0


def test_fitsdiff_wasp_cat(wasp_cat, pkg_wasp_cat):

    returncode = subprocess.call(
                     ['fitsdiff', '-k', 'CHECKSUM,DATASUM,DATETIME',
                      '-f', 'CNAME', wasp_cat, pkg_wasp_cat])
    
    assert returncode == 0

