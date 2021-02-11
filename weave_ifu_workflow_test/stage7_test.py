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


import subprocess

import pytest

import weave_ifu_workflow


@pytest.fixture(scope='module')
def combo_cat(pkg_mos_cat, pkg_ifu_cat, tmpdir_factory):

    output_dir = str(tmpdir_factory.mktemp('output'))

    file_path = weave_ifu_workflow.stage7.create_combo_fits_cat(
        pkg_mos_cat, pkg_ifu_cat, output_dir=output_dir)
    
    return file_path


def test_fitscheck_combo_cat(combo_cat):

    returncode = subprocess.call(['fitscheck', combo_cat])
    
    assert returncode == 0


def test_fitsdiff_combo_cat(combo_cat, pkg_combo_cat):

    returncode = subprocess.call(
                     ['fitsdiff', '-k', 'CHECKSUM,DATASUM,DATETIME',
                      combo_cat, pkg_combo_cat])
    
    assert returncode == 0

