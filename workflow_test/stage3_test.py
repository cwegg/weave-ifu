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

import workflow


@pytest.fixture(scope='module')
def tgc_xml_files(pkg_t_xml_files, tmpdir_factory):

    output_dir = str(tmpdir_factory.mktemp('output'))

    xml_filename_list = workflow.stage3.add_guide_and_calib_stars(
        pkg_t_xml_files, output_dir, mifu_num_guide_stars_request=None)
    
    return xml_filename_list


def test_diff_tgc_xml_files(tgc_xml_files, pkg_tgc_xml_files):
    
    assert len(tgc_xml_files) == len(pkg_tgc_xml_files)
    
    for ref_file, copy_file in zip(tgc_xml_files, pkg_tgc_xml_files):
    
        assert os.path.basename(ref_file) == os.path.basename(copy_file)

        returncode = subprocess.call(['diff', '-q', ref_file, copy_file])
        
        assert returncode == 0

