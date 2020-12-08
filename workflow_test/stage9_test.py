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
def wasp_xml_files(pkg_wasp_cat, pkg_copy_tgcs_xml_files, tmpdir_factory):

    output_dir = str(tmpdir_factory.mktemp('output'))

    xml_filename_list = workflow.stage9.fill_xmls_with_fits_info(
        pkg_wasp_cat, pkg_copy_tgcs_xml_files, output_dir)
    
    # Modify the output files
    
    for xml_filename in xml_filename_list:
    
        with open(xml_filename, 'r') as xml_file:
            contents = xml_file.read()
        
        contents = contents.replace('%%%', '')
        
        with open(xml_filename, 'w') as xml_file:
            xml_file.write(contents)
    
    return xml_filename_list


def test_diff_t_xml_files(wasp_xml_files, pkg_wasp_xml_files):
    
    assert len(wasp_xml_files) == len(pkg_wasp_xml_files)
    
    for ref_file, copy_file in zip(wasp_xml_files, pkg_wasp_xml_files):
    
        assert os.path.basename(ref_file) == os.path.basename(copy_file)

        returncode = subprocess.call(['diff', '-q', ref_file, copy_file])
        
        assert returncode == 0

