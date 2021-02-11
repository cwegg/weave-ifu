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


# These tests have not been implemented, as they depend on the configure tools

#@pytest.fixture(scope='module')
#def tgcs_xml_files(pkg_tgc_xml_files, tmpdir_factory):

#    output_dir = str(tmpdir_factory.mktemp('output'))

#    raise NotImplementedError
#    
#    return xml_filename_list


#def test_diff_tgcs_xml_files(tgcs_xml_files, pkg_tgcs_xml_files):
#    
#    assert len(tgcs_xml_files) == len(pkg_tgcs_xml_files)
#    
#    for ref_file, copy_file in zip(tgcs_xml_files, pkg_tgcs_xml_files):
#    
#        assert os.path.basename(ref_file) == os.path.basename(copy_file)

#        returncode = subprocess.call(['diff', '-q', ref_file, copy_file])
#        
#        assert returncode == 0

