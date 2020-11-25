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


import pytest

import workflow


def test_plot_data_from_xml(pkg_wasp_xml_files, tmpdir_factory):

    output_dir = str(tmpdir_factory.mktemp('output'))
    
    output_file_list = []
    
    for xml_file in pkg_wasp_xml_files:
    
        output_filename = workflow.stage10.plot_data_from_xml(
            xml_file, output_dir=output_dir, sim=False, verbose=False)
        
        output_file_list.append(output_filename)
    
    return output_file_list


def test_plot_data_from_xml_verbose(pkg_wasp_xml_files, tmpdir_factory):

    output_dir = str(tmpdir_factory.mktemp('output'))
    
    output_file_list = []
    
    for xml_file in pkg_wasp_xml_files:
    
        output_filename = workflow.stage10.plot_data_from_xml(
            xml_file, output_dir=output_dir, sim=False, verbose=True)
        
        output_file_list.append(output_filename)
    
    return output_file_list

