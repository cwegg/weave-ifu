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

import workflow


@pytest.fixture(scope='module')
def ifu_cat_from_xmls(pkg_tgcs_xml_files, pkg_wc_cat, tmpdir_factory):

    output_dir = str(tmpdir_factory.mktemp('output'))

    file_path = workflow.stage5.create_ifu_fits_cat(
        pkg_wc_cat, pkg_tgcs_xml_files,
        cat_nme1='First Name', cat_nme2='Surname', output_dir=output_dir)
    
    return file_path


def test_fitscheck_ifu_cat_from_xmls(ifu_cat_from_xmls):

    returncode = subprocess.call(['fitscheck', ifu_cat_from_xmls])
    
    assert returncode == 0


def test_fitsdiff_ifu_cat_from_xmls(ifu_cat_from_xmls, pkg_ifu_cat_from_xmls):

    returncode = subprocess.call(
                     ['fitsdiff', '-k', 'CHECKSUM,DATASUM,DATETIME',
                      ifu_cat_from_xmls, pkg_ifu_cat_from_xmls])
    
    assert returncode == 0

