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

from astropy.io import fits


def test_fitscheck_pkg_wc_cat(pkg_wc_cat):

    returncode = subprocess.call(['fitscheck', pkg_wc_cat])
    
    assert returncode == 0


def test_datamver_pkg_wc_cat(pkg_wc_cat, master_cat):

    with fits.open(pkg_wc_cat) as hdu_list:
        pkg_wc_cat_datamver = hdu_list[0].header['DATAMVER']

    with fits.open(master_cat) as hdu_list:
        master_cat_datamver = hdu_list[0].header['DATAMVER']
    
    assert pkg_wc_cat_datamver == master_cat_datamver


def test_fitscheck_pkg_ifu_driver_template(pkg_ifu_driver_template):

    returncode = subprocess.call(['fitscheck', pkg_ifu_driver_template])
    
    assert returncode == 0


def test_fitscheck_pkg_ifu_driver_cat(pkg_ifu_driver_cat):

    returncode = subprocess.call(['fitscheck', pkg_ifu_driver_cat])
    
    assert returncode == 0


def test_pkg_t_xml_files(pkg_t_xml_files):

    # We simply pass, as the existence of the files is checked in the fixture
    pass


def test_pkg_tgc_xml_files(pkg_tgc_xml_files):

    # We simply pass, as the existence of the files is checked in the fixture
    pass


def test_pkg_tgcs_xml_files(pkg_tgcs_xml_files):

    # We simply pass, as the existence of the files is checked in the fixture
    pass


def test_fitscheck_pkg_ifu_cat_from_xmls(pkg_ifu_cat_from_xmls):

    returncode = subprocess.call(['fitscheck', pkg_ifu_cat_from_xmls])
    
    assert returncode == 0


def test_fitscheck_pkg_ifu_cat(pkg_ifu_cat):

    returncode = subprocess.call(['fitscheck', pkg_ifu_cat])
    
    assert returncode == 0


def test_fitscheck_pkg_mos_cat(pkg_mos_cat):

    returncode = subprocess.call(['fitscheck', pkg_mos_cat])
    
    assert returncode == 0


def test_datamver_pkg_mos_cat(pkg_mos_cat, master_cat):

    with fits.open(pkg_mos_cat) as hdu_list:
        pkg_mos_cat_datamver = hdu_list[0].header['DATAMVER']

    with fits.open(master_cat) as hdu_list:
        master_cat_datamver = hdu_list[0].header['DATAMVER']
    
    assert pkg_mos_cat_datamver == master_cat_datamver


def test_fitscheck_pkg_combo_cat(pkg_combo_cat):

    returncode = subprocess.call(['fitscheck', pkg_combo_cat])
    
    assert returncode == 0


def test_fitscheck_pkg_wasp_cat(pkg_wasp_cat):

    returncode = subprocess.call(['fitscheck', pkg_wasp_cat])
    
    assert returncode == 0


def test_pkg_copy_tgcs_xml_files(pkg_copy_tgcs_xml_files):

    # We simply pass, as the existence of the files is checked in the fixture
    pass


def test_diff_pkg_copy_tgcs_xml_files(pkg_tgcs_xml_files,
                                      pkg_copy_tgcs_xml_files):
    
    assert len(pkg_tgcs_xml_files) == len(pkg_copy_tgcs_xml_files)
    
    for ref_file, copy_file in zip(pkg_tgcs_xml_files, pkg_copy_tgcs_xml_files):
    
        assert os.path.basename(ref_file) == os.path.basename(copy_file)

        returncode = subprocess.call(['diff', '-q', ref_file, copy_file])
        
        assert returncode == 0


def test_pkg_wasp_xml_files(pkg_wasp_xml_files):

    # We simply pass, as the existence of the files is checked in the fixture
    pass

