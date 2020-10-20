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


def test_aladin_jar(aladin_jar):

    # We simply pass, as the existence of the downloaded file is checked in the
    # fixture
    pass


def test_fitscheck_master_cat(master_cat):

    returncode = subprocess.call(['fitscheck', master_cat])
    
    assert returncode == 0


def test_blank_xml_template(blank_xml_template):

    # We simply pass, as the existence of the downloaded file is checked in the
    # fixture
    pass


def test_progtemp_file(progtemp_file):

    # We simply pass, as the existence of the downloaded file is checked in the
    # fixture
    pass


def test_obstemp_file(obstemp_file):

    # We simply pass, as the existence of the downloaded file is checked in the
    # fixture
    pass

