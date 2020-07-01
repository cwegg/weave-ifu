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


import glob
import os

from workflow.utils.classes import OBXML


if __name__ == '__main__':

    
    input_xml_loc = './input/'
    
    try:
        assert os.path.isdir(input_xml_loc)
    except:
        print('Please supply valid XML location')
        print('Usage: ./stage3.py location')
        raise SystemExit(0)

    output_dir = '../stage3/output/'

    if 1:
        cmd = 'rm -Rf %s/*.xml'%(output_dir)
        os.system(cmd)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print('Reading in XMLs from %s'%(input_xml_loc))

    filename_list = glob.glob(input_xml_loc+'/*.xml')
    filename_list.sort()

    if len(filename_list) == 0:
        raise SystemExit('No XML files found in %s'%(input_loc))

    for filename in filename_list:
        output_file = output_dir+os.path.basename(filename)
        
        stage3_ifu = OBXML(filename)
        stage3_ifu.add_guide_and_calib_stars()
        stage3_ifu.write_xml(output_file)
    
    print('IFU XMLs written to: {0}'.format(output_dir))

