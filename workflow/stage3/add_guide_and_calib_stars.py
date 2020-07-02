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


import argparse
import logging
import os

from workflow.utils.classes import OBXML


def add_guide_and_calib_stars(xml_file_list, output_dir, mifu_num_calibs=2,
                              all_guide_stars=False, all_calib_stars=False,
                              overwrite=False):
    """
    Add guide and calib stars to XML files.
    
    Parameters
    ----------
    xml_file_list : list of str
        A list of input OB XML files.
    output_dir : str
        Name of the directory which will containe the output XML files.
    mifu_num_calibs : int, optional
        Number of mIFU calibration stars to be added (it will be limited by
        max_guides attribute too).
    all_guide_stars : bool, optional
        Add all available guide stars to the XML file.
    all_calib_stars : bool, optional
        Add all available calib stars to the XML file.
    overwrite : bool, optional
        Overwrite the output FITS file.

    Returns
    -------
    output_file_list : list of str
        A list with the output XML files.
    """

    output_file_list = []

    xml_file_list.sort()

    for xml_file in xml_file_list:

        # Check that the input XML exists and is a file

        assert os.path.isfile(xml_file)

        # Choose the output filename depedending on the input filename

        input_basename_wo_ext = os.path.splitext(os.path.basename(xml_file))[0]

        if (input_basename_wo_ext.endswith('-t') or
            input_basename_wo_ext.endswith('-')):
            output_basename_wo_ext = input_basename_wo_ext + 'gc'
        else:
            output_basename_wo_ext = input_basename_wo_ext + '-gc'

        output_file = os.path.join(output_dir, output_basename_wo_ext + '.xml')

        # Save the output filename for the result

        output_file_list.append(output_file)

        # If the output file already exists, delete it or continue with the next
        # one

        if os.path.exists(output_file):
            if overwrite == True:
                logging.info('Removing previous file: {}'.format(output_file))
                os.remove(output_file)
            else:
                logging.info(
                    'Skipping file {} as its output already exists: {}'.format(
                        xml_file, output_file))
                continue

        # Read the input file, add the guide and calib stars and write it to the
        # output file
        
        ob_xml = OBXML(xml_file)

        ob_xml.add_guide_and_calib_stars(mifu_num_calibs=mifu_num_calibs,
                                         all_guide_stars=all_guide_stars,
                                         all_calib_stars=all_calib_stars)

        ob_xml.write_xml(output_file)

    return output_file_list


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
             description='Add guide and calib stars to XML files')

    parser.add_argument('xml_file', nargs='+',
                        help="""an input OB XML file""")

    parser.add_argument('--mifu_num_calibs', dest='mifu_num_calibs', default=2,
                        choices=range(20), type=int,
                        help="""number of mIFU calibration stars to be added
                        (it will be limited by max_guides attribute too)""")

    parser.add_argument('--all_guide_stars', dest='all_guide_stars',
                        action='store_true',
                        help='add all available guide stars to the XML file')

    parser.add_argument('--all_calib_stars', dest='all_calib_stars',
                        action='store_true',
                        help='add all available calib stars to the XML file')

    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name of the directory which will contain the
                        output XML files""")

    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='overwrite the output files')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}

    logging.basicConfig(level=level_dict[args.log_level])

    if not os.path.exists(args.output_dir):
        logging.info('Creating the output directory')
        os.mkdir(args.output_dir)

    add_guide_and_calib_stars(args.xml_file, args.output_dir,
                              mifu_num_calibs=args.mifu_num_calibs,
                              all_guide_stars=args.all_guide_stars,
                              all_calib_stars=args.all_calib_stars,
                              overwrite=args.overwrite)
 
