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


def add_guide_and_calib_stars(
        xml_file_list, output_dir, lifu_num_guide_stars_request=1,
        mifu_num_guide_stars_request=8, mifu_num_central_guide_stars=1,
        mifu_min_guide_cut=0.9, mifu_max_guide_cut=1.0,
        num_calib_stars_request=2, num_central_calib_stars=0,
        min_calib_cut=0.2, max_calib_cut=0.4, overwrite=False):
    """
    Add guide and calib stars to XML files.
    
    Parameters
    ----------
    xml_file_list : list of str
        A list of input OB XML files.
    output_dir : str
        Name of the directory which will containe the output XML files.
    lifu_num_guide_stars_request : int, optional
        Maximum number of LIFU guide stars in the output. None means no limit.
    mifu_num_guide_stars_request : int, optional
        Maximum number of mIFU guide stars in the output. None means no limit.
    mifu_num_central_guide_stars : int, optional
        Number of mIFU guide stars near to centre to be selected.
    mifu_min_guide_cut : float, optional
        Minimum cut factor to be used for the non-central mIFU guide stars.
    mifu_max_guide_cut : float, optional
        Maximum cut factor to be used for the non-central mIFU guide stars.
    num_calib_stars_request : int, optional
        Maximum number of calib stars in the output. None means no limit.
    num_central_calib_stars : int, optional
        Number of calib stars near to centre to be selected.
    min_calib_cut : float, optional
        Minimum cut factor to be used for the non-central calib stars.
    max_calib_cut : float, optional
        Maximum cut factor to be used for the non-central calib stars.
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
        guide_plot_filename = os.path.join(
            output_dir,output_basename_wo_ext + '-guide_stars.png')
        calib_plot_filename = os.path.join(
            output_dir, output_basename_wo_ext + '-calib_stars.png')

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

        ob_xml.add_guide_and_calib_stars(
            guide_plot_filename=guide_plot_filename,
            lifu_num_guide_stars_request=lifu_num_guide_stars_request,
            mifu_num_guide_stars_request=mifu_num_guide_stars_request,
            mifu_num_central_guide_stars=mifu_num_central_guide_stars,
            mifu_min_guide_cut=mifu_min_guide_cut,
            mifu_max_guide_cut=mifu_max_guide_cut,
            calib_plot_filename=calib_plot_filename,
            num_calib_stars_request=num_calib_stars_request,
            num_central_calib_stars=num_central_calib_stars,
            min_calib_cut=min_calib_cut,
            max_calib_cut=max_calib_cut)

        ob_xml.write_xml(output_file)

    return output_file_list


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
             description='Add guide and calib stars to XML files')

    parser.add_argument('xml_file', nargs='+',
                        help="""an input OB XML file""")

    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name of the directory which will contain the
                        output XML files""")

    parser.add_argument('--lifu_num_guide_stars_request', default=1,
                        help="""maximum number of LIFU guide stars in the
                        output; None means no limit""")

    parser.add_argument('--mifu_num_guide_stars_request', default=8,
                        help="""maximum number of mIFU guide stars in the
                        output; None means no limit""")

    parser.add_argument('--mifu_num_central_guide_stars', default=1, type=int,
                        help="""number of mIFU guide stars near to centre to be
                        selected""")

    parser.add_argument('--mifu_min_guide_cut', default=0.9, type=float,
                        help="""minimum cut factor to be used for the
                        non-central mIFU guide stars""")

    parser.add_argument('--mifu_max_guide_cut', default=1.0, type=float,
                        help="""maximum cut factor to be used for the
                        non-central mIFU guide stars""")

    parser.add_argument('--num_calib_stars_request', default=2, type=int,
                        help="""maximum number of calib stars in the output;
                        None means no limit""")

    parser.add_argument('--num_central_calib_stars', default=0, type=int,
                        help="""number of calib stars near to centre to be
                        selected""")

    parser.add_argument('--min_calib_cut', default=0.2, type=float,
                        help="""minimum cut factor to be used for the
                        non-central calib stars""")

    parser.add_argument('--max_calib_cut', default=0.4, type=float,
                        help="""maximum cut factor to be used for the
                        non-central calib stars""")

    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='overwrite the output files')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.output_dir):
        logging.info('Creating the output directory')
        os.mkdir(args.output_dir)
    
    if args.lifu_num_guide_stars_request != 'None':
        lifu_num_guide_stars_request = int(args.lifu_num_guide_stars_request)
    else:
        lifu_num_guide_stars_request = None
    
    if args.mifu_num_guide_stars_request != 'None':
        mifu_num_guide_stars_request = int(args.mifu_num_guide_stars_request)
    else:
        mifu_num_guide_stars_request = None
    
    if args.num_calib_stars_request != 'None':
        num_calib_stars_request = args.num_calib_stars_request
        assert type(num_calib_stars_request) == int
    else:
        num_calib_stars_request = None

    add_guide_and_calib_stars(
        args.xml_file, args.output_dir,
        lifu_num_guide_stars_request=lifu_num_guide_stars_request,
        mifu_num_guide_stars_request=mifu_num_guide_stars_request,
        mifu_num_central_guide_stars=args.mifu_num_central_guide_stars,
        mifu_min_guide_cut=args.mifu_min_guide_cut,
        mifu_max_guide_cut=args.mifu_max_guide_cut,
        num_calib_stars_request=num_calib_stars_request,
        num_central_calib_stars=args.num_central_calib_stars,
        min_calib_cut=args.min_calib_cut,
        max_calib_cut=args.max_calib_cut,
        overwrite=args.overwrite)
 
