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

from astropy.io import fits

from workflow.utils.classes import OBXML


def fill_xmls_with_fits_info(fits_cat, xml_files, output_dir, add_sim=False,
                             overwrite=False):
    """
    Combine MOS and IFU catalogues to create a combo catalogue.
    
    Parameters
    ----------
    fits_cat : str
        Name of a FITS file containing a combo FITS catalogue.
    xml_files : list of str
        A list of strings with the filenames of the XML files.
    output_dir : str
        Name of the directory which will contain the output XML files.
    add_sim : bool, optional
        Add simulation elements to the targets.
    overwrite : bool, optional
        Overwrite the output XML files.
    """
        
    # Open the combo FITS catalogue and read the needed data

    with fits.open(fits_cat) as hdu_list:
        
        fits_data = hdu_list[1].data
        
        if add_sim == True:
            sim_data = hdu_list[2].data
        else:
            sim_data = None

    # For each XML file

    for xml_file in xml_files:

        # Create an OB from the XML file
        
        ob_xml = OBXML(xml_file)
        
        # Choose an output filename from the input filename
        
        input_basename = os.path.basename(xml_file)

        output_basename_wo_ext = '-'.join(input_basename.split('-')[:-1])
        
        output_basename = '{}.xml'.format(output_basename_wo_ext)
        
        output_path = os.path.join(output_dir, output_basename)
        
        # If the output file exists, remove or skip it

        if os.path.exists(output_path):
            
            if overwrite == True:
                logging.info('Removing previous file: {}'.format(output_path))
                os.remove(output_path)
            else:
                logging.info(
                    'Skipping file because it already exists: {}'.format(
                        output_path))
                continue

        # Populate the OB XML with the FITS data

        ob_xml.populate_targets_with_fits_data(fits_data, sim_data=sim_data)

        # Write the OB XML to a file
        
        ob_xml.write_xml(output_path)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Fill XMLs with FITS info')

    parser.add_argument('fits_cat',
                        help='a FITS file containing a combo FITS catalogue')

    parser.add_argument('xml_file', nargs='+', help='name of a XML file')

    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name of the directory which will contain the
                        output XML files""")

    parser.add_argument('--add_sim', action='store_true',
                        help='add simulation elements to the targets')

    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite the output files')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.output_dir):
        logging.info('Creating the output directory')
        os.mkdir(args.output_dir)

    fill_xmls_with_fits_info(args.fits_cat, args.xml_file, args.output_dir,
                             add_sim=args.add_sim, overwrite=args.overwrite)

