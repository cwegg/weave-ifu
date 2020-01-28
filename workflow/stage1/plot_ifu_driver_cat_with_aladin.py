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
import os.path
import logging
import urllib.request
import subprocess

from astropy.io import fits


def _get_aladin_jar(aladin_jar_path='Aladin.jar',
                  aladin_jar_url='https://aladin.u-strasbg.fr/java/Aladin.jar'):

    urllib.request.urlretrieve(aladin_jar_url, aladin_jar_path)

    return aladin_jar_path


def _send_command(cmd, p, option='write', encoding='utf-8'):

    assert option in ['write', 'communicate']

    if option == 'write':
        send_func = p.stdin.write
        send_func_str = 'write'
    elif option == 'communicate':
        send_func = p.communicate
        send_func_str = 'communicate'

    logging.debug('{}({})'.format(send_func_str, cmd))

    send_func('{}\n'.format(cmd).encode(encoding))


def _get_obsmode_from_progtemp(progtemp):

    first_char = progtemp[0]

    if first_char in ['1', '2', '3']:
        result = 'MOS'
    elif first_char in ['4', '5', '6']:
        result = 'LIFU'
    elif first_char in ['7', '8', '9']:
        result = 'mIFU'
    else:
        raise ValueError

    return result


def plot_ifu_driver_cat_with_aladin(cat_filename, output_dir='output/',
                                    aladin_jar='Aladin.jar'):
    """
    Plot targets contained in an IFU driver catalogue with Aladin.

    Parameters
    ----------
    cat_filename : str
        A FITS file with an IFU driver catalogue.
    output_dir : str, optional
        The directory which will contain the plots generated with Aladin.
    aladin_jar : str, optional
        The location of the Java JAR file of Aladin.
    """

    logging.info(
        'This plotting tool needs time. If you want to monitor its progress, '
        'see the directory which will contain its ouput images.')

    # Open Aladin
    
    aladin_starting_cmd = 'java -jar {} -nogui'.format(aladin_jar)

    logging.debug('Popen({})'.format(aladin_starting_cmd))
    p = subprocess.Popen([aladin_starting_cmd], shell=True,
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    # Set the grid one

    _send_command('grid on', p)

    # # Load the catalogue

    # _send_command('load {}'.format(cat_filename), p)

    # Open the catalogue

    with fits.open(cat_filename) as hdu_list:

        # Get the data from the catalogue

        data = hdu_list[1].data

        # Set the format string to build the numbers of the images
        
        num_digits = len(str(len(data)))
        num_digits_fmt = '{:0' + str(num_digits) + '}'

        # For each target in the catalogue

        for i in range(len(data)):

            # Get the needed information of the target

            ra = data['GAIA_RA'][i]
            dec = data['GAIA_DEC'][i]
            progtemp = data['PROGTEMP'][i]
            targname = data['TARGNAME'][i]
            targid = data['TARGID'][i]

            # Set a filename for the output image

            cat_basename_wo_ext = os.path.splitext(
                                      os.path.basename(cat_filename))[0]

            img_filename = (output_dir + cat_basename_wo_ext + '-' +
                            num_digits_fmt.format(i + 1) + '-' +
                            targname + '-' + targid + '.png')

            # Create a the string with the coordinates which will be used in the
            # Aladin commands

            coord_str = '{:.10f} {:.10f}'.format(ra, dec)

            # Guess the obsmode and set some parameters acording to it

            obsmode = _get_obsmode_from_progtemp(progtemp)

            if obsmode == 'LIFU':
                zoom_size_str = '30arcmin'
                get_radius_str = '1deg'
            elif obsmode == 'mIFU':
                zoom_size_str = '60arcsec'
                get_radius_str = '2arcmin'
            else:
                logging.warning(
                    'Skipping plotting of target {} of {}'.format(
                        i + 1, cat_filename))

                continue

            # Create a list of commands to be submitted to Aladin

            cmd_list = []

            # Reset and load the catalogue

            cmd_list.append('reset')
            cmd_list.append('load {}'.format(os.path.abspath(cat_filename)))

            # Set the coordinates

            cmd_list.append(coord_str)

            # Get the image

            cmd_list.append('get aladin {} {}'.format(coord_str,
                                                      get_radius_str))

            # Draw some circles to show the field of view of the instrument:
            # - For LIFU, draw circles to show the central bundle and the sky
            #   bundles.
            # - For mIFU, draw a circle for each bundle

            if obsmode == 'LIFU':
                cmd_list.append('draw yellow circle({} 1.51arcmin)'.format(
                    coord_str))
                cmd_list.append('draw green circle({} 8.6arcmin)'.format(
                    coord_str))
                cmd_list.append('draw green circle({} 8.28arcmin)'.format(
                    coord_str))
            elif obsmode == 'mIFU':
                cmd_list.append('draw red circle({} 8.7arcsec)'.format(
                    coord_str))

            # Set the desired zoom

            cmd_list.append('zoom {}'.format(zoom_size_str))

            # Save the image

            cmd_list.append('save {}'.format(os.path.abspath(img_filename)))

            # Join the commands for this target and send them to Aladin

            cmd = '; '.join(cmd_list)
            _send_command(cmd, p)

        # Ask Aladin to quit and wait

        _send_command('quit', p, option='communicate')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
        """
        Plot targets contained in an IFU driver catalogue with Aladin.
        """)

    parser.add_argument('catalogue',
                        help='A FITS file with an IFU driver catalogue')

    parser.add_argument('--dir', default='output/', help=
                        """
                        The directory which will contain the plots generated
                        with Aladin
                        """)

    parser.add_argument('--aladin', default='aux/Aladin.jar',
                        help='The location of the Java JAR file of Aladin')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='The level for the logging messages')

    args = parser.parse_args()
    
    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}
    
    logging.basicConfig(level=level_dict[args.log_level])

    if not os.path.exists(args.aladin):
        logging.info('Downloading the Java JAR file of Aladin')
        _get_aladin_jar(aladin_jar_path=args.aladin)
    
    plot_ifu_driver_cat_with_aladin(args.catalogue, output_dir=args.dir,
                                    aladin_jar=args.aladin)

