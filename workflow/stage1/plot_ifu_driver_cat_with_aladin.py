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


def plot_ifu_driver_cat_with_aladin(filename, output_dir='output/',
                                    aladin_jar='Aladin.jar'):

    logging.info(
        """
        This plotting tool need time. If you want to monitor its progress, see
        the directory which will contain its ouput images.
        """)

    aladin_starting_cmd = 'aladin'.format(aladin_jar)
    aladin_starting_cmd = 'java -jar {}'.format(aladin_jar)
    aladin_starting_cmd = 'java -jar {} -nogui'.format(aladin_jar)

    logging.debug('Popen({})'.format(aladin_starting_cmd))
    p = subprocess.Popen([aladin_starting_cmd], shell=True,
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    _send_command('grid on', p)

    # _send_command('load {}'.format(filename), p)
    
    with fits.open(filename) as hdu_list:

        # Get the data

        data = hdu_list[1].data

        # Set the format string for the numbers of the images
        
        num_digits = len(str(len(data)))
        num_digits_fmt = '{:0' + str(num_digits) + '}'

        # For each row of the data

        for i in range(len(data)):

            # Get the needed information

            ra = data['GAIA_RA'][i]
            dec = data['GAIA_DEC'][i]
            progtemp = data['PROGTEMP'][i]
            targname = data['TARGNAME'][i]
            targid = data['TARGID'][i]

            # Set a filename

            img_filename = (output_dir + num_digits_fmt.format(i + 1) +
                            '-' + targname + '-' + targid + '.png')

            # Create a the string with the coordinates which will be used in the
            # Aladin commands

            coord_str = '{} {}'.format(ra, dec)

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
                    'Skipping plotting of row {} of {}'.format(
                        i + 1, filename))

                continue

            # Create a list of commands to be submitted to Aladin

            cmd_list = []

            # Reset and load the catalogue

            cmd_list.append('reset')
            # cmd_list.append('grid on')
            cmd_list.append('load {}'.format(os.path.abspath(filename)))

            # Set the coordinates and the desired zoom

            cmd_list.append(coord_str)

            cmd_list.append('zoom {}'.format(zoom_size_str))

            # Get the image

            # cmd_list.append('get aladin {}'.format(coord_str))
            # cmd_list.append('get aladin {} {}'.format(coord_str,
            #                                           get_radius_str))
            cmd_list.append('get ESO(DSS2/color) {} {}'.format(coord_str,
                                                         get_radius_str))
            # cmd_list.append('get hips(CDS/P/DSS2/color)'.format(coord_str,
            #                                           get_radius_str))

            # Draw some circles to show the field of view of the instrument

            if obsmode == 'LIFU':

                # Draw a circle for the inner LIFU bundle

                cmd_list.append('draw yellow circle({} 1.51arcmin)'.format(
                    coord_str))

                # Draw two circles for the sky LIFU bundles

                cmd_list.append('draw green circle({} 8.6arcmin)'.format(
                    coord_str))
                cmd_list.append('draw green circle({} 8.28arcmin)'.format(
                    coord_str))

            if obsmode == 'mIFU':

                # Draw a circle for a mIFU bundle

                cmd_list.append('draw red circle({} 8.7arcsec)'.format(
                    coord_str))

            # Set the desired zoom

            cmd_list.append('zoom {}'.format(zoom_size_str))

            # Save the image

            cmd_list.append('save {}'.format(os.path.abspath(img_filename)))

            cmd = '; '.join(cmd_list)

            _send_command(cmd, p)

        # Ask Aladin to quit and wait

        _send_command('quit', p, option='communicate')


if __name__ == '__main__':

    # logging.basicConfig(level=logging.DEBUG)

    filename = 'output/WC_IFU.fits'
    
    output_dir = 'output/'

    aladin_jar_dir = 'aux/'
    aladin_jar_filename = 'Aladin.jar'

    if aladin_jar_dir[-1] != os.path.sep:
        aladin_jar_dir = aladin_jar_dir + os.path.sep

    aladin_jar_path = aladin_jar_dir + aladin_jar_filename

    if not os.path.exists(aladin_jar_path):
        _get_aladin_jar(aladin_jar_path=aladin_jar_path)
    
    plot_ifu_driver_cat_with_aladin(filename, output_dir=output_dir,
                                    aladin_jar=aladin_jar_path)

