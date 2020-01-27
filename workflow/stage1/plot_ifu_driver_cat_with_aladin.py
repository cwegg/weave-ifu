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


def get_aladin_jar(aladin_jar_path='Aladin.jar',
                   aladin_jar_url='https://aladin.u-strasbg.fr/java/Aladin.jar'):

    urllib.request.urlretrieve(aladin_jar_url, aladin_jar_path)

    return aladin_jar_path


def _write_command(cmd, p):

    cmd = cmd + '\n'

    logging.debug('p.stdin.write({})'.format(cmd))

    p.stdin.write(cmd.encode())


def plot_ifu_driver_cat_with_aladin(filename, output_dir='img/',
                                    aladin_jar='Aladin.jar'):

    aladin_starting_cmd = ['aladin'.format(aladin_jar)]
    aladin_starting_cmd = ['java -jar {}'.format(aladin_jar)]
    aladin_starting_cmd = ['java -jar {} -nogui'.format(aladin_jar)]

    logging.debug(aladin_starting_cmd)
    p = subprocess.Popen(aladin_starting_cmd, shell=True, stdin=subprocess.PIPE)

    # p.wait()

    _write_command('grid on', p)

    _write_command('load {}'.format(filename), p)
    
    with fits.open(filename) as hdu_list:
    
        data = hdu_list[1].data

        for i in range(len(data)):

            img_filename = output_dir + '{}.png'.format(i + 1)

            ra = data['GAIA_RA'][i]
            dec = data['GAIA_DEC'][i]

            # _write_command('reset', p)
            # _write_command('grid on', p)
            # _write_command('load {}'.format(filename), p)
            
            _write_command('get aladin {} {}'.format(ra, dec), p)
            _write_command('zoom 10arcmin', p)
            _write_command('save {}'.format(img_filename), p)

            break

        _write_command('quit', p)

    p.wait()


if __name__ == '__main__':

    filename = './output/WC_IFU.fits'
    
    output_dir = 'img/'

    aladin_jar_dir = 'aux/'
    aladin_jar_filename = 'Aladin.jar'

    if aladin_jar_dir[-1] != os.path.sep:
        aladin_jar_dir = aladin_jar_dir + os.path.sep

    aladin_jar_path = aladin_jar_dir + aladin_jar_filename

    if not os.path.exists(aladin_jar_path):
        get_aladin_jar(aladin_jar_path=aladin_jar_path)
    
    plot_ifu_driver_cat_with_aladin(filename, output_dir=output_dir,
                                    aladin_jar=aladin_jar_path)

