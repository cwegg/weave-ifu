#!/usr/bin/env python3

#
# Copyright (C) 2019 Cambridge Astronomical Survey Unit
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

from astropy.io import fits


def plot_ifu_driver_cat_with_aladin(filename, output_dir='./img/',
                                    aladin_jar='/usr/local/aladin/Aladin.jar'):
    
#    p = subprocess.Popen(['aladin'.format(aladin_jar)], shell=True,
#                         stdin=subprocess.PIPE)
#    p = subprocess.Popen(['java -jar {}'.format(aladin_jar)], shell=True,
#                         stdin=subprocess.PIPE)
# p = subprocess.Popen(['java -jar {} -nogui'.format(aladin_jar)], shell=True,

#    p.wait()
    print('grid on\n')
#    p.stdin.write('grid on\n')
    print('load {}\n'.format(filename))
#    p.stdin.write('load {}\n'.format(filename))
    
    with fits.open(filename) as hdu_list:
    
        data = hdu_list[1].data

        for i in range(len(data)):
    #        p.stdin.write('reset\n')
            ra = data['GAIA_RA'][i]
            dec = data['GAIA_DEC'][i]
            
            img_filename = output_dir + '{}.png'.format(i + 1)
            
            print('get aladin {} {}\n'.format(ra, dec))
#            p.stdin.write('get aladin {} {}\n'.format(ra, dec))
            print('zoom 10arcmin\n')
#            p.stdin.write('zoom 10arcmin\n')
            print('save {}\n'.format(img_filename))
#            p.stdin.write('save {}\n'.format(img_filename))

#        p.stdin.write('quit\n')

#    p.wait()


if __name__ == '__main__':

    filename = './output/WC_IFU.fits'
    
    output_dir = './img/'
    
    aladin_jar='/usr/local/aladin/Aladin.jar'
    
    plot_ifu_driver_cat_with_aladin(filename, output_dir=output_dir,
                                    aladin_jar=aladin_jar)

