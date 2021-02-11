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

from astropy.io import fits
        
    
def create_empty_mos_cat(catalogue_template, output_filename, overwrite=False):

    with fits.open(catalogue_template) as hdu_list:
        
        hdu_list[0].header['TRIMESTE'] = '2020A1'
        hdu_list[0].header['CAT_NME1'] = 'First Name'
        hdu_list[0].header['CAT_NME2'] = 'Surname'
        hdu_list[0].header['CAT_MAIL'] = 'a@domain.com'
        hdu_list[0].header['CAT_CC'] = 'b@domain.com,c@domain.com'

        hdu_list.writeto(output_filename, checksum=True, overwrite=overwrite)


if __name__ == '__main__':

    catalogue_template = os.path.join('..', 'weave_ifu_workflow', 'stage5', 'aux',
                                      'WC_CatalogueTemplate.fits')

    output_filename = os.path.join('..', 'weave_ifu_workflow', 'stage7', 'input',
                                   'WC_2020A1-mos.fits')

    create_empty_mos_cat(catalogue_template, output_filename)

