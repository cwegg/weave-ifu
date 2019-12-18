#!/usr/bin/env python3

#
# Copyright (C) 2018 Cambridge Astronomical Survey Unit
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


import xml.dom.minidom

from astropy.io import fits


def get_lookup():

    lookup = {}

    lookup[''] = 'TARGSRVY'
    lookup[''] = 'TARGPROG'
    #lookup['target:targcat'] = 'TARGCAT'
    lookup['target:targid'] = 'TARGID'
    lookup['target:targname'] = 'TARGNAME'
    lookup[''] = 'TARGPRIO'
    lookup[''] = 'TARGCLASS'
    lookup['observation:progtemp'] = 'PROGTEMP'
    lookup['obsconstraints:obstemp'] = 'OBSTEMP'
    lookup[''] = 'GAIA_ID'
    lookup[''] = 'GAIA_DR'
    lookup['target:targra'] = 'GAIA_RA'
    lookup['target:targdec'] = 'GAIA_DEC'
    lookup['target:targepoch'] = 'GAIA_EPOCH'
    lookup[''] = 'GAIA_PMRA'
    lookup[''] = 'GAIA_PMDEC'
    lookup[''] = 'GAIA_PARAL'
    # lookup[''] = 'GAIA_PARAL_ERR'
    lookup[''] = 'IFU_SPAXEL'
    lookup['observation:pa'] = 'IFU_PA'
    lookup['dithering:apply_dither'] = 'IFU_DITHER'
    lookup[''] = 'MAG_G'
    lookup[''] = 'MAG_G_ERR'
    lookup[''] = 'MAG_R'
    lookup[''] = 'MAG_R_ERR'
    lookup[''] = 'MAG_I'
    lookup[''] = 'MAG_I_ERR'
    lookup[''] = 'GAIA_MAG_G'
    lookup[''] = 'GAIA_MAG_G_ERR'
    lookup[''] = 'GAIA_MAG_BP'
    lookup[''] = 'GAIA_MAG_BP_ERR'
    lookup[''] = 'GAIA_MAG_RP'
    lookup[''] = 'GAIA_MAG_RP_ERR'

    lookup.pop('')

    return lookup


def get_formats():

    formats = {}

    formats['target:targprio'] = int
    formats['target:targra'] = float
    formats['target:targdec'] = float
    formats['target:targepoch'] = float
    formats['target:targpmra'] = float
    formats['target:targpmdec'] = float
    formats['target:targparal'] = float
    formats['observation:pa'] = float
    formats['dithering:apply_dither'] = int
    formats['target:mag_g'] = float
    formats['target:emag_g'] = float
    formats['target:mag_r'] = float
    formats['target:emag_r'] = float
    formats['target:mag_i'] = float
    formats['target:emag_i'] = float
    formats['target:mag_gg'] = float
    formats['target:emag_gg'] = float
    formats['target:mag_bp'] = float
    formats['target:emag_bp'] = float
    formats['target:mag_rp'] = float
    formats['target:emag_rp'] = float

    return formats


def get_xml_data(xml_filename):

    xml_data = {}

    dom = xml.dom.minidom.parse(xml_filename)

    root = dom.childNodes[0]

    programme = root.childNodes[3]
    observation = root.childNodes[5]

    obsconstraints = dom.getElementsByTagName('obsconstraints')[0]
    dithering = dom.getElementsByTagName('dithering')[0]
    target = dom.getElementsByTagName('target')

    # configure = dom.getElementsByTagName('configure')[0]
    # field = dom.getElementsByTagName('field')[0]

    # offset = observation.getElementsByTagName('offsets')[0]

    # base_target = field.getElementsByTagName('target')[0]
    # targets_base = field.getElementsByTagName('target')

    xml_data['dom'] = dom
    xml_data['root'] = root
    xml_data['programme'] = programme
    xml_data['observation'] = observation
    xml_data['obsconstraints'] = obsconstraints
    xml_data['dithering'] = dithering
    xml_data['target'] = target

    return xml_data


def get_data(xml_filename_list):

    # Get dictionaries with the lookup information and the formats

    lookup = get_lookup()
    formats = get_formats()

    # Create an dictionary with the desired keywords

    data = {lookup[key]: [] for key in lookup.keys()}

    # For each XML file
    
    for xml_filename in xml_filename_list:

        # Get the data from the XML file

        xml_data = get_xml_data(xml_filename)

        # For each target in the XML

        for target in xml_data['target']:

            # Skip the target if it is a guide or sky fibre

            if str(target.getAttribute('targuse')) != 'T':
                continue

            # For each key in the lookup
            
            for key in lookup.keys():
                element_name, attribute_name = key.split(':')

                if element_name == 'target':
                    element = target
                else:
                    element = xml_data[element_name]

                raw_value = element.getAttribute(attribute_name)

                if key in formats.keys():
                    value = formats[key](raw_value)
                else:
                    value = str(raw_value)

                data[lookup[key]].append(value)

    return data


def write_data_as_in_template(data, fits_template, outname):
    columns = []

    template = fits.open(fits_template)

    for key in data.keys():
        c_base = template[1].columns[key]

        name = c_base.name
        if c_base.name == 'IFU_PA':
            name = 'IFU_PA_REQUEST'

        if 'A' in c_base.format:
            fmt = c_base.disp
            if fmt == 'A':
                col_format = fits.column._ColumnFormat('A')
            else:
                fno = fmt.split('A')[-1]
                col_format = fits.column._ColumnFormat('%sA'%(fno))
        else:
            col_format = c_base.format

        column = fits.Column(name=name, format=col_format,
                             disp=c_base.disp, unit=c_base.unit,
                             null=c_base.null, array=data[key])

        columns.append(column)

    new_coldef = fits.ColDefs(columns)
    new_table = fits.BinTableHDU.from_columns(new_coldef)
    new_table.update()
    template[1] = new_table
    template.writeto(outname,overwrite=True)


if __name__ == '__main__':

    import os
    import glob

    xml_filename_list = glob.glob('../stage5/input/*.xml')

    fits_template = '../../test_data/stage0_base.fits'
    output_dir = './output/'
    outname = output_dir + 'WC_IFU.fits'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    data = get_data(xml_filename_list)

    write_data_as_in_template(data, fits_template, outname)

