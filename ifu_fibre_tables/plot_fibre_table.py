#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# <one line to give the program's name and a brief idea of what it does.>
# Copyright (C) 2020  Luis Peralta de Arriba
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
from astropy.io import ascii


def _read_table(filename):

    # Get the data using astropy

    data = ascii.read(filename, header_start=1)

    # Read the file to get manually some values from the header afterwards

    with open(filename) as f:
        lines = f.readlines()

    # Get the datamver
    # The first line should like that:
    # '#DATAMVER 7.00'

    first_line = lines[0]

    datamver = first_line.split()[1]

    # Get the units of x and y
    # They should be the same and appear in the third line:
    # '# mm     mm       --         --        --'

    third_line = lines[2]

    units_x = third_line.split()[1]
    units_y = third_line.split()[2]

    assert units_x == units_y
    
    units = units_x

    return datamver, units, data


def _get_value_from_ifu_spaxel(data, col_name, ifu_spaxel):

    mask = (data['ifu_spaxel'] == ifu_spaxel)
    value = data[col_name][mask][0]

    return value


def _mv_lifu_sky_fibres(input_data, shift_factor=4):

    data = input_data.copy()

    # Get the IFU spaxels of the fibres in the centres of the sky bundles

    central_ifu_spaxel_sky_bundles = ['S{:02d}'.format(4 + i * 7)
                                      for i in range(8)]

    # Get the coordinates of the centres of sky bundles

    sky_bundles_x_list = np.array(
        [_get_value_from_ifu_spaxel(data, 'x', ifu_spaxel)
         for ifu_spaxel in central_ifu_spaxel_sky_bundles])

    sky_bundles_y_list = np.array(
        [_get_value_from_ifu_spaxel(data, 'y', ifu_spaxel)
         for ifu_spaxel in central_ifu_spaxel_sky_bundles])

    # Set the shifts

    shift_x_list = [(1 - 1.0 / shift_factor) * x for x in sky_bundles_x_list]
    shift_y_list = [(1 - 1.0 / shift_factor) * y for y in sky_bundles_y_list]

    # Get the coordinates of the centres of sky bundles

    for row in data:
        if row['ifu_spaxel'][0] == 'S':

           # Guess the sky bundle depending on its ifu_spaxel

           sky_num = int(row['ifu_spaxel'][1:])
           sky_bundle_i = (sky_num - 1) // 7

           # Move the fibre

           row['x'] -= shift_x_list[sky_bundle_i]
           row['y'] -= shift_y_list[sky_bundle_i]

    return data


def _plot_fibres(ax, x, y, fibre_size, id_array=None, color=(0.9, 0.9, 1.0),
                 alpha=0.25, edge=True, for_printing=False):

    radius = fibre_size / 2.0
    if for_printing is True:
        radius *= 1.7
    edge_outer_radius = 1.05 * radius
    edge_width = 0.05 * radius
    
    if for_printing is True:
        text_size = 3
    else:
        text_size = 2.5
    
    for i in range(len(x)):

        # Add a circle with the fibre size

        circle = Circle((x[i], y[i]), radius=radius, color=color)

        ax.add_patch(circle)

        # Add an edge above the circle

        if edge == True:
            wedge = Wedge((x[i], y[i]), edge_outer_radius, 0.0, 360.0,
                          width=edge_width, alpha=alpha)

            ax.add_patch(wedge)

        # Add a text with the ID

        if id_array is not None:
            ax.text(x[i], y[i], id_array[i], size=text_size,
                    ha='center', va='center')


def plot_lifu_table(lifu_filename, output_filename, id_name,
                    sky_shift_factor=4, fibre_size_mm=0.146, add_mini_ax=True,
                    for_printing=False):

    datamver, units, data = _read_table(lifu_filename)

    assert units == 'mm'

    fig, ax = plt.subplots(figsize=(11.69, 8.27))

    fig.suptitle('{} (version {}): {}'.format(
        lifu_filename, datamver, id_name.upper()), y=0.93)

    mv_data = _mv_lifu_sky_fibres(data, shift_factor=sky_shift_factor)

    _plot_fibres(ax, mv_data['x'], mv_data['y'], fibre_size_mm,
                 id_array=mv_data[id_name], for_printing=for_printing)

    ax.set_xlim([np.min(mv_data['x']) - 2.5 * fibre_size_mm,
                 np.max(mv_data['x']) + 2.5 * fibre_size_mm])
    ax.set_ylim([np.min(mv_data['y']) - 1.5 * fibre_size_mm,
                 np.max(mv_data['y']) + 1.5 * fibre_size_mm])

    ax.set_aspect('equal')

    ax.tick_params(direction='in', bottom=True, top=True, left=True, right=True)

    ax.set_xlabel('x ({})'.format(units))
    ax.set_ylabel('y ({})'.format(units))

    if add_mini_ax is True:

        mini_ax = fig.add_axes((0.7, 0.12, 0.1, 0.1))

        _plot_fibres(mini_ax, data['x'], data['y'], fibre_size_mm,
                     color=(0.5, 0.5, 0.75), alpha=1.0, edge=False,
                     for_printing=for_printing)

        sky_ifu_spaxel = 'S18'

        sky_x = _get_value_from_ifu_spaxel(data, 'x', sky_ifu_spaxel)
        sky_y = _get_value_from_ifu_spaxel(data, 'y', sky_ifu_spaxel)

        mini_ax.set_xlim([np.min(data['x']) - 20 * fibre_size_mm,
                          np.max(data['x']) + 20 * fibre_size_mm])
        mini_ax.set_ylim([np.min(data['y']) - 20 * fibre_size_mm,
                          np.max(data['y']) + 20 * fibre_size_mm])

        mini_ax.set_aspect('equal')

        mini_ax.set_xticks([])
        mini_ax.set_yticks([])

    fig.savefig(output_filename)


def plot_mifu_table(mifu_filename, output_filename, id_name, num_mifu=20,
                    fibre_size_mm=0.073, for_printing=False):

    datamver, units, data = _read_table(mifu_filename)

    assert units == 'mm'

    mifu_id_array = np.array([int(ifu_spaxel[1:3])
                              for ifu_spaxel in data['ifu_spaxel']])

    assert (min(set(mifu_id_array)) == 1)
    assert (max(set(mifu_id_array)) == num_mifu)

    nrows = 4
    ncols = 5

    fig, ax_array = plt.subplots(nrows=nrows, ncols=ncols,
                                 figsize=(11.69, 8.27))

    ax_vect = ax_array.flatten()

    fig.suptitle('{} (version {}): {}'.format(
        mifu_filename, datamver, id_name.upper()), y=0.96)

    for i in range(num_mifu):

        mifu_i = i + 1

        ax_vect[i].set_title('mIFU #{}'.format(mifu_i))

        mask = (mifu_id_array == mifu_i)

        mifu_i_data = data[mask]
        x = mifu_i_data['x']
        y = mifu_i_data['y']
        id_array = mifu_i_data[id_name]

        if id_name == 'ifu_spaxel':
            id_array = np.array(['{}\n{}'.format(value[:3], value[3:])
                                 for value in id_array])

        _plot_fibres(ax_vect[i], x, y, fibre_size_mm, id_array=id_array,
                     for_printing=for_printing)
        
        assert np.min(mifu_i_data['x']) > -0.5
        assert np.max(mifu_i_data['x']) < 0.5
        assert np.min(mifu_i_data['y']) > -0.5
        assert np.max(mifu_i_data['y']) < 0.5

        ax_vect[i].set_xticks([-0.5, 0, 0.5])
        ax_vect[i].set_yticks([-0.5, 0, 0.5])

        ax_vect[i].set_xlim([-0.5, 0.5])
        ax_vect[i].set_ylim([-0.5, 0.5])

        ax_vect[i].set_aspect('equal')

        ax_vect[i].tick_params(direction='in',
                               bottom=True, top=True, left=True, right=True)

        if ((i // ncols) != (nrows - 1)):
            ax_vect[i].set_xticklabels(['', '', ''])
        if ((i % ncols) != 0):
            ax_vect[i].set_yticklabels(['', '', ''])
        
        if (((i // ncols) == (nrows - 1)) and
            ((i % ncols) == (ncols - ncols // 2 - 1))):
            ax_vect[i].set_xlabel('x ({})'.format(units))
        # if ((i % ncols) == 0):
        #     ax_vect[i].set_ylabel('y ({})'.format(units))

    fig.text(0.085, 0.495, 'y ({})'.format(units),
             ha='center', va='center', rotation=90)

    fig.savefig(output_filename)

    
def plot_fibre_table(lifu_filename, mifu_filename, for_printing=False):

    # Set the patterng for the filenames    
    
    if for_printing is True:
        filename_pattern = '{}-{}-for_printing.pdf'
    else:
        filename_pattern = '{}-{}.pdf'

    # For each kind of ID

    for id_name in ['ifu_spaxel', 'ifu_slits', 'fibreid']:

        # Plot it for the LIFU
        
        lifu_output_filename = filename_pattern.format(
            os.path.splitext(lifu_filename)[0], id_name)

        plot_lifu_table(lifu_filename, lifu_output_filename, id_name,
                        for_printing=for_printing)
        
        # Plot it for the mIFU

        mifu_output_filename = filename_pattern.format(
            os.path.splitext(mifu_filename)[0], id_name)

        plot_mifu_table(mifu_filename, mifu_output_filename, id_name,
                        for_printing=for_printing)

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description='Create plots of the WEAVE IFU tables')

    parser.add_argument('--in_lifu', dest='lifu_filename',
                        default='LIFUfibreTable.dat',
                        help='LIFU fibre table')

    parser.add_argument('--in_mifu', dest='mifu_filename',
                        default='mIFUfibreTable.dat',
                        help='mIFU fibre table')

    parser.add_argument('--for_printing', dest='for_printing',
                        action='store_true',
                        help='scale the fibre sizes for easy printing')

    args = parser.parse_args()

    plot_fibre_table(args.lifu_filename, args.mifu_filename,
                     for_printing=args.for_printing)

