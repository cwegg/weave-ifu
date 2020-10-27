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
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
from astropy import coordinates

from workflow.utils.get_data_from_xmls import \
    get_obs_mode, get_coord_of_first_field, get_spa_data_per_field_of_one_xml


class _IFUPlot():

    def __init__(self, xml_file, output_dir='', verbose=False):

        # Save the input XML filename

        self.xml_file = xml_file

        # Save the verbose mode

        self.verbose = verbose

        # Set the filename for the PDF file and create its object

        xml_basename_wo_ext = os.path.splitext(os.path.basename(xml_file))[0]

        self.pdf_filename = os.path.join(output_dir,
                                         xml_basename_wo_ext + '.pdf')

        self.pdf = PdfPages(self.pdf_filename)

        # Get the obs_mode and assert that it is an IFU XML

        self.obs_mode = get_obs_mode(xml_file)

        assert self.obs_mode in ['LIFU', 'mIFU']
    
        # Get the data from the XML file

        self.data_dict_list = get_spa_data_per_field_of_one_xml(self.xml_file)

        # Set the aspect to be used in the figures

        self.figaspect = 1 / ((1 + 5 ** 0.5) / 2)

        # Set the aspect ratio

        self.field_ra, self.field_dec = get_coord_of_first_field(xml_file)

        self.axaspect = np.cos(np.deg2rad(self.field_dec))

        # Set the fibre size, the list with the central spaxels and the field
        # of views

        self._set_fibre_size()
        self._set_central_spaxel_list()
        self._set_fov_list()

        # Make the plots

        self._make_plots()

    def _set_central_spaxel_list(self):
        raise NotImplementedError

    def _set_fov_list(self):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.pdf.__exit__(exc_type, exc_val, exc_tb)

    def close(self):
        self.pdf.close()

    def get_pdf_filename(self):
        return self.pdf_filename

    def _get_fig_ax_vect(self):
        raise NotImplementedError

    def _set_aspect_on_axes(self, ax_vect):
        for ax in ax_vect:
            ax.set_aspect(self.axaspect)

    def _set_ax_lims(self, ax_vect):

        data_dict = self.data_dict_list[0]

        for i, ax in enumerate(ax_vect):
            central_spaxel = self.central_spaxel_list[i]
            fov = self.fov_list[i]

            try:
                index = data_dict['IFU_SPAXEL'].index(central_spaxel)
                ra = data_dict['GAIA_RA'][index]
                dec = data_dict['GAIA_DEC'][index]
            except ValueError:
                # If the central bundle is not present, we put the limits out
                # of the range of the coordinates
                ra = 1000
                dec = self.field_dec

            xlim = [ra + (fov / 2) / self.axaspect,
                    ra - (fov / 2) / self.axaspect]
            ylim = [dec - (fov / 2), dec + (fov / 2)]

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

    def _clean_ticks_on_axes(self, ax_vect):
        for ax in ax_vect:
            ax.set_xticks([])
            ax.set_yticks([])

    def _get_setted_fig_ax_vect(self):

        fig, ax_vect = self._get_fig_ax_vect()

        self._set_aspect_on_axes(ax_vect)
        self._set_ax_lims(ax_vect)
        self._clean_ticks_on_axes(ax_vect)

        return fig, ax_vect

    def _set_title_on_fig(self, fig, suffix_title='', x=0.425, y=0.91):
        title = os.path.basename(self.xml_file) + suffix_title
        fig.text(x, y, title, ha='center', va='bottom')

    def _set_title_on_axes(self, ax_vect):
        raise NotImplementedError

    def _set_fig_legend(self, fig, text_list, title=True,
                        x=0.82, y=0.87, ha='left', va='top'):

        for i, element in enumerate(text_list):

            if (i == 0) or (title is False):
                fontsize = 10
            else:
                fontsize = 8

            if type(element) == str:
                text_line = element
                color = 'k'
            else:
                text_line, color = element
            text = i * '\n' + text_line

            fig.text(x, y, text, color=color, ha=ha, va=va, fontsize=fontsize)

    def _set_fig_colorbar(self, fig, mappable, left=0.85, height=0.8,
                          fraction=1):

        aux_ax = fig.add_axes([left, 0.07, 0.05, height])

        fig.colorbar(mappable, ax=aux_ax, fraction=fraction)

        aux_ax.remove()

    def _plot_empty_fig(self):

        fig, ax_vect = self._get_setted_fig_ax_vect()

        self.pdf.savefig(fig)
        plt.close(fig)

    def _assign_color_pointings(self, i, legend=False):

        color_list = ['r', 'g', 'b', 'orange', 'brown', 'pink', 'purple']

        if legend is False:
            facecolor = color_list[i]
            edgecolor = None
            result = (facecolor, edgecolor)
        else:
            legend_text_list = [
                ('Field #{}'.format(j + 1),
                 color_list[j % (len(color_list) - 1)])
                for j in range(i)]

            result = legend_text_list

        return result

    def _assign_color_attrib(self, type_attrib, attrib=None, set_of_values=None,
                             legend=False):

        if legend is False:
            assert attrib is not None

        color_dict = OrderedDict()

        if type_attrib == 'TARGUSE':
            color_dict['T'] = ('b', None)
            color_dict['S'] = ('g', None)
            color_dict['C'] = ('r', None)
            color_dict['R'] = ('orange', None)
        elif type_attrib == 'TARGCLASS':
            color_dict['']          = ('w', 'k')
            color_dict['UNKNOWN']   = ('brown', None)
            color_dict['SKY']       = ('g', None)
            color_dict['GALAXY']    = ('r', None)
            color_dict['NEBULA']    = ('orange', None)
            color_dict['QSO']       = ('pink', None)
            color_dict['STAR']      = ('b', None)
            color_dict['STAR_BHB']  = ('LightBlue', None)
            color_dict['STAR_CEP']  = ('LightSteelBlue', None)
            color_dict['STAR_EM']   = ('LightSkyBlue', None)
            color_dict['STAR_EMP']  = ('Turquoise', None)
            color_dict['STAR_FGK']  = ('Cyan', None)
            color_dict['STAR_IB']   = ('DeepSkyBlue', None)
            color_dict['STAR_MLT']  = ('DodgerBlue', None)
            color_dict['STAR_MLUM'] = ('SteelBlue', None)
            color_dict['STAR_OB']   = ('CornflowerBlue', None)
            color_dict['STAR_BA']   = ('SlateBlue', None)
            color_dict['STAR_RRL']  = ('MediumSlateBlue', None)
            color_dict['STAR_VAR']  = ('BlueViolet', None)
            color_dict['STAR_WD']   = ('MidnightBlue', None)
            color_dict['STAR_YSO']  = ('DarkBlue', None)
        else:
            color_list = [
                'r', 'g', 'b', 'orange', 'deeppink',
                'darkred', 'olivedrab', 'mediumslateblue', 'orangered', 'darkviolet',
                'magenta', 'lime', 'teal', 'gold', 'indigo',
                'salmon', 'purple', 'deepskyblue', 'yellow', 'pink'
            ]

            if '' in set_of_values:
                color_dict[''] = ('w', 'k')

            copy_set_of_values = set_of_values[:]
            if '' in copy_set_of_values:
                copy_set_of_values.remove('')

            for i, value in enumerate(copy_set_of_values):
                color_dict[value] = (color_list[i % len(color_list)], None)

        if legend is False:
            facecolor, edgecolor = color_dict[attrib]
            result = (facecolor, edgecolor)
        else:
            legend_text_list = [type_attrib, '']

            if '' in color_dict.keys():
                color_dict['-Empty-'] = color_dict[''][1]
                color_dict.pop('')

            for key in color_dict.keys():
                legend_text_list.append((key, color_dict[key][0]))

            result = legend_text_list

        return result

    def _plot_on_ax_vect(self, ax_vect, contents='pointings', field=None):

        width = self.fibre_size
        height = self.fibre_size / self.axaspect

        if contents == 'pointings':
            set_of_values = None
            legend_text_list = self._assign_color_pointings(
                len(self.data_dict_list), legend=True)
            mappable = None
            alpha = 0.5
        elif contents in ['TARGUSE', 'TARGCLASS']:
            set_of_values = None
            legend_text_list = self._assign_color_attrib(contents, legend=True)
            mappable = None
            alpha = 1
        elif contents in ['TARGCAT', 'TARGID', 'TARGNAME', 'TARGPROG',
                          'TARGSRVY']:
            values_array = np.concatenate(
                [data_dict[contents] for data_dict in self.data_dict_list])
            set_of_values = list(set(values_array))
            set_of_values.sort()

            legend_text_list = self._assign_color_attrib(
                contents, set_of_values=set_of_values, legend=True)
            mappable = None
            alpha = 1
        else:
            set_of_values = None
            legend_text_list = None

            values_array = np.concatenate(
                [data_dict[contents] for data_dict in self.data_dict_list])

            if np.all(np.isnan(values_array)):
                vmin = -1
                vmax =  1
            else:
                vmin = np.nanmin(values_array)
                vmax = np.nanmax(values_array)

                if vmin == vmax:
                    vmin -= 0.5
                    vmax += 0.5

            mappable = cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax),
                                         cmap=cm.rainbow)
            alpha = 1

        for i, ax in enumerate(ax_vect):

            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            for j, data_dict in enumerate(self.data_dict_list):

                # Skip this field if a specific one has been requested

                if field is not None:
                    if j != field:
                        continue

                # Set the color of the contents are pointings
                
                if contents == 'pointings':
                    facecolor, edgecolor = self._assign_color_pointings(j)

                # Plot each fibre in the field of view of the axis

                for k, (ra, dec) in enumerate(zip(data_dict['GAIA_RA'],
                                                  data_dict['GAIA_DEC'])):

                    if ((ra >= xlim[1]) and (ra <= xlim[0]) and
                        (dec >= ylim[0]) and (dec <= ylim[1])):

                        if contents in ['TARGUSE', 'TARGCLASS', 'TARGCAT',
                                        'TARGID', 'TARGNAME', 'TARGPROG',
                                        'TARGSRVY']:
                            attrib = data_dict[contents][k]
                            facecolor, edgecolor = self._assign_color_attrib(
                                contents, attrib=attrib,
                                set_of_values=set_of_values)
                        elif contents != 'pointings':
                            attrib = data_dict[contents][k]

                            if np.isnan(attrib):
                                facecolor = 'w'
                                edgecolor = 'k'
                            else:
                                facecolor = mappable.to_rgba(attrib)
                                edgecolor = None

                        ax.add_artist(
                            Ellipse((ra, dec), width=width, height=height,
                                    facecolor=facecolor, edgecolor=edgecolor,
                                    alpha=alpha))

        return legend_text_list, mappable

    def _plot_pointings(self):

        fig, ax_vect = self._get_setted_fig_ax_vect()

        self._set_title_on_fig(fig)
        self._set_title_on_axes(ax_vect)

        legend_text_list, mappable = self._plot_on_ax_vect(ax_vect,
                                                           contents='pointings')

        self._set_fig_legend(fig, legend_text_list, title=False)

        self.pdf.savefig(fig)
        plt.close(fig)

    def _plot_attrib(self, type_attrib, split=True):

        if split is True:

            for i in range(len(self.data_dict_list)):

                fig, ax_vect = self._get_setted_fig_ax_vect()

                suffix_title = ' - {} - Field #{}'.format(type_attrib, i + 1)

                self._set_title_on_fig(fig, suffix_title=suffix_title)
                self._set_title_on_axes(ax_vect)

                legend_text_list, mappable = self._plot_on_ax_vect(
                    ax_vect, field=i, contents=type_attrib)

                if legend_text_list is not None:
                    self._set_fig_legend(fig, legend_text_list)

                if mappable is not None:
                    self._set_fig_colorbar(fig, mappable)

                self.pdf.savefig(fig)
                plt.close(fig)

        else:

            fig, ax_vect = self._get_setted_fig_ax_vect()

            suffix_title = ' - {}'.format(type_attrib)

            self._set_title_on_fig(fig, suffix_title=suffix_title)
            self._set_title_on_axes(ax_vect)

            legend_text_list, mappable = self._plot_on_ax_vect(
                ax_vect, contents=type_attrib)

            if legend_text_list is not None:
                self._set_fig_legend(fig, legend_text_list)

            if mappable is not None:
                self._set_fig_colorbar(fig, mappable)

            self.pdf.savefig(fig)
            plt.close(fig)

    def _make_plots(self):

        logging.debug('plotting pointings of {}'.format(self.xml_file))
        self._plot_pointings()

        mag_list = ['MAG_G', 'MAG_R', 'MAG_I',
                    'GAIA_MAG_G', 'GAIA_MAG_BP', 'GAIA_MAG_RP']

        err_mag_list = [mag + '_ERR' for mag in mag_list]

        non_verbose_attrib_list = mag_list + ['TARGUSE', 'TARGCLASS']

        verbose_attrib_list = (
            ['TARGCAT', 'TARGID', 'TARGNAME', 'TARGPROG', 'TARGSRVY'] +
             err_mag_list +
            ['GAIA_EPOCH', 'GAIA_PARAL', 'GAIA_PARAL', 'GAIA_PMRA',
             'GAIA_PMDEC', 'IFU_PA', 'TARGPRIO'])

        for attrib in non_verbose_attrib_list:
            logging.debug('plotting {} of {}'.format(attrib, self.xml_file))
            self._plot_attrib(attrib)

        if self.verbose is True:

            for attrib in verbose_attrib_list:
                logging.debug('plotting {} of {}'.format(attrib, self.xml_file))
                self._plot_attrib(attrib)

        self._plot_empty_fig()

    
class _LIFUPlot(_IFUPlot):

    def _set_fibre_size(self):
        self.fibre_size = 2.6 / 3600

    def _set_central_spaxel_list(self):
        self.central_spaxel_list = [
            'C14', 'S25', 'S32', 'S39', 'S46', 'S18', 'S11', 'S04', 'S53']

    def _set_fov_list(self):
        self.fov_list = [110/3600] + [25/3600 for i in range(8)]

    def _get_fig_ax_vect(self, nrows=4, ncols=6,
                         left=0.05, right=0.80, bottom=0.07, top=0.87):

        # Get the size of the fibure and create it

        figsize = plt.figaspect(self.figaspect)

        fig = plt.figure(figsize=figsize)

        # Get a grid for the axes

        gridspec = fig.add_gridspec(
            nrows, ncols, left=left, bottom=bottom, right=right, top=top)

        # Create the axes of the central bundle
        
        central_ax = fig.add_subplot(gridspec[:, 1:(ncols - 1)])

        # Create the axes of the sky bundles

        sky_ax_vect = []

        for i in range(4):
            sky_ax_vect.append(fig.add_subplot(gridspec[i, 0]))

        for i in range(4):
            sky_ax_vect.append(fig.add_subplot(gridspec[i, ncols - 1]))

        # Set some properties of the axes

        ax_vect = [central_ax] + sky_ax_vect

        return fig, ax_vect

    def _set_title_on_axes(self, ax_vect):

        ax_title_list = ['Ul', 'uL', 'dL', 'Dl', 'Ur', 'uR', 'dR', 'Dr']

        for i, ax_title in enumerate(ax_title_list):
            ax_vect[i + 1].set_title(ax_title, pad=0, fontdict={'fontsize': 10})

        data_dict = self.data_dict_list[0]

        for i, central_spaxel in enumerate(self.central_spaxel_list):

            index = data_dict['IFU_SPAXEL'].index(central_spaxel)
            ra = data_dict['GAIA_RA'][index]
            dec = data_dict['GAIA_DEC'][index]

            if i == 0:
                centre_coord = coordinates.SkyCoord(ra, dec, unit='deg')
            else:
                sky_bundle_coord = coordinates.SkyCoord(ra, dec, unit='deg')

                dist = centre_coord.separation(sky_bundle_coord)
                pos_angle = centre_coord.position_angle(sky_bundle_coord)

                text_coord = centre_coord.directional_offset_by(pos_angle,
                                                                dist / 5)

                ax_vect[0].text(text_coord.ra.deg, text_coord.dec.deg,
                                ax_title_list[i - 1])


class _MIFUPlot(_IFUPlot):

    def _set_fibre_size(self):
        self.fibre_size = 1.3 / 3600

    def _set_central_spaxel_list(self):
        self.central_spaxel_list = ['m{:02d}C04'.format(i + 1)
                                    for i in range(20)]

    def _set_fov_list(self):
        self.fov_list = [15/3600 for i in range(20)]

    def _get_fig_ax_vect(self, nrows=4, ncols=5,
                         left=0.10, right=0.75, bottom=0.07, top=0.87):

        # Get the size of the fibure and create it

        figsize = plt.figaspect(self.figaspect)

        fig = plt.figure(figsize=figsize)
    
        # Create a figure and the axes

        gridspec_kw = {'left': left, 'right': right,
                       'bottom': bottom, 'top': top}

        fig, ax_array = plt.subplots(
            figsize=figsize, nrows=nrows, ncols=ncols, gridspec_kw=gridspec_kw)

        ax_vect = ax_array.flatten()

        return fig, ax_vect

    def _set_title_on_axes(self, ax_vect):
        for i, ax in enumerate(ax_vect):
            ax.set_title('m{:02d}'.format(i + 1), pad=0,
                         fontdict={'fontsize': 10})


def plot_data_from_xml(xml_file, output_dir='output', verbose=False):
    """
    Plot targets contained in an IFU driver catalogue with Aladin.

    Parameters xml_file, output_dir=args.dir)
    ----------
    xml_file : str
        The name of a XML file.
    output_dir : str, optional
        The directory which will contain the plot.
    verbose : bool, optional
        Option for verbose mode.

    Returns
    -------
    output_filename : str
        The filename of the output plot.
    """

    # Get rid of a matplotlib warning
    
    plt.rcParams.update({'figure.max_open_warning': 0})
    
    # Get the obs_mode of the XML
    
    obs_mode = get_obs_mode(xml_file)

    # Use a plotting class depending on the obs_mode

    if obs_mode == 'LIFU':

        with _LIFUPlot(xml_file, output_dir=output_dir, verbose=verbose) as plot:
            output_filename = plot.get_pdf_filename()
        
    elif obs_mode == 'mIFU':

        with _MIFUPlot(xml_file, output_dir=output_dir, verbose=verbose) as plot:
            output_filename = plot.get_pdf_filename()

    else:

        logging.error(
            'unexpected obs_mode in {}: {}'.format(xml_file, obs_mode))
    
        output_filename = None
        
    return output_filename


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Plot the data of one or more XML files.')

    parser.add_argument('xml_file', nargs='+', help='name of a XML file')

    parser.add_argument('--dir', default='output',
                        help='the directory which will contain the plots')

    parser.add_argument('--verbose', action='store_true',
                        help='option for verbose mode')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()
    
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.dir):
        logging.info('Creating the output directory')
        os.mkdir(args.dir)
    
    for xml_file in args.xml_file:
        plot_data_from_xml(xml_file, output_dir=args.dir, verbose=args.verbose)

