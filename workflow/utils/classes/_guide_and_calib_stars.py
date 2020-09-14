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


"""
Retrieve of WEAVE guide and calib stars.

The authors of this module are:

- David Murphy (dmurphy@ast.cam.ac.uk),
  Cambridge Astronomical Survey Unit (CASU, IoA).
- Luis Peralta de Arriba (lperalta@ast.cam.ac.uk),
  Cambridge Astronomical Survey Unit (CASU, IoA).

Notes
----- 
The dependencies of this module are:

- Python core
- numpy
- matplotlib
- astropy
- astropy_healpix
"""


import logging
import tempfile as _tempfile

import numpy as _np
from astropy import coordinates as _coordinates
from astropy.io import fits as _fits
from astropy.table import Column as _Column
from astropy.table import Table as _Table
from astropy.table import vstack as _vstack
from astropy_healpix import HEALPix as _HEALPix
import matplotlib as _mpl
_mpl.use('Agg')
import matplotlib.pyplot as _plt
from matplotlib.patches import Ellipse as _Ellipse

from workflow.utils.get_resources import get_guide_cat as _get_guide_cat
from workflow.utils.get_resources import get_calib_cat as _get_calib_cat


class _AuxStars:


    def __init__(self, ra, dec, obsmode, dither_size_arcsec=10.0, nside=32):

        self.star_type = ''

        self.ra = ra
        self.dec = dec
        self.obsmode = obsmode
        self.dither_size_arcsec = dither_size_arcsec
        self.nside = nside
        
        self._assert_args()

        self._set_geometry()

        self.healpix_indices = None
        self.full_table = None
        self.useful_table = None
        self.selected_pa = None
        self.selected_table = None

        
    def _assert_args(self):
        
        if not ((self.ra >= 0.0) and (self.ra < 360.0)):
            logging.error('invalid ra value: {}'.format(self.ra))
            raise ValueError
        
        if not ((self.dec >= -90.0) and (self.dec < 90.0)):
            logging.error('invalid dec value: {}'.format(self.dec))
            raise ValueError

        if self.obsmode not in ['MOS', 'LIFU', 'mIFU']:
            logging.error('invalid obsmode: {}'.format(obsmode))
            raise ValueError
        
        
    def _set_geometry(self):

        # NB: The dither size is substracted twice to avoid edge issues

        if (self.obsmode == 'LIFU'):
            
            logging.critical(
                'Do not trust on the output of this code by now: ' +
                'This software is great, but the assumed position of the LIFU '
                'guide camera must to be confirmed')

            # Set the position of the guide camera

            self.cam_x_offset =   0.0 / 60.0
            self.cam_y_offset = -27.7 / 60.0

            cam_offset_dist = _np.hypot(self.cam_x_offset, self.cam_y_offset)

            # Set the effective radious of the field of view of the guide camera
            
            cam_width  =  4.0 / 60.0
            cam_height = 3.75 / 60.0
            
            self.fov_radious = (min(cam_width, cam_height) / 2 -
                                2 * self.dither_size_arcsec / 3600)

            # Set the limits for looking for stars
            
            self.ra_min = (self.ra - ((cam_offset_dist + self.fov_radious) /
                                      _np.cos(_np.deg2rad(self.dec))))
            self.ra_max = (self.ra + ((cam_offset_dist + self.fov_radious) /
                                      _np.cos(_np.deg2rad(self.dec))))
            self.dec_min = self.dec - (cam_offset_dist + self.fov_radious)
            self.dec_max = self.dec + (cam_offset_dist + self.fov_radious)

        else:

            # Set the geometry of the guiding camera to None

            self.cam_x_offset = None
            self.cam_y_offset = None

            # Set the radious of the field of view

            self.fov_radious = 1.0 - 2 * self.dither_size_arcsec / 3600
            
            # Set the limits for looking for stars

            self.ra_min = self.ra - (self.fov_radious / 
                                     _np.cos(_np.deg2rad(self.dec)))
            self.ra_max = self.ra + (self.fov_radious / 
                                     _np.cos(_np.deg2rad(self.dec)))
            self.dec_min = self.dec - self.fov_radious
            self.dec_max = self.dec + self.fov_radious

        # Fix edge issues for declination (for right ascension is not needed)
        
        if self.dec_min < -90.0:
            self.dec_min = -90.0

        if self.dec_max > 90.0:
            self.dec_max = 90.0


    def _set_healpix_indices(self, num_sample=20):

        if self.healpix_indices is None:

            hp = _HEALPix(nside=self.nside, order='nested',
                          frame=_coordinates.ICRS())

            healpix_indices = []

            dec_min = self.dec_min
            dec_max = self.dec_max

            if dec_min < -90.0:
                dec_min = -90.0
            if dec_max > 90.0:
                dec_max = 90.0

            for ra in _np.linspace(self.ra_min, self.ra_max, num_sample):
                for dec in _np.linspace(dec_min, dec_max, num_sample):

                    healpix_index = hp.skycoord_to_healpix(
                        _coordinates.SkyCoord(ra, dec, unit='deg'))

                    if healpix_index not in healpix_indices:
                        healpix_indices.append(healpix_index)

            healpix_indices.sort()

        self.healpix_indices = healpix_indices


    def _get_cat(self, healpix_index, directory):

        raise NotImplementedError


    def _retrieve_cats(self):

        assert self.healpix_indices is not None

        table_list = []

        for healpix_index in self.healpix_indices:

            with _tempfile.TemporaryDirectory() as tmp_dir:

                file_path = self._get_cat(healpix_index, tmp_dir)

                with _fits.open(file_path) as hdu_list:
                    aux_table = _Table(hdu_list[1].data)

            table_list.append(aux_table)

        full_table = _vstack(table_list)

        return full_table


    def _add_distance_and_angle_to_table(self, full_table):

        num_rows = len(full_table)

        centre_coord = _coordinates.SkyCoord(self.ra, self.dec, unit='deg')

        distance_list = []
        angle_list = []

        for row in full_table:

            row_coord = _coordinates.SkyCoord(row['GAIA_RA'], row['GAIA_DEC'],
                                              unit='deg')

            distance_list.append(centre_coord.separation(row_coord).deg)
            angle_list.append(centre_coord.position_angle(row_coord).wrap_at(
                _coordinates.Angle(360, unit='deg')).deg)

        distance_column = _Column(distance_list, name='DISTANCE')
        angle_column = _Column(angle_list, name='ANGLE')

        full_table.add_column(distance_column)
        full_table.add_column(angle_column)

        return full_table


    def _set_full_table(self):

        if self.full_table is None:

            self._set_healpix_indices()
            full_table = self._retrieve_cats()

            full_table = self._add_distance_and_angle_to_table(full_table)

            self.full_table = full_table


    def _get_stars_in_plate_a_b(self, table):

        mask = (table['DISTANCE'] <= self.fov_radious)

        fov_table = table[mask]

        return fov_table

    
    def _get_stars_in_guide_cam_ring(self, table):

        cam_offset_dist = _np.hypot(self.cam_x_offset, self.cam_y_offset)

        min_radious = cam_offset_dist - self.fov_radious 
        max_radious = cam_offset_dist + self.fov_radious 

        mask = ((table['DISTANCE'] >= min_radious) *
                (table['DISTANCE'] <= max_radious))

        ring_table = table[mask]

        return ring_table


    def _set_useful_table(self):

        if self.obsmode != 'LIFU':

            useful_table = self._get_stars_in_plate_a_b(self.full_table)

        else:
            
            useful_table = self._get_stars_in_guide_cam_ring(self.full_table)

        self.useful_table = useful_table

        
    def _select_central_stars_in_plate_a_b(self, table, num_central_stars=1):

        argsort = _np.argsort(table['DISTANCE'])

        selected_indices = argsort[0:num_central_stars]

        mask = _np.array([(i in selected_indices)
                          for i in range(len(table))])

        central_star_table = table[mask]
        non_selected_table = table[_np.logical_not(mask)]

        return central_star_table, non_selected_table


    def _get_index_of_nearest_angle(self, ref_angle, angle_array,
                                    selected_indices=None):

        index = None

        wrapped_angle_list = [
            _coordinates.Angle(angle , unit='deg').wrap_at(
                _coordinates.Angle(ref_angle + 180, unit='deg')).deg
            for angle in angle_array]

        for i, angle_i in enumerate(wrapped_angle_list):

            if (selected_indices is None) or (i not in selected_indices):

                angle_diff = _np.abs(angle_i - ref_angle)

                if (index is None) or (angle_diff < min_angle_diff):
                    index = i
                    min_angle_diff = angle_diff

        return index

    
    def _select_stars_in_ring_of_plate_a_b(self, table, num_stars,
                                           min_cut=0.9, max_cut=1.0):

        min_radious = min_cut * self.fov_radious
        max_radious = max_cut * self.fov_radious

        mask = ((table['DISTANCE'] >= min_radious) *
                (table['DISTANCE'] <= max_radious))
        
        ring_table = table[mask]

        if num_stars is None:

            selected_table = ring_table

        elif num_stars > 0:

            selected_indices = []

            for i in range(num_stars):

                angle = i * 360 / num_stars

                index = self._get_index_of_nearest_angle(
                    angle, ring_table['ANGLE'],
                    selected_indices=selected_indices)

                if index is not None:

                    selected_indices.append(index)

            mask = _np.array([(i in selected_indices)
                              for i in range(len(ring_table))])

            selected_table = ring_table[mask]

        else:

            selected_table = ring_table[0:0]

        return selected_table


    def _get_guide_cam_coord(self, pa):

        cam_offset_dist = _coordinates.Angle(
            _np.hypot(self.cam_x_offset, self.cam_y_offset), unit='deg')
        cam_pa_angle = _coordinates.Angle(
            _np.rad2deg(_np.arctan2(self.cam_y_offset, self.cam_x_offset)) +
            pa - 90, unit='deg')

        centre_coord = _coordinates.SkyCoord(self.ra, self.dec, unit='deg')
        cam_coord = centre_coord.directional_offset_by(cam_pa_angle,
                                                       cam_offset_dist)

        return cam_coord


    def _select_stars_in_guide_cam(self, table, pa_request, num_stars_request):

        # Get the coordinate of the guide cam for the requested PA

        cam_coord = self._get_guide_cam_coord(pa_request)

        # Compute the distances to the centre of the guide cam

        cam_distance_list = []

        for row in table:

            row_coord = _coordinates.SkyCoord(row['GAIA_RA'], row['GAIA_DEC'],
                                              unit='deg')

            cam_distance_list.append(cam_coord.separation(row_coord).deg)

        cam_distance_array = _np.array(cam_distance_list)

        # Get a table with the guide stars in the field of view of the guide
        # camera

        mask = (cam_distance_array <= self.fov_radious)

        pre_selected_table = table[mask]

        # Select the nearest star to the centre of the guide camera if it is in
        # the field of view. Otherwise, return an empty table

        if _np.count_nonzero(mask) > 0:

            argsort = _np.argsort(cam_distance_array[mask])

            if num_stars_request is not None:
                num_stars = num_stars_request
            else:
                num_stars = len(pre_selected_table)

            selected_table = pre_selected_table[argsort][0:num_stars]

        else:

            selected_table = pre_selected_table[0:0]

        return selected_table


    def _select_stars_in_near_guide_cam_position(self, table, pa_request,
                                                 num_stars_request):

        if len(table) > 0:
            
            cam_pa_at_request = _coordinates.Angle(
                _np.rad2deg(_np.arctan2(self.cam_y_offset, self.cam_x_offset)) +
                pa_request - 90, unit='deg').wrap_at(
                    _coordinates.Angle(360, unit='deg')).deg

            argmin = self._get_index_of_nearest_angle(cam_pa_at_request,
                                                      table['ANGLE'])
            
            cam_pa_at_zero = _coordinates.Angle(
                _np.rad2deg(_np.arctan2(self.cam_y_offset, self.cam_x_offset)) -
                90, unit='deg').wrap_at(
                    _coordinates.Angle(360, unit='deg')).deg

            selected_pa = _coordinates.Angle(
                table['ANGLE'][argmin] - cam_pa_at_zero, unit='deg').wrap_at(
                    _coordinates.Angle(360, unit='deg')).deg

            selected_table = self._select_stars_in_guide_cam(table, selected_pa,
                                                             num_stars_request)

        else:

            selected_pa = pa_request
            selected_table = table
            
        return selected_pa, selected_table


    def _set_selected_table(self, num_stars_request, pa_request=None,
                            num_central_stars=1, min_cut=0.9, max_cut=1.0):

        useful_table = self.useful_table

        # If there is not request, let's start with a position angle of 0.0

        if pa_request not in [None, _np.nan]:
            pa_request = pa_request
        else:
            pa_request = 0.0

        if self.obsmode != 'LIFU':

            # For MOS/mIFU, get central stars and other ones in the outer ring
            # of the plate

            selected_pa = pa_request

            central_star_table, non_selected_table = (
                self._select_central_stars_in_plate_a_b(
                    useful_table, num_central_stars=num_central_stars))

            if num_stars_request is not None:
                num_ring_stars = num_stars_request - num_central_stars
            else:
                num_ring_stars = None

            ring_stars_table = self._select_stars_in_ring_of_plate_a_b(
                non_selected_table, num_ring_stars,
                min_cut=min_cut, max_cut=max_cut)

            selected_table = _vstack([central_star_table, ring_stars_table])

        else:

            # For LIFU, try to get a star in the central in the guide cam, and
            # if not possible, choose the nearest available star to the desired
            # position

            selected_table = self._select_stars_in_guide_cam(useful_table,
                                                             pa_request,
                                                             num_stars_request)

            if len(selected_table) > 0:
                selected_pa = pa_request
            else:
                selected_pa, selected_table = (
                    self._select_stars_in_near_guide_cam_position(
                        useful_table, pa_request, num_stars_request))

        self.selected_table = selected_table
        self.selected_pa = float(selected_pa)

        
    def _print_summary(self):

        if self.full_table is not None:
            logging.info(
                '{} {} stars have been retrieved for ({}, {})'.format(
                    len(self.full_table), self.star_type, self.ra, self.dec))

        if self.selected_table is not None:
            logging.info(
                '{} {} stars are useful for ({}, {})'.format(
                    len(self.useful_table), self.star_type, self.ra, self.dec))

        if self.selected_table is not None:
            if self.obsmode != 'LIFU':
                logging.info(
                    '{} {} stars has been selected for ({}, {})'.format(
                        len(self.selected_table), self.star_type,
                        self.ra, self.dec))
            else:
                logging.info(
                    '{} {} stars has been selected for ({}, {})'.format(
                        len(self.selected_table), self.star_type,
                        self.ra, self.dec) +
                    ' PA = {} deg'.format(self.selected_pa))


    def _fix_ra_for_plotting(self, ra_array):

        mask_below_ra_min = (ra_array < self.ra_min)

        if _np.count_nonzero(mask_below_ra_min) > 0:
            ra_array[mask_below_ra_min] = ra_array[mask_below_ra_min] + 360

        mask_above_ra_max = (ra_array > self.ra_max)

        if _np.count_nonzero(mask_above_ra_max) > 0:
            ra_array[mask_above_ra_max] = ra_array[mask_above_ra_max] - 360

        return ra_array

                
    def _plot_table(self, table, ax, color):

        if (table is not None) and (len(table) > 0):

            # Get the coordinates from the table

            ra_array = table['GAIA_RA']
            dec_array = table['GAIA_DEC']

            # Fix coordinates where there are edge problems

            ra_array = self._fix_ra_for_plotting(ra_array)

            # Plot the coordinates

            ax.scatter(ra_array, dec_array, c=color)


    def _plot_useful_region_in_plate(self, ax):
    
        center = (self.ra, self.dec)

        width = self.ra_max - self.ra_min
        height = self.dec_max - self.dec_min

        ellipse = _Ellipse(center, width, height,
                           facecolor='lightgrey', edgecolor='none',
                           linewidth=3, zorder=-10)

        ax.add_patch(ellipse)

                
    def _plot_cut_region_in_plate(self, ax, cut):

        if cut is not None:

            center = (self.ra, self.dec)

            width = (self.ra_max - self.ra_min) * cut
            height = (self.dec_max - self.dec_min) * cut

            ellipse = _Ellipse(center, width, height,
                            facecolor='none', edgecolor='k',
                            linewidth=1, linestyle='--', zorder=-8)

            ax.add_patch(ellipse)

                
    def _plot_useful_region_for_guide_cam(self, ax):
    
        center = (self.ra, self.dec)

        out_width = self.ra_max - self.ra_min
        out_height = self.dec_max - self.dec_min

        out_ellipse = _Ellipse(center, out_width, out_height,
                               facecolor='lightgrey', edgecolor='none',
                               linewidth=0, zorder=-10)

        ax.add_patch(out_ellipse)

        in_width = ((self.ra_max - self.ra_min) -
                    2 * 2 * self.fov_radious / _np.cos(_np.deg2rad(self.dec)))
        in_height = ((self.dec_max - self.dec_min) -
                     2 * 2 * self.fov_radious)

        _ellipse = _Ellipse(center, in_width, in_height,
                               facecolor='w', edgecolor='none',
                               linewidth=0, zorder=-9)

        ax.add_patch(_ellipse)


    def _plot_guide_cam_fov(self, ax, pa, linestyle='-', color='k', zorder=-5):

        if pa not in [None, _np.nan]:

            guide_coord = self._get_guide_cam_coord(pa)

            ra = guide_coord.ra.deg
            dec = guide_coord.dec.deg

            if ra < self.ra_min:
                ra += 360

            if ra > self.ra_max:
                ra -= 360

            center = (ra, dec)

            width = 2 * self.fov_radious / _np.cos(_np.deg2rad(self.dec))
            height = 2 * self.fov_radious

            ellipse = _Ellipse(center, width, height,
                               facecolor='none', edgecolor=color,
                               linewidth=1, linestyle=linestyle, zorder=zorder)

            ax.add_patch(ellipse)

            ax.scatter([ra], [dec], c=color, marker='x', zorder=zorder)
    
                
    def _plot(self, plot_filename, pa_request='',
              num_stars='', num_stars_request='',
              num_central_stars='', min_cut='', max_cut=''):

        fig, ax = _plt.subplots(
            subplot_kw={'position': [0.03, 0.13, 0.75, 0.75]})

        ax.set_title(
            '{}\n{} stars at ({:.5f}, {:.5f})'.format(
                self.obsmode, self.star_type.capitalize(), self.ra, self.dec))

        ax.set_xlabel('RA (deg)')
        ax.set_ylabel('Dec (deg)')

        text_str = ('num_stars:\n' +
                    '   {}\n'.format(num_stars) +
                    'num_stars_request:\n' +
                    '   {}\n\n'.format(num_stars_request))

        if self.obsmode != 'LIFU':
            self._plot_useful_region_in_plate(ax)
            self._plot_cut_region_in_plate(ax, min_cut)
            self._plot_cut_region_in_plate(ax, max_cut)

            text_str += ('num_central_stars:\n' +
                         '   {}\n'.format(num_central_stars) +
                         'min_cut:\n' +
                         '   {}\n'.format(min_cut) +
                         'max_cut:\n' +
                         '   {}'.format(max_cut))
        else:
            self._plot_useful_region_for_guide_cam(ax)
            self._plot_guide_cam_fov(ax, pa_request,
                                     linestyle=':', color='blue')
            self._plot_guide_cam_fov(ax, self.selected_pa,
                                     linestyle='--', color='k')

            if type(self.selected_pa) == float:
                selected_pa_str = '{:.3f}'.format(self.selected_pa)
            else:
                selected_pa_str = '{}'.format(self.selected_pa)

            if type(pa_request) == float:
                pa_request_str = '{:.3f}'.format(pa_request)
            else:
                pa_request_str = '{}'.format(pa_request)

            text_str += ('selected_pa:\n' +
                         '   {} deg\n'.format(selected_pa_str) +
                         'pa_request:\n' +
                         '   {} deg'.format(pa_request_str))

        ax.text(1.05, 1.0, text_str, transform=ax.transAxes, fontsize=12,
                va='top', ha='left')

        self._plot_table(self.full_table, ax, 'red')
        self._plot_table(self.useful_table, ax, 'darkorange')
        self._plot_table(self.selected_table, ax, 'green')

        ax.scatter([self.ra], [self.dec], c='k', marker='+')

        ax.set_xlim([self.ra_max, self.ra_min])
        ax.set_ylim([self.dec_min, self.dec_max])
        ax.set_aspect(1. / _np.cos(_np.deg2rad(self.dec)))

        fig.savefig(plot_filename)

                
    def get_table(self, selection='filter', verbose=True, plot_filename=None,
                  pa_request=None, num_stars_request=None,
                  num_central_stars=1, min_cut=0.9, max_cut=1.0):

        # Assert the input values

        assert selection in ['filter', 'useful', 'full']

        if (num_stars_request) is not None:
            assert num_stars_request > 0
        
        if pa_request not in [None, _np.nan]:
            assert (pa_request >= 0.0) or (pa_request < 360.0)

            if self.obsmode != 'LIFU':

                if pa_request != 0.0:
                    logging.error('{} obsmode requires PA=0.0'.format(
                        self.obsmode))
                    raise ValueError

        # Set the needed tables

        self._set_full_table()

        if selection in ['useful', 'filter']:

            self._set_useful_table()

            if selection == 'filter':
                self._set_selected_table(num_stars_request,
                                         pa_request=pa_request,
                                         num_central_stars=num_central_stars,
                                         min_cut=min_cut, max_cut=max_cut)

        # Choose the output depending on the selection method

        if selection == 'full':
            pa = pa_request
            table = self.full_table
        elif selection == 'useful':
            pa = pa_request
            table = self.full_table
        elif selection == 'filter':
            pa = self.selected_pa
            table = self.selected_table

        # Print a summary and plot a figure if requested

        if verbose == True:
            self._print_summary()

        if plot_filename is not None:

            num_stars = len(table)

            self._plot(plot_filename, pa_request=pa_request,
                       num_stars=num_stars, num_stars_request=num_stars_request,
                       num_central_stars=num_central_stars,
                       min_cut=min_cut, max_cut=max_cut)

        return pa, table

    
class GuideStars(_AuxStars):
    """
    Handle the retrieval of WEAVE guide stars.
    
    This class provides a mechanism for querying and retrieving guide star
    targets for the construction of WEAVE 'protofields'.
    
    Parameters
    ----------
    ra : float
        The right ascension (decimal degrees) of either the central spaxel or
        central FoV.
    dec : float
        The declination (decimal degrees) of either the central spaxel or
        central FoV.
    obsmode : str
        Either LIFU, mIFU or MOS.
    dither_size_arcsec : float, optional
        Maximum allowed size for dithering patterns.
    nside : int, optional
        Override the default HEALPix nside value. Will likely end in tears.
    """

    def __init__(self, ra, dec, obsmode, dither_size_arcsec=10.0, nside=32):

        super().__init__(ra, dec, obsmode,
                         dither_size_arcsec=dither_size_arcsec, nside=nside)

        self.star_type = 'guide'


    def _get_cat(self, healpix_index, directory):

        file_path = _get_guide_cat(healpix_index, directory)

        return file_path


    def get_table(self, selection='filter', verbose=True, plot_filename=None,
                  pa_request=None, num_stars_request=None,
                  num_central_stars=1, min_cut=0.9, max_cut=1.0):
        """
        Get a table with a set of selected stars with the requested conditions.

        Parameters
        ----------
        selection : {'filter', 'useful', 'full'}
            Selection to be made in the output table: 'full' returns all the
            stars near to the coordinates, 'useful' returns all the stars which
            could be used for guiding at that coordinates (ignoring any request
            for the position angle), 'filter' returns a table following the
            prescriptions given in the other keywords of this method.
        verbose : bool
            Print a summary or not.
        plot_filename : str
            Name of the file to save a figure with the full table, the useful
            table and the selected table.
        pa_request : float, optional
            The position angle (degrees) of rotation (it can be different to
            zero only for LIFU). np.nan and None are also valid values.
        num_stars_request : int, optional
            Maximum number of guide stars in the output. None means no limit.
        num_central_stars : int, optional
            Number of stars near to centre to be selected (only used when
            obsmode is not LIFU).
        min_cut : float, optional
            Minimum cut factor to be used for the non-central stars (only used
            when obsmode is not LIFU).
        max_cut : float, optional
            Maximum cut factor to be used for the non-central stars (only used
            when obsmode is not LIFU).

        Returns
        -------
        pa : float
            The selected position angle.
        table : astropy.table.Table
            A table containing the guide stars.
        """

        pa, table = super().get_table(selection=selection, verbose=verbose,
                                      plot_filename=plot_filename,
                                      pa_request=pa_request,
                                      num_stars_request=num_stars_request,
                                      num_central_stars=num_central_stars,
                                      min_cut=min_cut, max_cut=max_cut)

        return pa, table


class CalibStars(_AuxStars):
    """
    Handle the retrieval of WEAVE calib stars.
    
    This class provides a mechanism for querying and retrieving calib star
    targets for the construction of WEAVE 'protofields'.
    
    Parameters
    ----------
    ra : float
        The right ascension (decimal degrees) of either the central spaxel or
        central FoV.
    dec : float
        The declination (decimal degrees) of either the central spaxel or
        central FoV.
    obsmode : str
        Either LIFU, mIFU or MOS.
    dither_size_arcsec : float, optional
        Maximum allowed size for dithering patterns.
    nside : int, optional
        Override the default HEALPix nside value. Will likely end in tears.
    """


    def __init__(self, ra, dec, obsmode, dither_size_arcsec=10.0, nside=32):

        assert obsmode != 'LIFU'

        super().__init__(ra, dec, obsmode,
                         dither_size_arcsec=dither_size_arcsec, nside=nside)

        self.star_type = 'calib'


    def _get_cat(self, healpix_index, directory):

        file_path = _get_calib_cat(healpix_index, directory)

        return file_path


    def get_table(self, selection='filter', verbose=True, plot_filename=None,
                  num_stars_request=None,
                  num_central_stars=0, min_cut=0.2, max_cut=0.4):
        """
        Get a table with a set of selected stars with the requested conditions.

        Parameters
        ----------
        selection : {'filter', 'useful', 'full'}
            Selection to be made in the output table: 'full' returns all the
            stars near to the coordinates, 'useful' returns all the stars which
            which are in the field of view at that coordinates, 'filter' returns
            a table following the prescriptions given in the other keywords of
            this method.
        verbose : bool
            Print a summary or not.
        plot_filename : str
            Name of the file to save a figure with the full table, the useful
            table and the selected table.
        num_stars_request : int, optional
            Maximum number of guide stars in the output. None means no limit.
        num_central_stars : int, optional
            Number of stars near to centre to be selected.
        min_cut : float, optional
            Minimum cut factor to be used for the non-central stars.
        max_cut : float, optional
            Maximum cut factor to be used for the non-central stars.

        Returns
        -------
        table : astropy.table.Table
            A table containing the guide stars.
        """

        pa, table = super().get_table(selection=selection, verbose=verbose,
                                      plot_filename=plot_filename,
                                      pa_request=None,
                                      num_stars_request=num_stars_request,
                                      num_central_stars=num_central_stars,
                                      min_cut=min_cut, max_cut=max_cut)

        return table


if __name__ == '__main__':

    # Get a parser

    import argparse

    parser = argparse.ArgumentParser(
        description='Get WEAVE guide and calib stars')

    parser.add_argument('type_stars', choices=['guide', 'calib'],
                        help='type of stars to be searched')
    parser.add_argument('obsmode', choices=['MOS', 'LIFU', 'mIFU'],
                        help='obsmode to be considered in the search')
    parser.add_argument('central_ra', type=float,
                        help='RA in degrees of the centre of the FoV')
    parser.add_argument('central_dec', type=float,
                        help='Dec in degrees of the centre of the FoV')
    parser.add_argument('--pa_request', default='None',
                        help="""Requested PA in degress for LIFU observations;
                        None means no limit""")
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
    parser.add_argument('--num_calib_stars_request', default=2,
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
    parser.add_argument('--no_verbose', dest='verbose', action='store_false',
                        help='deactivate verbose output')

    # Get the arguments from the parser

    args = parser.parse_args()

    # In arguments which accept None value, do casting or transform None values
    
    if args.pa_request != 'None':
        pa_request = float(args.pa_request)
    else:
        pa_request = None

    if args.lifu_num_guide_stars_request != 'None':
        lifu_num_guide_stars_request = int(args.lifu_num_guide_stars_request)
    else:
        lifu_num_guide_stars_request = None
    
    if args.mifu_num_guide_stars_request != 'None':
        mifu_num_guide_stars_request = int(args.mifu_num_guide_stars_request)
    else:
        mifu_num_guide_stars_request = None
    
    if args.num_calib_stars_request != 'None':
        num_calib_stars_request = int(args.num_calib_stars_request)
    else:
        num_calib_stars_request = None

    # Set the logging level
    
    logging.basicConfig(level=logging.INFO)
    
    # Choose a filename for the output plot
    
    if (args.obsmode == 'LIFU') and (pa_request is not None):
    
        if pa_request is not None:
            pa_request_str = '{:.3f}'.format(pa_request)
        else:
            pa_request_str = 'none'
    
        plot_filename = '{}s_lifu_{:.5f}_{:.5f}_{}.png'.format(
            args.type_stars, args.central_ra, args.central_dec, pa_request_str)
    else:
        plot_filename = '{}s_{:.5f}_{:.5f}.png'.format(
            args.type_stars, args.central_ra, args.central_dec)

    # Get the requested typo of stars

    if args.type_stars == 'guide':

        if args.obsmode == 'LIFU':
            num_stars_request = args.lifu_num_guide_stars_request
        else:
            num_stars_request = args.mifu_num_guide_stars_request
    
        guide_stars = GuideStars(args.central_ra, args.central_dec,
                                 args.obsmode)

        actual_pa, stars_table = guide_stars.get_table(
            verbose=args.verbose, plot_filename=plot_filename,
            pa_request=pa_request, num_stars_request=num_stars_request,
            num_central_stars=args.mifu_num_central_guide_stars,
            min_cut=args.mifu_min_guide_cut, max_cut=args.mifu_max_guide_cut)

    elif args.type_stars == 'calib':

        assert obsmode != 'LIFU'

        calib_stars = CalibStars(args.central_ra, args.central_dec,
                                 args.obsmode)

        stars_table = calib_stars.get_table(
            verbose=args.verbose, plot_filename=plot_filename,
            num_stars_request=args.num_calib_stars_request,
            num_central_stars=args.num_central_calib_stars,
            min_cut=args.min_calib_cut, max_cut=args.max_calib_cut)

    else:

        raise ValueError
    
    # Print the most interesting columns of the returned table
    
    stars_table.keep_columns(
        ['CNAME', 'GAIA_RA', 'GAIA_DEC', 'DISTANCE', 'ANGLE'])
    
    for line in str(stars_table).split('\n'):
        logging.info(line)
    
    # Report the plot file created
    
    logging.info('Plot available at file {}'.format(plot_filename))

