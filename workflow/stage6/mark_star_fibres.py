#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For a given IFU FITS catalogue, query the Gaia eDR3 database for
stars near individual fibres and then mark these fibres as
TARGCLASS="STAR" and enter the GAIA_MAG_X into the catalogues.

@author: sctrager@astro.rug.nl
"""

import argparse
import logging
import numpy as np
import os
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord, ICRS
from astropy.table import Table, vstack
from astroquery.gaia import Gaia

LIFU_FIBRE_RADIUS = 1.305 * u.arcsec  # arcsec (2.61 arcsec diameter)
LIFU_LONG_AXIS = 90.39 * u.arcsec
mIFU_FIBRE_RADIUS = 1.305/2. * u.arcsec # arcsec (1.305 arcsec diameter)
mIFU_LONG_AXIS = 12.35 * u.arcsec

def mark_star_fibres(cat_filename, match_dist=2.5, ifu='LIFU',
                      parallax_over_error=3.):
    """
    For a given IFU FITS catalogue, query the Gaia eDR3 database for
    stars near individual fibres and then mark these fibres as
    TARGCLASS="STAR" and enter the GAIA_MAG_X into the catalogues.

    Parameters
    ----------
    cat_filename : str
        A FITS file containing an IFU catalogues

    match_dist : float
        The distance in (m/L)IFU fibre radii to search for nearby Gaia stars

    ifu : str
        IFU mode used (LIFU/mIFU), to determine fibre radius
    """
    logging.info(
        'Searching for Gaia eDR stars near fibres in {}'.format(cat_filename))

    with fits.open(cat_filename, mode='update') as hdu_list:

        Gaia_stars=Table()
        # loop over each target in the catalogue
        for galaxy in np.unique(hdu_list[1].data['TARGID']):
            rows=np.where(hdu_list[1].data['TARGID'] == galaxy)
            logging.info('Searching for stars in the field of {}'.format(galaxy))

            central_fibre=np.where(np.char.equal(hdu_list[1].data['TARGID'],galaxy) &
                                   np.char.equal(hdu_list[1].data['IFU_SPAXEL'],'C14'))
            central_ra=hdu_list[1].data['GAIA_RA'][central_fibre][0]
            central_dec=hdu_list[1].data['GAIA_DEC'][central_fibre][0]
            logging.info('Central fibre of {} located at {},{}'.format(galaxy,
                                                                       central_ra,
                                                                       central_dec))

            central_coord = SkyCoord(ra=central_ra, dec=central_dec,
                                     unit=(u.degree,u.degree), frame='icrs')
            # search radius for ALL stars in field: use IFU size + 25% for dithers
            # match radius is match_dist*fibre radius
            if ifu=='LIFU':
                radius = 1.25 * LIFU_LONG_AXIS / 2.
                match_radius=LIFU_FIBRE_RADIUS * match_dist
            elif ifu=='mIFU':
                radius = 1.25 * mIFU_LONG_AXIS / 2.
                match_radius=mIFU_FIBRE_RADIUS * match_dist
            else:
                print('wrong IFU mode, exiting')
                sys.exit()
            # serach Gaia eDR3
            cone_search = Gaia.cone_search_async(central_coord, radius,
                                                 table_name='gaiaedr3.gaia_source')
            Gaia_sources = cone_search.get_results()
            Gaia_stars = vstack([Gaia_stars,
                                 Gaia_sources[Gaia_sources['parallax_over_error'] >
                                              parallax_over_error]])

        logging.info(Gaia_stars.pprint())
        Gaia_stars.add_index('source_id')
        logging.info('Building coordinate masks')
        # construct star coordinates
        Gaia_star_coords = SkyCoord(ra=Gaia_stars['ra'], dec=Gaia_stars['dec'],
                                    unit=(u.degree,u.degree), frame='icrs')
        # construct fibre coordinates
        fibre_coords = SkyCoord(ra=hdu_list[1].data['GAIA_RA'],
                                dec=hdu_list[1].data['GAIA_DEC'],
                                unit=(u.degree,u.degree), frame='icrs')
        # compute separation of each fibre from each star found
        fibre_separations = np.array([fibre_coords.separation(Gaia_star_coord)
                                      for Gaia_star_coord in Gaia_star_coords]) * u.degree
        # is there a star within match_radis?
        nearby_star = np.logical_or.reduce(np.less_equal(fibre_separations,
                                                         match_radius))
        # which is the nearest star (even if it isn't within match_radius)?
        which_nearby_star = np.argmin(fibre_separations,axis=0)

        # change TARGCLASS and TARGNAME
        newtarguse = np.where(nearby_star,'T',hdu_list[1].data['TARGUSE'])
        hdu_list[1].data['TARGUSE']=newtarguse
        newtargclass = np.where(nearby_star,'STAR',hdu_list[1].data['TARGCLASS'])
        hdu_list[1].data['TARGCLASS']=newtargclass
        logging.info('Updated TARGCLASS')
        # get Gaia photometry for each star
        gaia_g = np.where(nearby_star,Gaia_sources['phot_g_mean_mag'][which_nearby_star],None)
        gaia_bp = np.where(nearby_star,Gaia_sources['phot_bp_mean_mag'][which_nearby_star],None)
        gaia_rp = np.where(nearby_star,Gaia_sources['phot_rp_mean_mag'][which_nearby_star],None)
        hdu_list[1].data['GAIA_MAG_G'] = gaia_g
        hdu_list[1].data['GAIA_MAG_BP'] = gaia_bp
        hdu_list[1].data['GAIA_MAG_RP'] = gaia_rp
        logging.info('Updated Gaia photometry')

        logging.info('Saving {}'.format(cat_filename))
        hdu_list.flush()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""Add Gaia eDR3 stars to IFU FITS catalogue.""")

    parser.add_argument('catalogue',
                        help='a FITS file with an IFU driver catalogue')

    parser.add_argument('--match_dist', default=2.5,
                        help='Star matching distance in units of IFU fibre radii')

    parser.add_argument('--ifu_mode', default='LIFU',
                        choices=['LIFU','mIFU'],
                        help='IFU mode')

    parser.add_argument('--parallax_over_error', default=3.,
                        help='Gaia sources are assumed to be real when parallax/error > parallax_over_error')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    mark_star_fibres(args.catalogue, match_dist=args.match_dist,
                      ifu=args.ifu_mode, parallax_over_error=args.parallax_over_error)
