#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Add broadband photometry from SDSS (g,r,i) imaging to an existing IFU catalogue.

For every source in the IFU catalogue, fill the corresponding columns
with circular aperture photometric measurements within each
fibre. Note that unlike the PanSTARRS1 equivalent code by Yago, this
code does NOT download SDSS mosaics, as the web interface for these is
quite clunky. See https://dr12.sdss.org/mosaics for instructions and
access to these mosaics.

Created on Tue 1 Dec 2020 09:35

@authors: yago.ascasibar@uam.es, sctrager@astro.rug.nl

"""

import argparse
import logging
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, FK5, ICRS
from astropy.wcs import WCS, utils
from photutils import CircularAperture, aperture_photometry

SDSS_PIXEL = 0.396  # arcsec
LIFU_FIBRE_RADIUS = 1.305  # arcsec (2.61 arcsec diameter)
LIFU_FIBRE_AREA = np.pi*LIFU_FIBRE_RADIUS**2

# <codecell> following is taken from https://github.com/esheldon/sdsspy/blob/master/sdsspy/util.py

_bands = ['u','g','r','i','z']
_bvalues=[1.4, 0.9, 1.2, 1.8, 7.4]
_log10 = np.log(10.0)

def _nmgy2lups_1band(nmgy, band):
    # luptitudes = asinh mags
    b=_bvalues[band]
    lups = 2.5*(10.0-np.log10(b)) - 2.5*np.arcsinh(5.0*nmgy/b)/_log10
    # pogson magnitudes
    #lups = 22.5 - 2.5 * np.log10(nmgy)
    return lups

# <codecell> Specific WEAVE (not necessarily LIFU?) code

def photometry_from_SDSS(cat_filename,input_dir='output'):
    """
    Creates a directory where one file per band is saved for every
    unique TARGID in the catalog.

    Parameters
    ----------
    cat_filename : str
        A FITS file with an IFU catalogue.
    """
    logging.info(
        'Reading SDSS broadmand imaging for {}'.format(cat_filename))

    with fits.open(cat_filename) as hdu_list:

        # Get the data from the catalogue

        data = hdu_list[1].data

        for galaxy in np.unique(data['TARGID']):
            rows = np.where(data['TARGID'] == galaxy)
            ra_Gaia = data['GAIA_RA'][rows]
            dec_Gaia = data['GAIA_DEC'][rows]
            coords_Gaia = SkyCoord(ra=ra_Gaia, dec=dec_Gaia,
                                   unit='deg', frame=ICRS())
            coords_SDSS = coords_Gaia.transform_to(FK5(equinox='J2000'))
            ra = coords_SDSS.ra.value
            dec = coords_SDSS.dec.value

            # Center and size of the images:
            ra_min = np.min(ra)
            ra_max = np.max(ra)
            dec_min = np.min(dec)
            dec_max = np.max(dec)
            center_ra = (ra_min + ra_max)/2
            center_dec = (dec_min + dec_max)/2
            delta_ra = ra_max - ra_min
            delta_dec = dec_max - dec_min
            size_deg = np.max([delta_ra*np.sin(np.radians(center_dec)),
                               delta_dec]) * 1.1  # add 10% buffer
            size_pix = int(size_deg * 3600. / SDSS_PIXEL)
            logging.info('{}: center=({}, {}); size=({} deg, {} pix)'.format(
                galaxy, center_ra, center_dec, size_deg, size_pix))

            # Photometer images:
            for band in "gri":
                filename = '{}_{}.fits'.format(galaxy, band)
                full_path = os.path.join(input_dir, filename)
                if os.path.isfile(full_path):
                    logging.info('{} exists - good'.format(filename))
                    hdu = fits.open(full_path)
                else:
                    logging.info('Please ensure that SDSS fits file with name {} exits'.format(full_path))
                    sys.exit()
 
                wcs_info = WCS(hdu[0].header)
                img = hdu[0].data

                coords_pix = utils.skycoord_to_pixel(coords_SDSS, wcs_info,
                                                     origin=0, mode='all')
                apertures = CircularAperture(np.array(coords_pix).transpose(),
                                             r=LIFU_FIBRE_RADIUS/SDSS_PIXEL)
                phot_table = aperture_photometry(img, apertures)
                bad_photo = np.where(phot_table['aperture_sum'] <= 0)
                photo = _nmgy2lups_1band(phot_table['aperture_sum'].data,
                                         _bands.index(band))

                # Update table:

                colname = 'MAG_'+band.upper()
                hdu_list[1].data[colname][rows] = photo
                
                # Scott's addition, 28.11.2020:
                # set core spaxels with r>=25. to TARGUSE="S" -- 
                # anything this faint in PanSTARRS
                # imaging is likely to be sky
#                if band=='r':
#                    logging.info('Changing TARGUSE')
#                    targuse=hdu_list[1].data['TARGUSE'][rows]
#                    newtarguse=np.where(np.less_equal(hdu_list[1].data[colname][rows],25.) & np.char.equal(targuse,'T'), 'T', 'S')
#                    hdu_list[1].data['TARGUSE'][rows]=newtarguse

                # xx = hdu_list[1].data[colname]
                # print('\n\n')
                # print(len(rows[0]), rows)
                # # print(len(good_photo[0]), good_photo)
                # print(len(photo), photo)
                # # print(len(xx), len(xx[rows][good_photo]), xx[rows][good_photo])
                # print(xx[rows][0:10])

                # Diagnostic plots:

                # mu=zero_point-2.5*np.log10(img*LIFU_FIBRE_AREA/SDSS_PIXEL**2)
                # # plt.figure()
                # plt.figure(figsize=(8, 8))
                # plt.imshow(mu, origin='lower',
                #            vmin=16, vmax=24, cmap='inferno_r')
                # # plt.xlim(size_pix/2-100, size_pix/2+100)
                # # plt.ylim(size_pix/2-100, size_pix/2+100)
                # plt.colorbar()
                # plt.show()
                # plt.figure(figsize=(8, 8))
                # # plt.imshow(mu, origin='lower',
                #             # vmin=0, vmax=2, cmap='inferno_r')
                # # plt.scatter(coords_pix[0], coords_pix[1], s=10)
                # plt.scatter(coords_pix[0][good_photo],
                #             coords_pix[1][good_photo], s=100,
                #             c=photo, vmin=16, vmax=24, cmap='inferno_r')
                # plt.xlim(size_pix/2-150, size_pix/2+150)
                # plt.ylim(size_pix/2-150, size_pix/2+150)
                # plt.colorbar().set_label('fibre '+band+' mag')
                # plt.title(filename)
                # plt.show()
                fig = plt.figure(figsize=(12, 6))
                ax = fig.add_axes([.1, .1, .4, .8])  # , projection=wcs_info)
                sc = ax.scatter(ra, dec, s=4, c=photo,
                                vmin=16, vmax=24, cmap='inferno_r')
                ax.scatter(ra[bad_photo], dec[bad_photo], s=4, c='grey')
                ax.grid('both')
                ax.set_xlabel('RA (J2000)')
                ax.set_ylabel('DEC (J2000)')
                ax = fig.add_axes([.5, 0, .5, 1])  # , projection=wcs_info)
                ax.axis('off')
                sc = ax.scatter(ra, dec, s=40, c=photo,
                                vmin=16, vmax=24, cmap='inferno_r')
                ax.scatter(ra[bad_photo], dec[bad_photo], s=40, c='grey')
                ax.set_xlim(center_ra-.02, center_ra+.02)
                ax.set_ylim(center_dec-.02, center_dec+.02)
                fig.suptitle(filename)
                ax = fig.add_axes([.55, .1, .4, .025])
                cb = fig.colorbar(sc, orientation="horizontal", cax=ax)
                cb.set_label('fibre '+band+' mag')
                plt.savefig(full_path[:-4]+'png')
                # plt.show()
                plt.close()

        # Save a new table with the updated columns

        new_filename = cat_filename[:-15]+'.fits'
        logging.info('Saving '+new_filename)
        hdu_list.writeto(new_filename, overwrite=True)

        # g = hdu_list[1].data['MAG_G']
        # print(len(rows[0]), rows)
        # print(len(good_photo[0]), good_photo)
        # print(len(photo), photo)
        # print(len(g), len(g[rows][good_photo]), g[rows][good_photo])


# photometry_from_SDSS('../output/WA_2020A1-ifu_from_xmls.fits')


# <codecell> When called from the command line

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""Add SDSS photometry to FITS catalogue.""")

    parser.add_argument('catalogue',
                        help='a FITS file with an IFU driver catalogue')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    parser.add_argument('--input_dir',default='output',help='Directory where SDSS mosaic images correspond to the individual objects can be found')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    photometry_from_SDSS(args.catalogue, args.input_dir)


# <codecell> Bye
# -----------------------------------------------------------------------------
#                                                           ... Paranoy@ Rulz!
# -----------------------------------------------------------------------------
