#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download PS1 images and add broadband photometry to an existing IFU catalogue.

For every source in the IFU catalogue, download Pan-STARRS (g, r, i) images
and fill the corresponding columns with circular aperture photometric
measurements within each fibre.

Created on Wed Nov  4 17:33:49 2020

@author: yago.ascasibar@uam.es, sctrager@astro.rug.nl
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

PS1_PIXEL = 0.25  # arcsec
LIFU_FIBRE_RADIUS = 1.305  # arcsec (2.6 arcsec diameter)
LIFU_FIBRE_AREA = np.pi*LIFU_FIBRE_RADIUS**2


# <codecell> Code from PanSTARRS https://ps1images.stsci.edu/ps1image.html

def getimages(ra, dec, size=240, filters="grizy"):
    """Query ps1filenames.py service to get a list of images.

    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table


def geturl(ra, dec, size=240, output_size=None, filters="grizy",
           format="jpg", color=False):
    """
    Get URL for images in the table.

    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URL for single-filter grayscale images.
    Returns a string with the URL
    """
    if color and format == "fits":
        raise ValueError(
            "color images are available only for jpg or png formats")
    if format not in ("jpg", "png", "fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra, dec, size=size, filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0, len(table)//2, len(table)-1]]
        for i, param in enumerate(["red", "green", "blue"]):
            url = url + "&{}={}".format(param, table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url


# <codecell> Specific WEAVE (not necessarily LIFU?) code

def photometry_from_PS1(cat_filename):
    """
    Download broadband imaging for a given catalogue.

    Creates a directory where one file per band is saved for every
    unique TARGID in the catalog.

    Parameters
    ----------
    cat_filename : str
        A FITS file with an IFU catalogue.
    """
    logging.info(
        'Downloading PS1 broadband imaging for {}'.format(cat_filename))

    output_dir = os.path.splitext(cat_filename)[0]+'_PS1'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        logging.info('Creating {}'.format(output_dir))

    with fits.open(cat_filename) as hdu_list:

        # Get the data from the catalogue

        data = hdu_list[1].data

        for galaxy in np.unique(data['TARGID']):
            rows = np.where(data['TARGID'] == galaxy)
            ra_Gaia = data['GAIA_RA'][rows]
            dec_Gaia = data['GAIA_DEC'][rows]
            coords_Gaia = SkyCoord(ra=ra_Gaia, dec=dec_Gaia,
                                   unit='deg', frame=ICRS())
            coords_PS1 = coords_Gaia.transform_to(FK5(equinox='J2000'))
            ra = coords_PS1.ra.value
            dec = coords_PS1.dec.value

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
            size_pix = int(size_deg * 3600. / PS1_PIXEL)
            logging.info('{}: center=({}, {}); size=({} deg, {} pix)'.format(
                galaxy, center_ra, center_dec, size_deg, size_pix))

            # Retrieve images:
            for band in "gri":
                filename = '{}_{}.fits'.format(galaxy, band)
                full_path = os.path.join(output_dir, filename)
                if os.path.isfile(full_path):
                    logging.info('{} already exists'.format(filename))
                    hdu = fits.open(full_path)
                else:
                    logging.info('Downloading {}'.format(filename))
                    fitsurl = geturl(center_ra, center_dec, size=size_pix,
                                     filters=band, format="fits")
                    for url in fitsurl:
                        band = url[-13]
                        # print(filename, band)
                        hdu = fits.open(url)
                        hdu.writeto(full_path, overwrite=True)

                zero_point = 25 + 2.5*np.log10(hdu[0].header['exptime'])
                wcs_info = WCS(hdu[0].header)
                img = hdu[0].data

                coords_pix = utils.skycoord_to_pixel(coords_PS1, wcs_info,
                                                     origin=0, mode='all')
                apertures = CircularAperture(np.array(coords_pix).transpose(),
                                             r=LIFU_FIBRE_RADIUS/PS1_PIXEL)
                phot_table = aperture_photometry(img, apertures)
                # good_photo = np.where(phot_table['aperture_sum'] > 0)
                bad_photo = np.where(phot_table['aperture_sum'] <= 0)
                photo = zero_point - 2.5*np.log10(
                    phot_table['aperture_sum'].data)

                # Update table:

                colname = 'MAG_'+band.upper()
                hdu_list[1].data[colname][rows] = photo
                # Scott's addition, 28.11.2020:
                # set core spaxels with r>=25. to TARGUSE="S" -- 
                # anything this faint in PanSTARRS imaging
                # is indistinguishable from sky
                if band=='r':
                    logging.info('Changing TARGUSE, TARGCLASS, APS keys')
                    targuse=hdu_list[1].data['TARGUSE'][rows]
                    newtarguse=np.where(np.less_equal(hdu_list[1].data[colname][rows],25.) & np.char.equal(targuse,'T'), 'T', 'S')
                    hdu_list[1].data['TARGUSE'][rows]=newtarguse
                    targclass=np.where(np.char.equal(newtarguse,'T') & np.char.not_equal(hdu_list[1].data['IFU_SPAXEL'][0],'S'),'GALAXY','SKY')
                    hdu_list[1].data['TARGCLASS'][rows]=targclass
                    aps_keys=['APS_WL_MIN','APS_WL_MAX','APS_Z','APS_SIGMA','APS_IFU_TSSL_TARG_SNR']
                    aps_defaults=[3700.,9550.,0.09,200.,20.]
                    for i,aps_key in enumerate(aps_keys):
                        aps_default=np.where(np.char.equal(newtarguse,'T') & \
                                             np.char.not_equal(hdu_list[1].data['IFU_SPAXEL'][0],'S'), \
                                             aps_defaults[i],None)
                        hdu_list[1].data[aps_key][rows]=aps_default
                # xx = hdu_list[1].data[colname]
                # print('\n\n')
                # print(len(rows[0]), rows)
                # # print(len(good_photo[0]), good_photo)
                # print(len(photo), photo)
                # # print(len(xx), len(xx[rows][good_photo]), xx[rows][good_photo])
                # print(xx[rows][0:10])

                # Diagnostic plots:

                # mu=zero_point-2.5*np.log10(img*LIFU_FIBRE_AREA/PS1_PIXEL**2)
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
        logging.info('Saving'+new_filename)
        hdu_list.writeto(new_filename, overwrite=True)

        # g = hdu_list[1].data['MAG_G']
        # print(len(rows[0]), rows)
        # print(len(good_photo[0]), good_photo)
        # print(len(photo), photo)
        # print(len(g), len(g[rows][good_photo]), g[rows][good_photo])


# photometry_from_PS1('../output/WA_2020A1-ifu_from_xmls.fits')


# <codecell> When called from the command line

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""Add PS1 photometry to FITS catalogue.""")

    parser.add_argument('catalogue',
                        help='a FITS file with an IFU driver catalogue')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    photometry_from_PS1(args.catalogue)


# <codecell> Bye
# -----------------------------------------------------------------------------
#                                                           ... Paranoy@ Rulz!
# -----------------------------------------------------------------------------
