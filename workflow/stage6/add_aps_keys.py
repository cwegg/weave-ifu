#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Add APS_ parameters to IFU FITS catalogue. Assumes each target has a single set of
APS_ parameters and ignores all fibres marked as TARGUSE="S".

@author: sctrager@astro.rug.nl
"""

import argparse
import logging
import numpy as np
from astropy.io import fits

def add_aps_keys(cat_filename, aps_dict):
    """
    Add APS_ parameters for a given catalogue.

    Parameters
    ----------
    cat_filename : str
        A FITS file with an IFU catalogue.
        
    aps_dict : dict
        A dictionary containing {aps_key : value} flags
    """
    logging.info(
        'Adding APS_ parameters to {}'.format(cat_filename))

    with fits.open(cat_filename, mode='update') as hdu_list:

        # get the data from the catalogue

        for i,galaxy in enumerate(np.unique(hdu_list[1].data['TARGID'])):
            rows = np.where(hdu_list[1].data['TARGID'] == galaxy)
            logging.info('Modifying APS_ values for {}'.format(galaxy))

            targclass=hdu_list[1].data['TARGCLASS'][rows]
            for aps_key in aps_dict.keys():
                aps_default=np.where(np.char.equal(targclass,'GALAXY'),
                                     aps_dict[aps_key][i],None)
                hdu_list[1].data[aps_key][rows]=aps_default

        logging.info('Saving {}'.format(cat_filename))
        hdu_list.flush()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""Add APS_ defaults to IFU FITS catalogue.""")

    parser.add_argument('catalogue',
                        help='a FITS file with an IFU driver catalogue')

    parser.add_argument('aps_dict',
                        help='a dictionary of APS_ keywords and default values')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    add_aps_keys(args.catalogue, args.aps_dict)
