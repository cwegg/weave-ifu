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
Retrieve of WEAVE calib stars for the mIFU or MOS.

The authors of this module are:

- David Murphy (dmurphy@ast.cam.ac.uk), Cambridge Astronomical Survey Unit
  (CASU, IoA).

Notes
-----
The dependencies of this module are:

- Python core
- numpy
- astropy
- astropy_healpix

Log:

- v0.1: Still work in progress.
"""


import os
from copy import deepcopy
from operator import indexOf
from urllib.request import urlretrieve

import numpy as np
from astropy import coordinates
from astropy.io import fits
from astropy.table import Table
from astropy_healpix import HEALPix


class CalibStars:
    """
    Handle the retrieval of WEAVE calib stars for the mIFU or MOS.
    
    This class provides a mechanism for querying and retrieving calibration star
    targets for the construction of WEAVE 'protofields'.
    
    Parameters
    ----------
    ra : float
        The Right Ascension (decimal degrees) - of either the central spaxel or
        central FOV.
    dec : float
        The Declination (decimal degrees) - of either the central spaxel or
        central FOV.
    pa : float
        The position angle (degrees) of rotation (disallowed for calibration
        selection).
    obsmode : str
        Either MOS or mIFU.
    nside : int, optional
        Override the default HEALPix nside value. Will likely end in tears.
    annular : bool, optional
        Search only within an annular radius at the edge of the WEAVE FOV.
    plot : bool, optional
        Plot the results?
    """


    def __init__(self,ra,dec,pa,obsmode,nside=32,annular=False,plot=False):
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.obsmode = obsmode
        self.nside = nside
        self.annular = annular
        self.plot = plot
        
        self.calibs = None
        self.calib_url = "http://casu.ast.cam.ac.uk/~dmurphy/opr3/swg/resources/calibs/"
        self.calib_filename = "WD_S4_<HEALPIX>.fits"
        self.cache_dir = './aux/'
        
        if not self.obsmode in ['MOS','mIFU']:
            print('ERROR: must specify MOS or mIFU obsmode')
            raise SystemExit(0)
        
        if self.pa != 0.0:
            print('ERROR: CalibStar requires PA=0')
            raise SystemExit(0)

        if not os.path.isdir(self.cache_dir):
            os.makedirs(os.path.abspath(self.cache_dir), exist_ok=True)
        
        
    def _set_geometry(self,healpix=True):

        self.ra_min = self.ra - 1.0
        self.ra_max = self.ra + 1.0
        self.dec_min = self.dec - 1.0
        self.dec_max = self.dec + 1.0

        if healpix:
            hp = HEALPix(nside=self.nside, order='nested', frame=coordinates.ICRS())
            self.healpix_indices = []
            ra = [self.ra_min,self.ra_min,self.ra_max,self.ra_max]
            dec = [self.dec_min,self.dec_max,self.dec_min,self.dec_max]

            dx = self.ra_max - self.ra_min
            dy = self.dec_max - self.dec_min

            for i in np.arange(500):
                ra.append(self.ra_min+(np.random.random()*dx))
                dec.append(self.dec_min+(np.random.random()*dy))

            for r,d in zip(ra,dec):
                self.healpix_indices.append(hp.skycoord_to_healpix(coordinates.SkyCoord(r,d,unit='deg')))
            self.healpix_indices = np.unique(self.healpix_indices)


    def _retrieve_calibcats(self,clobber=False):
        if len(self.healpix_indices) == 0:
            self._set_geometry()
            
        self.calib_files = []

        for hp in self.healpix_indices:
            fn = self.calib_filename.replace('<HEALPIX>',str(hp))
            url = self.calib_url+fn
            fn = self.cache_dir+fn
            if os.path.isfile(fn):
                if clobber == False:
                    print('Using existing file %s'%(fn))
                    self.calib_files.append(fn)
                    continue
            print(url)
            urlretrieve(url, fn)
            self.calib_files.append(fn)
            
        
        tabs = []
        for cat in self.calib_files:
            hdu_list = fits.open(cat)
            tabs.append(Table(hdu_list[1].data))
            
        if len(tabs) == 1:
            self.calibs = tabs[0]
            
        else:
            t0 = tabs[0]
            for t in tabs[1:]:
                for g in t:
                    t0.add_row(g)
            self.calibs = t0


    def _select_target(self):
        self.ra_c = self.calibs['GAIA_RA']
        self.dec_c = self.calibs['GAIA_DEC']
        
#        if (not self.lifu):
#            annular = True
            #fig = plt.figure(); sp = plt.subplot(aspect='equal'); plt.plot(self.ra_g,self.dec_g,'ko'); cir = plt.Circle((gs.ra,gs.dec),radius=1.0,lw=2.0,color='red'); sp.add_patch(cir); plt.show()

            
        if self.annular == False:        
            r_min = 0.0
            r_max = 1.0
            
        else:
            r_min = 0.95
            r_max = 1.0

                
        self.radii = np.array([(((_ra-self.ra)**2)+((_dec-self.dec)**2))**0.5 for _ra,_dec in zip(self.ra_c,self.dec_c)])
        filter = (self.radii > r_min) & (self.radii < r_max)
        self.calibs_filter = self.calibs[filter]
        self.ra_x = self.calibs_filter['GAIA_RA']
        self.dec_x = self.calibs_filter['GAIA_DEC']

        radii = self.radii[filter]

        r0 = r_min + (0.5*(r_max-r_min))
        self.dist = [abs(d-r0) for d in radii]

        self.dist = radii
        minval = min(self.dist)

        #do a quick report
        dist_sort = deepcopy(self.dist)
        dist_sort.sort()
        if self.annular:
            print('Annular search summary:')
        else:
            print('Search summary:')
        print("#\t Dist (')  CNAME\t\t RA\t    Dec\t      angle\t Gaia_G mag")
        i = 0
        self.calibs_filter['dist'] = self.dist

        for d in dist_sort:
            i = i + 1
            index = indexOf(self.dist,d)
            calib_candidate = self.calibs_filter[index]
            ra_trans = self.ra - calib_candidate['GAIA_RA']
            dec_trans = self.dec - calib_candidate['GAIA_DEC']
            ang = np.arctan2(dec_trans,ra_trans)
            ang = (ang*180.0) / np.pi
            if ang < 0:
                ang += 360.0
            print('#%d\t %1.2f\t   %s\t %1.4f   %1.4f   %1.3f\t %1.3f'%(i,d*60.0,calib_candidate['CNAME'],calib_candidate['GAIA_RA'],calib_candidate['GAIA_DEC'],ang,calib_candidate['GAIA_MAG_GG']))

        self.calibs_filter.sort('dist')
            
        if self.plot:
            self.ra_c = self.calibs_filter['GAIA_RA']
            self.dec_c = self.calibs_filter['GAIA_DEC']
            fig = plt.figure(); sp = plt.subplot(aspect='equal'); plt.plot(self.ra_c,self.dec_c,'bo'); cir = plt.Circle((self.ra,self.dec),radius=1.0,lw=2.0,color='red'); sp.add_patch(cir)
            if self.annular:
                cir = plt.Circle((self.ra,self.dec),radius=0.95,lw=1.0,color='blue',fill=False,ls='--')
                sp.add_patch(cir)
            plt.xlim(self.ra-1.2,self.ra+1.2)
            plt.ylim(self.dec-1.2,self.dec+1.2)
            plt.show()
        
        return self.calibs_filter
        

    def _annular_search(self):
        hp = HEALPix(nside=self.nside, order='nested', frame=coordinates.ICRS())
        
        self._set_geometry(healpix=False)
        r_min = self.ra_min-self.ra
        r_max = self.ra_max-self.ra
        radius = self.ra_max-self.ra

        #get the healpix IDs covering an annulus centred on self.ra,dec

        in_annulus = []
        self.healpix_indices = []
        print('Populating annulus and determining HEALPix coverage...')
        while(len(in_annulus)) < 500:
            rnd_ra = self.ra+(2*(np.random.random()-0.5)*radius)
            rnd_dec = self.dec+(2*(np.random.random()-0.5)*radius)
            rnd_dist = (((rnd_ra-self.ra)**2)+((rnd_dec-self.dec)**2))**0.5
            if rnd_dist > r_min:
                if rnd_dist < r_max:
                    self.healpix_indices.append(hp.skycoord_to_healpix(coordinates.SkyCoord(rnd_ra,rnd_dec,unit='deg')))
                    in_annulus.append([rnd_ra,rnd_dec])
                #print(len(in_annulus))
        print('....done')
        self.healpix_indices = np.unique(self.healpix_indices)
        print(self.healpix_indices)
        self._retrieve_calibcats()
        target = self._select_target(annular=True)
        return target

        
    def get_calib(self, annular_fail=True):
        """
        Master function to return a calib star once the object is instantiated.
        
        Parameters
        ----------
        annular_fail : bool, optional
             (for mIFU) If there are no calibrators in the annulus, search the
             whole field.
        
        Returns
        -------
        calibs : astropy.Table
             Row(s) from the WD calibration star catalogue.
        """


        self._set_geometry()
        self._retrieve_calibcats()
        calibs = self._select_target()
        
        if (type(calibs) == type(None)) and (annular_fail == True) and (self.obsmode != 'MOS'):
            print('No calibration stars found in annulus - performing full search')
            calibs = self._annular_search(annular=False)

        if type(calibs) == type(None):
            print('No calibration stars found...')
            return None
            
        return calibs

