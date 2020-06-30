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
Retrieve of WEAVE guide stars for the LIFU or MOS.

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


class GuideStars:
    """
    Handle the retrieval of WEAVE guide stars for the LIFU or MOS.
    
    This class provides a mechanism for querying and retrieving guide star
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
        The position angle (degrees) of rotation (LIFU only).
    obsmode : str
        Either LIFU or mIFU.
    nside : int, optional
        Override the default HEALPix nside value. Will likely end in tears.
    """


    def __init__(self,ra,dec,pa,obsmode,nside=32,dither_pattern=None,max_guides=None):
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.obsmode = obsmode
        self.nside = nside
        self.dither_pattern = dither_pattern
        self.max_guides = max_guides
        
        self.plate_scale = 17.8   ##  "/mm
        
        self.guides = None
        self.guide_url = "http://casu.ast.cam.ac.uk/~dmurphy/opr3/swg/resources/guides/"
        self.guide_filename = "Guides_S4_<HEALPIX>.fits"
        self.cache_dir = './cache/'

        if not self.obsmode in ['LIFU','mIFU']:
            print('ERROR: must specify LIFU or mIFU obsmode')
            raise SystemExit(0)

        if (self.obsmode != 'LIFU') and (self.pa != 0):
            print('ERROR: mIFU obsmode requires PA=0')
            raise SystemExit(0)

        if not os.path.isdir(self.cache_dir):
            os.makedirs(os.path.abspath(self.cache_dir), exist_ok=True)
        
        
    def _set_geometry(self,healpix=True):

        if (self.obsmode == 'LIFU'):
            self.g_dx = 3.75
            self.g_dy = 4.0

            #testing - let's make it bigger
    #        self.g_dx = 17.0
    #        self.g_dy = 20.0

            #needs to be about 20x larger in area
    #        self.g_dx = 16.8
    #        self.g_dy = 17.9

            self.ra_max = self.ra + ((27.7+(0.5*self.g_dx))/60.0)
            self.ra_min = self.ra + ((27.7-(0.5*self.g_dx))/60.0)
            self.dec_max = self.dec + ((0.5*self.g_dy)/60.0)
            self.dec_min = self.dec - ((0.5*self.g_dy)/60.0)

            self.ra_gc0 = self.ra + (27.7/60.0)
            self.dec_gc0 = self.dec

        else:
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


    def _retrieve_guidecats(self,clobber=False):

        if len(self.healpix_indices) == 0:
            self._set_geometry()
            
        self.guide_files = []

        for hp in self.healpix_indices:
            fn = self.guide_filename.replace('<HEALPIX>',str(hp))
            url = self.guide_url+fn
            fn = self.cache_dir+fn
            if os.path.isfile(fn):
                if clobber == False:
                    print('Using existing file %s'%(fn))
                    self.guide_files.append(fn)
                    continue
            print(url)
            urlretrieve(url, fn)
            self.guide_files.append(fn)
            
        #self.guide_files = ['/scratch/Guides_S4_cname.fits']
        tabs = []
        for cat in self.guide_files:
            hdu_list = fits.open(cat)
            tabs.append(Table(hdu_list[1].data))
            
        if len(tabs) == 1:
            self.guides = tabs[0]
            
        else:
            t0 = tabs[0]
            for t in tabs[1:]:
                for g in t:
                    t0.add_row(g)
            self.guides = t0


    def _dither_check(self,guide):
        # This is a check to ensure that, for the defined pointing and specified dither pattern,
        # the guidestar will not fall outside the guidecam FOV upon dithering

        # We will need access to configure.cfg for this (to access the dither patterns themselves)
        print('WARNING: LIFU dither pattern checks not implemented!')
        return True


    def _select_target(self,annular=False):
        self.ra_g = self.guides['GAIA_RA']
        self.dec_g = self.guides['GAIA_DEC']
        mifu_subset = []
        
        if (self.obsmode != 'LIFU'):
            annular = True
            #fig = plt.figure(); sp = plt.subplot(aspect='equal'); plt.plot(self.ra_g,self.dec_g,'ko'); cir = plt.Circle((gs.ra,gs.dec),radius=1.0,lw=2.0,color='red'); sp.add_patch(cir); plt.show()

            
        if annular ==False:        
            filter1 = (self.ra_g > self.ra_min) & (self.ra_g < self.ra_max)
            filter2 = (self.dec_g > self.dec_min) & (self.dec_g < self.dec_max)
            filter = filter1 & filter2

            if (True in filter) == False:
                print('No guide stars within the GC FOV!!')
                return None

            self.guides_filter = self.guides[np.where(filter)[0]]
            self.ra_g = self.guides_filter['GAIA_RA']
            self.dec_g = self.guides_filter['GAIA_DEC']
        
            self.dist = [((abs(r-self.ra_gc0)**2)+(abs(d-self.dec_gc0)**2))**0.5 for r,d in zip(self.ra_g,self.dec_g)]
        else:
            if (self.obsmode == 'LIFU'):
                #in annulus, want closest to the central radius of the annulus
                r_min = self.ra_min-self.ra
                r_max = self.ra_max-self.ra
            else:
                r_min = 0.95
                r_max = 1.0

                
            self.radii = np.array([(((_ra-self.ra)**2)+((_dec-self.dec)**2))**0.5 for _ra,_dec in zip(self.ra_g,self.dec_g)])
            filter = (self.radii > r_min) & (self.radii < r_max)
            self.guides_filter = self.guides[filter]
            self.ra_g = self.guides_filter['GAIA_RA']
            self.dec_g = self.guides_filter['GAIA_DEC']

            radii = self.radii[filter]
            if (self.obsmode == 'LIFU'):
                r0 = r_min + (0.5*(r_max-r_min))
                self.dist = [abs(d-r0) for d in radii]
            else:
                self.dist = radii
                
        minval = min(self.dist)
        g_index = indexOf(self.dist,minval)

        guide_sel = self.guides_filter[g_index]
        if (annular == False):
            if (guide_sel['GAIA_RA'] > self.ra_min) and (guide_sel['GAIA_RA'] < self.ra_max) and (guide_sel['GAIA_DEC'] > self.dec_min) and (guide_sel['GAIA_DEC'] < self.dec_max):
                print('Guidestar candidate is %f arcmin from GC centre'%(minval*60.0))
                # print(self.guides[g_index])
                self.guide_sel = self.guides[g_index]

                if 1:
                    print("#\t Dist (')  CNAME\t\t RA\t    Dec\t      angle\t Gaia_G mag")
                    i = 0
                    angles = []
                    for d in [minval]:
                        i = i + 1
                        sel = ''
                        # if (i == 1) and (self.obsmode == 'LIFU'):
                        #     sel = ' <-----'
                        index = 1
                        guide_candidate = self.guide_sel
                        ra_trans = self.ra - guide_candidate['GAIA_RA']
                        dec_trans = self.dec - guide_candidate['GAIA_DEC']
                        ang = np.arctan2(dec_trans,ra_trans)
                        ang = (ang*180.0) / np.pi
                        if ang < 0:
                            ang += 360.0
                        angles.append(ang)

                        print('#%d\t %1.2f\t   %s\t %1.4f   %1.4f   %1.3f\t %1.3f%s'%(i,d*60.0,guide_candidate['CNAME'],guide_candidate['GAIA_RA'],guide_candidate['GAIA_DEC'],ang,guide_candidate['GAIA_MAG_GG'],sel))
                    
                    a = 1

                
            else:
                print('Closest candidate still lies outside of the GC FOV!!')
                self.guide_sel = self.guides[g_index]
                return None

        else:
            #do a quick report
            dist_sort = deepcopy(self.dist)
            dist_sort.sort()
            angles = {}
            for d in dist_sort:
                index = indexOf(self.dist,d)
                guide_candidate = self.guides_filter[index]
                ra_trans = self.ra - guide_candidate['GAIA_RA']
                dec_trans = self.dec - guide_candidate['GAIA_DEC']
                ang = np.arctan2(dec_trans,ra_trans)
                ang = (ang*180.0) / np.pi
                if ang < 0:
                    ang += 360.0
                # angles.append(ang)
                angles[index] = ang

            _angles = [angles[ii] for ii in range(len(angles.keys()))]
            
            self.guides_filter['ANGLE'] = _angles
            


            angular_diffs = []
            
            if (self.max_guides) and (self.obsmode != 'LIFU'):
                # need to select the "best" guides up to max_guides, azimuthally distributed
                # sampling = 360.0 / float(self.max_guides)
                sampling_points = np.arange(0,360.0,360.0 / float(self.max_guides))
                for s in sampling_points:
                    diffs = [abs(s-ang) for ang in _angles]
                    minval = min(diffs)
                    if s == 0.0:
                        diffs360 = [abs(360-ang) for ang in _angles]
                        if min(diffs360) < minval:
                            diffs = diffs360
                            minval = min(diffs360)
                    jj = indexOf(diffs,minval)
                    diff_sort = deepcopy(diffs)
                    diff_lookup = {}
                    for i in range(len(diff_sort)):
                        diff_lookup[diff_sort[i]] = i
                    diff_sort.sort()
                    dist_iter = iter(diff_sort)
                    while jj in mifu_subset:
                        jj = diff_lookup[dist_iter.next()]
                    mifu_subset.append(jj)
                    angular_diffs.append(diffs[jj])
                    
                
            if (self.obsmode == 'LIFU'):
                print('Annular search summary (selected closest to centre of a rotated guidecam):')
            else:
                print('Annular search summary:')
            print("#\t Dist (')  CNAME\t\t RA\t    Dec\t      angle\t Gaia_G mag")
            i = 0
            self.guides_filter['dist'] = self.dist

            for index in range(len(self.guides_filter)):
            # for d in dist_sort:
                # i = i + 1
                sel = ''
                # index = indexOf(self.dist,d)
                if (index == 0) and (self.obsmode == 'LIFU'):
                    sel = ' <-----'
                elif (len(mifu_subset) > 0) and (index in mifu_subset):
                    sel = ' <-----'
                
                guide_candidate = self.guides_filter[index]
                # ra_trans = self.ra - guide_candidate['GAIA_RA']
                # dec_trans = self.dec - guide_candidate['GAIA_DEC']
                # ang = np.arctan2(dec_trans,ra_trans)
                # ang = (ang*180.0) / np.pi
                # if ang < 0:
                #     ang += 360.0
                # angles.append(ang)

                print('#%d\t %1.2f\t   %s\t %1.4f   %1.4f   %1.3f\t %1.3f%s'%(index+1,guide_candidate['dist']*60.0,guide_candidate['CNAME'],guide_candidate['GAIA_RA'],guide_candidate['GAIA_DEC'],guide_candidate['ANGLE'],guide_candidate['GAIA_MAG_GG'],sel))

                
            # self.guides_filter.sort('dist')

        if len(mifu_subset) > 0:
            self.guides_filter = self.guides_filter[mifu_subset]
            
        if (self.obsmode == 'LIFU'):
            if self.dither_pattern:
                dither_OK = seff._dither_check(guide_sel)
                if not dither_OK:
                    return None
            return [guide_sel]
        return self.guides_filter


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
        self._retrieve_guidecats()
        target = self._select_target(annular=True)
        return target

        
    def get_guide(self,annular_fail=True,as_xml=True):
        """
        Master function to return a guide star once the object is instantiated.
        
        Parameters
        ----------
        annular_fail : bool, optional
             If there is no guidestar in GC FOV, search an annulus and define
             the PA required to get a guidestar. Return most centralised
             candidate.
        as_xml : bool, optional
             Returns the result as an XML <target> element that can be added to
             a <field> element.
        
        Returns
        -------
        guide : astropy.Table
             Row from the Guide star catalogue.
        """


        self._set_geometry()
        self._retrieve_guidecats()
        guides = self._select_target()
        if (type(guides) == type(None)) and (annular_fail == True):
            print('No guide(s) found at fixed position - performing annular search')
            guides = self._annular_search()

        if type(guides) == type(None):
            print('No guide star(s) found...')
            return None
            
        return guides

