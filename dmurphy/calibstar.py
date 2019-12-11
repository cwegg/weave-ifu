#!/usr/bin/env python3

#
# Copyright (C) 2018 Cambridge Astronomy Survey Unit
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
Notes:
    v0.1 - Still work in progress

Dependencies:
    Python core, numpy, astropy, astropy_healpix

Authors:
    David Murphy, Cambridge Astronomical Survey Unit (CASU, IoA)
                  dmurphy@ast.cam.ac.uk
"""


from astropy.io import fits as pyfits
import os
from numpy import unique, arange, random, where
from copy import deepcopy
import xml.dom.minidom
from xml.dom.minidom import Node
from math import radians,cos
import numpy


class CalibStar:
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
    mode : str
        Either MOS or mIFU.
    nside : int, optional
        Override the default HEALPix nside value. Will likely end in tears.
    annular : bool
        Search only within an annular radius at the edge of the WEAVE FOV.
    plot : bool
        Plot the results?
    """

    import xml.dom.minidom
    from xml.dom.minidom import Node

    def __init__(self,ra,dec,pa,mode,nside=32,annular=False,plot=False):
        self.xml_template_url = 'http://casu.ast.cam.ac.uk/~dmurphy/opr3/swg/resources/BlankXMLTemplate.xml'
        self.xml_template = 'BlankXMLTemplate.xml'
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.nside = nside
        self.calibs = None
        self.plate_scale = 17.8   ##  "/mm
        self.calib_url = "http://casu.ast.cam.ac.uk/~dmurphy/opr3/swg/resources/calibs/"
        self.calib_filename = "WD_S4_<HEALPIX>.fits"
        self.cache_dir = './cache/'
        self.annular = annular
        self.plot = plot
        
        if not mode in ['MOS','mIFU']:
            print 'ERROR: must specify MOS or mIFU mode'
            raise SystemExit(0)
        self.mos = (mode == 'MOS')
        
        if self.pa != 0:
            print 'ERROR: CalibStar requires PA=0'
            raise SystemExit(0)

        if not os.path.isdir(self.cache_dir):
            cmd = 'mkdir -p %s'%(self.cache_dir)
            os.system(cmd)
        
        if not os.path.isfile(self.xml_template):
            self.wget(self.xml_template_url)

        self.cspax_id = 'C03'
        if self.mos:
            self.cspax_id = ""

        
    def get_calib(self,annular_fail=True,as_xml=True,print_xml=False):
        """
        Master function to return a calib star once the object is instantiated.
        
        Parameters
        ----------
        annular_fail : bool, optional
             (for mIFU) If there are no calibrators in the annulus, search the
             whole field.
        as_xml : bool, optional
             Returns the result as an XML <target> element that can be added to
             a <field> element.
        print_xml : bool, optional
             Prints the XML results if as_xml=True.
        
        Returns
        -------
        calibs : astropy.Table
             row(s) from the WD calibration star catalogue.
        calibs (if as_xml=True) : xml.dom.minidom.Element
             XML <target> element that can be inserted into a field XML.
        """


        self.set_geometry()
        self.retrieve_calibcats()
        calibs = self.select_target()
        if (type(calibs) == type(None)) and (annular_fail == True) and (not self.mos):
            print 'No calibration stars found in annulus - performing full search'
            calibs = self.annular_search(annular=False)

        if type(calibs) == type(None):
            print 'No calibration stars found...'
            return None
            
        if as_xml:
            xmls = [self.to_xml(calib) for calib in calibs]
            if print_xml:
                for x in xmls:
                    print x.toxml()
            return xmls
        else:
            return calibs
        

    def wget(self,url,outname=None):
        import os
        import time
        cmd = 'wget -q -t 1 -T 5 %s'%(url)
        if outname != None:
            print 'Downloading URL %s to %s'%(url,outname)
            cmd += ' -O %s'%(outname)
        os.system(cmd)
        return

    def annular_search(self):
        from astropy_healpix import HEALPix
        from astropy.coordinates import SkyCoord
        import astropy.coordinates as cc
        hp = HEALPix(nside=self.nside, order='nested', frame=cc.ICRS())
        
        self.set_geometry(healpix=False)
        r_min = self.ra_min-self.ra
        r_max = self.ra_max-self.ra
        radius = self.ra_max-self.ra

        #get the healpix IDs covering an annulus centred on self.ra,dec

        in_annulus = []
        self.healpix_indices = []
        print 'Populating annulus and determining HEALPix coverage...'
        while(len(in_annulus)) < 500:
            rnd_ra = self.ra+(2*(random.random()-0.5)*radius)
            rnd_dec = self.dec+(2*(random.random()-0.5)*radius)
            rnd_dist = (((rnd_ra-self.ra)**2)+((rnd_dec-self.dec)**2))**0.5
            if rnd_dist > r_min:
                if rnd_dist < r_max:
                    self.healpix_indices.append(hp.skycoord_to_healpix(SkyCoord(rnd_ra,rnd_dec,unit='deg')))
                    in_annulus.append([rnd_ra,rnd_dec])
                #print len(in_annulus)
        print '....done'
        self.healpix_indices = unique(self.healpix_indices)
        print self.healpix_indices
        self.retrieve_calibcats()
        target = self.select_target(annular=True)
        return target
        
    def set_geometry(self,healpix=True):
        from astropy_healpix import HEALPix
        from astropy.coordinates import SkyCoord
        import astropy.coordinates as cc
        from numpy import unique


        self.ra_min = self.ra - 1.0
        self.ra_max = self.ra + 1.0
        self.dec_min = self.dec - 1.0
        self.dec_max = self.dec + 1.0

            
        if healpix:
            hp = HEALPix(nside=self.nside, order='nested', frame=cc.ICRS())
            self.healpix_indices = []
            ra = [self.ra_min,self.ra_min,self.ra_max,self.ra_max]
            dec = [self.dec_min,self.dec_max,self.dec_min,self.dec_max]

            dx = self.ra_max - self.ra_min
            dy = self.dec_max - self.dec_min

            for i in arange(500):
                ra.append(self.ra_min+(random.random()*dx))
                dec.append(self.dec_min+(random.random()*dy))

            for r,d in zip(ra,dec):
                self.healpix_indices.append(hp.skycoord_to_healpix(SkyCoord(r,d,unit='deg')))
            self.healpix_indices = unique(self.healpix_indices)

        return

    def retrieve_calibcats(self,clobber=False):
        from astropy.table import Table
        import astropy.io.fits as pyfits
        if len(self.healpix_indices) == 0:
            self.set_geometry()
            
        self.calib_files = []

        for hp in self.healpix_indices:
            fn = self.calib_filename.replace('<HEALPIX>',str(hp))
            url = self.calib_url+fn
            fn = self.cache_dir+fn
            if os.path.isfile(fn):
                if clobber == False:
                    print 'Using existing file %s'%(fn)
                    self.calib_files.append(fn)
                    continue
            print url
            self.wget(url,outname=fn)
            self.calib_files.append(fn)
            
        
        tabs = []
        for cat in self.calib_files:
            fits = pyfits.open(cat)
            tabs.append(Table(fits[1].data))
            
        if len(tabs) == 1:
            self.calibs = tabs[0]
            
        else:
            t0 = tabs[0]
            for t in tabs[1:]:
                for g in t:
                    t0.add_row(g)
            self.calibs = t0
                
        return

    def select_target(self):
        import numpy
        from operator import indexOf
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

                
        self.radii = numpy.array([(((_ra-self.ra)**2)+((_dec-self.dec)**2))**0.5 for _ra,_dec in zip(self.ra_c,self.dec_c)])
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
        from copy import deepcopy
        dist_sort = deepcopy(self.dist)
        dist_sort.sort()
        if self.annular:
            print 'Annular search summary:'
        else:
            print 'Search summary:'
        print "#\t Dist (')  CNAME\t\t RA\t    Dec\t      angle\t Gaia_G mag"
        i = 0
        self.calibs_filter['dist'] = self.dist

        for d in dist_sort:
            i = i + 1
            index = indexOf(self.dist,d)
            calib_candidate = self.calibs_filter[index]
            ra_trans = self.ra - calib_candidate['GAIA_RA']
            dec_trans = self.dec - calib_candidate['GAIA_DEC']
            ang = numpy.arctan2(dec_trans,ra_trans)
            ang = (ang*180.0) / numpy.pi
            if ang < 0:
                ang += 360.0
            print '#%d\t %1.2f\t   %s\t %1.4f   %1.4f   %1.3f\t %1.3f'%(i,d*60.0,calib_candidate['CNAME'],calib_candidate['GAIA_RA'],calib_candidate['GAIA_DEC'],ang,calib_candidate['GAIA_MAG_GG'])

        self.calibs_filter.sort('dist')
            
        if self.plot:
            import matplotlib.pyplot as plt
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
    
    def ingest_xml(self,dom):
        self.dom = dom
        self.root = dom.childNodes[0]
        self.programme = self.root.childNodes[3]
        self.observation = self.root.childNodes[5]
        self.configure = dom.getElementsByTagName('configure')[0]
        self.field = dom.getElementsByTagName('field')[0]
        self.base_target = self.field.getElementsByTagName('target')[0]
        self.offset = self.observation.getElementsByTagName('offsets')[0]
        self.targets_base = self.field.getElementsByTagName('target')

    def new_xml(self):
        #init the new XML
        try:
            dom = xml.dom.minidom.parse(self.xml_template)
        except xml.parsers.expat.ExpatError:
            print("File {0} would not parse".format(self.xml_template))
            raise SystemExit(0)
        self.ingest_xml(dom)


    
    def to_xml(self,calib):
        self.new_xml()
        xml_target = self.targets_base[2].cloneNode(True)
        calib_ra = calib['GAIA_RA']
        calib_dec = calib['GAIA_DEC']

        xml_target.setAttribute('targx','%%%')
        xml_target.setAttribute('targy','%%%')
        
        xml_target.setAttribute('fibreid','%%%')
        xml_target.setAttribute('configid',"%%%")


        ## NB: this will FAIL when we move from OpR3b catalogues, due to new column names for errors etc
        
        xml_target.setAttribute('cname',str(calib['CNAME']))
        xml_target.setAttribute('targid',str(calib['TARGID']))
        xml_target.setAttribute('targra',str(calib['GAIA_RA']))
        xml_target.setAttribute('targdec',str(calib['GAIA_DEC']))
        xml_target.setAttribute('targpmra',str(calib['GAIA_PMRA']))
        xml_target.setAttribute('targpmdec',str(calib['GAIA_PMDEC']))
        xml_target.setAttribute('targprio',str(calib['TARGPRIO']))
        xml_target.setAttribute('targuse',str(calib['TARGUSE']))
        xml_target.setAttribute('targsrvy',str(calib['TARGSRVY']))
        xml_target.setAttribute('targname',str(calib['TARGNAME']))
        xml_target.setAttribute('targprog',str(calib['TARGPROG']))
        xml_target.setAttribute('targclass',str(calib['TARGCLASS']))
        xml_target.setAttribute('targcat',str(calib['TARGCAT']))
        xml_target.setAttribute('targepoch',str(calib['GAIA_EPOCH']))
        xml_target.setAttribute('targparal',str(calib['GAIA_PARAL']))
        xml_target.setAttribute('targprio',"10")
        xml_target.setAttribute('ifu_spaxel',self.cspax_id)

        xml_photom = xml_target.getElementsByTagName('photometry')[0]
        bands = ['g','r','i']
        for b in bands:
            xml_photom.setAttribute('mag_%s'%(b),"")
            xml_photom.setAttribute('emag_%s'%(b),"")
        xml_photom.setAttribute('mag_gg',str(calib['GAIA_MAG_GG']))
        xml_photom.setAttribute('emag_gg',str(calib['GAIA_EMAG_GG']))
        xml_photom.setAttribute('mag_bp',str(calib['GAIA_MAG_BP']))
        xml_photom.setAttribute('emag_bp',str(calib['GAIA_EMAG_BP']))
        xml_photom.setAttribute('mag_rp',str(calib['GAIA_MAG_RP']))
        xml_photom.setAttribute('emag_rp',str(calib['GAIA_EMAG_RP']))

#        xml_target.appendChild(xml_photom)

        #now do photometry
        return xml_target




def warning(self, message, *args, **kws):
    if self.isEnabledFor(logging.WARNING):
        self._log(logging.WARNING, message, args, **kws)
        raise Exception(message)


if __name__ =='__main__':

#    import logging
#    logging.Logger.warning = warning

    import warnings

    warnings.simplefilter('error', FutureWarning)
    
    
    cs = CalibStar(316.369609537,-4.71060356792,0,'mIFU',annular=False,plot=True)
    calibs = cs.get_calib(print_xml=True)
    print 'Fin'





