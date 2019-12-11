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


class GuideStar:
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
    mode : str
        Either LIFU or mIFU.
    nside : int, optional
        Override the default HEALPix nside value. Will likely end in tears.
    """

    import xml.dom.minidom
    from xml.dom.minidom import Node

    def __init__(self,ra,dec,pa,mode,nside=32):
        self.xml_template_url = 'http://casu.ast.cam.ac.uk/~dmurphy/opr3/swg/resources/BlankXMLTemplate.xml'
        self.xml_template = 'BlankXMLTemplate.xml'
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.nside = nside
        self.guides = None
        self.plate_scale = 17.8   ##  "/mm
        self.guide_url = "http://casu.ast.cam.ac.uk/~dmurphy/opr3/swg/resources/guides/"
        self.guide_filename = "Guides_S4_<HEALPIX>.fits"
        self.cache_dir = './cache/'

        if not mode in ['LIFU','mIFU']:
            print 'ERROR: must specify LIFU or mIFU mode'
            raise SystemExit(0)
        self.lifu = (mode == 'LIFU')

        if (not self.lifu) and (self.pa != 0):
            print 'ERROR: mIFU mode requires PA=0'
            raise SystemExit(0)

        if not os.path.isdir(self.cache_dir):
            cmd = 'mkdir -p %s'%(self.cache_dir)
            os.system(cmd)
        
        if not os.path.isfile(self.xml_template):
            self.wget(self.xml_template_url)
        
        
    def get_guide(self,annular_fail=True,as_xml=True,print_xml=False):
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
        print_xml : bool, optional
             Prints the XML results if as_xml=True.
        
        Returns
        -------
        guide : astropy.Table
             row from the Guide star catalogue.
        guide (if as_xml=True) : xml.dom.minidom.Element
             XML <target> element that can be inserted into a field XML.
        """


        self.set_geometry()
        self.retrieve_guidecats()
        guides = self.select_target()
        if (type(guides) == type(None)) and (annular_fail == True):
            print 'No guide(s) found at fixed position - performing annular search'
            guides = self.annular_search()

        if type(guides) == type(None):
            print 'No guide star(s) found...'
            return None
            
        if as_xml:
            if self.lifu:
                return self.to_xml(guides)
            else:
                xmls = [self.to_xml(guide) for guide in guides]
                if print_xml:
                    for x in xmls:
                        print x.toxml()

                return xmls
        else:
            return guides
        

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
        self.retrieve_guidecats()
        target = self.select_target(annular=True)
        return target
        
    def set_geometry(self,healpix=True):
        from astropy_healpix import HEALPix
        from astropy.coordinates import SkyCoord
        import astropy.coordinates as cc
        from numpy import unique

        if self.lifu:
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

    def retrieve_guidecats(self,clobber=False):

        from astropy.table import Table
        import astropy.io.fits as pyfits
        if len(self.healpix_indices) == 0:
            self.set_geometry()
            
        self.guide_files = []

        for hp in self.healpix_indices:
            fn = self.guide_filename.replace('<HEALPIX>',str(hp))
            url = self.guide_url+fn
            fn = self.cache_dir+fn
            if os.path.isfile(fn):
                if clobber == False:
                    print 'Using existing file %s'%(fn)
                    self.guide_files.append(fn)
                    continue
            print url
            self.wget(url,outname=fn)
            self.guide_files.append(fn)
            
        #self.guide_files = ['/scratch/Guides_S4_cname.fits']
        tabs = []
        for cat in self.guide_files:
            fits = pyfits.open(cat)
            tabs.append(Table(fits[1].data))
            
        if len(tabs) == 1:
            self.guides = tabs[0]
            
        else:
            t0 = tabs[0]
            for t in tabs[1:]:
                for g in t:
                    t0.add_row(g)
            self.guides = t0
                
        return

    def select_target(self,annular=False):
        import numpy
        from operator import indexOf
        from astropy.table import Table
        self.ra_g = self.guides['GAIA_RA']
        self.dec_g = self.guides['GAIA_DEC']
        
        if (not self.lifu):
            annular = True
            #fig = plt.figure(); sp = plt.subplot(aspect='equal'); plt.plot(self.ra_g,self.dec_g,'ko'); cir = plt.Circle((gs.ra,gs.dec),radius=1.0,lw=2.0,color='red'); sp.add_patch(cir); plt.show()

            
        if annular ==False:        
            filter1 = (self.ra_g > self.ra_min) & (self.ra_g < self.ra_max)
            filter2 = (self.dec_g > self.dec_min) & (self.dec_g < self.dec_max)
            filter = filter1 & filter2

            if (True in filter) == False:
                print 'No guide stars within the GC FOV!!'
                return None

            self.guides_filter = self.guides[where(filter)[0]]
            self.ra_g = self.guides_filter['GAIA_RA']
            self.dec_g = self.guides_filter['GAIA_DEC']
        
            self.dist = [((abs(r-self.ra_gc0)**2)+(abs(d-self.dec_gc0)**2))**0.5 for r,d in zip(self.ra_g,self.dec_g)]
        else:
            if self.lifu:
                #in annulus, want closest to the central radius of the annulus
                r_min = self.ra_min-self.ra
                r_max = self.ra_max-self.ra
            else:
                r_min = 0.95
                r_max = 1.0

                
            self.radii = numpy.array([(((_ra-self.ra)**2)+((_dec-self.dec)**2))**0.5 for _ra,_dec in zip(self.ra_g,self.dec_g)])
            filter = (self.radii > r_min) & (self.radii < r_max)
            self.guides_filter = self.guides[filter]
            self.ra_g = self.guides_filter['GAIA_RA']
            self.dec_g = self.guides_filter['GAIA_DEC']

            radii = self.radii[filter]
            if self.lifu:
                r0 = r_min + (0.5*(r_max-r_min))
                self.dist = [abs(d-r0) for d in radii]
            else:
                self.dist = radii
                
        minval = min(self.dist)
        g_index = indexOf(self.dist,minval)

        guide_sel = self.guides_filter[g_index]
        if (annular == False):
            if (guide_sel['GAIA_RA'] > self.ra_min) and (guide_sel['GAIA_RA'] < self.ra_max) and (guide_sel['GAIA_DEC'] > self.dec_min) and (guide_sel['GAIA_DEC'] < self.dec_max):
                print 'Guidestar candidate is %f arcmin from GC centre'%(minval*60.0)
                print self.guides[g_index]
                self.guide_sel = self.guides[g_index]
            else:
                print 'Closest candidate still lies outside of the GC FOV!!'
                self.guide_sel = self.guides[g_index]
                return None

        else:
            #do a quick report
            from copy import deepcopy
            dist_sort = deepcopy(self.dist)
            dist_sort.sort()
            if self.lifu:
                print 'Annular search summary (selected closest to centre of a rotated guidecam):'
            else:
                print 'Annular search summary:'
            print "#\t Dist (')  CNAME\t\t RA\t    Dec\t      angle\t Gaia_G mag"
            i = 0
            self.guides_filter['dist'] = self.dist

            for d in dist_sort:
                i = i + 1
                sel = ''
                if (i == 1) and (self.lifu):
                    sel = ' <-----'
                index = indexOf(self.dist,d)
                guide_candidate = self.guides_filter[index]
                ra_trans = self.ra - guide_candidate['GAIA_RA']
                dec_trans = self.dec - guide_candidate['GAIA_DEC']
                ang = numpy.arctan2(dec_trans,ra_trans)
                ang = (ang*180.0) / numpy.pi
                if ang < 0:
                    ang += 360.0
                print '#%d\t %1.2f\t   %s\t %1.4f   %1.4f   %1.3f\t %1.3f%s'%(i,d*60.0,guide_candidate['CNAME'],guide_candidate['GAIA_RA'],guide_candidate['GAIA_DEC'],ang,guide_candidate['GAIA_MAG_GG'],sel)

            self.guides_filter.sort('dist')

        if self.lifu:
            return guide_sel
        return self.guides_filter
    
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


    
    def to_xml(self,guide):
        self.new_xml()
        xml_target = self.targets_base[0].cloneNode(True)
        guide_ra = guide['GAIA_RA']
        guide_dec = guide['GAIA_DEC']

        dx = (self.ra - guide_ra)*self.plate_scale
        dy = (self.dec - guide_dec)*self.plate_scale
        xml_target.setAttribute('targx',str(dx))
        xml_target.setAttribute('targy',str(dy))

#        print 'WARNING - overriding targx, targy for now!'
        #manual override for the moment, position of targx,y
#        guide_targ.setAttribute('targx',"-110.0")
#        guide_targ.setAttribute('targy',"-500.55")

        #xml_target.setAttribute('fibreid',"9999")
        #xml_target.setAttribute('configid',"9999")
        xml_target.setAttribute('fibreid',"")
        xml_target.setAttribute('configid',"")


        xml_target.setAttribute('cname',str(guide['CNAME']))
        xml_target.setAttribute('targid',str(guide['TARGID']))
        xml_target.setAttribute('targra',str(guide['GAIA_RA']))
        xml_target.setAttribute('targdec',str(guide['GAIA_DEC']))
        xml_target.setAttribute('targpmra',str(guide['GAIA_PMRA']))
        xml_target.setAttribute('targpmdec',str(guide['GAIA_PMDEC']))
        xml_target.setAttribute('targprio',str(guide['TARGPRIO']))
        xml_target.setAttribute('targuse',str(guide['TARGUSE']))
        xml_target.setAttribute('targsrvy',str(guide['TARGSRVY']))
        xml_target.setAttribute('targname',str(guide['TARGNAME']))
        xml_target.setAttribute('targprog',str(guide['TARGPROG']))
        xml_target.setAttribute('targclass',str(guide['TARGCLASS']))
        xml_target.setAttribute('targcat',str(guide['TARGCAT']))
        xml_target.setAttribute('targepoch',str(guide['GAIA_EPOCH']))
        xml_target.setAttribute('targparal',str(guide['GAIA_PARAL']))
        xml_target.setAttribute('targprio',"10")
#        xml_target.setAttribute('ifu_spaxel',"")
        #xml_photom = self.targets_base[0].getElementsByTagName('photometry')[0]
        xml_photom = xml_target.getElementsByTagName('photometry')[0]
        bands = ['g','r','i']
        for b in bands:
            xml_photom.setAttribute('mag_%s'%(b),"")
            xml_photom.setAttribute('emag_%s'%(b),"")
        xml_photom.setAttribute('mag_gg',str(guide['GAIA_MAG_GG']))
        xml_photom.setAttribute('emag_gg',str(guide['GAIA_EMAG_GG']))
        xml_photom.setAttribute('mag_bp',str(guide['GAIA_MAG_BP']))
        xml_photom.setAttribute('emag_bp',str(guide['GAIA_EMAG_BP']))
        xml_photom.setAttribute('mag_rp',str(guide['GAIA_MAG_RP']))
        xml_photom.setAttribute('emag_rp',str(guide['GAIA_EMAG_RP']))

#        xml_target.appendChild(xml_photom)

        #now do photometry
        return xml_target





if __name__ =='__main__':


    if 0:
        import ifu
        ra = 178.835488822
        dec = 58.2835493041
        pa = 0.0
        gs = ifu.guidestar(ra,dec,pa)
        gs.set_geometry()
        guide = gs.get_guide(annular_fail=True,as_xml=True)

        #gs.retrieve_guidecats()
        #gs.select_target()

    if 0:
        gs = GuideStar(316.369609537,-4.71060356792,0,'mIFU')
        guides = gs.get_guide()
        for g in guides[:8]:
            print g.toxml()

    if 1:
        gs = GuideStar(316.369609537,-4.71060356792,0,'LIFU')
        guide = gs.get_guide()
        print guide.toxml()

            
    # gs.set_geometry()
    # gs.retrieve_guidecats()
    # gs.select_target(annular=True)
    
    print 'Fin'





