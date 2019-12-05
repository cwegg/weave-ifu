#!/usr/bin/env python
from astropy.io import fits as pyfits
import os
from numpy import unique, arange, random, where
from copy import deepcopy
import xml.dom.minidom
from xml.dom.minidom import Node
from math import radians,cos
import numpy
from guidestar import GuideStar

class lifu:
    """
 *+
 *  Name:
 *      lifu
 *
 *  Purpose:
 *      Convert the target-level data from a FITS input catalogue to a set of XMLs
 *
 *  Description:
 *      This class provides code to take an input FITS catalogue containing IFU spaxel
 *      data and convert it into XML targets that are written into either a supplied XML
 *      or the BlankTemplate.xml file
 *
 *  Arguments:
 *      input_fits : str
 *          Filename of the input FITS catlalogue
 *      res : str, optional
 *          The resolution you wish to select from this input catalogue
 *      output : str, optional
 *          Override the output filename, or leave as 'auto' to auto-generate 
 *          names based on the pointing ID
 *      version : str, optional
 *          Special actions are performed if this is set to 'OpR3b', including 
 *          reading a 2nd extension containing <simulation> data
 *
 *  Returned values:
 *      
 *  Notes:
 *      v0.1 - Still work in progress. This has been integrated into a Jupyter Notebook
 *             so that targets can be generated for an XML file that already has the 
 *             required elements populated (apart from the <targets>)
 *             
 *  Dependencies:
 *      Python core, numpy, astropy, astropy_healpix, GuideStar
 *
 *  Authors:
 *      David Murphy, Cambridge Astronomical Survey Unit (CASU, IoA)
 *                    dmurphy@ast.cam.ac.uk
 *
 *  Copyright:
 *      Copyright (C) 2018 Cambridge Astronomy Survey Unit.
 *
 *      This program is free software: you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation, either version 3 of the License, or
 *      (at your option) any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      You should have received a copy of the GNU General Public License
 *      along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
+*  
    """

    def __init__(self,input_fits,res='LR',output='auto',version='OpR3b',binning='1'):
        self.version = 0.6
        self.plate_scale = 17.8   ##  "/mm
        print 'Reading in %s'%(input_fits)
        self.input_fits = input_fits
        fits = pyfits.open(input_fits)
        self.data = fits[1].data
        if version == 'OpR3b':
            try:
                self.sim_data = fits[2].data
            except:
                print 'No 2nd extension supplied in FITS catalogue - cannot continue'
                raise SystemExit(0)
            
            self.sim_lookup = {}
            for i in xrange(len(self.data)):
                key = '%s %s %s %s'%(self.data[i]['TARGID'],self.data[i]['GAIA_RA'],self.data[i]['GAIA_DEC'],self.data[i]['IFU_SPAXEL'])
                bail = False
                try:
                    aa = self.sim_lookup[self.data[i]]
                    bail = True
                except KeyError:
                    self.sim_lookup[key] = self.sim_data[i]
                if bail:
                    print 'There was a duplicate entry for:'
                    print key
                    raise SystemExit(0)


        else:
            self.sim_data = None
            self.sim_lookup = {}
        self.binning = binning
        self.res = res
        self.xml_out = output
        self.xml_template = 'BlankXMLTemplate.xml'
        ifu_lookup = './LIFUfibreTable.dat'


        file = open(ifu_lookup,'r')
        ifu_tab = file.readlines()
        file.close()
        self.spax_ids = [i.split()[2] for i in ifu_tab]
        self.spax_lookup = {}
        for line in ifu_tab:
            x = float(line.split()[0])
            y = float(line.split()[1])
            id = line.split()[2]
            fibreid = line.split()[4]
            self.spax_lookup[id] = {'x':x,'y':y,'fibreid':fibreid}


        if len(self.spax_ids) != len(unique(self.spax_ids)):
            raise SystemExit('List of spaxel ids not unique')

        self.cspax_id = 'C14'
        return
    

    def get_central_spaxels(self,dither_group=False,verbose=False):
        res_lookup = {}
        res_lookup['1'] = 'LR'
        res_lookup['2'] = 'HR'
        res_lookup['3'] = 'HR'
        res_lookup['4'] = 'LR'
        res_lookup['5'] = 'HR'
        res_lookup['6'] = 'HR'
        res_lookup['7'] = 'LR'
        res_lookup['8'] = 'HR'
        res_lookup['9'] = 'HR'

        centrals = []
        #NORBI.X
        data_filter = []
        for d in self.data:
            #filter for IFU
            if int(d['PROGTEMP'][0]) > 3:
                #filter for resolution
                if res_lookup[d['PROGTEMP'][0]] == self.res:
                    #filter for binning
                    if d['PROGTEMP'][4] == self.binning:
                        #ok, add!
                        data_filter.append(d)

        if len(data_filter) == 0:
            print 'The supplied catalogue does not provide target data for the following requested configuration:'
            print 'PROGTEMP (N) > 3 (L/mIFU)'
            print 'Resolution = %s'%(self.res)
            print 'Binning = %sx'%(self.binning)
            raise SystemExit(0)

        for d in data_filter:
            if (d['IFU_SPAXEL'] == self.cspax_id):
#            if (d['IFU_SPAXEL'] == self.cspax_id) and (res_lookup[d['PROGTEMP'][0]] == self.res) and (int(d['PROGTEMP']) > 3):
#            if (d['TARGID'] not in self.spax_ids) and (d['TARGPROG'] == self.res):

                centrals.append(d)
                if verbose:
                    print d


        if (len(centrals) > 0) and dither_group:
            if verbose:
                print ''
                print 'Groupings:'

            groups = {}
            for c in centrals:
                radec = '%1.10f %1.10f'%(c['GAIA_RA'],c['GAIA_DEC'])
                try:
                    groups[c['TARGID']].append(c)
                except KeyError:
                    groups[c['TARGID']] = [c]

            if len(groups.keys()) > 0:
                if verbose:
                    for key in groups.keys():
                        print 'Group %s'%(key)
                        for d in groups[key]:
                            print d
                    print '\n'
            else:
                if verbose:
                    print 'No dithers... sorry'

        if ((len(groups)) == 0) and len(centrals) == 0:
            raise SystemExit('No centrals of groups found')

        if len(groups) == 0:
            ##??
            groups = [centrals[0]]

        self.groups = groups
        self.centrals = centrals
        if dither_group:
            return centrals,groups
        return centrals


    def ingest_xml(self,dom):
        self.dom = dom
        self.root = dom.childNodes[0]
        self.programme = self.root.getElementsByTagName('programme')[0]
        self.observation = self.root.getElementsByTagName('observation')[0]
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


    def generate_targets(self,field,targ_base,cross_match=None):
        ra0 = float(field.attributes['RA_d'].value)
        dec0 = float(field.attributes['Dec_d'].value)
        targets = []

        if cross_match:
            user_lookup = {}
            #get the central spaxel
            cspax = None
            for spax in cross_match:
                if spax['IFU_SPAXEL'] == self.cspax_id:
                    cspax = spax
                    break
            if not cspax:
                raise SystemExit('Could not identify central spaxel')

            user_lookup[self.cspax_id] = cspax

            
            for spax in cross_match:
                if spax['IFU_SPAXEL'] in self.spax_ids:
                    #this is not the central spaxel - put in the lookup table, overwrite the targid
                    user_lookup[spax['IFU_SPAXEL']] = spax
                    #spax['TARGID'] = cspax['TARGID']    ##no need to do this now

            #what stuff should we definitely not overwrite? spatial data mainly
            exclusions = ['targx','targy','targra','targdec','targprog','fibreid','configid']
            exclusions = ['targx','targy','targra','targdec','fibreid','configid']

        keyword_translate = {}
        keyword_translate['targepoch'] = 'GAIA_EPOCH'
        keyword_translate['targpmra'] = 'GAIA_PMRA'
        keyword_translate['targpmdec'] = 'GAIA_PMDEC'


        mandatory = ['targsrvy','targcat','targid','targprio','targuse','targclass','progtemp','obstemp','','']
        
        #need to set a fibreID as well (as D.Terrett won't do this in configure for IFU)
        fibreID = 0
        non_compliant = False
        nc_messages = []
        for spaxID in self.spax_lookup.keys():
            fibreID = self.spax_lookup[spaxID]['fibreid']
            target = targ_base.cloneNode(True)
            targ_dec = dec0+((float(self.spax_lookup[spaxID]['y'])/3600))
            targ_ra = ra0-((float(self.spax_lookup[spaxID]['x'])/3600)/cos(radians(targ_dec)))


            target.setAttribute('targra', value=str(targ_ra)) 
            target.setAttribute('targdec', value=str(targ_dec)) 
            target.setAttribute('targx', value=str(float(self.spax_lookup[spaxID]['x'])/self.plate_scale))
            target.setAttribute('targy', value=str(float(self.spax_lookup[spaxID]['y'])/self.plate_scale))


            target.setAttribute('ifu_spaxel', value=spaxID) 

            target.setAttribute('fibreid', value=str(fibreID)) 
            target.setAttribute('configid', value=str(fibreID)) 

            if cross_match:
                try:
                    match = user_lookup[spaxID]
                except KeyError:
                    print 'Missing spaxel %s - cannot continue'%(spaxID)
                    raise SystemExit(0)
                for key in target.attributes.keys():
                    if key not in exclusions:
                        key_fits = key
                        if key in keyword_translate.keys():
                            key_fits = keyword_translate[key]
                            
                        value = match[key_fits]

                        if value == 'GALAX':
                            #grr..
                            value = 'GALAXY'
                    
                        
                        if (type(value) != type('')) and  numpy.isnan(value):
                            value = ""

                        if (type(value) != type('')) and  numpy.isnan(value) and (False):  ##don't eval this one...
                            non_compliant = True
                            err_string = 'NaN entry in the input FITS catalogue for %s - this catalogue is NOT compliant'%(key)
                            nc_messages.append(err_string)
                            value = "0.0"
                        
                        elif type(value) != type(''):
                            if numpy.isnan(match[key_fits]):
                                value = ''
                            else:
                                value = str(match[key_fits])
                        target.setAttribute(key, value=value) 


                targ_photom = target.getElementsByTagName('photometry')[0]


                #now for the photometry

                

                for key in targ_photom.attributes.keys():
                    try:
                        photom_data = match[key]
                        if numpy.isnan(photom_data):
                            photom_data = ''
                    except:
                        photom_data = ''
                    
                    if photom_data != '':
                        photom_data = '%1.3f'%(photom_data)
                    targ_photom.setAttribute(key, value=photom_data)
                if type(self.sim_data) != type(None):
                    #now add the <simulation tag>
                    match_key = '%s %s %s %s'%(match['TARGID'],match['GAIA_RA'],match['GAIA_DEC'],match['IFU_SPAXEL'])
                    sim_match = self.sim_lookup[match_key]
                    targ_sim = target.getElementsByTagName('simulation')[0]
                    for key in targ_sim.attributes.keys():
                        sim_data = sim_match[key.upper()]
                        
                        if (type(sim_data) != type('')) and (numpy.isnan(sim_data)):
                            sim_data = ''

                        if key == 'template':
                            if '/' in sim_data:
                                sim_data = sim_data.split('/')[-1]
#                        try:
#                            sim_data = sim_match[key.upper()]
#                            if numpy.isnan(sim_data):
#                                sim_data = ''
#                        except:
#                            sim_data = ''
                        targ_sim.setAttribute(key, value=str(sim_data))
                    a = 1                        
            #if this is a sky spaxel, remove the photometry subelements
            if match['TARGUSE'] == 'S':
                target.removeChild(targ_photom)
                if type(self.sim_data) != type(None):
                    targ_sim = target.getElementsByTagName('simulation')[0]
                    target.removeChild(targ_sim)

            targets.append(target)

        if len(nc_messages) > 0:
            nc_messages = unique(nc_messages)
            print 'WARNING: input FITS file target data was NON-COMPLIANT:'
            for nc in nc_messages:
                print '\t %s'%(nc)
            
        return targets

    def generate(self):
       print 'Generating XML OBs for %d IFU pointings'%(len(self.groups.keys()))
       self.root.setAttribute('comment', value='XMLIFU version %s'%(str(self.version)))
        
       #set configure elements for configure
       if 1:
           #self.configure.setAttribute('min_guide', value='1')
           #self.configure.setAttribute('min_sky', value='1')
           #self.configure.setAttribute('num_guide_fibres', value='1')
           self.configure.setAttribute('num_sky_fibres', value='603')
           self.configure.setAttribute('plate', value='LIFU')
           
       
       for groupID in self.groups.keys():
            #these need to be put into one OB, with distinct field entries
            group = self.groups[groupID]

            order = 2


            #create new observation instance
            _observation = self.observation.cloneNode(True)  #deepcopy(observation)
            #wipe any field entries

            _Fields = _observation.getElementsByTagName('fields')[0]
            _field_list = _Fields.getElementsByTagName('field')
            for _field in _field_list:
                _Fields.removeChild(_field)

            _observation.removeChild(_Fields)
                
            offset_ra = "0.0"
            offset_dec = "0.0"
            base_ra = float(group[0]['GAIA_RA'])
            base_dec = float(group[0]['GAIA_DEC'])
            fields = []
            children = []
            j = 0
            for dither in group:
                j += 1
                print 'Pointing %s dither %d'%(groupID,j)
                order += 1
                ra = float(dither['GAIA_RA'])
                dec = float(dither['GAIA_DEC'])
                pa = dither['IFU_PA']


                _field = self.field.cloneNode(True)#deepcopy(field)
                _targ_tmps = _field.getElementsByTagName('target')
                for _t in _targ_tmps:
                    _field.removeChild(_t)
                _field.setAttribute('RA_d', value=str(ra)) 
                _field.setAttribute('Dec_d', value=str(dec)) 
#                _field.setAttribute('pa', value=str(pa)) 
                if len(group) > 1:
                    _field.setAttribute('order', value=str(order))
                else:
                    _field.setAttribute('order', value="")
                fields.append(_field)

#                user_spaxel_data = [dither] + list(get_spaxel_children(dither,self.data))
                user_spaxel_data = list(get_spaxel_children(dither,self.data))
                
                targets = self.generate_targets(_field,self.targets_base[1],cross_match=user_spaxel_data)

                #now add a guide star
                print 'WARNING - using PA=0 for now!'
                gs = GuideStar(ra,dec,0,'LIFU',nside=32)
                guide_targ = gs.get_guide(as_xml=True)
                self.guide_targ = guide_targ

                
                _field.appendChild(guide_targ)
                
                for t in targets:
                    _field.appendChild(t)

                _Fields.appendChild(_field)

                if dither == group[0]:
                    #skip the rest - this is the first dither
                    continue

                delta_ra = (ra-base_ra)*cos(radians(dec))

                delta_ra = ra-base_ra
                delta_dec = dec-base_dec

                offset_ra = offset_ra + ' %s'%(delta_ra*3600)
                offset_dec = offset_dec + ' %s'%(delta_dec*3600)

            _observation.appendChild(_Fields)
            
            _observation.setAttribute('name', value=dither['TARGID'])
            _observation.setAttribute('pa', value=str(dither['IFU_PA']))

            #create the offsets
            _offsets = _observation.getElementsByTagName('offsets')[0]
            _offsets.setAttribute('offset_step_ra', value=offset_ra) 
            _offsets.setAttribute('offset_step_dec', value=offset_dec) 
            self.root.removeChild(self.observation)
            self.root.appendChild(_observation)
            
            #set the "new" observation to self.observation (for the next pointing)
            self.observation = _observation


            newxml = remove_empty_lines(self.dom.toprettyxml())
            finalxml = remove_xml_declaration(newxml)
            import os
            full_path = os.getcwd()

            if self.xml_out == 'auto':
                xml_out = '%s.xml'%(dither['TARGID'])
                print 'Writing to %s/%s'%(full_path,xml_out)
                
            else:
                xml_out = self.xml_out

            with open(xml_out, 'w') as f:
                f.write(finalxml)    


       return



def remove_xml_declaration(xml_text):
    doc = xml.dom.minidom.parseString(xml_text)
    root = doc.documentElement
    xml_text_without_declaration = root.toxml(doc.encoding)
    return xml_text_without_declaration

def remove_empty_lines(xml_text):
    reparsed = xml.dom.minidom.parseString(xml_text)
    return '\n'.join([line for line in reparsed.toprettyxml(indent=' ').split('\n') if line.strip()])



def get_spaxel_children(cspax,data):
    radec_hash = '%s %s'%(cspax['GAIA_RA'],cspax['GAIA_DEC'])
    children = []
    for d in data:
#        if d['IFU_SPAXEL'] == 'C14':
#            if (_radec_hash == radec_hash):
#                a = 1
        _radec_hash = '%s %s'%(d['GAIA_RA'],d['GAIA_DEC'])
        if (_radec_hash == radec_hash):
#        if (_radec_hash == radec_hash) and (d != cspax):
            children.append(d)
    return children




if __name__ =='__main__':
    input_fits = './WA_FITS_catalogue_test.fits'
    input_fits = './WA_FITS_catalogue_180523-132206ex2.fits'
    input_fits = './WA_FITS_catalogue_180608-143105ex2.fits'
    input_fits = './scipLFG108.2-00.6.fits'
    input_fits = './WC_S4.fits'
    input_fits = '../test_data/WA_FITS_catalogue_180608-143105ex2.fits'
    import sys
    try:
        assert os.path.isfile(input_fits)
    except:
        print 'Please supply valid FITS catalogue name'
        print 'Usage: ./ifu.py filename'
        raise SystemExit(0)


    
    mylifu = lifu(input_fits,res='HR',binning='1')
    centrals,groups = mylifu.get_central_spaxels(dither_group=True)
    mylifu.new_xml()
    mylifu.generate()
    print 'Fin'





