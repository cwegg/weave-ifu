#!/usr/bin/env python

#
# Copyright (C) 2018 Cambridge Astronomical Survey Unit
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


import xml.dom.minidom
import os
from astropy.io import fits


class xml_data:
    def __init__(self):
        self.xml_template_url = 'http://casu.ast.cam.ac.uk/~dmurphy/opr3/swg/resources/BlankXMLTemplate.xml'
        self.xml_template = 'BlankXMLTemplate.xml'
        if not os.path.isfile(self.xml_template):
            self.wget(self.xml_template_url)
        return

    def new_xml(self,root_data=None):
        #init the new XML
        try:
            dom = xml.dom.minidom.parse(self.xml_template)
        except xml.parsers.expat.ExpatError:
            print("File {0} would not parse".format(self.xml_template))
            raise SystemExit(0)
        self.ingest_xml(dom)
        if root_data:
            for key in root_data.keys():
                self.root.setAttribute(key, value=str(root_data[key]))


    def wget(self,url,outname=None):

        ##
        ## Put this into a generic tools module
        ##
        
        import os
        import time
        cmd = 'wget -q -t 1 -T 5 %s'%(url)
        if outname != None:
            print 'Downloading URL %s to %s'%(url,outname)
            cmd += ' -O %s'%(outname)
        os.system(cmd)
        return

    def ingest_xml(self,dom):
        self.dom = dom
        self.root = dom.childNodes[0]
        self.exposures = self.root.getElementsByTagName('exposures')[0]
        self.observation = self.root.getElementsByTagName('observation')[0]
        self.configure = dom.getElementsByTagName('configure')[0]
        self.surveys = dom.getElementsByTagName('surveys')[0]
        self.obsconstraints = dom.getElementsByTagName('obsconstraints')[0]
        self.dithering = dom.getElementsByTagName('dithering')[0]
        self.fields = dom.getElementsByTagName('fields')[0]
        self.base_target = self.fields.getElementsByTagName('target')[0]
        self.offset = self.observation.getElementsByTagName('offsets')[0]
        self.targets_base = self.fields.getElementsByTagName('target')

    def write_xml(self,filename):
        newxml = self.remove_empty_lines(self.dom.toprettyxml())
        finalxml = self.remove_xml_declaration(newxml)
        print 'Writing to %s'%(filename)
        with open(filename, 'w') as f:
            f.write(finalxml)    

        

    def remove_xml_declaration(self,xml_text):
        doc = xml.dom.minidom.parseString(xml_text)
        root = doc.documentElement
        xml_text_without_declaration = root.toxml(doc.encoding)
        return xml_text_without_declaration

    def remove_empty_lines(self,xml_text):
        reparsed = xml.dom.minidom.parseString(xml_text)
        return '\n'.join([line for line in reparsed.toprettyxml(indent=' ').split('\n') if line.strip()])

        



class ifu:
    """
    Convert the target-level data from a FITS input catalogue to a set of XMLs.
    
    This class provides code to take an input FITS catalogue containing IFU
    spaxel data and convert it into XML targets that are written into either a
    supplied XML or the BlankTemplate.xml file.
    
    Parameters
    ----------
    input_fits : str
        Filename of the input FITS catalogue.
    res : str, optional
        The resolution you wish to select from this input catalogue.
    output : str, optional
        Override the output filename, or leave as 'auto' to auto-generate names
        based on the pointing ID.
    version : str, optional
        Special actions are performed if this is set to 'OpR3b', including
        reading a 2nd extension containing <simulation> data.
    binning : str, optional
        The binning.
    """

    def __init__(self,input_fits,trimester,res='LR',output='auto',version='OpR3b',binning='1'):
        self.version = 0.6
        self.plate_scale = 17.8   ##  "/mm
        print 'Reading in %s'%(input_fits)
        self.input_fits = input_fits
        input_data = fits.open(input_fits)
        self.data = input_data[1].data
        self.author = input_data[0].header['AUTHOR']
        self.cc_report = input_data[0].header['CCREPORT']
        self.report_verbosity = input_data[0].header['VERBOSE']
        self.targsrvy = input_data[0].header['TARGSRVY']
        self.root_data = {}
        self.root_data['author'] = self.author
        self.root_data['cc_report'] = self.cc_report
        self.root_data['report_verbosity'] = str(self.report_verbosity)
        self.trimester = trimester
        
        self.binning = binning
        self.res = res
        self.xml_out = output
        self.xml_template = 'BlankXMLTemplate.xml'

        # ifu_lookup = './LIFUfibreTable.dat'
        # file = open(ifu_lookup,'r')
        # ifu_tab = file.readlines()
        # file.close()
        # self.spax_ids = [i.split()[2] for i in ifu_tab]
        # self.spax_lookup = {}
        # for line in ifu_tab:
        #     x = float(line.split()[0])
        #     y = float(line.split()[1])
        #     id = line.split()[2]
        #     fibreid = line.split()[4]
        #     self.spax_lookup[id] = {'x':x,'y':y,'fibreid':fibreid}


        # if len(self.spax_ids) != len(unique(self.spax_ids)):
        #     raise SystemExit('List of spaxel ids not unique')

        # self.cspax_id = 'C14'
        return

    def row_validator(self,row):
        assert row['IFU_DITHER'] in [-1,3,5]
        
        
        return False
    
    def generate_xmls(self):
       
        centrals = []
        #NORBI.X
        data_filter = []
        for d in self.data:
            #filter for IFU
            if int(d['PROGTEMP'][0]) > 3:
                data_filter.append(d)

        if len(data_filter) == 0:
            print 'The supplied catalogue does not provide IFU target data'

        lifu = []
        mifu = []
        for d in data_filter:

            # validation checks

            if d['PROGTEMP'][0] in ['4','5','6']:
                lifu.append(d)
            else:
                mifu.append(d)                


        # import pdb
        # pdb.set_trace()

        
        #how do you group entries belonging to the same OB (for custom dithers)?
        group_id = 'TARGID:TARGNAME:PROGTEMP:OBSTEMP'
        
        custom_dithers = {}
        for lifu_entry in lifu:
            if lifu_entry['IFU_DITHER'] != -1:
                #process this right away
                self.process_rows([lifu_entry])
            else:
                key = ':'.join([lifu_entry[subkey] for subkey in group_id.split(':')])
                try:
                    custom_dithers[key].append(lifu_entry)
                except KeyError:
                    custom_dithers[key] = [lifu_entry]

        if len(custom_dithers.keys()) > 0:
            #what happens if there are (eg) 10 pointings in a given key?
            #something like:
            # assert len(custom_dithers[key]) == custom_dithers[key]['PROGTEMP'][2]
            
            print 'Processing %d custom dither pointing groupings for LIFU'%(len(custom_dithers.keys()))
            for key in custom_dithers.keys():
                lifu_entry = custom_dithers[key]
                self.process_rows(lifu_entry)
                    
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

    def process_rows(self,rows):
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
        
        mode_lookup = {}
        mode_lookup['4'] = 'LIFU'
        mode_lookup['5'] = 'LIFU'
        mode_lookup['6'] = 'LIFU'
        mode_lookup['7'] = 'mIFU'
        mode_lookup['8'] = 'mIFU'
        mode_lookup['9'] = 'mIFU'



        ndither = len(rows)
        if len(rows) == 1:
            ndither = rows[0]['IFU_DITHER']

        this_xml = xml_data()
        this_xml.new_xml(root_data=self.root_data)

        #get a clone of the exposures element and wipe it of everything bar the initial calibs
        #make no assumptions about the OB calibration strategy, just take from the template
        exposures = this_xml.exposures.cloneNode(True)
        calibs_after_science = []
        pre_sci = True
        order = 0

        #remove the "ORDER undefined" comment nodes (as we will define them here)
        del_next = False
        ii = 0
        nodes = []
        exposures_filter = exposures.cloneNode(True)
        while len(exposures_filter.childNodes) != 0:
            for node in exposures_filter.childNodes:
                exposures_filter.removeChild(node)
        
        include_node = []
        for node in exposures.childNodes:
            if del_next:
                del_next = False
                continue
            try:
                if 'ORDER undefined' in str(node.data):
                    #remove carraige return text node as well!
                    del_next = True
                    continue
                else:
                    include_node.append(node)
            except AttributeError:
                ii += 1
                #add the node to the exposures_filter
                include_node.append(node)
                continue
            ii += 1

        for node in include_node:
            exposures_filter.appendChild(node)

        exposures = exposures_filter

        #identify the first comment node after the set of <exposure> elements
        #this allows us to insert the new elements there, rather than at the end
        
        comment_ref_node = None
        cr_node = exposures.childNodes[-1]
        pre_exp = True
        previous_node = exposures.childNodes[0]
        ii = 0
        for node in exposures.childNodes:
            if str(node.nodeName)[0] != '#':
                pre_exp = False

            if str(node.nodeName)[0] == '#':
                if pre_exp == False:
                    #must be a comment node
                    #previous must be a text node (<DOM Text node "u'\n      \n\n '...">)
                    if (str(previous_node.nodeName) == '#text') and (str(node.nodeName) == '#comment'):
                        comment_ref_node = node
                        #want to insert new exposures before this node
                        break
            previous_node = node
            ii += 1
                
        for exposure in exposures.getElementsByTagName('exposure'):
            if exposure.getAttribute('type') == 'science':
                pre_sci = False
                sibling = exposure.nextSibling
                exposures.removeChild(exposure)
                exposures.removeChild(sibling)
            else:
                sibling = exposure.nextSibling
                exposures.removeChild(sibling)                
                if pre_sci == False:
                    calibs_after_science.append(exposure)
                    exposures.removeChild(exposure)
                else:
                    order = int(exposure.getAttribute('order'))
        
        sci_dummy = this_xml.exposures.getElementsByTagName('exposure')[4]
        sci_exps = []
        order += 1
        first_sci_order = order
        for i in xrange(ndither):
            sci = sci_dummy.cloneNode(True)
            sci.setAttribute('order',value=str(order))
            sci.setAttribute('arm',value='both')
            sci_exps.append(sci)
            order += 1
            #exposures.appendChild(sci)
            exposures.insertBefore(sci,comment_ref_node)
            cr_node_clone = cr_node.cloneNode(True)
#            exposures.insertBefore(cr_node_clone,comment_ref_node)

        #now add the trailing calibration exposure elements
        for calib in calibs_after_science:
            calib.setAttribute('order',value=str(order))
            #exposures.appendChild(calib)
            exposures.insertBefore(calib,comment_ref_node)
            cr_node_clone = cr_node.cloneNode(True)
#            exposures.insertBefore(cr_node_clone,comment_ref_node)
            order += 1

        cr_node_clone = cr_node.cloneNode(True)
        exposures.insertBefore(cr_node_clone,comment_ref_node)


        #now the <observation> element:
        this_xml.observation.setAttribute('name',str(rows[0]['TARGID']))
        this_xml.observation.setAttribute('progtemp',str(rows[0]['PROGTEMP']))
        this_xml.observation.setAttribute('obs_type',str(mode_lookup[rows[0]['PROGTEMP'][0]]))
        this_xml.observation.setAttribute('trimester',str(self.trimester))
        this_xml.observation.setAttribute('pa',str(rows[0]['IFU_PA_REQUEST']))

        #the <configure> amd <survey> elements:
        survey = this_xml.surveys.getElementsByTagName('survey')[0]
        survey.setAttribute('name',value=str(self.targsrvy))
        survey.setAttribute('priority',value='1.0')

        if str(mode_lookup[rows[0]['PROGTEMP'][0]]) == 'LIFU':
            mode = 'LIFU'
            this_xml.configure.setAttribute('plate','LIFU')
            survey.setAttribute('max_fibres',value='603')
        else:
            #mIFU
            mode = 'mIFU'
            this_xml.configure.setAttribute('plate','PLATE_B')
            survey.setAttribute('max_fibres',value='740')
            
        this_xml.obsconstraints.setAttribute('obstemp',str(rows[0]['OBSTEMP']))

        #now generate field template
        fields_clone = this_xml.fields.cloneNode(True)

        #clear it
        field_list = fields_clone.getElementsByTagName('field')
        field_clone = field_list[0].cloneNode(True)
        for field in field_list:
            fields_clone.removeChild(field)

        target_list = field_clone.getElementsByTagName('target')
        node_list = field_clone.childNodes

        while len(field_clone.childNodes) != 0:
            for target in field_clone.childNodes:
                field_clone.removeChild(target)

        #now remove the existing <fields> element
        this_xml.observation.removeChild(this_xml.fields)

        this_xml.dithering.setAttribute('apply_dither',str(rows[0]['IFU_DITHER']))
        if len(rows) == 1:

            #this is a fixed dither-pattern - only provide the initial pointing

            #assembly order:
            #1.target
            #2.field
            #3.fields

            target = this_xml.base_target.cloneNode(True)
            print 'What about fibreid ifu_spaxel??' 
            target.setAttribute('targra',value=str(rows[0]['GAIA_RA']))
            target.setAttribute('targdec',value=str(rows[0]['GAIA_DEC']))
            target.setAttribute('targpmra',value=str(rows[0]['GAIA_PMRA']))
            target.setAttribute('targpmdec',value=str(rows[0]['GAIA_PMDEC']))
            target.setAttribute('targepoch',value=str(rows[0]['GAIA_EPOCH']))
            target.setAttribute('targparal',value=str(rows[0]['GAIA_PARAL']))
            target.setAttribute('targid',value=str(rows[0]['TARGID']))
            target.setAttribute('targname',value=str(rows[0]['TARGNAME']))
            target.setAttribute('targprog',value=str(rows[0]['TARGPROG']))
            target.setAttribute('targcat',value=self.input_fits.split('/')[-1])
            target.setAttribute('targprio',value=str(rows[0]['TARGPRIO']))
            target.setAttribute('targuse',value='T')

            #now add field data, and add the target to the field
            field_clone.setAttribute('order',value=str(first_sci_order))
            field_clone.setAttribute('RA_d',value=str(rows[0]['GAIA_RA']))
            field_clone.setAttribute('Dec_d',value=str(rows[0]['GAIA_DEC']))
            field_clone.appendChild(target)

            #now add this field to the fields element
            fields_clone.appendChild(field_clone)

            #...and finally to observation
            this_xml.observation.appendChild(fields_clone)


            
        index = 1
        output_name = '%s%s%d.xml'%(self.xml_out,mode.lower(),index)
        while os.path.isfile(output_name):
            index += 1
            output_name = '%s%s%d.xml'%(self.xml_out,mode.lower(),index)

        this_xml.write_xml(output_name)

        a = 1

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




if __name__ == '__main__':

    import os
    import glob

    input_fits = './input/WC_IFU.fits'
    trimester = '2016B2'
    
    try:
        assert os.path.isfile(input_fits)
    except:
        print 'Please supply valid FITS IFU catalogue driver file'
        print 'Usage: ./stage2.py filename'
        raise SystemExit(0)

    output_dir = '../stage3/input_tmp/'

    if 1:
        cmd = 'rm -Rf ../stage3/input_tmp/*.xml'
        os.system(cmd)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    stage2_ifu = ifu(input_fits,trimester,output=output_dir)   ## add specifiers for mIFU grouping / overloading fields /etc behaviour here
    stage2_ifu.generate_xmls()
    
    print('IFU XMLs written to: {0}'.format(output_dir))
