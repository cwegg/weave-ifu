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
import numpy


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
        self.root_data = {}
        self.root_data['author'] = self.author
        self.root_data['cc_report'] = self.cc_report
        self.root_data['report_verbosity'] = str(self.report_verbosity)
        self.trimester = trimester
        
        self.binning = binning
        self.res = res
        self.xml_out = output
        self.xml_template = 'BlankXMLTemplate.xml'

        return

    def row_validator(self,row):
        print 'WARNING: row validator not implemented!'
        #assert row['IFU_DITHER'] in [-1,3,5]
        #check length of cc_report != max length (70? 69?) issue warning


        return True
    
    def generate_xmls(self,mifu_mode=1,mifu_ncalibs=2):
        assert mifu_mode in [1,2,3]
        #what happens if there are (eg) 100 bundles in a given key?
        #user needs to decide 3 options
        #1. Include all
        #2. Aufbau principle - fill up to N (18? ..leaving 2 for calibration bundles) then make a new XML
        #3. Equipartition - divide up the bundles equally. In the case of 100:
        #                   100/18. = 5.5 --> 6 fields
        #                   100/6 = 17 bundles per field (with remainder added to a final XML), ie:
        #                   5 x 17 bundles + 1 x 15 bundles


        
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
            self.row_validator(d)
            if d['PROGTEMP'][0] in ['4','5','6']:
                lifu.append(d)
            else:
                mifu.append(d)           

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
            print ''
            print 'Processing %d custom dither pointing groupings for LIFU'%(len(custom_dithers.keys()))
            for key in custom_dithers.keys():
                lifu_entry = custom_dithers[key]
                self.process_rows(lifu_entry)


        #how do you group bundles belonging to the same field?
        group_id = 'TARGNAME:PROGTEMP:OBSTEMP:IFU_DITHER'
        

        fields = {}
        for mifu_entry in mifu:
            key = ':'.join([str(mifu_entry[subkey]) for subkey in group_id.split(':')])
            try:
                fields[key].append(mifu_entry)
            except KeyError:
                fields[key] = [mifu_entry]

        if len(fields.keys()) > 0:
            print ''
            print 'Processing %d mIFU bundle groupings'%(len(fields.keys()))
            # for the moment, just go with (1), but implement (2) and (3)
            for key in fields.keys():
                mifu_entry = fields[key]
                if mifu_mode == 1:
                    if len(mifu_entry) > (20 - mifu_ncalibs):
                        print 'WARNING: Bundles in field %s (%d) exceeds maximum when including calibration bundles (%d). Change mifu_mode or remove excess bundles downstream'%(key,len(mifu_entry),mifu_ncalibs)
                        self.process_rows(mifu_entry)

                elif len(mifu_entry) > (20 - mifu_ncalibs):
                    print 'Group %s: filling XML files according to mifu_mode=%d'%(key,mifu_mode)
                    if mifu_mode == 2:
                        #2. Aufbau principle - fill up to N (18? ..leaving 2 for calibration bundles) then make a new XML
                        max_sci_bundles = 20 - mifu_ncalibs
                        print 'Will create XML files with %d science bundles inside'%(max_sci_bundles)
                        
                    if mifu_mode == 3:
                        max_sci_bundles = 20 - mifu_ncalibs
                        if (len(mifu_entry)/float(max_sci_bundles)).is_integer():
                            nxml = (len(mifu_entry) / (max_sci_bundles))
                        else:
                            nxml = ((len(mifu_entry) / (max_sci_bundles))) + 1

                        max_sci_bundles = float(len(mifu_entry)) / float(nxml)
                        if (max_sci_bundles).is_integer():
                            max_sci_bundles = int(max_sci_bundles)
                        else:
                            max_sci_bundles = int(max_sci_bundles) - 1
                        print 'Will create %d XML files each with %d science bundles inside'%(nxml,max_sci_bundles)

                    added_rows = []
                    rows = iter(mifu_entry)
                    while len(added_rows) != len(mifu_entry):
                        this_xml = []
                        while len(this_xml) != max_sci_bundles:
                            try:
                                _row = rows.next()
                                added_rows.append(_row)
                                this_xml.append(_row)
                            except StopIteration:
                                break
                        self.process_rows(this_xml)
                        
                        
                else:
                    self.process_rows(mifu_entry)
                    
                    
                        


            
        

        if 1:
            return


                
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

        if str(mode_lookup[rows[0]['PROGTEMP'][0]]) == 'LIFU':
            mode = 'LIFU'
            ndither = len(rows)
            if len(rows) == 1:
                ndither = rows[0]['IFU_DITHER']
            
        else:
            #mIFU
            mode = 'mIFU'
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
            if (calib.getAttribute('arm') in ['both','blue']):
                #assumes sequence is always red,blue for arm-independent <exposure> entries!
                order += 1

        cr_node_clone = cr_node.cloneNode(True)
        exposures.insertBefore(cr_node_clone,comment_ref_node)

        #now remove the old exposures element and replace with this one
        this_xml.exposures.parentNode.replaceChild(exposures,this_xml.exposures)
        this_xml.exposures = exposures
        #this_dom.insertBefore()
        

        
        #the <configure> amd <survey> elements:
        survey = this_xml.surveys.getElementsByTagName('survey')[0]
        #make a clone and remove the placeholder <survey>
        survey_clone = survey.cloneNode(True)
        this_xml.surveys.removeChild(survey)
        #get initial comment in <surveys>, to allow an insertBefore
        sruveys_comment = this_xml.surveys.childNodes[0]
        a = 1
        
        all_surveys = numpy.unique([r['TARGSRVY'] for r in rows])
        for s in all_surveys:
            this_survey = survey_clone.cloneNode(True)
            this_survey.setAttribute('name',value=str(s))
            this_survey.setAttribute('priority',value='1.0')
            if mode == 'LIFU':
                this_xml.configure.setAttribute('plate','LIFU')
                this_survey.setAttribute('max_fibres',value='603')
            elif mode == 'mIFU':
                this_xml.configure.setAttribute('plate','PLATE_B')
                this_survey.setAttribute('max_fibres',value='740')
            this_xml.surveys.insertBefore(this_survey,sruveys_comment)
            
        this_xml.obsconstraints.setAttribute('obstemp',str(rows[0]['OBSTEMP']))

        #now the <observation> element:
        if mode == 'LIFU':
            this_xml.observation.setAttribute('name',str(rows[0]['TARGID']))
        elif mode == 'mIFU':
            this_xml.observation.setAttribute('name',str(rows[0]['TARGNAME']))
        this_xml.observation.setAttribute('progtemp',str(rows[0]['PROGTEMP']))
        this_xml.observation.setAttribute('obs_type',str(mode_lookup[rows[0]['PROGTEMP'][0]]))
        this_xml.observation.setAttribute('trimester',str(self.trimester))
        this_xml.observation.setAttribute('pa',str(rows[0]['IFU_PA_REQUEST']))


        
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

        #assembly order:
        #1.target
        #2.field
        #3.fields

        order = first_sci_order
        field = None
        col_names = rows[0].array.columns.names
        for i in xrange(len(rows)):
            row = rows[i]
            target = this_xml.base_target.cloneNode(True)

            # remove things we shouldn't have as the xml gets passed into Configure
            to_remove = ['configid','fibreid']
            for rem in to_remove:
                target.removeAttribute(rem)

            _row = {}
            for col in col_names:
                _row[col] = row[col]
                if (str(row[col]) == 'nan') and numpy.isnan(row[col]):
                    _row[col] = ""
            row = _row
             
            target.setAttribute('targra',value=str(row['GAIA_RA']))
            target.setAttribute('targdec',value=str(row['GAIA_DEC']))
            target.setAttribute('targpmra',value=str(row['GAIA_PMRA']))
            target.setAttribute('targpmdec',value=str(row['GAIA_PMDEC']))
            target.setAttribute('targepoch',value=str(row['GAIA_EPOCH']))
            target.setAttribute('targparal',value=str(row['GAIA_PARAL']))
            target.setAttribute('targid',value=str(row['TARGID']))
            target.setAttribute('targname',value=str(row['TARGNAME']))
            target.setAttribute('targprog',value=str(row['TARGPROG']))
            target.setAttribute('targcat',value=self.input_fits.split('/')[-1])
            target.setAttribute('targprio',value=str(row['TARGPRIO']))
            target.setAttribute('targuse',value='T')
            target.setAttribute('targsrvy',value=str(row['TARGSRVY']))

            if (mode == 'LIFU'):
                #generate a <field> element to add this target to
                field = field_clone.cloneNode(True)
                field.setAttribute('order',value=str(order))
                field.setAttribute('RA_d',value=str(row['GAIA_RA']))
                field.setAttribute('Dec_d',value=str(row['GAIA_DEC']))
                field.appendChild(target)
                #now add this field to the fields element
                fields_clone.appendChild(field)

                order += 1

            elif (mode == 'mIFU'):
                #if this is the first of the targets, then generate the field, otherwise use the existing <field>
                #remember - we require fixed IFU_DITHER patterns for mIFU, so multiple <field> entries not needed
                if i == 0:
                    field = field_clone.cloneNode(True)
                    field.setAttribute('order',value=str(order))
                    field.setAttribute('RA_d',value=str(row['GAIA_RA']))
                    field.setAttribute('Dec_d',value=str(row['GAIA_DEC']))
                field.appendChild(target)

        if (mode == 'mIFU'):
            #now add this field to the fields element
            fields_clone.appendChild(field)
            

        #...and finally to observation
        this_xml.observation.appendChild(fields_clone)


            
        index = 1
        output_name = '%s%s%d.xml'%(self.xml_out,mode.lower(),index)
        while os.path.isfile(output_name):
            index += 1
            output_name = '%s%s%d.xml'%(self.xml_out,mode.lower(),index)

        this_xml.write_xml(output_name)


if __name__ == '__main__':

    import os
    import glob

    input_fits = './input/WC_2020A1-ifu_driver_cat.fits'
    trimester = '2020A1'
    
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
    stage2_ifu.generate_xmls(mifu_mode=1)
    print 
    print('IFU XMLs written to: {0}'.format(output_dir))
