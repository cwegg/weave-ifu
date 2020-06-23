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


import argparse
import glob
import logging
import os
import re
import xml.dom.minidom
from collections import OrderedDict

import numpy as np
from astropy.io import fits

from workflow.utils.get_progtemp_info import get_obsmode_from_progtemp


class _XMLdata:

    
    def __init__(self):
        self.xml_template_url = 'http://casu.ast.cam.ac.uk/~dmurphy/opr3/swg/resources/BlankXMLTemplate.xml'
        self.xml_template = 'BlankXMLTemplate.xml'
        if not os.path.isfile(self.xml_template):
            self._wget(self.xml_template_url)
        return


    def new_xml(self,root_data=None):
        #init the new XML
        try:
            dom = xml.dom.minidom.parse(self.xml_template)
        except xml.parsers.expat.ExpatError:
            logging.error('File {} would not parse'.format(self.xml_template))
            raise SystemExit(1)
        self._ingest_xml(dom)
        if root_data:
            for key in root_data.keys():
                self.root.setAttribute(key, value=str(root_data[key]))


    def _wget(self,url,outname=None):

        ##
        ## Put this into a generic tools module
        ##
        
        cmd = 'wget -q -t 1 -T 5 {}'.format(url)
        if outname != None:
            logging.info('Downloading URL {} to {}'.format(url, outname))
            cmd += ' -O {}'.format(outname)
        os.system(cmd)
        return

    
    def _ingest_xml(self,dom):
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
        newxml = self._remove_empty_lines(self.dom.toprettyxml())
        finalxml = self._remove_xml_declaration(newxml)
        logging.info('Writing to {}'.format(filename))
        with open(filename, 'w') as f:
            f.write(finalxml)    
        

    def _remove_xml_declaration(self,xml_text):
        doc = xml.dom.minidom.parseString(xml_text)
        root = doc.documentElement
        xml_text_without_declaration = root.toxml(doc.encoding)
        return xml_text_without_declaration

    
    def _remove_empty_lines(self,xml_text):
        reparsed = xml.dom.minidom.parseString(xml_text)
        return '\n'.join([line for line in reparsed.toprettyxml(indent=' ').split('\n') if line.strip()])


class _IFU:
    """
    Convert the target-level data from an IFU driver catalogue to a set of XMLs.
    
    This class provides code to take an IFU driver catalogue containing
    information of the central fibres and convert it into XML targets that are
    written into a supplied XML.
    
    Parameters
    ----------
    ifu_driver_cat : str
        A FITS file containing an IFU driver cat.
    output_dir : str
        Name of the directory which will containe the output XML files.
    prefix : str, optional
        Prefix to be used in the output files (it will be derived from
        ifu_driver_cat if None is provided).
    suffix : str, optional
        Suffix to be used in the output files.
    """

    
    def __init__(self, ifu_driver_cat, output_dir='.', prefix='', suffix=''):

        # Save the input parameters

        self.ifu_driver_cat = ifu_driver_cat
        self.output_dir = output_dir
        self.prefix = prefix 
        self.suffix = suffix

        # Read and save the information into the IFU driver catalogue
        
        logging.info('Reading {}'.format(self.ifu_driver_cat))

        with fits.open(self.ifu_driver_cat) as hdu_list:

            self.ifu_driver_cat_datamver = hdu_list[0].header['DATAMVER']
            self.trimester = hdu_list[0].header['TRIMESTE']
            self.report_verbosity = hdu_list[0].header['VERBOSE']
            self.author = hdu_list[0].header['AUTHOR']
            self.cc_report = hdu_list[0].header['CCREPORT']

            self.data = hdu_list[1].data

        # Save the data needed for the root element of the XML file
        
        self.root_data = {}
        self.root_data['author'] = self.author
        self.root_data['cc_report'] = self.cc_report
        self.root_data['report_verbosity'] = str(self.report_verbosity)
        
        self.xml_template = 'BlankXMLTemplate.xml'

    
    def _process_ob(self, entry_group):
        
        mode_lookup = {}
        mode_lookup['4'] = 'LIFU'
        mode_lookup['5'] = 'LIFU'
        mode_lookup['6'] = 'LIFU'
        mode_lookup['7'] = 'mIFU'
        mode_lookup['8'] = 'mIFU'
        mode_lookup['9'] = 'mIFU'

        if str(mode_lookup[entry_group[0]['PROGTEMP'][0]]) == 'LIFU':
            obsmode = 'LIFU'
            ndither = len(entry_group)
            if len(entry_group) == 1:
                ndither = entry_group[0]['IFU_DITHER']
            
        else:
            #mIFU
            obsmode = 'mIFU'
            ndither = entry_group[0]['IFU_DITHER']


        this_xml = _XMLdata()
        this_xml.new_xml(root_data=self.root_data)

        #get a clone of the exposures element and wipe it of everything bar the initial calibs
        #make no assumptions about the OB calibration strategy, just take from the template
        exposures = this_xml.exposures.cloneNode(True)
        calibs_after_science = []
        pre_sci = True
        order = 0

        #remove the 'ORDER undefined' comment nodes (as we will define them here)
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
                    #previous must be a text node (<DOM Text node 'u'\n      \n\n '...'>)
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
        for i in range(ndither):
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
        
        all_surveys = np.unique([r['TARGSRVY'] for r in entry_group])
        for s in all_surveys:
            this_survey = survey_clone.cloneNode(True)
            this_survey.setAttribute('name',value=str(s))
            this_survey.setAttribute('priority',value='1.0')
            if obsmode == 'LIFU':
                this_xml.configure.setAttribute('plate','LIFU')
                this_survey.setAttribute('max_fibres',value='603')
            elif obsmode == 'mIFU':
                this_xml.configure.setAttribute('plate','PLATE_B')
                this_survey.setAttribute('max_fibres',value='740')
            this_xml.surveys.insertBefore(this_survey,sruveys_comment)
            
        this_xml.obsconstraints.setAttribute('obstemp',str(entry_group[0]['OBSTEMP']))

        #now the <observation> element:
        if obsmode == 'LIFU':
            this_xml.observation.setAttribute('name',str(entry_group[0]['TARGID']))
        elif obsmode == 'mIFU':
            this_xml.observation.setAttribute('name',str(entry_group[0]['TARGNAME']))
        this_xml.observation.setAttribute('progtemp',str(entry_group[0]['PROGTEMP']))
        this_xml.observation.setAttribute('obs_type',str(mode_lookup[entry_group[0]['PROGTEMP'][0]]))
        this_xml.observation.setAttribute('trimester',str(self.trimester))
        this_xml.observation.setAttribute('pa',str(entry_group[0]['IFU_PA_REQUEST']))


        
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

        this_xml.dithering.setAttribute('apply_dither',str(entry_group[0]['IFU_DITHER']))

        #assembly order:
        #1.target
        #2.field
        #3.fields

        order = first_sci_order
        field = None
        col_names = entry_group[0].array.columns.names
        for i in range(len(entry_group)):
            row = entry_group[i]
            target = this_xml.base_target.cloneNode(True)

            # remove things we shouldn't have as the xml gets passed into Configure
            to_remove = ['configid','fibreid']
            for rem in to_remove:
                target.removeAttribute(rem)

            _row = {}
            for col in col_names:
                _row[col] = row[col]
                if (str(row[col]) == 'nan') and np.isnan(row[col]):
                    _row[col] = ''
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
            target.setAttribute('targcat',value=self.ifu_driver_cat.split('/')[-1])
            target.setAttribute('targprio',value=str(row['TARGPRIO']))
            target.setAttribute('targuse',value='T')
            target.setAttribute('targsrvy',value=str(row['TARGSRVY']))

            if (obsmode == 'LIFU'):
                #generate a <field> element to add this target to
                field = field_clone.cloneNode(True)
                field.setAttribute('order',value=str(order))
                field.setAttribute('RA_d',value=str(row['GAIA_RA']))
                field.setAttribute('Dec_d',value=str(row['GAIA_DEC']))
                field.appendChild(target)
                #now add this field to the fields element
                fields_clone.appendChild(field)

                order += 1

            elif (obsmode == 'mIFU'):
                #if this is the first of the targets, then generate the field, otherwise use the existing <field>
                #remember - we require fixed IFU_DITHER patterns for mIFU, so multiple <field> entries not needed
                if i == 0:
                    field = field_clone.cloneNode(True)
                    field.setAttribute('order',value=str(order))
                    field.setAttribute('RA_d',value=str(row['GAIA_RA']))
                    field.setAttribute('Dec_d',value=str(row['GAIA_DEC']))
                field.appendChild(target)

        if (obsmode == 'mIFU'):
            #now add this field to the fields element
            fields_clone.appendChild(field)
            

        #...and finally to observation
        this_xml.observation.appendChild(fields_clone)


            
        index = 0
        output_path = ''
        
        while (index < 1) or os.path.exists(output_path):
            index += 1
            output_basename = '{}{}_{:02d}{}.xml'.format(
                self.prefix, obsmode.lower(), index, self.suffix)
            output_path = os.path.join(self.output_dir, output_basename)

        this_xml.write_xml(output_path)

        
    def _generate_lifu_xmls(self, lifu_entry_list):

        # How do you group entries belonging to the same OB?
        # (for custom dithers)

        group_id = ('TARGID', 'TARGNAME', 'PROGTEMP', 'OBSTEMP')

        # Detect the non-custom dithers and save them into a list

        non_custom_dithers_list = []
        
        for lifu_entry in lifu_entry_list:

            if lifu_entry['IFU_DITHER'] != -1:

                non_custom_dithers_list.append(lifu_entry)

        # Generate the IFU XMLs of the non-custom dithers

        logging.info(
            'Processing {} non-custom dither OBs for LIFU'.format(
                len(non_custom_dithers_list)))

        for entry in non_custom_dithers_list:

            entry_group = [lifu_entry]

            self._process_ob(entry_group)

        # Detect the custom dithers and save them into a dict

        custom_dithers_dict = OrderedDict()

        for lifu_entry in lifu_entry_list:

            if lifu_entry['IFU_DITHER'] == -1:

                key = tuple(lifu_entry[col] for col in group_id)

                if key not in custom_dithers_dict.keys():
                    custom_dithers_dict[key] = []
                
                custom_dithers_dict[key].append(lifu_entry)

        # Process the custom dither pointing groupings
        
        logging.info(
            'Processing {} custom dither OBs for LIFU'.format(
                len(custom_dithers_dict)))

        for entry_group in custom_dithers_dict.values():
            self._process_ob(entry_group)

                
    def _generate_mifu_xmls(self, mifu_entry_list, mifu_mode='aufbau',
                            mifu_num_calibs=2, mifu_num_extra=0):

        # What happens if there are (e.g.) 100 bundles in a given key and
        # mifu_num_calibs is 2?
        #
        # User needs to decide 3 options:
        # - all: Include all
        # - aufbau: Aufbau principle, i.e. fill up to 18, leaving 2 for
        #           calibration bundles) then make a new XML
        # - equipartition: divide up the bundles equally:
        #                  100 / 18 = 5.5 --> 6 OBs
        #                  100 / 6 = 17 bundles per OB
        #                            (with remainder added to a final OB)
        #                  i.e. 5 OBs x 17 bundles + 1 OB x 15 bundles

        assert mifu_mode in ['aufbau', 'equipartition', 'all']
        assert (mifu_num_calibs >= 0) and (mifu_num_calibs < 20)

        # How do you group bundles belonging to the same field?

        group_id = ('TARGNAME', 'PROGTEMP', 'OBSTEMP', 'IFU_DITHER')

        # Group the mIFU entries per field

        fields_dict = OrderedDict()

        for mifu_entry in mifu_entry_list:

            key = tuple(mifu_entry[col] for col in group_id)

            if key not in fields_dict.keys():
                fields_dict[key] = []

            fields_dict[key].append(mifu_entry)

        # Group the targets in each field per OB

        ob_nested_list = []

        for key in fields_dict.keys():

            # Get the entries in the field

            field_entry_list = fields_dict[key]

            # Guess the number of entries in the field

            num_entries_in_field = len(field_entry_list)

            # Create the group of mIFU entries per OB depending on the mode

            if (mifu_mode == 'aufbau') or (mifu_mode == 'equipartition'):

                # Maximum number of target entries per OB

                max_entries_per_ob = 20 - mifu_num_calibs

                if mifu_mode == 'aufbau':

                    # All the OBs (except the last one) will contain the maximum
                    # number of target entries plus the amount of extra targets
                    
                    num_entries_per_ob = max_entries_per_ob + mifu_num_extra

                elif mifu_mode == 'equipartition':

                    # Compute the number of OBs

                    num_obs = int(np.ceil(num_entries_in_field /
                                          max_entries_per_ob))

                    # Divide the targets to be the same (except the last one,
                    # due to rounding effects)

                    num_entries_per_ob = int(np.ceil(num_entries_in_field /
                                                     num_obs))

                else:

                    raise ValueError

                # Save the targets of each OB using Python slices
                # (do not worry, if stop is higher than the lenght of the list,
                #  this is allowed and the trick)
                
                start = 0

                while start < num_entries_in_field:

                    stop = start + num_entries_per_ob

                    ob_nested_list.append(field_entry_list[start:stop])

                    start = stop
                    
            elif mifu_mode == 'all':

                ob_nested_list.append(field_entry_list)

            else:

                raise ValueError

        # Proccess OB grouping

        logging.info(
            'Processing {} mIFU fields into {} OBs according to {} mode'.format(
                len(fields_dict), len(ob_nested_list), mifu_mode))

        for entry_group in ob_nested_list:
            self._process_ob(entry_group)

                    
    def generate_xmls(self, mifu_mode='aufbau', mifu_num_calibs=2,
                      mifu_num_extra=0):

        # Classify the entries into LIFU and mIFU
        
        lifu_entry_list = []
        mifu_entry_list = []

        for i, entry in enumerate(self.data):

            progtemp = entry['PROGTEMP']

            try:
                obsmode = get_obsmode_from_progtemp(progtemp)
            except:
                obsmode = None

            if obsmode == 'LIFU':
                lifu_entry_list.append(entry)
            elif obsmode == 'mIFU':
                mifu_entry_list.append(entry)
            else:
                logging.warning('unexpected PROGTEMP in row {}: {}'.format(
                    i + 1, progtemp))

        # Generate the LIFU XMLs
                
        self._generate_lifu_xmls(lifu_entry_list)

        # Generate the mIFU XMLs
                
        self._generate_mifu_xmls(mifu_entry_list, mifu_mode=mifu_mode,
                                 mifu_num_calibs=mifu_num_calibs,
                                 mifu_num_extra=mifu_num_extra)


def _get_previous_files(output_dir, prefix, suffix):

    glob_pattern = os.path.join(output_dir, '{}*{}.xml'.format(prefix, suffix))
    regex = '^{}[lm]ifu_[0-9]+{}.xml$'.format(prefix, suffix)

    candidate_list = glob.glob(glob_pattern)
    candidate_list.sort()

    filename_list = [filename for filename in candidate_list
                     if re.match(regex, os.path.basename(filename))]

    return filename_list


def create_xml_files(ifu_driver_cat, output_dir, prefix=None, suffix='-t',
                     mifu_mode='aufbau', mifu_num_calibs=2, mifu_num_extra=0,
                     overwrite=False):
    """
    Create XML files with targets from an IFU driver cat.
    
    Parameters
    ----------
    ifu_driver_cat : str
        A FITS file containing an IFU driver cat.
    output_dir : str
        Name of the directory which will containe the output XML files.
    prefix : str, optional
        Prefix to be used in the output files (it will be derived from
        ifu_driver_cat if None is provided).
    suffix : str, optional
        Suffix to be used in the output files.
    mifu_mode : {'aufbau', 'equipartition', 'all'}, optional
        Grouping mode for mIFU targets.
    mifu_num_calibs : int, optional
        Number of mIFU calibration stars and sky bundles to be considered in
        the mIFU grouping.
    mifu_num_extra : int, optional
        Number of extra mIFU targets to be included in the OBs when 'aufbau'
        grouping mode is used.
    overwrite : bool, optional
        Overwrite the output FITS file.
    """

    # Check that the input IFU driver cat exists and is a file

    assert os.path.isfile(ifu_driver_cat)

    # Set the prefix of the output file if it has not been provided
    
    if prefix is None:
        basename_wo_ext = os.path.splitext(os.path.basename(ifu_driver_cat))[0]

        prefix = basename_wo_ext.replace('-ifu_driver_cat', '')

        if prefix[-1] != '-':
            prefix = prefix + '-'

    # Remove the previous files if overwriting has been requested
    
    if overwrite == True:

        prev_filename_list = _get_previous_files(output_dir, prefix, suffix)

        if len(prev_filename_list) > 0:

            logging.info('Removing previous files: {}'.format(
                prev_filename_list))

            for filename in prev_filename_list:
                os.remove(filename)

    # Create the XML files

    stage2_ifu = _IFU(ifu_driver_cat, output_dir=output_dir,
                      prefix=prefix, suffix=suffix)

    stage2_ifu.generate_xmls(mifu_mode=mifu_mode,
                             mifu_num_calibs=mifu_num_calibs,
                             mifu_num_extra=mifu_num_extra)
    
    logging.info('IFU XMLs written to: {}'.format(output_dir))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
             description='Create XML files with targets from an IFU driver cat')

    parser.add_argument('ifu_driver_cat',
                        help="""a FITS file containing an IFU driver cat""")

    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name of the directory which will containe the
                        output XML files""")

    parser.add_argument('--mifu_mode', dest='mifu_mode', default='aufbau',
                        choices=['aufbau', 'equipartition', 'all'],
                        help="""grouping mode for mIFU targets: 'aufbau'
                        follows the aufbau principle, 'equipartition' divides up
                        the bundles equally, 'all' includes all""")

    parser.add_argument('--mifu_num_calibs', dest='mifu_num_calibs', default=2,
                        choices=range(20), type=int,
                        help="""number of mIFU calibration stars and sky bundles
                        to be considered in the mIFU grouping""")

    parser.add_argument('--mifu_num_extra', dest='mifu_num_extra', default=0,
                        type=int,
                        help="""number of extra mIFU targets to be included in
                        the OBs when 'aufbau' grouping mode is used""")

    parser.add_argument('--prefix', dest='prefix', default=None,
                        help="""prefix to be used in the output files (it will
                        be derived from ifu_driver_cat if non provided)""")

    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='overwrite the output files')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    level_dict = {'debug': logging.DEBUG, 'info': logging.INFO,
                  'warning': logging.WARNING, 'error': logging.ERROR}

    logging.basicConfig(level=level_dict[args.log_level])

    if not os.path.exists(args.output_dir):
        logging.info('Creating the output directory')
        os.mkdir(args.output_dir)

    create_xml_files(args.ifu_driver_cat, args.output_dir, prefix=args.prefix,
                     mifu_mode=args.mifu_mode,
                     mifu_num_calibs=args.mifu_num_calibs,
                     mifu_num_extra=args.mifu_num_extra,
                     overwrite=args.overwrite)

