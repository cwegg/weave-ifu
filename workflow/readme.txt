Documentation for IFU workflow
==============================

v. 20191205

This document details the workflow shown in the relevant html file
included within this repository. Please use this guide in conjunction
with the html file. Hyperlinks therein (depicted as globes) provide
access to the resources described here.

Please note: during development, not all resources are available, and
links will be generated when sample files, code is committed.

Disclaimer: this software is designed to aide in the preparation of
WEAVE IFU observations. We provide no warranty for this software, and
offer support on a best-effort basis.

Definitions:
End-user - the user that is preparing IFU observations
SPA-columns - the mandatory WEAVE SPA columns that must be included in all input FITS catalogues

Stage 0
The end-user must compile an initial list of IFU
targets. Consideration must be made of the instrument mode the
observations are made in, what dither strategy is to be used, what
targets are to be provided.

Stage 1
The information from Stage 0 needs to be standardised into a
preliminary FITS catalogue containing "input IFU drivers": IFU target
positions in each row alongside other critical columns that need to be
propagated through into the XML-based OB files. The example data in
Stage1 provides a template for what quantities must be included, but
should ideally encompass all SPA-columns.

Stage 2
The code at stage 2 analyses the input IFU driver catalogue and
converts this list into a series of distinct pointings (grouping mIFU
pointings into fields where required). This process outputs the
pointings as XML files (using the base XML template), filling the
details as interpreted from the PROGTEMP and OBSTEMP codes.


Stage 3
Each XML file is further developed by injecting the non-science target
information. This includes addition of the guidestar(s) for the LIFU
(mIFU) as well as calibration target options (more than one, to allow
for flexibility when allocating mIFU bundles). The selection of a
guidestar for LIFU also sets the required position angle (PA) at this
point. Users have the option to specify a PA and determine if this is
permitted based on availability of a suitable guiding candidate. For
LIFU pointings with specified preset dither patterns (LIFU_1.xml in
the stage3 example data), we require that the selected guidestar be
visible in each dither. In principle for a custom dither pattern
(LIFU_2.xml) end-users are able to specify a guide star for each
pointing. We recommend LIFU xmls use the same guide star for each
exposure.

Stage 4
The output of Stage 3 can be seen in this container. The files here
are ready to be passed to Configure.

Stage 5
End-users are required to have a copy of Configure installed and
running for this stage. XML files are passed to the Configure
tool. For LIFU observations, this serves solely to populate the XML
with the LIFU array, including correct positions for the fibres. For
the mIFU, this involves a more interatcive method of dropping mIFU
bundles onto input targets. Configure then calculates the bundle
rotation and consequently fibre positions. For the mIFU, this can be
an iterative process, whereby un-allocated targets may be re-submitted
to Configure in order to generate a new XML file (and hence
observation) with these additional targets. For both IFU modes, where
required Configure will implement the pre-set dither pattern. This is
controlled by the apply_dither' attribute in the <dithering> element,
which in turn is inherited from the ifu_dither column in stage 0
catalogue.

Stage 6



Stage 7
Stage 8
Stage 9
Stage 10
Stage 11
Stage 12
