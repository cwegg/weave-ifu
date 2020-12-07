Tutorial for IFU workflow
=========================

Installation
------------

The core code of this workflow has been written in Python 3 using the following
non-standard packages:

- numpy
- matplotlib
- astropy
- astropy_healpix

Additionally, for the plotting tool of stage 1 a Java Runtime Environment is
needed.

In Ubuntu 20.04, these dependencies can be solved installing the following
packages:

- `python3-numpy`
- `python3-matplotlib`
- `python3-astropy` and `astropy-utils`
- `python3-astropy-healpix`
- `openjdk-11-jre`

Once you have the dependencies installed, you simple have to add the directory
which contains the `workflow` folder to the Python path:

```
export PYTHONPATH=`pwd`
cd workflow
```

Stage 1: Creation of IFU driver catalogue
-----------------------------------------

### Stage 1a: Creation of IFU driver catalogue

```
stage1/_create_ifu_driver_cat_example.py
```

### Stage 1b: Checking the IFU driver catalogue

```
stage1/check_ifu_driver_cat.py output/WC_2020A1-ifu_driver_cat.fits
```

### Stage 1c: Plotting the IFU driver catalogue

```
stage1/plot_ifu_driver_cat.py output/WC_2020A1-ifu_driver_cat.fits
```

Stage 2: Creation of the XML files with the targets
---------------------------------------------------

```
stage2/create_xml_files.py output/WC_2020A1-ifu_driver_cat.fits
```

Stage 3: Adding guide and calibration stars to the XML files
------------------------------------------------------------

```
stage3/add_guide_and_calib_stars.py output/WC_2020A1-*-t.xml
```

Stage 4: Configuring the XML files
----------------------------------

```
lifu_configure --epoch 2020 -i output/WC_2020A1-lifu_01-tgc.xml -o output/WC_2020A1-lifu_01-tgcs-tmp-orig.xml
dither -i output/WC_2020A1-lifu_01-tgcs-tmp-orig.xml -o output/WC_2020A1-lifu_01-tgcs-orig.xml

lifu_configure --epoch 2020 -i output/WC_2020A1-lifu_02-tgc.xml -o output/WC_2020A1-lifu_02-tgcs-tmp-orig.xml
dither -i output/WC_2020A1-lifu_02-tgcs-tmp-orig.xml -o output/WC_2020A1-lifu_02-tgcs-orig.xml

configure --epoch 2020 -f output/WC_2020A1-mifu_01-tgc.xml -o output/WC_2020A1-mifu_01-tgcs-tmp-orig.xml
dither -i output/WC_2020A1-mifu_01-tgcs-tmp-orig.xml -o output/WC_2020A1-mifu_01-tgcs-orig.xml

configure --epoch 2020 -f output/WC_2020A1-mifu_02-tgc.xml -o output/WC_2020A1-mifu_02-tgcs-tmp-orig.xml
dither -i output/WC_2020A1-mifu_02-tgcs-tmp-orig.xml -o output/WC_2020A1-mifu_02-tgcs-orig.xml
```

### Extra trick

The other stages of the IFU worflow produce XML files preserving the W3C XML
Canonicalisation (C14N) inherited from the blank XML template and making usage
of a pretty indentation.

If you want to convert the outputs from configure to such convention (which
would be useful for diffing), you could find useful the following command:

```
xmllint --c14n input.xml | xmllint --format --encode utf-8 --output output.xml -
```

Users of Ubuntu 20.04 could be interested in knowing that `xmllint` is available
in the package `libxml2-utils`.

Therefore, this trick applied to the files of the example followed in this
tutorial would lead to the following commands:

```
xmllint --c14n output/WC_2020A1-lifu_01-tgcs-orig.xml | xmllint --format --encode utf-8 --output output/WC_2020A1-lifu_01-tgcs.xml -

xmllint --c14n output/WC_2020A1-lifu_02-tgcs-orig.xml | xmllint --format --encode utf-8 --output output/WC_2020A1-lifu_02-tgcs.xml -

xmllint --c14n output/WC_2020A1-mifu_01-tgcs-orig.xml | xmllint --format --encode utf-8 --output output/WC_2020A1-mifu_01-tgcs.xml -

xmllint --c14n output/WC_2020A1-mifu_02-tgcs-orig.xml | xmllint --format --encode utf-8 --output output/WC_2020A1-mifu_02-tgcs.xml -
```

Stage 5: Creation the IFU FITS catalogue
----------------------------------------

The FITS catalogue template should be downloaded from the WEAVE Operational
Repository ( http://casu.ast.cam.ac.uk/weave/ ). For this tutorial, we have
included the template for the WEAVE Clusters survey:

```
stage5/create_ifu_fits_cat.py aux/WC_CatalogueTemplate.fits output/WC_2020A1-*-tgcs.xml --cat_nme1 "First Name" --cat_nme2 "Surname" --gaia_dr 2
```

Stage 6: Improving the IFU FITS catalogue (optional, but recommended)
---------------------------------------------------------------------

In the most simple case (although it is not recommented at all), you can simply
copy the catalogue of the previous step:

```
cp output/WC_2020A1-ifu_from_xmls.fits output/WC_2020A1-ifu.fits
stage6/check_populated_ifu_fits_cat.py output/WC_2020A1-ifu_from_xmls.fits output/WC_2020A1-ifu.fits
```

Stage 7: Creation of the combo FITS catalogue
---------------------------------------------

```
stage7/create_combo_fits_cat.py stage7/input/WC_2020A1-mos.fits output/WC_2020A1-ifu.fits
```

Stage 8: Submission to WASP of the combo FITS catalogue
-------------------------------------------------------

This step is performed via the WEAVE Automated Submission Platform (WASP,
http://wasp.ast.cam.ac.uk/ ).

Stage 9: Finishing the XML files
--------------------------------

```
stage9/fill_xmls_with_fits_info.py output/WC_2020A1.fits output/WC_2020A1-*-tgcs.xml
```

Stage 10: Submission to WASP of the XML files
---------------------------------------------

This step is performed via the WEAVE Automated Submission Platform (WASP,
http://wasp.ast.cam.ac.uk/ ).

Make sure that your submissions are OK: it is your responsability.

For that goal, you can find useful the figures generated with the following
script:

```
stage10/plot_data_from_xml.py output/WC_2020A1-?ifu_??.xml
```

Extra tips
----------

- You can find more help for each script running it with the option `-h`.

