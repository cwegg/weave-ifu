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
lifu_configure --epoch 2020 -i output/WC_2020A1-lifu_01-tgc.xml -o output/WC_2020A1-lifu_01-tgcs-tmp.xml
dither -i output/WC_2020A1-lifu_01-tgcs-tmp.xml -o output/WC_2020A1-lifu_01-tgcs.xml

lifu_configure --epoch 2020 -i output/WC_2020A1-lifu_02-tgc.xml -o output/WC_2020A1-lifu_02-tgcs-tmp.xml
dither -i output/WC_2020A1-lifu_02-tgcs-tmp.xml -o output/WC_2020A1-lifu_02-tgcs.xml

configure --epoch 2020 -f output/WC_2020A1-mifu_01-tgc.xml -o output/WC_2020A1-mifu_01-tgcs-tmp.xml
dither -i output/WC_2020A1-mifu_01-tgcs-tmp.xml -o output/WC_2020A1-mifu_01-tgcs.xml

configure --epoch 2020 -f output/WC_2020A1-mifu_02-tgc.xml -o output/WC_2020A1-mifu_02-tgcs-tmp.xml
dither -i output/WC_2020A1-mifu_02-tgcs-tmp.xml -o output/WC_2020A1-mifu_02-tgcs.xml
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

Stage 5: Creation the IFU FITS catalogue
----------------------------------------

```
stage5/create_ifu_fits_cat.py aux/WC_CatalogueTemplate.fits output/WC_2020A1-*-tgcs.xml --out output/WC_2020A1-ifu_from_xmls.fits --cat_nme1="First Name" --cat_nme2="Surname"
```

Stage 6: Improving the IFU FITS catalogue (optional, but recommended)
---------------------------------------------------------------------

**TBW**

```
cp output/WC_2020A1-ifu_from_xmls.fits output/WC_2020A1-ifu.fits
```

Stage 7: Creation of the combo FITS catalogue
---------------------------------------------

```
stage7/create_combo_fits_cat.py stage7/input/WC_2020A1-mos.fits output/WC_2020A1-ifu.fits --out output/WC_2020A1.fits
```

Stage 8: Submission to WASP of the combo FITS catalogue
-------------------------------------------------------

**TBW**

Stage 9: Finishing the XML files
--------------------------------

```
stage9/fill_xmls_with_fits_info.py output/WC_2020A1.fits output/WC_2020A1-*-tgcs.xml
```

Extra tips
----------

- You can find more help for each script running it with the option `-h`.

