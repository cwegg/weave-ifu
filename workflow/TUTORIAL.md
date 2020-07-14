Tutorial for IFU workflow
=========================

Stage 1: Creation of IFU driver catalogue
-----------------------------------------

```
stage1/_create_ifu_driver_cat_example.py
stage1/check_ifu_driver_cat.py output/WC_2020A1-ifu_driver_cat.fits
stage1/plot_ifu_driver_cat.py output/WC_2020A1-ifu_driver_cat.fits
```

Stage 2: Creation of the XML files with the targets
---------------------------------------------------

```
stage2/create_xml_files.py input/WC_2020A1-ifu_driver_cat.fits
```

Stage 3: Adding guide and calibration stars to the XML files
------------------------------------------------------------

```
stage3/add_guide_and_calib_stars.py input/WC_2020A1-*-t.xml
```

Stage 4: Configuring the XML files
----------------------------------

TBW

Stage 5: Creation the IFU FITS catalogue
----------------------------------------

```
stage5/create_ifu_fits_cat.py aux/WC_CatalogueTemplate.fits input/WC_2020A1-*-tgcs.xml --out output/WC_2020A1-ifu_from_xmls.fits
```

Stage 6: Improving the IFU FITS catalogue (optional, but recommended)
---------------------------------------------------------------------

TBW

Stage 7: Creation of the combo FITS catalogue
---------------------------------------------

```
stage7/create_combo_fits_cat.py input/WC_2020A1-mos.fits input/WC_2020A1-ifu.fits --out output/WC_2020A1.fits
```

Stage 8: Submission to WASP of the combo FITS catalogue
-------------------------------------------------------

TBW

Stage 9: Finishing the XML files
--------------------------------

TBW

