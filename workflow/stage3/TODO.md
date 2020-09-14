- Should new survey elements be added for guide/calib stars? To be
  discussed between DM, DT and LPdA.
- Should targprio be overwritten for guide/calib stars? To be discussed
  between DM, DT and LPdA.
- Check the position of the LIFU guide camera is right (for PA=0 is enough,
  rotation should work fine; there is logging.critical call in
  _guide_and_calib_stars.py where the information should be edited).
- LIFU guide camera is currentily circularised to make things much more simple.
- Consider to update the maximum size for dithering (c.f. stage 1).
- Write the tests.

