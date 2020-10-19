- Should new survey elements be added for guide/calib stars? To be
  discussed between DM, DT and LPdA.
- Should targprio be overwritten for guide/calib stars? To be discussed
  between DM, DT and LPdA.
- Write the tests.
- Check the position of the LIFU guide camera is right (if it is OK for PA=0,
  rotation should work fine; there is logging.critical call in
  _guide_and_calib_stars.py where the information should be edited).
  See WEAVESPA-537 and WEAVESPA-540 for more details.
- LIFU guide camera is currently circularised to make things much more simple.
  Depending on the orientation in the sky of the guide camera: (i) we could be
  covering all the stars if the major axis of the camera is pointing towards
  the central fibre, (ii) we could be missing small annulus in the sky in the
  search of the guiding star. Anyway, this is only important if we do not find
  a star to know if there could be an useful star in these tiny annuli.
  See WEAVESPA-537 and WEAVESPA-540 for more details.
- Consider to update the maximum size for dithering (c.f. stage 1).

