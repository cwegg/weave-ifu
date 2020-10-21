- To modify the XML template to contain binning_X="1" as default value.
  Otherwise, decide an alternative to assign this value. By now, this has to be
  fixed by hand at this stage.
- To develop an OB linking code. In case of no linking, we can simply replace
  "%%%" in the XML files by "". It is worth noting that this should be done once
  the value of binning_X has been assigned.

