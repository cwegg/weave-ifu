from . import utils

from . import stage1
from . import stage2
from . import stage3
# Stage 4 is configure
from . import stage5
# Stage 6 has to be developed
from . import stage7
# Stage 8 is WASP


from .stage1 import create_ifu_driver_cat as stage1_main
from .stage1 import check_ifu_driver_cat as stage1_check
from .stage1 import plot_ifu_driver_cat as stage1_plot
from .stage2 import create_xml_files as stage2_main
from .stage3 import add_guide_and_calib_stars as stage3_main
# Stage 4 is configure
from .stage5 import create_ifu_fits_cat as stage5_main
# Stage 6 has to be developed
from .stage7 import create_combo_fits_cat as stage7_main
# Stage 8 is WASP

