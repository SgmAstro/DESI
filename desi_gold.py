import os
import runtime
import numpy as np
import astropy.units as u
import argparse

from   config              import Configuration
from   findfile            import findfile
from   astropy.coordinates import SkyCoord
from   astropy.table       import Table, vstack, hstack, unique, join
from   ros_tools           import tile2rosette, calc_rosr, ros_limits
from   gama_limits         import gama_field
from   cartesian           import cartesian, rotate
from   cosmo               import cosmo, distmod
from   lss                 import fetch_lss
from   bitmask             import lumfn_mask
from   ddp_zlimits         import ddp_zlimits


def desi_gold(args, survey='sv3', release='fuji'):
    from   desiutil.dust                 import mwdust_transmission
    from   desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask


    raise  