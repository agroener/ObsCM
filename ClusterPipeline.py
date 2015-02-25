import numpy as np
from xlrd import open_workbook

import ipdb

## The general conversion process will take data in one convention (either 200 or virial), and
## convert these measurements to a common cosmology and uncertainty convention.

## Functions for converting between overdensity convention (if needed) come from MConvert_Personal.py.
import MConvert_Personal as MC
## Functions for converting between cosmologies come from CosmoConvert.py.
import CosmoConvert as CC
## Functions for converting confidence intervals to symmetric uncertainties come from Uncertainties.py.
import Uncertainties as UN

# Get the data form the excel sheet cm_data.xlsx
import DataHandler as DH
clusters,redshift,methods,c200,c200_plus,c200_minus,m200,m200_plus,m200_minus,cvir,cvir_plus,cvir_minus,mvir,mvir_plus,mvir_minus,short_refs,orig_convention,cosmology=DH.startup()
ipdb.set_trace()
