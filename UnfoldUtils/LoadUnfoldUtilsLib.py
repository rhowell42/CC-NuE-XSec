"""
  LoadUnfoldUtilsLib.py:
   The code necessary to load the libUnfoldUtils.so library
   so that UnfoldUtils C++ objects are available.
   
   This is cribbed directly from the corresponding module
   in Ana/PlotUtils/python/PlotUtils.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    June 2014
"""

import sys
import os

import ROOT
import os,sys
#import PyCintex # needed to properly load the Reflex dictionary

import UnfoldUtils

# the ROOT interpreter gobbles up any arguments
# that might be in sys.argv, whether you like it or not.
# we need to tuck them away for later.
args = sys.argv[:]
sys.argv = sys.argv[:0]

if "UNFOLDUTILSROOT" in os.environ:
	ROOT.gSystem.Load("libUnfoldUtils")

	# copy the classes from the Reflex library
	# into the "UnfoldUtils" namespace for more
	# straightforward access.
	for cls in [ "MnvUnfold", "MnvResponse" ]:
		setattr(UnfoldUtils, cls, getattr(ROOT.MinervaUnfold, cls))
else:
	print >> sys.stderr, "Note: $UNFOLDUTILSROOT is not defined in the current environment.  UnfoldUtils libraries were not loaded."
	
# now that ROOT has done its thing,
# we can restore the arguments...
sys.argv = args
