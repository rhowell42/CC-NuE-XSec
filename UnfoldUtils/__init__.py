# This file, and the fact that the other files here are in the subdirectory UnfoldUtils,
# exist only so that the line 'import UnfoldUtils' will work in other packages.
#
# See http://docs.python.org/2/tutorial/modules.html#packages if you're curious
# how this works.

# load the C++ objects and bind them into the namespace.
from . import LoadUnfoldUtilsLib
